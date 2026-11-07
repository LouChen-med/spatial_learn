#!/usr/bin/env python3
"""Read-only bridge from Jupyter kernels to a Codex-readable event log.

The bridge attaches a second Jupyter client to an existing ipykernel process.
It never sends code to the kernel. Kernel IOPub messages are normalized and
written under the project's tmp/ directory, while the same events are echoed to
stdout so a Codex command session can consume them as they arrive.
"""

from __future__ import annotations

# %%
import argparse
import base64
import fcntl
import json
import os
import queue
import re
import signal
import sys
import threading
import time
import uuid
from datetime import datetime
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Sequence, Tuple

from jupyter_client import BlockingKernelClient


# %%
ANSI_ESCAPE = re.compile(r"\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])")
DEFAULT_LOG_DIR = "tmp/codex_kernel_bridge"
IMAGE_EXTENSIONS = {
    "image/png": "png",
    "image/jpeg": "jpg",
    "image/svg+xml": "svg",
    "application/pdf": "pdf",
}
OUTPUT_MESSAGE_TYPES = {
    "execute_input",
    "execute_result",
    "display_data",
    "update_display_data",
    "stream",
    "error",
    "clear_output",
}


class KernelInfo(NamedTuple):
    pid: int
    cwd: Path
    kernel_name: str
    connection_file: Path

    def as_dict(self) -> Dict[str, object]:
        return {
            "pid": self.pid,
            "cwd": str(self.cwd),
            "kernel_name": self.kernel_name,
            "connection_file": str(self.connection_file),
        }


# %%
def utc_timestamp() -> str:
    return datetime.now().astimezone().isoformat(timespec="milliseconds")


def clean_text(value: object, limit: int) -> str:
    text = value if isinstance(value, str) else repr(value)
    text = ANSI_ESCAPE.sub("", text)
    if len(text) <= limit:
        return text
    omitted = len(text) - limit
    return "{}... [truncated {} chars]".format(text[:limit], omitted)


def connection_file_from_argv(argv: Sequence[str]) -> Optional[Path]:
    for index, argument in enumerate(argv):
        for prefix in ("--f=", "-f=", "--connection-file="):
            if argument.startswith(prefix):
                return Path(argument[len(prefix) :])
        if argument in {"--f", "-f", "--connection-file"}:
            if index + 1 < len(argv):
                return Path(argv[index + 1])
    return None


def read_process_argv(pid: int) -> List[str]:
    raw = Path("/proc/{}/cmdline".format(pid)).read_bytes()
    return [part.decode("utf-8", errors="replace") for part in raw.split(b"\0") if part]


def connection_kernel_name(path: Path) -> str:
    try:
        with path.open("r", encoding="utf-8") as handle:
            metadata = json.load(handle)
        return str(metadata.get("kernel_name") or "unknown")
    except (OSError, ValueError, TypeError):
        return "unknown"


def is_within(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def discover_kernels(project_root: Path) -> List[KernelInfo]:
    kernels: List[KernelInfo] = []
    for process_dir in Path("/proc").iterdir():
        if not process_dir.name.isdigit():
            continue
        pid = int(process_dir.name)
        try:
            argv = read_process_argv(pid)
        except OSError:
            continue
        if not any("ipykernel_launcher" in argument for argument in argv):
            continue

        connection_file = connection_file_from_argv(argv)
        if connection_file is None or not connection_file.is_file():
            continue
        try:
            cwd = Path(os.readlink("/proc/{}/cwd".format(pid))).resolve()
        except OSError:
            continue
        if not is_within(cwd, project_root):
            continue

        kernels.append(
            KernelInfo(
                pid=pid,
                cwd=cwd,
                kernel_name=connection_kernel_name(connection_file),
                connection_file=connection_file.resolve(),
            )
        )
    return sorted(kernels, key=lambda item: item.pid)


def resolve_log_dir(project_root: Path, requested: str) -> Path:
    candidate = Path(requested)
    if not candidate.is_absolute():
        candidate = project_root / candidate
    candidate = candidate.resolve()
    tmp_root = (project_root / "tmp").resolve()
    if not tmp_root.is_dir():
        raise ValueError("project tmp directory does not exist: {}".format(tmp_root))
    try:
        candidate.relative_to(tmp_root)
    except ValueError:
        raise ValueError("log directory must be inside {}".format(tmp_root))
    candidate.mkdir(parents=True, exist_ok=True)
    return candidate


# %%
class EventWriter:
    def __init__(
        self,
        project_root: Path,
        log_dir: Path,
        max_text: int,
        echo: bool = True,
    ) -> None:
        self.project_root = project_root
        self.log_dir = log_dir
        self.images_dir = log_dir / "images"
        self.images_dir.mkdir(parents=True, exist_ok=True)
        self.events_path = log_dir / "events.jsonl"
        self.readable_path = log_dir / "readable.log"
        self.latest_json_path = log_dir / "latest_event.json"
        self.latest_text_path = log_dir / "latest_event.txt"
        self.max_text = max_text
        self.echo = echo
        self.run_id = datetime.now().strftime("%Y%m%dT%H%M%S")
        self.sequence = 0
        self.lock = threading.Lock()
        self.file_lock = (log_dir / ".bridge.lock").open("a+")
        try:
            fcntl.flock(self.file_lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError(
                "another bridge is already writing to {}".format(log_dir)
            )
        self.events_file = self.events_path.open("a", encoding="utf-8", buffering=1)
        self.readable_file = self.readable_path.open(
            "a", encoding="utf-8", buffering=1
        )

    def close(self) -> None:
        with self.lock:
            self.events_file.close()
            self.readable_file.close()
            fcntl.flock(self.file_lock.fileno(), fcntl.LOCK_UN)
            self.file_lock.close()

    def write(self, event: Dict[str, object], human: str) -> None:
        with self.lock:
            self.sequence += 1
            event["run_id"] = self.run_id
            event["sequence"] = self.sequence
            event["timestamp"] = utc_timestamp()
            human = "[{}] seq={} pid={} {}\n{}".format(
                event["timestamp"],
                self.sequence,
                event.get("kernel_pid", "bridge"),
                event.get("msg_type", "event"),
                human,
            )
            serialized = json.dumps(
                event, ensure_ascii=False, separators=(",", ":"), default=str
            )
            self.events_file.write(serialized + "\n")
            self.readable_file.write(human.rstrip() + "\n\n")
            self._atomic_write(self.latest_json_path, serialized + "\n")
            self._atomic_write(self.latest_text_path, human.rstrip() + "\n")
            if self.echo:
                print(human.rstrip(), flush=True)

    def system_event(self, message: str, **fields: object) -> None:
        event = {
            "kernel_pid": fields.pop("kernel_pid", "bridge"),
            "msg_type": fields.pop("msg_type", "bridge"),
            "content": {"message": message, **fields},
        }
        human = message
        if fields:
            human += "\n" + json.dumps(fields, ensure_ascii=False, indent=2)
        self.write(event, human)

    @staticmethod
    def _atomic_write(path: Path, text: str) -> None:
        temporary = path.with_name(path.name + ".tmp")
        temporary.write_text(text, encoding="utf-8")
        os.replace(str(temporary), str(path))


# %%
def save_rich_data(
    writer: EventWriter,
    data: Dict[str, object],
    sequence_hint: str,
) -> Dict[str, object]:
    normalized: Dict[str, object] = {}
    for mime, value in data.items():
        if mime in IMAGE_EXTENSIONS:
            extension = IMAGE_EXTENSIONS[mime]
            image_path = writer.images_dir / "{}-{}.{}".format(
                sequence_hint, re.sub(r"[^A-Za-z0-9_.-]", "_", mime), extension
            )
            try:
                if mime == "image/svg+xml":
                    image_path.write_text(str(value), encoding="utf-8")
                else:
                    image_path.write_bytes(base64.b64decode(str(value)))
            except (OSError, ValueError) as exc:
                normalized[mime] = {
                    "error": "could not save image: {}".format(exc),
                }
            else:
                normalized[mime] = {
                    "path": str(image_path.resolve()),
                    "bytes": image_path.stat().st_size,
                }
            continue
        normalized[mime] = clean_text(value, writer.max_text)
    return normalized


def normalize_message(
    message: Dict[str, object],
    writer: EventWriter,
    kernel: KernelInfo,
) -> Optional[Tuple[Dict[str, object], str]]:
    header = message.get("header") or {}
    parent_header = message.get("parent_header") or {}
    content = message.get("content") or {}
    message_type = str(header.get("msg_type") or "")
    if message_type not in OUTPUT_MESSAGE_TYPES:
        return None

    message_id = str(header.get("msg_id") or uuid.uuid4().hex)
    parent_id = str(parent_header.get("msg_id") or "")
    base = {
        "kernel_pid": kernel.pid,
        "kernel_name": kernel.kernel_name,
        "kernel_cwd": str(kernel.cwd),
        "msg_type": message_type,
        "msg_id": message_id,
        "parent_msg_id": parent_id,
    }

    if message_type == "execute_input":
        code = clean_text(content.get("code", ""), writer.max_text)
        event = {
            **base,
            "content": {
                "execution_count": content.get("execution_count"),
                "code": code,
            },
        }
        return event, "code:\n{}".format(code)

    if message_type == "stream":
        text = clean_text(content.get("text", ""), writer.max_text)
        stream_name = str(content.get("name") or "stdout")
        event = {
            **base,
            "content": {"name": stream_name, "text": text},
        }
        return event, "{}:\n{}".format(stream_name, text)

    if message_type in {"execute_result", "display_data", "update_display_data"}:
        execution_count = content.get("execution_count")
        data = save_rich_data(
            writer,
            content.get("data") or {},
            "{}x{}".format(kernel.pid, message_id[:16]),
        )
        event = {
            **base,
            "content": {
                "execution_count": execution_count,
                "data": data,
                "metadata": content.get("metadata") or {},
            },
        }
        preferred = data.get("text/plain")
        lines = ["{}:".format(message_type)]
        if preferred:
            lines.append(str(preferred))
        for mime, value in data.items():
            if mime == "text/plain":
                continue
            if isinstance(value, dict) and "path" in value:
                lines.append("{} -> {}".format(mime, value["path"]))
            else:
                lines.append("{}: {}".format(mime, clean_text(value, writer.max_text)))
        return event, "\n".join(lines)

    if message_type == "error":
        traceback_text = clean_text(
            "\n".join(str(line) for line in content.get("traceback") or []),
            writer.max_text,
        )
        ename = str(content.get("ename") or "Error")
        evalue = clean_text(content.get("evalue", ""), writer.max_text)
        event = {
            **base,
            "content": {
                "ename": ename,
                "evalue": evalue,
                "traceback": traceback_text,
            },
        }
        return event, "error: {}: {}\n{}".format(ename, evalue, traceback_text)

    event = {**base, "content": content}
    return event, "{}: {}".format(message_type, clean_text(content, writer.max_text))


def stream_kernel(
    kernel: KernelInfo,
    writer: EventWriter,
    stop_event: threading.Event,
    reconnect: bool,
) -> None:
    connected_once = False
    while not stop_event.is_set():
        if not kernel.connection_file.is_file():
            writer.system_event(
                "kernel connection file disappeared; bridge stopped",
                kernel_pid=kernel.pid,
                msg_type="bridge_error",
                connection_file=str(kernel.connection_file),
            )
            return

        client = BlockingKernelClient()
        try:
            client.load_connection_file(str(kernel.connection_file))
            client.start_channels(
                shell=False,
                iopub=True,
                stdin=False,
                hb=False,
                control=False,
            )
            connected_once = True
            writer.system_event(
                "attached to kernel",
                kernel_pid=kernel.pid,
                msg_type="bridge_attached",
                kernel_name=kernel.kernel_name,
                kernel_cwd=str(kernel.cwd),
            )
            next_liveness_check = time.monotonic() + 2.0
            while not stop_event.is_set():
                if time.monotonic() >= next_liveness_check:
                    if not kernel.connection_file.is_file():
                        raise RuntimeError("kernel connection file disappeared")
                    if kernel.pid and not Path("/proc/{}".format(kernel.pid)).exists():
                        raise RuntimeError("kernel process exited")
                    next_liveness_check = time.monotonic() + 2.0
                try:
                    message = client.get_iopub_msg(timeout=0.5)
                except queue.Empty:
                    continue
                normalized = normalize_message(message, writer, kernel)
                if normalized is not None:
                    event, human = normalized
                    writer.write(event, human)
        except Exception as exc:
            writer.system_event(
                "kernel stream error: {}: {}".format(type(exc).__name__, exc),
                kernel_pid=kernel.pid,
                msg_type="bridge_error",
            )
        finally:
            try:
                client.stop_channels()
            except Exception:
                pass

        if not reconnect:
            return
        if stop_event.wait(1.0):
            return
        if connected_once:
            writer.system_event(
                "reconnecting to kernel",
                kernel_pid=kernel.pid,
                msg_type="bridge_reconnect",
            )


# %%
def create_writer(args: argparse.Namespace, project_root: Path) -> EventWriter:
    log_dir = resolve_log_dir(project_root, args.log_dir)
    return EventWriter(
        project_root=project_root,
        log_dir=log_dir,
        max_text=args.max_text,
        echo=not args.quiet,
    )


def selected_kernels(args: argparse.Namespace, project_root: Path) -> List[KernelInfo]:
    if args.connection_file:
        connection_file = Path(args.connection_file).expanduser().resolve()
        if not connection_file.is_file():
            raise ValueError("connection file not found: {}".format(connection_file))
        return [
            KernelInfo(
                pid=args.kernel or 0,
                cwd=project_root,
                kernel_name=connection_kernel_name(connection_file),
                connection_file=connection_file,
            )
        ]

    kernels = discover_kernels(project_root)
    requested = set(args.kernel or [])
    if requested:
        kernels = [kernel for kernel in kernels if kernel.pid in requested]
        missing = requested - {kernel.pid for kernel in kernels}
        if missing:
            raise ValueError(
                "kernel pid(s) not found in project: {}".format(
                    ", ".join(str(pid) for pid in sorted(missing))
                )
            )
    if not kernels:
        raise ValueError(
            "no ipykernel process found under {}".format(project_root)
        )
    return kernels


def command_list(args: argparse.Namespace, project_root: Path) -> int:
    kernels = discover_kernels(project_root)
    if args.json:
        print(json.dumps([kernel.as_dict() for kernel in kernels], indent=2))
        return 0
    if not kernels:
        print("No project ipykernel processes found.")
        return 1
    for kernel in kernels:
        print(
            "pid={pid} kernel={kernel_name} cwd={cwd} connection={connection_file}".format(
                **kernel.as_dict()
            )
        )
    return 0


def command_follow(args: argparse.Namespace, project_root: Path) -> int:
    kernels = selected_kernels(args, project_root)
    writer = create_writer(args, project_root)
    stop_event = threading.Event()
    watch_new = not (
        args.no_reconnect or args.connection_file or args.kernel
    )

    def stop_handler(signum: int, frame: object) -> None:
        stop_event.set()

    previous_sigint = signal.signal(signal.SIGINT, stop_handler)
    previous_sigterm = signal.signal(signal.SIGTERM, stop_handler)
    active: Dict[int, threading.Thread] = {}

    def start_kernel(kernel: KernelInfo) -> None:
        thread = threading.Thread(
            target=stream_kernel,
            name="kernel-{}".format(kernel.pid),
            args=(kernel, writer, stop_event, not args.no_reconnect),
            daemon=True,
        )
        thread.start()
        active[kernel.pid] = thread

    try:
        writer.system_event(
            "bridge started",
            msg_type="bridge_started",
            kernels=[kernel.as_dict() for kernel in kernels],
            events_path=str(writer.events_path),
            readable_path=str(writer.readable_path),
            latest_text_path=str(writer.latest_text_path),
        )
        for kernel in kernels:
            start_kernel(kernel)

        deadline = (
            time.monotonic() + args.duration if args.duration is not None else None
        )
        next_scan = time.monotonic() + 2.0
        while not stop_event.is_set():
            if deadline is not None and time.monotonic() >= deadline:
                break
            for pid, thread in list(active.items()):
                if not thread.is_alive():
                    active.pop(pid, None)
            if watch_new and time.monotonic() >= next_scan:
                for kernel in discover_kernels(project_root):
                    if kernel.pid not in active:
                        start_kernel(kernel)
                next_scan = time.monotonic() + 2.0
            if not active and not watch_new:
                break
            time.sleep(0.2)
    finally:
        stop_event.set()
        for thread in active.values():
            thread.join(timeout=2.0)
        writer.system_event("bridge stopped", msg_type="bridge_stopped")
        writer.close()
        signal.signal(signal.SIGINT, previous_sigint)
        signal.signal(signal.SIGTERM, previous_sigterm)
    return 0


def command_read(args: argparse.Namespace, project_root: Path) -> int:
    log_dir = resolve_log_dir(project_root, args.log_dir)
    events_path = log_dir / "events.jsonl"
    if not events_path.is_file():
        print("No bridge event log found at {}".format(events_path))
        return 1

    recent = queue.Queue(maxsize=args.count)
    with events_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            try:
                event = json.loads(line)
            except ValueError:
                continue
            if recent.full():
                recent.get()
            recent.put(event)

    events = list(recent.queue)
    if not events:
        print("Bridge event log is empty.")
        return 1
    for event in events:
        print(format_stored_event(event))
    return 0


def format_stored_event(event: Dict[str, object]) -> str:
    content = event.get("content") or {}
    message_type = event.get("msg_type")
    if message_type == "execute_input":
        body = "code:\n{}".format(content.get("code", ""))
    elif message_type == "stream":
        body = "{}:\n{}".format(
            content.get("name", "stdout"), content.get("text", "")
        )
    elif message_type in {"execute_result", "display_data", "update_display_data"}:
        lines = ["{}:".format(message_type)]
        for mime, value in (content.get("data") or {}).items():
            if isinstance(value, dict) and "path" in value:
                lines.append("{} -> {}".format(mime, value["path"]))
            else:
                lines.append("{}: {}".format(mime, value))
        body = "\n".join(lines)
    elif message_type == "error":
        body = "error: {}: {}\n{}".format(
            content.get("ename", "Error"),
            content.get("evalue", ""),
            content.get("traceback", ""),
        )
    else:
        body = str(content.get("message") or json.dumps(content, ensure_ascii=False))
    return "[{}] seq={} pid={} {}\n{}".format(
        event.get("timestamp"),
        event.get("sequence"),
        event.get("kernel_pid"),
        message_type,
        body,
    )


# %%
def build_parser() -> argparse.ArgumentParser:
    default_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description=(
            "Read-only Jupyter IOPub bridge for Codex. "
            "All generated files stay under <project>/tmp/codex_kernel_bridge."
        )
    )
    parser.add_argument(
        "--project-root",
        default=str(default_root),
        help="only attach to kernels whose cwd is inside this directory",
    )
    subparsers = parser.add_subparsers(dest="command")
    subparsers.required = True

    list_parser = subparsers.add_parser("list", help="list project kernel processes")
    list_parser.add_argument("--json", action="store_true")

    follow_parser = subparsers.add_parser(
        "follow", help="follow one or more existing kernel IOPub streams"
    )
    follow_parser.add_argument(
        "--kernel",
        action="append",
        type=int,
        help="kernel PID; repeat to select multiple (default: all in project)",
    )
    follow_parser.add_argument(
        "--connection-file",
        help="explicit kernel connection JSON (primarily for isolated tests)",
    )
    follow_parser.add_argument(
        "--log-dir",
        default=DEFAULT_LOG_DIR,
        help="relative to project root and must remain under tmp/",
    )
    follow_parser.add_argument(
        "--max-text",
        type=int,
        default=20000,
        help="maximum characters retained for each text field",
    )
    follow_parser.add_argument(
        "--duration",
        type=float,
        help="stop after this many seconds (default: run until interrupted)",
    )
    follow_parser.add_argument(
        "--no-reconnect",
        action="store_true",
        help="exit a kernel stream after its first channel error",
    )
    follow_parser.add_argument(
        "--quiet",
        action="store_true",
        help="write logs without echoing events to stdout",
    )

    read_parser = subparsers.add_parser(
        "read", help="print the most recent normalized events"
    )
    read_parser.add_argument("--count", type=int, default=20)
    read_parser.add_argument(
        "--log-dir",
        default=DEFAULT_LOG_DIR,
        help="relative to project root and must remain under tmp/",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    project_root = Path(args.project_root).expanduser().resolve()
    if not project_root.is_dir():
        parser.error("project root does not exist: {}".format(project_root))
    max_text = getattr(args, "max_text", None)
    if max_text is not None and max_text <= 0:
        parser.error("--max-text must be positive")
    duration = getattr(args, "duration", None)
    if duration is not None and duration <= 0:
        parser.error("--duration must be positive")
    count = getattr(args, "count", None)
    if count is not None and count <= 0:
        parser.error("--count must be positive")

    try:
        if args.command == "list":
            return command_list(args, project_root)
        if args.command == "follow":
            return command_follow(args, project_root)
        if args.command == "read":
            return command_read(args, project_root)
    except (ValueError, RuntimeError) as exc:
        print("error: {}".format(exc), file=sys.stderr)
        return 2
    parser.error("unknown command")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
