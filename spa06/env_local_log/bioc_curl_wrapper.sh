#!/usr/bin/env bash
# Serve verified Bioconductor archives from the project tmp cache.
curl() {
    local cached=""
    case "$*" in
        *GenomeInfoDbData_1.2.13.tar.gz*)
            cached=/home/lcpyspa/spa/tmp/GenomeInfoDbData_1.2.13.tar.gz
            ;;
        *GO.db_3.20.0.tar.gz*)
            cached=/home/lcpyspa/spa/tmp/GO.db_3.20.0.tar.gz
            ;;
        *org.Mm.eg.db_3.20.0.tar.gz*)
            cached=/home/lcpyspa/spa/tmp/org.Mm.eg.db_3.20.0.tar.gz
            ;;
        *org.Hs.eg.db_3.20.0.tar.gz*)
            cached=/home/lcpyspa/spa/tmp/org.Hs.eg.db_3.20.0.tar.gz
            ;;
    esac

    if [[ -n "$cached" ]]; then
        cat "$cached"
        return 0
    fi

    command /usr/bin/curl "$@"
}
