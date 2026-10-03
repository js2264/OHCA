#!/bin/sh
## TEMPORARY, throwaway branch `claude/ohca-devel-test`: render each chapter on
## its own (quarto renders chapters in separate R sessions anyway), carry on past
## failures, and list the chapters that fail, with the end of their log
set -u
cp -r /opt/pkg/inst /tmp/book
mkdir -p /tmp/render-logs
cd /tmp/book
fails=""
## CHAPTERS: chapters to render (default: all, in the book's order)
for f in ${CHAPTERS:-$(grep -oE '(index|pages/[a-z-]+)\.qmd' assets/_book.yml)}; do
    log=/tmp/render-logs/$(basename "${f}" .qmd).log
    start=$(date +%s)
    status=0 ; quarto render "${f}" --to html > "${log}" 2>&1 || status=$?
    grep -h "HARFBUZZ-DIAG" "${log}" | sed 's/^/    | /'
    if [ "${status}" -eq 0 ]; then
        echo "PASS ${f} ($(( $(date +%s) - start )) s)"
    else
        echo "FAIL ${f} ($(( $(date +%s) - start )) s)"
        tail -n 40 "${log}" | sed 's/^/    | /'
        fails="${fails} ${f}"
    fi
done
echo "Chapters that failed:${fails:- none}"
cd / && rm -rf /tmp/book /tmp/render-logs
[ -z "${fails}" ]
