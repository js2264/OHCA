#!/bin/sh
## TEMPORARY, throwaway branch `claude/ohca-devel-test`: install HiCExperiment,
## HiContacts and plyinteractions from their Bioconductor devel source packages,
## with the fixes of their `claude/new-session-craecj` branches applied, as
## they would be pushed to Bioconductor
set -eu
command -v patch > /dev/null || { apt-get update -qq && apt-get install -y -qq patch > /dev/null ; }
cd /tmp

Rscript -e '
suppressPackageStartupMessages(library(InteractionSet))
cat("S4Vectors", format(packageVersion("S4Vectors")), "| InteractionSet", format(packageVersion("InteractionSet")), "\n")
gi <- GInteractions(GRanges("chr1:1-10"), GRanges("chr1:20-30"))
r <- tryCatch({ base::as.data.frame(gi) ; "works" }, error = function(e) paste("ERROR:", conditionMessage(e)))
cat("Without the fixes, base::as.data.frame(<GInteractions>):", r, "\n")
bioc <- sprintf("https://bioconductor.org/packages/%s/bioc", BiocManager::version())
invisible(download.packages(c("HiCExperiment", "HiContacts", "plyinteractions"), destdir = "/tmp", repos = bioc, type = "source"))
'

for pkg in HiCExperiment HiContacts plyinteractions; do
    tar -xzf ${pkg}_*.tar.gz
    echo "== ${pkg} $(grep -m1 '^Version:' ${pkg}/DESCRIPTION) on Bioconductor devel"
    for p in /opt/pkg/ci-devel-test/${pkg}/*.patch; do
        echo "-- applying $(basename "${p}")"
        patch -d ${pkg} -p1 --forward --no-backup-if-mismatch < "${p}"
    done
    R CMD INSTALL ${pkg} > /tmp/install-${pkg}.log 2>&1 || { tail -n 40 /tmp/install-${pkg}.log ; exit 1 ; }
    echo "-- installed $(Rscript -e "cat(format(packageVersion('${pkg}')))")"
done
rm -rf /tmp/HiCExperiment* /tmp/HiContacts* /tmp/plyinteractions* /tmp/install-*.log

Rscript -e '
stopifnot(identical(get("as.data.frame", asNamespace("HiCExperiment")), BiocGenerics::as.data.frame))
suppressPackageStartupMessages(library(plyinteractions))
gi <- GInteractions(GRanges("chr1:1-10"), GRanges("chr1:20-30"))
print(as_tibble(gi))
'
