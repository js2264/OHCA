## TEMPORARY, throwaway branch `claude/ohca-devel-test`: base::as.data.frame()
## and its callers, on a GInteractions, on Bioconductor devel, with whichever
## InteractionSet is installed. Usage: Rscript interactionset.R <before|after>
stage <- commandArgs(TRUE)[[1]]
suppressPackageStartupMessages(library(InteractionSet))
grDevices::pdf(NULL)
cat("==", stage, "| R", format(getRversion()), 
    "| S4Vectors", format(packageVersion("S4Vectors")), 
    "| InteractionSet", format(packageVersion("InteractionSet")), 
    "| HiCExperiment", format(packageVersion("HiCExperiment")), 
    "| HiContacts", format(packageVersion("HiContacts")), "\n")
bins <- GRanges("chr1", IRanges(seq(1, 91, by = 10), width = 10))
pairs <- expand.grid(i = seq_along(bins), j = seq_along(bins))
pairs <- pairs[pairs$i <= pairs$j, ]
gis <- GInteractions(bins[pairs$i], bins[pairs$j])
gis$count <- seq_along(gis)
calls <- list(
    "base::as.data.frame(gis)" = function() base::as.data.frame(gis),
    "data.frame(gis)" = function() data.frame(gis),
    "tibble::as_tibble(gis)" = function() tibble::as_tibble(gis),
    "as.data.frame(gis), HiCExperiment's generic" = function() {
        get("as.data.frame", asNamespace("HiCExperiment"))(gis)
    },
    "HiContacts::plotMatrix(gis)" = function() ggplot2::ggplotGrob(
        HiContacts::plotMatrix(gis, use.scores = "count", rasterize = FALSE)
    )
)
ok <- vapply(names(calls), function(n) {
    r <- tryCatch({ calls[[n]]() ; "works" }, error = function(e) paste("ERROR:", conditionMessage(e)))
    cat(sprintf("%-45s %s\n", n, r))
    identical(r, "works")
}, logical(1))
if (stage == "after") {
    res <- as.data.frame(testthat::test_dir(
        "/tmp/InteractionSet/tests/testthat", package = "InteractionSet", 
        load_package = "installed", reporter = "summary", stop_on_failure = FALSE
    ))
    cat(sprintf("InteractionSet tests: %d failed, %d errors, %d expectations passed\n", 
        sum(res$failed), sum(res$error), sum(res$passed)))
    if (!all(ok) || sum(res$failed) > 0 || any(res$error)) quit(status = 1)
}
