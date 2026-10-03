# syntax=docker/dockerfile:1
ARG BIOC_VERSION
FROM bioconductor/bioconductor_docker:${BIOC_VERSION}
COPY . /opt/pkg

## Keep R's user cache (the book's `micromamba` and its conda env) inside the image
ENV R_USER_CACHE_DIR=/opt/R-cache

## Optionally pin `quarto`, e.g. to a release that renders `llms.txt` for books
## (quarto >= 1.11). Empty: keep the quarto shipped with the Bioconductor image
ARG QUARTO_VERSION=""
ARG TARGETARCH
RUN if [ -n "${QUARTO_VERSION}" ]; then \
      curl -fsSL -o /tmp/quarto.tar.gz "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-${TARGETARCH:-amd64}.tar.gz" && \
      mkdir -p "/opt/quarto-${QUARTO_VERSION}" && \
      tar -xzf /tmp/quarto.tar.gz -C "/opt/quarto-${QUARTO_VERSION}" --strip-components=1 && \
      ln -sf "/opt/quarto-${QUARTO_VERSION}/bin/quarto" /usr/local/bin/quarto && \
      rm /tmp/quarto.tar.gz ; \
    fi && \
    quarto --version

## Install book package. GITHUB_PAT arrives as a BuildKit secret, so that 
## GitHub-only dependencies do not hit the anonymous API rate limit, and is 
## never stored in the image
RUN --mount=type=secret,id=GITHUB_PAT \
    GITHUB_PAT="$(cat /run/secrets/GITHUB_PAT 2>/dev/null || true)" \
    Rscript -e 'install.packages("remotes") ; repos <- BiocManager::repositories() ; remotes::install_github("js2264/HiContactsData") ; remotes::install_local(path = "/opt/pkg/", repos=repos, dependencies=TRUE, build_vignettes=FALSE, upgrade=TRUE) ; sessioninfo::session_info(installed.packages()[,"Package"], include_base = TRUE)'

## TEMPORARY (throwaway branch `claude/ohca-devel-test`): HiCExperiment, 
## HiContacts and plyinteractions with the fixes for S4Vectors >= 0.51.10
RUN sh /opt/pkg/ci-devel-test/patch-deps.sh

## Install the book's python environment, if any page runs python, and point
## every R session in the image at it
RUN if grep -rqsE '^```\{python' /opt/pkg/inst; then \
      Rscript -e 'BiocBook::setup_python(requirements = "/opt/pkg/inst/requirements.yml"); cat(sprintf("RETICULATE_PYTHON=%s\n", reticulate::py_config()$python), file = file.path(R.home("etc"), "Renviron.site"), append = TRUE)' && \
      chmod -R a+rX "${R_USER_CACHE_DIR}" && \
      Rscript -e 'stopifnot(nzchar(Sys.getenv("RETICULATE_PYTHON"))); cat("book python:", reticulate::py_config()$python, "\n")' ; \
    fi

## Build/install using same approach than BBS
RUN R CMD INSTALL /opt/pkg

## TEMPORARY (throwaway branch `claude/ohca-devel-test`): render every chapter, 
## and list those that fail
RUN sh /opt/pkg/ci-devel-test/render-chapters.sh
RUN R CMD build --keep-empty-dirs --no-resave-data /opt/pkg
