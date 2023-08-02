RCMD := Rscript -e
CONDA := /root/.cache/R/basilisk/1.12.1/0/bin/conda

.PHONY: quarto
quarto: ## Update to latest available quarto
	@echo "💾 Installing latest Quarto CLI"
	wget https://quarto.org/docs/download/_download.json;\
	ver=$$(grep -o '"version": "[^"]*' _download.json | grep -o '[^"]*$$');\
	wget https://github.com/quarto-dev/quarto-cli/releases/download/v"$${ver}"/quarto-"$${ver}"-linux-amd64.deb;\
	sudo dpkg -i quarto-"$${ver}"-linux-amd64.deb;\
	rm quarto-"$${ver}"-linux-amd64.deb;\
	rm _download.json;\
	quarto --version;\

.PHONY: setup
setup: ## Install HiCExperiment & co packages with pak.
	@echo "📦 Installing OHCA core packages"
	$(RCMD) 'install.packages("pak", repos = "https://r-lib.github.io/p/pak/devel/")'
	$(RCMD) 'writeLines(paste0("R_BIOC_VERSION=", gsub(".[0-9]*$$", "", as.character(packageVersion("BiocVersion")))), ".Renviron")'
	$(RCMD) 'pak::repo_status()$$url'
	$(RCMD) 'pak::pkg_install("HiCExperiment", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("HiCool", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("HiContacts", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("HiContactsData", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("fourDNData", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("DNAZooData", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'

.PHONY: install
install: ## Install OHCA package and dependencies with pak.
	@echo "🚀 Installing OHCA book package"
	$(RCMD) 'remotes::install_local(".", ask = FALSE, upgrade = TRUE)'

.PHONY: info
info: ## list installed packages
	@echo "🗂 Listing installed packages"
	$(RCMD) 'pkgs <- installed.packages()[,"Package"] ; sessioninfo::session_info(pkgs, include_base = TRUE)'

.PHONY: render
render: ## Render OHCA book
	@echo "📖 Rendering OHCA book"
	quarto render --to html

.PHONY: cleanup
cleanup: ## Removing temp files before pushing to ghcr.io
	@echo "🧹 Cleaning up"
	$(RCMD) 'pak::cache_clean()'
	rm -rf /tmp/*
	rm -rf *_cache
	df -h 

.PHONY: .serve
serve: ## serve local static site
	$(RCMD) 'servr::httd("docs", port = 4444)'

.PHONY: .render-serve
render-serve: ## Test rendering locally
	@echo "📖 Rendering OHCA book locally"
	quarto render --to html
	$(RCMD) 'servr::httd("docs", port = 4444)'

.PHONY: .deps
deps: ## Install missing OHCA dependencies
	@echo "🔗 Installing missing OHCA dependencies"
	$(RCMD) 'devtools::install_dev_deps(".", dependencies = c("Depends", "Imports", "Suggests"))'

.PHONY: .git
git: ## Automated commit and pushing to github rpeo
	@echo "📨 Pushing to GH"
	git add .
	git commit -m "make git"
	git push

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.PHONY: all
all: quarto setup install info cleanup

.DEFAULT_GOAL := all
