RCMD := Rscript -e

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
	$(RCMD) 'write(paste0("R_BIOC_VERSION=", gsub(".[0-9]*$$", "", as.character(packageVersion("BiocVersion")))), paste0(Sys.getenv("R_HOME"), "/etc/Renviron.site"), append = TRUE)'
	$(RCMD) 'pak::repo_status()$$url'
	$(RCMD) 'pak::pkg_install("HiCExperiment", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'remove.packages(c("HiContacts"))'
	$(RCMD) 'pak::pkg_install("HiContacts", upgrade = TRUE, ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("HiCool", upgrade = TRUE, ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("HiContactsData", upgrade = TRUE, ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("fourDNData", upgrade = TRUE, ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("DNAZooData", upgrade = TRUE, ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'

.PHONY: install
install: ## Install OHCA package and dependencies with pak.
	@echo "🚀 Installing OHCA book package"
	$(RCMD) 'pak::pkg_install(".", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'

.PHONY: info
info: ## list installed packages
	@echo "🗂 Listing installed packages"
	$(RCMD) 'pkgs <- installed.packages()[,"Package"] ; sessioninfo::session_info(pkgs, include_base = TRUE)'

.PHONY: render
render: ## Render OHCA book
	@echo "📖 Rendering OHCA book"
	$(RCMD) 'pkgs <- installed.packages()[,"Package"] ; sessioninfo::session_info(pkgs, include_base = TRUE)'
	quarto render --to html

.PHONY: cleanup
cleanup: ## Removing temp files before pushing to ghcr.io
	@echo "🧹 Cleaning up"
	$(RCMD) 'pak::cache_clean()'
	rm -rf /tmp/*
	rm -rf *_cache
	df -h 

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.PHONY: all
all: quarto setup install info cleanup

.DEFAULT_GOAL := all
