RCMD := Rscript -e
SHCMD := /bin/bash

.PHONY: render
render: ## Render OHCA book
	@echo "📖 Rendering OHCA book"
	$(RCMD) 'bookdown::render_book("index.Rmd", "bookdown::gitbook")'

.PHONY: install
install: ## Install OHCA package and dependencies.
	@echo "🚀 Installing package"
	$(RCMD) 'BiocManager::install("js2264/HiCExperiment", dependencies = TRUE)'
	$(RCMD) 'BiocManager::install("js2264/HiContacts", dependencies = TRUE)'
	$(RCMD) 'BiocManager::install("js2264/HiContactsData", dependencies = TRUE)'
	# $(RCMD) 'BiocManager::install("js2264/HiCool", dependencies = TRUE)'
	# $(RCMD) 'BiocManager::install("js2264/fourDNData", dependencies = TRUE)'
	# $(RCMD) 'BiocManager::install("js2264/DNAZooData", dependencies = TRUE)'
	$(RCMD) 'devtools::install()'

.PHONY: deps
deps: ## Install missing OHCA dependencies
	@echo "🔗 Installing missing dependencies"
	$(RCMD) 'devtools::install_dev_deps(".")'

.PHONY: git
git: ## Automated commit and pushing to github rpeo
	@echo "📨 Pushing to GH"
	git add .
	git commit -m "automated commit"
	git push

.PHONY: render-test
render-test: ## Test rendering locally
	@echo "📖 Rendering OHCA book locally"
	$(RCMD) 'bookdown::render_book(".")'
	$(RCMD) 'file.copy("docs/", "/data/scratch/OHCA/", recursive=TRUE)'
	$(RCMD) 'unlink("docs/", recursive=TRUE)'
	$(RCMD) 'unlink("OHCA: Orchestrating Hi-C analysis with Bioconductor.rds")'
	$(RCMD) 'unlink(list.files(".", pattern = ".*rds"))'
	$(RCMD) 'servr::httd("/data/scratch/OHCA/")'

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.DEFAULT_GOAL := help
