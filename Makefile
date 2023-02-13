RCMD := Rscript -e

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

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.DEFAULT_GOAL := help
