RCMD := Rscript -e

.PHONY: render
render: ## Render OHCA book
	@echo "📖 Rendering OHCA book"
	quarto render --to html

.PHONY: serve
serve: ## serve local static site
	$(RCMD) 'servr::httd("docs", port = 4444)'

.PHONY: render-serve
render-serve: ## Test rendering locally
	@echo "📖 Rendering OHCA book locally"
	quarto render --to html
	$(RCMD) 'servr::httd("docs", port = 4444)'

.PHONY: install
install: ## Install OHCA package and dependencies with pak.
	@echo "🚀 Installing OHCA package"
	$(RCMD) 'install.packages("pak", repos = "https://r-lib.github.io/p/pak/devel/")'
	$(RCMD) 'pak::pkg_install("js2264/HiCExperiment", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("js2264/HiCool", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("js2264/HiContacts", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("js2264/HiContactsData", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("js2264/fourDNData", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install("js2264/DNAZooData", ask = FALSE, dependencies = c("Depends", "Imports", "Suggests"))'
	$(RCMD) 'pak::pkg_install(".", ask = FALSE)'

.PHONY: deps
deps: ## Install missing OHCA dependencies
	@echo "🔗 Installing missing OHCA dependencies"
	$(RCMD) 'devtools::install_dev_deps(".", dependencies = c("Depends", "Imports"))'

.PHONY: git
git: ## Automated commit and pushing to github rpeo
	@echo "📨 Pushing to GH"
	git add .
	git commit -m "make git"
	git push

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.DEFAULT_GOAL := help

