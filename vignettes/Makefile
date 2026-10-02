## Books need quarto >= 1.11 to render `llms.txt`, and quarto <= 1.8 rejects
## the option, so the `llms` profile (inst/_quarto-llms.yml) is only switched
## on when the installed quarto supports it.
LLMS_PROFILE := $(shell quarto --version | awk -F. '($$1 > 1) || ($$1 == 1 && $$2 >= 11) { print "llms" }')

all: render

render: 
	QUARTO_PROFILE=$(LLMS_PROFILE) quarto render ../inst/
	mkdir -p ../inst/doc && mv ../inst/docs ../inst/doc/book
