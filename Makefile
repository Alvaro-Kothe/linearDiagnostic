RSCRIPT ?= Rscript

check:
	$(RSCRIPT) -e 'devtools::check()'

test:
	$(RSCRIPT) -e 'devtools::test()'

lint:
	$(RSCRIPT) -e 'lintr::lint_package()'

style:
	$(RSCRIPT) -e 'styler::style_pkg()'

.PHONY: test lint style
