RSCRIPT ?= Rscript

test:
	$(RSCRIPT) -e 'devtools::test()'
