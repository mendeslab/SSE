LATEXMK ?= latexmk
DOCS := guide theory explorations
OUT := ../build/docs

.PHONY: all clean $(DOCS)
all: $(DOCS)

# Separate auxiliary directories make independent or parallel document builds safe.
$(DOCS):
	mkdir -p $(OUT)/$@
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error -outdir=$(OUT)/$@ quasse/$@.tex
	cp $(OUT)/$@/$@.pdf $(OUT)/$@.pdf

clean:
	$(foreach name,$(DOCS),$(LATEXMK) -C -outdir=$(OUT)/$(name) quasse/$(name).tex;)
	rm -f $(addprefix $(OUT)/,$(addsuffix .pdf,$(DOCS)))
