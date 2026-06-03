PAPER = diversityanddissimilarity
BUILD = build

.PHONY: all clean

all: $(PAPER).pdf

$(PAPER).pdf: $(PAPER).tex references.bib Figs/framework.tex Figs/validation.tex Figs/benchmark-jaccard.tex benchmark-results.tex case-study-results.tex
	mkdir -p $(BUILD)
	pdflatex -interaction=nonstopmode -halt-on-error -output-directory=$(BUILD) $(PAPER).tex
	bibtex $(BUILD)/$(PAPER)
	pdflatex -interaction=nonstopmode -halt-on-error -output-directory=$(BUILD) $(PAPER).tex
	pdflatex -interaction=nonstopmode -halt-on-error -output-directory=$(BUILD) $(PAPER).tex
	cp $(BUILD)/$(PAPER).pdf $(PAPER).pdf

clean:
	rm -rf $(BUILD) $(PAPER).pdf
