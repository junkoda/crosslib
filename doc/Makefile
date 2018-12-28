all: cross.pdf

%.pdf: %.tex 
	pdflatex $<
	@grep -e 'may have changed' $(basename $<).log && rm $@; $(MAKE) $@ || echo "done."
	test -f copy.sh && sh copy.sh || :

.PHONY: clean check

clean:
	rm *.pdf *.aux *.log *.out

check:
	aspell --lang=en_GB -t check cross.tex
