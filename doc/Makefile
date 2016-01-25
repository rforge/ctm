
all:	mlt.pdf

mlt.tex:	mlt.Rnw
	echo "library(\"knitr\"); knit(\"mlt.Rnw\")" | R --vanilla 

mlt.R:	mlt.Rnw
	echo "library(\"knitr\"); knit(\"mlt.Rnw\", tangle = TRUE)" | R --vanilla 
	grep -v "^opts" mlt.R > tmp.R
	mv tmp.R mlt.R

	
mlt.pdf:	mlt.tex
	texi2dvi --pdf mlt.tex
	pdflatex mlt.tex

clean:	
	rm -rf mlt.aux
	rm -rf mlt.tex
	rm -rf mlt.pdf
	rm -rf mlt.R
	rm -rf *~
	rm -rf packages.bib
	
dist-clean:
	make clean
	rm -rf cache
	rm -rf figure
	