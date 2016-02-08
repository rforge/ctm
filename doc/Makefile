
all:	mlt.pdf

mlt.tex:	mlt.Rnw
	echo "library(\"knitr\"); knit(\"mlt.Rnw\")" | R --vanilla 

mlt.R:	mlt.Rnw
	echo "library(\"knitr\"); knit(\"mlt.Rnw\", tangle = TRUE)" | R --vanilla 
	grep -v "^opts" mlt.R > tmp.R
	mv tmp.R mlt.R

mlt.Rout:	mlt.R
	R CMD BATCH mlt.R
	
mlt.pdf:	mlt.tex todo.tex mlt.Rnw
	texi2dvi --pdf mlt.tex
	pdflatex mlt.tex

clean:	
	rm -rf mlt.aux
	rm -rf mlt.tex
	rm -rf mlt.pdf
	rm -rf mlt.R
	rm -rf *aux
	rm -rf *out
	rm -rf *~
	rm -rf packages.bib
	
dist-clean:
	make clean
	rm -rf analysis
	rm -rf cache
	rm -rf figure
	