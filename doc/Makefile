
all:	mlt.pdf

mlt.tex:	mlt.Rnw
	echo "library(\"knitr\"); knit(\"mlt.Rnw\")" | R --vanilla 

mlt.R:	mlt.Rnw
	echo "library(\"knitr\"); knit(\"mlt.Rnw\", tangle = TRUE)" | R --vanilla 
	grep -v "^opts" mlt.R > tmp.R
	mv tmp.R mlt.R

mlt.Rout:	mlt.R
	R CMD BATCH mlt.R
	
mlt.pdf:	mlt.tex todo.tex mlt.bib
	texi2dvi --pdf mlt.tex
	pdflatex mlt.tex

texclean:
	rm -rf mlt.aux
	rm -rf mlt.pdf
	rm -rf *.aux
	rm -rf *.out
	rm -rf *.blg  
	rm -rf *.log  
	rm -rf *.tpt
	rm -rf *.bbl
	rm -rf *.fff
	rm -rf *~

clean:	
	make texclean
	rm -rf mlt.tex
	rm -rf mlt.R
	rm -rf mlt.Rout
	rm -rf packages.bib
	rm -rf .Rdata
	rm -rf Rplots.pdf
	
dist-clean:
	make clean
	rm -rf analysis
	rm -rf cache
	rm -rf figure
	rm -rf DVC.tgz 
	