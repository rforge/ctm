
all:	mlt.pdf tram.pdf

mlt.tex:	mlt.Rnw
	echo "library(\"knitr\"); knit(\"mlt.Rnw\")" | R --vanilla 
#	diff mlt.tex mlt.tex.save

mlt.R:	mlt.Rnw
	echo "library(\"knitr\"); knit(\"mlt.Rnw\", tangle = TRUE)" | R --vanilla 
	grep -v "^opts" mlt.R > tmp.R
	mv tmp.R mlt.R

mlt.Rout:	mlt.R
	R CMD BATCH mlt.R
#	diff mlt.Rout mlt.Rout.save
	
mlt.pdf:	mlt.tex mlt.bib
	texi2dvi --pdf mlt.tex
	pdflatex mlt.tex

tram.tex:	tram.Rnw
	echo "library(\"knitr\"); knit(\"tram.Rnw\")" | R --vanilla 
#	diff mlt.tex mlt.tex.save

tram.R:	tram.Rnw
	echo "library(\"knitr\"); knit(\"tram.Rnw\", tangle = TRUE)" | R --vanilla 
	grep -v "^opts" tram.R > tmp.R
	mv tmp.R tram.R

tram.Rout:	tram.R
	R CMD BATCH tram.R
#	diff mlt.Rout mlt.Rout.save
	
tram.pdf:	tram.tex mlt.bib
	texi2dvi --pdf tram.tex
	pdflatex tram.tex


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
	rm -rf tram.aux
	rm -rf tram.pdf
	rm -rf mlt.tex
	rm -rf mlt.R
	rm -rf mlt.Rout
	rm -rf packages.bib
	rm -rf .Rdata
	rm -rf Rplots.pdf
	rm -rf .RData
	rm -rf tram.tex
	rm -rf tram.R
	rm -rf tram.Rout
	rm -rf packages.bib
	rm -rf figure
	
dist-clean:
	make clean
	rm -rf analysis
	rm -rf cache
	rm -rf figure
	rm -rf DVC.tgz 
	
