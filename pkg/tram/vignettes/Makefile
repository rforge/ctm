
all:	tram.pdf

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
	rm -rf tram.aux
	rm -rf tram.pdf
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
	rm -rf tram.tex
	rm -rf tram.R
	rm -rf tram.Rout
	rm -rf packages.bib
	rm -rf .Rdata
	rm -rf Rplots.pdf
	rm -rf .RData
	
