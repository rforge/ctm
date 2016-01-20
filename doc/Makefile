
all:	mlt.pdf

mlt.tex:	mlt.Rnw
	echo "library(\"knitr\"); knit(\"mlt.Rnw\")" | R --vanilla 
	
mlt.pdf:	mlt.tex
	texi2dvi --pdf mlt.tex
	pdflatex mlt.tex

clean:	
	rm -rf mlt.aux
	rm -rf mlt.tex
	rm -rf mlt.pdf