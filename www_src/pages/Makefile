
Rmd = $(shell ls *.Rmd)
md = $(Rmd:.Rmd=.md)

.SUFFIXES:
.SUFFIXES: .Rmd .md

.Rmd.md:
	echo "download.file(\"http://user.math.uzh.ch/hothorn/TH.bib\", dest = \"TH.bib\"); library(\"knitr\"); knit(\"$<\");" | R --vanilla 
	
md: $(md)

clean:
	rm -rf $(md)
	rm -rf TH.bib
	rm -rf *~
	
