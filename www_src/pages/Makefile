
Rmd = $(shell ls *.Rmd)
md = $(Rmd:.Rmd=.md)

.SUFFIXES:
.SUFFIXES: .Rmd .md

.Rmd.md:
	rm -rf figure/*
	echo "download.file(\"http://user.math.uzh.ch/hothorn/TH.bib\", dest = \"TH.bib\"); library(\"knitr\"); knit(\"$<\");" | R --vanilla 
	cp -ra figure/* ../assets/img/
	
md: $(md)

