
Rnw = $(shell ls ex_*.Rnw)
R = $(Rnw:.Rnw=.R)
Rda = $(Rnw:.Rnw=.Rda)
tex = $(Rnw:.Rnw=.tex)

.SUFFIXES:
.SUFFIXES: .Rnw .tex .R .Rda

.R.Rda:
	echo "library("utils"); source(\"$<\");" | R --vanilla

.Rnw.tex:
	echo "library("utils"); Sweave(\"$<\");" | R --vanilla 

all: Rdas texs 
	
Rdas: $(Rda)
texs: $(tex) Rdas 
 
clean:
	rm -f $(tex)
	rm -f *pdf
