
ex = $(shell ls  ex_*.R)
Rda = $(ex:.R=.rda)

.SUFFIXES:
.SUFFIXES: .R .rda

.R.rda:	setup.R
	R CMD BATCH $< &

Rdas: $(Rda)

all: Rdas

summary.pdf:	summary.R
	R CMD BATCH summary.R

dist-clean:	
	rm -rf $(Rda)
	rm -rf *~
	rm -rf *Rout
	
clean:
	rm -rf .RData
	rm -rf *~
	rm -rf summary.Rout
	rm -rf summary.pdf
	