
tramSrc = $(shell ls  tram*.R)
tramRda = $(tramSrc:.R=.rda)

ctmSrc = $(shell ls  ctm*.R)
ctmRda = $(ctmSrc:.R=.rda)


.SUFFIXES:
.SUFFIXES: .R .rda

.R.rda:	dgp.rda setup.R
	R CMD BATCH $< &

dgp.rda:	dgp.R
	R CMD BATCH dgp.R
        
Rdas: $(tramRda) $(ctmRda)

all: dgp.rda Rdas

stab.tex:	summary.R
	R CMD BATCH summary.R

dist-clean:	
	rm -rf $(tramRda)
	rm -rf $(ctmRda)
	rm -rf dgp.rda
	rm -rf *~
	rm -rf *Rout
	
clean:
	rm -rf *~
	rm -rf stab.tex
	