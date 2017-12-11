
	Transformation Forests
	Torsten Hothorn and Achim Zeileis

	Simulation Experiments

Install packages

### from R-forge
### svn checkout -r1638 svn://scm.r-forge.r-project.org/svnroot/coin/pkg/
### R CMD build libcoin
### R CMD INSTALL libcoin_1.0-0.tar.gz
library("libcoin")
### R CMD build inum
### R CMD INSTALL inum_9.3-0.tar.gz
library("inum")
### svn checkout -r1814 svn://scm.r-forge.r-project.org/svnroot/partykit/pkg/
### R CMD build partykit
### R CMD INSTALL partykit_1.2-0.tar.gz
library("partykit")
### R CMD build ATR
### R CMD INSTALL ATR_0.1-0.tar.gz
library("ATR")
### svn checkout -r771 svn://scm.r-forge.r-project.org/svnroot/ctm/pkg/
### R CMD build trtf
### R CMD INSTALL trtf_0.3-0.tar.gz
library("trtf")
### R CMD build mlt
### R CMD INSTALL mlt_0.2-1.tar.gz
library("mlt")

To reproduce the figures, run ``runsim.R''.

