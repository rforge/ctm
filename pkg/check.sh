
#!/bin/bash

R CMD build variables
R CMD check variables
R CMD INSTALL variables

R CMD build basefun
R CMD check basefun
R CMD INSTALL basefun

R CMD build mlt
R CMD check mlt
R CMD INSTALL mlt

R CMD build ltm
R CMD check ltm
R CMD INSTALL ltm






