
#!/bin/bash

R CMD BATCH ctm_glm.R 
R CMD BATCH tram_glm.R 

R CMD BATCH ctm_gam.R 
R CMD BATCH tram_gam.R 

R CMD BATCH ctm_tree.R 
R CMD BATCH tram_tree.R 

