################################################################################
# Description: Setup directory structure for obesity data set   
# - creates directory and subdirectories under ./obesity_results
################################################################################
#
rm(list = ls())

#set location
loc = "obesity"

#
dir.create(paste0(loc, "/obesity_results/NMA"),showWarnings = FALSE,recursive = TRUE)
dir.create(paste0(loc, "/obesity_results/Time"),showWarnings = FALSE,recursive = TRUE)
dir.create(paste0(loc, "/obesity_results/DoseTime"),showWarnings = FALSE,recursive = TRUE)
dir.create(paste0(loc, "/obesity_results/DoseTimeRCOV"),showWarnings = FALSE,recursive = TRUE)
#
