args <- commandArgs(trailingOnly=T)
prefix <- args[1]
genome <- args[2]

library(BuenRTools)
library(data.table)

chromSizesFile <- paste0('/mnt/users/snagaraja/ref/ref_sizes/',genome,'genome.bed')

if (genome == 'hg19'){
	myBlackListFile <- '/mnt/Apps/JDB_tools/01_additionalData/BlacklistFiles/hg19_blacklist.JDB.bed'
} else if (genome == 'mm10'){
	myBlackListFile <- '/mnt/Apps/JDB_tools/01_additionalData/BlacklistFiles/mm10_blacklist.JDB.bed'
} else { myBlackListFile <- '~/comp_tools/Split_seq/null_blacklist.bed' }
summitsToCleanPeaks(summitFiles = paste0(prefix,'.macs2_summits.bed'),
                    peakWidth = 800,
                    blackList = myBlackListFile,
                    chromSizes =  chromSizesFile,
                    topNPeaks = NULL, # Filter top N peaks after?
                    useQuantileRanks = FALSE, # Use normalized ranks instead of raw MACS score?
                    FDRcutoff = 0.01, # MACS peak FDR cut-off used for all peaks
                    outDir = "./", # Make sure this directory exists first
                    expName= paste0(prefix,".cleanPeaks"), # Name to give as prefix to resulting filtered peak bed file
                    resizeTo = 300 # Re-size peaks after filtering?
                    )
