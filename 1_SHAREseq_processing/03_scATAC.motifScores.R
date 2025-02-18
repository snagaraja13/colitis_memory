output <- ''
motif_match_file <- '' # peaks x motifs, 0 or 1 annotation for match of each motif in a givenbn peak
unfilt_se_file <- '' # singlet filtered SummarizedExperiment from 02c_scATAC.getSinglet_cluster
bag_file <- '' # grouping of motifs into families (see Supplemental Data), two column matrix: leader motif1,motif2,motif3...
genome <- '' #mm10 or hg38

cores <- 4
nBg <- 250

library(Matrix)
library(SummarizedExperiment)
library(motifmatchr)
library(chromVAR)
library(stringi)
BiocParallel::register(BiocParallel::MulticoreParam(cores, progressbar = TRUE))
if (genome == 'hg19'){
  library(BSgenome.Hsapiens.UCSC.hg19)
  BSg <-BSgenome.Hsapiens.UCSC.hg19
} else if (genome == 'mm10'){
  library(BSgenome.Mmusculus.UCSC.mm10)
  BSg <- BSgenome.Mmusculus.UCSC.mm10
} else if (genome == 'hg38'){
  library(BSgenome.Hsapiens.UCSC.hg38)
  BSg <- BSgenome.Hsapiens.UCSC.hg38
} else {
  cat('Genome not recognized')
}


## read in data and reformat
motif_match <- read.table(motif_match_file,header=T)
dim(motif_match)
motif_ix <- SummarizedExperiment(assays = list(motifMatches = as(as.matrix(motif_match > 0),'sparseMatrix'))) 
saveRDS(motif_ix,paste0(output,'.motif_ix.rds'))

## read in and process count data
se <- readRDS(unfilt_se_file)
se <- addGCBias(se,genome=BSg)
keep_ind <- (Matrix::rowSums(assays(se)$counts) > 0)
se.filt <- se[keep_ind,]
motif_ix.filt <- motif_ix[keep_ind,]

## rescore on motif annotation
set.seed(123)
bg <- getBackgroundPeaks(se.filt,niterations=nBg)
dev_motif <- computeDeviations(object = se.filt,annotations = motif_ix.filt,background_peaks=bg)
saveRDS(dev_motif,paste0(output,'.dev_motif.unbagged.rds'))

assays(dev_motif)$z[1:5,1:5]
z <- assays(dev_motif)$z

## bag deviations / scores
bag_data <- read.table(bag_file,header=T,sep='\t')
leader_list <- lapply(X=1:nrow(bag_data),FUN=function(i){
  bag_members <- c(stri_split_fixed(bag_data[i,2],',',simplify=T))
  if (length(bag_members) > 1){
    SDs <- apply(z[bag_members,],1,sd,na.rm=T)
    bag_leaders <- names(which(SDs == max(SDs)))
    return(bag_leaders[1])
  } else {
    return(bag_members)
  }
})
bagged <- dev_motif[unlist(leader_list),]
saveRDS(bagged,paste0(output,'.dev_motif.bagged.rds'))
