output <- ''
#raw filtered output from 02_spatial
raw_rna_files <- list('sample1' = '/path/raw.rds',
                     'sample2' = '/path/raw.rds'
)
tumor_call_file <- '' # tumor calls from 03_spatial

ncores <- 16

library(Matrix)
library(dplyr)
library(stringi)
library(pbmcapply)

##
d_tumor <- read.table(tumor_call_file,sep='\t',header=T)
uniq_samp <- unique(d_tumor$sample)

## pseudobulk counts
ps_list <- list()
for (samp in uniq_samp){
  d <- d_tumor %>% filter(sample == samp) %>% filter(isTumor == 'yes') %>% mutate(tumorID = paste0(sample,'_tumor',cluster))
  uniq_tumors <- unique(d$tumorID)
  rna_all <- readRDS(raw_rna_files[[samp]])
  raw_counts <- (rna_all$counts)[,stri_split_fixed(rownames(d),'.',simplify=T)[,3]]
  rm(rna_all)
  gc()
  identical(colnames(raw_counts),stri_split_fixed(rownames(d),'.',simplify=T)[,3])
  
  ps_samp_list <- pbmclapply(X=1:nrow(raw_counts),FUN=function(i){
    return(cbind(d,counts = raw_counts[i,]) %>% group_by(tumorID) %>% summarise(total_counts = sum(counts)) %>% arrange(factor(tumorID,levels = uniq_tumors)) %>% pull(total_counts))
  },mc.cores=ncores)
  ps_samp <- do.call('rbind',ps_samp_list)
  colnames(ps_samp) <- uniq_tumors
  rownames(ps_samp) <- rownames(raw_counts)
  
  ps_list[[samp]] <- as(ps_samp,'sparseMatrix')
}
ps_raw <- do.call('cbind',ps_list)
saveRDS(ps_raw,paste0(output,'.tumor_pseudobulk.raw.rds'))

## this pseudobulk matrix can now be scored for AP-1/P20 program expression using same code as single cells
## a t-test was then run on colitis recovered vs control tumors
