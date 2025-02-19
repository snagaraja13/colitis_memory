## There is a low-level transcript mixing between cells that happens in SHARE-seq, as with many split pool methods. 
## This code is meant to identify and remove these events
## The cleaned clone call data has also been provided in /ref/ in mouse_organoids.cleaned_clone_calls.txt

output <- ''
clone_call_file <- '' #cell_clone_calls.txt file from 02_SHARE_TRACE, two colums: cell_barcode  cloneA
atac_meta_file <- '' # metadata file for scATAC-seq barcodes, needs columns: "barcode" = scATAC cell barcode, "replicate" = organoid line

## Ad barcode matching, c('P1.barcode','P1.ATAC')
adList <- list()

cores <- 4
clone_thresh <- 5      # minimum cells to keep clone

library(stringi)
library(BuenRTools)
library(pbmcapply)
library(reshape2)
library(ggplot2)
library(dplyr)

##
matchShareSamples <- function(atac_bc, # ATAC barcodes ('R1.XXX,R2.YYY,R3.ZZZ,P1.AAA')
                              rna_bc, # RNA barcodes ('R1.XXX,R2.YYY,R3.ZZZ,P1.BBB')
                              adList # list of P1 matches (list(c('P1.AAA','P1.BBB'),c()...))
){
  d.a <- data.frame(stri_split_fixed(atac_bc,',',simplify=T))
  d.r <- data.frame(stri_split_fixed(rna_bc,',',simplify=T))
  match_list <- lapply(X=adList,FUN=function(pair){
    rna_sub <- d.r %>% filter(X4 == pair[[1]]) %>% mutate(cell = paste0(X1,',',X2,',',X3)) %>% mutate(barcode = paste0(cell,',',X4))
    atac_sub <- d.a %>% filter(X4 == pair[[2]]) %>% mutate(cell = paste0(X1,',',X2,',',X3)) %>% mutate(barcode = paste0(cell,',',X4))
    match_cells <- intersect(rna_sub$cell,atac_sub$cell)
    rna_match <- rna_sub$barcode[match(match_cells,rna_sub$cell)]
    atac_match <- atac_sub$barcode[match(match_cells,atac_sub$cell)]
    return(list(rna_match,atac_match))
  })
  
  r.barcodes <- unlist(lapply(X=match_list,FUN=function(x){return(x[[1]])}))
  a.barcodes <- unlist(lapply(X=match_list,FUN=function(x){return(x[[2]])}))
  
  return(list(a.barcodes,r.barcodes))
}



# load data
clone_calls <- read.table(clone_call_file,sep='\t',header=F)
meta <- read.table(atac_meta_file,sep='\t',header=T)

# filter for barcodes with matching ATAC barcode
matchBC <- matchShareSamples( atac_bc = meta$barcode,
                              rna_bc = clone_calls[,1], # not actually RNA, but giving clonal barcode P1.xx cells here
                              adList = adList
                              )
cat('Total barcode cells:',nrow(clone_calls),'\n')
cat('Total ATAC barcodes:',nrow(meta),'\n')
cat('Shared barcodes between barcode and ATAC:',length(matchBC[[2]]),'\n')
cat(paste0(100*round(length(matchBC[[2]])/nrow(clone_calls),3),'% of barcode cells\n'))
cat(paste0(100*round(length(matchBC[[2]])/nrow(meta),3),'% of ATAC cells\n'))

clone_calls <- clone_calls[match(matchBC[[2]],clone_calls[,1]),]
colnames(clone_calls) <- c('barcode_cell','clone')
clone_calls$atac_cell <- matchBC[[1]]
clone_calls$replicate <- meta$replicate[match(clone_calls$atac_cell,meta$barcode)]
write.table(clone_calls,paste0(output,'.unfiltered_clone_calls.txt'),sep='\t',row.names=F,quote=F)

# clean out bad match clones
# Given the complexity of the barcode space (>10^6 unique sequeces), it is virtually impossible for two cells from different replicates (organoid lines) to get the same clone
# Therefore, we assign each clone to one mouse/organoid line - the one with the majority of cells (usually >95%)
uniq_clones <- unique(unlist(lapply(X=clone_calls$clone,FUN=function(x){return(stri_split_fixed(x,',',simplify=T))})))
clone_calls_cleaned <- clone_calls
for (x  in uniq_clones){
  clone_ind <- (Matrix::rowSums(stri_split_fixed(clone_calls$clone,',',simplify=T) == x) > 0)
  if (sum(clone_ind) == 0){
    cat(x,'\n')
  }
  rep_counts <- table(clone_calls$replicate[clone_ind])
  top_rep <- names(rep_counts)[which(rep_counts == max(rep_counts))]
  if (length(top_rep) > 1){
    top_rep <- top_rep[1]
  }
  remove_clone_ind <- (clone_ind) & (clone_calls$replicate != top_rep)
  clone_calls_cleaned$clone[remove_clone_ind] <- gsub(paste0(x,','),'',clone_calls_cleaned$clone[remove_clone_ind])
  clone_calls_cleaned$clone[remove_clone_ind] <- gsub(paste0(',',x),'',clone_calls_cleaned$clone[remove_clone_ind])
  clone_calls_cleaned$clone[remove_clone_ind & clone_calls_cleaned$clone == x] <- ''
}
clone_calls_cleaned <- clone_calls_cleaned[grepl('clone',clone_calls_cleaned$clone),]
write.table(clone_calls_cleaned,paste0(output,'.cleaned_clone_calls.txt'),sep='\t',row.names=F,quote=F)


# subset to clones with minimum number of cells
clone_counts <- table(unlist(lapply(X=clone_calls_cleaned$clone,FUN=function(x){return(stri_split_fixed(x,',',simplify=T))})))
keep_clones <- names(clone_counts[clone_counts >= clone_thresh])
clone_split <- stri_split_fixed(clone_calls_cleaned$clone,',',simplify=T)
keep_ind <- unlist(lapply(X=1:nrow(clone_split),FUN=function(x){
  return(sum(clone_split[x,] %in% keep_clones) > 0)
}))
clone_calls_filt <- clone_calls_cleaned[keep_ind,]
write.table(clone_calls_filt,paste0(output,'.cleaned_clone_calls.thresh_',clone_thresh,'.txt'),sep='\t',row.names=F,quote=F)
