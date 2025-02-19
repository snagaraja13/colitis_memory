fastq <- '' # R1 fastq for barcode sequences extracted by 01_SHARE_TRACE
ncores <- 4
min_bc_length <- 40   # exclude barcodes shorter than this, will depend on sequencing forma
max_bc_length <- 48   # exclude barcodes longer than this after trimming
min_read <- 5         # minimum raw sequencing reads to include a putative barcode sequence, will depend on sequencing depth but meant to exclude sequencing error
lv_dist_max <- 4     # max Levenshtein distance to collapse two barcodes to the same clone
removePG <- T         # remove poly-G UMIs

## options for using a mixture model to UMI number for cutoffs
## This was not used for analysis in this paper but can be used for more stringent clone calls for very deep datasets with high rates of clone capture
useModel <- F         # use mixture model to call clones from barcodes
Nfit <- 100         # number of iterations for mixture model
logUMI <- T         # whether model off of logUMIs to raw UMIs



source('../tools/static_barcode.functions.v1.3.R')

cat('Reading in barcode file...\n')
split_prefix <- stri_split_fixed(fastq,'/',simplify=T)[length(split_prefix <- stri_split_fixed(fastq,'/',simplify=T))]
prefix <- stri_split_fixed(split_prefix,'.matched',simplify=T)[1]

f1 <- gzfile(fastq,'rt')
d_all <- readLines(f1)
close(f1)
tags  <- d_all[seq(1,length(d_all),by=4)]
reads <- d_all[seq(2,length(d_all),by=4)]
num.reads <- length(tags)
cat('> Number of unfiltered barcode reads detected:',num.reads,'\n')

d <- stri_split_fixed(tags,'_',simplify=T)
d <- d[,(ncol(d)-1):ncol(d)]
d <- data.frame(cbind(d,reads))
colnames(d) <- c('cell','UMI','read')
if (removePG){
  d <- d %>% filter(!(UMI %in% c('GGGGGGG','GGGGGGGG','GGGGGGGGG','GGGGGGGGGG')))
}
rm(d_all,tags,reads)
gc()

####
#### 01 - Filtering raw reads
####

cat('Filtering barcodes\n')
max_read_char <- max(nchar(d$read))
if (max_read_char < min_bc_length){
  cat('-- Max barcode sequence length detected as',max_read_char,'bases. This is shorter than the minimum length set',paste0('(',min_bc_length,')'),'\n')
  cat('-- Resetting minimum barcode length to',max_read_char,'\n')
  min_bc_length <- max_read_char
}
cat('- Validating barcode sequences\n')
ind_pf1 <- validateBarcodeSeq(d$read, min_length = min_bc_length, max_length = max_bc_length, ncores = ncores)
pf1 <- d[ind_pf1,]
pf1$read <- substr(pf1$read,1,max_bc_length)
rm(d,ind_pf1)
gc()

cat('- Removing cell-UMI-barcode triples with low read counts\n')
bc_trips <- apply(pf1,1,paste,collapse='_')
read_counts <- data.frame(table(bc_trips))
hist_reads <- hist(log10(read_counts$Freq),plot=F)

pdf(paste0(prefix,'.static_bc.all_valid.reads_per_umi.pdf'),width=6,height=4)
with(hist_reads,plot(mids,log10(counts),xlab='log10( # of reads )',ylab='log10( # of doubles )',main='Reads per UMI-Barcode Double'))
with(hist_reads,lines(mids,log10(counts)))
abline(v=log10(min_read),lty=2,col='red')
dev.off()

read_filt <- stri_split_fixed(as.character(read_counts$bc_trips[read_counts$Freq>=min_read]),'_',simplify=T)[,3]
pf2 <- pf1[(pf1$read %in% read_filt),]
cat('> Barcodes with minimum',min_read,'reads:',nrow(pf2),'/',nrow(pf1),paste('(',round(100*nrow(pf2)/nrow(pf1),1),'%)\n',sep=''))
rm(read_counts,read_filt,pf1,bc_trips)
gc()

bc_counts <- data.frame(table(pf2$read))
bc_counts <- bc_counts[order(bc_counts$Freq,decreasing=T),]
write.table(bc_counts,paste0(prefix,'.static_bc.min_reads_',min_read,'.barcode_counts.csv'),quote=F,sep=',',row.names=F)
pdf(paste0(prefix,'.static_bc.min_reads_',min_read,'.rep_dist.pdf'),width=6,height=6)
plot(1:nrow(bc_counts),cumsum(bc_counts$Freq)/sum(bc_counts$Freq),pch=16,cex=0.5,
     main='Raw barcode distribution',
     xlab='Barcodes ranked by abundance',
     ylab='Fraction of total reads')
dev.off()

####
#### 02 - Collapsing barcodes to clones
####

cat('Collapsing barcode sequences to clones...\n')
uniq_bc <- unique(bc_counts[,1])
cat('> Number of static barcode sequences:',length(uniq_bc),'\n')
dist_list <- pbmclapply(X=1:length(uniq_bc),FUN=function(x){return(t(adist(uniq_bc[x],uniq_bc)))},mc.cores=ncores)
all_dist <- unlist(dist_list)
all_dist <- all_dist[all_dist>0]
pdf(paste0(prefix,'.static_bc.pairwise_dist.pdf'),width=6,height=4)
hist(all_dist,main='Pairwise barcode distances',xlab='Levenshtein distance',freq=F)
abline(v=lv_dist_max,lty=2,col='red')
dev.off()

pdf(paste0(prefix,'.static_bc.lev_collapse_',lv_dist_max,'.pdf'),width=6,height=4)
bc_clones <- levBarcodeCollapse(uniq_bc,lv_dist_max,ncores,plot=T,prompt=F)
dev.off()
write.table(bc_clones,paste0(prefix,'.static_bc.lv_dist_',lv_dist_max,'.collapsed_clones.txt'),sep=',',row.names=F,quote=F)

####
#### 03 - Filter clone calls
####

cat('- Filtering by UMIs per cell-clone pairing...\n')
pf2$clone <- paste0('clone',bc_clones$clone[match(pf2$read,bc_clones$bc)])
clone_trips2 <- apply(pf2[,c(1,2,4)],1,paste,collapse='_')
clone_collap <- unique(clone_trips2)
rm(clone_trips2)
cat('> Number of unique cell-UMI-clone triples:',length(clone_collap),'\n')

trip_collap_split <- stri_split_fixed(clone_collap,'_',simplify=T)
doub_collap <- paste(trip_collap_split[,1],trip_collap_split[,3],sep='_')
rm(trip_collap_split,clone_collap)
umi_counts <- data.frame(table(doub_collap))

if (useModel){
  cat('-- Running mixture model...\n')
  pdf(paste0(prefix,'.static_bc.model_fit.pdf'),width=7,height=5)
  fit <- fitCompModel(umi_counts$Freq,Nfit,logUMI)
  dev.off()
  write.table(fit[[1]],paste0(prefix,'.model_iter.csv'),quote=F,sep=',',row.names=F)
  min_umi <- 2^fit[[3]]
} else {
  min_umi <- 1
}
cat('> Minimum UMIs for cell-barcode double:',min_umi,'\n')
write.table(umi_counts[umi_counts$Freq >= min_umi,],paste0(prefix,'.static_bc.min_umi_',min_umi,'.barcode_counts.csv'),quote=F,sep=',',row.names=F)
hist_umis <- hist(log10(umi_counts$Freq),plot=F)

pdf(paste0(prefix,'.static_bc.umi_filter.umi_per_double.pdf'),width=6,height=4)
with(hist_umis,plot(mids,counts,xlab='log10( # of UMIs )',ylab='# of doubles',main='UMIs per Cell-Barcode Double'))
with(hist_umis,lines(mids,counts))
abline(v=log10(min_umi),lty=2,col='red')
dev.off()

doub_filt <- doub_collap[(doub_collap %in% umi_counts$doub_collap[(umi_counts$Freq >= min_umi)])]
pairs_filt <- unique(doub_filt)
cat('> Cell-barcode doubles passing threshold:',length(pairs_filt),'/',length(unique(doub_collap)),paste('(',round(100*length(pairs_filt)/length(unique(doub_collap)),1),'%)\n',sep=''),'\n')
rm(pf2)
gc()

counts_filt <- umi_counts[(umi_counts$Freq >= min_umi),]

####
#### 03 - Call final clones
####

cat('- Creating cell - clone assignments\n')
pairs <- data.frame(cell = stri_split_fixed(counts_filt[,1],'_',simplify=T)[,1],
                    clone = stri_split_fixed(counts_filt[,1],'_',simplify=T)[,2]
                    )
pairs$UMIs <- counts_filt$Freq
pairs$clone_pair <- paste(pairs$cell,pairs$clone,sep='_')
umis_by_clone <- data.frame(pairs %>% group_by(clone_pair) %>% summarise(UMIs = sum(UMIs)))
umis_by_clone$cell <- stri_split_fixed(umis_by_clone$clone_pair,'_',simplify=T)[,1]
umis_by_clone$clone <- stri_split_fixed(umis_by_clone$clone_pair,'_',simplify=T)[,2]
rm(pairs)

uniq_clones <- c(unique(umis_by_clone[,4]))
uniq_cells <- c(unique(umis_by_clone[,3]))
clone_count_mat <- sparseMatrix(i = match(umis_by_clone[,4],uniq_clones),
                                j = match(umis_by_clone[,3],uniq_cells),
                                x = umis_by_clone$UMI
                                )
dimnames(clone_count_mat) <- list(uniq_clones, uniq_cells)
static_call_list <- lapply(X=1:ncol(clone_count_mat),FUN=function(x){
  cell <- colnames(clone_count_mat)[x]
  clones <- paste(rownames(clone_count_mat)[(clone_count_mat[,x] > 0)],collapse=',')
  return(c(cell,clones))
})
static_calls <- do.call('rbind',static_call_list)

cat('> Cells with static barcode assignments:',ncol(clone_count_mat),'\n')
cat('> Cells with greater than 1 barcode assigned:',sum(Matrix::colSums(clone_count_mat>0) > 1),paste('(',round(100*sum(Matrix::colSums(clone_count_mat>0) > 1)/ncol(clone_count_mat),1),'%)\n',sep=''),'\n')
cat('> Clones with more than 5 cells assigned:',sum(Matrix::rowSums(clone_count_mat > 0) > 5),'/',nrow(clone_count_mat),'\n')
write.table(static_calls,paste0(prefix,'.static_bc.min_umi_',min_umi,'.lv_dist_',lv_dist_max,'.cell_clone_calls.txt'),sep='\t',col.names=F,row.names=F,quote=F)
saveRDS(clone_count_mat,paste0(prefix,'.static_bc.min_umi_',min_umi,'.lv_dist_',lv_dist_max,'.cell_clone_calls.umi_counts.rds'))

#clone sizes
cells_per_static <- data.frame(clone=rownames(clone_count_mat),cells=Matrix::rowSums(clone_count_mat > 0))
pdf(paste0(prefix,'.static_bc.min_umi_',min_umi,'.lv_dist_',lv_dist_max,'.cells_per_clone.pdf'),width=6,height=4)
hist(log10(cells_per_static[,2]),breaks=30,main='Cells per clone',xlab='log10(cells per static barcode)')
plotCloneSizes(cells_per_static)
dev.off()

cat('Done!\n')

