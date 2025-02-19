## This code is for differential gene expression in the in vivo memory dataset
## The adenoma vs stem differential was performed analogously between cell types rather than between colitis stage and control

output <- ''
raw_file <- '' # genes x cells, raw UMI counts
meta_file <- '' # metadata file, needs columns: "barcode" for each cell barcode, "celltype" for stem_prog, etc, "replicate" for unique identifier for each mouse

## ATAC-RNA mapping to find overlapping cells
## c('P1.RNA','P1.ATAC'), etc
adList <- list()

ncores <- 16

library(dplyr)
library(Matrix)
library(pbmcapply)
library(stringi)
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

################################################
#### 1) Pseudobulk by cell type and replicate
################################################

## filter data
all_rna <- readRDS(raw_file)
meta <- read.table(meta_file,sep='\t',header=T) %>% filter(celltype %in% c('stem_prog','AE_inter','AE_diff')) ## subset to cell types you want to test

matchBC <- matchShareSamples(atac_bc = meta$barcode,
                             rna_bc = colnames(all_rna),
                             adList = adList
                             )
length(matchBC[[1]])
rna_sub <- all_rna[,matchBC[[2]]]
meta_sub <- meta[match(matchBC[[1]],meta$barcode),]
dim(rna_sub)
dim(meta_sub)

rm(all_rna)
gc()

## pseudobulk - can take a few hours 
uniq_samples <- meta_sub %>% group_by(replicate,celltype) %>% mutate(sample2 = paste0(replicate,'+',celltype)) %>% pull(sample2) %>% unique()
ps_list <- pbmclapply(X=1:nrow(rna_sub),FUN=function(i){
  d <- cbind(meta_sub,counts = rna_sub[i,]) %>% group_by(replicate,celltype) %>% summarise(total_counts = sum(counts)) %>% mutate(sample2 = paste0(replicate,'+',celltype))
  return(d$total_counts[match(uniq_samples,d$sample2)])
},mc.cores=ncores)
ps <- do.call('rbind',ps_list)
dim(ps)
colnames(ps) <- uniq_samples
rownames(ps) <- rownames(rna_sub)
saveRDS(as(ps,'sparseMatrix'),paste0(output,'.pseudobulk.raw.rds'))


################################################
#### 2) Differential testing
################################################
library(DESeq2)

### organize data
## example is given for in vivo tissue memory, testing each stage vs control. 
## For other datasets (e.g. adenoma vs stem), will need to set up data table, conditions, etc accordingly

sample_meta <- meta %>% group_by(replicate) %>% summarise(condition = unique(sample))
info <- data.frame(stri_split_fixed(colnames(ps),'+',simplify=T))
colnames(info) <- c('replicate','cellType')
info$condition <- meta$condition[match(info$replicate,meta$replicate)]

# test by condition and cell type
conds_to_test <- c('Colitis_acute','Colitis_chronic','Colitis_recovered')
uniq_groups <- unique(info$cellType)

stat_list <- list()
for (i in 1:length(uniq_groups)){
  gp <- uniq_groups[i]
  keep_ind <- which(info$cellType == gp)
  counts_gp <- ps[,keep_ind]
  info_gp <- info[keep_ind,]
  
  diff_list <- pbmclapply(X=conds_to_test,FUN=function(cond){
    keep_ind2 <- which(info_gp$condition %in% c(cond,'Control'))
    counts_sub <- counts_gp[,keep_ind2]
    info_sub <- info_gp[keep_ind2,]
    
    # filter for min RPM 10 - can also test all genes if preferred, like was done for adenoma vs stem
    total_counts <- Matrix::colSums(counts_sub)
    rpm_scale <- 1000000/total_counts
    rpm_list <- lapply(X=1:ncol(counts_sub),FUN=function(j){return(counts_sub[,j]*rpm_scale[j])})
    rpm <- do.call('cbind',rpm_list)
    
    keep_ind3 <- apply(rpm,1,max,na.rm=T) > 10
    
    # test
    dds <- DESeqDataSetFromMatrix(countData = counts_sub[keep_ind3,],
                                  colData = info_sub,
                                  design = ~ condition
    )
    dds <- DESeq(dds)
    res <- data.frame(results(dds,contrast=c("condition", cond, "Control"))) %>% mutate(celltype = gp, timepoint = cond) %>% arrange(pvalue)
    res$gene <- rownames(res)
    diff_list[[cond]] <- res
  },mc.cores=ncores)
  stat_list[[i]] <- do.call('rbind',diff_list)
}
stats <- do.call('rbind',stat_list)
write.table(stats,paste0(output,'.pseudobulk_diff.deseq2.txt'),sep='\t',quote=F,row.names=F)

## black list multimappers - we have found certain genes to get very high counts in all samples (>100x higher than other genes).
## We have found these genes to contain repeat elements within their gene bodies and the vast majority of those counts to come from those elements
## As such, we routinely exclude them from analysis

counts_per_gene <- data.frame(gene = rownames(counts), raw_counts = Matrix::rowSums(counts)) %>% 
  arrange(desc(raw_counts)) %>% mutate(rank = 1:nrow(counts))
hist(log10(counts_per_gene$raw_counts))
sum(counts_per_gene$raw_counts > 150000)
out_genes <- counts_per_gene %>% filter(raw_counts > 150000) %>% pull(gene)
dev.off()

## find genes activated during colitis in stem cells
genes1 <- stats %>% filter(!(gene %in% out_genes)) %>% filter(celltype == 'stem_prog') %>% 
  filter(timepoint == 'Colitis_acute') %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>% pull(gene) %>% unique()
genes2 <- stats %>% filter(!(gene %in% out_genes)) %>% filter(celltype == 'stem_prog') %>% 
  filter(timepoint == 'Colitis_chronic') %>%
  filter(padj < 0.05 & log2FoldChange > 0) %>% pull(gene) %>% unique()
length(genes1)
length(genes2)
mem <- stats %>% filter(celltype == 'stem_prog') %>% filter(gene %in% genes1) %>% filter(timepoint == 'Colitis_recovered') %>% filter(padj < 0.05 & log2FoldChange > 0)

# plot
library(ggplot2)
library(cowplot)
dsum2 <- stats %>% filter(gene %in% c(genes1,genes2)) %>% group_by(timepoint,celltype) %>%
  summarise(nUp = sum(padj < 0.05 & log2FoldChange > 0,na.rm = T),
            nDown = sum(padj < 0.05 & log2FoldChange < 0,na.rm = T)
  )
dsum2$timepoint <- factor(dsum2$timepoint,levels=c('Colitis_acute','Colitis_chronic','Colitis_recovered'))
dsum2$celltype <- factor(dsum2$celltype,levels=c('stem_prog','AE_inter','AE_diff'))

gup2 <- ggplot(dsum2, aes(x=timepoint,y=nUp,fill=celltype)) + geom_bar(stat='identity',position='dodge') + 
  theme_bw() + xlab('') + ylab('Number of up-regulated genes') + ggtitle('Colitis activated genes')
gdown2 <- ggplot(dsum2, aes(x=timepoint,y=nDown,fill=celltype)) + geom_bar(stat='identity',position='dodge') + 
  theme_bw() + xlab('') + ylab('Number of down-regulated genes') + ggtitle('Colitis downregulated genes')

p3 <- plot_grid(gup2,gdown2,nrow=1)
print(p3)





