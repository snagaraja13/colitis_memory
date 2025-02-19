## This example is for testing differential motif accessibility for in vivo memory
## The analogous procedure was used for adenomas, human organoids etc, just swapping colitis stages for the corresponding conditions in each dataset

output <- ''
meta_file <- '' # metdata with columns: "barcode","celltype", "replicate", "sample" (for colitis stage). See 01_diff_genes for more details
motif_file <- '' # chromVAR motif deviation object, see 1_SHAREseq_processing
bag_file <- '' # bag file tsv with structure: leader  motif1,motif2,motif3

groups_to_test <- c('stem_prog','AE_inter','AE_diff','immune_other','secretory')

ncores <- 4

library(Matrix)
library(pbmcapply)
library(pheatmap)
library(dplyr)
library(BuenColors)
library(stringi)
library(BuenRTools)
library(chromVAR)

# read in data
meta <- read.table(meta_file,sep='\t',header=T)

bag_data <- read.table(bag_file,header=F,sep='\t')
dev <- readRDS(motif_file)
bagged <- dev[bag_data[,1],]
bagZ <- assays(bagged)$z
rm(dev,bagged)
gc()

if (identical(groups$barcode,colnames(bagZ))){
  cat('Score matrix cells match group matrix\n')
} else {
  cat('Subsetting to scored cells\n')
  groups <- groups[match(colnames(bagZ),groups$barcode),]
}

# find average motif score per replicate
meta$rep_group <- paste0(meta$replicate,'+',meta$celltype)
uniq_rep <- unique(meta$rep_group)
r_list <- pbmclapply(X=uniq_rep,FUN=function(x){
  rep_cells <- meta %>% filter(rep_group == x) %>% pull(barcode)
  bagZ_sub <- bagZ[,rep_cells]
  return(Matrix::rowMeans(bagZ_sub,na.rm=T))
},mc.cores=ncores)
brd <- data.frame(do.call('rbind',r_list))
rownames(brd) <- uniq_rep
write.table(brd,paste0(output,'.motifs_bagged.rep_average.txt'),sep='\t',quote=F)

# now gather data across samples
uniq_samples <- unique(meta$sample)
uniq_cellTypes <- unique(meta$celltype)
rep_data <- stri_split_fixed(rownames(brd),'+',simplify=T)

samp_list <- vector("list",length=length(uniq_samples))
for (i in 1:length(uniq_samples)){
  samp <- uniq_samples[i]
  all_reps <- unique(meta$replicate[meta$sample == samp])
  sum_data <- pbmclapply(X=uniq_cellTypes,FUN=function(x){
    brd_sub <- brd[(rep_data[,1] %in% all_reps)&(rep_data[,2] == x),]
    means <- Matrix::colMeans(brd_sub)
    SDs <- apply(brd_sub,2,sd)
    SEs <- SDs/sqrt(nrow(brd_sub))
    return(data.frame(sample = samp,
                      cell_type = x,
                      motif = names(means),
                      mean = means, 
                      sd = SDs,
                      se = SEs )
    )
  },mc.cores = ncores)
  samp_list[[i]] <- data.frame(do.call('rbind',sum_data))
}
bd <- do.call('rbind',samp_list)
saveRDS(bd,paste0(output,'.motifs_bagged.rep_average.sample_stats.rds'))

#filter bags to test
motif_var <- data.frame(motif = colnames(brd),
                        sd = apply(brd[grepl('stem_prog',rownames(brd)),],2,sd,na.rm=T)
) %>% arrange(desc(sd))
motif_var$rank <- 1:nrow(motif_var)
motifs_to_keep <- motif_var$motif[1:50]

# p-values
pair_list <- list(c('Colitis_acute','Control'),
                  c('Colitis_chronic','Control'),
                  c('Colitis_recovered','Control')
)
test_list <- vector("list",length=length(pair_list))

for (i in 1:length(pair_list)){
  pair <- pair_list[[i]]
  cat(pair,'\n')
  test_reps <- meta %>% filter(sample == pair[1]) %>% pull(replicate) %>% unique()
  control_reps <- meta %>% filter(sample == pair[2]) %>% pull(replicate) %>% unique()
  toMatch <- c(test_reps,control_reps)
  rd_sub <- brd[grep(paste(toMatch,collapse='|'),rownames(brd)),motifs_to_keep]
  
  rd_sub1 <- rd_sub[grep(paste(groups_to_test,collapse="|"),rownames(rd_sub)),]
  stat_list <- vector("list",length=length(groups_to_test))
  for(j in 1:length(groups_to_test)){
    G <- groups_to_test[j]
    sub_G <- rd_sub1[grep(G,rownames(rd_sub1)),]
    t.list <- lapply(X=1:ncol(sub_G),FUN=function(x){
      delta <- mean(sub_G[grep(paste(test_reps,collapse='|'),rownames(sub_G)),x]) - mean(sub_G[grep(paste(control_reps,collapse='|'),rownames(sub_G)),x])
      test <- t.test(x = sub_G[grep(paste(test_reps,collapse='|'),rownames(sub_G)),x],
                     y = sub_G[grep(paste(control_reps,collapse='|'),rownames(sub_G)),x])
      return(c(pair[1],G,colnames(sub_G[x]),delta,test$p.value))
    })
    gp_stats <- data.frame(do.call('rbind',t.list))
    gp_stats$padj <- p.adjust(gp_stats[,5],method='BH')
    stat_list[[j]] <- gp_stats
    
  }
  test_list[[i]] <- data.frame(do.call('rbind',stat_list))
}
stat_data <- do.call('rbind',test_list)
colnames(stat_data) <- c('sample','cell_type','bag_motif','delta','pval','padj')
write.table(stat_data,paste0(output,'.motifs_bagged.by_rep.group_stats.by_timepoint.txt'),row.names=F,sep='\t',quote=F)

## plot
sample_order <- rev(c('Colitis_acute','Colitis_chronic','Colitis_recovered'))
ds <- stat_data %>% filter(cell_type == 'stem_prog')
motifs_to_keep2 <- ds %>% group_by(bag_motif) %>% summarise(minP = min(padj,na.rm=T)) %>% filter(minP < 0.05) %>% pull(bag_motif) # motif significant in any comparison
dplot <- ds %>% filter(bag_motif %in% motifs_to_keep2)

dplot$sample <- factor(dplot$sample,levels = sample_order)
dplot$pval <- as.numeric(dplot$pval)
dplot$delta <- as.numeric(dplot$delta)
dplot$bag_motif <- factor(dplot$bag_motif,levels=(dplot %>% filter(sample == 'Colitis_recovered') %>% arrange(delta) %>% pull(bag_motif)))
max_delta <- max(abs(dplot$delta),na.rm=T)

g2 <- ggplot(dplot,aes(y=sample,x=bag_motif,size=-log10(padj),color=delta)) + geom_point() + 
  theme_bw() + xlab('') + ylab('Motif Family') + ggtitle('Motif Differences - Stem/Progenitor') + 
  theme(axis.text.x = element_text(angle = 45)) +
  scale_color_gradient2(low='navy',mid='white',high='firebrick') + 
  geom_point(data=dplot,aes(y=sample,x=bag_motif,size=-log10(padj)),shape=1,color='gray25') + 
  geom_point(data = (dplot %>% filter(padj < 0.05)), aes(y=sample,x=bag_motif), shape = 8, color = 'black', size = 0.75)

pdf(paste0(output,'.motifs_bagged.by_rep.group_stats.by_timepoint.top_bags.pdf'),width=12,height=6)
print(g2)
dev.off()
