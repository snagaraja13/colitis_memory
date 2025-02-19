## For testing whether a given feature displays "clonality" - whether variance in the feature amongst clonal cells is less than chance
## The code here is for motif accessibility but may be generalized to any single cell feature

output <- ''
type <- '' # "motif" or some other type of feature
clean_clone_call_file <- '../ref/mouse_organoids.cleaned_clone_calls.txt' # see 03_SHARE_TRACE.clean_clones
score_file <- '' # chromVAR motif deviation object
bag_file <- '' # motif family file tsv with structure: leader  motif1,motif2,motif3
meta_file <- '' # scATAC-seq metdata file  

scores_to_plot <- c("Fos","Fosl1","Hnf4g","Ctcf",
                    "Mafk","Etv2","Hoxd10","Hhex","Sox9",
                    "Tead4","Tead2","Tcf7","Tcf7l2")

cores <- 8
nPerm <- 1000

library(BuenRTools)
library(stringi)
library(pbmcapply)
library(reshape2)
library(dplyr)
library(cowplot)
if (type %in% c('motif')){
  library(chromVAR)
}
library(ggrepel)

findCloneStats <- function(d,features,ncores=4,seed=123){
  set.seed(seed)
  feat_list <- pbmclapply(X=features,FUN=function(x){
    d_sub <- d %>% filter(feature == x)
    
    # find clonal distribution
    clone_var <- data.frame(d_sub %>% group_by(clone) %>% summarise(sd(score,na.rm = T)))
    avg_clone_var <- median(clone_var[,2]^2)
    
    if (sum(is.na(clone_var[,2])) == nrow(clone_var)){
      return( data.frame(feature = x,
                         clone_median_var = avg_clone_var,
                         shuf_median_var = NA,
                         shuf_sd_var = NA,
                         pval = NA
                         )
              )
    } else {
      # find permuted distribution
      perm_list <- lapply(X=1:nPerm,FUN=function(i){
        d_shuf <- d_sub
        d_shuf$clone <- sample(d_shuf$clone)
        shuf_data <- data.frame(d_shuf %>% group_by(clone) %>% summarise(sd(score,na.rm = T)))
        return(median(shuf_data[,2]))
      })
      shuf_data <- unlist(perm_list)
      shuf_mean_var <- mean(shuf_data^2)
      shuf_sd_var <- sd(shuf_data^2)
      clone_z <- (avg_clone_var - shuf_mean_var)/(shuf_sd_var)
      # find stats
      feat_stats <- data.frame(feature = x,
                               clone_median_var = avg_clone_var,
                               shuf_median_var = shuf_mean_var,
                               shuf_sd_var = shuf_sd_var,
                               pval = 2*pnorm(-abs(clone_z))
                               )
      return(feat_stats)
    }
  },mc.cores = ncores)
  return(data.frame(do.call('rbind',feat_list)))
}


# read in data
d <- readRDS(score_file)
meta <- read.table(meta_file,sep='\t',header=T)
if(!identical(meta$barcode,colnames(d))){
  cat('Subsetting metadata to barcodes in data\n')
  meta <- meta %>% filter(barcode %in% colnames(d))
}
clone_calls <- read.table(clean_clone_call_file,sep='\t',header=T)

bag_data <- read.table(bag_file,sep='\t')
z <- assays(d[bag_data[,1],])$z

# score by clone for all motifs
if (sum(clone_calls$atac_cell %in% colnames(d)) != length(clone_calls$atac_cell %in% colnames(d))){
  cat('Subsetting clone data to barcodes in data\n')
  clone_calls <- clone_calls %>% filter(atac_cell %in% meta$barcode)
}
z_bc <- z[,clone_calls$atac_cell]
clone_counts <- table(unlist(lapply(X=clone_calls$clone,FUN=function(x){return(stri_split_fixed(x,',',simplify=T))})))
keep_clones <- names(clone_counts)[clone_counts > 1]
z_by_clone <- pbmclapply(X=keep_clones,FUN=function(x){
  clone_ind <- (Matrix::rowSums(stri_split_fixed(clone_calls$clone,',',simplify=T) == x) > 0)
  z_clone <- melt(z_bc[,clone_calls$atac_cell[clone_ind]])
  colnames(z_clone) <- c('feature','atac_cell','score')
  z_clone$clone <- x
  return(z_clone)
},mc.cores=cores)
z_by_clone <- data.frame(do.call('rbind',z_by_clone))
z_by_clone$sample <- meta$sample[match(z_by_clone$atac_cell,meta$barcode)]
write.table(z_by_clone,paste0(output,'.scores_by_clone.txt'),sep='\t',row.names=F,quote=F)

# find variance for each TF by clone and permute
uniq_feat <- as.character(unique(z_by_clone$feature))
all_stats <- findCloneStats(d = z_by_clone,
                            features = uniq_feat,
                            ncores = cores,
                            seed = 123
                            )
all_stats$fdr <- p.adjust(all_stats$pval,method='BH')
write.table(all_stats,paste0(output,'.clone_var_stats.txt'),sep='\t',quote=F,row.names=F)

# plot summary stats
g1 <- ggplot(all_stats,aes(x=(clone_median_var - shuf_median_var),y=-log10(fdr))) + geom_point(shape=1) + 
  theme_bw() + xlab('Clonal SD - Random SD') + ylab('-log10(FDR)') + 
  geom_hline(yintercept=-log10(0.05),linetype='dashed') + 
  geom_point(data = all_stats %>% filter(feature %in% scores_to_plot),aes(x=(clone_median_var - shuf_median_var),y=-log10(fdr)),shape=16,col='red') +
  geom_text_repel(data=(all_stats %>% filter(fdr < 0.05 | feature %in% scores_to_plot)),aes(x=(clone_median_var - shuf_median_var),y=-log10(fdr),label=feature),color='red',max.overlaps=50)

pdf(paste0(output,'.clone_stats.by_sample.pdf'),width=5,height=5)
print(g1)
dev.off()