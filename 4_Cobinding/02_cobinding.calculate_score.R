## Calculating co-binding scores from seq2PRINT output
## Example here is for in vivo memory

output <- ''
all_tfbs_file <- '' # filtered and quantile normalized TF footprint score file from seq2PRINT
all_peak_file <- '../ref/mouse_colitis_tissue.clean_peaks.bed' # peaks associated with the seq2PRINT output
bag_file <- ''  # motif family file tsv with structure: leader  motif1,motif2,motif3
## pairwise differences using corresponding controls from each batch
## for other datasets, change differentials accordingly below
pair_list <- list('acute' = c('Colitis_acute','Control_3'),
                  'chronic' = c('Colitis_chronic','Control_1'),
                  'recovered' = c('Colitis_recovered','Control_2')
                  )

nMotif <- 20
tfbs_thresh <- 0.2
ncores <- 6

library(Matrix)
library(dplyr)
library(SummarizedExperiment)
library(BuenRTools)
library(stringi)
library(GenomicRanges)
library(pbmcapply)

###
### 1) Define TF families
###

bag_data <- read.table(bag_file,sep='\t')
getBagMembers <- function(TF){return(c(stri_split_fixed(bag_data[bag_data[,1] == TF,2],',',simplify = T)))}
grouped_list <- list('AP1' = getBagMembers('Fos'),
                     'NR1/2' = getBagMembers('Hnf4g'),
                     'ETS' = c('Spi1','Spib','Spic',getBagMembers('Etv2')),
                     'BCL11A' = 'Bcl11a',
                     'CTCF' = getBagMembers('Ctcf'),
                     'NFI' = getBagMembers('Nfix'),
                     'HOX/CDX' = getBagMembers('Hoxd10'),
                     'RUNX' = c('Runx1','Runx2','Runx3'),
                     'FOX' = getBagMembers('Foxa1'),
                     'NFkB/REL' = getBagMembers('Nfkb1'),
                     'SNAI/MESP' = c('Mesp1','Mesp2','Snai2','Snai3','Tcf3','Tcf4'),
                     'ESRR' = c('Esrra','Esrrb','Esrrg'),
                     'Retinoid' = c('Rara','Rarb','Rxra','Rxrb','Rxrg'),
                     'KLF' = getBagMembers('Klf1'),
                     'SOX' = getBagMembers('Sox9')
                     )

rm(dev)
gc()

###
### 2) Find all memory TFBS motifs
###
tfbs_all <- read.table(all_tfbs_file,sep='\t',header=T)
# adjust this step accordingly for different datasets, eg (T5224 - DMSO) or (IBD - healthy)
delta_list <- lapply(X=pair_list,FUN=function(pair){
  delta <- tfbs_all[,grepl(pair[1],colnames(tfbs_all))] - tfbs_all[,grepl(pair[2],colnames(tfbs_all))]
  return(delta)
})
dmat <- do.call('cbind',delta_list)
colnames(dmat) <- paste0(names(pair_list),'_delta')
head(dmat)
filt_ind <- (apply(dmat,1,max) > tfbs_thresh) # for T5224 experiments, we are looking for loss in signal so this would be < (-tfbs_thresh)
tfbs_filt <- cbind(tfbs_all,dmat)[filt_ind,]
dim(tfbs_filt)
rm(tfbs_all)
gc()

# all peaks with any memory TFBS
all_peaks <- makeGRangesFromDataFrame(read.table(all_peak_file,col.names=c('chr','start','end')))
filt_ov <- findOverlaps(all_peaks,makeGRangesFromDataFrame(tfbs_filt[,1:3]))
filt_peaks <- all_peaks[unique(filt_ov@from)]

###
### 3) Find overlaps for all bag pairs and enrichment
###
all_bag_tfs <- unlist(grouped_list)
tfbs_sub <- tfbs_filt %>% filter(TF %in% all_bag_tfs)
dim(tfbs_sub)

# rename with bag
tf2_list <- lapply(X=tfbs_sub$TF,FUN=function(tf){
  for (i in 1:length(grouped_list)){
    if (tf %in% grouped_list[[i]]){
      return(names(grouped_list)[i])
    }
  }
})
tfbs_sub$bag <- unlist(tf2_list)

# iterate across all pairs and find enrichment
uniq_bags <- unique(tfbs_sub$bag)
n_total <- length(filt_peaks)

ed_list <- list()
for (i in 1:length(pair_list)){
  tp <- names(pair_list)[i]
  cat('--',tp,'--\n')
  e_list <- pbmclapply(X=uniq_bags,FUN=function(tf1){
    tp_delta <- tfbs_sub[,paste0(tp,'_delta')]
    tf1_ov <- findOverlaps(all_peaks,makeGRangesFromDataFrame(tfbs_sub %>% filter(tp_delta > tfbs_thresh) %>% filter(bag == tf1)))
    peaks_tf1 <- all_peaks[unique(tf1_ov@from)]
    n_tf1 <- length(peaks_tf1)
    tf2_list <- lapply(X=uniq_bags,FUN=function(tf2){
      tf2_ov <- findOverlaps(all_peaks,makeGRangesFromDataFrame(tfbs_sub %>% filter(tp_delta > tfbs_thresh) %>% filter(bag == tf2)))
      peaks_tf2 <- all_peaks[unique(tf2_ov@from)]
      
      shared_ov <- findOverlaps(peaks_tf1,peaks_tf2)
      n_both <- length(unique(shared_ov@from))
      n_only1 <- n_tf1 - n_both
      n_only2 <- length(peaks_tf2) - n_both
      n_neither <- n_total - (n_both + n_only1 + n_only2)
      
      d_cont <- data.frame(hasTF1 = c(n_both, n_only1),
                           noTF1 = c(n_only2,n_neither)
      )
      
      
      f.test <- fisher.test(d_cont, alternative = 'greater')
      #f.test <- fisher.test(d_cont)
      return(data.frame(TF1 = tf1, TF2 = tf2, condition = tp,OR = f.test$estimate, pval = f.test$p.value))
    })
    enrich_data <- do.call('rbind',tf2_list)
    return(enrich_data)
  },mc.cores=ncores)
  ed_list[[i]] <- do.call('rbind',e_list)
}
ed <- do.call('rbind',ed_list)
write.table(ed,paste0(output,'.seq2PRINT_TFBS.',tfbs_thresh,'.grouped.enrichment_stats.csv'),sep=',',quote=F,row.names=F)

###
### 4) Convert enrichment to co-binding score
###
# reformat as matrix
mat_list <- list()
for (i in 1:length(pair_list)){
  tp <- names(pair_list)[[i]]
  tmat_list <- lapply(X=uniq_bags,FUN=function(motif){
    e_sub <- ed %>% filter(condition == tp) %>% filter(TF1 == motif) 
    e_sub <- e_sub[match(uniq_bags,e_sub$TF2),]
    return(log2(e_sub$OR))
  })
  te_mat <- do.call('rbind',tmat_list)
  colnames(te_mat) <- uniq_bags
  rownames(te_mat) <- uniq_bags
  for (i in 1:nrow(te_mat)){
    te_mat[i,i] <- NA
  }
  mat_list[[tp]] <- te_mat
}


# scale -1 to 1
all_vals <- c(do.call('cbind',mat_list))
# scale positive values 0 to 1 and negative 0 to -1
probs_range <- c(0.1,0.9)
pos_perc <- quantile(all_vals[all_vals > 0],probs=probs_range,na.rm=T)
neg_perc <- quantile(all_vals[all_vals < 0],probs=probs_range,na.rm=T)
apos <- 1/(pos_perc[2] - pos_perc[1])
bpos <- -1*(pos_perc[1]*apos)
aneg <- 1/(neg_perc[2] - neg_perc[1])
bneg <- -1*(neg_perc[2]*aneg)
scaleEnrichment <- function(vals){
  scale_list <- lapply(X=vals,FUN=function(x){
    if (is.na(x)){
      return(NA)
    } else if (x > 0){
      return(apos*x + bpos)
    } else if (x < 0) {
      return(aneg*x + bneg)
    } else {
      return(0)
    }
  })
  scale_vals <- unlist(scale_list)
  names(scale_vals) <- names(vals)
  return(scale_vals)
}

cmat_list <- list()
for (i in 1:length(mat_list)){
  e_mat <- mat_list[[i]]
  c_list <- lapply(X=1:ncol(e_mat),FUN=function(j){
    e_vals <- scaleEnrichment(e_mat[,j])
    e_vals[e_vals > 1] <- 1
    e_vals[e_vals < -1] <- -1
    return(e_vals)
  })
  c_mat <- do.call('cbind',c_list)
  colnames(c_mat) <- colnames(e_mat)
  cmat_list[[names(mat_list)[i]]] <- c_mat
}
saveRDS(cmat_list,paste0(output,'.seq2PRINT_TFBS.',tfbs_thresh,'.grouped.cobinding_matrix.rds'))

# plot
library(ggplot2)
library(cowplot)

pdf(paste0(output,'.seq2PRINT_TFBS.',tfbs_thresh,'.grouped_cooperativity.barplots.pdf'),width=8,height=5)
for (TF in uniq_bags){
  g_list <- lapply(X=1:length(cmat_list),FUN=function(i){
    c_mat <- cmat_list[[i]]
    dc <- data.frame(tf = rownames(c_mat), coop = c_mat[,TF]) %>% arrange(desc(coop)) %>% mutate(rank = 1:nrow(c_mat))
    dc$tf <- factor(dc$tf,levels=dc$tf)
    g <- ggplot(dc %>% filter(!is.na(coop)),aes(x=tf,y=coop, fill = coop)) + geom_bar(stat='identity',color='black') +
      theme_bw() + xlab('Ranked motifs') + ylab('Co-binding') + ggtitle(paste0(TF,' - ',names(cmat_list)[i])) +
      geom_hline(yintercept = 0, linetype = 'dashed') + 
      scale_fill_gradient2(low = 'navy', mid = 'beige', high = 'firebrick', midpoint = 0, limits = c(-1,1)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
    return(g)
  })
  p <- do.call('plot_grid',c(g_list,nrow=1))
  print(p)
}
dev.off()