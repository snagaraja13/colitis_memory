output <- ''
dev_file <- '' # chromVAR motif deviation object, see 1_SHAREseq_processing
meta_file <- '' # metadata with columns: 

motifs_to_plot <- c('Fos','Etv2','Spib','Spi1','Rorc','Ctcf')

library(chromVAR)
library(ggplot2)
library(BuenRTools)
library(cowplot)
library(dplyr)
library(ggpubr)

### read in data
dev <- readRDS(dev_file)
z <- assays(dev)$z
meta <- read.table(meta_file,header=T) %>% filter(celltype == 'stem_prog')
shared_cells <- intersect(colnames(z),meta$barcode)
meta <- meta[match(shared_cells,meta$barcode),]
z <- z[,shared_cells]
identical(meta$barcode,colnames(z))

### find high cells per replicate
# Motif scores should not be considered comparable between datasets that were not scored together.
# The deviations are dependent on background peaks, which are generated to control for bias within each dataset, as well as average accessibility within the scored cells
# As such, cutoffs should be varied for each dataset

r_list <- lapply(X=motifs_to_plot,FUN=function(motif){
  d <- cbind(meta,score = z[motif,]) %>% filter(sample %in% c('Control','Colitis_recovered'))
  dr <- d %>% group_by(replicate) %>% summarise(Sample = unique(sample), n_cells = n(),
                                                TF = motif,
                                                n_gt1.5 = sum(score > 1.5)
                                                ) %>% mutate(frac_gt1.5 = n_gt1.5/n_cells)
  
  return(dr)
})

rd <- do.call('rbind',r_list)

# get sig stats
stat_list <- lapply(X=motifs_to_plot,FUN=function(motif){
  rd_sub <- rd %>% filter(TF == motif)
  ttest <- t.test(x = (rd_sub %>% filter(Sample == 'Colitis_recovered') %>% pull(frac_gt1.5)),
                  y = (rd_sub %>% filter(Sample == 'Control') %>% pull(frac_gt1.5)),
                  var.equal = T
                  )
  return(data.frame(TF = motif, pval1.5 = ttest1.5$p.value))
})
stats <- do.call('rbind',stat_list)

rd$Sample <- factor(rd$Sample,levels = c('Control','Colitis_recovered'))

# plot
g1 <- ggbarplot(rd,x='TF',y='frac_gt1.5',fill='Sample', add = c('mean_se'),
                position = position_dodge()
) + 
  geom_point(data = rd, aes(x=TF,y=frac_gt1.5,fill = Sample), shape = 1, size=2,stat = 'identity', position = position_jitterdodge(jitter.width=0.2,jitter.height = 0)) + 
  geom_point(data = stats %>% filter(pval1.5 < 0.05), aes(x=TF,y=pval1.5),shape=8, y = max(rd$frac_gt1.5))


pdf(paste0(output,'.top_memory_TFs.high_barplot.pdf'),width=8,height=5)
print(g1)
dev.off()

### Plot single cell data
library(ggbeeswarm)
library(dplyr)
nDown <- 500 # each condition was downsampled to 500 cells to account for there being many more control mice than colitis

pdf(paste0(output,'.diff_motifs.memory_beeswarm.pdf'),width=8,height=8)
for (motif in motifs_to_plot){
  cat(motif,'\n')
  d <- cbind(meta,score = z[motif,])
  dg <- d[d$sample %in% c('Control','Colitis_recovered'),]
  
  control_ind <- which(dg$sample == 'Control')
  recovered_ind <- which(dg$sample == 'Colitis_recovered')
  if (length(control_ind) < nDown | length(recovered_ind) < nDown){
    nDown <- min(length(control_ind),length(recovered_ind))
  }
  
  control_down <- sample(control_ind,size = nDown)
  recovered_down <- sample(recovered_ind,size = nDown)
  dg <- dg[c(control_down,recovered_down),]
  dg$sample <- factor(dg$sample,levels=c('Control','Colitis_recovered'))
  
  g2 <- ggplot(dg,aes(x=sample,y=score,color=sample)) + 
    geom_beeswarm(size = 0.1) + 
    geom_crossbar(data = (dg %>% group_by(sample) %>% summarise(avg_score = mean(score, na.rm=T))),aes(x=sample,y = avg_score,ymin=avg_score, ymax=avg_score), width = 0.2, color = 'black') + 
    theme_bw() + xlab('') + ylab('Motif accessibility score') + ggtitle(motif) + guides(color = 'none')
  print(g2)
  #cont <- readline(prompt = 'Enter')
}
dev.off()
