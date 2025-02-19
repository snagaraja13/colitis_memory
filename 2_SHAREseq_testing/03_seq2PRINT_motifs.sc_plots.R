output <- 'DSS-prim-merge10-ATAC.rescore.cisBP_only'
dev_file <- '../DSS-prim-merge10-ATAC.modisco.dev_motif.unbagged.rds'
meta_file <- '/mnt/users/snagaraja/2022-12-20_DSS-share/DSS-prim-merge10/DSS-prim-merge10-ATAC/cell_groups/DSS-prim-merge10-ATAC.cell_groups.txt'

motifs_to_plot <- c('Fos','Etv2','Spib','Spi1','Rorc','Ctcf')

library(chromVAR)
library(ggplot2)
library(BuenRTools)
library(cowplot)

###
dev <- readRDS(dev_file)
z <- assays(dev)$z
meta <- read.table(meta_file,header=T)
identical(meta$barcode,colnames(z))

uniq_groups <- unique(meta$group)
## all timepoints
pdf(paste0(output,'.diff_motifs.all_sample_violin.pdf'),width=8,height=8)
for (motif in motifs_to_plot){
  cat(motif,'\n')
  d <- cbind(meta,score = z[motif,])
  g_list1 <- lapply(X=uniq_groups,FUN=function(gp){
    dg <- d[d$group == gp,]
    dg$sample <- factor(dg$sample,levels=samp_order)
    
    g1 <- ggplot(dg,aes(x=sample,y=score,fill=sample)) + 
      geom_violin(trim=T) + 
      #geom_jitter(width = 0.3, height = 0, size = 0.1) + 
      geom_boxplot(outlier.shape = NA, width = 0.2) + 
      theme_bw() + xlab('') + ylab('Motif Score') + ggtitle(paste0(gp,' - ',motif)) + 
      guides(fill = 'none')
    return(g1)
  })
  p1 <- do.call('plot_grid',g_list1)
  print(p1)
  #cont <- readline(prompt = 'Enter')
}
dev.off()

pdf(paste0(output,'.diff_motifs.memory_violin.pdf'),width=8,height=8)
for (motif in motifs_to_plot){
  cat(motif,'\n')
  d <- cbind(meta,score = z[motif,])
  d <- d[d$sample %in% c('Control','Colitis_recovered'),]
  g_list2 <- lapply(X=uniq_groups,FUN=function(gp){
    dg <- d[d$group == gp,]
    dg$sample <- factor(dg$sample,levels=samp_order)
    g2 <- ggplot(dg,aes(x=sample,y=score,fill=sample)) + 
      geom_violin(trim=T) + 
      #geom_jitter(width = 0.3, height = 0, size = 0.1) + 
      geom_boxplot(outlier.shape = NA, width = 0.2) + 
      theme_bw() + xlab('') + ylab('Motif Score') + ggtitle(paste0(gp,' - ',motif)) + 
      guides(fill = 'none')
    return(g2)
  })
  p2 <- do.call('plot_grid',g_list2)
  print(p2)
  #cont <- readline(prompt = 'Enter')
}
dev.off()

## dot plots
library(ggbeeswarm)
library(dplyr)
nDown <- 500

pdf(paste0(output,'.diff_motifs.memory_beeswarm.pdf'),width=8,height=8)
for (motif in motifs_to_plot){
  cat(motif,'\n')
  d <- cbind(meta,score = z[motif,])
  d <- d[d$sample %in% c('Control','Colitis_recovered'),]
  g_list2 <- lapply(X=uniq_groups,FUN=function(gp){
    dg <- d[d$group == gp,]
    control_ind <- which(dg$sample == 'Control')
    recovered_ind <- which(dg$sample == 'Colitis_recovered')
    if (length(control_ind) < nDown | length(recovered_ind) < nDown){
      nDown <- min(length(control_ind),length(recovered_ind))
    }
    
    control_down <- sample(control_ind,size = nDown)
    recovered_down <- sample(recovered_ind,size = nDown)
    dg <- dg[c(control_down,recovered_down),]
    dg$sample <- factor(dg$sample,levels=samp_order)
    
    g2 <- ggplot(dg,aes(x=sample,y=score,color=sample)) + 
      geom_beeswarm(size = 0.1) + 
      geom_crossbar(data = (dg %>% group_by(sample) %>% summarise(avg_score = mean(score, na.rm=T))),aes(x=sample,y = avg_score,ymin=avg_score, ymax=avg_score), width = 0.2, color = 'black') + 
      theme_bw() + xlab('') + ylab('Motif accessibility score') + ggtitle(paste0(gp,' - ',motif)) + guides(color = 'none')
    
    return(g2)
  })
  p2 <- do.call('plot_grid',g_list2)
  print(p2)
  #cont <- readline(prompt = 'Enter')
}
dev.off()
