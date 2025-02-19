output <- 'DSS-prim-merge10-ATAC.rescore.cisBP_only'
dev_file <- '../DSS-prim-merge10-ATAC.modisco.dev_motif.unbagged.rds'
meta_file <- '../DSS-prim-merge10-ATAC.metadata.txt'

motifs_to_plot <- c('Fos','Etv2','Bcl11a','Spib','Spi1','Rorc')

library(chromVAR)
library(ggplot2)
library(BuenRTools)
library(cowplot)
library(dplyr)
library(ggpubr)

###
dev <- readRDS(dev_file)
z <- assays(dev)$z
meta <- read.table(meta_file,header=T) %>% filter(group == 'stem_prog')
shared_cells <- intersect(colnames(z),meta$barcode)
meta <- meta[match(shared_cells,meta$barcode),]
z <- z[,shared_cells]
identical(meta$barcode,colnames(z))

r_list <- lapply(X=motifs_to_plot,FUN=function(motif){
  d <- cbind(meta,score = z[motif,]) %>% filter(sample %in% c('Control','Colitis_recovered'))
  dr <- d %>% group_by(replicate) %>% summarise(Sample = unique(sample), n_cells = n(),
                                                TF = motif,
                                                n_gt1 = sum(score > 1),
                                                n_gt1.5 = sum(score > 1.5),
                                                n_gt2 = sum(score > 2),
                                                n_gt2.5 = sum(score > 2.5)
  ) %>%
    mutate(frac_gt1 = n_gt1/n_cells,
           frac_gt1.5 = n_gt1.5/n_cells,
           frac_gt2 = n_gt2/n_cells,
           frac_gt2.5 = n_gt2.5/n_cells)
  
  return(dr)
})

rd <- do.call('rbind',r_list)

# get sig stats
stat_list <- lapply(X=motifs_to_plot,FUN=function(motif){
  rd_sub <- rd %>% filter(TF == motif)
  ttest1 <- t.test(x = (rd_sub %>% filter(Sample == 'Colitis_recovered') %>% pull(frac_gt1)),
                   y = (rd_sub %>% filter(Sample == 'Control') %>% pull(frac_gt1)),
                   var.equal = T
                   )
  ttest1.5 <- t.test(x = (rd_sub %>% filter(Sample == 'Colitis_recovered') %>% pull(frac_gt1.5)),
                   y = (rd_sub %>% filter(Sample == 'Control') %>% pull(frac_gt1.5)),
                   var.equal = T
  )
  
  ttest2 <- t.test(x = (rd_sub %>% filter(Sample == 'Colitis_recovered') %>% pull(frac_gt2)),
                     y = (rd_sub %>% filter(Sample == 'Control') %>% pull(frac_gt2)),
                     var.equal = T
  )
  
  return(data.frame(TF = motif, pval1 = ttest1$p.value, pval1.5 = ttest1.5$p.value, pval2 = ttest2$p.value))
})
stats <- do.call('rbind',stat_list)

rd$Sample <- factor(rd$Sample,levels = c('Control','Colitis_recovered'))

g1 <- ggbarplot(rd,x='TF',y='frac_gt1',fill='Sample', add = c('mean_se'),
          position = position_dodge()
          ) + 
  geom_point(data = rd, aes(x=TF,y=frac_gt1,fill = Sample), shape = 1, size=2,stat = 'identity', position = position_jitterdodge(jitter.width=0.2,jitter.height = 0)) + 
  geom_point(data = stats %>% filter(pval1 < 0.05), aes(x=TF,y=pval1),shape=8, y = max(rd$frac_gt1))

g1.5 <- ggbarplot(rd,x='TF',y='frac_gt1.5',fill='Sample', add = c('mean_se'),
                position = position_dodge()
) + 
  geom_point(data = rd, aes(x=TF,y=frac_gt1.5,fill = Sample), shape = 1, size=2,stat = 'identity', position = position_jitterdodge(jitter.width=0.2,jitter.height = 0)) + 
  geom_point(data = stats %>% filter(pval1.5 < 0.05), aes(x=TF,y=pval1.5),shape=8, y = max(rd$frac_gt1.5))


g2 <- ggbarplot(rd,x='TF',y='frac_gt2',fill='Sample', add = c('mean_se'),
                position = position_dodge()
) + 
  geom_point(data = rd, aes(x=TF,y=frac_gt2,fill = Sample), shape = 1, size=2,stat = 'identity', position = position_jitterdodge(jitter.width=0.2,jitter.height = 0)) + 
  geom_point(data = stats %>% filter(pval2 < 0.05), aes(x=TF,y=pval2),shape=8, y = max(rd$frac_gt2))


pdf(paste0(output,'.top_memory_TFs.high_barplot.pdf'),width=8,height=5)
print(g1)
print(g1.5)
print(g2)
dev.off()
