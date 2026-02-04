output <- ''
# score file with columns for feature/motif, atac_cell_barcode, feature_score (e.g chromvar), clone_call
score_by_clone_file <- '' 

# metdata file with columns barcode (ATAC barcode), sample (exposure to DSS colitis or not)
meta_file <- '' 

# plot favorite ones
top_features <- c('Fos','Hnf4g','Hhex','Hoxd10','Tcf7l2','Sox9','Runx1')

cores <- 4

library(dplyr)
library(Matrix)
library(pbmcapply)
library(variancePartition)
library(stringi)

## make metadata matrix
scores <- read.table(score_by_clone_file,header=T,sep='\t')
all_meta <- read.table(meta_file,header=T,sep='\t')
shared_cells <- intersect(scores$atac_cell,all_meta$barcode)
length(shared_cells)
d <- data.frame(is_colitis = all_meta$sample[match(shared_cells,all_meta$barcode)],
                clone = scores$clone[match(shared_cells,scores$atac_cell)],
                row.names = shared_cells
                )
head(d)

## make score matrix, columns = atac cell barcodes, rows = features (e.g. motifs)
uniq_feat <- unique(scores$feature)
score_list <- lapply(X=uniq_feat,FUN=function(feat){
  score_sub <- scores %>% filter(feature == feat) 
  return(score_sub$score[match(shared_cells,score_sub$atac_cell)])
})
score_data <- do.call('rbind',score_list)
colnames(score_data) <- shared_cells
rownames(score_data) <- uniq_feat
score_data[is.na(score_data)] <- 0
head(score_data[,1:5])

## run linear mixed model
form <- ~ (1|is_colitis) + (1|clone)
varPart <- fitExtractVarPartModel(score_data, form, d) 
var.df <- as.data.frame(varPart)
# reformat
stat_list <- lapply(X=1:ncol(var.df),FUN=function(i){
  return(data.frame(feature = rownames(var.df),
                    stat = 'variance',
                    variable = colnames(var.df)[i],
                    value = var.df[,i]
                    )
         )
})
stats <- do.call('rbind',stat_list)
write.table(stats,paste0(output,'.variance.stats.txt'),sep='\t',quote=F,row.names=F)

## plots
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggrepel)

# select TFs
ds <- stats %>% filter(feature %in% top_features) %>% filter(stat == 'variance')
feat_order <- ds %>% filter(variable == 'clone') %>% arrange(desc(value)) %>% pull(feature)
ds$feature <- factor(ds$feature,levels=feat_order)

g1 <- ggbarplot(data = (ds %>% filter(variable != 'Residuals')),
          x = 'feature',y='value',fill='variable',add=c('mean_se'),
          position = position_dodge(width = 0.75),
          xlab = '', ylab = 'Fraction of variance explained', title = output
          )

pdf(paste0(output,'.top_feature.variance_stats.pdf'),width=8,height=4)
print(g1)
dev.off()

# all features
d1 <- stats %>% filter(stat == 'variance') %>% filter(variable == 'clone') %>% arrange(desc(value)) %>% mutate(rank_clone = 1:length(value))
feat_to_plot1 <- unique(c(d1$feature[1:20],top_features))
d2 <- stats %>% filter(stat == 'variance') %>% filter(variable == 'is_colitis') %>% arrange(desc(value)) %>% mutate(rank_colitis = 1:length(value))
feat_to_plot2 <- unique(c(d2$feature[1:20],top_features))

g3 <- ggplot(d1, aes(x=rank_clone,y=100*value)) + geom_point() + 
  theme_bw() + xlab('Rank') + ylab('Variance explained by clone (%)') + ggtitle('Clone') + 
  geom_point(data = d1 %>% filter(feature %in% feat_to_plot1),aes(x=rank_clone,y=100*value),color='red') + 
  geom_text_repel(data = d1 %>% filter(feature %in% feat_to_plot1),aes(x=rank_clone,y=100*value,label=feature),max.overlaps = 50,color='red')

g4 <- ggplot(d2, aes(x=rank_colitis,y=100*value)) + geom_point() + 
  theme_bw() + xlab('Rank') + ylab('Variance explained by colitis (%)') + ggtitle('Colitis') +
  geom_point(data = d2 %>% filter(feature %in% feat_to_plot2),aes(x=rank_colitis,y=100*value),color='red') + 
  geom_text_repel(data = d2 %>% filter(feature %in% feat_to_plot2),aes(x=rank_colitis,y=100*value,label=feature),max.overlaps = 50,color='red')

d3 <- data.frame(feature = d1$feature,
            value_clone = d1$value,
            rank_clone = d1$rank_clone,
            value_colitis = d2$value[match(d1$feature,d2$feature)],
            rank_colitis = d2$rank_colitis[match(d1$feature,d2$feature)]
            )
head(d3)
feat_to_plot3 <- unique(feat_to_plot1,feat_to_plot2)

g5 <- ggplot(d3,aes(x=100*value_clone,y=100*value_colitis)) + geom_point() + 
  geom_abline(slope=1,intercept=0,linetype = 'dashed') +
  geom_point(data = d3 %>% filter(feature %in% feat_to_plot3),aes(x=100*value_clone,y=100*value_colitis),color='red') + 
  theme_bw() + xlab('Variance explained by clone (%)') + ylab('Variance explained by colitis (%)') + ggtitle('Clone vs Colitis') +
  geom_text_repel(data = d3 %>% filter(feature %in% feat_to_plot3),aes(x=100*value_clone,y=100*value_colitis,label=feature),max.overlaps=50)


p2 <- plot_grid(g3,g4,g5,nrow=1)

pdf(paste0(output,'.all_feature.variance_stats.pdf'),width=12,height=4)
print(p2)
dev.off()
