output <- ''
frac_file <- '' # file with methylation fraction per peak across all CpGs
meta_file <- '' # metadata file with columns: replicate, condition (DSS colitis vs control), treatment (DMSO vs T5224)
annot_file <- '' # motif annotatin file, peaks x motifs
annot_features <- c('Fos') # testing AP-1 sites here
bg_file <- '' # previously generated background peaks for peak set in annotatin file, from chromVAR

nvar <- 500
cores <- 4

library(dplyr)
library(Matrix)
library(pbmcapply)
library(ggplot2)
library(SummarizedExperiment)
library(cowplot)

meta <- read.table(meta_file,header=T,sep='\t') %>% 
  mutate(sample = paste0(condition,'_',treatment))
d <- readRDS(frac_file)
identical(colnames(d),meta$replicate)

# filter and subset to annotated peaks
peak_stats <- data.frame(region  = 1:nrow(d),
                         peak_avg = apply(d,1,mean,na.rm=T),
                         peak_max = apply(d,1,max,na.rm=T),
                         peak_sd = apply(d,1,sd,na.rm=T)
)
keep_ind <- peak_stats %>% filter(peak_max >= 0.1 & peak_sd > 0.05) %>% pull(region) # at least 1 peak with at least 10% methylation, S.D. at least 0.05
length(keep_ind)

# annotated peaks
annot <- readRDS(annot_file)
if (length(annot_features) > 1){
  feat_peaks <- Matrix::rowSums(assays(annot)$motifMatches[,annot_features]) > 0
} else {
  feat_peaks <- assays(annot)$motifMatches[,annot_features] > 0
}
sum(feat_peaks)
feat_ind <- intersect(which(feat_peaks),keep_ind) # intersect variable peaks with AP1 peaks
length(feat_ind)
pf1 <- d[feat_ind,]

# most variable features
feat_var <- data.frame(mean = Matrix::rowMeans(pf1), sd = apply(pf1,1,sd)) %>% arrange(desc(sd)) %>% mutate(rank = 1:nrow(pf1))
ggplot(feat_var,aes(x=rank,y=sd)) + geom_point()
top_regions <- rownames(feat_var[1:nvar,])
pf2 <- pf1[rownames(pf1) %in% top_regions,]

# reformat for plotting
samp_order <- c('Control_DMSO','Control_T5224','DSS_DMSO','DSS_T5224')
rep_order <- lapply(X=samp_order,FUN=function(samp){return(meta %>% filter(sample == samp) %>% pull(replicate) %>% sample())})
rep_order <- unlist(rep_order)
pf2 <- pf2[,rep_order]
meta <- meta[match(rep_order,meta$replicate),]
dim(pf2)

# find change in methylation
identical(meta$replicate,colnames(d))
d2 <- d[,rep_order]
delta_list <- pbmclapply(X=1:nrow(d2),FUN=function(i){
  dg <- data.frame(meta,methyl = d2[i,]) %>% group_by(sample) %>% summarise(avg_meth = mean(methyl))
  return(data.frame(region = rownames(d2)[i],
                    dssDelta = dg$avg_meth[dg$sample == 'DSS_DMSO'] - dg$avg_meth[dg$sample == 'Control_DMSO'], # change in methylation following DSS colitis and recovery (over control)
                    t5224Delta = dg$avg_meth[dg$sample == 'DSS_T5224'] - dg$avg_meth[dg$sample == 'DSS_DMSO']) # change in methylation following T-5224 AP-1 inhibition, within colitis recovered mice
         )
},mc.cores=cores)
df <- do.call('rbind',delta_list)
head(df)
write.table(df,paste0(output,'.delta_stats.all_peaks.txt'),sep='\t',quote=F,row.names=F)

# find background peaks
bg <- as.matrix(read.table(bg_file,header=F,sep=','))
dim(bg)
dim(df) # rows should match

df$is_var <- ifelse(df$region %in% rownames(pf2),'var','not_var')
df$is_bg <- 'not_bg'
df$is_bg[unique(c(bg[(df$is_var == 'var'),]))] <- 'bg'
head(df)

dp <- rbind(df %>% filter(is_var == 'var') %>% select(region,dssDelta,t5224Delta) %>% mutate(type = 'var'),
            df %>% filter(is_bg == 'bg') %>% select(region,dssDelta,t5224Delta) %>% mutate(type = 'bg')
            )

dp$type <- factor(dp$type,levels=c('var','bg'))

# Wilcox test for change in selected peaks relative to background set
wtest1 <- wilcox.test(x = dp %>% filter(type == 'var') %>% pull(dssDelta),
                      y = dp %>% filter(type == 'bg') %>% pull(dssDelta),
                      alternative = 'two.sided'
                      )
wtest2 <- wilcox.test(x = dp %>% filter(type == 'var') %>% pull(t5224Delta),
                      y = dp %>% filter(type == 'bg') %>% pull(t5224Delta),
                      alternative = 'two.sided'
)

# plot

maxY <- max(dp[,c('dssDelta','t5224Delta')],na.rm=T)
minY <- min(dp[,c('dssDelta','t5224Delta')],na.rm=T)

g1 <- ggplot(dp,aes(x=type,y=dssDelta,fill=type)) + geom_violin(trim=T) + geom_boxplot(outlier.shape=NA,width=0.25) + 
  theme_bw() + guides(fill='none') + ylim(minY,maxY) +
  xlab('') + ylab('Change in methylation fraction') + ggtitle('DSS DMSO vs Control DMSO') + 
  annotate('text',x=1.5,y=maxY*0.9,label=paste0('pval=',formatC(wtest1$p.value,digits=2,format='e')))

g2 <- ggplot(dp,aes(x=type,y=t5224Delta,fill=type)) + geom_violin(trim=T) + geom_boxplot(outlier.shape=NA,width=0.25) +  
  theme_bw() + guides(fill='none') + 
  xlab('') + ylab('Change in methylation fraction') + ggtitle('DSS T5224 vs DSS DMSO') + 
  annotate('text',x=1.5,y=maxY*0.9,label=paste0('pval=',formatC(wtest2$p.value,digits=2,format='e')))

p1 <- plot_grid(g1,g2,nrow=1)

pdf(paste0(output,'.AP1_sites.methyl_change.nvar_',nvar,'.pdf'),width=5,height=4)
print(p1)
dev.off()

