output <- 'APC-Visium-02'
#smoothened output from 02_spatial
smooth_files <- list('sample1' = '/path/norm.smooth.rds',
                     'sample2' = '/path/norm.smooth.rds'
                     )

feature_to_use <- 'Axin2'
start_thresh <- 2
ncores <- 4

library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(BuenColors)
library(FNN)
library(igraph)
library(pbmcapply)
library(viridis)
#######

#######
smooth_list <- lapply(X=smooth_files,FUN=function(file_name){return(readRDS(file_name))})
names(smooth_list) <- names(smooth_files)

#### Call individual tumor cells 
# iterate across each sample, call tumor cells - should use same threshold across samples
d_list <- list()
for (i in 1:length(smooth_list)){
  samp_name <- names(smooth_list)[i]
  cat('--',samp_name,'--\n')
  
  obj <- smooth_list[[i]]
  exp.s <- obj$pc_smooth
  coord <- obj$spatial
  df <- data.frame(X = coord$X, Y = coord$Y, exp = scale(exp.s[feature_to_use,])[,1])
  
  # plot scaled value
  thresh <- start_thresh
  cont <- 'y'
  while (cont == 'Y' | cont == 'y'){
    cat('Threshold:',thresh,'\n')
    d <- df %>% mutate(isTumor = ifelse(exp > thresh,'yes','no'))
    cat('Percent of cells called as tumor:',100*(round(mean(d$isTumor == 'yes'),3)),'\n')
    g1 <- ggplot(d,aes(x=exp)) + geom_histogram(bins = 50) + theme_bw() + 
      xlab('Scaled expression') + ggtitle('Expression Range') +
      geom_vline(xintercept = thresh, linetype='dashed', color = 'red')
    g2 <- ggplot(d,aes(x=X,y=Y,color=exp)) + geom_point(size = 0.1) + 
      theme_bw() + xlab('') + ylab('') + ggtitle('Spatial Expression') +
      coord_fixed() + labs(color = feature_to_use) + 
      scale_color_gradientn(colors = jdb_palette('brewer_heat'))
    g3 <- ggplot(d,aes(x=X,y=Y,color=isTumor)) + geom_point(size = 0.1) + 
      theme_bw() + xlab('') + ylab('') + ggtitle('Tumor Calls') + 
      coord_fixed() + 
      scale_color_manual(values = c('yes' = 'black', 'no' = 'gray75'))
    p1 <- plot_grid(g1,g2,g3,nrow=1)
    print(p1)
    cont <- readline(prompt = 'Set new threshold (y/n)?: ')
    if (cont == 'y' | cont == 'Y'){
      thresh <- as.numeric(readline(prompt = 'New threshold to use: '))
    }
  }
  
  d_list[[samp_name]] <- d %>% mutate(sample = samp_name)  
}
df <- do.call('rbind',d_list)
write.table(df,paste0(output,'.tumor_cell_calls.txt'),sep='\t',quote=F)

### Cluster cells into tumors
uniq_samp <- unique(df$sample)
K <- 5 
clust_list <- list()
for (samp in uniq_samp){
  cat(samp,'\n')
  dsub <- df %>% filter(sample == samp)
  dtumor <- dsub %>% filter(isTumor == 'yes')
  
  knn_data.tumor <- get.knn(dtumor[,c("X",'Y')],k=K)
  knn.tumor <- knn_data.tumor$nn.index
  
  ig <- graph.empty(nrow(knn.tumor))
  edge_list <- pbmclapply(X=1:nrow(knn.tumor),FUN=function(x){
    return(unlist(mapply(c,rep(x,ncol(knn.tumor)),knn.tumor[x,],SIMPLIFY=F)))
  },mc.cores=ncores)
  ig <- add_edges(ig,unlist(edge_list))
  comm <- cluster_louvain(as.undirected(ig))
  clust <- cbind(dtumor,cluster = as.vector(membership(comm)))
  
  clust <- rbind(clust,(dsub %>% filter(isTumor == 'no') %>% mutate(cluster = 'none')))
  
  ge <- ggplot(dsub,aes(x=X,y=Y,color=exp)) + geom_point(size = 0.1) + 
    theme_bw() + xlab('') + ylab('') + ggtitle(paste0(samp,'- Spatial Expression')) +
    coord_fixed() + labs(color = feature_to_use) + 
    scale_color_gradientn(colors = jdb_palette('brewer_heat'))
  gt <- ggplot(dsub %>% arrange(isTumor),aes(x=X,y=Y,color=isTumor)) + geom_point(size = 0.1) + 
    theme_bw() + xlab('') + ylab('') + ggtitle(paste0(samp,'- Tumor Cells')) + 
    coord_fixed() + 
    scale_color_manual(values = c('yes' = 'black', 'no' = 'gray75'))
  gc <- ggplot(clust %>% filter(isTumor == 'no'),aes(x=X,y=Y,color=cluster)) + geom_point(size = 0.1, col = 'gray75') + 
    theme_bw() + xlab('') + ylab('') + ggtitle(paste0(samp,'- Tumor Clusters')) + 
    coord_fixed() +
    geom_point(data = clust %>% filter(isTumor == 'yes'), aes(x=X,y=Y, color=cluster),size=0.1) + 
    scale_color_viridis_d() + guides(color = 'none')
  
  p2 <- plot_grid(ge,gt,gc,nrow=1)
  
  pdf(paste0(output,'.tumor_clusters.k_',K,'.',samp,'.pdf'),width=15,height=5)
  print(p2)
  #cont <- readline(prompt = 'Enter')
  dev.off()
  
  clust_list[[samp]] <- clust
}
dc <- do.call('rbind',clust_list)
write.table(dc,paste0(output,'.tumor_clusters.k_',K,'.txt'),sep='\t',quote=F)



