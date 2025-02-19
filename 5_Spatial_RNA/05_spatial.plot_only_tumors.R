## for plotting the relative expression of single genes in only tumor cells
## AP-1/P20 program expression was not smoothened or re-normalized for plotting

output <- ''
#output from 02_spatial
raw_files <- list('sample1' = '/path/raw.rds',
                  'sample2' = '/path/raw.rds'
                  )

pc_files <- c('sample1' = '/path/pc_scores.rds',
              'sample2' = '/path/pc_scores.rds'
              )

coord_files <- c('sample1' = '/path/norm.filtered.rds',
                 'sample2' = '/path/norm.filtered.rds'
                 )

tumor_file <- '' # tumor calls from 03_spatial

ncores <- 8
Kp <- 10 
Ks <- 10

library(Matrix)
library(FNN)
library(pbmcapply)
library(stringi)
library(dplyr)

#######
smoothCellsNN <- function( E, #expression matrix, genes x cells
                           K, # KNN matrix, each row = index of KNN
                           cellThresh = 20000, #threshold for splitting up large matrices
                           nChunks = 4, #for splitting up large matrices
                           cores = 4
){
  nCells <- ncol(E)
  if (nCells > cellThresh){
    chunkSize <- ceiling(nCells/nChunks)
    cat('Large matrix detected (greater than 20k cells), smoothening in chunks of',chunkSize,'cells\n')
    starts <- seq(1,nCells,by=chunkSize)
    ends <- starts + chunkSize - 1 
    ends[length(ends)] <- nCells
    chunkList <- mapply(c,starts,ends,SIMPLIFY=F)
    E.smooth.list <- vector('list',length=nChunks)
    for (i in 1:length(chunkList)){
      cat('Chunk',i,'\n')
      chunk <- chunkList[[i]]
      smooth.list.sub <- pbmclapply(X=chunk[1]:chunk[2],FUN=function(x){
        return(Matrix::rowMeans(cbind(E[,x],E[,K[x,]])))
      },mc.cores=cores)
      E.smooth.list[[i]] <- do.call('cbind',smooth.list.sub)
      E.smooth.list[[i]] <- as(E.smooth.list[[i]],'sparseMatrix')
      rm(smooth.list.sub)
      gc(verbose=F)
    }
    cat('Combining chunks...\n')
    E.smooth <- do.call('cbind',E.smooth.list)
    colnames(E.smooth) <- colnames(E)
    return(E.smooth)
  } else {
    smooth.list <- pbmclapply(X=1:ncol(E),FUN=function(x){
      return(Matrix::rowMeans(cbind(E[,x],E[,K[x,]])))
    },mc.cores=cores)
    E.smooth <- do.call('cbind',smooth.list)
    E.smooth <- as(E.smooth,'sparseMatrix')
    colnames(E.smooth) <- colnames(E)
    return(E.smooth)
  }
}
##########

## read data, subset cells
raw_list <- lapply(X=raw_files,FUN=function(x){return(readRDS(x))})
count_list <- lapply(X=raw_list,FUN=function(x){return(x$counts)})
names(count_list) <- c('A1-13','A1-7')
rm(raw_list)
gc()

tumor_calls <- read.table(tumor_file,header=T,sep='\t')
tumor_calls$barcode <- stri_split_fixed(rownames(tumor_calls),'.',simplify=T)[,3]
tumor_counts <- lapply(X=1:length(count_list),FUN=function(i){
  tumor_barcodes <- tumor_calls %>% filter(sample == names(count_list)[i]) %>% filter(isTumor == 'yes') %>% pull(barcode)
  counts_sub <- count_list[[i]]
  return(counts_sub[,colnames(counts_sub) %in% tumor_barcodes])
})
names(tumor_counts) <- names(count_list)
saveRDS(tumor_counts,paste0(output,'.raw.filtered.rds'))

avg_counts <- mean(Matrix::colSums(do.call('cbind',tumor_counts)))

## renormalize across only tumor cells
renorm_list <- lapply(X=tumor_counts,FUN=function(counts){
  renorm_counts_list <- pbmclapply(X=1:ncol(counts),FUN=function(i){
    return(counts[,i]*avg_counts/sum(counts[,i]))
  },mc.cores=ncores)
  renorm_counts <- do.call('cbind',renorm_counts_list)
  colnames(renorm_counts) <- colnames(counts)
  renorm_counts <- as(renorm_counts,'sparseMatrix')
  return(renorm_counts)
})
saveRDS(renorm_list,paste0(output,'.joint_norm.filtered.rds'))
rm(tumor_counts)
gc()

## find new PC KNN
pc_list <- lapply(X=1:length(pc_files),FUN=function(i){
  pc_all <- readRDS(pc_files[i])
  return(pc_all[match(colnames(renorm_list[[i]]),rownames(pc_all)),])
})

knn_list <- pbmclapply(X=pc_list,FUN=function(d){
  knn_data <- get.knn(d[,1:20],k=Kp)
  return(knn_data$nn.index)
},mc.cores=ncores)
names(knn_list) <- names(renorm_list)
saveRDS(knn_list,paste0(output,'.pc_knn_',Kp,'.rds'))

# smoothen
smooth_list <- list()
for (i in 1:length(renorm_list)){
  d.s <- smoothCellsNN(E = renorm_list[[i]], K = knn_list[[i]], cores = ncores)
  smooth_list[[names(renorm_list)[i]]] <- d.s
}
saveRDS(smooth_list,paste0(output,'.smooth.pc20_knn_',Kp,'.rds'))

## find new spatial KNN
coord_list <- lapply(X=1:length(coord_files),FUN=function(i){
  obj <- readRDS(coord_files[[i]])
  return(obj$spatial[match(colnames(smooth_list[[i]]),obj$spatial$bc),])
})
knn_list2 <- pbmclapply(X=coord_list,FUN=function(d){
  knn_data2 <- get.knn(d[,2:3],k=Ks)
  return(knn_data2$nn.index)
},mc.cores=ncores)
names(knn_list2) <- names(smooth_list)
saveRDS(knn_list2,paste0(output,'.spatial_knn_',Ks,'.rds'))


# smoothen
smooth_list2 <- list()
for (i in 1:length(smooth_list)){
  d2.s <- smoothCellsNN(E = smooth_list[[i]], K = knn_list2[[i]], cores = ncores)
  smooth_list2[[names(smooth_list)[i]]] <- d2.s
}
saveRDS(smooth_list2,paste0(output,'.smooth.pc20_knn_',Kp,'.spatial_knn_',Ks,'.rds'))


# plot 
library(ggplot2)
library(cowplot)
library(BuenColors)
gene_list <- c('Tnfrsf11b','Macf1','Cdh13','Fn1','Pmepa1','Thbs1','Mllt3','Dock4','Runx2','Clu','Srgap1')
cap <- 3

coord_list2 <- lapply(X=1:length(coord_files),FUN=function(i){
  obj <- readRDS(coord_files[[i]])
  return(obj$spatial)
})

for (gene in gene_list){
  pdf(paste0(output,'.',gene,'.pdf'),width=12,height=5)
  cat(gene,'\n')
  d_list <- lapply(X=1:length(smooth_list2),FUN=function(i){
    samp <- names(smooth_list2)[i]
    exp.t <- smooth_list2[[i]]
    coord_all <- coord_list2[[i]]
    df.t <- data.frame(sample = samp,
                       X = coord_all$X[match(colnames(exp.t),coord_all$bc)],
                       Y = coord_all$Y[match(colnames(exp.t),coord_all$bc)],
                       exp = exp.t[gene,],
                       isTumor = 'yes'
    )
    df.t$barcode <- rownames(df.t)
    coord_rest <- coord_all %>% filter(!(bc %in% colnames(exp.t)))
    df.rest <- data.frame(sample = samp,
                          X = coord_rest$X,
                          Y = coord_rest$Y,
                          exp = 0,
                          isTumor = 'no',
                          barcode = coord_rest$bc
    )
    return(rbind(df.t,df.rest))
  })
  d <- do.call('rbind',d_list)
  #d$norm <- scale(d$exp)[,1]
  dtumor <- d %>% filter(isTumor == 'yes')
  dtumor$norm <- scale(dtumor$exp)[,1]
  dnot <- d %>% filter(isTumor == 'no')
  
  uniq_samp <- unique(d$sample)
  g_list <- lapply(X=uniq_samp,FUN=function(samp){
    dt.s <- dtumor %>% filter(sample == samp)
    dn.s <- dnot %>% filter(sample == samp)
    
    dt.s$norm[dt.s$norm > cap] <- cap
    dt.s$norm[dt.s$norm < -cap] <- -cap
    
    g <- ggplot(dn.s,aes(x=X,y=Y)) + geom_point(size = 0.1,color = 'gray95') +
      theme_bw() + xlab('') + ylab('') + ggtitle(paste0(samp,' - Tumor Only')) +
      coord_fixed() + 
      geom_point(data = dt.s,aes(x=X,y=Y,color=norm),size=0.1) + 
      labs(color = gene) +
      scale_color_gradient2(low = 'dodgerblue',mid='beige', high = 'firebrick',midpoint=0)
    return(g)
  })
  p <- do.call('plot_grid',c(g_list,nrow=1))
  print(p)
  dev.off()
}

