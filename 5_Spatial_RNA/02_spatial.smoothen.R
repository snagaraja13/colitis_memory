output <- ''
# 16um binds were used for all the below
barcode_file <- 'barcodes.tsv.gz' # barcodes.tsv.gz from filtered spaceranger output
feature_file <- 'features.tsv.gz' # features.tsv.gz from filtered spaceranger output
count_file <- 'matrix.mtx.gz' # matrix.mtx.gz from filtered spaceranger output
spatial_file <- 'tissue_positions.parquet' # spatial/tissue_positions.parquet from filtered spaceranger output

k_spatial <- 20 # spatial smoothening window 
k_pca <- 20 # k-NN smoothening window 

cores <- 4

# favorite genes
marker_list <- c('Epcam','Acta2',
                 'Lgr5','Lrig1','Mki67',
                 'Ctnnb1','Axin2',
                 'Car1','Car4',
                 'Ptprc'
                 )

library(Matrix)
library(dplyr)
library(nanoparquet)
library(pbmcapply)
library(FNN)
library(ggplot2)
library(BuenRTools)
library(cowplot)

#####################
#### FUNCTIONS
#####################

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

plotSmoothVisium <- function(obj,gene,cap = T){
  d.pc <- obj$pc_smooth
  d.both <- obj$pc_spatial_smooth
  coord <- obj$spatial
  
  if (!identical(coord$bc,colnames(d.pc)) | !identical(coord$bc,colnames(d.both))){
    cat('Coordinate barcodes do not match count matrix rows, exiting...\n')
    return()
  }
  
  if (!(gene %in% rownames(d.pc))){
    cat(paste0(gene,'not present in expression matrix\n'))
  }
  
  df <- data.frame(X = coord$X, Y = coord$Y,
                   pc = scale(d.pc[gene,])[,1],
                   both = scale(d.pc[gene,])[,1]
  )
  
  if (cap){
    df$pc[df$pc > 3] <- 3
    df$pc[df$pc < (-3)] <- (-3)
    df$both[df$both > 3] <- 3
    df$both[df$both < (-3)] <- (-3)
  }
  
  g1 <- ggplot(df,aes(x=X,y=Y,color=pc)) + geom_point(size = 0.25) + 
    scale_color_gradient2(low='navy',mid='beige',high='firebrick',midpoint=0) + 
    labs(color = gene) +
    theme_bw() + xlab('') + ylab('') + ggtitle('PC Smoothened')
  
  g2 <- ggplot(df,aes(x=X,y=Y,color=both)) + geom_point(size = 0.25) + 
    scale_color_gradient2(low='navy',mid='beige',high='firebrick',midpoint=0) + 
    labs(color = gene) +
    theme_bw() + xlab('') + ylab('') + ggtitle('PC and Spatial Smoothened')
  
  p <- plot_grid(g1,g2,nrow=1)
  print(p)
}
#####################
#### END FUNCTIONS
#####################

## read in data, reformat
counts <- Matrix::readMM(count_file)
rownames(counts) <- read.table(feature_file,header=F,sep='\t') %>% pull(V2)
colnames(counts) <- readLines(barcode_file)
counts[1:5,1:5]

spatial_raw <- nanoparquet::read_parquet(spatial_file)
head(spatial_raw)
spatial <- spatial_raw[match(colnames(counts),spatial_raw$barcode),] %>% summarise(bc = barcode, X = pxl_row_in_fullres, Y = pxl_col_in_fullres) 
identical(spatial$bc,colnames(counts))
saveRDS(list(counts = counts, spatial = spatial),paste0(output,'.raw.all.rds'))

## filter cells based on total reads
counts_per_bead <- Matrix::colSums(counts)
cont <- 'Y'

while(cont == 'Y' | cont == 'y'){
  hist(log10(counts_per_bead+1))
  thresh <- as.numeric(readline(prompt = 'Minimum read cutoff: '))
  abline(v=log10(thresh+1),lty=2,col='red')
  cont <-readline(prompt = 'Change cut off? (y/n): ')
}
dev.off()

counts_filt <- counts[,(Matrix::colSums(counts) > thresh)]
spatial_filt <- spatial %>% filter(bc %in% colnames(counts_filt))
identical(spatial_filt$bc,colnames(counts_filt))
saveRDS(list(counts = counts_filt, spatial = spatial_filt),paste0(output,'.raw.filtered.rds'))
rm(counts,spatial)
gc()

## KNN smoothening
if (ncol(counts_filt) > 50000){
  starts <- seq(1,ncol(counts_filt),by=10000)
  ends <- starts + 9999
  ends[length(ends)] <- ncol(counts_filt)
  chunkList <- mapply(c,starts,ends,SIMPLIFY=F)
  
  norm_list <- pbmclapply(X=chunkList,FUN=function(x){
    chunk <- c(x[1]:x[2])
    norm_sub <- t(t(counts_filt[,chunk])/(Matrix::colSums(counts_filt[,chunk])))
    return(norm_sub)
  },mc.cores=cores)
  norm1 <- do.call('cbind',norm_list)
  rm(norm_list)
  gc()
} else {
  norm1 <- t(t(counts_filt)/(Matrix::colSums(counts_filt))) # normalized to total transcripts per cell
}
exp <- as(norm1*mean(Matrix::colSums(counts_filt)),'sparseMatrix') # scale to average transcripts across all cells
saveRDS(list(counts = exp, spatial = spatial_filt),paste0(output,'.norm.filtered.rds'))
l2e <- as(log2(exp[(Matrix::rowSums(exp,na.rm = T) > 0),]+1),'sparseMatrix')
rm(counts_filt, norm1)
gc()

PCA <- cachePCA(cachePath = './pcaCache',
                dataSet = as.matrix(t(l2e)),
                center = F, scale = F)
pc.scores <- data.frame(PCA$x[,1:50])
saveRDS(pc.scores,paste0(output,'.pc_scores.rds'))
rm(PCA,l2e)
gc()

knn_data.p <- get.knn(pc.scores[,1:20],k=k_pca)
knn.p <- knn_data.p$nn.index
write.table(knn.p,paste0(output,'.pc20_knn_',k_pca,'.csv'),col.names=F,row.names=F,quote=F,sep=',')

n.smooth <- smoothCellsNN(E = exp, K = knn.p, cores = cores)

## Smooth spatially
identical(spatial_filt$bc,rownames(exp))
knn_data.s <- get.knn(spatial_filt[,c(2,3)],k=k_spatial)
knn.s <- knn_data.s$nn.index
write.table(knn.s,paste0(output,'.spatial_',k_spatial,'_knn_',k_pca,'.csv'),col.names=F,row.names=F,quote=F,sep=',')

ns.smooth <- smoothCellsNN(E = n.smooth, K = knn.s, cores = cores)
exp.s <- list(pc_smooth = n.smooth, pc_spatial_smooth = ns.smooth, spatial = spatial_filt)

saveRDS(exp.s,paste0(output,'.norm.smooth_',k_spatial,'_',k_pca,'.rds'))

## ## Plot some genes
pdf(paste0(output,'.smooth_',k_spatial,'_',k_pca,'.marker_genes.pdf'),width=15,height=6)
for (gene in marker_list){
  cat(gene,'\n')
  plotSmoothVisium(exp.s,gene)
  #cont <- readline(prompt = 'Enter')
}
dev.off()













