output <- '' #prefix for output files
umiFile <- '' #filtered UMI counts from 01a_refilter_barcodes.RNA
qc_file <- '' #counts.csv.gz file from SHARE-seq V2 workflow

two_thresh <- F #use two size thresholds for filtering
# this was used for tissue datasets since we see significantly lower RNA counts from immune cells than epithelial
if (two_thresh){
  thresh_marker <- 'Ptprc' # marker to split cells for different thresholding
  umi_thresh1 <- 300
  umi_thresh2 <- 800
}

cores <- 4
umi_thresh <- 500 # 300 for mouse colitis tissue, 500 for mouse organoids, mouse adenoma tissue and human organoids
doub_est <- 0.06
umap_pt_size <- 0.5
zCap <- 3
nTop <- 50
markerList <- readLines('../ref/gut_primary.ms.txt') # convenient lists of marker genes to plot

source('../tools/singlet_clustering.shared.v1.7.R')
source('../tools/singlet_clustering.RNA.v1.4.R')

#############################################################
###################### 01 - DECLUMPING ######################
#############################################################
## We have found that large clumps of cells (>10) often do not get called properly by doublet callers.
## However, they consistently demonstrate markedly higher estimated library sizes than singlets or doublets
## This is most easily detected by a k-NN smoothening approach 

umi <- readRDS(umiFile)
cat('Mean transcripts per cell:',mean(Matrix::colSums(umi)),'\n')
umi_filt <- umi[,(Matrix::colSums(umi) > umi_thresh)]
cat('Number of cells passing UMI threshold',paste0('(',umi_thresh,')'),':',ncol(umi_filt),'out of',ncol(umi),'\n')

qc <- read.csv(qc_file,header=T)
colnames(qc)[grepl('_libsize',colnames(qc))] <- 'libsize'
qc$barcode <- with(qc,paste(R1,R2,R3,P5,sep=','))
if (sum(!(colnames(umi_filt) %in% qc$barcode)) > 0){
  cat('There are cells in the UMI matrix not found in QC files!\n')
}
qc_sub1 <- qc[(qc$barcode %in% colnames(umi_filt)),]
qc_sub1 <- qc_sub1[match(colnames(umi_filt),qc_sub1$barcode),]

if (two_thresh){
  keep.cells <- splitThreshold( counts = umi_filt,
                                sizes = qc_sub1$libsize,
                                thresh1 = umi_thresh1,
                                thresh2 = umi_thresh2,
                                marker = thresh_marker
                                )
  dev.off()
  cat(sum(keep.cells),'cells of',length(keep.cells),'total retained after split thresholding\n')
  umi_filt <- umi_filt[,keep.cells]
  qc_sub1 <- qc_sub1[keep.cells,]
}

if (ncol(umi_filt) > 50000){
  starts <- seq(1,ncol(umi_filt),by=10000)
  ends <- starts + 9999
  ends[length(ends)] <- ncol(umi_filt)
  chunkList <- mapply(c,starts,ends,SIMPLIFY=F)
  
  norm_list <- pbmclapply(X=chunkList,FUN=function(x){
    chunk <- c(x[1]:x[2])
    norm_sub <- t(t(umi_filt[,chunk])/(Matrix::colSums(umi_filt[,chunk])))
    return(norm_sub)
  },mc.cores=cores)
  norm1 <- do.call('cbind',norm_list)
  rm(norm_list)
  gc()
} else {
  norm1 <- t(t(umi_filt)/(Matrix::colSums(umi_filt))) # normalized to total transcripts per cell
}
exp <- as(norm1*mean(Matrix::colSums(umi_filt)),'sparseMatrix') # scale to average transcripts across all cells
saveRDS(exp,paste0(output,'.min_umi_',umi_thresh,'.norm_UMI.all.rds'))

# filter for variant / expressed genes
exp.sd <- rowSds(exp)
names(exp.sd) <- rownames(exp)
var <- exp.sd[order(exp.sd,decreasing=T)]
l2e <- log2(t(exp[names(var)[1:5000],])+1)

# PCA
PCA <- cachePCA(	cachePath = './pcaCache',
					dataSet = as.matrix(l2e),
					center = F, scale = F)
pc.scores <- data.frame(PCA$x[,1:50])
saveRDS(pc.scores,paste0(output,'.min_umi_',umi_thresh,'.var_5000.l2e_pc_scores.all.rds'))
saveRDS(data.frame(PCA$rotation[,1:50]),paste0(output,'.min_umi_',umi_thresh,'.var_5000.l2e_pc_rotations.all.rds'))

# top 20 PC UMAP
UMAP <- cacheUMAP(	cachePath = './umapCache',
					dataSet = pc.scores[,1:20],
					seed = 123)
umap.p <- data.frame(UMAP$layout)
colnames(umap.p) <- c('UMAP1','UMAP2')
saveRDS(umap.p,paste0(output,'.min_umi_',umi_thresh,'.pc_20.umap_coord.all.rds'))
rm(umi,PCA,UMAP)
gc(verbose=FALSE)

## Finding KNN
knn_data <- get.knn(pc.scores[,1:20],k=20)
knn <- knn_data$nn.index
write.table(knn,paste0(output,'.pc20_knn_20.all.csv'),col.names=F,row.names=F,quote=F,sep=',')

# QC filtering
qc_sub2 <- qc_sub1[(qc_sub1$barcode %in% rownames(pc.scores)),]
qc_sub2 <- qc_sub2[match(rownames(pc.scores),qc_sub2$barcode),]
qc_sub2$log10UMI <- log10(Matrix::colSums(umi_filt))

# identifying clumps
pdf(paste0(output,'.umap_all.log10UMI.pdf'),width=4.5,height=4)
g <- ggplot(cbind(qc_sub2,umap.p),aes(x=UMAP1,y=UMAP2,color=log10UMI)) + geom_point(size=umap_pt_size,stroke=0) + 
	scale_color_gradient(low='navy',high='yellow') +
		theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
print(g)
dev.off()
clump.ind <- getShareClumps(qc_sub2$libsize,knn,umap.p,ptSize=1)
writeLines(rownames(umap.p)[clump.ind],paste0(output,'.RNA_clump_barcodes.txt'))
writeLines(rownames(umap.p)[!(clump.ind)],paste0(output,'.RNA_declump_barcodes.txt'))
umi_declump <- umi_filt[,!(clump.ind)]
saveRDS(umi_declump,paste0(output,'.min_umi_',umi_thresh,'.raw_UMI.declump.rds'))

if (ncol(umi_declump) > 50000){
  starts <- seq(1,ncol(umi_declump),by=10000)
  ends <- starts + 9999
  ends[length(ends)] <- ncol(umi_declump)
  chunkList <- mapply(c,starts,ends,SIMPLIFY=F)
  
  norm_list <- pbmclapply(X=chunkList,FUN=function(x){
    chunk <- c(x[1]:x[2])
    norm_sub <- t(t(umi_declump[,chunk])/(Matrix::colSums(umi_declump[,chunk])))
    return(norm_sub)
  },mc.cores=cores)
  norm_declump <- do.call('cbind',norm_list)
  rm(norm_list)
  gc()
} else {
  norm_declump <- t(t(umi_declump)/(Matrix::colSums(umi_declump))) # normalized to total transcripts per cell
}
exp_declump <- as(norm_declump*mean(Matrix::colSums(umi_declump)),'sparseMatrix') # scale to average transcripts across all cells
saveRDS(exp_declump,paste0(output,'.min_umi_',umi_thresh,'.norm_UMI.declump.rds'))
rm(umi_filt,exp,norm1,l2e)
gc(verbose=FALSE)

# variant genes, PCA, UMAP on declumped data
exp_declump.sd <- rowSds(exp_declump)
names(exp_declump.sd) <- rownames(exp_declump)
var_declump <- exp_declump.sd[order(exp_declump.sd,decreasing=T)]
l2e_declump <- log2(t(exp_declump[names(var_declump)[1:5000],])+1)
PCA.d <- cachePCA(	cachePath = './pcaCache',
                   dataSet = as.matrix(l2e_declump),
                   center = F, scale = F)
pc.scores.d <- data.frame(PCA.d$x[,1:50])
saveRDS(pc.scores.d,paste0(output,'.min_umi_',umi_thresh,'.var_5000.l2e_pc_scores.declump.rds'))
saveRDS(data.frame(PCA.d$rotation[,1:50]),paste0(output,'.min_umi_',umi_thresh,'.var_5000.l2e_pc_rotations.declump.rds'))

UMAP.d <- cacheUMAP(	cachePath = './umapCache',
                     dataSet = pc.scores.d[,1:20],
                     seed = 123)
umap.d <- data.frame(UMAP.d$layout)
colnames(umap.d) <- c('UMAP1','UMAP2')
saveRDS(umap.d,paste0(output,'.min_umi_',umi_thresh,'.pc_20.umap_coord.declump.rds'))
rm(PCA.d,UMAP.d)
gc(verbose=FALSE)

#############################################################
################ 02 - DOUBLETS / COLLISIONS #################
#############################################################
## Once clumps are removed, doublet callers are better able to identify multiple cells

# scrublet to identify doublets
mtx_file <- paste0(output,'.min_umi_',umi_thresh,'.raw_UMI.declump.mtx')
writeMM(t(umi_declump),file=mtx_file)
doublet_info <- runPyScrublet(mtx_file,doub_est)
rownames(doublet_info) <- colnames(umi_declump)
write.table(doublet_info,paste0(output,'.min_umi_',umi_thresh,'.scrublet_est_',doub_est,'.doublet_info.txt'),col.names=F,row.names=T,quote=F,sep='\t')
doub <- findDoublets(doublet_info[,1],umap.d,ptSize=1)
doub.index <- doub[[1]]
doub.cutoff <- doub[[2]]
writeLines(rownames(umap.d)[doub.index],paste0(output,'.RNA_doublets_',round(doub.cutoff,2),'.barcodes.txt'))

#############################################################
################### 03 - SINGLET ANALYSIS ###################
#############################################################

# Re-run PCA, UMAP on singlets
umi.s <- umi_declump[,!(doub.index)]
saveRDS(umi.s,paste0(output,'.min_umi_',umi_thresh,'.raw_UMI.singlet.rds'))

if (ncol(umi.s) > 50000){
  starts <- seq(1,ncol(umi.s),by=10000)
  ends <- starts + 9999
  ends[length(ends)] <- ncol(umi.s)
  chunkList <- mapply(c,starts,ends,SIMPLIFY=F)
  
  norm_list <- pbmclapply(X=chunkList,FUN=function(x){
    chunk <- c(x[1]:x[2])
    norm_sub <- t(t(umi.s[,chunk])/(Matrix::colSums(umi.s[,chunk])))
    return(norm_sub)
  },mc.cores=cores)
  norm.s <- do.call('cbind',norm_list)
  rm(norm_list)
  gc()
} else {
  norm.s <- t(t(umi.s)/(Matrix::colSums(umi.s))) # normalized to total transcripts per cell
}
exp.s <- as(norm.s*mean(Matrix::colSums(umi.s)),'sparseMatrix') # scale to average transcripts across all cells
saveRDS(exp.s,paste0(output,'.min_umi_',umi_thresh,'.norm_UMI.singlet.rds'))
rm(umi_declump,norm_declump,exp_declump,l2e_declump,norm.s)
gc(verbose=FALSE)

exp.s.sd <- rowSds(exp.s)
names(exp.s.sd) <- rownames(exp.s)
var.s <- exp.s.sd[order(exp.s.sd,decreasing=T)]
l2e.s <- log2(t(exp.s[names(var.s)[1:5000],])+1)

PCA.s <- cachePCA(	cachePath = './pcaCache',
                 dataSet = as.matrix(l2e.s),
                 center = F, scale = F)
pc.scores.s <- data.frame(PCA.s$x[,1:50])
saveRDS(pc.scores.s,paste0(output,'.min_umi_',umi_thresh,'.var_5000.l2e_pc_scores.singlet.rds'))
saveRDS(data.frame(PCA.s$rotation[,1:50]),paste0(output,'.min_umi_',umi_thresh,'.var_5000.l2e_pc_rotations.singlet.rds'))
rm(PCA.s,l2e.s)
gc(verbose=FALSE)

qc_sub.s <- qc_sub2[(qc_sub2$barcode %in% rownames(pc.scores.s)),]
qc_sub.s <- qc_sub.s[match(rownames(pc.scores.s),qc_sub.s$barcode),]
colnames(qc_sub.s)[grepl('_unique',colnames(qc_sub.s))] <- 'unique'
## plot singlet QC stats
q1 <- ggplot(qc_sub.s,aes(x=1,y=log10(libsize)))+geom_violin(fill='dodgerblue')+xlab('')
detected_genes <- Matrix::colSums(umi.s)
detected_genes <- detected_genes[match(qc_sub.s$barcode,names(detected_genes))]
qc_sub.s$det_genes <- detected_genes
q2 <- ggplot(qc_sub.s,aes(x=1,y=log10(det_genes)))+geom_violin(fill='red')+xlab('')
q3 <- ggplot(qc_sub.s,aes(x=unique,y=det_genes))+geom_point() + xlab('UMIs detected')+ylab('Genes detected')
clump_stats <- data.frame(n_cells=c(sum(clump.ind),sum(doub.index),ncol(umi.s)),type=c('clumps','doublets','singlets'))
clump_stats$frac_cells <- clump_stats$n_cells/sum(clump_stats$n_cells)
q4 <- ggplot(clump_stats,aes(y=n_cells,x=1,fill=type))+geom_bar(stat='identity',position='stack')+coord_flip()+xlim(0,2)+xlab('')+ylab('Number of cells')
q5 <- ggplot(clump_stats,aes(y=frac_cells,x=1,fill=type))+geom_bar(stat='identity',position='stack')+coord_flip()+xlim(0,2)+xlab('')+ylab('Fraction of cells')
pdf(paste0(output,'.min_umi_',umi_thresh,'.qc_stats.singlet.pdf'),width=8,height=5)
plot_grid(q1,q2,q3,nrow=1)
with(qc_sub.s,plot(log10(qc_sub.s$libsize)[order(qc_sub.s$libsize,decreasing=T)],col='firebrick',ylab='log10(library size)',xlab='Rank',main='Singlets - RNA'))
plot_grid(q4,q5,ncol=1)
dev.off()

## This is to generate metdata, can just load from publicly available files if preferred
r1.s <- as.numeric(stri_split_fixed(qc_sub.s$R1,'.',simplify=T)[,2])
qc_sub.s$sample <- 'blank' # e.g., Control, Acute, Chronic, Recovered
qc_sub.s$replicate <- 'blank' # e.g., DSS04-DSS-m1...
qc_sub.s$batch <- 'blank' # e.g., DSS04

## e.g. 
## qc_sub.s$sample <- 'blank'
## qc_sub.s$sample[(r1.s<=16)] <- 'd35.1'
## qc_sub.s$sample[(r1.s>16) & (r1.s<=48)] <- 'Control-P1'

table(qc_sub.s$batch)
table(qc_sub.s$sample)
table(qc_sub.s$replicate)

write.table(qc_sub.s,paste0(output,'.min_umi_',umi_thresh,'.qc_stats.singlet.txt'),quote=F,sep='\t',row.names=F)

### check for depth in PCs
depth.check <- pcaDepthCor(pc.scores.s,qc_sub.s$libsize)

cors <- depth.check[[1]]
usePC1 <- depth.check[[2]]

pdf(paste0(output,'.min_umi_',umi_thresh,'.var_5000.l2e_pc_scores.depth_cor.singlet.pdf'),width=5,height=4)
plot(cors[1:20],xlab='Principal Component',ylab='Correlation with log10 library size')
lines(cors[1:20])
text(2,0.9*max(cors[1:20]),labels=paste0('Use PC1: ',usePC1))
dev.off()

if (usePC1 == 'n' | usePC1 == 'N'){
  cat('Excluding PC1 from clustering\n')
  p <- pc.scores.s[,2:21]
} else {
  p <- pc.scores.s[,1:20]
}

UMAP.s <- cacheUMAP(	cachePath = './umapCache',
                   dataSet = p,
                   seed = 123)
umap.s <- data.frame(UMAP.s$layout)
colnames(umap.s) <- c('UMAP1','UMAP2')
saveRDS(umap.s,paste0(output,'.min_umi_',umi_thresh,'.pc_20.umap_coord.singlet.rds'))
rm(UMAP.s)
gc(verbose=FALSE)

# Smoothening KNN
knn_data.s20 <- get.knn(p,k=20)
knn.s20 <- knn_data.s20$nn.index
write.table(knn.s20,paste0(output,'.pc20_knn_20.singlet.csv'),col.names=F,row.names=F,quote=F,sep=',')

knn_data.s <- get.knn(p,k=5)
knn.s <- knn_data.s$nn.index
write.table(knn.s,paste0(output,'.pc20_knn_5.singlet.csv'),col.names=F,row.names=F,quote=F,sep=',')

# plot depth
qc_sub.s <- cbind(qc_sub.s,umap.s)
pdf(paste0(output,'.umap_singlet.libsize.pdf'),width=4.5,height=4)
g <- ggplot((qc_sub.s),aes(x=UMAP1,y=UMAP2,color=log10(libsize))) + geom_point(size=umap_pt_size,stroke=0) + 
  scale_color_gradient(low='navy',high='yellow') +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
print(g)
dev.off()

# plot markers
exp.s.smooth <- smoothCellsNN(E = exp.s, K = knn.s, cores = cores)
saveRDS(exp.s.smooth,paste0(output,'.min_umi_',umi_thresh,'.norm_UMI.knn_5_smooth.singlet.rds'))
gc(verbose=FALSE)

pdf(paste0(output,'.umap_singlet.marker_genes.pdf'),width=4.5,height=4)
for (marker in markerList){
  cat(marker,'\n')
  plotUMIgene(umap.s,exp.s.smooth,marker,zCap,umap_pt_size,palette='brewer_heat')
  #cont <- readline(prompt = 'Enter')
}
dev.off()


# louvain clustering
ig <- graph.empty(nrow(knn.s))
edge_list <- pbmclapply(X=1:nrow(knn.s),FUN=function(x){
  return(unlist(mapply(c,rep(x,ncol(knn.s)),knn.s[x,],SIMPLIFY=F)))
},mc.cores=cores)
ig <- add_edges(ig,unlist(edge_list))
comm <- cluster_louvain(as.undirected(ig))
clust <- data.frame(sample=rownames(p),community=as.vector(membership(comm)),stringsAsFactors = F)
write.table(clust,paste(output,'pc20_singlet.louvain_clusters.txt',sep='.'),quote=F,sep='\t',row.names=F)

u <- cbind(umap.s,as.character(clust$community))
colnames(u)[3] <- 'community'

pdf(paste0(output,'.umap_singlet.louvain_clusters.pdf'),width=5,height=4)
g <- ggplot(u,aes(x=UMAP1,y=UMAP2,color=community)) + geom_point(size=umap_pt_size,stroke=0) + 
  scale_color_viridis(discrete = T, option = "D")+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
print(g)
for (comm in unique(u$community)){
  g <- ggplot(u,aes(x=UMAP1,y=UMAP2,color=(community==comm))) + geom_point(size=umap_pt_size,stroke=0) +
    guides(color=guide_legend(title=paste0('Community ',comm))) +
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
  print(g)
}
dev.off()