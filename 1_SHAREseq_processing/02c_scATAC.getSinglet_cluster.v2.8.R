output <- '' # output prefix
genome <- 'mm10'
fragBed <- '' #combined fragments.tsv.gz file
bcFile <- '' #bc_pf from 01b_refilter_barcodes.ATAC
peakFile <- '' #fixed width and resized peaks - see /ref/ for peaks used in manuscript
qcFile <- '' #counts.csv.gz from SHARE V2 Workflow

# both declump and doublet cell barcodes to exclude from 02a_scRNA.getSinglet_cluster
rna_out_files <- c('',
                   ''
)

#for plotting genes on UMAP
rna_smooth_rds <- '' #smoothened RNA from  02a_scRNA.getSinglet_cluster
markerList <- readLines('../ref/')

## Ad barcode matching, c('P1.RNA','P1.ATAC'), see barcode mapping metadata files
adList <- list( c('','')
)

size.thresh <- 500 #starting threshold for filtering barcodes, can adjust manually later based on actual data
frip.thresh <- 0.2 #starting threshold for filtering barcodes, can adjust manually later based on actual data
# note - FRIP cutoffs are highly dependent on peak calls and should be carefully considered for each dataset
num_topics <- 80
umap_pt_size <- 0.5
zCap <- 3
cores <- 4
nBg <- 250
nTFplot <- 30
nTop <- 50
doGeneScores <- F

###

source('../tools/singlet_clustering.shared.v1.6.R')
source('../tools/singlet_clustering.ATAC.v1.6.R')

### General QC, filtering 

## Getting reads in peaks
peaks <- read.table(peakFile,stringsAsFactors=F)
colnames(peaks) <- c('chr','start','end')
peaks <- makeGRangesFromDataFrame(peaks)
bcList <- readLines(bcFile)
if (ncol(stri_split_fixed(bcList,'.',simplify=T)) == 8){
  cell_split <- data.frame(stri_split_fixed(bcList,'.',simplify=T))
  bcList <- with(cell_split,paste0(X1,'.',X2,',',X3,'.',X4,',',X5,'.',X6,',',X7,'.',X8))
}

fragments <- Signac::CreateFragmentObject(path = fragBed, cells = bcList, validate.fragments = FALSE)
fm <- Signac::FeatureMatrix(
  fragments,
  peaks,
  cells = bcList,
  verbose = TRUE
)
se <- SummarizedExperiment(list(counts=fm))
rowRanges(se) <- peaks

## Filtering cells by QC
qc <- read.csv(qcFile)

qc$barcode <- with(qc,paste(R1,R2,R3,P5,sep=','))
qc_sub <- qc[(qc$barcode %in% colnames(se)),]
colnames(qc_sub)[grepl('_libsize',colnames(qc_sub))] <- 'libsize'
colData(se)$libsize <- qc_sub$libsize[match(colnames(se),qc_sub$barcode)]
colData(se)$tot_reads <- qc_sub[,grepl('_unique',colnames(qc_sub))]
colData(se)$depth <- colSums(assays(se)$counts)
colData(se)$FRIP <- with(colData(se),(depth/tot_reads))

saveRDS(se,paste(output,'all_counts.rds',sep='.'))

pf <- plotAtacQC(se,var1='libsize',thresh1=size.thresh,frip.thresh=frip.thresh)
cat('Library size cutoff:',size.thresh,'\n')
cat('FRIP cutoff:',frip.thresh,'\n')
cont <- readline(prompt='Change cutoffs? (y/n): ')

## See methods section of manuscript for specific cutoffs used for each dataset
while (cont=='y' | cont=='Y'){
  size.thresh <- as.numeric(readline(prompt='New size threshold: '))
  frip.thresh <- as.numeric(readline(prompt='New FRIP threshold: '))
  pf <- plotAtacQC(se,var1='libsize',thresh1=size.thresh,frip.thresh=frip.thresh)
  cont <- readline(prompt='Change cutoffs? (y/n): ')
}
dev.off()
pdf(paste0(output,'.qc_plot.pdf'),width=4,height=4)
pf <- plotAtacQC(se,var1='libsize',thresh1=size.thresh,frip.thresh=frip.thresh)
dev.off()

se.qc <- se[,pf]
saveRDS(se.qc,paste0(output,'.lib_',size.thresh,'.frip_',frip.thresh,'.counts.all.rds'))
rm(se,fm)
gc()

#############################################################
################# 01 - ATAC ONLY DECLUMPING #################
#############################################################
### See notes in 02a_scRNA file regarding large clumps

## all cell cisTopics
cis.all <- getCisTopics(se.qc,'all',num_topics) # this can take a while to run (many hours) depending on number of cells being used
topics <- t(modelMatSelection(cis.all,target='cell',method='Z-score'))
write.table(topics,paste0(output,'.cisTopics_topics_zCell.all.txt'),sep='\t',quote=F)
write.table(cis.all@region.data,paste0(output,'.cisTopics_region_data.all.txt'),sep=',',quote=F)

## UMAP
UMAP <- cacheUMAP(	cachePath = './umapCache',
                   dataSet = topics,
                   seed = 123)
umap.t <- data.frame(UMAP$layout)
colnames(umap.t) <- c('UMAP1','UMAP2')
saveRDS(umap.t,paste0(output,'.cisTopic_',ncol(topics),'.umap_embed.all.rds'))
rm(UMAP)
gc()

## Finding KNN
knn_data <- get.knn(topics,k=20)
knn <- knn_data$nn.index

# QC filtering
qc_sub.all <- qc_sub[(qc_sub$barcode %in% rownames(topics)),]
qc_sub.all <- cbind(qc_sub.all,umap.t)

pdf(paste0(output,'.umap.log10libsize.all.pdf'),width=4.5,height=4)
g <- ggplot(qc_sub.all,aes(x=UMAP1,y=UMAP2,color=log10(libsize))) + geom_point(size=1,stroke=0) + 
  scale_color_gradient(low='navy',high='yellow') +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
print(g)
dev.off()

clump.ind <- getShareClumps(qc_sub.all$libsize,knn,umap.t,ptSize=1)
writeLines(rownames(umap.t)[clump.ind],paste0(output,'.ATAC_clump_barcodes.txt'))
se.d <- se.qc[,!(clump.ind)]
saveRDS(se.d,paste0(output,'.lib_',size.thresh,'.frip_',frip.thresh,'.counts.declump.rds'))

#############################################################
################ 02 - FILTERING RNA BARCODES ################
#############################################################
## Remove clumps and doublets identified in scRNA data

rna_bc_out_list <- vector("list",length=length(rna_out_files))
for (i in 1:length(rna_out_files)){
  rna_bc_out_list[[i]] <- readLines(rna_out_files[i])
}
rna_bc_out <- unique(unlist(rna_bc_out_list))

matchBC <- matchShareSamples( atac_bc = colnames(se.d),
                              rna_bc = rna_bc_out,
                              isMulti = isMultiSample,
                              adList = adList
                              )
keep_ind <- !(colnames(se.d) %in% matchBC[[1]])
sum(keep_ind)
se.s <- se.d[,keep_ind]
saveRDS(se.s,paste0(output,'.lib_',size.thresh,'.frip_',frip.thresh,'.counts.singlet.rds'))

#############################################################
################### 03 - SINGLET ANALYSIS ###################
#############################################################

cis.s <- getCisTopics(se.s,'singlet',nTopics = seq(ncol(topics)-10,ncol(topics)+10,by=10)) ## This can also take many hours to run
#topics.s <- cis.s[[1]]
topics.s <- t(modelMatSelection(cis.s,target='cell',method='Z-score'))
write.table(topics.s,paste0(output,'.cisTopics_topics_zCell.singlet.txt'),sep='\t',quote=F)
write.table(cis.s@region.data,paste0(output,'.cisTopics_region_data.singlet.txt'),sep=',',quote=F)

## UMAP
UMAP.s <- cacheUMAP(	cachePath = './umapCache',
                   dataSet = topics.s,
                   seed = 123)
umap.s <- data.frame(UMAP.s$layout)
colnames(umap.s) <- c('UMAP1','UMAP2')
saveRDS(umap.s,paste0(output,'.cisTopic_',ncol(topics.s),'.umap_embed.singlet.rds'))

## Finding KNN
knn_data.s20 <- get.knn(topics.s,k=20)
knn.s20 <- knn_data.s20$nn.index
write.table(knn.s20,paste0(output,'.cisTopics_knn_20.singlet.csv'),col.names=F,row.names=F,quote=F,sep=',')

knn_data.s <- get.knn(topics.s,k=5)
knn.s <- knn_data.s$nn.index
write.table(knn.s,paste0(output,'.cisTopics_knn_5.singlet.csv'),col.names=F,row.names=F,quote=F,sep=',')

# plot depth
qc_sub.s <- qc_sub[(qc_sub$barcode %in% rownames(topics.s)),]
qc_sub.s$FRIP <- colData(se.s)$FRIP

## Similar to RNA data, can important metdata if preferred. 
r1.s <- as.numeric(stri_split_fixed(qc_sub.s$R1,'.',simplify=T)[,2])
qc_sub.s$sample <- 'blank'
qc_sub.s$replicate <- 'blank'
qc_sub.s$batch <- 'blank'

# e.g.
#qc_sub.s$sample <- 'blank'
#qc_sub.s$rep[(r1.s<=24)] <- 'DSS04-DSS-m1'
# etc

table(qc_sub.s$batch)
table(qc_sub.s$sample)
table(qc_sub.s$replicate)

write.table(qc_sub.s,paste0(output,'.qc_stats.singlet.txt'),quote=F,sep='\t',row.names=F)

qc_sub.s <- cbind(qc_sub.s,umap.s)
pdf(paste0(output,'.umap.log10libsize.singlet.pdf'),width=4.5,height=4)
g <- ggplot(qc_sub.s,aes(x=UMAP1,y=UMAP2,color=log10(libsize))) + geom_point(size=1,stroke=0) + 
  scale_color_gradient(low='navy',high='yellow') +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank())
print(g)
dev.off()

rm(UMAP.s,knn_data.s,g)
gc()

#louvain
ig <- graph.empty(nrow(knn.s))
edge_list <- pbmclapply(X=1:nrow(knn.s),FUN=function(x){
  return(unlist(mapply(c,rep(x,ncol(knn.s)),knn.s[x,],SIMPLIFY=F)))
},mc.cores=cores)
ig <- add_edges(ig,unlist(edge_list))
comm <- cluster_louvain(as.undirected(ig))
clust <- data.frame(sample=rownames(topics.s),community=as.vector(membership(comm)),stringsAsFactors = F)
write.table(clust,paste(output,'cisTopics_singlet.louvain_clusters.txt',sep='.'),quote=F,sep='\t',row.names=F)

u <- cbind(umap.s,as.character(clust$community))
colnames(u)[3] <- 'community'

pdf(paste0(output,'.cisTopics_umap_singlet.louvain_clusters.pdf'),width=5,height=4)
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

rm(edge_list,g,ig)
gc()

# plot genes on UMAP
rna.smooth <- readRDS(rna_smooth_rds)
matchBC.s <- matchShareSamples( atac_bc = colnames(se.s),
                              rna_bc = colnames(rna.smooth),
                              isMulti = isMultiSample,
                              adList = adList
                              )
cat('Total ATAC barcodes:',ncol(se.s),'\n')
cat('Total RNA barcodes:',ncol(rna.smooth),'\n')
cat('Shared barcodes between ATAC and RNA:',length(matchBC.s[[1]]),
    paste0('(',100*round(length(matchBC.s[[1]])/ncol(se.s),3),'% of ATAC)\n'))

rna.sub <- rna.smooth[,matchBC.s[[2]]]
umap.sub <- umap.s[matchBC.s[[1]],]

pdf(paste0(output,'.umap_singlet.RNA.marker_genes.pdf'),width=4.5,height=4)
for (marker in markerList){
  plotUMIgene(umap.sub,rna.sub,marker,zCap,umap_pt_size,palette='brewer_heat')
  #cont <- readline(prompt = 'Enter')
}
dev.off()
