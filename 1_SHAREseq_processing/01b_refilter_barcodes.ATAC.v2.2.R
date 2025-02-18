output <- '' # prefix for output
# output from SHARE v2 workflow
bcFile <- '' #barcodes.txt
readFile <- '' #filtered.counts.csv.gz
scQCfile <- '' #counts.csv.gz

read_thresh <- 500 # 300 for mouse colitis tissue, 500 for mouse organoids and mouse adenoma tissue, 800 for human organoids
lib_thresh <- 500 # 500 for mouse colitis tissue, mouse adenoma tissue, mouse organoids and human organoids


source('../tools/barcode_filtering.v1.1.R')

bc <- read.csv(bcFile,header=F)
readCounts <- read.csv(readFile)
scQC <- read.csv(scQCfile)

bc$bcID <- with(bc,paste(V1,V2,V3,V4,sep=','))
readCounts$bcID <- with(readCounts,paste(R1,R2,R3,P5,sep=','))
readCounts <- readCounts[(readCounts$bcID %in% bc$bcID),]
scQC$bcID <- with(scQC,paste(R1,R2,R3,P5,sep=','))
scQC <- scQC[(scQC$bcID %in% bc$bcID),]

## filter by total fragments per barcode
pdf(paste0(output,'.frag_distribution.pdf'),width=5,height=4)
with(readCounts,hist(log10(fragments+1),main='Fragments per barcode'))
abline(v=log10(read_thresh),lty=2,col='red')
dev.off()
bc_pf1 <- readCounts$bcID[(readCounts$fragments >= read_thresh)]

## check barcode wells
wells <- wellReadPlot(bc_pf1,readCounts,3)
pdf(paste0(output,'.frag_well_distribution.pf1.pdf'),width=5,height=3)
plotShareWells(wells[[1]],1)
plotShareWells(wells[[2]],2)
plotShareWells(wells[[3]],3)
dev.off()

## filter by lib size
qc_pf1 <- scQC[(scQC$bcID %in% bc_pf1),]
colnames(qc_pf1)[grepl('_libsize',colnames(qc_pf1))] <- 'libsize'
pdf(paste0(output,'.libsize_pf1.pdf'),width=5,height=5)
lib <- qc_pf1[order(qc_pf1$libsize,decreasing=T),]
plot(log10(lib$libsize),pch=16,col='navy',xlab='Barcode rank',ylab='log10(Estimated library size)',main='Library Size PF1')
abline(h=log10(lib_thresh),lty=2,col='red')
avg <- floor(median(qc_pf1$libsize))
text(0.9*nrow(qc_pf1),0.9*max(log10(qc_pf1$libsize)),labels=paste0('Median lib size: ',avg),pos=2)
dev.off()
qc_pf2 <- qc_pf1[(qc_pf1$libsize > lib_thresh),]


writeLines(bc_pf1,paste0(output,'.bc_pf.txt'))

