output <- ''
methyl_file <- '' # file with methylation info per CpG (consolidated across strands). columns: chr, start, end, sample1_frac_methyl, sample2_frac_methyl ...
meta_file <- '' # metadata file where replicate column  matches sample order in methyl file

rep_order <- c() # order for plotting replicates in heatmap

ncores <- 4

library(dplyr)
library(stringi)
library(Matrix)matrixStats
library(ComplexHeatmap)
library(circlize)
library(pbmcapply)
library(data.table)

###################

getByPos <- function(interval,bed,cores = 4){
  split_pos <- stri_split_fixed(interval,':',simplify=T)
  range <- stri_split_fixed(split_pos[2],'-',simplify = T)
  df <- bed %>% filter(chr == split_pos[1]) %>% filter(start >= as.numeric(range[1])) %>% filter(end <= as.numeric(range[2]))
  
  dshift <- data.frame(startX = df$start - df$start[1],
                       endX = df$end - df$start[1])
  # extend each CpG to halfway to next CpG
  pad_list <- list()
  for (i in 1:nrow(dshift)){
    if (i == 1){
      padstart <- dshift$startX[i]
    } else {
      padstart <- (dshift$endX[i-1] + dshift$startX[i])/2
    }
    if (i == nrow(dshift)){
      padend <- dshift$endX[i]
    } else {
      padend <- (dshift$endX[i] + dshift$startX[i+1])/2
    }
    pad_list[[i]] <- data.frame(padStart = padstart,padEnd = padend)
  }
  padded <- cbind(dshift,do.call('rbind',pad_list))
  
  # make a position matrix
  pos_list <- pbmclapply(X=1:nrow(padded),FUN=function(i){
    n_bases <- floor(padded$padEnd[i]) - floor(padded$padStart[i])
    d_pos <- df[rep(i,n_bases),4:ncol(df)]
    rownames(d_pos) <- (floor(padded$padStart[i])+1):floor(padded$padEnd[i])
    return(d_pos)
  })
  mat <- do.call('rbind',pos_list)
  colnames(mat) <- colnames(df)[4:ncol(df)]
  
  return(t(mat))
}

plotByPos <- function(interval,bed,colors=NULL,breaks=NULL,rowOrder=NULL,cores=4){
  dmat <- getByPos(interval,bed,cores)
  if (!is.null(rowOrder)){
    dmat <- dmat[rowOrder,]
  }
  
  if (is.null(colors)){
    ht <- Heatmap(dmat,show_column_names = F,cluster_rows = F, cluster_columns = F,
                  column_title = interval)
  } else {
    col_fun <- colorRamp2(breaks = breaks, colors = colors)
    ht <- Heatmap(dmat,show_column_names = F,cluster_rows = F, cluster_columns = F,
                  column_title = interval,
                  col = col_fun,
                  heatmap_legend_param = list(title = "methylation",
                                              at = breaks,
                                              col = col_fun
                  ))
  }
  draw(ht)
}

###################

d <- fread(methyl_file)
meta <- read.table(meta_file,header=T,sep='\t')
d <- data.frame(d)
colnames(d) <- c('chr','start','end',meta$replicate)
head(d)
meta <- meta[match(rep_order,meta$replicate),]


pdf(paste0(output,'.EM_by_CpG.region_set1.pdf'),width=8,height=2)
# region 1
plotByPos(interval = "chr2:118115000-118118000", bed=d,
          rowOrder = meta$replicate,
          colors=c('white','gray50','black'),
          breaks=c(0,0.33,1)
)

# region 2
plotByPos(interval = "chr3:30019000-30022000", bed=d,
          rowOrder = meta$replicate,
          colors=c('white','gray75','gray50','gray25','black'),
          breaks=c(0,0.25,0.5,0.75,1)
)
dev.off()
