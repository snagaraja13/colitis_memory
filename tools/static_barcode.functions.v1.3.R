library(parallel)
library(pbmcapply)
library(igraph)
library(stringi)
library(Matrix)
library(dplyr)
library(ggplot2)
library(pheatmap)

checkSpacerSeq <- function(	bc,				# LARRY barcode sequence following search sequence
                            spacer_tab # matrix with start index of spacer and spacer bases 
                            ){	
  num_spacer <- floor(nchar(bc)/6)
  starts <- spacer_tab$ind[1:num_spacer]
  spacer_inds <- mapply(c,starts,starts+1,SIMPLIFY=F)
  check = T
  for (i in 1:length(spacer_inds)){
    inds = spacer_inds[[i]]
    if (as.character(spacer_tab$seq[i]) != substr(bc,inds[1],inds[2])){
      check = F
      break
    }
  }
  return(check)
}

validateBarcodeSeq <- function( bc,
                                min_length = 24,
                                max_length = 48,
                                spacerSeq = c('CT','AC','TC','GT','TG','CA','AT','GC'),
                                spacerPos = c(5,11,17,23,29,35,41,47),
				ncores = 4
){
  cat('-- Trimming excess bases\n')
  bc <- substr(bc,1,max_length)
  
  cat('-- Removing barcodes shorter than',min_length,'\n')
  bc_pf1 <- bc[(nchar(bc) >= min_length)]
  
  cat('-- Checking spacer sequence in reads...\n')
  spacers <- data.frame(ind=spacerPos,seq=spacerSeq)
  valid_bc <- unlist(pbmclapply(X=bc_pf1,FUN=function(x){checkSpacerSeq(x,spacers)},mc.cores=ncores))
  bc_pf2 <- bc_pf1[valid_bc]
  n_pf2 <- sum(valid_bc)
  cat('> Barcode reads with valid spacers:',n_pf2,'/',length(bc_pf1),paste('(',round(100*n_pf2/length(bc_pf1),1),'%)',sep=''),'\n')
  
  return(bc %in% bc_pf2)
}

# base prob density functions
pois <- function(lambda,k){return(exp(-lambda)*(lambda^k)/factorial(k))}
gaus <- function(mu,sigma,k){return(exp(-0.5*((k-mu)/sigma)^2)/(sigma*sqrt(2*pi)))}
mixmodel <- function(alpha1,alpha2,L,M,S,k){return(alpha1*pois(L,k)+alpha2*gaus(M,S,k))}
postProbGaus <- function(alpha1,alpha2,lambda,mu,sigma,k){
  return(gaus(mu,sigma,k)*alpha2/mixmodel(alpha1,alpha2,lambda,mu,sigma,k))
}

reEstimateParam <- function(old.wt1,old.wt2,old.pois.m,old.gaus.m,old.gaus.sd,data){
  #postProbGaus <- function(k){return(gaus(old.gaus.m,old.gaus.sd,k)*old.wt2/mixmodel(old.wt1,old.wt2,old.pois.m,old.gaus.m,old.gaus.sd,k))}
  pp <- data.frame(       x=data,
                          pGaus=unlist(lapply(X=data,FUN=function(x){
                            postProbGaus(old.wt1,old.wt2,old.pois.m,old.gaus.m,old.gaus.sd,x)
                          }))
  )
  pp$pPois <- 1-pp$pGaus
  new.wt1 <- sum(pp$pPois)/nrow(pp)
  new.wt2 <- sum(pp$pGaus)/nrow(pp)
  new.pois.m <- with(pp,sum(x*pPois)/sum(pPois))
  new.gaus.m <- with(pp,sum(x*pGaus)/sum(pGaus))
  new.gaus.sd <- sqrt(with(pp,sum(pGaus*(x-new.gaus.m)^2)/sum(pGaus)))
  log.likelihood <- sum(log(mixmodel(new.wt1,new.wt2,new.pois.m,new.gaus.m,new.gaus.sd,data)))
  return(c(new.wt1,new.wt2,new.pois.m,new.gaus.m,new.gaus.sd,log.likelihood))
}

fitCompModel <- function( vector, #numeric vector of values to model
                          N = 100, #number of iternations to fit model
                          dir.name='./', #directory to write files in plot = TRUE
                          log=TRUE, #model log transformed values of vector
                          plot=TRUE #write plots to file
                          ){
  if (log){
    d <- log2(vector)
  } else {
    d <- vector
  }
  
  ## estimate initial parameters by clustering data
  clust <- kmeans(d,2)$cluster
  mean1 <- mean(d[clust==1])
  sd1 <- sd(d[clust==1])
  mean2 <- mean(d[clust==2])
  sd2 <- sd(d[clust==2])
  
  if (mean1 < mean2){
    #clust1 is poisson
    pois.m <- mean1
    gaus.m <- mean2
    gaus.sd <- sd2
    wt1 <- sum(clust==1)/length(clust)
    wt2 <- sum(clust==2)/length(clust)
  } else {
    #clust1 is gaussian
    gaus.m <- mean1
    gaus.sd <- sd1
    pois.m <- mean2
    wt1 <- sum(clust==2)/length(clust)
    wt2 <- sum(clust==1)/length(clust)
  }
  log.likelihood <- sum(log(mixmodel(wt1,wt2,pois.m,gaus.m,gaus.sd,d)))
  
  ## Plot initial estimate
  xmax <- ceiling(max(d,na.rm=T))
  if (plot){
    cat('Plotting initial estimate of fit\n')
    hist(d,freq=F,xlim=c(1,ceiling(max(d,na.rm=T))))
    x <- seq(1,xmax,by=1)
    lines(x,mixmodel(wt1,wt2,pois.m,gaus.m,gaus.sd,x),col='red',lwd=2)
  }
  
  mle.table <- data.frame(matrix(data=0,ncol=7,nrow=(N+1)))
  colnames(mle.table) <- c('iter','pois.alpha','gaus.alpha',"pois.lambda",'gaus.mu','gaus.sd','log.likelihood')
  mle.table$iter <- 0:N
  mle.table[1,] <- c(0,wt1,wt2,pois.m,gaus.m,gaus.sd,log.likelihood)
  
  cat('Refitting model for',N,'iterations...\n')
  p <- c(wt1,wt2,pois.m,gaus.m,gaus.sd)
  for (i in 1:N){
    mle.table[(i+1),2:7] <- reEstimateParam(p[1],p[2],p[3],p[4],p[5],d)
    p <- as.numeric(mle.table[(i+1),2:6])
    if ((i %% 10) == 0){
      cat('-',i,'iterations complete...\n')
    }
  }
  
  gaus.prob <- postProbGaus(p[1],p[2],p[3],p[4],p[5],d)
  gaus.thresh <- min(d[gaus.prob > 0.5])
  
  if (plot){
    cat('Plotting convergence of solution and final fit\n')
    
    with(mle.table,
         plot(   iter,
                 log.likelihood,
                 pch=16,xlab='iteration',ylab='ln(likelihood)',main='Parameter estimation')
    )
    hist(   d,freq=F,
            ylab='Probability density',
            xlim=c(1,ceiling(max(d,na.rm=T)))
    )
    x <- seq(1,xmax,by=1)
    lines(x,mixmodel(p[1],p[2],p[3],p[4],p[5],x),col='gray',lwd=2)
    lines(x,p[1]*pois(p[3],x),col='red',lwd=1)
    lines(x,p[2]*gaus(p[4],p[5],x),col='blue',lwd=1)
    abline(v=gaus.thresh,lty=2,lwd=2,col='black')
    legend( 'topright',
            legend=c('Gaussian','Poisson','Mixture'),
            lty=1,col=c('blue','red','gray'),cex=0.75
    )
  }
  return(list(mle.table,p,gaus.thresh))
}


levBarcodeCollapse <- function(	bcList,		# vector of barcode sequences
                                bcCounts = NULL, # corresponding vector of counts per barcode
                                dist_max,	# maximum Levenshtein distance between barcodes to link
                                num_cores,
                                graph = F, #whether to use graph based collapsing
                                plot = F,
                                prompt = F
                                ){	# cores for parallelization
  if (graph){
    cat('Collapsing with graph-based clustering\n')
    unique_bc <- unique(bcList)
    dist_list <- mclapply(X=1:length(unique_bc),FUN=function(x){
      return(t(adist(unique_bc[x],unique_bc)))
    },mc.cores=num_cores)
    lv_dist <- do.call('cbind',dist_list)
    g <- graph.empty(nrow(lv_dist))
    for (i in 1:(nrow(lv_dist)-1)){
      sub <- lv_dist[i,(i+1):ncol(lv_dist)]
      comm <- sum(sub <= dist_max)
      if (length(comm) > 0){
        comm_shift <- comm + i
        edges <- vector()
        for (j in 1:length(comm_shift)){
          edge <- c(i,comm_shift[j])
          edges <- append(edges,edge)
        }
      }
      g <- add_edges(g,edges)
    }
    bc_clust <- cluster_louvain(as.undirected(g))
    bc_collapsed <- data.frame(bc=unique_bc,clone=as.vector(membership(bc_clust)))
    cat('-- Collapsed',length(unique_bc),'distinct barcode sequences into',length(unique(bc_collapsed$clone)),'barcode families\n')
    if (plot){
      plot.igraph(g,vertex.label=NA,vertex.size=0.2,edge.arrow.size=0,edge.width=1)
    }
  } else {
    cat('- Collapsing by sequence abundance\n')
    cat('-- Finding initial sequence similarity\n')
    d <- data.frame(bc=bcList,counts=bcCounts)
    dist_list <- pbmclapply(X=1:nrow(d),FUN=function(x){return(t(adist(d[x,1],d[,1])))},mc.cores=ncores)
    dist_mat <- do.call('cbind',dist_list)
    rownames(dist_mat) <- d$bc
    if(plot){
      pheatmap(dist_mat,main='Uncollapsed barcodes',show_rownames = F,show_colnames = F)
      if (prompt){
        cont <- readline(prompt = 'Enter to continue')
      }
    }
    
    # find initial match
    cat('-- Matching to most abundant sequences\n')
    d <- d[order(d$counts,decreasing=T),]
    match_list <- pbmclapply(X=d$bc,FUN=function(x){
      match_vect <- which(t(adist(x,d$bc)) < dist_max)
      if (length(match_vect) > 0){
        return(d$bc[min(match_vect)])
      } else {
        return(0)
      }
    })
    bc_collapsed1 <- data.frame(bc = d$bc, matched_seq = unlist(match_list))
    cat('-- Initially collapsed',length(bcList),'distinct barcode sequences into',length(unique(bc_collapsed1$matched_seq)),'barcode families\n')
    top_bc <- unique(bc_collapsed1$matched_seq)
    top_dist_list <- pbmclapply(X=top_bc,FUN=function(x){return(t(adist(x,top_bc)))},mc.cores=ncores)
    top_dist_mat <- do.call('cbind',top_dist_list)
    rownames(top_dist_mat) <- top_bc
    if(plot){
      pheatmap(top_dist_mat,main='Consensus barcodes - Round 1',show_rownames = F,show_colnames = F)
      if (prompt){
        cont <- readline(prompt = 'Enter to continue')
      }
    }
    
    # collapse consensus sequences if top bc < 2 from another
    cat('-- Creating final clone calls\n')
    d2 <- d[match(top_bc,d$bc),]
    match_list2 <- pbmclapply(X=1:nrow(d2),FUN=function(x){
      match_vect <- which(t(adist(d2$bc[x],d2$bc)) < dist_max)
      if (length(match_vect) > 0){
        return(min(match_vect))
      } else {
        return(0)
      }
    })
    bc_collapsed2 <- data.frame(bc = top_bc, clone = unlist(match_list2))
    
    top_bc2 <- d2$bc[unique(unlist(match_list2))]
    top_dist_list2 <- pbmclapply(X=top_bc2,FUN=function(x){return(t(adist(x,top_bc2)))},mc.cores=ncores)
    top_dist_mat2 <- do.call('cbind',top_dist_list2)
    rownames(top_dist_mat2) <- top_bc2
    if(plot){
      pheatmap(top_dist_mat2,main='Consensus barcodes - Round 2',show_rownames = F,show_colnames = F)
      if (prompt){
        cont <- readline(prompt = 'Enter to continue')
      }
    }
    bc_collapsed <- data.frame( bc = bc_collapsed1$bc, 
                                clone = bc_collapsed2$clone[match(bc_collapsed1$matched_seq,bc_collapsed2$bc)]
                                )
    cat('-- Final: Collapsed',length(bcList),'distinct barcode sequences into',length(unique(bc_collapsed$clone)),'barcode families\n')
  }
  return(bc_collapsed)
}

plotCloneSizes <- function( clone_mat, # data frame with clone names in first column and counts in second
                            N=1000,
                            log=T
                            ){
  if (nrow(clone_mat) > N){
    clone_mat <- clone_mat[1:N,]
    title <- paste('Top',N,'Clones',sep=' ')
  } else {
    title <- 'All Clones'
  }
  colnames(clone_mat)[1:2] <- c('clone','counts')
  clone_mat <- clone_mat[order(clone_mat$counts,decreasing=F),]
  clone_mat$clone <- factor(clone_mat$clone,levels = clone_mat$clone)
  if (log){
    if (min(clone_mat$counts,na.rm=T) > 0){
      g <- ggplot(clone_mat,aes(y=clone,x=log10(counts)))
    } else {
      g <- ggplot(clone_mat,aes(y=clone,x=log10(counts+1)))
    }
  } else {
    g <- ggplot(clone_mat,aes(y=clone,x=counts))
  }
  print(g + geom_bar(stat='identity',fill='blue') + ggtitle(title))
}
