output <- ''
by_clone_file <- '' # scores_by_clone.txt file from 04_SHARE_TRACE.variance_testing

features_to_plot <- c('Hnf4g','Fos','Hoxd10','Tead2','Hhex','Sox9','Ctcf')

n_plot <- 50
cap <- 5
nPerm <- 100
cores <- 4

myPath <- .libPaths()
myPath <- c(myPath,'/mnt/bin/R/library/')
.libPaths(myPath)
library(dplyr)
library(ggplot2)
library(pbmcapply)

z_by_clone <- read.table(by_clone_file,sep='\t',header=T)

# violin plots
pdf(paste0(output,'.top_violin_plots.pdf'),width=10,height=5)
for (feat in features_to_plot){
  d <- z_by_clone %>% filter(feature == feat)
  #make matrix for all cells
  d_all <- d
  d_all$clone <- 'all_cells'
  d_all$sample <- 'all_cells'
  
  # combine and cap values at +/-3
  d2 <- data.frame(rbind(d,d_all))
  d2$score[d2$score > cap] <- cap
  d2$score[d2$score < (-1*cap)] <- (-1*cap)
  
  # filter clones to most abundant
  clone_data <- d %>% group_by(clone) %>% summarise(n_cells = n(),
                                                    sd = sd(score,na.rm=T),
                                                    mean = mean(score,na.rm=T)
                                                    )
  top_clones <- clone_data %>% arrange(desc(n_cells)) %>% head(.,n_plot) %>% pull(clone)
  clone_order <- clone_data %>% filter(clone %in% top_clones) %>% arrange(desc(mean)) %>% pull(clone)
  
  d2$clone <- factor(d2$clone,levels = c('all_cells',clone_order))
  dsum <- (d2 %>% filter(clone %in% c('all_cells',clone_order)) %>% group_by(clone) %>% summarise(medScore = mean(score))) 
  dsum$clone <- factor(dsum$clone,levels=c('all_cells',clone_order))
  g1 <- ggplot((d2 %>% filter(clone %in% c('all_cells',clone_order))),aes(x=clone,y=score,fill=sample)) + geom_violin() +
    scale_fill_manual(values = c('all_cells' = 'blue', 'DSS' = 'red', 'Control' = 'gray75')) + 
    ggtitle(feat) + xlab('Clones ordered by average score') + ylab('Chromatin activation score') +
    geom_hline(yintercept = 0, linetype = 'dashed') + 
    theme(axis.text.x = element_text(angle=90, vjust=.5, hjust=1),
        panel.background = element_rect(fill = 'white', color = 'black')
        ) + 
    geom_point(data = dsum,aes(x=clone,y=medScore),size=1,fill='black')
  print(g1)
  #cont <- readline(prompt = 'Enter')
}
dev.off()

