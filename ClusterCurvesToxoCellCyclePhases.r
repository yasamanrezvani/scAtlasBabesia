library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(tidytext)
library(RColorBrewer)
library(parallel)
library(ComplexHeatmap)
library(circlize)

source('./util_funcs.R')

getCurvePeakLoc <- function(t, y){
  
  ## Fitting the estimated kernel with smooth splines
  spline.fit <- smooth.spline(x = t, y = y)
  
  ## Compute the derivatives of the fitted splines
  s.0 <- predict(spline.fit, spline.fit$x, deriv=0)
  s.1 <- predict(spline.fit, spline.fit$x, deriv=1)
  s.derv <- data.frame(s0=s.0$y, s1=s.1$y)
  
  ## Get the location of the extrema
  locs <- rle(den.sign <- sign(s.derv$s1))
  
  ## Maxima
  inc.ind <- which(locs$values == 1)
  if(length(inc.ind) > 1){
    maxima.ind = {}
    for(i in inc.ind){
      maxima.ind = c(maxima.ind, sum(locs$lengths[1:i]))
    }
    ## Interpolate a point between the location where derivative changes sign
    maxima = (spline.fit$x[maxima.ind] + spline.fit$x[(maxima.ind + 1)]) / 2
    maxima = maxima[!is.na(maxima)]
    ## Get the maximum values
    maxval = predict(spline.fit, maxima)
    
    ## Get the outliers
    maxima.outliers = which(maxval$y >= quantile(maxval$y, prob = 0.9))
    
    ## Peaks for entities of interest
    entity.x = maxval$x[maxima.outliers]
    entity.y = maxval$y[maxima.outliers]
  }else{
    entity.x <- spline.fit$x[which.max(spline.fit$y)]
    entity.y <- spline.fit$y[which.max(spline.fit$y)]
  }
  
  return(entity.x)
}

# 
# sc.tc.df.adj <- readRDS('../Input/compScBabesia/RData/all_sc_tc_df_adj.RData')
# sc.tc.df.adj <- sc.tc.df.adj[-5]
# # # 
# sc.tc.fits <- readRDS('../Input/compScBabesia/RData/all_sme_fits_sc_tc_20min.RData')
# sc.tc.fits <- sc.tc.fits[-5]

sc.tc.df.adj <- readRDS('~/Downloads/all_sc_tc_df_adj.RData')
sc.tc.fits <- readRDS('~/Downloads/all_sme_fits_sc_tc_20min.RData')

#markers.sig <- readRDS('../Input/compScBabesia/RData/all.markers.sig.RData')
markers.sig <- readRDS('../Input/compScBabesia/RData/all.markers.sig.cell.cycle.phase.RData')

sc.tc.mus <- lapply(1:length(sc.tc.fits), function(i){
  sc.tc.mu <- smoothSplineSmeFits(sc.tc.fits[[i]], unique(sc.tc.df.adj[[i]]$variable), extend = F)
  colnames(sc.tc.mu) <- c('GeneID', 't', 'y')
  sc.tc.mu <- sc.tc.mu %>% dplyr::filter(GeneID %in% markers.sig[[i]]$markers.sig$GeneID) %>%
    pivot_wider(names_from = 'GeneID', values_from = 'y') %>%
    as.data.frame()
  
  return(sc.tc.mu)
})

names(sc.tc.mus) <- names(sc.tc.fits)


# sync.tc.df <- readRDS('../Input/compScBabesia/RData/bd_sync_tc_df.RData')
# sync.tc.fit <- readRDS('../Input/compScBabesia/RData/bd_sme_fits_sync_tc_20min.RData')

sync.tc.df <- readRDS('~/Downloads/bd_sync_tc_df.RData')
sync.tc.fit <- readRDS('~/Downloads/bd_sme_fits_sync_tc_20min.RData')


sync.tc.mu <- smoothSplineSmeFits(sync.tc.fit, unique(sync.tc.df$variable), extend = F)
colnames(sync.tc.mu) <- c('GeneID', 't', 'y')
sync.tc.mu <- sync.tc.mu %>% dplyr::filter(GeneID %in% markers.sig[[4]]$markers.sig$GeneID) %>%
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>%
  as.data.frame()


## Generate the clusters
num.clust <- 5L

sc.hc_dtws <- lapply(sc.tc.mus, function(x){
  sc.hc_dtw <- dtwClustCurves(x[,2:ncol(x)], nclust = num.clust)
  return(sc.hc_dtw)
})


# toxo top marker info 

toxo.markers <- readRDS('../Input/compScBabesia/RData/TG.markers.sig.RData')
toxo.bdiv.orth <- read.xlsx('../Input/compScBabesia/Orthologs/GT1_BDiv.xlsx')


toxo.markers <- inner_join(toxo.markers, toxo.bdiv.orth, by = c('GeneID' = 'GT1'))
top.toxo.markers <- toxo.markers %>% dplyr::filter(avg_log2FC > 1) %>% 
  group_by(cluster) %>% slice_min(order_by = p_val_adj, n = 20)

# ## Extract the top markers in Each babesia spp. Filter for ones that appear in Toxo markers

top.markers <- lapply(1:length(markers.sig), function(i){
  tmp <- markers.sig[[i]]$markers.sig %>% dplyr::filter(avg_log2FC > 1) %>% 
    group_by(cluster) %>% slice_min(order_by = p_val_adj, n = 100) %>%
    ungroup() %>% transmute(GeneID = GeneID, spp = names(markers.sig)[i])
  tmp <- inner_join(tmp, top.toxo.markers, by = c('GeneID' = 'BDiv'))
  tmp <- tmp %>% transmute(GeneID = GeneID, spp = spp, Phase = cluster, ProductDescription = ProductDescription)
})
names(top.markers) <- names(markers.sig)

phase.boundaries <- lapply(1:length(sc.tc.mus), function(i){
  sc.tc.mu.scale <- sc.tc.mus[[i]]
  sc.tc.mu.scale[,2:ncol(sc.tc.mu.scale)] <- scale(sc.tc.mu.scale[,2:ncol(sc.tc.mu.scale)],
                                                   center = T,scale = T)
  sc.tc.mu.scale <- sc.tc.mu.scale %>%  as.data.frame() %>% 
    pivot_longer(cols = -t, names_to = 'GeneID', values_to = 'y')
  
  sc.peak.order <- sc.tc.mu.scale %>% group_by(GeneID) %>% summarise(peak.ord = getCurvePeakLoc(t, y))
  sc.tc.mu.scale <- left_join(sc.tc.mu.scale, sc.peak.order, by = 'GeneID')
  
  sc.tc.mu.scale <- inner_join(sc.tc.mu.scale, top.markers[[i]], by = 'GeneID')
  
  tmp <- sc.tc.mu.scale %>% ungroup() %>% group_by(Phase) %>% mutate(mean.peak = mean(peak.ord), sd.peak = sd(peak.ord)) %>%
    mutate(lb = mean.peak - sd.peak, ub = mean.peak + sd.peak)
  tmp$spp <- names(sc.tc.fits)[i]
  
  return(tmp)
  
})
names(phase.boundaries) <- names(sc.tc.fits)

phase.boundaries <- do.call('rbind', phase.boundaries)

phase.boundaries <- phase.boundaries %>% dplyr::filter(Phase != 'G1.a')
phase.boundaries$Phase <- as.character(phase.boundaries$Phase)
phase.boundaries$Phase[phase.boundaries$Phase == 'G1.b'] <- 'G1'
phase.boundaries <- phase.boundaries %>% dplyr::filter(spp != 'bmic')
phase.boundaries$Phase <- factor(phase.boundaries$Phase, levels = c('G1', 'S', 'M', 'C'))

p <- ggplot(phase.boundaries, aes(x=mean.peak, color = Phase)) +
  geom_density(size = 1) + theme_bw() + 
  geom_vline(data=phase.boundaries, aes(xintercept=mean.peak, color=Phase),
             linetype="dashed", size = 0.8) +
  scale_x_continuous(breaks = round(seq(min(phase.boundaries$t), max(phase.boundaries$t), by = 0.5),1)) + 
  facet_grid(spp~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE)) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )+ expand_limits(x=c(0,12))
plot(p)

ggsave(filename="../Output/compScBabsSpecies/figs/phase_boundaries.png",
       plot=p,
       width = 7, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

# add transition boundaries line to data and figure 
phase.boundaries.list <- split(phase.boundaries, phase.boundaries$spp)

# manual bounds
B.big.ph.bnd <- c(0, 2.5, 6.75, 8.25, 11.35, 12)
B.bov.ph.bond <- c(0, 2.25, 7.5 , 8.25, 11.35, 12)
B.div.cow.ph.bond <- c(0, 2.75, 7.15, 8.60, 11.35, 12)
B.div.hum.ph.bond <- c(0, 2.25, 6.75, 8.60, 11.35 ,12)

ph.trans <- list(B.big.ph.bnd, B.bov.ph.bond, B.div.cow.ph.bond, B.div.hum.ph.bond)
names(ph.trans) <- names(phase.boundaries.list)

cell.cycle.boundaries <- lapply(1:length(ph.trans), function(i){
  df <- phase.boundaries.list[[i]]
  df <- phase.boundaries.list[[i]] %>% 
    mutate(Phase.transition = ifelse(Phase == "G1", ph.trans[[i]][2], 
                                     ifelse(Phase == "S", ph.trans[[i]][3], 
                                            ifelse(Phase == "M", ph.trans[[i]][4], ph.trans[[i]][5]))))
})


names(cell.cycle.boundaries) <-  names(ph.trans)

cell.cycle.boundaries <-  do.call("rbind",cell.cycle.boundaries)

p <- ggplot(cell.cycle.boundaries, aes(x=mean.peak, color = Phase)) +
  geom_density(size = 1) + theme_bw() + 
  geom_vline(data=cell.cycle.boundaries, aes(xintercept=Phase.transition , 
                                             #color=Phase
  ),
  linetype="dashed", size = 0.8) +
  scale_x_continuous(breaks = round(seq(min(phase.boundaries$t), max(phase.boundaries$t), by = 0.5),1)) + 
  facet_grid(spp~., scales = "free",  space='free',
             labeller=label_wrap_gen(multi_line = TRUE)) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )+
  expand_limits(x=c(0,12))
plot(p)

ggsave(filename="../Output/compScBabsSpecies/figs/phase_boundaries_all_spp.png",
       plot=p,
       width = 7, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



## phase boundries figure 3

SM: "#c44237"

MC: "#ad8842"

G: "#1b5878"

C: "#e99093"

saveRDS(phase.boundaries, "../Input/compScBabesia/RData/phase.boundaries.RData")

phase.boundries.Bdiv_human <- phase.boundaries %>% filter(spp == 'bdiv_human')
p.bdiv <- ggplot(phase.boundries.Bdiv_human, aes(x=mean.peak, color = Phase)) +
  geom_density(aes(color = Phase), size = 1.5) + 
  scale_color_manual(values=c("G1" = "#1b5878", "S" = "#c44237", "M" = "#ad8842","C" = "#e99093"))+
  theme_bw() +
  geom_vline(xintercept=B.div.hum.ph.bond[-c(1,6)] ,linetype="dashed", size = 0.8) +
  scale_x_continuous(breaks = round(seq(min(phase.boundaries$t), max(phase.boundaries$t), by = 0.5),1)) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold"), 
    legend.title = element_text(size = "14", face = "bold"), 
    axis.text = element_text(size = 14)
  ) +
  expand_limits(x=c(0,12))

plot(p.bdiv)

ggsave(filename="../Output/compScBabsSpecies/figs/Bdiv_transition_boundaries.png",
       plot=p.bdiv,
       width = 12, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)




# # mid point of peaks in each phase
# mid.ph.trans <- lapply(phase.boundaries.list, function(ph){
#   
#   ph.mean.peak <- sort(unique(ph$mean.peak))
#   mids <- c()
#   for(i in 2:length(ph.mean.peak)){
#     mids[i-1] <- (ph.mean.peak[i-1] + ph.mean.peak[i]) / 2
#   }
#   
#   mids.phases <- c(mids)
#   
# })
# 
# mid.ph.trans


## transition phases identified by vision

B.big.ph.bnd <- c(0, 2.5, 6.75, 8.25, 11.35, 12)
B.bov.ph.bond <- c(0, 2.25, 7.5 , 8.25, 11.35, 12)
B.div.cow.ph.bond <- c(0, 2.75, 7.30, 8.75, 11.35, 12)
B.div.hum.ph.bond <- c(0, 2.25, 6.75, 8.75, 11.35 ,12)

ph.trans <- list(B.big.ph.bnd, B.bov.ph.bond, B.div.cow.ph.bond, B.div.hum.ph.bond)
names(ph.trans) <- names(sc.tc.fits)

#  reorder and scale genes for heatmap

sc.tc.mus.scale <- lapply(1:length(sc.tc.mus), function(i){
  
  sc.tc.mus.scale <- sc.tc.mus[[i]]
  
  sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)] <- scale(sc.tc.mus.scale[,2:ncol(sc.tc.mus.scale)],
                                                     center = T,scale = T)
  sc.tc.mus.scale <- sc.tc.mus.scale %>%  as.data.frame() %>% 
    pivot_longer(cols = -t, names_to = 'GeneID', values_to = 'y')
  
  
  ## Add curve cluster info
  sc.hc_dtw.df <- data.frame(GeneID = unique(sc.tc.mus.scale$GeneID), 
                             order = as.numeric(sc.hc_dtws[[i]]$order),
                             cluster = cutree(sc.hc_dtws[[i]],k = num.clust))
  
  sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.hc_dtw.df, by = 'GeneID')
  
  ## Reorder the genes within each cluster.
  sc.peak.order <- sc.tc.mus.scale %>% group_by(GeneID) %>% summarise(peak.ord = getCurvePeakLoc(t, y))
  sc.tc.mus.scale <- left_join(sc.tc.mus.scale, sc.peak.order, by = 'GeneID')
                                                             
  

  tmp <- sc.tc.mus.scale %>% ungroup() %>% group_by(cluster) %>% mutate(mean.peak = mean(peak.ord))
  tmp <- data.frame(cluster = unique(tmp$cluster[order(tmp$mean.peak)]), cluster.order = 1:length(unique(tmp$cluster)))
  tmp <- sc.tc.mus.scale %>% ungroup() %>% group_by(cluster) %>% mutate(mean.peak = mean(peak.ord)) %>%
    transmute(mean.peak = mean.peak) %>% distinct()
  tmp2 <- data.frame(cluster = unique(tmp$cluster[order(tmp$mean.peak)]), cluster.order = 1:length(unique(tmp$cluster)))
  tmp <- left_join(tmp, tmp2, by = 'cluster')
  
  
  sc.tc.mus.scale <- left_join(sc.tc.mus.scale, tmp, by = 'cluster')
  sc.tc.mus.scale$stage <- paste('dtw', sc.tc.mus.scale$cluster.order, sep = '')
  

  return(sc.tc.mus.scale)
})
  
names(sc.tc.mus.scale) <- names(sc.tc.fits)


# add  toxo phases 
titles <- c("B. big", "B. bov", "B. div (cow)", "B. div (human)")

sc.tc.mus.scale.phase <- lapply(1:length(titles), function(i){
  
  sc.tc.mus.scale <- sc.tc.mus.scale[[i]]
  ph <- ph.trans[[i]]
  
  # sc.tc.mus.scale <- sc.tc.mus.scale %>%
  #   dplyr::mutate(cell.cycle.phase = ifelse((t >=  ph[2] & t < ph[3]), 'SM',
  #                                           ifelse((t >= ph[3] & t < ph[4]), 'MC',
  #                                                  ifelse(t >= ph[4] & t < ph[5] , 'C', 'G'))))

  sc.tc.mus.scale <- sc.tc.mus.scale %>%
    dplyr::mutate(cell.cycle.phase = ifelse((t >=  ph[1] & t < ph[2]), 'G1.L',
                                     ifelse((t >= ph[2] & t < ph[3]), 'SM',
                                            ifelse(t >= ph[3] & t < ph[4] , 'MC',
                                                   ifelse(t >= ph[4] & t < ph[5], 'C', 'G1.E')))))

  sc.tc.mus.scale <- sc.tc.mus.scale %>%
    mutate(trans.time = ifelse((t >=  ph[1] & t < ph[2]), ph[2],
                                ifelse((t >= ph[2] & t < ph[3]), ph[3],
                                       ifelse(t >= ph[3] & t < ph[4] , ph[4],
                                              ifelse(t >= ph[4] & t < ph[5], ph[5], ph[6])))))
  # 
  
  sc.tc.mus.scale$cell.cycle.phase <- factor(sc.tc.mus.scale$cell.cycle.phase, levels = c('G1.L',  'SM', 'MC', 'C','G1.E'))
  #sc.tc.mus.scale$stage   <- gsub("^\\D+", "", sc.tc.mus.scale$stage)
  
  return(sc.tc.mus.scale)
})

names(sc.tc.mus.scale.phase) <- names(sc.tc.mus.scale)

## add Marker phase info 

sc.tc.mus.scale.df <- lapply(1:length(sc.tc.mus.scale.phase), function(i){
  
  markers <- markers.sig[[i]]$markers.sig
  names(markers) <- gsub('cluster', 'marker.phase', names(markers))
  markers <- markers %>% group_by(marker.phase) %>% mutate(num.degs  = n())
  
  sc.tc.mus.tab <- sc.tc.mus.scale.phase[[i]]

  sc.tc.mus.scale.markers <- left_join(sc.tc.mus.tab, markers, by = 'GeneID')
  sc.tc.mus.scale.markers$cell.cycle.phase <-  factor( sc.tc.mus.scale.markers$cell.cycle.phase ,
                                                      levels = c('G1.L', 'SM', 'MC', 'C', 'G1.E'))
  sc.tc.mus.scale.markers$marker.phase <- factor( sc.tc.mus.scale.markers$marker.phase, 
                                                 levels = c('G', 'SM', 'MC', 'C'))
  
  sc.tc.mus.scale.markers$num.deg.phase <- paste(sc.tc.mus.scale.markers$num.degs, sc.tc.mus.scale.markers$marker.phase, sep = " ")
  
  sc.tc.mus.scale.markers <- sc.tc.mus.scale.markers %>%
    mutate(num.deg.phase = factor(num.deg.phase, unique(arrange(sc.tc.mus.scale.markers, marker.phase)$num.deg.phase)))

  # sc.tc.mus.scale.markers <- sc.tc.mus.scale.markers %>% 
  #   mutate(num.degs = factor(num.degs, unique(arrange(sc.tc.mus.scale.markers, marker.phase)$num.degs)))
  # 
  
  return(sc.tc.mus.scale.markers)
  
})

names(sc.tc.mus.scale.df) <- names(sc.tc.mus.scale)



saveRDS(sc.tc.mus.scale.df,'../Input/compScBabesia/RData/all_sc_tc_mus_scale_toxo_inferred_cell_cycle_phases_Marker_phases_Progression_heatmap.RData')
sc.tc.mus.scale.df <- readRDS('../Input/compScBabesia/RData/all_sc_tc_mus_scale_toxo_inferred_cell_cycle_phases_Marker_phases_Progression_heatmap.RData')

#sc.tc.mus.scale.df <- readRDS('../Input/compScBabesia/RData/all_sc_tc_mus_scale_toxo_inferred_cell_cycle_phases_Marker_phases_Progression_heatmap.RData')
titles <- c("B. big", "B. bov", "B. div (cow)", "B. div (human)")
title.phases <- c("G1.L", "SM", "MC", "C", "G1.E")
names(title.phases) <- c('G1.L',  'SM', 'MC', 'C','G1.E')



ps <- lapply(1:length(titles), function(i){
  
  p <- ggplot(sc.tc.mus.scale.df[[i]], aes(x = t, reorder_within(GeneID, -peak.ord, stage ), fill = y)) +
    
    geom_tile() +
    scale_x_discrete(expand=c(0,0)) +
    
    ylab("Genes") + xlab("time/cells") +
    scale_fill_gradientn(colours = viridis::inferno(10)) +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 12, face = "bold"),
      #axis.title.y = element_text(size = 12, face = "bold"),
      axis.title.y = element_blank(),
      legend.position = "none") +
    facet_grid(num.deg.phase ~ cell.cycle.phase  , scales = 'free', space = 'free',
               labeller = labeller(cell.cycle.phase = title.phases)) +
    theme(strip.background=element_rect(fill='white', color = 'black'))+
    theme(panel.spacing = unit(0.1, "lines")) +
    # theme(strip.text.x=element_text(angle=0, hjust=0.5, size = 6,vjust=0.5, face = 'bold'), 
    #       strip.text.y = element_text( hjust = 0.5, size = 6, vjust = 0.5, face = 'bold'))+
    theme(strip.text = element_text(size = 14, face = 'bold'))+
    ggtitle(titles[i])+
    theme(
      plot.title = element_text(size=18, face = "bold.italic", color = 'black', hjust = 0.5),
      axis.title.x = element_text(size=14, face="bold"),
      #axis.title.y = element_text(size=14, face="bold")
    )
  
})
p <- grid.arrange(ps[[1]], ps[[2]], ps[[3]], ps[[4]], ncol = 4)


ggsave(filename="../Output/compScBabsSpecies/figs/gene_cluster_curves_cinferred_cell_cycle_phase_marker_phase_num_degs.png",
       plot=p,
       width = 32, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

# Bdiv_human 
ggsave(filename="../Output/compScBabsSpecies/figs/Bdiv_human_gene_cluster_curves_cinferred_cell_cycle_phase_marker_phase_num_degs.png",
       plot=ps[[4]],
       width = 8, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

## comparative  needs  work 
# ggplot 
sc.tc.mus.scale.df <- readRDS('../Input/compScBabesia/RData/all_sc_tc_mus_scale_toxo_inferred_cell_cycle_phases_Marker_phases_Progression_heatmap.RData')

spps <- c('bbig', 'bbov', 'bdiv_c', 'bdiv_h')
for(i in 1:length(titles)){
  sc.tc.mus.scale.df[[i]]$spp <- spps[i]
}

com.genes  <- map(sc.tc.mus.scale.df, ~.$GeneID) %>% purrr::reduce(intersect) # common cyclic genes across spp

#cons.markers.phases.sig.cc <- readRDS("../Input/compScBabesia/RData/conserved_markers_cross_spp_phase_based.RData")



all.sc.tc.mus.scale <- do.call("rbind", sc.tc.mus.scale.df)
all.sc.tc.mus.scale.com <- all.sc.tc.mus.scale %>% filter(GeneID %in% com.genes)

# get  the order of Bdiv human to  order the genes 
Bdiv_hum <- all.sc.tc.mus.scale %>% filter(spp == "bdiv_h") %>% select(GeneID, peak.ord) %>% distinct()
names(Bdiv_hum) <- gsub("peak.ord", "bdiv_peak.ord", names(Bdiv_hum)) 

all.sc.tc.mus.scale.com  <- inner_join(all.sc.tc.mus.scale.com, Bdiv_hum, by  = "GeneID")
all.sc.tc.mus.scale.com$aug.time <- all.sc.tc.mus.scale.com$t

for(i in 1:length(spps)){
  
  ind <- all.sc.tc.mus.scale.com$spp == spps[i]
  
  all.sc.tc.mus.scale.com$aug.time[ind] <- all.sc.tc.mus.scale.com$aug.time[ind] + (12 + 1/3) * (i - 1)
  
}

p <- ggplot(all.sc.tc.mus.scale.com, aes(x = aug.time, reorder(GeneID, -bdiv_peak.ord),  fill = y)) +
  geom_tile() +
  scale_x_discrete(expand=c(0,0)) +
  ylab("Genes") + xlab("time/cells") +
  facet_grid(. ~ spp , scales = 'free', space = 'free')+
  scale_fill_gradientn(colours = viridis::inferno(10)) +
  theme(
    axis.text.x = element_blank(),
    #axis.text.y = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.position = "none")+
  theme(strip.background=element_rect(fill='white', color = 'black'))+
  theme(panel.spacing = unit(1.5, "lines")) +
  theme(strip.text.x=element_text(angle=0, hjust=0.5,vjust=0.5, size = "14",face = 'bold'))+
  theme(
    plot.title = element_text(size=18, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold")
  )

p

ggsave(filename="../Output/compScBabsSpecies/figs/gene_cluster_curves_common_genes_cross_spp.png",
       plot=p,
       width = 28, height = 3,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

##########


## add toxo marker phases to S.O objects

addCellCyclePhaseManual <- function(S.O, phase.trans){
  
  ph <- phase.trans
  
  S.O@meta.data <- S.O@meta.data %>%
    mutate(cell.cycle.phase = ifelse((adj.time >=  ph[2] & adj.time < ph[3]), 'SM',
                                     ifelse((adj.time >= ph[3] & adj.time < ph[4]), 'MC',
                                            ifelse(adj.time >= ph[4] & adj.time < ph[5] , 'C', "G"))))
  
  
  S.O@meta.data <- S.O@meta.data %>%
    mutate(trans.time = ifelse((adj.time >=  ph[2] & adj.time < ph[3]), ph[3],
                               ifelse((adj.time >= ph[3] & adj.time < ph[4]), ph[4],
                                      ifelse(adj.time >= ph[4] & adj.time < ph[5] , ph[5], ph[6]))))
  return(S.O)
}



S.O.orth <- readRDS('../Input/compScBabesia/RData/S.O.integrated.list.pstime.GAM.indiv.Rdata')
S.O.orth <- S.O.orth[-5]


S.O.updated <- lapply(names(S.O.orth), function(nn){
  addCellCyclePhaseManual(S.O.orth[[nn]], ph.trans[[nn]])
})
names(S.O.updated) <- names(S.O.orth)

saveRDS(S.O.updated, "../Input/compScBabesia/RData/S.O.integrated.list.pstime.GAM.indiv.cell.cycle.phase.Rdata")

for(i in  1:length(S.O.updated)){
  Idents(S.O.updated[[i]]) <- "cell.cycle.phase"
}

p.list <- list()
p.list <- lapply(S.O.updated, function(S){
  
  p.list <- DimPlot(S, reduction = "pca", 
                    #group.by = "cells", 
                    split.by = 'spp',
                    pt.size = 1,
                    shape.by='spp',
                    label = TRUE, label.size = 6) + NoLegend() + 
    theme(panel.spacing = unit(0.5, "lines")) + 
    theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
    theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
    theme(
      axis.title.x = element_text(size=14, face="bold"),
      axis.title.y = element_text(size=14, face="bold")
    )
  
})

pp <- grid.arrange(p.list[[1]], p.list[[2]], p.list[[3]], p.list[[4]], ncol = 4)

ggsave(filename="../Output/compScBabsSpecies/figs/integrated_pca_toxo_phases_updated.png",
       plot=pp,
       width = 12, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



# Bdiv pca plot

Idents(S.O.updated$bdiv_human) <- "cell.cycle.phase"
S.O.updated$bdiv_human$cell.cycle.phase <- factor(S.O.updated$bdiv_human$cell.cycle.phase, 
                                                  levels = c('G', 'SM', 'MC', 'C'))


SM: "#c44237"

MC: "#ad8842"

G: "#1b5878"

C: "#e99093"

p <- DimPlot(S.O.updated$bdiv_human, reduction = "pca", 
             #group.by = "cells", 
             cols = c("#1b5878", "#c44237", "#ad8842","#e99093" ),
             split.by = 'spp',
             pt.size = 1,
             shape.by='spp',
             label = TRUE, label.size = 8)  + 
  NoLegend()+
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

plot(p)

ggsave("../Output/compScBabsSpecies/figs/Bdiv_human_pca.png", 
       plot = p, 
       width = 6, height = 6, units = "in",
       dpi = 300)
