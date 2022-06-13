library(Seurat)
library(gam)
library(princurve)
library(parallel)
library(tidyverse)
library(MyEllipsefit)
library(sctransform)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidytext)
library(RColorBrewer)
library(parallel)
library(ComplexHeatmap)
library(circlize)

source('./util_funcs.R')

prod.desc <- read.csv('../Input/compScBabesia/genes/BDiv_Prod_desc.csv')
prod.desc <- prod.desc %>% transmute(GeneID = Gene.ID, ProductDescriptionBdiv = Product.Description) %>% distinct()
toxo.bdiv.orth <- read.xlsx('../Input/compScBabesia/Orthologs/GT1_BDiv.xlsx')
colnames(toxo.bdiv.orth) <- c('GeneID_Toxo', 'GeneID', 'ProductDescriptionToxo')

num.cores <- detectCores()

makeMatchedContrasts <- function(S.O.integrated){
  
  objs <- as.character(unique(S.O.integrated@meta.data$phase.spp))
  
  contrasts <- data.frame(ref = objs, dummy = 1) %>% full_join( data.frame(query = objs, dummy = 1), by = 'dummy') %>% 
    mutate(ref.spp = gsub(":.*", "", ref), ref.phase = gsub(".*:", "", ref), 
           query.spp = gsub(":.*", "", query), query.phase = gsub(".*:", "", query))
  my.contrasts <- contrasts %>% dplyr::filter(ref.phase == query.phase & ref.spp != query.spp)
  
  return(my.contrasts)
  
}


getCellCyclePhaseMarkers <- function(S.O.integrated){
  S.O.integrated.merge.list <- SplitObject(S.O.integrated, split.by = "spp")
  spps <- names(S.O.integrated.merge.list)
  all.markers.list <- mclapply(1:length(spps), function(i){
    S.O <- S.O.integrated.merge.list[[i]]
    DefaultAssay(S.O) <- "RNA"
    Idents(S.O) <- 'cell.cycle.phase'
    markers <- FindAllMarkers(object = S.O, only.pos = TRUE)
    markers$GeneID = gsub('-', '_', markers$gene)
    colnames(markers)[colnames(markers) == 'cluster'] <- 'phase'
    markers$spp <- spps[i]
    return(markers)
  }, mc.cores = num.cores) 
  all.markers <- bind_rows(all.markers.list)
  
  all.markers.sig <- all.markers %>% dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.01)
  
  return(all.markers.sig)
}

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

exprExprData <- function(S.O, genes){
  expr <- data.frame(as.matrix(S.O@assays$RNA@data))
  expr$GeneID <- gsub('-', '_', rownames(as.matrix(S.O@assays$RNA@data)))
  expr <- expr %>% pivot_longer(-GeneID, names_to = 'Sample', values_to = 'expr')
  pc <- S.O@reductions$pca@cell.embeddings
  pc <- data.frame(pc) %>% dplyr::mutate(Sample = rownames(pc)) %>% 
    transmute(Sample = Sample, PC_1 = PC_1, PC_2 = PC_2, PC_3 = PC_3)
  umap <- S.O@reductions$umap@cell.embeddings
  umap <- data.frame(umap) %>% dplyr::mutate(Sample = rownames(umap)) %>% 
    transmute(Sample = Sample, UMAP_1 = UMAP_1, UMAP_2 = UMAP_2)
  
  meta.data <- data.frame(Sample = rownames(S.O@meta.data), Phase = S.O@meta.data$cell.cycle.phase)
  meta.data <- left_join(meta.data,
                         pc, by = 'Sample')
  meta.data <- left_join(meta.data, umap, by = 'Sample')
  
  expr.filt <- expr %>% dplyr::filter(GeneID %in% genes)
  
  expr.filt <- inner_join(meta.data, expr.filt, by = 'Sample')
  return(expr.filt)  
}


###cell cycle marker analysis done per spp (no common genes)
markers.sig <- readRDS('../Input/compScBabesia/RData/all.markers.sig.cell.cycle.phase.RData')
clusters <-  c('G', 'SM', 'MC', 'C')
spps <- c('Bbig', 'Bbov', 'Bdiv_cow', 'Bdiv_human')
titles <- c("B. big", "B. bov", "B. div (cow)", "B. div (human)")

for(i in 1:length(titles)){
  markers.sig[[i]]$markers.sig$spp <- spps[i]
}

all.markers.sig <- do.call(Map, c(f = rbind, markers.sig))
all.markers.sig <- all.markers.sig$markers.sig
all.markers.sig$spp <-  factor(all.markers.sig$spp, levels = c('Bbig', 'Bbov', 'Bdiv_cow', 'Bdiv_human'))
names(all.markers.sig) <- gsub("cluster", "phase", names(all.markers.sig))
all.markers.sig$phase <- factor(all.markers.sig$phase, levels = c('G', 'SM', 'MC', 'C'))
all.markers.sig.stat <- all.markers.sig %>% group_by(spp, phase) %>% summarise(num.deg = n()) 

####  bar plot figure

p <- ggplot(all.markers.sig.stat, aes(x=phase, y=num.deg, fill = phase)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=num.deg), vjust=1.5, color="black", size=8, fontface = 'bold')+
  facet_grid(. ~ spp, scales = "free", space='free',
             labeller=label_wrap_gen(multi_line = TRUE))+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    #axis.text.y = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"))+
  theme(strip.background=element_rect(fill='white', color = 'black'),
        panel.spacing = unit(1.5, "lines"), 
        strip.text.x=element_text(angle=0, hjust=0.5,vjust=0.5, size = "14",face = 'bold'),
        plot.title = element_text(size=18, face = "bold.italic", color = 'black'),
        axis.title.x = element_text(size=16, face="bold"),
        axis.title.y = element_text(size=16, face="bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(face = "bold", size = 12,  angle = 0), 
        strip.placement = "outside") +
  theme(legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 14))+
  theme(panel.spacing = unit(1.5, "lines")) +
  theme(legend.position = "none")

p


ggsave(filename="../Output/compScBabsSpecies/figs/cell_cycle_phase_per_spp_bar_plot.png",
       plot=p,
       width = 20, height = 3,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

out.dir <- "../Output/compScBabsSpecies/GoTerm/"
spp.markers <- names(markers.sig)

for(i in 1:length(names(markers.sig))){
  file.name <- names(markers.sig)[i]
  cat(paste('processing file', file.name))
  cat('\n')
  
  df <- markers.sig[[i]]$markers.sig
  file.name.xlsx <- paste(file.name, "xlsx", sep =".")
  write.xlsx(df, paste(out.dir, file.name.xlsx, sep = ""))


  }

##############################

##### Doing diff exp using anchored datasets
# Orthologs S.Os after  integration
S.O.integrated <- readRDS('../Input/compScBabesia/RData/S.O.integrated.list.pstime.GAM.indiv.cell.cycle.phase.Rdata')
spps <- names(S.O.integrated)
## Differential expression must be done on genes present in all datasets
comm.genes <- lapply(S.O.integrated, function(S.O){
  genes <- rownames(S.O@assays$RNA@data)
})

comm.genes <- Reduce(intersect, comm.genes)

S.O.integrated <- lapply(S.O.integrated, function(S.O){
  S.O <- subset(S.O, features = comm.genes)
})
names(S.O.integrated) <- spps

## We nned to re-merge and integrate the data for cross spp marker analysis
S.O.integrated.merge <- merge(S.O.integrated[[1]], S.O.integrated[2:4], 
                              add.cell.ids = names(S.O.integrated))
S.O.integrated.merge.list <- SplitObject(S.O.integrated.merge, split.by = "spp")
S.O.integrated.merge.list <- lapply(X = S.O.integrated.merge.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
features <- SelectIntegrationFeatures(object.list = S.O.integrated.merge.list)
anchors <- FindIntegrationAnchors(object.list = S.O.integrated.merge.list, 
                                  anchor.features = features, reference = 4)

S.O.integrated <- IntegrateData(anchorset = anchors)
# switch to integrated assay for PCA. Make sure to set to RNA for Differential Expression
DefaultAssay(S.O.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
S.O.integrated <- ScaleData(S.O.integrated, verbose = FALSE)
S.O.integrated <- RunPCA(S.O.integrated, npcs = 30, verbose = FALSE)
S.O.integrated <- RunUMAP(S.O.integrated, reduction = "pca", dims = 1:30)
S.O.integrated <- FindNeighbors(S.O.integrated, reduction = "pca", dims = 1:30)
S.O.integrated <- FindClusters(S.O.integrated, resolution = 0.12)

S.O.integrated@meta.data$spp <- factor(S.O.integrated.merge@meta.data$spp, 
                                       levels = unique(S.O.integrated.merge@meta.data$spp))

S.O.integrated@meta.data$phase.spp <- paste(S.O.integrated@meta.data$spp, S.O.integrated@meta.data$cell.cycle.phase, sep = ':')

saveRDS(S.O.integrated, '../Input/compScBabesia/RData/S.O.integrated.list.pstime.GAM.indiv.cell.cycle.phase_sinlge_obj.Rdata')

S.O.integrated <- readRDS('../Input/compScBabesia/RData/S.O.integrated.list.pstime.GAM.indiv.cell.cycle.phase_sinlge_obj.Rdata')

Idents(S.O.integrated) <- "cell.cycle.phase"
p <- DimPlot(S.O.integrated, reduction = "pca", 
             #group.by = "cells", 
             split.by = 'spp',
             pt.size = 1,
             shape.by='spp',
             label = TRUE, label.size = 4, ncol = 2) + NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold"))

plot(p)

ggsave(filename="../Output/compScBabsSpecies/figs/integrated_pca_comm_gene_inferred_phases.png",
       plot=p,
       width = 13, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


## Cell cycle marker: Independently done per spp (common genes)
cell.cycle.markers.sig <- getCellCyclePhaseMarkers(S.O.integrated)

cell.cycle.markers.sig <- left_join(cell.cycle.markers.sig, prod.desc, by = 'GeneID') 
cell.cycle.markers.sig <- left_join(cell.cycle.markers.sig, toxo.bdiv.orth, by='GeneID') %>% arrange(spp, desc(avg_log2FC))
write.xlsx(cell.cycle.markers.sig, '../Output/compScBabsSpecies/tables/per_spp_inferred_cell_cycle_markers_sig_int_data_comm_genes.xlsx')

cell.cycle.markers.sig <- read.xlsx('../Output/compScBabsSpecies/tables/per_spp_inferred_cell_cycle_markers_sig_int_data_comm_genes.xlsx')
cell.cycle.markers.top <- cell.cycle.markers.sig %>% group_by(spp, phase) %>% slice_max(n = 2, order_by = avg_log2FC)
cell.cycle.markers.stat <- cell.cycle.markers.sig %>% group_by(spp, phase) %>% summarise(num.deg = n())



# bar plot 
cell.cycle.markers.stat$phase <- factor(cell.cycle.markers.stat$phase, levels = c("G", "SM", "MC", "C"))
p <-  ggplot(cell.cycle.markers.stat, aes(x=phase, y=num.deg)) +
  geom_bar(stat="identity", fill = "steelblue")+
  geom_text(aes(label=num.deg), vjust=2, color="black", size=6, fontface = 'bold')+
  theme_bw()+
  facet_wrap(.~spp, scales = 'free',labeller=label_wrap_gen(multi_line = TRUE))+
  theme(axis.text.x = element_text(face="bold", size=12, angle=0),
        axis.text.y = element_text(face="bold", size=12, angle=0),
        axis.title.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12),
        strip.placement = "outside")+
  theme(axis.line = element_line(color = 'black'), 
        axis.ticks = element_blank())

p 

ggsave(filename="../Output/compScBabsSpecies/figs/per_spp_inferred_cell_cycle_markers_sig_int_data_comm_genes_bar_plot.png",
       plot=p,
       width = 10, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



#################################################################
## Cross spp differential expression: Phase-based unique markers
## spp specific cell cycle markers 
######## Identity based comparison. EX: compare C to C accross all spps
#ident.1 case, ident.2 is control
DefaultAssay(S.O.integrated) <- "RNA"
Idents(S.O.integrated) <- "phase.spp"


contrasts <- makeMatchedContrasts(S.O.integrated)
contrasts.groups <- contrasts %>% group_by(ref) %>% summarise(query = list(query))
contrasts.groups$phase <- gsub('.*:', '', contrasts.groups$ref)
contrasts.groups$ref.spp<- gsub(':.*', '', contrasts.groups$ref)
matched.DEGs <- mclapply(1:nrow(contrasts.groups), function(i){
  tmp <- FindMarkers(S.O.integrated, ident.1 = contrasts.groups$ref[i],
                     ident.2 = c(unlist(contrasts.groups$query[i])), 
                     only.pos = T, verbose = T)
  tmp$ref <- contrasts.groups$ref[i]
  tmp$ref.spp <- contrasts.groups$ref.spp[i]
  tmp$phase <- contrasts.groups$phase[i]
  tmp$gene <- rownames(tmp)
  tmp$GeneID <- gsub('-', '_', tmp$gene)
  return(tmp)
})
markers.int.local <- bind_rows(matched.DEGs)
markers.int.local.sig <- markers.int.local %>% dplyr::filter(avg_log2FC > 0.58 & p_val_adj < 0.01)
markers.int.local.sig <- markers.int.local.sig  %>% mutate(pct.ratio = (pct.1 / pct.2))

## Now consider markers that are phase specific
markers.int.local.sig.cc <- inner_join(markers.int.local.sig, cell.cycle.markers.sig, 
                                       by = c('GeneID', 'ref.spp' = 'spp', 'phase'))

markers.int.local.sig.specific <- markers.int.local.sig.cc %>% 
  arrange(GeneID, ref.spp, desc(pct.ratio)) %>% group_by(GeneID) %>% slice_max(n = 1, order_by = pct.ratio)
markers.int.local.sig.specific.df <- markers.int.local.sig.specific %>% dplyr::select(!contains(".y"))

#write.xlsx(markers.int.local.sig.specific.df, '../Output/compScBabsSpecies/tables/cross_spp_phase_based_markers_sig_fc_1_5_rev.xlsx')

# markers.int.local.top <- markers.int.local.sig.specific %>% group_by(phase, ref.spp) %>% 
#   slice_max(n = 2, order_by = avg_log2FC)


#### Here we need to generate a bar plot of spp specific local markers total numbers.
markers.int.local.sig.stats <- markers.int.local.sig.specific.df %>% group_by(phase, ref.spp) %>% 
  summarise(genes = list(GeneID), num.deg = n())

#saveRDS(markers.int.local.sig.stats, "../Input/compScBabesia/RData/spp_specific_unique_markers_fc_1_5_rev.RData")
markers.int.local.sig.stats$phase <- factor(markers.int.local.sig.stats$phase, levels = c('G', 'SM', 'MC', 'C'))

tmp <- unlist(markers.int.local.sig.stats$genes[markers.int.local.sig.stats$ref.spp == "Bdiv_human" & markers.int.local.sig.stats$phase=="C"])

colnames(markers.int.local.sig.stats) <- gsub("ref.spp", "spp", colnames(markers.int.local.sig.stats))
markers.int.local.sig.stats$phase <- factor(markers.int.local.sig.stats$phase, levels = c('G', "SM", "MC", "C"))
p <-  ggplot(markers.int.local.sig.stats, aes(x=num.deg, y=spp, fill = spp)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("Bbig" = "firebrick","Bbov" ="darkorchid3",
                               'Bdiv_cow' = 'darkslateblue', 'Bdiv_human' = 'darkolivegreen4'))+
  geom_text(aes(label=num.deg), hjust= 1.15, color="black", size=6, fontface = 'bold')+
  theme_bw()+
  facet_grid(phase~., scales = "free", space='free',switch = "y",
             labeller=label_wrap_gen(multi_line = TRUE))+
  theme(#axis.text.x = element_blank(),
    axis.text.x = element_text(face="bold", size=12, angle=0, hjust = 1),
    axis.title.y = element_text(size=14, face="bold"),
    axis.title.x = element_text(size=14, face="bold"),
    #axis.text.y = element_text(face="bold", size=12, angle=0),
    axis.text.y = element_blank())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12,  angle = 0), 
        strip.placement = "outside")+
  theme(axis.line = element_line(color = 'black'), 
        axis.ticks = element_blank())+
  theme(legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 14))+
  theme(strip.text.y.left = element_text(angle = 0))

p



ggsave(filename="../Output/compScBabsSpecies/figs/cross_spp_uniq_phase_based_markers.png",
       plot=p,
       width = 10, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


# p <- FeaturePlot(object = S.O.integrated, features = markers.int.local.top$gene[3], label = T,repel = T,
#                  #shape.by  = 'spp',
#                  split.by = 'spp',
#                  cols = c("grey", "blue"), reduction = "pca")
# 
# plot(p)

# table to do GO term

cross.spp.phase.based <- markers.int.local.sig.specific.df %>% group_by(phase, ref.spp) %>% na.omit() %>%
  summarise(genes = list(GeneID_Toxo), num.deg = n())
#write.xlsx(cross.spp.phase.based, "../Input/compScBabesia/GO/cross_spp_phase_based_markers/cross_spp_uniq_phase_based_markers_TGGT1_orth_summarized.xlsx")
write.xlsx(cross.spp.phase.based, "../Output/compScBabsSpecies/tables/cross_spp_phase_based_markers_TGGT1_orth_summarized_1_5_fc_rev.xlsx")


############# Find conserved markers in matched inferred phases

DefaultAssay(S.O.integrated) <- 'RNA'
Idents(S.O.integrated) <- 'cell.cycle.phase'
phases  <- c('G', 'SM', 'MC', 'C')

cons.markers <- lapply(phases, function(cc){
  cons.marker <- FindConservedMarkers(S.O.integrated, ident.1 = cc, grouping.var = "spp", verbose = T)
  cons.marker$gene <- rownames(cons.marker)
  cons.marker$GeneID <- gsub('-', '_', cons.marker$gene)
  cons.marker$phase <- cc
  return(cons.marker)
})



cons.markers.all.phases <- do.call("rbind", cons.markers)

cons.FC <- cons.markers.all.phases[grepl('GeneID|log2FC|phase', names(cons.markers.all.phases))]
colnames(cons.FC) <- gsub("_avg_log2FC", "", colnames(cons.FC))
cons.FC.lng <- cons.FC %>% pivot_longer(!c(GeneID , phase),names_to = "spp" , values_to = 'max_avg_log2FC')

cons.FC.max <- cons.FC.lng %>%  arrange(GeneID,  desc(max_avg_log2FC)) %>% group_by(GeneID) %>% slice_max(n= 1, order_by= max_avg_log2FC)
cons.FC.max <- cons.FC.max %>% dplyr::select(GeneID, phase,spp, max_avg_log2FC)

cons.markers.phases <- left_join(cons.FC.max, cons.markers.all.phases, by = c('GeneID',"phase")) 
cons.markers.phases.sig <- cons.markers.phases %>% dplyr::filter(max_pval < 0.01 & max_avg_log2FC > 1)
cons.markers.phases.sig.df <- cons.markers.phases.sig %>%
  dplyr::select(GeneID, phase, spp, gene, max_avg_log2FC, minimump_p_val, max_pval )


# ## Now consider markers that are phase specific
cons.markers.phases.sig.cc <- inner_join(cons.markers.phases.sig.df, cell.cycle.markers.sig, 
                                         by = c("gene" ,"GeneID", "phase", "spp"))

write.xlsx(cons.markers.phases.sig.cc, '../Output/compScBabsSpecies/tables/conserved_markers_cross_spp_phase_based_sig.xlsx')
saveRDS(cons.markers.phases.sig.cc, "../Input/compScBabesia/RData/conserved_markers_cross_spp_phase_based.RData")

## GO term tab
cons.markers.GT1.orth <- cons.markers.phases.sig.cc %>% na.omit() %>%
  group_by(phase) %>% summarise(genes= list(GeneID_Toxo),  num.genes = n())
write.xlsx(cons.markers.GT1.orth, "../Output/compScBabsSpecies/tables/conserved_markers_cross_spp_phase_based_TGGT1_orth_summarized.xlsx")



cons.markers.phases.sig.cc.stat <- cons.markers.phases.sig.cc %>% group_by(phase) %>% 
  summarise(genes = list(GeneID), num.deg = n())

cons.markers.phases.sig.cc.stat$phase <- factor(cons.markers.phases.sig.cc.stat$phase, levels = c('G', 'SM', 'MC', 'C'))
# 
# SM: "#c44237"
# 
# MC: "#ad8842"
# 
# G: "#1b5878"
# 
# C: "#e99093"

p <-  ggplot(cons.markers.phases.sig.cc.stat, aes(x=num.deg, y=phase, fill = phase)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = c('G' = "#1b5878","SM" ="#c44237",
                               'MC' = '#ad8842', 'C' = '#e99093'))+
  geom_text(aes(label=num.deg), hjust=1.25, color="black", size=6, fontface = 'bold')+
  theme_bw()+
  theme(axis.text.x = element_text(face="bold", size=12, angle=0, hjust = 1),
        axis.text.y = element_text(face="bold", size=12, angle=0),
        axis.title.y = element_text(size=14, face="bold"),
        axis.title.x = element_text(size=14, face="bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"))+
  theme(axis.line = element_line(color = 'black'), 
        axis.ticks = element_blank())+
  theme(legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 14))+
  theme(legend.justification = "top")
p


ggsave(filename="../Output/compScBabsSpecies/figs/unique_conserved_markers_bar_plot_updated.png",
       plot=p,
       width = 6, height = 4,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)




### conserved heatmap
sc.tc.df.adj <- readRDS('../Input/compScBabesia/RData/all_sc_tc_df_adj.RData')
sc.tc.df.adj <- sc.tc.df.adj[-5]

sc.tc.fits <- readRDS('../Input/compScBabesia/RData/all_sme_fits_sc_tc_20min.RData')
sc.tc.fits <- sc.tc.fits[-5]


cons.markers.phases.sig <- readRDS("../Input/compScBabesia/RData/Unique_conserved_markers_phases_based.RData")
conserved.genes <- cons.markers.phases.sig$GeneID

sc.tc.mus <- lapply(1:length(sc.tc.fits), function(i){
  sc.tc.mu <- smoothSplineSmeFits(sc.tc.fits[[i]], unique(sc.tc.df.adj[[i]]$variable), extend = F)
  colnames(sc.tc.mu) <- c('GeneID', 't', 'y')
  sc.tc.mu <- sc.tc.mu %>% dplyr::filter(GeneID %in% cons.markers.phases.sig$GeneID) %>%
    pivot_wider(names_from = 'GeneID', values_from = 'y') %>%
    as.data.frame()
  
  return(sc.tc.mu)
})

names(sc.tc.mus) <- names(sc.tc.fits)

## Generate the clusters
num.clust <- 5L

sc.hc_dtws <- lapply(sc.tc.mus, function(x){
  sc.hc_dtw <- dtwClustCurves(x[,2:ncol(x)], nclust = num.clust)
  return(sc.hc_dtw)
})

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
saveRDS(sc.tc.mus.scale,'../Input/compScBabesia/RData/all_sc_tc_mus_scale_conserved_phase_based_markers.RData')

sc.tc.mus.scale <- readRDS('../Input/compScBabesia/RData/all_sc_tc_mus_scale_conserved_phase_based_markers.RData')
titles <- c("B. big", "B. bov", "B. div (cow)", "B. div (human)")
spps <- c('Bbig', 'Bbov', 'Bdiv_cow', 'Bdiv_human')

for(i in 1:length(titles)){
  sc.tc.mus.scale[[i]]$spp <- spps[i]
}


all.sc.tc.mus.scale <- do.call("rbind", sc.tc.mus.scale)
#all.sc.tc.mus.scale.com <- all.sc.tc.mus.scale %>% filter(GeneID %in% com.genes)
#all.sc.tc.mus.scale.com <- all.sc.tc.mus.scale %>% filter(GeneID %in% conserved.genes)

# get  the order of Bdiv human to  order the genes 
Bdiv_hum <- all.sc.tc.mus.scale %>% filter(spp == "Bdiv_human") %>% dplyr::select(GeneID, peak.ord) %>% distinct()
names(Bdiv_hum) <- gsub("peak.ord", "bdiv_peak.ord", names(Bdiv_hum)) 

all.sc.tc.mus.scale  <- inner_join(all.sc.tc.mus.scale, Bdiv_hum, by  = "GeneID")
all.sc.tc.mus.scale$aug.time <- all.sc.tc.mus.scale$t

for(i in 1:length(spps)){
  ind <- all.sc.tc.mus.scale$spp == spps[i]
  all.sc.tc.mus.scale$aug.time[ind] <- all.sc.tc.mus.scale$aug.time[ind] + (12 + 1/3) * (i - 1)
}

p <- ggplot(all.sc.tc.mus.scale, aes(x = aug.time, reorder(GeneID, -bdiv_peak.ord),  fill = y)) +
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
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.position = "none")+
  theme(strip.background=element_rect(fill='white', color = 'black'))+
  theme(panel.spacing = unit(1.5, "lines")) +
  theme(strip.text.x=element_text(angle=0, hjust=0.5,vjust=0.5, size = "14",face = 'bold'))+
  theme(
    plot.title = element_text(size=18, face = "bold.italic", color = 'black'),
    axis.title.x = element_text(size=16, face="bold"),
    axis.title.y = element_text(size=18, face="bold")
  )

p

ggsave(filename="../Output/compScBabsSpecies/figs/gene_cluster_curves_common_conserved_phase_based_markers.png",
       plot=p,
       width = 20, height = 3,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)




################################################################
## second approach of finding conserved and spp specific markers
## looking at the iintersection and differences 


markers.sig <- readRDS('../Input/compScBabesia/RData/all.markers.sig.cell.cycle.phase.RData')
clusters <-  c('G', 'SM', 'MC', 'C')
spps <- c('Bbig', 'Bbov', 'Bdiv_cow', 'Bdiv_human')
titles <- c("B. big", "B. bov", "B. div (cow)", "B. div (human)")

for(i in 1:length(titles)){
  markers.sig[[i]]$markers.sig$spp <- spps[i]
}

all.markers.sig <- do.call(Map, c(f = rbind, markers.sig))
all.markers.sig <- all.markers.sig$markers.sig
all.markers.sig$spp <-  factor(all.markers.sig$spp, levels = c('Bbig', 'Bbov', 'Bdiv_cow', 'Bdiv_human'))
names(all.markers.sig) <- gsub("cluster", "phase", names(all.markers.sig))
all.markers.sig$phase <- factor(all.markers.sig$phase, levels = c('G', 'SM', 'MC', 'C'))
all.markers.sig.stat <- all.markers.sig %>% group_by(spp, phase) %>% summarise(num.deg = n()) 

p <- ggplot(all.markers.sig.stat, aes(x=phase, y=num.deg, fill = phase, color = phase)) +
  geom_bar(stat="identity",size = 2, width = 0.75, aes(fill = phase))+
  scale_color_manual(values = c('G' = "#1b5878","SM" ="#c44237",'MC' = '#ad8842', 'C' = '#e99093'))+
  scale_fill_manual(values = c('G' = "#758FA2","SM" = "#D4807F",'MC' = "#C3AE85", 'C' = "#EAADAF"))+
  geom_text(aes(label=num.deg), vjust=1.15, color="black", size=10, fontface = 'bold')+
  facet_grid(. ~ spp, scales = "free", space='free',labeller=label_wrap_gen(multi_line = TRUE))+
  theme_bw()+
  theme(
    axis.text.x = element_blank(),
    #axis.text.y = element_text(size = 10),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_text(size = 24, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold"))+
  theme(strip.background=element_rect(fill='white', color = 'black'),
        panel.spacing = unit(1.5, "lines"), 
        strip.text.x=element_text(angle=0, hjust=0.5,vjust=0.5, size = 24,face = 'bold'),
        plot.title = element_text(size=18, face = "bold.italic", color = 'black'),
        axis.title.x = element_text(size=22, face="bold"),
        axis.title.y = element_text(size=22, face="bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.text = element_text(face = "bold", size = 24,  angle = 0), 
        strip.placement = "outside") +
  theme(legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 14))+
  theme(panel.spacing = unit(1.5, "lines")) 
#theme(legend.position = "none")

p


all.markers <- all.markers.sig %>% group_by(spp, phase) %>% summarise(genes = list(GeneID), total = n()) %>%
  group_by(phase) %>% 
  mutate(intersect = list(reduce(genes, intersect))) %>% ungroup() %>% mutate(intersect.length = lengths(intersect))

########3 shared markers in phases across all spp

all.markers.shared.stat <- all.markers.sig %>% group_by(spp, phase) %>% summarise(genes = list(GeneID), total = n()) %>%
  group_by(phase) %>% 
  summarise(intersect = list(reduce(genes, intersect)), intersect.length = lengths(intersect)) 

write.xlsx(all.markers.shared.stat, "../Output/compScBabsSpecies/tables/shared_phase_markers_across_speciies_intersect_summ.xlsx")

all.markers.shared.stat$phase <- factor(all.markers.shared.stat$phase, levels = c('C', 'MC', 'SM', 'G'))
p <-  ggplot(all.markers.shared.stat, aes(x=intersect.length, y=phase, fill = phase, color = phase)) +
  geom_bar(stat="identity",  size = 2, width = 0.75, aes(fill = phase))+
  scale_color_manual(values = c('G' = "#1b5878","SM" ="#c44237",
                                'MC' = '#ad8842', 'C' = '#e99093'))+
  scale_fill_manual(values = c('G' = "#758FA2","SM" = "#D4807F",
                               'MC' = "#C3AE85", 'C' = "#EAADAF"))+
  
  geom_text(aes(label=intersect.length), hjust=1.15, color="black", size=12, fontface = 'bold')+
  theme_bw()+
  theme(axis.text.x = element_text(face="bold", size=30, color = "black",angle=0, hjust = 1),
        axis.text.y = element_text(face="bold", size=30, color = "black", angle=0),
        axis.title.y = element_text(size=28, face="bold"),
        axis.title.x = element_text(size=28, face="bold"))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"))+
  theme(axis.line = element_line(color = 'black'), 
        axis.ticks = element_blank())+
  theme(legend.text = element_text(face = "bold", size = 18),
        legend.title = element_text(face = "bold", size = 20))+
  theme(legend.justification = "top")+ 
  scale_x_continuous(expand = c(0,0))  +
  xlab("num.deg")+ 
  ylab("")
p

ggsave(filename="../Output/compScBabsSpecies/figs/spp_specific_markers_shared.png",
       plot=p,
       width = 8, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


####### speciesiec specific markers in each phase subtracting the interct 



all.markers.stat <- all.markers.sig %>% group_by(spp, phase) %>% summarise(genes = list(GeneID), total = n())

XX <- full_join(all.markers.stat, all.markers.shared.stat, by = "phase")
XX.diff <- XX %>% rowwise() %>% 
  mutate(spp.genes = list(setdiff(genes, intersect)), num.spp.genes = length(unlist(setdiff(genes, intersect))))

spp.specific.markers <- XX.diff %>% select(spp, phase, spp.genes,num.spp.genes)

spp.specific.markers$phase <- factor(spp.specific.markers$phase, levels = c('G', 'SM', 'MC', 'C'))

write.xlsx(spp.specific.markers, "../Output/compScBabsSpecies/tables/spp_specific_markers_set_diff_summ.xlsx")
# Bbov = "#D09FE9"
# Bbig = "#E8B5B4"
# Bdiv_cow = "#BFB9D6"
# Bdiv_human = "#CED5BD"

p <-
  ggplot(spp.specific.markers, aes(x=num.spp.genes, y=spp, fill = spp, color = spp)) +
  geom_bar(stat="identity", size = 2, width = 0.75, aes(fill=spp)) +
  scale_color_manual(values = c("Bbig" = "firebrick","Bbov" ="darkorchid3",
                                'Bdiv_cow' = 'darkslateblue', 'Bdiv_human' = 'darkolivegreen4'))+
  scale_fill_manual(values =c("Bbig" = "#E8B5B4","Bbov" ="#D09FE9",
                              'Bdiv_cow' = "#BFB9D6", 'Bdiv_human' = "#CED5BD"))+
  geom_text(aes(label=num.spp.genes), hjust= 1, color="black", size=10, fontface = 'bold')+
  theme_bw()+
  facet_grid(phase~., scales = "free", space='free',switch = "y",
             labeller=label_wrap_gen(multi_line = TRUE)) +
  theme(#axis.text.x = element_blank(),
    axis.text.x = element_text(face="bold", color = "black", size=20, angle=0, hjust = 1),
    axis.title.y = element_text(size=30, face="bold"),
    axis.title.x = element_text(size=30, face="bold"),
    #axis.text.y = element_text(face="bold", size=12, angle=0),
    axis.text.y = element_blank())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 28,  angle = 0),
        strip.placement = "outside")+
  theme(axis.line = element_line(color = 'black'),
        axis.ticks = element_blank())+
  theme(legend.text = element_text(face = "bold", size = 20),
        legend.title = element_text(face = "bold", size = 18))+
  theme(strip.text.y.left = element_text(angle = 0)) +
  scale_x_continuous(expand = c(0,0))+
  xlab("num.deg") + 
  ylab("")

p


ggsave(filename="../Output/compScBabsSpecies/figs/spp_specific_markers_set_diff.png",
       plot=p,
       width = 10, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



################################################################
# 
##  Bdiv_human vs Bdiv_cow: Phase_based

### B div human vs cow: Global unique
DefaultAssay(S.O.integrated) <- 'RNA'
Idents(S.O.integrated) <- 'spp'


S.O.int.bdivs <- subset(S.O.integrated, ident = c('Bdiv_human', 'Bdiv_cow'))
DefaultAssay(S.O.int.bdivs) <- 'RNA'
Idents(S.O.int.bdivs) <- 'phase.spp'

contrasts <- makeMatchedContrasts(S.O.integrated)
contrasts.groups <- contrasts %>% group_by(ref) %>% summarise(query = list(query))
contrasts.groups$phase <- gsub('.*:', '', contrasts.groups$ref)
contrasts.groups$ref.spp<- gsub(':.*', '', contrasts.groups$ref)

contrasts.bdivs <- contrasts %>% dplyr::filter(grepl('Bdiv', ref.spp) & grepl('Bdiv', query.spp))
contrasts.bdivs$phase <- gsub('.*:', '', contrasts.bdivs$ref)
contrasts.bdivs$ref.spp<- gsub(':.*', '', contrasts.bdivs$ref)

matched.DEGs.bdivs <- mclapply(1:nrow(contrasts.bdivs), function(i){
  tmp <- FindMarkers(S.O.int.bdivs, ident.1 = contrasts.bdivs$ref[i],
                     ident.2 = contrasts.bdivs$query[i], 
                     only.pos = T, verbose = T)
  tmp$ref <- contrasts.bdivs$ref[i]
  tmp$ref.spp <- contrasts.bdivs$ref.spp[i]
  tmp$phase <- contrasts.bdivs$phase[i]
  tmp$gene <- rownames(tmp)
  tmp$GeneID <- gsub('-', '_', tmp$gene)
  return(tmp)
})
markers.int.local.bdivs <- bind_rows(matched.DEGs.bdivs)
markers.int.local.bdivs.sig <- markers.int.local.bdivs %>% dplyr::filter(avg_log2FC > 0.58 & p_val_adj < 0.01)
markers.int.local.bdivs.sig <- markers.int.local.bdivs.sig  %>% mutate(pct.ratio = (pct.1 / pct.2))

## Now consider markers that are phase specific
markers.int.local.bdivs.sig.cc <- inner_join(markers.int.local.bdivs.sig, cell.cycle.markers.sig, 
                  by = c('GeneID', 'ref.spp' = 'spp', 'phase'))
markers.int.local.bdivs.sig.cc <- markers.int.local.bdivs.sig.cc %>% dplyr::select(!contains(".y"))

markers.int.local.bdivs.sig.df <- markers.int.local.bdivs.sig.cc %>% 
  arrange(GeneID, desc(pct.ratio)) %>% group_by(GeneID, ref.spp) %>% slice_max(n = 1, order_by = pct.ratio)


write.xlsx(markers.int.local.bdivs.sig.df, '../Output/compScBabsSpecies/tables/bdiv_human_vs_bdiv_cow_phase_based_markers.xlsx')
saveRDS(markers.int.local.bdivs.sig.df, "../Input/compScBabesia/RData/bdiv_human_vs_bdiv_cow_phase_based_markers.RData")

markers.int.local.bdivs.orth <- markers.int.local.bdivs.sig.df %>% group_by(ref.spp, phase) %>% 
  na.omit() %>% summarise(genes = list(GeneID_Toxo), num.genes = n())
write.xlsx(markers.int.local.bdivs.orth, '../Output/compScBabsSpecies/tables/bdiv_human_vs_bdiv_cow_phase_based_TGGT1_orth_summarized.xlsx')

## bar plot 
markers.int.local.bdivs.sig.stat <- markers.int.local.bdivs.sig.df %>% group_by(ref.spp, phase) %>% 
  summarise(genes = list(GeneID), num.deg = n())

colnames(markers.int.local.bdivs.sig.stat) <- gsub("ref.spp", "spp", colnames(markers.int.local.bdivs.sig.stat))
markers.int.local.bdivs.sig.stat$phase <- factor(markers.int.local.bdivs.sig.stat$phase, levels = c('G', 'SM', 'MC', 'C'))

p <-  ggplot(markers.int.local.bdivs.sig.stat, aes(x=num.deg, y=spp, fill = spp)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = c('Bdiv_cow' = 'darkslateblue', 'Bdiv_human' = 'darkolivegreen4'))+
  geom_text(aes(label=num.deg), hjust= 1.75, color="black", size=6, fontface = 'bold')+
  theme_bw()+
  facet_grid(phase~., scales = "free", space='free',switch = "y",
             labeller=label_wrap_gen(multi_line = TRUE))+
  theme(#axis.text.x = element_blank(),
    axis.text.x = element_text(face="bold", size=12, angle=0, hjust = 1),
    axis.title.y = element_text(size=14, face="bold"),
    axis.title.x = element_text(size=14, face="bold"),
    #axis.text.y = element_text(face="bold", size=12, angle=0),
    axis.text.y = element_blank())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12,  angle = 0), 
        strip.placement = "outside")+
  theme(axis.line = element_line(color = 'black'), 
        axis.ticks = element_blank())+
  theme(legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 14))+
  theme(strip.text.y.left = element_text(angle = 0))

p

ggsave(filename = "../Output/compScBabsSpecies/figs/Bdiv_human_vs_Bdiv_cow_phase_based_markers_bar_plot.png", 
       plot = p, 
       width = 8, height = 4, 
       dpi = 300)

## expression plot 
markers.int.local.bdivs.sig.top <- markers.int.local.bdivs.sig.df %>% group_by(ref.spp) %>% slice_max(n = 1, order_by = avg_log2FC.x)
S.O.int.bdivs@meta.data <- S.O.int.bdivs@meta.data %>% mutate(host = ifelse( spp == "Bdiv_human", "Bdiv_human", "Bdiv_cow"))
S.O.int.bdivs[['pca']]@cell.embeddings[,2] = -S.O.int.bdivs[['pca']]@cell.embeddings[,2]
S.O.int.bdivs[['pca']]@cell.embeddings[,1] = -S.O.int.bdivs[['pca']]@cell.embeddings[,1]

p <- FeaturePlot(object = S.O.int.bdivs, 
                 features = rev.default(markers.int.local.bdivs.sig.top$gene.x), label = F,repel = F,
                 #shape.by  = 'spp',
                 split.by = 'host',
                 cols = c("grey", "blue"), reduction = "pca")

plot(p)
num.plt <- 4
p2 <- lapply(1:num.plt, function(i){
  
  plt <- p[[i]] + 
    theme(axis.text.x  = element_text(face = "bold", size = 16, angle = 0, hjust = 0.5),
          axis.text.y  = element_text(face = "bold", size = 18),
          axis.title = element_blank(),
          #axis.title.x = element_text(face = "bold", size = 20), 
          #axis.title.y = element_text(face = "bold", size = 20),
          plot.title = element_text(face = "bold", size = 20)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #panel.border = element_rect(color = "black", size = 0.5),
          axis.line.x = element_line(colour = "black", size = 0.5),
          axis.line.y = element_line(colour = "black", size = 0.5)
          
    ) 
  
  
})

p2 <- grid.arrange(p2[[1]], p2[[2]], p2[[3]], p2[[4]], nrow = 2)


ggsave(filename="../Output/compScBabsSpecies/figs/bdiv_human_vs_bdiv_cow_featureplot_V1.png",
       plot=p2,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)




# p <- FeaturePlot(object = S.O.int.bdivs, features = markers.int.local.bdivs.sig.top$gene.x, label = T,repel = F,
#                  #shape.by  = 'spp',
#                  split.by = 'spp',
#                  cols = c("grey", "blue"), reduction = "pca")
# 
# 
# plot(p)

exprs_bdivs <- exprExprData(S.O.int.bdivs, markers.int.local.bdivs.sig.top$GeneID)
exprs_bdivs$spp[grepl('cow', exprs_bdivs$Sample)] <- 'Bdiv_cow'
exprs_bdivs$spp[grepl('human', exprs_bdivs$Sample)] <- 'Bdiv_human'

col_range <- colorRampPalette(c('white', 'blue'))
p <- ggplot(exprs_bdivs, aes(x=-PC_1,y=-PC_2)) +
  geom_point(aes(fill = expr,
                 #color = ""
  ), color = 'lightgray',
  shape=21, size = 1.5)+ 
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  ylab('PC2') + xlab('PC1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 20, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 20, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid"),
        strip.text = element_text(size = 28, face = "bold", color = "black")) +
  facet_grid(GeneID ~ spp) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=20, face="bold", hjust = 1),
    axis.title.y = element_text(size=20, face="bold"), 
    legend.title = element_text(size = 20, face = "bold"), 
    legend.text = element_text(size = 18, color = "black"))+
  theme(
    #axis.text.x = element_blank(),
    axis.text.x = element_text(face="bold", size=20, angle=0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(face="bold", size=20, color = "black"),
    axis.title.y = element_text(size=20, face="bold"),
    axis.title.x = element_text(size=20, face="bold"),
    #axis.text.y = element_text(face="bold", size=12, angle=0),
    #axis.text.y = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.border = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.y = element_blank()
  )+
  theme(panel.spacing = unit(0.5, "lines"))

plot(p)
ggsave(filename="../Output/compScBabsSpecies/figs/expr_Bdiv_human_vs_Bdiv_cow_top_1_phase_based_V2_rev.png",
       plot=p,
       width = 10, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

## violin plot of top markers
DefaultAssay(S.O.int.bdivs) <- 'RNA'
Idents(S.O.int.bdivs) <- 'spp'
S.O.int.bdivs@meta.data$spp <- factor(S.O.int.bdivs@meta.data$spp, levels = c("Bdiv_cow", "Bdiv_human"))
top.markers <- markers.int.local.bdivs.sig.top$GeneID

#VlnPlot(subset(S.O.int.bdivs, idents = c('Brady', 'Intra')), features = gsub('_', '-', show.lab))
p <- VlnPlot(S.O.int.bdivs, features = gsub('_', '-', top.markers))
p + theme(axis.text  = element_text(face = "bold", size = 14))

pp <- lapply(1:2, function(i){
  
  plt <- p[[i]] + 
    scale_fill_manual(values = c('Bdiv_cow' = 'darkslateblue', 'Bdiv_human' = 'darkolivegreen4'))+
    theme(axis.text.x  = element_text(face = "bold", size = 16, angle = 0, hjust = 0.5),
          axis.text.y  = element_text(face = "bold", size = 18),
          axis.title.x = element_blank(), 
          axis.title.y = element_text(face = "bold", size = 20),
          plot.title = element_text(face = "bold", size = 20)) + 
    ylab("expr")
  
  
})

p2 <- grid.arrange(pp[[2]], pp[[1]], nrow = 1)

ggsave(filename="../Output/compScBabsSpecies/figs/vln_expr_Bdiv_human_vs_Bdiv_cow_top_1_phase_based_rev.png",
       plot=p2,
       width = 10, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

#########################
##############################################
###  human vs cow: phase based
DefaultAssay(S.O.integrated) <- "RNA"
S.O.integrated@meta.data <- S.O.integrated@meta.data %>% mutate(host = ifelse( spp == "Bdiv_human", "human", "cow"))
S.O.integrated@meta.data$host.phase <- paste(S.O.integrated@meta.data$host, S.O.integrated@meta.data$cell.cycle.phase, sep = ":")

Idents(S.O.integrated) <- 'host.phase'

makeHostMatchedContrasts <- function(S.O.integrated){
  
  objs <- as.character(unique(S.O.integrated@meta.data$host.phase))
  
  contrasts <- data.frame(ref = objs, dummy = 1) %>% full_join( data.frame(query = objs, dummy = 1), by = 'dummy') %>% 
    mutate(ref.host = gsub(":.*", "", ref), ref.phase = gsub(".*:", "", ref), 
           query.host = gsub(":.*", "", query), query.phase = gsub(".*:", "", query))
  my.contrasts <- contrasts %>% dplyr::filter(ref.phase == query.phase & ref.host != query.host)
  
  return(my.contrasts)
  
}

contrasts <- makeHostMatchedContrasts(S.O.integrated)
contrasts.groups <- contrasts %>% group_by(ref) %>% summarise(query = list(query))
contrasts.groups$phase <- gsub('.*:', '', contrasts.groups$ref)
contrasts.groups$ref.host <- gsub(':.*', '', contrasts.groups$ref)

matched.DEGs <- mclapply(1:nrow(contrasts.groups), function(i){
  tmp <- FindMarkers(S.O.integrated, ident.1 = contrasts.groups$ref[i],
                     ident.2 = c(unlist(contrasts.groups$query[i])), 
                     only.pos = T, verbose = T)
  tmp$ref <- contrasts.groups$ref[i]
  tmp$ref.host <- contrasts.groups$ref.host[i]
  tmp$phase <- contrasts.groups$phase[i]
  tmp$gene <- rownames(tmp)
  tmp$GeneID <- gsub('-', '_', tmp$gene)
  return(tmp)
})


host.markers.int.local <- bind_rows(matched.DEGs)
host.markers.int.local.sig <- host.markers.int.local %>% dplyr::filter(avg_log2FC > 0.58 & p_val_adj < 0.01)
host.markers.int.local.sig  <- host.markers.int.local.sig %>% mutate(pct.ratio = (pct.1 / pct.2))

cell.cycle.markers.sig <- cell.cycle.markers.sig %>% mutate(host = ifelse( spp == "Bdiv_human", "human", "cow"))
host.markers.int.local.sig.cc <- inner_join(host.markers.int.local.sig, cell.cycle.markers.sig, 
                  by = c('GeneID', 'phase', 'ref.host' = 'host'))
host.markers.int.local.sig.cc <- host.markers.int.local.sig.cc %>% dplyr::select(!contains(".y"))
# unique to spp and phase
host.markers.int.local.sig.specific <- host.markers.int.local.sig.cc  %>%  distinct(GeneID, .keep_all = T) %>%
  arrange(GeneID, desc(pct.ratio)) %>% group_by(GeneID) %>% slice_max(n = 1, order_by = pct.ratio)


write.xlsx(host.markers.int.local.sig.specific, '../Output/compScBabsSpecies/tables/bdiv_human_vs_all_spp_in_cow_phase_based_markers_sig.xlsx')
saveRDS(host.markers.int.local.sig.specific, '../Input/compScBabesia/RData/bdiv_human_vs_all_spp_in_cow_phase_based_markers_sig.RData')
## GO term table
host.int.local.markers.orth <- host.markers.int.local.sig.specific %>% group_by(ref.host, phase) %>% na.omit() %>%
  summarise(genes = list(GeneID_Toxo), total = n())
write.xlsx(host.int.local.markers.orth, "../Output/compScBabsSpecies/tables/bdiv_human_vs_all_spp_in_cow_phase_based_markers_TGGT1_orth_summarized.xlsx")


## bar plot 
host.markers.int.local.sig.stats <- host.markers.int.local.sig.specific %>% group_by(phase, ref.host) %>% 
  summarise(genes = list(GeneID), num.deg = n())

host.markers.int.local.sig.stats$phase <- factor(host.markers.int.local.sig.stats$phase, levels = c('G', 'SM', 'MC', 'C'))
host.markers.int.local.sig.stats$ref.host <- factor(host.markers.int.local.sig.stats$ref.host, levels = c('cow', 'human'))
colnames(host.markers.int.local.sig.stats) <- gsub("ref.host", "host", colnames(host.markers.int.local.sig.stats))


p <-  ggplot(host.markers.int.local.sig.stats, aes(x=num.deg, y=host, fill = host)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = c('cow' = 'darkgoldenrod', 'human' = 'darkolivegreen4'))+
  geom_text(aes(label=num.deg), hjust= 1.15, color="black", size=6, fontface = 'bold')+
  theme_bw()+
  facet_grid(phase~., scales = "free", space='free',switch = "y",
             labeller=label_wrap_gen(multi_line = TRUE))+
  theme(#axis.text.x = element_blank(),
    axis.text.x = element_text(face="bold", size=12, angle=0, hjust = 1),
    axis.title.y = element_text(size=14, face="bold"),
    axis.title.x = element_text(size=14, face="bold"),
    #axis.text.y = element_text(face="bold", size=12, angle=0),
    axis.text.y = element_blank())+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12,  angle = 0), 
        strip.placement = "outside")+
  theme(axis.line = element_line(color = 'black'), 
        axis.ticks = element_blank())+
  theme(legend.text = element_text(face = "bold", size = 10),
        legend.title = element_text(face = "bold", size = 14))+
  theme(strip.text.y.left = element_text(angle = 0))

p

ggsave(filename = "../Output/compScBabsSpecies/figs/Bdiv_human_vs_all_spp_in_cow_phase_based_markers_bar_plot.png", 
       plot = p, 
       width = 8, height = 4, 
       dpi = 300)


## top marker expression
host.markers.int.local.top <- host.markers.int.local.sig.specific %>% group_by( ref.host) %>% 
  slice_max(n = 1, order_by = avg_log2FC.x)
# 
# p <- FeaturePlot(object = S.O.integrated, features = host.markers.int.local.top$gene.x, label = F,repel = F,
#                  #shape.by  = 'spp',
#                  split.by = 'host',
#                  cols = c("grey", "blue"), reduction = "pca")
# 
# 
# plot(p)

exprs_hosts.phase <- exprExprData(S.O.integrated, host.markers.int.local.top$GeneID)
exprs_hosts.phase <- exprs_hosts.phase %>% mutate(host = ifelse(str_detect(Sample,"human"), "human", "cow"))


col_range <- colorRampPalette(c('white', 'blue'))
p <- ggplot(exprs_hosts.phase, aes(x=-PC_1,y=-PC_2)) +
  geom_point(aes(fill = expr,
                 #color = ""
  ), color = 'lightgray',
  shape=21, size = 1.5)+ 
  theme_bw(base_size = 14) +
  theme(legend.position = "right") +
  #scale_fill_gradientn(colours = viridis::inferno(10)) +
  scale_fill_gradientn(colours = col_range(10)) +
  #scale_fill_gradient(low = "gray66", high = "blue", midpoint = mid_point) + 
  #scale_fill_brewer(palette = "BuPu") +
  ylab('PC2') + xlab('PC1') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 16, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 16, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid"),
        strip.text = element_text(size = 16, face = "bold")) +
  facet_grid(GeneID~host) + 
  #ggtitle(titles[i]) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold"), 
    legend.title = element_text(size = 14, face = "bold"))+
  theme(
    #axis.text.x = element_blank(),
    axis.text.x = element_text(face="bold", size=16, angle=0, hjust = 1),
    axis.title.y = element_text(size=16, face="bold"),
    axis.title.x = element_text(size=16, face="bold"),
    #axis.text.y = element_text(face="bold", size=12, angle=0),
    #axis.text.y = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    #panel.border = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.y = element_blank()
  )

plot(p)


ggsave(filename="../Output/compScBabsSpecies/figs/expr_Bdiv_human_vs_all_spp_in_cow_top_1_phase_based_markers.png",
       plot=p,
       width = 10, height = 10,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

## please run sppDifferences.V3 whiich corrects for spp effect on host markers


#####

## this is not the right way of dooing it 
# ## species specific markers in cow
# 
# DefaultAssay(S.O.integrated) <- 'RNA'
# Idents(S.O.integrated) <- 'spp'
# S.O.spp.cow <- subset(S.O.integrated, ident = c('Bbig', 'Bbov', 'Bdiv_cow'))
# S.O.spp.cow@meta.data$host.spp <- paste(S.O.spp.cow@meta.data$host, S.O.spp.cow@meta.data$spp, sep = ":")
# unique(S.O.spp.cow$host.spp)
# 
# 
# DefaultAssay(S.O.spp.cow) <- 'RNA'
# Idents(S.O.spp.cow) <- 'host.spp'
# 
# makeCowSppMatchedContrasts <- function(S.O.integrated){
#   
#   objs <- as.character(unique(S.O.integrated@meta.data$host.spp))
#   
#   contrasts <- data.frame(ref = objs, dummy = 1) %>% full_join( data.frame(query = objs, dummy = 1), by = 'dummy') %>% 
#     mutate(ref.host = gsub(":.*", "", ref), ref.spp = gsub(".*:", "", ref), 
#            query.host = gsub(":.*", "", query), query.spp = gsub(".*:", "", query))
#   my.contrasts <- contrasts %>% filter(ref.spp != query.spp)
#   
#   return(my.contrasts)
#   
# }
# 
# contrasts <- makeCowSppMatchedContrasts(S.O.spp.cow)
# contrasts.groups <- contrasts %>% group_by(ref) %>% summarise(query = list(query))
# contrasts.groups$ref.spp <- gsub('.*:', '', contrasts.groups$ref)
# contrasts.groups$ref.host <- gsub(':.*', '', contrasts.groups$ref)
# 
# 
# cow.spp.DEGs <- mclapply(1:nrow(contrasts.groups), function(i){
#   tmp <- FindMarkers(S.O.spp.cow, ident.1 = contrasts.groups$ref[i],
#                      ident.2 = c(unlist(contrasts.groups$query[i])), 
#                      only.pos = T, verbose = T)
#   tmp$ref <- contrasts.groups$ref[i]
#   tmp$ref.spp <- contrasts.groups$ref.spp[i]
#   tmp$ref.host <- contrasts.groups$ref.host[i]
#   tmp$gene <- rownames(tmp)
#   tmp$GeneID <- gsub('-', '_', tmp$gene)
#   return(tmp)
# })
# cow.spp.DEGs.all <- bind_rows(cow.spp.DEGs)
# 
# cow.spp.DEGs.sig <- cow.spp.DEGs.all %>% dplyr::filter(avg_log2FC > 0.58 & p_val_adj < 0.01)
# #cow.spp.DEGs.sig  <- cow.spp.DEGs.sig %>% mutate(pct.ratio = (pct.1 / pct.2))
# 
# cow.spp.DEGs.sig.cc <- inner_join(cow.spp.DEGs.sig, cell.cycle.markers.sig, 
#                                   by = c('GeneID','ref.spp' = 'spp'))
# cow.spp.DEGs.sig.cc <- cow.spp.DEGs.sig.cc %>% dplyr::select(!contains(".y")) %>% distinct(GeneID, .keep_all  = T)
# comm <-  host.markers.int.local.sig.specific$GeneID[unique(host.markers.int.local.sig.specific$GeneID) %in% cow.spp.DEGs.sig.cc$GeneID]
# 
# gene_ids <- c("Bdiv_010630", "Bdiv_020470c", "Bdiiv_035020", "Bdiv_003910", 
#               "Bdiv_005170", "Bdiv_009600c", "Bdiv_018340", "Bdiv_022710", "Bdiv_000950", "Bdiv_037500")
# 
# gene_ids %in% comm
# 
# # # unique to spp and phase
# # host.markers.int.local.sig.specific <- host.markers.int.local.sig.cc %>% distinct(GeneID, .keep_all = T) %>%
# #   arrange(GeneID, desc(pct.ratio)) %>% group_by(GeneID) %>% slice_max(n = 1, order_by = pct.ratio)
# 
# 
# write.xlsx(host.markers.int.local.sig.specific, '../Output/compScBabsSpecies/tables/bdiv_human_vs_all_spp_in_cow_phase_based_markers_sig.xlsx')
# saveRDS(host.markers.int.local.sig.specific, '../Input/compScBabesia/RData/bdiv_human_vs_all_spp_in_cow_phase_based_markers_sig.RData')


#########################################################

## GO plot 

## GO trm 
in.dir <- '../Output/compScBabsSpecies/GO/GOToxoDB_Cow_vs_Human/'
#in.dir <- '../Output/compScBabsSpecies/GO/GOToxoDB_Species_Specific/'
all.files <- list.files(in.dir)



all.clust.items <- list()
for(f in all.files){
  nn <- gsub('\\.tsv', '', f)
  GF <- strsplit(nn, split = '_')[[1]][1]
  phase <- strsplit(nn, split = '_')[[1]][2]
  tmp <- read_tsv(paste(in.dir, f, sep = ''))
  tmp$GF <- GF
  tmp$phase <- phase
  all.clust.items <- c(all.clust.items, list(tmp))
}

all.clust.items <- do.call(rbind, all.clust.items)
names(all.clust.items) <- gsub(" ", ".", names(all.clust.items))

#saveRDS(all.clust.items, '../Input/compScBabesia/RData/GO_Conserved_Markers_TGGT1.RData')

filtered.Go <- all.clust.items %>% arrange(phase, Benjamini) %>% distinct() %>%
  group_by(phase) %>% mutate(rank = row_number()) %>%
  dplyr::filter(Benjamini < 0.1 & rank < 30 ) %>% 
  arrange(phase, Benjamini) 

#top <- filtered.Go %>% filter(phase == "HumG") %>% slice_min(n = 10, order_by = Benjamini)
# 

filtered.Go.lng <- filtered.Go %>% 
  mutate(Result.gene.list = strsplit(as.character(Result.gene.list), ",")) %>% 
  unnest(Result.gene.list)


Bdiv_TGGT1_orth <-  read.xlsx("../Input/compScBabesia/Orthologs/Bdiv_GT1_orth.xlsx")


GO_tab_orth <- inner_join(Bdiv_TGGT1_orth, filtered.Go.lng , by = c("TGGT1" = "Result.gene.list"))
GO_tab_df <- GO_tab_orth %>% group_by(ID, Name, GF, phase) %>% summarise(Bdiv_ID = list(Bdiv), Toxo_ID = list(TGGT1), total = n())
Go.Bdiv.TG1.IDs <- left_join(GO_tab_df, filtered.Go, by = c("ID", "Name", "GF", "phase"))

write.xlsx(Go.Bdiv.TG1.IDs, "../Output/compScBabsSpecies/GO/GO_Term_bdiv_human_vs_all_spp_in_cow_BdivIDs.xlsx")



# 
# p <- ggplot(filtered.Go, aes(x = phase, y = Name)) + 
#   geom_point(aes(colour = "red", size = -log(Benjamini))) +
#   theme_bw(base_size = 14) +
#   #scale_colour_gradient(limits=c(0, 0.01), low="red") +
#   ylab(NULL) + xlab(NULL) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face="bold")) + 
#   theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 14, face="bold")) +
#   theme(legend.position="none") +
#   theme(strip.background = element_rect(colour="black", fill="white", 
#                                         size=0.5, linetype="solid")) + guides(color = FALSE)+
#   ggtitle(strsplit(in.dir,split = "/")[[1]][5]) +
#   theme(plot.title = element_text(size = 10))
# 
# plot(p)
# 
# ggsave(filename="../Output/compScBabsSpecies/figs/GO_Human_vs_Cow.png", 
#        plot=p,
#        width = 12, height = 12, 
#        units = "in", # other options are "in", "cm", "mm" 
#        dpi = 300
# )


