library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(RColorBrewer)
library(parallel)

source('./util_funcs.R')

## This script reads the expr.csv files for all Babesia species
## It generates and saves the Seurat object. 

mito.genes <- read.xlsx("../Input/compScBabesia/genes/mitochondiral_genes.xlsx")
input.dir <- "../Input/compScBabesia/scRNAseqBabesia/"
count.files <- list.files(input.dir)
num.total.files <- length(count.files)

file.info <- data.frame(Species = rep("Babesia", num.total.files))

## These were determined with visual inspection
bdiv.human.cutoff <- c(200, 1200)
bdiv.cow.cutoff <- c(200, 1300)
bbig.cutoff <- c(80, 1100)
bbov.cutoff <- c(100, 1300)
bmic.cutoff <- c(50, 900)

bdiv.human.mt.cutoff <- 0.5
bdiv.cow.mt.cutoff <- 0.5
bbig.mt.cutoff <- 0.6
bbov.mt.cutoff <- 2.2
bmic.mt.cutoff <- 0.7

filter.mit <- F

processCountBabs <- function(input.dir, filename, down.sample = T){
  file.dir <- paste(input.dir, filename, sep = '')
  file.csv <- paste(filename, ".expr.csv", sep = "")
  
  print(file.csv)
  
  file.counts <- read.csv(paste(file.dir, file.csv, sep = '/'))
  genes <- file.counts$X
  expr <- file.counts[,-1]
  rownames(expr) <- genes

  
  spp <- filename
  print(spp)
  if(grepl('BBIG', toupper(filename))){
    cutoff <- bbig.cutoff
    mt.genes <- mito.genes %>% dplyr::filter(Spp == 'BBIG') %>% select(GeneID)
    mt.genes <- gsub('_', '-', mt.genes$GeneID)
    mt.cutoff <- bbig.mt.cutoff
  }else if(grepl('BBOV', toupper(filename))){
    cutoff <- bbov.cutoff
    mt.genes <- mito.genes %>% dplyr::filter(Spp == 'BBOV') %>% select(GeneID)
    mt.genes <- gsub('_', '-', mt.genes$GeneID)
    mt.cutoff <- bbov.mt.cutoff
  }else if(grepl('HUMAN', toupper(filename))){
    cutoff <- bdiv.human.cutoff
    mt.genes <- mito.genes %>% dplyr::filter(Spp == 'BDIV') %>% select(GeneID)
    mt.genes <- gsub('_', '-', mt.genes$GeneID)[c(1,2,3,5,6)] ## 4th one DNE
    mt.cutoff <- bdiv.human.mt.cutoff
  }else if(grepl('COW', toupper(filename))){
    cutoff <- bdiv.cow.cutoff
    mt.genes <- mito.genes %>% dplyr::filter(Spp == 'BDIV') %>% select(GeneID)
    mt.genes <- gsub('_', '-', mt.genes$GeneID)[c(1,2,3,5,6)] ## 4th one DNE
    mt.cutoff <- bdiv.cow.mt.cutoff
  }else if(grepl('BMIC', toupper(filename))){
    cutoff <- bmic.cutoff
    mt.genes <- mito.genes %>% dplyr::filter(Spp == 'BMIC') %>% select(GeneID)
    mt.genes <- gsub('_', '-', mt.genes$GeneID)[c(1,4,5,8,9,10)]
    mt.cutoff <- bmic.mt.cutoff
  }
  
  
  print(cutoff)
  print(mt.cutoff)
   
  S.O <- CreateSeuratObject(counts = expr, min.cells = 10, min.features = 100)
  S.O[["percent.mt"]] <- PercentageFeatureSet(object = S.O, features = mt.genes)
  
  if(filter.mit){
    S.O <- subset(S.O, subset = nFeature_RNA > cutoff[1] & nFeature_RNA < cutoff[2] & percent.mt < mt.cutoff)
  }else{
    S.O <- subset(S.O, subset = nFeature_RNA > cutoff[1] & nFeature_RNA < cutoff[2])
  }
  
  
  #VlnPlot(S.O, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  #FeatureScatter(S.O, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  #cutoffs <- quantile(S.O$nCount_RNA, probs = c(0.01, 0.9))
  #print(cutoffs)
  #S.O <- subset(S.O, subset = nFeature_RNA > cutoffs[1] & nFeature_RNA < cutoffs[2] )
  
 
  S.O <- prep_S.O(S.O)
  

  cluster <- as.character(S.O$seurat_clusters)
  
  pheno <- data.frame(Sample = names(S.O$orig.ident))
  pheno$spp <- spp
  pheno$cluster <- cluster
  pheno$cell <- paste(pheno$spp, pheno$cluster , sep = "")
  pheno$NAME <- paste(pheno$spp, pheno$Sample, sep = "_")
  
  ind <- match(pheno$Sample,  colnames(expr))
  colnames(expr) <- pheno$NAME[ind]
  
  if(down.sample){
    set.seed(100)
    S.O <- subset(x = S.O, downsample = 800)
  }
  
  S.O.obj <- paste(paste("S.O", spp, sep = "."), "RData", sep = ".")
  S.O.dir <- paste('../Input/compScBabesia/RData/' , S.O.obj, sep = "")
  saveRDS(S.O, S.O.dir)
  
  L <- list(pheno = pheno, S.O = S.O)
  return(L)
  
  return(L)
}

S.O.list <- list()

for(i in 1:num.total.files){
  file.info$filename[i] <- count.files[i]
  cat(paste('processing file', file.info$filename[i]))
  cat('\n')
  L <- processCountBabs(input.dir, file.info$filename[i])
  S.O.list <- c(S.O.list, list(L))
}

names(S.O.list) <- c('bbig', 'bbov', 'bdiv_cow', 'bdiv_human', 'bmic')
saveRDS(S.O.list, '../Input/compScBabesia/RData/S.O.list.Rdata')


### Working with orthologous genes only
orthologs <- read.xlsx("../Input/compScBabesia/Orthologs/Bdiv_Bmic_Bbov_Bbig_orth.xlsx")
#orthologs <- read.xlsx("../Input/compScBabesia/Orthologs/Bdiv_Bbig_Bbov_orth.xlsx")


num.objs <- length(S.O.list)

phenos <- lapply(S.O.list, `[[`, 1)
S.Os <-  lapply(S.O.list, `[[`, 2)

spps <- unlist(lapply(phenos, function(x) x$spp[1]))


ortho.genes.bdiv.id <- lapply(1:num.objs, function(i){
  pheno <- phenos[[i]]
  spp <- gsub('_.*', '', pheno$spp[1])
  col.ind <- match(toupper(spp), toupper(colnames(orthologs)))
  ind.orth <- which(gsub('-', '_', rownames(S.Os[[i]]@assays$RNA@data)) %in% orthologs[, col.ind])
  genes.with.orth <- gsub('-', '_', rownames(S.Os[[i]]@assays$RNA@data))[ind.orth]
  genes.bdiv.id <- orthologs$Bdiv[match(genes.with.orth,orthologs[, col.ind])]
} )

common.genes <- Reduce(intersect, ortho.genes.bdiv.id)
orthologs.common <- orthologs[which(orthologs$Bdiv %in% common.genes), ]


processCountBabsBdivOrth <- function(input.dir, filename, orthologs.common, down.sample = T){
  file.dir <- paste(input.dir, filename, sep = '')
  file.csv <- paste(filename, ".expr.csv", sep = "")
  
  file.counts <- read.csv(paste(file.dir, file.csv, sep = '/'))
  genes <- file.counts$X
  expr <- file.counts[,-1]
  rownames(expr) <- genes
  
  expr.mat <- as.matrix(S.Os[[i]]@assays$RNA@data)
  
  spp <- file.info$filename[i]
  
  col.ind <- match(toupper(gsub('_.*', '', spp)), toupper(colnames(orthologs.common)))
  expr <- expr[which(gsub('-', '_',rownames(expr)) %in% orthologs.common[, col.ind]),]
  rownames(expr) <- orthologs.common$Bdiv[match(gsub('-', '_',rownames(expr)), orthologs.common[, col.ind])]
  
  
  mt.genes <- mito.genes %>% dplyr::filter(Spp == 'BDIV') %>% select(GeneID)
  mt.genes <- gsub('_', '-', mt.genes$GeneID)[c(1,2,3,5,6)] ## 4th one DNE
  
  spp <- filename
  print(spp)
 
  if(grepl('BBIG', toupper(filename))){
    cutoff <- bbig.cutoff
    mt.cutoff <- bbig.mt.cutoff
  }else if(grepl('BBOV', toupper(filename))){
    cutoff <- bbov.cutoff
    mt.cutoff <- bbov.mt.cutoff
  }else if(grepl('HUMAN', toupper(filename))){
    cutoff <- bdiv.human.cutoff
    mt.cutoff <- bdiv.human.mt.cutoff
  }else if(grepl('COW', toupper(filename))){
    cutoff <- bdiv.cow.cutoff
    mt.cutoff <- bdiv.cow.mt.cutoff
  }else if(grepl('BMIC', toupper(filename))){
    cutoff <- bmic.cutoff
    mt.cutoff <- bmic.mt.cutoff
  }
  
  
  print(cutoff)
  
  S.O <- CreateSeuratObject(counts = expr, min.cells = 10, min.features = 100)
  S.O[["percent.mt"]] <- PercentageFeatureSet(object = S.O, features = mt.genes)
  
  if(filter.mit){
    S.O <- subset(S.O, subset = nFeature_RNA > cutoff[1] & nFeature_RNA < cutoff[2] & percent.mt < mt.cutoff)
  }else{
    S.O <- subset(S.O, subset = nFeature_RNA > cutoff[1] & nFeature_RNA < cutoff[2])
  }
  
  #FeatureScatter(S.O, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  #cutoffs <- quantile(S.O$nCount_RNA, probs = c(0.01, 0.9))
  #print(cutoffs)
  #S.O <- subset(S.O, subset = nFeature_RNA > cutoffs[1] & nFeature_RNA < cutoffs[2] )
  

  S.O <- prep_S.O(S.O)
  
  
  cluster <- as.character(S.O$seurat_clusters)
  
  pheno <- data.frame(Sample = names(S.O$orig.ident))
  pheno$spp <- spp
  pheno$cluster <- cluster
  pheno$cell <- paste(pheno$spp, pheno$cluster , sep = "")
  pheno$NAME <- paste(pheno$spp, pheno$Sample, sep = "_")
  
  ind <- match(pheno$Sample,  colnames(expr))
  colnames(expr) <- pheno$NAME[ind]
  
  if(down.sample){
    set.seed(100)
    S.O <- subset(x = S.O, downsample = 800)
  }
  
  S.O.obj <- paste(paste("S.O", spp, "ortholog", sep = "."), "RData", sep = ".")
  S.O.dir <- paste('../Input/compScBabesia/RData/' , S.O.obj, sep = "")
  saveRDS(S.O, S.O.dir)
  
  L <- list(pheno = pheno, S.O = S.O)
  return(L)
}

S.O.list <- list()

for(i in 1:num.total.files){
  file.info$filename[i] <- count.files[i]
  cat(paste('processing file', file.info$filename[i]))
  cat('\n')
  L <- processCountBabsBdivOrth(input.dir, file.info$filename[i], orthologs.common)
  S.O.list <- c(S.O.list, list(L))
}

names(S.O.list) <- c('bbig', 'bbov', 'bdiv_cow', 'bdiv_human', 'bmic')
saveRDS(S.O.list, '../Input/compScBabesia/RData/S.O.list.ortholog.Rdata')


