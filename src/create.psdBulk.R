#' ---
#' title: "create psd bulk dataset"
#' output: 
#'  html_document:
#'    code_folding: hide
#' ---

#' ## DEG analysis 

#R3.6
rm(list=ls())

condaENV <- "/home/chenzh/miniconda3/envs/R3.6"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(pheatmap))
suppressMessages(library(Seurat))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)
suppressMessages(library(ggvenn))
suppressMessages(library(DESeq))

# working directory
DIR <- "/home/chenzh/My_project/SHI_Glu"
setwd(DIR)

#' Loading R functions
source("~/PC/R_code/functions.R")
source("~/PC/SnkM/SgCell.R")
#source("src/local.quick.function.R")

options(digits = 4)
options(future.globals.maxSize= 3901289600)
TD="Aug_2022"


#pt.cutoff <- 0.4
#' loading data
meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds"))
load("tmp_data/gene.meta.Rdata",verbose=T)
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/whole.IT.UMAP.coord.rds")) %>% rows_update(readRDS(paste0("tmp_data/",TD,"/Fibryo.IT.UMAP.coord.rds")) %>% select(cell,EML,cluster_EML),by="cell")



#sel.expG <- rownames(lognormExp.mBN)
data.ob.umap <- data.ob.umap %>% mutate(cluster=EML) 


counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
rownames(counts.filter) <- gsub("_","-",rownames(counts.filter))

cells.list <- data.ob.umap %>%mutate(SID=EML)%>% split(.$SID) %>% lapply(function(x){x$cell})
psd.merge.meta <- data.ob.umap %>%mutate(SID=EML) %>% select(SID,sample,EML) %>% unique()
#' get the psd-bulk normalization results
counts.merge <- list()
for (n in names(cells.list)) {
  temp.cells <- cells.list[[n]]
  counts.merge[[n]] <- as.data.frame(rowSums(counts.filter[,temp.cells]) ) %>% setNames(n)
}
counts.merge <- counts.merge %>% do.call("cbind",.)
temp.counts <- counts.merge[,colnames(counts.merge)%>% setdiff(NULL)]#,"JPF2019_Tsw_AMLC"
temp.sel.gene <- rownames(temp.counts)[rowSums(temp.counts >= 2) >0] ## at least expressed with 2 counts in one dataset
temp.counts <- temp.counts[temp.sel.gene,]
temp.cds=newCountDataSet(temp.counts,colnames(temp.counts) ) %>% estimateSizeFactors()
counts.merge.norm <- counts(temp.cds, normalized=TRUE)

saveRDS(counts.merge.norm ,paste0("tmp_data/",TD,"/psd.EML.counts.merge.norm.rds"))


#' split by samples
cells.list <- data.ob.umap %>%mutate(SID=paste(sample,EML,sep=":"))%>% split(.$SID) %>% lapply(function(x){x$cell})
psd.merge.meta <- data.ob.umap %>%mutate(SID=paste(sample,EML,sep=":")) %>% select(SID,sample,EML) %>% unique()
#' get the psd-bulk normalization results
counts.merge <- list()
for (n in names(cells.list)) {
  temp.cells <- cells.list[[n]]
  counts.merge[[n]] <- as.data.frame(rowSums(counts.filter[,temp.cells]) ) %>% setNames(n)
}
counts.merge <- counts.merge %>% do.call("cbind",.)
temp.counts <- counts.merge[,colnames(counts.merge)%>% setdiff(NULL)]#,"JPF2019_Tsw_AMLC"
temp.sel.gene <- rownames(temp.counts)[rowSums(temp.counts >= 2) >0] ## at least expressed with 2 counts in one dataset
temp.counts <- temp.counts[temp.sel.gene,]
temp.cds=newCountDataSet(temp.counts,colnames(temp.counts) ) %>% estimateSizeFactors()
counts.merge.norm <- counts(temp.cds, normalized=TRUE)

saveRDS(counts.merge.norm ,paste0("tmp_data/",TD,"/psd.counts.merge.norm.rds"))





#' check macro cells only
cells.list <- data.ob.umap %>% filter(EML %in% c("M0Macro","M1Macro","M2Macro","Spp1_hi_Macro")) %>%mutate(SID=paste(sample,EML,sep=":"))%>% split(.$SID) %>% lapply(function(x){x$cell})
psd.merge.meta <- data.ob.umap %>%mutate(SID=paste(sample,EML,sep=":")) %>% select(SID,sample,EML) %>% unique()
#' get the psd-bulk normalization results
counts.merge <- list()
for (n in names(cells.list)) {
  temp.cells <- cells.list[[n]]
  counts.merge[[n]] <- as.data.frame(rowSums(counts.filter[,temp.cells]) ) %>% setNames(n)
}
counts.merge <- counts.merge %>% do.call("cbind",.)
temp.counts <- counts.merge[,colnames(counts.merge)%>% setdiff(NULL)]#,"JPF2019_Tsw_AMLC"
temp.sel.gene <- rownames(temp.counts)[rowSums(temp.counts >= 3) >3] ## at least expressed with 2 counts in one dataset
temp.counts <- temp.counts[temp.sel.gene,]
temp.cds=newCountDataSet(temp.counts,colnames(temp.counts) ) %>% estimateSizeFactors()
counts.merge.norm <- counts(temp.cds, normalized=TRUE)


temp.pca <-  prcomp(log2(t(counts.merge.norm+0.001)), retx=TRUE)
c=round(100*summary(temp.pca)$importance[2,1],digits=2)
d=round(100*summary(temp.pca)$importance[2,2],digits=2)
temp.pca <- temp.pca$x[,1:2] %>% as.data.frame() %>% tibble::rownames_to_column("SID") %>% tbl_df %>% inner_join(psd.merge.meta ,by="SID")
print(
  temp.pca  %>% ggplot(mapping=aes(x=PC1,y=PC2,col=EML,shape=sample))+geom_point()+geom_text(mapping=aes(label=SID))+theme_classic()#+xlim(-150,250)+ylim(-75,150)+ggtitle("PCA plot of Savana cells")+FunTitle()
)
pcol=colorRampPalette(c("#0E7DBE","white","#B41D22"))(100)
PH=pheatmap(t(cor(log2(counts.merge.norm+0.001),method = "spearman")),border="white",scale="none",col=pcol,legend=T,display_numbers = F,number_format="%.3f",annotation_col=(psd.merge.meta  %>% tibble::column_to_rownames("SID")))
plot(PH$tree_col,hang=-1)


#' check Fibro cells only
cells.list <- data.ob.umap %>% filter(big_EML %in% c("Fibro")) %>% filter(EML!="Prolifer_Fibro")  %>%mutate(SID=paste(EML,sep=":"))%>% split(.$SID) %>% lapply(function(x){x$cell})
psd.merge.meta <- data.ob.umap %>%mutate(SID=paste(sample,EML,sep=":")) %>% select(SID,sample,EML) %>% unique()
#' get the psd-bulk normalization results
counts.merge <- list()
for (n in names(cells.list)) {
  temp.cells <- cells.list[[n]]
  counts.merge[[n]] <- as.data.frame(rowSums(counts.filter[,temp.cells]) ) %>% setNames(n)
}
counts.merge <- counts.merge %>% do.call("cbind",.)
temp.counts <- counts.merge[,colnames(counts.merge)%>% setdiff(NULL)]
temp.sel.gene <- rownames(temp.counts)[rowSums(temp.counts >= 1) >4] ## at least expressed with 2 counts in one dataset
temp.counts <- temp.counts[temp.sel.gene,]
temp.cds=newCountDataSet(temp.counts,colnames(temp.counts) ) %>% estimateSizeFactors()
counts.merge.norm <- counts(temp.cds, normalized=TRUE)


temp.pca <-  prcomp(log2(t(counts.merge.norm+0.001)), retx=TRUE)
c=round(100*summary(temp.pca)$importance[2,1],digits=2)
d=round(100*summary(temp.pca)$importance[2,2],digits=2)
temp.pca <- temp.pca$x[,1:2] %>% as.data.frame() %>% tibble::rownames_to_column("SID") %>% tbl_df #%>% inner_join(psd.merge.meta ,by="SID")
print(
  temp.pca  %>% ggplot(mapping=aes(x=PC1,y=PC2))+geom_point()+geom_text(mapping=aes(label=SID))+theme_classic()#+xlim(-150,250)+ylim(-75,150)+ggtitle("PCA plot of Savana cells")+FunTitle()
)
pcol=colorRampPalette(c("#0E7DBE","white","#B41D22"))(100)
PH=pheatmap(t(cor(log2(counts.merge.norm+0.001),method = "spearman")),border="white",scale="none",col=pcol,legend=T,display_numbers = F,number_format="%.3f")
plot(PH$tree_col,hang=-1)

