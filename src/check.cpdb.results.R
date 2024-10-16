#' ---
#' title: ""
#' output: 
#'  html_document:
#'    code_folding: hide
#' ---

#' ## DEG analysis 

#R4.0
rm(list=ls())

condaENV <- "/home/chenzh/miniconda3/envs/R4.0"
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

load("tmp_data/gene.meta.Rdata",verbose=T)
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/whole.IT.UMAP.coord.rds"))
data.ob.umap <- data.ob.umap %>% mutate(cluster=EML) 
#' overlap of DEGs
load(paste0("tmp_data/",TD,"/DEG.cluster.ov.Rdata"),verbose = T)

cpdb.list <- list()
#' loading cell-cell interaction results
for (sa in c("Sample_C","Sample_H","Sample_M")) {
  cpdb.list[[sa]] <- list()
  cpdb.list[[sa]]$pv <-  read.delim(paste0("tmp_data/",TD,"/",sa,"_cpdb/",sa,"/pvalues.txt"),head=T,stringsAsFactors = F) %>% tbl_df() %>% select(-c(interacting_pair:is_integrin)) 
  cpdb.list[[sa]]$mean_value <- read.delim(paste0("tmp_data/",TD,"/",sa,"_cpdb/",sa,"/means.txt"),head=T,stringsAsFactors = F) %>% tbl_df() 
  cpdb.list[[sa]]$decon <- read.delim(paste0("tmp_data/",TD,"/",sa,"_cpdb/",sa,"/deconvoluted.txt"),head=T,stringsAsFactors = F) %>% tbl_df() 
  temp <-  cpdb.list[[sa]]$pv%>% gather(cluster_pair,pvalue,-id_cp_interaction) %>% separate(cluster_pair,c("cluster_pair1","cluster_pair2"),sep="\\.",remove=F) %>% filter(pvalue < 0.05) 
  cpdb.list[[sa]]$sig_id_cp <- temp %>% pull(id_cp_interaction)
}

pair="Spp1_hi_Macro.Fibro2"
temp.list <- list()
for (pair in setdiff(colnames(cpdb.list$Sample_M$pv),"id_cp_interaction")) {
  temp.pv <- cpdb.list$Sample_C$pv %>% select(id_cp_interaction,pair) %>% mutate(sample="Sample_C") %>% bind_rows(cpdb.list$Sample_H$pv %>% select(id_cp_interaction,pair) %>% mutate(sample="Sample_H")) %>% bind_rows(cpdb.list$Sample_M$pv %>% select(id_cp_interaction,pair) %>% mutate(sample="Sample_M")) %>% setNames(c("id_cp_interaction","pvalue","sample"))
  temp.mean <- cpdb.list$Sample_C$mean_value %>% select(id_cp_interaction,pair) %>% mutate(sample="Sample_C") %>% bind_rows(cpdb.list$Sample_H$mean_value %>% select(id_cp_interaction,pair) %>% mutate(sample="Sample_H")) %>% bind_rows(cpdb.list$Sample_M$mean_value %>% select(id_cp_interaction,pair) %>% mutate(sample="Sample_M")) %>% setNames(c("id_cp_interaction","mean_value","sample"))
  temp.decon <-  cpdb.list$Sample_C$decon %>% select(gene_name,id_cp_interaction,complex_name) %>% group_by(id_cp_interaction,complex_name) %>% summarise(gene=paste(gene_name,collapse=":"))%>% bind_rows(cpdb.list$Sample_H$decon %>% select(gene_name,id_cp_interaction,complex_name) %>% group_by(id_cp_interaction,complex_name) %>% summarise(gene=paste(gene_name,collapse=":")))%>% bind_rows(cpdb.list$Sample_M$decon %>% select(gene_name,id_cp_interaction,complex_name) %>% group_by(id_cp_interaction,complex_name) %>% summarise(gene=paste(gene_name,collapse=":"))) %>% unique()
  temp.sel_id <-  temp.pv  %>% filter(pvalue < 0.05) %>% pull(id_cp_interaction) %>% unique()
  temp.out <- temp.pv %>% filter(id_cp_interaction %in% temp.sel_id)%>% left_join(temp.mean,by=c("id_cp_interaction","sample")) %>% left_join(temp.decon  %>% unique(),by="id_cp_interaction")
  temp.list[[pair]] <- temp.out %>% ggplot+geom_point(mapping=aes(x=sample,y=gene,col=(-1*log10(pvalue)),size=log(mean_value)))+ggtitle(pair)+ theme(plot.title = element_text(hjust=0.5))
}
temp.pv <- read.delim(paste0("tmp_data/",TD,"/",sa,"_cpdb/",sa,"/pvalues.txt"),head=T,stringsAsFactors = F) %>% tbl_df() %>% select(-c(interacting_pair:is_integrin)) %>% gather(cluster_pair,pvalue,-id_cp_interaction) %>% separate(cluster_pair,c("cluster_pair1","cluster_pair2"),sep="\\.",remove=F) %>% filter(pvalue < 0.05)
temp.mean <- read.delim(paste0("tmp_data/",TD,"/",sa,"_cpdb/",sa,"/means.txt"),head=T,stringsAsFactors = F) %>% tbl_df() %>% filter(id_cp_interaction %in% temp.pv$id_cp_interaction)%>% select(-c(interacting_pair:is_integrin)) %>% gather(cluster_pair,mean_exp,-id_cp_interaction) %>% separate(cluster_pair,c("cluster_pair1","cluster_pair2"),sep="\\.",remove=F)
temp.decon <- read.delim(paste0("tmp_data/",TD,"/",sa,"_cpdb/",sa,"/deconvoluted.txt"),head=T,stringsAsFactors = F) %>% tbl_df() %>% filter(id_cp_interaction %in% temp.pv$id_cp_interaction)

temp <- temp.pv %>% filter(cluster_pair1 %in% c("Spp1_hi_Macro","Fibro2")) %>% filter(cluster_pair2 %in% c("Spp1_hi_Macro","Fibro2")) %>% filter(cluster_pair1!=cluster_pair2) %>% left_join(temp.decon %>% select(gene_name,complex_name,id_cp_interaction) %>% unique(),by="id_cp_interaction")
