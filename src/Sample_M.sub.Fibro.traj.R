#' ---
#' title: "define the lineage traj of sample M sub Fibro clusters"
#' output: 
#'  html_document:
#'    code_folding: hide
#'    
#' ---


rm(list=ls())
rewrite=FALSE

#' check whether in local computer
if (grepl("KI-",Sys.info()['nodename'])) {
  print("local computer")
  source("/Users/cheng.zhao/chzhao_bioinfo/PC/SnkM/SgCell.R")
  base_dir <- "/Users/cheng.zhao/Documents"
} else {
  print("On server")
  condaENV <- "/home/chenzh/miniconda3/envs/R4.0"
  LBpath <- paste0(condaENV ,"/lib/R/library")
  .libPaths(LBpath)
  base_dir="/home/chenzh"
}


suppressMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(tibble)
  #library(scran)
  #library(batchelor)
  library(SeuratWrappers)
  
  #library(clustree)
  library(monocle)
})



#' option 
set.seed(1)
DIR=paste0(base_dir,"/My_project/SHI_Glu")
setwd(DIR)
knitr::opts_knit$set(root.dir=DIR)
para_cal <- FALSE
rewrite=FALSE

#' loading local function 
source("~/PC/SnkM/SgCell.R")
#source("src/local.quick.fun.R")
filter <- dplyr::filter
rename <- dplyr::rename
select <- dplyr::select

if (para_cal) {
  suppressMessages(library(foreach))
  suppressMessages(library(doParallel))
  numCores <- 10
  registerDoParallel(numCores)
}


options(digits = 4)
options(future.globals.maxSize= 3901289600)
TD="Aug_2022"


# #### loading dataset
meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds"))
counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
data.all.ob.umap <- readRDS(paste0("tmp_data/",TD,"/whole.IT.UMAP.coord.rds")) %>% select(cell:phase_seurat,EML,cluster_EML,big_EML)   %>% rows_update(readRDS(paste0("tmp_data/",TD,"/Fibryo.IT.UMAP.coord.rds")) %>% select(cell,EML,cluster_EML),by="cell") 
data.all.ob.umap <- data.all.ob.umap %>% mutate(cluster=EML) 
TF.genes=data.table::fread("~/Genome_new/Rats/TF/Rattus_norvegicus_TF",stringsAsFactors = F,head=T) %>% tbl_df() %>% pull(Symbol) %>% intersect(rownames(counts.filter))


load("tmp_data/gene.meta.Rdata",verbose=T)
nGene=2000
npc=25


#' check all cells
#temp.M <- data.all.ob.umap %>% filter(sample=="Sample_M")
# temp.sel.expG <- rownames(counts.filter)[rowSums(counts.filter[,c(temp.M$cell)] >=1) >4]
# 
# data.ob <- CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) %>% FindVariableFeatures( selection.method = "vst", nfeatures = nGene, verbose = FALSE) %>% ScaleData(verbose=F,vars.to.regress=c("nGene"))%>% RunPCA(verbose=F) %>% RunUMAP(dims=1:npc,verbose=F) %>% FindNeighbors( dims = 1:npc,verbose = FALSE) %>%  FindClusters(resolution = 0.4,verbose = FALSE)
# cowplot::plot_grid(
#   DimPlot(data.ob,label=T)+theme_void()+NoLegend()+ggtitle("raw_cluster")+FunTitle(),
#   DimPlot(data.ob,label=T,group.by="sample")+theme_void()+NoLegend()+ggtitle("Sample")+FunTitle(),
#   DimPlot(data.ob,label=T,group.by="EML")+theme_void()+NoLegend()+ggtitle("Sample")+FunTitle(),
#   DimPlot(data.ob,label=F,group.by="EML")+theme_void()+NoLegend()+ggtitle("Sample")+FunTitle(),
#   #FeaturePlot(data.ob,"mt.perc"),
#   #FeaturePlot(data.ob,"nGene")
#   nrow=2,ncol=2
# )
# cowplot::plot_grid(plotlist=FunFP_plot(data.ob,c("Mki67","Top2a")))

#' check Fibro cells only
temp.M <- data.all.ob.umap %>% filter(sample=="Sample_M") %>% filter(big_EML=="Fibro") %>% filter(EML !="Prolifer_Fibro") # %>% filter(cell %in% (colnames(data.ob)[Idents(data.ob) %in% c(0,3,7)]))%>% left_join(readRDS(paste0("tmp_data/",TD,"/Fibryo.IT.UMAP.coord.rds")) %>% mutate(sub_EML=EML) %>% select(cell,sub_EML) ,by="cell")
temp.sel.expG <- rownames(counts.filter)[rowSums(counts.filter[,c(temp.M$cell)] >=1) >4]

data.ob <- CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) %>% FindVariableFeatures( selection.method = "vst", nfeatures = nGene, verbose = FALSE) %>% ScaleData(verbose=F,vars.to.regress="nGene")%>% RunPCA(verbose=F) %>% RunUMAP(dims=1:npc,verbose=F) %>% FindNeighbors( dims = 1:npc,verbose = FALSE) %>%  FindClusters(resolution = 0.7,verbose = FALSE)#,


cowplot::plot_grid(
  DimPlot(data.ob,label=T)+theme_void()+NoLegend()+ggtitle("raw_cluster")+FunTitle(),
  DimPlot(data.ob,label=T,group.by="sample")+theme_void()+NoLegend()+ggtitle("Sample")+FunTitle(),
  # DimPlot(data.ob,label=T,group.by="sub_EML")+theme_void()+NoLegend()+ggtitle("Sample")+FunTitle(),
  DimPlot(data.ob,label=T,group.by="EML")+theme_void()+NoLegend()+ggtitle("Sample")+FunTitle(),
  nrow=2,ncol=2
)

if (file.exists(paste0("tmp_data/",TD,"/Sample_M.subFibro.trajectory.dm.Rdata"))) {
  load(paste0("tmp_data/",TD,"/Sample_M.subFibro.trajectory.dm.Rdata"),verbose = T)
}else{
  #'data.temp <- data.ob
  data.temp <- subset(data.ob,cells=(temp.M %>% filter(EML %in% c("Fib_SC3","Fib_SC4","Fib_SC5","Fib_SC6","Fib_SC7"))) %>% pull(cell))
  #data.temp=data.ob
  Idents(data.temp) <- factor(data.temp@meta.data$EML)
  #temp.mk <- FunRF_FindAllMarkers_para(data.temp)
  
  #' check the trajecory analysis
  cds <- newimport(data.temp)
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)
  
  #1000 nGene not good , 500 seems good 
  cds.traj.gene <- data.ob %>% FindVariableFeatures( selection.method = "vst", nfeatures = 500, verbose = FALSE) %>% VariableFeatures()
  #cds.traj.gene <-  temp.mk$sig %>% filter(power > 0.4) %>% pull(gene)
  #cds.traj.gene <- temp.mk$detail %>% filter(power > 0.4)  %>% group_by(BG,EG) %>% top_n(50,power) %>% pull(gene) %>% unique()
  length( cds.traj.gene )
  cds.out <- setOrderingFilter(cds, cds.traj.gene  )
  cds.out <- reduceDimension(cds.out, max_components = 2, num_dim = npc,reduction_method = 'tSNE', verbose = T,residualModelFormulaStr = "~nGene")#,)
  cds.out<- clusterCells(cds.out, num_clusters = 6)
  cds.out <- reduceDimension(cds.out, max_components = 2,, num_dim = npc,reduction_method = 'DDRTree', verbose = F,residualModelFormulaStr = "~nGene")# work with 500
  cds.out <- orderCells(cds.out)
  cds.out <- orderCells(cds.out,root_state=9)
  
  
  #' get the genes changed along the pseudo-time 
  diff_TF_pse <- differentialGeneTest(cds.out[TF.genes %>% intersect(rownames(cds.out)),], fullModelFormulaStr="~sm.ns(Pseudotime)")%>% tbl_df() %>% filter(qval < 0.05) %>% arrange(qval)
  diff_TF_pse.EML <- differentialGeneTest(cds.out[TF.genes %>% intersect(rownames(cds.out)),], fullModelFormulaStr="~EML")%>% tbl_df() %>% filter(qval < 0.05) %>% arrange(qval)
  
  cds.out.DM <- t(cds.out@reducedDimS) %>% as.data.frame()%>% tbl_df() %>% mutate(cell=colnames(cds.out)) %>% rename(traj_Dim1=V1,traj_Dim2=V2)%>% left_join(pData(cds.out) %>% tibble::rownames_to_column("cell") %>% mutate(monocle_Cluster=Cluster) %>% select(cell,Size_Factor,Pseudotime,State,monocle_Cluster) %>% mutate(State=paste0("St",as.vector(State)),monocle_Cluster=paste0("MC",monocle_Cluster)) %>% tbl_df() ,by="cell")   %>% left_join(temp.M  ,by="cell")
  temp.traj.plot <- plot_cell_trajectory(cds.out, show_branch_points=F,cell_size=0.75)+NoAxes()+NoLegend()
  branch.layer <- temp.traj.plot$layer[[1]]
  
  #' save object
  saveRDS(cds.out.DM, file=paste0("tmp_data/",TD,"/Sample_M.subFibro.trajectory.dm.rds"))
  save(cds,  cds.out, cds.traj.gene, cds.out.DM,diff_TF_pse,diff_TF_pse.EML,file=paste0("tmp_data/",TD,"/Sample_M.subFibro.trajectory.dm.Rdata"))
  
}


plot_genes_in_pseudotime(cds.out[as.vector(head(diff_TF_pse$gene_short_name)),], color_by = "EML")
plot_genes_in_pseudotime(cds.out[as.vector(diff_TF_pse$gene_short_name[16:25]),], color_by = "EML",ncol = 3)
c("Cebpd","Cebpb","Fosb","Egr1")

pdf("tmp_data/fig_pdf/temp.test.pdf",6,4.5)
plot_genes_in_pseudotime(cds.out[c("Cebpd","Cebpb","Fosb","Egr1"),], color_by = "EML",ncol = 2)+scale_color_manual(values=EML.col[c("Fib_SC3","Fib_SC4","Fib_SC5","Fib_SC6","Fib_SC7")])
dev.off()

plot_genes_branched_heatmap(cds.out[as.vector(diff_TF_pse$gene_short_name),],branch_point = 1,num_clusters = 5,cores = 1,use_gene_short_name = T, show_rownames = T)

# BEAM_res <- BEAM(cds.out, branch_point = 1, cores = 10)
# BEAM_res <- BEAM_res[order(BEAM_res$qval),]
# BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
#plot_genes_branched_heatmap(cds.out[as.vector(head(BEAM_res$gene_short_name,50)),],branch_point = 1,num_clusters = 4,cores = 1,use_gene_short_name = T, show_rownames = T)
#plot_genes_branched_heatmap(cds.out[mk.list$all$sig %>% filter(NP=="pos") %>% group_by(set)  %>% top_n(15,power)%>% ungroup() %>% pull(gene),],branch_point = 1,num_clusters = 4,cores = 1,use_gene_short_name = T, show_rownames = T)

#' some check

cowplot::plot_grid(
  plot_cell_trajectory(cds.out, color_by = "State",show_branch_points=F)+NoAxes(),
  plot_cell_trajectory(cds.out, color_by = "EML")+NoAxes(),
  #plot_cell_trajectory(cds.out, color_by = "sub_EML")+NoAxes(),
  plot_cell_trajectory(cds.out, color_by = "Pseudotime")+NoAxes()
  
)
pdf("tmp_data/fig_pdf/sup.subM.traj.psd.pdf",4.5,4.5)
plot_cell_trajectory(cds.out, color_by = "Pseudotime",show_branch_points=F)+NoAxes()
dev.off()


EML.col <- c("Fib_SC1"="#4EA74A","Fib_SC2"="#8E4B99","Fib_SC3"="#E88F14","Fib_SC4"="#E57FB0","Fib_SC5"="#B699C6","Fib_SC6"="#A15428","Fib_SC7"="#55A7D6","Fib_SC8"="#202D70")

temp.traj.plot <- plot_cell_trajectory(cds.out, show_branch_points=F,cell_size=0.75)+NoAxes()+NoLegend()
branch.layer <- temp.traj.plot$layer[[1]]
temp.plot <- list()
for (n in unique(cds.out.DM$EML)) {
  temp.input <- cds.out.DM %>% mutate(sel=ifelse(EML==n , "sel","no")) %>% arrange(sel)
  temp.plot[[paste("traj.hl",n,sep=".")]] <-  ggplot()+branch.layer+geom_point(data=(temp.input %>% filter(sel=="no")),mapping=aes(x=traj_Dim1,y=traj_Dim2,col=selected),size=0.15,col="grey")+geom_point(data=(temp.input %>% filter(sel=="sel")),mapping=aes(x=traj_Dim1,y=traj_Dim2,col=selected),size=0.15,col=EML.col[n])+theme_classic()+NoAxes()+NoLegend()+ggtitle(n)+theme(plot.title = element_text(hjust=0.5))
}
pdf("tmp_data/fig_pdf/sup.subM.traj.hl.EML.pdf",8,5)
cowplot::plot_grid(plotlist=temp.plot[c(5,3,4,1,2)])
dev.off()

cds.out.DM <- cds.out.DM%>% mutate(path_branch=ifelse(State %in% c("St9"),"Path0",path_branch))  %>% mutate(path_branch=ifelse(State %in% c("St7","St8"),"Path1",path_branch)) %>% mutate(path_branch=ifelse(State %in% c("St1"),"Path2",path_branch))  %>% mutate(path_branch=ifelse(State %in% c("St2","St3","St4","St5","St6"),"Path3",path_branch))
print(
  # cds.out.DM  %>% group_by(EML,path_branch)  %>% summarise(nCell=n()) %>% group_by(EML) %>% mutate(prop=nCell/sum(nCell)) %>% group_by(path_branch) %>% mutate(prop=prop/sum(prop)) %>% ggplot+geom_bar(mapping=aes(x=path_branch,y=prop,fill=EML),stat="identity")
  cds.out.DM  %>% group_by(EML,path_branch)  %>% summarise(nCell=n())  %>% mutate(prop=nCell/sum(nCell)) %>% group_by(path_branch)%>% mutate(prop=nCell/sum(nCell)) %>% ggplot+geom_bar(mapping=aes(x=path_branch,y=prop,fill=EML),stat="identity")
)
print(
  cds.out.DM %>% ggplot()+geom_point(mapping=aes(x=traj_Dim1,y=traj_Dim2,col=path))+theme_classic()+NoAxes()+NoLegend()
)
