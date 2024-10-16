#' ---
#' title: "Defining the clusters"
#' output: 
#'  html_document:
#'    code_folding: hide
#' ---




# ### Loading R library

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
suppressMessages(library(clustree))

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


options(digits = 4)
options(future.globals.maxSize= 3901289600)
TD="Aug_2022"



# #### loading dataset
meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds"))
load("tmp_data/gene.meta.Rdata",verbose=T)



nGene=3000
npcs=40

gd.sel.mk <- list(fibroblast=c("Dpt","Dcn","Twist2"),sMCs=c("Des","Cox4i2","Tnn"),Adipocytes_Endothelial=c("Plvap","Rbp7","Fabp4"),Monocytes=c("Fcnb","Lyz2"),M1_Macrophages =c("Il1b"),NK_T=c("Ccl4"),Macrophages=c("C1qa","Spp1","Cd163"),M2a_Macrophages =c("Clec10a"),Erythroid=c("Hbb"),Neutrophils=c("G0s2","S100a9"),Tcell=c("Ccl5"),Keratinocytes_skins=c("Krt17","Ly6d","Dmkn"),T_reg=c("Stmn1","Ube2c"),prolifer=c("Top2a","Mki67","Stmn1"),Dendritic=c("RT1-Da","Cd74","Ccl17"),Endothelial=c("Apln","Krt15"),granulocyte=c("Cpa3","Cma1"),MAST=c("Tpsb2","Tpsab1"))


savefile <- paste0("tmp_data/",TD,"/whole.UMAP.IT.results.rds")
#saveRDS(data.ob,paste0("tmp_data/",TD,"/whole.UMAP.IT.results.nonGENE.rds"))
#' CCA for mouse dataset
if (file.exists(savefile)) {
  data.ob <- readRDS(savefile)
  data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/whole.IT.UMAP.coord.rds"))
}else{
  counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
  lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))
  rownames(lognormExp.mBN) <- gsub("_","-",rownames(lognormExp.mBN))
  
  
  # get the expressed genes ( expressed in at least 1 datatsets)
  expG.set <- list()
  for (b in unique(meta.filter$sample  %>% unique() %>% as.vector())) { 
    temp.cell <- meta.filter %>% filter(sample==b) %>% pull(cell)
    expG.set[[b]] <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]
  }
  sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()

  #
  temp.sel.expG <- sel.expG
  temp.M <- meta.filter 
  data.temp <-   CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) 
  data.list <- SplitObject(data.temp, split.by = "sample") %>% lapply( function(x) {SCTransform(x,,method = "glmGamPoi",verbose = FALSE,vars.to.regress = c("nGene"))}) #vars.to.regress = c("mt.perc"),#,vars.to.regress = c("nGene")
  #data.list <- SplitObject(data.temp, split.by = "sample") %>% lapply( function(x) {SCTransform(x,,method = "glmGamPoi",verbose = FALSE)}) #vars.to.regress = c("mt.perc"),#,vars.to.regress = c("nGene")
  
  anchor.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures =  nGene)
  
  #rm(counts.filter)
  data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = anchor.features ,verbose=F)
  data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT",anchor.features = anchor.features,verbose=F)
  
  data.ob <- IntegrateData(data.anchors, normalization.method = "SCT") %>%  RunPCA( verbose = FALSE,npcs=50)%>% RunUMAP(dims=1:npcs,verbose=F)   %>% FindNeighbors( dims = 1:npcs,verbose = FALSE)  #%>% FindNeighbors( dims = 1:npcs,verbose = FALSE,nn.method="annoy",annoy.metric="cosine") #45 not good
  
  data.ob@assays$RNA@data <- as.matrix(lognormExp.mBN[rownames(data.ob@assays$RNA@data),colnames(data.ob@assays$RNA@data)])
  
  # 50 ,regress out nGene , not working.
  
  # Annotating clusters (based on the WT cells only)
  if (file.exists(paste0("tmp_data/",TD,"/singleR.pred.rds"))) {
    wt.pred <- readRDS(paste0("tmp_data/",TD,"/singleR.pred.rds"))
  }else{
    suppressMessages(library(scran))
    temp.M <- meta.filter %>% filter(sample=="Sample_M") 
    temp.sel.expG <-expG.set$Sample_M
    temp.sce <-  SingleCellExperiment(list(counts=as.matrix(counts.filter[temp.sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors() 
    
    
    suppressMessages(library(SingleR))
    ref.data <- celldex::MouseRNAseqData()
    predictions <- SingleR(test=temp.sce, assay.type.test=1, ref=ref.data, labels=ref.data$label.main)
    wt.pred <- data.frame(cell=colnames(temp.sce),pred=predictions$labels) %>% tbl_df()
    saveRDS(wt.pred,file=paste0("tmp_data/",TD,"/singleR.pred.rds"))
  }
  data.ob@meta.data$pred <- NA
  data.ob@meta.data[wt.pred$cell,]$pred <- wt.pred$pred
  
  #' cell cycle from Seurat
  load("/home/chenzh/Genome_new/Mouse/seurat.mouse.phase.Rdata",verbose=T)
  temp <- CellCycleScoring(data.ob, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
  data.ob@meta.data$phase_seurat <- temp@meta.data$Phase
  
  data.temp <- data.ob  %>%  FindClusters(resolution = 0.3,verbose = FALSE)
  
  data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,nGene:sample,pred,phase_seurat)) %>% mutate(seurat_clusters=paste0("C",as.vector(Idents(data.temp))))  %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$pca@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:PC_40),by="cell") #%>% left_join(sn.result  %>% lapply( function(x) {return(x$sce.out)}) %>% do.call("bind_rows",.) %>% select(cell:AUC_tSNE2),by="cell")
  data.ob.umap <- data.ob.umap %>% mutate(EML=seurat_clusters)%>% mutate(EML=ifelse(seurat_clusters %in% c("C7"),"M2Macro",EML))%>% mutate(EML=ifelse(seurat_clusters %in% c("C1"),"Spp1_hi_Macro",EML)) %>% mutate(EML=ifelse(seurat_clusters %in% c("C2"),"M0Macro",EML))%>% mutate(EML=ifelse(seurat_clusters %in% c("C6"),"M1Macro",EML))%>% mutate(EML=ifelse(seurat_clusters %in% c("C13"),"Mait_Tcell",EML)) %>% mutate(EML=ifelse(seurat_clusters %in% c("C15"),"Dendritic",EML)) %>% mutate(EML=ifelse(seurat_clusters %in% c("C3"),"SMC",EML))%>% mutate(EML=ifelse(seurat_clusters %in% c("C16"),"Prolifer_SMC",EML))   %>% mutate(EML=ifelse(seurat_clusters %in% c("C0"),"Fibro2",EML))  %>% mutate(EML=ifelse(seurat_clusters %in% c("C4"),"Fibro1",EML))  %>% mutate(EML=ifelse(seurat_clusters %in% c("C14"),"Prolifer_Fibro",EML)) %>% mutate(EML=ifelse(seurat_clusters %in% c("C5"),"Endo",EML))  %>% mutate(EML=ifelse(seurat_clusters %in% c("C9"),"NKT",EML)) %>% mutate(EML=ifelse(seurat_clusters %in% c("C11"),"NK",EML))  %>% mutate(EML=ifelse(seurat_clusters %in% c("C10"),"Kera",EML)) %>% mutate(EML=ifelse(seurat_clusters %in% c("C12"),"MelScwann",EML)) %>% mutate(EML=ifelse(seurat_clusters %in% c("C17"),"Mast",EML))%>% mutate(EML=ifelse(seurat_clusters %in% c("C8"),"Neutrophil",EML)) %>% mutate(big_EML=EML) %>% mutate(big_EML=ifelse(EML %in% c("SMC","Prolifer_SMC"),"SMC",big_EML))  %>% mutate(big_EML=ifelse(EML %in% c("Fibro1","Fibro2","Prolifer_Fibro"),"Fibro",big_EML)) %>% mutate(big_EML=ifelse(EML %in% c("M0Macro","M1Macro","M2Macro","Mait_Tcell","Neutrophil","NK","NKT","Spp1_hi_Macro","Dendritic","Mast"),"Immun",big_EML))### Consult PMID: 35013299#
  data.ob.umap <- data.ob.umap %>% mutate(cluster_EML=EML)
  #Unk2 ->  Dendritic
  #Unk1 -> Mait T cells
  # Dendtritic -> M0 marcophage
  # Monocytes -> M1 macrophage
  #M1 macrophage -> Spp1_hi_Macrophage
  # data.ma.ob.umap <- data.ma.ob.umap %>% mutate(subCT=seurat_clusters) %>% mutate(subCT=ifelse(seurat_clusters %in% c("C3"),"M2Macro",subCT)) %>% mutate(subCT=ifelse(seurat_clusters %in% c("C1"),"Dendritic",subCT))%>% mutate(subCT=ifelse(seurat_clusters %in% c("C0"),"M1Macro",subCT))%>% mutate(subCT=ifelse(seurat_clusters %in% c("C2"),"Monocytes",subCT)) %>% mutate(subCT=ifelse(seurat_clusters %in% c("C5"),"Neutro_Macro",subCT)) %>% mutate(subCT=ifelse(seurat_clusters %in% c("C4"),"Unk3",subCT)) 
  #no Erytho
  
  # saving object
  saveRDS(data.ob,file=savefile)
  saveRDS(data.ob.umap,file=paste0("tmp_data/",TD,"/whole.IT.UMAP.coord.rds"))

}
data.temp <- data.ob  %>%  FindClusters(resolution = 0.3,verbose = FALSE)
data.temp@meta.data$EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"EML"]
data.temp@meta.data$big_EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"big_EML"]
cowplot::plot_grid(
  DimPlot(data.temp,label=T)+theme_void()+NoLegend()+ggtitle("raw_cluster")+FunTitle(),
  DimPlot(data.temp,label=T,group.by="sample")+theme_void()+NoLegend()+ggtitle("Sample")+FunTitle(),
  DimPlot(data.temp,label=T,group.by="EML")+theme_void()+NoLegend()+ggtitle("EML")+FunTitle(),
  DimPlot(data.temp,label=T,group.by="big_EML")+theme_void()+NoLegend()+ggtitle("big_EML")+FunTitle(),
  #FeaturePlot(data.ob,"mt.perc"),
  #FeaturePlot(data.ob,"nGene")
  nrow=2,ncol=2
)
DefaultAssay(data.temp) <- "RNA"
cowplot::plot_grid(plotlist=FunFP_plot(data.temp,gd.sel.mk$prolifer))


#  #### get the resolution
#data.ob  %>%  FindClusters(resolution = seq(0.3,0.8,by=0.1),verbose = FALSE) %>% clustree( prefix = "integrated_snn_res.") # 0.15 Or 0.1
#ac.cells <- colnames(data.ob)[(Idents(data.ob  %>%  FindClusters(resolution = 1.6,verbose = FALSE))==37)]
#data.temp <- data.ob %>% FindNeighbors( dims = 1:npcs,verbose = FALSE) %>%  FindClusters(resolution = 0.1,verbose = FALSE)
#sn.result <- readRDS("tmp_data/scenic.processs.rds")
#data.temp <- data.ob  %>%  FindClusters(resolution = 0.2,verbose = FALSE)
#' update data.ob
DefaultAssay(data.ob) <- "RNA"
Idents(data.ob) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[colnames(data.ob),"EML"])
data.ob@meta.data$EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"EML"]
data.ob@meta.data$cluster_EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"cluster_EML"]
data.ob@meta.data$big_EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"big_EML"]


print(
  DimPlot(data.ob,label=T,group.by="EML")+theme_void()+NoLegend()+ggtitle("EML")+FunTitle()
)
print(
  DimPlot(data.ob,label=T,split.by="sample")+theme_void()+NoLegend()+ggtitle("EML")+FunTitle()
)
# cowplot::plot_grid(
#   DimPlot(data.ob,label=T)+theme_void()+NoLegend()+ggtitle("raw_cluster")+FunTitle(),
#   DimPlot(data.ob,label=T,group.by="sample")+theme_void()+NoLegend()+ggtitle("Sample")+FunTitle(),
#   DimPlot(data.ob,label=T,group.by="big_EML")+theme_void()+NoLegend()+ggtitle("Sample")+FunTitle(),
#   #FeaturePlot(data.ob,"mt.perc"),
#   #FeaturePlot(data.ob,"nGene")
#   nrow=2,ncol=2
# )


# check cell proportation
cowplot::plot_grid(
  data.ob.umap %>% group_by(sample,cluster_EML) %>% summarise(nCell=n()) %>% group_by(sample) %>% mutate(prop=nCell/sum(nCell)) %>% ggplot+geom_bar(mapping=aes(x=cluster_EML,y=prop,fill=sample),stat="identity",position="dodge")+ theme(axis.text.x=element_text(angle = 90)),
  data.ob.umap %>% group_by(sample,big_EML) %>% summarise(nCell=n()) %>% group_by(sample) %>% mutate(prop=nCell/sum(nCell)) %>% ggplot+geom_bar(mapping=aes(x=big_EML,y=prop,fill=sample),stat="identity",position="dodge")+ theme(axis.text.x=element_text(angle = 90)),
  data.ob.umap %>% filter(big_EML=="Immun") %>% group_by(sample,cluster_EML) %>% summarise(nCell=n()) %>% group_by(sample) %>% mutate(prop=nCell/sum(nCell)) %>% ggplot+geom_bar(mapping=aes(x=cluster_EML,y=prop,fill=sample),stat="identity",position="dodge")+ theme(axis.text.x=element_text(angle = 90)),
  data.ob.umap %>% filter(big_EML=="Fibro") %>% group_by(sample,cluster_EML) %>% summarise(nCell=n()) %>% group_by(sample) %>% mutate(prop=nCell/sum(nCell)) %>% ggplot+geom_bar(mapping=aes(x=cluster_EML,y=prop,fill=sample),stat="identity",position="dodge")+ theme(axis.text.x=element_text(angle = 90))
)


#VlnPlot(data.ob,c("RGD1559482","Cd14","Ftl1"))
#  data.ob.umap %>% group_by(sample,EML) %>% summarise(nCell=n()) %>% group_by(sample) %>% mutate(prop=nCell/sum(nCell)) %>% ggplot+geom_bar(mapping=aes(x=EML,y=prop,fill=sample),stat="identity",position="dodge")+ theme(axis.text.x=element_text(angle = 90))
#  data.ma.ob.umap %>% group_by(sample,subCT) %>% summarise(nCell=n()) %>% group_by(sample) %>% mutate(prop=nCell/sum(nCell)) %>% ggplot+geom_bar(mapping=aes(x=subCT,y=prop,fill=sample),stat="identity",position="dodge")+ theme(axis.text.x=element_text(angle = 90))

cowplot::plot_grid(
  DimPlot(data.ob,label=T)+theme_void()+NoLegend()+ggtitle("raw_cluster")+FunTitle(),
  DimPlot(data.ob,label=T,group.by="sample")+theme_void()+NoLegend()+ggtitle("Sample")+FunTitle(),
  FeaturePlot(data.ob,"mt.perc"),
  FeaturePlot(data.ob,"nGene")
)

DimPlot(data.ob,label=T)+theme_void()+NoLegend()+ggtitle("raw_cluster")+FunTitle()
DimPlot(data.ob,label=T,split.by="sample")+theme_void()+NoLegend()+ggtitle("raw_cluster")+FunTitle()

#' get the lineage markers
if (file.exists(paste0("tmp_data/",TD,"/marker.rds"))) {
  total.mk <- readRDS(paste0("tmp_data/",TD,"/marker.rds"))
}else{
  total.mk <- list()
  total.mk$all <- FunRF_FindAllMarkers_para(data.ob)
  
  data.temp <- subset(data.ob,cells=(data.ob.umap %>% filter(!EML %in% c("Prolifer_SMC","Prolifer_Fibro")) %>% pull(cell)))
  data.temp=RenameIdents(data.temp,"Fibro1"="Fibro","Fibro2"="Fibro",M1Macro="big_monocytes",M2Macro="big_monocytes",Spp1_hi_Macro="big_monocytes",M0Macro="big_monocytes","Mait_Tcell"="big_monocytes")
  total.mk$merg_all <- FunRF_FindAllMarkers_para(data.temp)
  data.temp <- subset(data.ob,cells=(data.ob.umap %>% filter(EML %in% c("M1Macro","M2Macro","Spp1_hi_Macro","M0Macro")) %>% pull(cell)))#"M0Macro",
  total.mk$big_monocytes <- FunRF_FindAllMarkers_para(data.temp)
  data.temp <- subset(data.ob,cells=(data.ob.umap %>% filter(EML %in% c("Fibro1","Fibro2")) %>% pull(cell)))
  total.mk$Fibro <- FunRF_FindAllMarkers_para(data.temp)
  
  saveRDS(total.mk,file=paste0("tmp_data/",TD,"/marker.rds"))
}



#' check the marker gene expression
lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))
rownames(lognormExp.mBN) <- gsub("_","-",rownames(lognormExp.mBN))

#' #### check the marker genes expression
temp.DM.sel <-  c(total.mk$all$sig %>% group_by(set) %>% top_n(5,power) %>% group_by(gene) %>% top_n(1,power)  %>% split(.,.$set),(total.mk$merg_all$sig %>% group_by(set) %>% top_n(5,power) %>% group_by(gene) %>% top_n(1,power)  %>% split(.,.$set))[c("big_monocytes","Fibro")]) %>% lapply(function(x) {head(x,5)}) 
temp.DM.sel <- temp.DM.sel[c("Endo","Kera","Mast","MelScwann","SMC","Neutrophil","NK","NKT","Fibro1","Fibro2","Fibro","Dendritic","Mait_Tcell","M0Macro","M1Macro","M2Macro","Spp1_hi_Macro","big_monocytes")]
temp.M <-  data.ob.umap %>% select(cell,sample,EML) %>% split(.,.$EML) %>% lapply(function(x){FunMaSF(x,100)}) %>% do.call("bind_rows",.) %>% mutate(od=factor(EML,c(setdiff(names(temp.DM.sel),c("Fibro","big_monocytes")),"Prolifer_Fibro","Prolifer_SMC"),ordered = T)) %>% arrange(od,sample) %>% select(-EML) 
temp.exp <- lognormExp.mBN[unlist(lapply(temp.DM.sel,function(x){return(x$gene)})),temp.M$cell] %>% as.data.frame()#lognormExp.mBN
out.limit <- 2.5
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
temp.sel.exp[temp.sel.exp>out.limit] <- out.limit
temp.sel.exp[temp.sel.exp<  (-1*out.limit)] <- -1*out.limit
temp.anno <- temp.M %>% tibble::column_to_rownames("cell") 
#+ fig.width=15,fig.height=15
pheatmap::pheatmap(temp.sel.exp[,rownames(temp.anno)],cluster_rows=F,cluster_cols=F,scale="none",annotation_col=temp.anno,show_colnames=F,color=(colorRampPalette(c("royalblue3","white","firebrick4"))(50)), show_rownames = T,border_color="NA",main="Top markers",gaps_row = rep(unlist(lapply(temp.DM.sel,function(x){return(nrow(x))})) %>% cumsum(),each=2),fontsize_row=5.5)
print(temp.DM.sel %>% lapply(function(x)x$gene))

#' Fibro
temp.DM.sel <-  total.mk$Fibro$sig %>% group_by(set) %>% top_n(50,power)  %>% split(.,.$set) %>% lapply(function(x) {head(x,25)}) 
temp.DM.sel <- temp.DM.sel[c("Fibro1","Fibro2")]
temp.M <-  data.ob.umap %>% select(cell,sample,EML) %>% filter(EML %in% c("Fibro1","Fibro2","Fibro3")) %>% split(.,.$EML) %>% lapply(function(x){FunMaSF(x,300)}) %>% do.call("bind_rows",.) %>% mutate(od=factor(EML,setdiff(names(temp.DM.sel),c("Fibro","big_monocytes")),ordered = T)) %>% arrange(od,sample) %>% select(-EML) 
temp.exp <- lognormExp.mBN[unlist(lapply(temp.DM.sel,function(x){return(x$gene)})),temp.M$cell] %>% as.data.frame()#lognormExp.mBN
out.limit <- 2.5
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
temp.sel.exp[temp.sel.exp>out.limit] <- out.limit
temp.sel.exp[temp.sel.exp<  (-1*out.limit)] <- -1*out.limit
temp.anno <- temp.M %>% tibble::column_to_rownames("cell") 

Fibro.big.diff.mk <- total.mk$Fibro$detail %>% filter(power > 0.6) %>% filter(pct.1 < 0.25 | pct.2 < 0.25)
data.temp <- subset(data.ob,cells=(data.ob.umap %>% filter(EML %in% c("Fibro1","Fibro2")) %>% pull(cell)))
VlnPlot(data.temp,Fibro.big.diff.mk ,group.by="EML",ncol=6)
cowplot::plot_grid(plotlist=FunFP_plot(data.ob,Fibro.big.diff.mk$gene))
#+ fig.width=15,fig.height=15
pheatmap::pheatmap(temp.sel.exp[,rownames(temp.anno)],cluster_rows=F,cluster_cols=F,scale="none",annotation_col=temp.anno,show_colnames=F,color=(colorRampPalette(c("royalblue3","white","firebrick4"))(50)), show_rownames = T,border_color="NA",main="Top markers",gaps_row = rep(unlist(lapply(temp.DM.sel,function(x){return(nrow(x))})) %>% cumsum(),each=2),fontsize_row=7)
print(temp.DM.sel %>% lapply(function(x)x$gene))


#' big_monocytes
temp.DM.sel <-  total.mk$big_monocytes$sig %>% group_by(set) %>% top_n(15,power)  %>% split(.,.$set) %>% lapply(function(x) {head(x,15)}) 
temp.DM.sel <- temp.DM.sel[c("Dendritic","Mait_Tcell","M1Macro","M2Macro","Monocytes","Spp1_hi_Macro")]
temp.M <-  data.ob.umap %>% select(cell,sample,EML) %>% filter(EML %in% c("Dendritic","Mait_Tcell","M1Macro","M2Macro","Monocytes","Spp1_hi_Macro","M0Macro")) %>% split(.,.$EML) %>% lapply(function(x){FunMaSF(x,150)}) %>% do.call("bind_rows",.) %>% mutate(od=factor(EML,setdiff(c(names(temp.DM.sel),"M0Macro"),c("Fibro","big_monocytes")),ordered = T)) %>% arrange(od,sample) %>% select(-EML) 
temp.exp <- lognormExp.mBN[unlist(lapply(temp.DM.sel,function(x){return(x$gene)})),temp.M$cell] %>% as.data.frame()#lognormExp.mBN
out.limit <- 2.5
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
temp.sel.exp[temp.sel.exp>out.limit] <- out.limit
temp.sel.exp[temp.sel.exp<  (-1*out.limit)] <- -1*out.limit
temp.anno <- temp.M %>% tibble::column_to_rownames("cell") 
#+ fig.width=15,fig.height=15
pheatmap::pheatmap(temp.sel.exp[,rownames(temp.anno)],cluster_rows=F,cluster_cols=F,scale="none",annotation_col=temp.anno,show_colnames=F,color=(colorRampPalette(c("royalblue3","white","firebrick4"))(50)), show_rownames = T,border_color="NA",main="Top markers",gaps_row = rep(unlist(lapply(temp.DM.sel,function(x){return(nrow(x))})) %>% cumsum(),each=2),fontsize_row=5.5)
print(temp.DM.sel %>% lapply(function(x)x$gene))

# Gpr171 Cytip Klrb1b
# VlnPlot(data.ob,c("Ccl9","Cd68","Csf1r","Acp5","Mmp9")) for Dendritic

# for (n in c("M1Macro","M2Macro","Spp1_hi_Macro")) {
#   data.temp <- subset(data.ob,cells=(data.ob.umap %>% filter(EML==n) %>% pull(cell)))
#   DefaultAssay(data.temp)="RNA"
#   VlnPlot(data.temp,c("Dcn"),group.by="sample")
# }


#c("RGD1559482","Psap","G0s2")
total.mk$all$sig %>% bind_rows(total.mk$ma$sig) %>% select(-c(n,NP)) %>% mutate(power=round(power,3))%>% filter(power > 0.4) %>% group_by(set) %>% top_n(100,power) %>% write.table("tmp_data/temp.mk.tsv",quote=F,sep="\t",row.names = F,col.names = T)

# clear
#‘ loading human lineage markers 
#human.mk <- read.delim("doc/PMID35013299.sup2.mk.tsv",stringsAsFactors = F,head=F) %>% tbl_df() %>% inner_join(rats.human %>% mutate(V2=hsapiens_homolog_associated_gene_name) %>% select(external_gene_name,V2),by="V2")%>% split(.,.$V1) %>% lapply(function(x){unique(x$external_gene_name)})
#human.mk[["M2-Macro"]] <- human.mk[["M2-Macro"]] %>% setdiff("Cd74")
#human.mk.sub <- read.delim("doc/PMID35013299.sup3.mk.tsv",stringsAsFactors = F,head=F) %>% tbl_df() %>% inner_join(rats.human %>% mutate(V2=hsapiens_homolog_associated_gene_name) %>% select(external_gene_name,V2),by="V2")%>% split(.,.$V1) %>% lapply(function(x){unique(x$external_gene_name)})
#human.mk.sel <- list(SMCs=c("Tagln", "Acta2"),Fibro=c("Dcn", "Cfd"),T_lympho="Cd3d",CD14_Mono=c("Cd14","S100a9"),DiffKer=c("Krt1","Krt10"),BasalKera=c("Krt5","Krt14"),NK=c("Gzmb"),NKT=c("Cd3d"), CD16_Mono=c("Fcgr3a"), M1_Macro=c("Il1b"),M2_Macro="Cd163",Melano_Schwann=c("MMlana","Cdh19"),LymphEndo="Ccl21" ,DCs=c("Gzmb", "Irf8") ,B_lympho=c("Cd79a","Ms4a1"),Plasma=c("Mzb1"),Mast=c("Tpsab1"),VasEndo=c("Aqp1"))#Sweat_Seba="DCD",Erythro="HBB",
#shi.mk <- list(Myofibroblast=c("Sfrp2","Rgs5"),fibroblast=c("Dpt","Gsn","Cxcl14"),universal_fibroblast=c("Pi16"),Perimysial_fibroblast=c("Col12a1"),HE_Fibro=c("Mmp13"),unk=c("Igfbp5","RGD1564664","Acp5","Mmp9","Mt3"),sMCs=c("Des","Cox4i2","Tnn","Acta2"),Adipocytes_Endothelial=c("Plvap","Rbp7","Fabp4"),Monocytes=c("Fcnb","Lyz2"),M1_Macrophages =c("Il1b"),NK_T=c("Ccl4"),Macrophages=c("C1qa","Spp1","Cd163"),M2a_Macrophages =c("Clec10a"),Erythroid=c("Hbb"),Neutrophils=c("G0s2","S100a9"),Tcell=c("Ccl5"),Keratinocytes_skins=c("Krt17","Ly6d","Dmkn"),T_reg=c("Top2a","Stmn1","Ube2c"),prolifer=c("Mki67"),Bcell=c("Ighm","Cd24"),Dendritic=c("RT1-Da","Cd74","Ccl17"),Endothelial=c("Apln","Krt15"),granulocyte=c("Cpa3","Cma1"),MAST=c("Tpsb2","Tpsab1"),Cardiomyocytes=c("Acta1"))

# B-Lympho/T-Lympho not exists
#' sel.mk <- list()
#' sel.mk$Endo <- c("Plvap","Rbp7","Fabp4","Apln")
#' sel.mk$MyoFibro <- c("Sfrp2","Dpt","Igfbp5")
#' sel.mk$PeriFibro <- c("Col12a1")
#' sel.mk$HeFibro <- c("Mmp13")
#' sel.mk$Kera <- c("Krt17","Ly6d","Dmkn")
#' sel.mk$Mast <- c("Cpa3","Cma1","Tpsb2","Tpsab1")
#' sel.mk$MelScwann <- c("Cyb561a3")
#' sel.mk$Neutrophils=c("G0s2","S100a9")
#' sel.mk$NKT <- c("Nkg7","Gzmk","Ccl5")
#' sel.mk$SMC <- c("RGD1564664","Des","Cox4i2","Acta2")
#' sel.mk$Macro1 <- c("Spp1")
#' sel.mk$Macro2 <- c("Clec10a","Cd163","C1qa","C1qb","C1qc","Selenop")
#' sel.mk$Dendritic <- c("RT1-Da","Cd74")
#' sel.mk$Monocytes <- c("Fcnb","Lyz2")
#' #c("Il1b")
#' #sel.mk$NotSure <- c("Tnn")
#' #sel.mk$Unk2 <- c("Acp5","Mmp9","Mt3")
#' 
#' 
#' cowplot::plot_grid(plotlist=FunFP_plot(data.ob,c("Sparc","Col3a1","Dcn","Epcam","Alb","Pecam1","Msln","Rgs5","Myh11","Top2a","Mki67","Cd24a" )))
#' pdf("temp.mk.pdf",9,9)
#' for (n in names(sel.mk)) {
#'   print(
#'     cowplot::plot_grid(plotlist=FunFP_plot(data.ob,sel.mk[[n]]),nrow=2,ncol=2)
#'   )
#' }
#' dev.off()
#' 
#' pdf("temp.mk2.pdf",9,9)
#' for (n in c("Macro1","Macro2","Dendritic","Monocytes")) {
#'   print(
#'     cowplot::plot_grid(plotlist=FunFP_plot(data.list$ma,sel.mk[[n]]),nrow=2,ncol=2)
#'   )
#' }
#' dev.off()
#' 
#' print(
#   lognormExp.mBN[c("RGD1559482","Dcn"),data.ob.umap %>% filter(EML %in% c("M1Macro","M2Macro","Fibro1","Fibro2")) %>% pull(cell)] %>% tibble::rownames_to_column("Gene") %>% gather(cell,logExp,-Gene) %>% spread(Gene,logExp) %>% left_join(data.ob.umap,by="cell") %>% ggplot+geom_jitter(mapping=aes(x=Dcn,y=RGD1559482,col=EML))+facet_wrap(.~sample)
#' )
# sus.fp.Mf.cells <- lognormExp.mBN[c("RGD1559482","Dcn"),data.ob.umap %>% filter(EML %in% c("M1Macro","M2Macro","Fibro1","Fibro2","M0Macro","Spp1_hi_Macro")) %>% pull(cell)] %>% tibble::rownames_to_column("Gene") %>% gather(cell,logExp,-Gene) %>% spread(Gene,logExp)  %>% left_join(data.ob.umap,by="cell") %>% filter(EML  %in% c("M1Macro","M2Macro"))  %>% tbl_df() %>% filter(Dcn > 4 & RGD1559482 <4) %>% select(cell:RGD1559482) %>% pull(cell)
#' DimPlot(data.ob,cells.highlight=sus.fp.Mf.cells)
#' saveRDS(sus.fp.Mf.cells ,file="tmp_data/sus.fp.Mf.cells.rds")
#' #' show the marker
#' cowplot::plot_grid(plotlist=c(FunFP_plot(data.temp,unlist(gd.sel.mk))))

#temp.DM$sig %>% group_by(set) %>% top_n(5,power) %>% as.data.frame()
#temp.DM= FunRF_FindAllMarkers_para(data.temp)

#Idents(data.temp) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[colnames(data.temp),"EML"])
#Idents(data.ma.ob) <- factor((data.ma.ob.umap %>% tibble::column_to_rownames("cell"))[colnames(data.ma.ob),"subCT"])

#
#data.ob.umap %>% group_by(seurat_clusters, pred) %>% summarise(nCell=n()) %>% filter(!is.na(pred)) %>% group_by(seurat_clusters) %>% mutate(prop=nCell/sum(nCell)) %>% top_n(1,prop)

# for endocrine cell only integration
# if (file.exists(paste0("tmp_data/",TD,"/ma.UMAP.IT.results.rds"))) {
#   data.ma.ob <- readRDS(paste0("tmp_data/",TD,"/ma.UMAP.IT.results.rds"))
#   data.ma.ob.umap <- readRDS(paste0("tmp_data/",TD,"/ma.IT.UMAP.coord.rds"))
# }else{
#   counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
#   lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))
#   rownames(lognormExp.mBN) <- gsub("_","-",rownames(lognormExp.mBN))
#   
#   temp.sel.cells <- data.ob.umap %>% filter(EML %in% c("Monocyte")) %>% pull(cell)
#   temp.M <- meta.filter %>% filter(cell %in% temp.sel.cells)
#   
#   # get the expressed genes ( expressed in at least 1 datatsets)
#   expG.set <- list()
#   for (b in unique(temp.M$sample  %>% unique() %>% as.vector())) { 
#     expG.set[[b]] <- rownames(counts.filter )[rowSums(counts.filter[,(temp.M %>% filter(sample==b) %>% pull(cell))] >=1) >=5]
#   }
#   temp.sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()
#   
#   #
#   npcs=40
#   nGene=3000
#   
#   
#   data.end.temp <-   CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) 
#   data.list <- SplitObject(data.end.temp, split.by = "sample") %>% lapply( function(x) {SCTransform(x,,method = "glmGamPoi", verbose = FALSE,vars.to.regress = c("nGene"))}) #vars.to.regress = c("mt.perc")
#   anchor.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = nGene)
#   
#   #rm(counts.filter)
#   data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = anchor.features ,verbose=F)
#   data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT",anchor.features = anchor.features,verbose=F)
#   
#   data.ma.ob <- IntegrateData(data.anchors, normalization.method = "SCT") %>%  RunPCA( verbose = FALSE,npcs=50)%>% RunUMAP(dims=1:npcs,verbose=F) %>% FindNeighbors( dims = 1:npcs,verbose = FALSE) 
#   
#   data.ma.ob@assays$RNA@data <- as.matrix(lognormExp.mBN[rownames(data.ma.ob@assays$RNA@data),colnames(data.ma.ob@assays$RNA@data)])
#   
#   data.temp <-  data.ma.ob  %>%  FindClusters(resolution = 0.15,verbose = FALSE)
#   data.temp <- RenameIdents(data.temp,"6"="2")
#   #temp.DM= FunRF_FindAllMarkers_para(data.temp)
#   #temp.DM$sig %>% group_by(set) %>% top_n(5,power) %>% as.data.frame()
#   #cowplot::plot_grid(plotlist=FunFP_plot(data.temp,temp.DM$sig %>% group_by(set) %>% top_n(3,power) %>% as.data.frame() %>% pull(gene)))
#   #DimPlot(data.temp,label=T)+theme_void()+NoLegend()+ggtitle("raw_cluster")+FunTitle()
#   
#   data.ma.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,phases:sample,seurat_clusters)) %>% mutate(seurat_clusters=paste0("C",as.vector(Idents(data.temp))))  %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$pca@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:PC_10),by="cell") %>% mutate(EML="Macro")
#   data.ma.ob.umap <- data.ma.ob.umap %>% mutate(subCT=seurat_clusters) %>% mutate(subCT=ifelse(seurat_clusters %in% c("C3"),"M2Macro",subCT)) %>% mutate(subCT=ifelse(seurat_clusters %in% c("C1"),"Dendritic",subCT))%>% mutate(subCT=ifelse(seurat_clusters %in% c("C0"),"M1Macro",subCT))%>% mutate(subCT=ifelse(seurat_clusters %in% c("C2"),"Monocytes",subCT)) %>% mutate(subCT=ifelse(seurat_clusters %in% c("C5"),"Neutro_Macro",subCT)) %>% mutate(subCT=ifelse(seurat_clusters %in% c("C4"),"Unk3",subCT))  
#   
#   saveRDS(data.ma.ob ,file=paste0("tmp_data/",TD,"/ma.UMAP.IT.results.rds"))
#   saveRDS(data.ma.ob.umap,file=paste0("tmp_data/",TD,"/ma.IT.UMAP.coord.rds"))
#   
# }



