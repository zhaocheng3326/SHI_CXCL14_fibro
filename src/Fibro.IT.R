#' ---
#' title: "Defining the UMPA ( fibro cells)"
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

#https://www.nature.com/articles/s41467-018-08247-x fig1 
sel.big.fb.mk <- c("Mmp14", "Spon1", "Twist2", "Ccl2", "Il6", "Notch3", "Pdgfa", "Angptl4", "Col14a1", "Plagl1", "Cdk1",  "Ccl3", "Ccl4", "Trem1")

nGene=2000
npcs=30

data.all.ob.umap <- readRDS(paste0("tmp_data/",TD,"/whole.IT.UMAP.coord.rds")) %>% select(cell:phase_seurat,EML,cluster_EML,big_EML)
#saveRDS(data.ob,paste0("tmp_data/",TD,"/whole.UMAP.IT.results.nonGENE.rds"))
#' CCA for mouse dataset
if (file.exists(paste0("tmp_data/",TD,"/Fibryo.IT.ob.rds"))) {
  data.ob <- readRDS(paste0("tmp_data/",TD,"/Fibryo.IT.ob.rds"))
  data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/Fibryo.IT.UMAP.coord.rds"))
}else{
  counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
  rownames(counts.filter) <- gsub("_","-",rownames(counts.filter))
  lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))
  rownames(lognormExp.mBN) <- gsub("_","-",rownames(lognormExp.mBN))
  #
  
  expG.set <- list()
  for (b in unique(meta.filter$sample  %>% unique() %>% as.vector())) { 
    temp.cell <- data.all.ob.umap %>% filter(EML!="Prolifer_Fibro") %>% filter(big_EML=="Fibro") %>% filter(sample==b) %>% pull(cell)
    expG.set[[b]] <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]
  }
  sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()
  
  temp.sel.expG <-  rownames(lognormExp.mBN) %>% intersect(sel.expG )
  temp.M <-  data.all.ob.umap %>% filter(EML!="Prolifer_Fibro") %>% filter(big_EML=="Fibro")  %>% select(cell:sample,pred,phase_seurat)
  data.temp <-   CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) 
  data.list <- SplitObject(data.temp, split.by = "sample") %>% lapply( function(x) {SCTransform(x,,method = "glmGamPoi",verbose = FALSE)}) 
  #data.list <- data.list[c("Sample_M","Sample_C","Sample_H")]
  #data.list <- data.list[c("Sample_H","Sample_M","Sample_C")]
  #vars.to.regress = c("mt.perc"),#,vars.to.regress = c("nGene")
  #data.list <- SplitObject(data.temp, split.by = "sample") %>% lapply( function(x) {SCTransform(x,,method = "glmGamPoi",verbose = FALSE,vars.to.regress = c("nGene"))}) #vars.to.regress = c("mt.perc"),#,vars.to.regress = c("nGene")
  #data.list <- SplitObject(data.temp, split.by = "sample") %>% lapply( function(x) {SCTransform(x,,method = "glmGamPoi",verbose = FALSE)}) #vars.to.regress = c("mt.perc"),#,vars.to.regress = c("nGene")
  
  anchor.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures =  nGene)
  #rm(counts.filter)
  data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = anchor.features ,verbose=F)
  data.anchors <- FindIntegrationAnchors(object.list = data.list, normalization.method = "SCT",anchor.features = anchor.features,verbose=F)
  data.ob <- IntegrateData(data.anchors, normalization.method = "SCT") %>%  RunPCA( verbose = FALSE,npcs=50)%>% RunUMAP(dims=1:npcs,verbose=F) %>% FindNeighbors( dims = 1:npcs,verbose = FALSE)  #%>% FindNeighbors( dims = 1:npcs,verbose = FALSE,nn.method="annoy",annoy.metric="cosine") #45 not good
  #data.ob <- IntegrateData(data.anchors, normalization.method = "SCT") %>%  RunPCA( verbose = FALSE,npcs=50)%>% RunTSNE(dims=1:npcs,verbose=F)   %>% FindNeighbors( dims = 1:npcs,verbose = FALSE) 
  
  data.ob@assays$RNA@data <- as.matrix(lognormExp.mBN[rownames(data.ob@assays$RNA@data),colnames(data.ob@assays$RNA@data)])
  
  # if (recheck) {
  #   temp.mk.check.list <- list()
  #   for (reso in seq(0.2,0.8,by=0.05)) {
  #     print(reso)
  #     data.temp <- data.ob  %>%  FindClusters(resolution = 0.6,verbose = FALSE)
  #     DefaultAssay(data.temp) <- "RNA"
  #     temp.mk.check.list[[paste0("reso_",reso)]] <- FunRF_FindAllMarkers_para(data.temp)
  #     temp.mk.check.list[[paste0("reso_",reso)]]$n <- length(unique(Idents(data.temp)))
  #   }
  #   saveRDS(temp.mk.check.list,paste0("tmp_data/",TD,"/Fibro.diffReso.mk.rds"))
  # }
  
 
  
  data.temp <- data.ob  %>%  FindClusters(resolution = 0.6,verbose = FALSE)
  
  data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,nGene:sample,pred,phase_seurat))  %>% mutate(EML=paste0("FibroSub_",as.vector(Idents(data.temp)))) %>% mutate(EML=ifelse(EML=="FibroSub_0","FibroSub_3",EML)) %>% mutate(EML=recode(EML,"FibroSub_1"="Fib_SC1"))%>% mutate(EML=recode(EML,"FibroSub_2"="Fib_SC4"))%>% mutate(EML=recode(EML,"FibroSub_3"="Fib_SC3"))%>% mutate(EML=recode(EML,"FibroSub_4"="Fib_SC7"))%>% mutate(EML=recode(EML,"FibroSub_5"="Fib_SC2"))%>% mutate(EML=recode(EML,"FibroSub_6"="Fib_SC6"))%>% mutate(EML=recode(EML,"FibroSub_7"="Fib_SC5"))%>% mutate(EML=recode(EML,"FibroSub_8"="Fib_SC8")) %>% mutate(seurat_clusters=paste0("C",as.vector(Idents(data.temp)))) %>% mutate(big_EML="Fibro")  %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$pca@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:PC_40),by="cell") %>% mutate(cluster_EML=EML,big_EML="Fibro") #%>% left_join(sn.result  %>% lapply( function(x) {return(x$sce.out)}) %>% do.call("bind_rows",.) %>% select(cell:AUC_tSNE2),by="cell")
  # saving object
  saveRDS(data.ob,file=paste0("tmp_data/",TD,"/Fibryo.IT.ob.rds"))
  saveRDS(data.ob.umap,file=paste0("tmp_data/",TD,"/Fibryo.IT.UMAP.coord.rds"))
}

if (file.exists(paste0("tmp_data/",TD,"/marker.rds"))) {
  Fibro.mk <- readRDS(paste0("tmp_data/",TD,"/Fibro.marker.rds"))
}else{
  data.temp <- data.ob
  DefaultAssay(data.temp) <- "RNA"
  Idents(data.temp) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[colnames(data.ob),"EML"])
  data.temp@meta.data$EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[colnames(data.ob),"EML"]
  Fibro.mk <- FunRF_FindAllMarkers_para(data.temp)
  saveRDS(Fibro.mk,file=paste0("tmp_data/",TD,"/Fibro.marker.rds"))
}
lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))
rownames(lognormExp.mBN) <- gsub("_","-",rownames(lognormExp.mBN))

##data.ob  %>%  FindClusters(resolution = seq(0.2,1,by=0.1),verbose = FALSE) %>% clustree( prefix = "integrated_snn_res.")
#data.temp <- data.ob  %>%  FindClusters(resolution = 0.6,verbose = FALSE) # 0.4 , 1 ^ 0.6
#DefaultAssay(data.temp) <- "RNA"
#tp <- DimPlot(data.temp,label=T)+theme_void()+NoLegend()
#cowplot::plot_grid(plotlist=c(FunFP_plot(data.temp,c("Pdgfra","Rgs5","Cyp26a1","Crabp1","Acta2","Mest","Gpx3","Plac8","Creb3","Lyz2","Top2a")),list(tp=tp)))
#cowplot::plot_grid(plotlist=c(FunFP_plot(data.temp,c("Dpp4","Cd24","Prdm1","Acta2","Lrrc15","Pi16","Col15a1","Ccl19","Cxcl12","Comp","Npnt","Hhip","Adamdec1","Cxcl6")),list(tp=tp)))


DefaultAssay(data.ob) <- "RNA"
Idents(data.ob) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[colnames(data.ob),"EML"])
data.ob@meta.data$EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[colnames(data.ob),"EML"]
data.ob@meta.data$big_EML <- (data.all.ob.umap %>% tibble::column_to_rownames("cell"))[colnames(data.ob),"EML"]


cowplot::plot_grid(
  DimPlot(data.ob,label=T,group.by="big_EML")+theme_void()+NoLegend()+ggtitle("pre_EML")+FunTitle(),
  DimPlot(data.ob,label=T)+theme_void()+NoLegend()+ggtitle("raw_cluster")+FunTitle(),
  DimPlot(data.ob,label=T,group.by="sample")+theme_void()+NoLegend()+ggtitle("Sample")+FunTitle(),
  #FeaturePlot(data.ob,"mt.perc"),
  #FeaturePlot(data.ob,"nGene")
  nrow=2,ncol=2
)


#cowplot::plot_grid(plotlist=c(FunFP_plot(data.temp,c("Pdgfra","Rgs5","Cyp26a1","Crabp1","Acta2","Mest","Gpx3","Plac8","Creb3","Lyz2","Top2a")),list(tp=tp)))
#
data.temp=data.ob
tp <- DimPlot(data.temp,label=T)+theme_void()+NoLegend()
cowplot::plot_grid(plotlist=c(FunFP_plot(data.temp,c("Dpp4","Cd24","Prdm1","Acta2","Lrrc15","Pi16","Col15a1","Ccl19","Cxcl12","Comp","Npnt","Hhip","Adamdec1","Cxcl6","Mest")),list(tp=tp)))

# GS
cowplot::plot_grid(plotlist=c(FunFP_plot(data.temp,c("Col12a1","Acta2","Tagln","Cd14","Mgp","Cd24","Acan","Corin","Sox2","Lef1","Crabp1","Pi16","Dpp4","Col15a1","Dpt","Tgfbi","Mest","Plac8","Gpx3","Sparc","Dcn"))))
#"Acan","Corin" # lowly expressed

DimPlot(data.ob,label=T,group.by="EML",split.by="sample")+theme_void()+NoLegend()+ggtitle("Sample")+FunTitle()

print(
  data.ob.umap %>% group_by(sample,cluster_EML) %>% summarise(nCell=n()) %>% group_by(sample) %>% mutate(prop=nCell/sum(nCell)) %>% ggplot+geom_bar(mapping=aes(x=cluster_EML,y=prop,fill=sample),stat="identity",position="dodge")+ theme(axis.text.x=element_text(angle = 90))
)

#' check Fibro markers
temp.DM.sel <-  Fibro.mk$sig %>% filter(n==1) %>% group_by(set) %>% top_n(15,power)  %>% split(.,.$set) %>% lapply(function(x) {head(x,15)}) 
temp.M <-  data.ob.umap %>% select(cell,sample,EML) %>% split(.,.$EML) %>% lapply(function(x){FunMaSF(x,150)}) %>% do.call("bind_rows",.) %>% mutate(od=EML) %>% arrange(od,sample) %>% select(-EML) #%>% mutate(od=factor(EML,setdiff(names(temp.DM.sel),c("Fibro","big_monocytes")),ordered = T))
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


#' loading other marker
#' cowplot::plot_grid(plotlist=FunFP_plot(data.ob,c("Pdgfra","Cd34")))
#' cowplot::plot_grid(plotlist=FunFP_plot(data.ob,sel.big.fb.mk))
#' cowplot::plot_grid(
#'   plotlist=FunFP_plot(data.ob,c("Zo1","Mucin","Cd26","Lrig1","Sca1","Dlk1","Cd24","Igta5","Blimp1"))
#' )
#' cowplot::plot_grid(
#'   plotlist=FunFP_plot(data.ob,c("Tjp1","Muc1","Dpp4","Lig1","Atxn1","Dlk1","Cd24","Igta5","Prdm1"))
#' )
#' #Dpp4 Cd26+,  sub2 (Cd24+)
#' wound.sub12.mk <- read.delim("doc/PMID30737373.sup2.subFibro.mk.tsv",stringsAsFactors = F,head=T) %>% tbl_df() %>% filter(Gene %in% rownames(data.ob@assays$RNA@data))  %>% group_by(Sub.cluster) %>% top_n(3,P_val_adj*-1) %>% split(.,.$Sub.cluster) %>% lapply(function(x) {head(x$Gene,3)})
#' 
#' cowplot::plot_grid(plotlist=FunFP_plot(data.ob,unlist(wound.sub12.mk)))
#' cowplot::plot_grid(plotlist=FunFP_plot(data.ob,c("Lyz2","Pdgfra")))
#' 
#' #' check the co-expression of Lyz2,Ccl6, Pdgfra,  Sma
#' temp.M <- data.all.ob.umap %>% filter(big_EML %in% c("Fibro","Immun")) %>% filter(! EML %in% c("Dendritic","Neutrophil","Mast","NK","NKT")) %>% select(cell,sample,big_EML) %>% mutate(type=ifelse(cell %in% ( data.ob.umap %>% filter(EML=="FibroSub_5") %>% pull(cell)),"Lyz2_pos_Fibro",big_EML))%>% mutate(type=ifelse(type == "Fibro","Lyz2_neg_Fibro",type))
#' temp.out <- lognormExp.mBN[c("Lyz2","Pdgfra"),temp.M$cell] %>% tibble::rownames_to_column("gene")  %>%tbl_df()%>% gather(cell,logExp,-gene) %>% spread(gene,logExp)%>% inner_join(temp.M,by="cell")
#' print(
#'   temp.out %>% ggplot()+geom_point(mapping=aes(x=Lyz2,y=Pdgfra,col=type))+facet_grid(type~sample)
#' )
# cowplot::plot_grid(plotlist=FunFP_plot(data.ob,c("Pdgfra","Rgs5","Cyp26a1","Crabp1","Acta2","Gpx3","Plac8","Creb3","Top2a","Mest","Plac8")))
# sel.mk <- list()
# sel.mk$myofibroblast=c("Runx1","Tcf4","Zeb2")
# cowplot::plot_grid(
#   plotlist=FunFP_plot(data.ob,unlist(sel.mk))
# )
# cowplot::plot_grid(
#   plotlist=FunFP_plot(data.ob,c("Pdgfrb","Crabp1"))
# )


# data.ob  %>%  FindClusters(resolution = 0.35,verbose = FALSE) %>% DimPlot(label=T)+NoAxes()+NoLegend()
# data.temp= data.ob  %>%  FindClusters(resolution = 0.4,verbose = FALSE)
# DefaultAssay(data.temp)="RNA"
# a=FunRF_FindAllMarkers_para(data.temp)
# a$sig %>% filter(power > 0.6) %>% pull(set) %>% table() %>% length()
# unique(Idents(data.temp)) %>% length()
