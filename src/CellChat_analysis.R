#' ---
#' title: "CellChat analysis"
#' output: 
#'  html_document:
#'    code_folding: hide
#' ---



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
#suppressMessages(library(topGO))
#suppressMessages(library(clusterProfiler))
#suppressMessages(library(STRINGdb))
suppressMessages(library(CellChat))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
suppressMessages(library(patchwork))
suppressMessages(library(ComplexHeatmap))
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



#pt.cutoff <- 0.4
#' loading data
meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds"))
# DEG and overlap results
load(paste0("tmp_data/",TD,"/DEG.cluster.logFCtune.Rdata"),verbose=T)
#load(paste0("tmp_data/",TD,"/DEG.cluster.ov.Rdata"),verbose=T)
load("tmp_data/gene.meta.Rdata",verbose=T)
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/whole.IT.UMAP.coord.rds")) %>% rows_update(readRDS(paste0("tmp_data/",TD,"/Fibryo.IT.UMAP.coord.rds")) %>% select(cell,EML,cluster_EML),by="cell") %>% mutate(EML=recode(EML,"Spp1_hi_Macro"="M3Macro_Spp1_hi"))
data.ob.umap <- data.ob.umap %>% mutate(cluster=EML) 
#sel.expG <- rownames(lognormExp.mBN)
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB 
which(CellChatDB.use[["interaction"]]$ligand == "H2-BI") # 1887
CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-1887,]
which(CellChatDB.use[["interaction"]]$ligand == "H2-Ea-ps") #1900
CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-1900,]

temp <- CellChatDB.use$interaction %>% select(interaction_name,ligand,receptor) %>% tbl_df()
temp.complex <-  CellChatDB.use$complex %>% tibble::rownames_to_column("ID") %>% tbl_df() %>% gather(subunit,gene,-ID) %>% filter(gene!="")
temp <- temp %>% left_join(temp.complex %>% rename(ligand=ID),by="ligand") %>% left_join(temp.complex %>% rename(receptor=ID),by="receptor") 

Gene.IT <- temp %>% filter(is.na(gene.x) & is.na(gene.y)) %>% select(interaction_name,ligand,receptor) %>% bind_rows(temp %>% filter(!is.na(gene.x) & is.na(gene.y)) %>% mutate(ligand=gene.x) %>% select(interaction_name,ligand,receptor)) %>% bind_rows(temp %>% filter(is.na(gene.x) & !is.na(gene.y)) %>% mutate(receptor=gene.y) %>% select(interaction_name,ligand,receptor)) %>% bind_rows(temp %>% filter(!is.na(gene.x) & !is.na(gene.y)) %>% mutate(ligand=gene.x,receptor=gene.y) %>% select(interaction_name,ligand,receptor)) %>% unique() %>% inner_join(CellChatDB.use$interaction %>% select(interaction_name,pathway_name) %>% unique(),by="interaction_name")
Gene.IT.short <- Gene.IT %>% select(interaction_name,pathway_name,ligand) %>% rename(gene=ligand) %>% mutate(type="ligand") %>% unique() %>% bind_rows(Gene.IT %>% select(interaction_name,pathway_name,receptor) %>% rename(gene=receptor) %>% mutate(type="receptor") %>% unique()) %>% unique()



#group.cellType <- factor(c("FibroSub_0","FibroSub_1","FibroSub_2","FibroSub_3","FibroSub_4","FibroSub_5","FibroSub_6","FibroSub_7","FibroSub_8","M0Macro","M1Macro","M2Macro","Spp1_hi_Macro"),levels=c(rep("Fibro",9),rep("Macro",4)))



#showDatabaseCategory(CellChatDB)

#data.ob <- readRDS(paste0("tmp_data/",TD,"/whole.UMAP.IT.results.rds"))
#DefaultAssay(data.ob) <- "RNA"
#Idents(data.ob) <- factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[colnames(data.ob),"EML"])
#data.ob@meta.data$EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"EML"]
#data.ob@meta.data$cluster_EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"cluster_EML"]
#data.ob@meta.data$big_EML <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.ob@meta.data),"big_EML"]


if (file.exists(paste0("tmp_data/",TD,"/CellChat.ob.rds"))) {
  CellChat.ob <-  readRDS(paste0("tmp_data/",TD,"/CellChat.ob.rds"))
  #CellChat.merge.ob <- readRDS(paste0("tmp_data/",TD,"/CellChat.merge.ob.rds"))
  print("done with preprocessing")
}else{
  counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
  lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))
  #rownames(lognormExp.mBN) <- gsub("_","-",rownames(lognormExp.mBN))
  #rownames(counts.filter) <- gsub("_","-",rownames(counts.filter))
  
  CellChat.ob <- list()
  for ( sa in unique(data.ob.umap$sample)) {
    print(sa)
    temp.M <- data.ob.umap %>% select(cell:UMAP_2,EML:cluster)%>% filter(sample==sa) %>% filter(!EML %in% c("Prolifer_SMC","Prolifer_Fibro"))
    temp_log_counts <- lognormExp.mBN[,temp.M$cell]
    cellchat <- createCellChat(object = as.matrix(temp_log_counts), meta = (temp.M %>% tibble::column_to_rownames("cell")), group.by = "EML")
    
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    
    
    cellchat <- identifyOverExpressedGenes( cellchat)
    cellchat <- identifyOverExpressedInteractions( cellchat)
    
    #'  the projection reduces the dropout effects of signaling genes, 
    cellchat <- projectData(cellchat, PPI.mouse)
    
    cellchat <- computeCommunProb(cellchat)
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    
    #Infer the cell-cell communication at a signaling pathway level
    cellchat <- computeCommunProbPathway(cellchat)
    
    #Calculate the aggregated cell-cell communication network
    cellchat <- aggregateNet(cellchat)
    # Compute the network centrality scores
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
    CellChat.ob[[sa]] <-  cellchat
  }
  saveRDS(CellChat.ob,file=paste0("tmp_data/",TD,"/CellChat.ob.rds"))
  # CellChat.merge.ob <- list()
  # for ( sa in unique(data.ob.umap$sample)) {
  #   print(sa)
  #   temp.M <- data.ob.umap %>% select(cell:UMAP_2,EML:cluster)%>% filter(sample==sa) %>% filter(!EML %in% c("Prolifer_SMC","Prolifer_Fibro")) %>% mutate(EML=ifelse(EML %in% c("FibroSub_0","FibroSub_1","FibroSub_2","FibroSub_3","FibroSub_4","FibroSub_5","FibroSub_6","FibroSub_7","FibroSub_8"),"Fibro",EML)) %>% mutate(EML=ifelse(EML %in% c("M0Macro","M1Macro","M2Macro","Spp1_hi_Macro"),"Macro",EML))
  #   temp_log_counts <- lognormExp.mBN[,temp.M$cell]
  #   cellchat <- createCellChat(object = as.matrix(temp_log_counts), meta = (temp.M %>% tibble::column_to_rownames("cell")), group.by = "EML")
  # 
  #   cellchat@DB <- CellChatDB.use
  #   cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  # 
  # 
  #   cellchat <- identifyOverExpressedGenes( cellchat)
  #   cellchat <- identifyOverExpressedInteractions( cellchat)
  # 
  #   cellchat <- computeCommunProb(cellchat)
  #   # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  #   cellchat <- filterCommunication(cellchat, min.cells = 10)
  # 
  #   #Infer the cell-cell communication at a signaling pathway level
  #   cellchat <- computeCommunProbPathway(cellchat)
  # 
  #   #Calculate the aggregated cell-cell communication network
  #   cellchat <- aggregateNet(cellchat)
  #   # Compute the network centrality scores
  #   cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  #   CellChat.merge.ob[[sa]] <-  cellchat
  # 
  # }
    
   
  #cellchat <- mergeCellChat(CellChat.ob, add.names = names(CellChat.ob))
  #cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
  #cellchat <- netEmbedding(cellchat, type = "functional")
  #cellchat <- netClustering(cellchat, type = "functional")
  #cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
  #cellchat <- netEmbedding(cellchat, type = "structural")
  #cellchat <- netClustering(cellchat, type = "structural")
  #> Classification learning of the signaling networks for datasets 1 2
  
  #saveRDS(CellChat.merge.ob,file=paste0("tmp_data/",TD,"/CellChat.merge.ob.rds"))
}

#CellChat.ob <-  readRDS(paste0("tmp_data/",TD,"/CellChat.ob.rds"))
#CellChat.merge.ob <- readRDS(paste0("tmp_data/",TD,"/CellChat.merge.ob.rds"))

# rownames is source, colnames is target
#' CellChat.ob <-  readRDS(paste0("tmp_data/",TD,"/CellChat.merge.ob.rds"))
#' #' check total interactions counts which passed threshhold (merge the fibro and macro)
#' cc.out <- list()
#' for (sa in c("Sample_C","Sample_H","Sample_M")) {
#'   temp.list <- list()
#   temp.list$cluster <- CellChat.ob[[sa]]@net$count %>% as.data.frame() %>% tibble::rownames_to_column("pair1") %>% tbl_df() %>% gather(pair2,Valid_IT_count,-pair1) %>% full_join( CellChat.ob[[sa]]@net$weight %>% as.data.frame() %>% tibble::rownames_to_column("pair1") %>% tbl_df() %>% gather(pair2,weight,-pair1),by=c("pair1","pair2"))%>% mutate(Sample=sa) %>% rename(source=pair1,target=pair2)
#'   temp.list$pathway <- subsetCommunication(CellChat.ob[[sa]],slot.name="netP") %>% tbl_df()%>% mutate(Sample=sa)
#'   temp.list$it <-  subsetCommunication(CellChat.ob[[sa]]) %>% tbl_df() %>% select(-c(interaction_name_2,annotation,ligand,receptor,evidence))%>% mutate(Sample=sa)
#'   cc.out[[sa]] <- temp.list
#' }
#' 
#' #' check the total interaction weight 
#' temp.list <- lapply(cc.out,function(x){x$cluster})
#' #temp.out <- temp.list %>% do.call("bind_rows",.) %>% mutate(source=ifelse(source %in% c("FibroSub_0","FibroSub_1","FibroSub_2","FibroSub_3","FibroSub_4","FibroSub_5","FibroSub_6","FibroSub_7","FibroSub_8"),"Fibro",source)) %>% mutate(source=ifelse(source %in% c("M0Macro","M1Macro","M2Macro","Spp1_hi_Macro"),"Macro",source)) %>% mutate(target=ifelse(target %in% c("FibroSub_0","FibroSub_1","FibroSub_2","FibroSub_3","FibroSub_4","FibroSub_5","FibroSub_6","FibroSub_7","FibroSub_8"),"Fibro",target)) %>% mutate(target=ifelse(target %in% c("M0Macro","M1Macro","M2Macro","Spp1_hi_Macro"),"Macro",target)) %>% group_by(source,target,Sample) %>% summarise(Valid_IT_count=sum(Valid_IT_count),weight=sum(weight)) %>% mutate(SID=paste(source,target,sep=":"))%>% ungroup()#%>% spread(Sample,weight) %>% filter((Sample_C+Sample_H+Sample_M) >0 ) 
#' temp.out <- temp.list %>% do.call("bind_rows",.) %>% group_by(source,target,Sample) %>% summarise(Valid_IT_count=sum(Valid_IT_count),weight=sum(weight)) %>% mutate(SID=paste(source,target,sep=":"))%>% ungroup()#%>% spread(Sample,weight) %>% filter((Sample_C+Sample_H+Sample_M) >0 ) 
#' 
#' cowplot::plot_grid(
#'   temp.out %>% group_by(source,Sample)%>% summarise(Valid_IT_count=sum(Valid_IT_count),weight=sum(weight)) %>% ggplot()+geom_bar(mapping=aes(x=source,fill=Sample,y=weight),stat="identity",position = "dodge") + theme_classic() + theme(axis.text.x=element_text(angle = 90))+ggtitle("Source,(inclu same group)")+FunTitle(),
#'   temp.out %>% group_by(target,Sample)%>% summarise(Valid_IT_count=sum(Valid_IT_count),weight=sum(weight)) %>% ggplot()+geom_bar(mapping=aes(x=target,fill=Sample,y=weight),stat="identity",position = "dodge") + theme_classic() + theme(axis.text.x=element_text(angle = 90))+ggtitle("Target,(inclu same group)")+FunTitle(),
#'   temp.out %>% filter(source !=target) %>% group_by(source,Sample)%>% summarise(Valid_IT_count=sum(Valid_IT_count),weight=sum(weight)) %>% ggplot()+geom_bar(mapping=aes(x=source,fill=Sample,y=weight),stat="identity",position = "dodge") + theme_classic() + theme(axis.text.x=element_text(angle = 90))+ggtitle("Source,(exclu same group)")+FunTitle(),
#'   temp.out  %>% filter(source !=target)  %>% group_by(target,Sample)%>% summarise(Valid_IT_count=sum(Valid_IT_count),weight=sum(weight)) %>% ggplot()+geom_bar(mapping=aes(x=target,fill=Sample,y=weight),stat="identity",position = "dodge") + theme_classic() + theme(axis.text.x=element_text(angle = 90))+ggtitle("Target,(exclu same group)")+FunTitle(),
#'   nrow=2,ncol=2
#' )
#' 


netVisual_heatmap(cellchat, measure = "weight")
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
rankSimilarity(cellchat, type = "functional")
#' VISFATIN

#' check total interactions counts which passed threshhold (merge the fibro and macro)
cc.out <- list()
for (sa in c("Sample_C","Sample_H","Sample_M")) {
  temp.list <- list()
  temp.list$cluster <- CellChat.ob[[sa]]@net$count %>% as.data.frame() %>% tibble::rownames_to_column("pair1") %>% tbl_df() %>% gather(pair2,Valid_IT_count,-pair1) %>% full_join( CellChat.ob[[sa]]@net$weight %>% as.data.frame() %>% tibble::rownames_to_column("pair1") %>% tbl_df() %>% gather(pair2,weight,-pair1),by=c("pair1","pair2"))%>% mutate(Sample=sa) %>% rename(source=pair1,target=pair2)
  temp.list$pathway <- subsetCommunication(CellChat.ob[[sa]],slot.name="netP",thresh=0.05) %>% tbl_df()%>% mutate(Sample=sa)
  temp.list$it <-  subsetCommunication(CellChat.ob[[sa]],thresh=0.05) %>% tbl_df() %>% select(-c(interaction_name_2,annotation,ligand,receptor,evidence))%>% mutate(Sample=sa)
  cc.out[[sa]] <- temp.list
}

#' check the total interaction weight  (only subgroups)
temp.list <- lapply(cc.out,function(x){x$pathway})
temp.out <- temp.list %>% do.call("bind_rows",.) 

temp.sel.path <- "TGFb"
temp <- temp.out %>% filter(pathway_name==temp.sel.path)%>% group_by(source,Sample)%>% summarise(Outgoing=sum(prob)) %>% ungroup() %>% rename(EML=source) %>% full_join(temp.out %>% filter(pathway_name==temp.sel.path)%>% group_by(target,Sample)%>% summarise(Incoming=sum(prob)) %>% ungroup() %>% rename(EML=target),by=c("Sample","EML")) %>% replace(is.na(.), 0)

temp.plot <- list()
for (sa in c("Sample_C","Sample_H","Sample_M")) {
  temp.plot[[sa]] <- temp %>% filter(Sample==sa)%>% ggplot()+geom_point(mapping=aes(x=Outgoing,y=Incoming,col=EML,shape=Sample),size=3) +ggrepel::geom_text_repel(temp %>% filter(Sample==sa) %>% filter(EML %in% c("FibroSub_6","FibroSub_1","M1Macro")) %>% mutate(SID=EML),mapping=aes(x=Outgoing,y=Incoming,label=SID),box.padding = 0.5)+ theme_classic() +ggtitle(paste(temp.sel.path,sa ))+FunTitle()+xlab("Outgoing interaction strength")+ylab("Incoming interaction strength")+NoLegend()
}
cowplot::plot_grid(plotlist=temp.plot,ncol=3)



# num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
# weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
# gg <- list()
# for (i in 1:length(object.list)) {
#   gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
# }

levels(CellChat.ob$Sample_C@idents)
pdf("tmp_data/fig_pdf/temp.Interaction.Hier.pdf",9,9)
temp.sel.path <- "TGFb"
# netVisual_aggregate(subsetCellChat(CellChat.ob$Sample_C,idents.use= c("M0Macro","M1Macro","M2Macro","Spp1_hi_Macro",paste0("FibroSub_",1:8))), signaling = temp.sel.path, layout = "hierarchy",vertex.receiver=c(1:8))
# netVisual_aggregate(subsetCellChat(CellChat.ob$Sample_H,idents.use= c("M0Macro","M1Macro","M2Macro","Spp1_hi_Macro",paste0("FibroSub_",1:8))), signaling = temp.sel.path, layout = "hierarchy",vertex.receiver=c(1:8))
# netVisual_aggregate(subsetCellChat(CellChat.ob$Sample_M,idents.use= c("M0Macro","M1Macro","M2Macro","Spp1_hi_Macro",paste0("FibroSub_",1:8))), signaling = temp.sel.path, layout = "hierarchy",vertex.receiver=c(1:8))
netVisual_aggregate(CellChat.ob$Sample_C, signaling = temp.sel.path, layout = "hierarchy",vertex.receiver=c(3:10))
netVisual_aggregate(CellChat.ob$Sample_H, signaling = temp.sel.path, layout = "hierarchy",vertex.receiver=c(3:10))
netVisual_aggregate(CellChat.ob$Sample_M, signaling = temp.sel.path, layout = "hierarchy",vertex.receiver=c(3:10))
dev.off()
#netVisual_aggregate(CellChat.ob$Sample_H, sources.use=c(13,14,15,23),targets.use=c(3:11),idents.use=c(13,14,15,23,3:11), signaling = "TGFb", layout = "circle")
#netVisual_aggregate(CellChat.ob$Sample_M,idents.use=c(13,14,15,23,3:11), signaling = "TGFb", layout = "circle")

levels(CellChat.ob$Sample_C@idents)
pdf("temp.pdf",6,12)
temp.sel.path <- "TGFb"
netAnalysis_signalingRole_network(CellChat.ob$Sample_C, signaling = temp.sel.path)
netAnalysis_signalingRole_network(CellChat.ob$Sample_H, signaling = temp.sel.path)
netAnalysis_signalingRole_network(CellChat.ob$Sample_M, signaling = temp.sel.path)
dev.off()


object.list <- CellChat.ob[c("Sample_C","Sample_M")]
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for datasets 1
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for datasets 1 2
# Visualization in 2D-space
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")

pdf("temp.pdf",12,6)
gg1+gg2
dev.off()


temp.list <- lapply(cc.out,function(x){x$cluster})
temp.out <- temp.list %>% do.call("bind_rows",.)

p1=temp.out %>% group_by(Sample)%>% summarise(Valid_IT_count=sum(Valid_IT_count),weight=sum(weight)) %>% ggplot()+geom_bar(mapping=aes(x=Sample,y=Valid_IT_count,fill=Sample),stat="identity") + theme_classic() + theme(axis.text.x=element_text(angle = 90))+ggtitle("Number of Interaction")+FunTitle()

p2=temp.out %>% group_by(Sample)%>% summarise(Valid_IT_count=sum(Valid_IT_count),weight=sum(weight)) %>% ggplot()+geom_bar(mapping=aes(x=Sample,y=weight,fill=Sample),stat="identity") + theme_classic() + theme(axis.text.x=element_text(angle = 90))+ggtitle("Interaction Strength")+FunTitle()
p1+p2

cowplot::plot_grid(
  temp.out %>% group_by(source,Sample)%>% summarise(Valid_IT_count=sum(Valid_IT_count),weight=sum(weight)) %>% ggplot()+geom_bar(mapping=aes(x=source,fill=Sample,y=weight),stat="identity",position = "dodge") + theme_classic() + theme(axis.text.x=element_text(angle = 90))+ggtitle("Source,(inclu same group)")+FunTitle(),
  temp.out %>% group_by(target,Sample)%>% summarise(Valid_IT_count=sum(Valid_IT_count),weight=sum(weight)) %>% ggplot()+geom_bar(mapping=aes(x=target,fill=Sample,y=weight),stat="identity",position = "dodge") + theme_classic() + theme(axis.text.x=element_text(angle = 90))+ggtitle("Target,(inclu same group)")+FunTitle(),
  temp.out   %>% filter(source !=target) %>% group_by(source,Sample)%>% summarise(Valid_IT_count=sum(Valid_IT_count),weight=sum(weight)) %>% ggplot()+geom_bar(mapping=aes(x=source,fill=Sample,y=weight),stat="identity",position = "dodge") + theme_classic() + theme(axis.text.x=element_text(angle = 90))+ggtitle("Source,(exclu same group)")+FunTitle(),
  temp.out  %>% filter(source !=target)  %>% group_by(target,Sample)%>% summarise(Valid_IT_count=sum(Valid_IT_count),weight=sum(weight)) %>% ggplot()+geom_bar(mapping=aes(x=target,fill=Sample,y=weight),stat="identity",position = "dodge") + theme_classic() + theme(axis.text.x=element_text(angle = 90))+ggtitle("Target,(exclu same group)")+FunTitle(),
  nrow=2,ncol=2
)

#' interaction tables
temp.list <- lapply(cc.out,function(x){x$cluster})
temp.out <- temp.list %>% do.call("bind_rows",.)%>% mutate(SID=paste(source,target,sep=":"))%>% select(-Valid_IT_count)%>% spread(Sample,weight) %>% filter((Sample_C+Sample_H+Sample_M) >0 ) 
temp.out.sel <- temp.out %>% filter((Sample_C < Sample_H & Sample_H < Sample_M) | (Sample_C > Sample_H & Sample_H > Sample_M))%>% mutate(sum_value=Sample_C+Sample_H+Sample_M) %>% mutate(UpDown=ifelse(Sample_C > Sample_H,"DownRe","UpRe"),logFC_H_vs_C=log2((Sample_H+0.1)/(Sample_C+0.1)),logFC_M_vs_H=log2((Sample_M+0.1)/(Sample_H+0.1))) %>% arrange(desc(sum_value)) %>% arrange(desc(logFC_H_vs_C))

temp.list <- lapply(cc.out,function(x){x$pathway})
temp.out <- temp.list %>% do.call("bind_rows",.) %>% filter(pval < 0.05) %>% mutate(SID=paste(source,target,pathway_name,sep=":"))%>% select(-pval)%>% spread(Sample,prob) %>% filter((Sample_C+Sample_H+Sample_M) >0 ) 
temp.out.sel <- temp.out %>% filter((Sample_C < Sample_H & Sample_C < Sample_M) | (Sample_C > Sample_H & Sample_C > Sample_M))%>% mutate(sum_value=Sample_C+Sample_H+Sample_M) %>% mutate(UpDown=ifelse(Sample_C > Sample_H,"DownRe","UpRe"),logFC_H_vs_C=log2((Sample_H+0.1)/(Sample_C+0.1)),logFC_M_vs_C=log2((Sample_M+0.1)/(Sample_C+0.1))) %>% arrange(desc(sum_value)) %>% arrange(desc(logFC_H_vs_C))



temp.sel.path.set <- temp.list %>% do.call("bind_rows",.)  %>% filter(pval < 0.05 & prob > 0.025) %>% pull(pathway_name) %>% unique()

levels(CellChat.ob$Sample_C@idents)

#vertex.receiver = seq(3,11) # a numeric vector. 
vertex.receiver = c(12,13,14,22) # a numeric vector. 
pdf("temp.pdf")
for (temp.sel.path in temp.sel.path.set ) {
  netVisual_aggregate(CellChat.ob$Sample_C, signaling = temp.sel.path,  vertex.receiver = vertex.receiver,layout="hierarchy")
  netVisual_aggregate(CellChat.ob$Sample_H, signaling = temp.sel.path,  vertex.receiver = vertex.receiver,layout="hierarchy")
  netVisual_aggregate(CellChat.ob$Sample_M, signaling = temp.sel.path,  vertex.receiver = vertex.receiver,layout="hierarchy")
}


levels(CellChat.ob$Sample_C@idents)

netVisual_aggregate(subsetCellChat(CellChat.ob$Sample_C,idents.use= c("M0Macro","M1Macro","M2Macro","Spp1_hi_Macro",paste0("FibroSub_",0:8))), signaling = "TGFb", layout = "circle")
netVisual_aggregate(subsetCellChat(CellChat.ob$Sample_H,idents.use= c("M0Macro","M1Macro","M2Macro","Spp1_hi_Macro",paste0("FibroSub_",0:8))), signaling = "TGFb", layout = "circle")
netVisual_aggregate(subsetCellChat(CellChat.ob$Sample_M,idents.use= c("M0Macro","M1Macro","M2Macro","Spp1_hi_Macro",paste0("FibroSub_",0:8))), signaling = "TGFb", layout = "circle")
#netVisual_aggregate(CellChat.ob$Sample_H, sources.use=c(13,14,15,23),targets.use=c(3:11),idents.use=c(13,14,15,23,3:11), signaling = "TGFb", layout = "circle")
#netVisual_aggregate(CellChat.ob$Sample_M,idents.use=c(13,14,15,23,3:11), signaling = "TGFb", layout = "circle")
dev.off()

gd.sel.path <- c('TWEAK',"GALECTIN","GDF","SEMA4","TGFb","PTN","FN1","ncWNT","GRN","NCAM","VCAM","CD99","FN1","GAS","SEMA3")
pdf("temp.pdf")
temp.plot <- list()
for (temp.sel.path in gd.sel.path ) {
  temp.plot[[temp.sel.path]] <- list()
  netAnalysis_signalingRole_network(CellChat.ob$Sample_C, signaling = temp.sel.path)
  #temp.plot[[temp.sel.path]]$p1 <-  recordPlot()
  netAnalysis_signalingRole_network(CellChat.ob$Sample_H, signaling = temp.sel.path)
  #temp.plot[[temp.sel.path]]$p2 <-  recordPlot()
  netAnalysis_signalingRole_network(CellChat.ob$Sample_M, signaling = temp.sel.path)
  #temp.plot[[temp.sel.path]]$p3 <-  recordPlot()
  
}
dev.off()


temp.list <- lapply(cc.out,function(x){x$pathway})
temp.out <- temp.list %>% do.call("bind_rows",.) %>% filter(pval < 0.05) %>% mutate(SID=paste(source,target,pathway_name,sep=":"))
further.sel.path <- c("TGFb","ncWNT","NCAM","VCAM")
temp.plot <- list()
for (temp.sel.path in further.sel.path ) {
  temp <- temp.out %>% filter(pathway_name==temp.sel.path)%>% group_by(source,Sample)%>% summarise(Outgoing=sum(prob)) %>% ungroup() %>% rename(EML=source) %>% full_join(temp.out %>% filter(pathway_name==temp.sel.path)%>% group_by(target,Sample)%>% summarise(Incoming=sum(prob)) %>% ungroup() %>% rename(EML=target),by=c("Sample","EML")) %>% replace(is.na(.), 0)
  
  temp.plot[[temp.sel.path]] <- temp %>% ggplot()+geom_point(mapping=aes(x=Outgoing,y=Incoming,col=EML,shape=Sample),size=3) +ggrepel::geom_text_repel(temp %>% filter(EML %in% c("FibroSub_6","FibroSub_1","M1Macro")) %>% mutate(Sample=gsub("Sample_","",Sample)) %>% mutate(SID=paste(Sample,EML,sep=":")),mapping=aes(x=Outgoing,y=Incoming,label=SID),box.padding = 0.5)+ theme_classic() +ggtitle(temp.sel.path )+FunTitle()+xlab("Outgoing interaction strength")+ylab("Incoming interaction strength")+NoLegend()
}
cowplot::plot_grid(plotlist=temp.plot)
#' check the select pathway interaction strenth
#%>% spread(Sample,weight) %>% filter((Sample_C+Sample_H+Sample_M) >0 ) 

cowplot::plot_grid(
  temp.out %>% group_by(source,Sample)%>% summarise(Valid_IT_count=sum(Valid_IT_count),weight=sum(weight)) %>% ggplot()+geom_bar(mapping=aes(x=source,fill=Sample,y=weight),stat="identity",position = "dodge") + theme_classic() + theme(axis.text.x=element_text(angle = 90))+ggtitle("Source,(inclu same group)")+FunTitle(),
  temp.out %>% group_by(target,Sample)%>% summarise(Valid_IT_count=sum(Valid_IT_count),weight=sum(weight)) %>% ggplot()+geom_bar(mapping=aes(x=target,fill=Sample,y=weight),stat="identity",position = "dodge") + theme_classic() + theme(axis.text.x=element_text(angle = 90))+ggtitle("Target,(inclu same group)")+FunTitle(),
  temp.out   %>% filter(source !=target) %>% group_by(source,Sample)%>% summarise(Valid_IT_count=sum(Valid_IT_count),weight=sum(weight)) %>% ggplot()+geom_bar(mapping=aes(x=source,fill=Sample,y=weight),stat="identity",position = "dodge") + theme_classic() + theme(axis.text.x=element_text(angle = 90))+ggtitle("Source,(exclu same group)")+FunTitle(),
  temp.out  %>% filter(source !=target)  %>% group_by(target,Sample)%>% summarise(Valid_IT_count=sum(Valid_IT_count),weight=sum(weight)) %>% ggplot()+geom_bar(mapping=aes(x=target,fill=Sample,y=weight),stat="identity",position = "dodge") + theme_classic() + theme(axis.text.x=element_text(angle = 90))+ggtitle("Target,(exclu same group)")+FunTitle(),
  nrow=2,ncol=2
)


#temp.out.sel  %>% filter(grepl("Fibro",SID) & grepl("Macro",SID))

# select one cluster as example ( Spp1_hi_Macro)
#sel.cluster="Spp1_hi_Macro"

#' output the interaction pathway which has the lots of overlpapped DEGs
temp.list <- lapply(cc.out,function(x){x$pathway})
temp.out <- temp.list %>% do.call("bind_rows",.)%>% mutate(SID=paste(source,target,pathway_name,sep=":")) %>% rename(weight=prob)

temp.ov.genes <- data.frame(gene=DEG.ov[[sel.cluster]]$ov.upDEG,UpDown="up-re") %>% tbl_df() %>% bind_rows(data.frame(gene=DEG.ov[[sel.cluster]]$ov.downDEG,UpDown="down-re") %>% tbl_df())
temp.ov.genes.asso.pathway <- Gene.IT.short %>% inner_join(temp.ov.genes,by="gene")  %>% group_by(pathway_name,UpDown) %>% summarise(nGene=n_distinct(gene)) %>% arrange(desc(nGene)) %>% left_join(Gene.IT.short %>% group_by(pathway_name)  %>% summarise(nBGGene=n_distinct(gene)) ,by="pathway_name") %>% filter(nGene>1)

temp.sel.path <-  temp.out %>% filter((source==sel.cluster| target==sel.cluster ) & pathway_name%in% temp.ov.genes.asso.pathway$pathway_name) %>% filter(pval < 0.05) %>% pull(pathway_name) %>% unique()
temp.plot <- list()
for (sp in temp.sel.path ) {
  temp.sel.SID <- temp.out %>% filter((source==sel.cluster| target==sel.cluster ) & pathway_name==sp )%>% filter(pval < 0.05) %>% pull(SID) %>% unique()
  temp.plot[[paste(sel.cluster,sp,sep=":")]] <- temp.out %>% filter( SID %in% temp.sel.SID) %>% ggplot+geom_point(mapping=aes(x=Sample,y=SID,size=weight,col=-1*log10(pval+10^-50))) +ggtitle(paste0("sel.cluster:",sel.cluster,",sel.path:",sp))
}

for (sel.target.cluster in c("FibroSub_0","FibroSub_1","FibroSub_2","FibroSub_4")) {
  temp.sel.SID <- temp.out %>% filter(source==sel.cluster & target==sel.target.cluster)  %>% filter(pval < 0.05) %>% pull(SID) %>% unique()
  temp.plot[[paste(sel.cluster,sel.target.cluster)]] <- temp.out %>% filter(source==sel.cluster & target==sel.target.cluster & SID %in% temp.sel.SID) %>% ggplot+geom_point(mapping=aes(x=Sample,y=SID,size=weight,col=-1*log10(pval+10^-50))) +ggtitle(paste0("source:",sel.cluster,",target:",sel.target.cluster))
}



#' select one target cluster as example ( Spp1_hi_Macro)
sel.target.cluster="Dendritic"
temp.list <- lapply(cc.out,function(x){x$pathway})
temp.out <- temp.list %>% do.call("bind_rows",.)%>% mutate(SID=paste(source,target,pathway_name,sep=":")) %>% rename(weight=prob)
temp.plot <- list()
for (sel.path in c("TENASCIN","TGFb")) {
  temp.sel.SID <- temp.out %>% filter(target==sel.target.cluster & pathway_name==sel.path) %>% filter(pval < 0.05) %>% pull(SID) %>% unique()
  temp.plot[[paste(sel.target.cluster,sel.path)]] <- temp.out %>% filter(target==sel.target.cluster & pathway_name==sel.path & SID %in% temp.sel.SID) %>% ggplot+geom_point(mapping=aes(x=Sample,y=SID,size=weight,col=-1*log10(pval+10^-50))) +ggtitle(paste0("target:",sel.target.cluster,",sel.path:",sel.path))
}

for (sel.target.cluster in c("FibroSub_0","FibroSub_1","FibroSub_2","FibroSub_4")) {
  temp.sel.SID <- temp.out %>% filter(source==sel.cluster & target==sel.target.cluster)  %>% filter(pval < 0.05) %>% pull(SID) %>% unique()
  temp.plot[[paste(sel.cluster,sel.target.cluster)]] <- temp.out %>% filter(source==sel.cluster & target==sel.target.cluster & SID %in% temp.sel.SID) %>% ggplot+geom_point(mapping=aes(x=Sample,y=SID,size=weight,col=-1*log10(pval+10^-50))) +ggtitle(paste0("source:",sel.cluster,",target:",sel.target.cluster))
}





#weight.max <- getMaxWeight(CellChat.ob,slot.name = c("idents", "netP"), attribute = c("idents","weight"))
weight.max <-1
temp.plot <- list()
for (sa in c("Sample_C","Sample_H","Sample_M")) {
  mat <- CellChat.ob[[sa]]@net$weight
  i=which(rownames(mat)==sel.cluster)
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  temp.plot[[sa]] <- netVisual_circle(mat2, vertex.weight = length(unique(data.ob.umap$EML)), weight.scale = T, edge.weight.max = weight.max, title.name = paste(sa,sel.cluster))
}
#temp.list[[sa]] <-  netVisual_aggregate(CellChat.ob[[sa]], signaling = sel.path, layout = "circle")
#temp.list[[sa]] <-  netVisual_aggregate(CellChat.ob[[sa]], signaling = sel.path, layout = "chord")
#temp.list[[sa]] <- netVisual_heatmap(CellChat.ob[[sa]], signaling = sel.path, color.heatmap = "Reds")




#' select one pathway as example ( Spp1_hi_Macro)
#' check the COLLAGEN pathway
weight.max <- getMaxWeight(CellChat.ob, attribute = c("idents","count"))
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(CellChat.ob)) {
  netVisual_aggregate(CellChat.ob[[i]], signaling = sel.path, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, )#title.name = paste0("Number of interactions - ", names(CellChat.ob)[i])
}


sel.path="TGFb"
temp.list <- list()
for (sa in c("Sample_C","Sample_H","Sample_M")) {
  #temp.list[[sa]] <-  netVisual_aggregate(CellChat.ob[[sa]], signaling = sel.path, layout = "circle")
  #temp.list[[sa]] <-  netVisual_aggregate(CellChat.ob[[sa]], signaling = sel.path, layout = "chord")
  temp.list[[sa]] <- netVisual_heatmap(CellChat.ob[[sa]], signaling = sel.path, color.heatmap = "Reds")
}
temp.list$Sample_C+temp.list$Sample_H+temp.list$Sample_M

netAnalysis_contribution(CellChat.ob$Sample_C, signaling = sel.path)



num.link <- sapply(CellChat.ob, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
temp.list <- list()
for (sa in c("Sample_C","Sample_H","Sample_M")) {
  temp.list[[sa]] <-netAnalysis_signalingRole_scatter(CellChat.ob[[sa]],weight.MinMax = weight.MinMax)
}
temp.list$Sample_C+temp.list$Sample_H+temp.list$Sample_M

temp.plot <- list()
temp.plot$p1 <- rankSimilarity(cellchat, type = "functional",comparison2 = c(1,2))
temp.plot$p2 <- rankSimilarity(cellchat, type = "functional",comparison2 = c(1,3))
temp.plot$p3 <- rankSimilarity(cellchat, type = "functional",comparison2 = c(2,3))
temp.plot$p1+temp.plot$p2+temp.plot$p3



#' pathway reads (split by income and outcome for each cluster and pathway)
temp.list <- lapply(cc.out,function(x){x$pathway})
  temp.out <- temp.list %>% do.call("bind_rows",.) %>% filter(pval < 0.05) %>% mutate(SID=paste(source,target,pathway_name,sep=":"))%>% select(-c(pval)) %>% rename(weight=prob)%>% spread(Sample,weight) %>% filter((Sample_C+Sample_H+Sample_M) >0 ) 
temp.out.sel <- temp.out %>% filter((Sample_C < Sample_H & Sample_H < Sample_M) | (Sample_C > Sample_H & Sample_H > Sample_M))%>% mutate(sum_value=Sample_C+Sample_H+Sample_M) %>% mutate(UpDown=ifelse(Sample_C > Sample_H,"DownRe","UpRe"),logFC_H_vs_C=log2((Sample_H+0.1)/(Sample_C+0.1)),logFC_M_vs_H=log2((Sample_M+0.1)/(Sample_H+0.1))) %>% arrange(desc(sum_value)) %>% arrange(desc(abs(logFC_H_vs_C)))

temp.out %>% select(SID,Sample_C,Sample_H,Sample_M) %>% tibble::column_to_rownames("SID") %>% pheatmap::pheatmap(scale="row",cluster_cols=F,border_color="NA",show_rownames=F)
temp.out.sel %>% select(SID,Sample_C,Sample_H,Sample_M) %>% tibble::column_to_rownames("SID") %>% pheatmap::pheatmap(scale="row",cluster_cols=F,border_color="NA",show_rownames=F)

temp.out %>% write.table(paste0("tmp_data/",TD,"/CellChat.detail.tsv"),col.names = T,row.names = F,sep="\t",quote=F)
temp.out.sel %>% write.table(paste0("tmp_data/",TD,"/CellChat.filtered.tsv"),col.names = T,row.names = F,sep="\t",quote=F)


#temp.out.sel  %>% filter(grepl("Fibro",SID) & grepl("Macro",SID)) %>% head(10)
#temp.out.sel  %>% filter(grepl("Fibro",target) & grepl("Macro",source)) %>% head(10)



#' create cell chat object
cellchat <- mergeCellChat(CellChat.ob, add.names = names(CellChat.ob))
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")
#netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)



#####


temp.out %>% tibble::column_to_rownames("SID") %>% pheatmap(scale="row",cluster_cols=F,border_color="NA",show_rownames=F)
set.seed(123)
ph=temp.out %>% tibble::column_to_rownames("SID") %>% pheatmap::pheatmap(scale="row",cluster_cols=F,border_color="NA",kmeans_k=8)
sel.SID <- ph$kmeans$cluster %>% as.data.frame() %>% setNames("Kcluster")%>% tibble::rownames_to_column("SID") %>% tbl_df() %>% filter(Kcluster %in% c("2",4)) %>% arrange(Kcluster) %>%filter(grepl("Fibro",SID)) %>% filter(!grepl("Neutrophil",SID)) %>% pull(SID)
temp.out %>% filter(SID %in% sel.SID) %>% tibble::column_to_rownames("SID") %>% pheatmap::pheatmap(scale="row",cluster_cols=F,border_color="NA",display_numbers=T)

#' check total interactions prob
temp.list <- list()
for (sa in c("Sample_C","Sample_H","Sample_M")) {
  temp.list[[sa]] <- CellChat.ob[[sa]]@net$weight %>% as.data.frame() %>% tibble::rownames_to_column("pair1") %>% tbl_df() %>% gather(pair2,VC,-pair1) %>% mutate(Sample=sa)
}
temp.out <- temp.list %>% do.call("bind_rows",.)%>% mutate(SID=paste(pair1,pair2,sep=":")) %>% select(-c(pair1,pair2))%>% spread(Sample,VC) %>% filter((Sample_C+Sample_H+Sample_M) >0 & (Sample_C!=Sample_H))

temp.out %>% tibble::column_to_rownames("SID") %>% pheatmap(scale="row",cluster_cols=F,border_color="NA",show_rownames=F)
set.seed(123)
ph=temp.out %>% tibble::column_to_rownames("SID") %>% pheatmap::pheatmap(scale="row",cluster_cols=F,border_color="NA",kmeans_k=8)
sel.SID <- ph$kmeans$cluster %>% as.data.frame() %>% setNames("Kcluster")%>% tibble::rownames_to_column("SID") %>% tbl_df() %>% filter(Kcluster %in% c("2",8)) %>% arrange(Kcluster) %>%filter(grepl("Fibro",SID)) %>% filter(!grepl("Neutrophil",SID)) %>% pull(SID)
temp.out %>% filter(SID %in% sel.SID) %>% tibble::column_to_rownames("SID") %>% pheatmap::pheatmap(scale="row",cluster_cols=F,border_color="NA",display_numbers=T)


# CellChat.ob = lapply(CellChat.ob,function(x) {netAnalysis_computeCentrality(x, slot.name = "netP") })
# cellchat=CellChat.ob[[1]]
# groupSize <- as.numeric(table(cellchat@idents))
# par(mfrow = c(1,2), xpd=TRUE)
# netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
# netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#Compare the total number of interactions and interaction strength
cowplot::plot_grid(
  compareInteractions(cellchat, show.legend = F, group = c(1,2,3)),
  compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight"),
  nrow=2,ncol=2
)


#Compare the number of interactions and interaction strength among different cell populations
cowplot::plot_grid(
  netVisual_diffInteraction(cellchat, weight.scale = T,comparison = c(1, 2)),
  netVisual_diffInteraction(cellchat, weight.scale = T,comparison = c(2, 3)),
  netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",comparison = c(1, 2)),
  netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight",comparison = c(2, 3)),
  nrow=2,ncol=2
)

# heatmaps for visulizaiton of interactions
temp.plot <- list()
temp.plot$p1 <- netVisual_heatmap(cellchat,comparison = c(1, 2))
temp.plot$p2 <- netVisual_heatmap(cellchat,comparison = c(2, 3))
temp.plot$p3 <- netVisual_heatmap(cellchat,comparison = c(1, 2), measure = "weight")
temp.plot$p4 <- netVisual_heatmap(cellchat,comparison = c(2, 3), measure = "weight")

temp.plot$p1+temp.plot$p2+temp.plot$p3+temp.plot$p4



#Differential number of interactions or interaction strength among different cell types
weight.max <- getMaxWeight(CellChat.ob, attribute = c("idents","count"))
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(CellChat.ob)) {
  netVisual_circle(CellChat.ob[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(CellChat.ob)[i]))
}

#Compare the major sources and targets in 2D space
num.link <- sapply(CellChat.ob, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(CellChat.ob)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(CellChat.ob[[i]], title = names(CellChat.ob)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)

### Identify signaling changes associated with one cell group
temp.plot <- list()
for (n in setdiff(unique(data.ob.umap$EML),c("Prolifer_SMC","Prolifer_Fibro"))) {
  temp.plot[[paste(n,"Sample_H_vs_Sample_C",sep=".")]] <-netAnalysis_signalingChanges_scatter(cellchat, idents.use = n,comparison=c(2,1))#+ggtitle(paste(n,"Sample_H_vs_Sample_C",sep="."))+theme(plot.title = element_text(hjust=0.5))
  temp.plot[[paste(n,"Sample_M_vs_Sample_H",sep=".")]] <-netAnalysis_signalingChanges_scatter(cellchat, idents.use = n,comparison=c(3,2))#+ggtitle(paste(n,"Sample_M_vs_Sample_H",sep="."))+theme(plot.title = element_text(hjust=0.5))
  temp.plot[[paste(n,"Sample_M_vs_Sample_C",sep=".")]] <-netAnalysis_signalingChanges_scatter(cellchat, idents.use = n,comparison=c(3,2))#+ggtitle(paste(n,"Sample_M_vs_Sample_C",sep="."))+theme(plot.title = element_text(hjust=0.5))
}
for (n in names(temp.plot)) {
  print(temp.plot[[n]])
}


## Visualization in 2D-space
# not working
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)
#> 2D visualization of signaling networks from datasets 1 2


temp.plot <- list()
temp.plot$p1 <- rankSimilarity(cellchat, type = "functional",comparison2 = c(1,2))
temp.plot$p2 <- rankSimilarity(cellchat, type = "functional",comparison2 = c(1,3))
temp.plot$p3 <- rankSimilarity(cellchat, type = "functional",comparison2 = c(2,3))
temp.plot$p1+temp.plot$p2+temp.plot$p3

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,comparison = c(1,2))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,comparison = c(1,2))
gg1 + gg2

temp.plot <- list()
pathway.union <- union(CellChat.ob$Sample_C@netP$pathways, CellChat.ob$Sample_H@netP$pathways) %>% union(CellChat.ob$Sample_M@netP$pathways)
for (sa in c("Sample_C","Sample_H","Sample_M")) {
  temp.plot[[sa]] <-  netAnalysis_signalingRole_heatmap(CellChat.ob[[sa]], pattern = "outgoing", signaling = pathway.union, title = sa, width = 5, height = 6)
}
draw(temp.plot[[1]]+temp.plot[[2]]+temp.plot[[3]],ht_gap = unit(0.5, "cm"))

temp.plot <- list()
pathway.union <- union(CellChat.ob$Sample_C@netP$pathways, CellChat.ob$Sample_H@netP$pathways) %>% union(CellChat.ob$Sample_M@netP$pathways)
for (sa in c("Sample_C","Sample_H","Sample_M")) {
  temp.plot[[sa]] <-  netAnalysis_signalingRole_heatmap(CellChat.ob[[sa]], pattern = "incoming", signaling = pathway.union, title = sa, width = 5, height = 6,color.heatmap = "GnBu")
}
draw(temp.plot[[1]]+temp.plot[[2]]+temp.plot[[3]],ht_gap = unit(0.5, "cm"))


temp.plot <- list()
pathway.union <- union(CellChat.ob$Sample_C@netP$pathways, CellChat.ob$Sample_H@netP$pathways) %>% union(CellChat.ob$Sample_M@netP$pathways)
for (sa in c("Sample_C","Sample_H","Sample_M")) {
  temp.plot[[sa]] <-  netAnalysis_signalingRole_heatmap(CellChat.ob[[sa]], pattern = "all", signaling = pathway.union, title = sa, width = 5, height = 6,color.heatmap = "GnBu")
}
draw(temp.plot[[1]]+temp.plot[[2]]+temp.plot[[3]],ht_gap = unit(0.5, "cm"))

#Identify the upgulated and down-regulated signaling ligand-receptor pairs
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11),  angle.x = 45,comparison = c(1, 2))



pathways.show <- c("WNT")#,"ncWNT") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(CellChat.ob)) {
  ht[[i]] <- netVisual_heatmap(CellChat.ob[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(CellChat.ob)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]]+ ht[[2]], ht_gap = unit(0.5, "cm"))

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(CellChat.ob)) {
  netVisual_aggregate(CellChat.ob[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(CellChat.ob)[i]))
}

plotGeneExpression(cellchat, signaling = "CXCL", split.by = "sample", colors.ggplot = T)