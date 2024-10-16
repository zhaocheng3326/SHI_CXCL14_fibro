#' ---
#' title: "CellChat analysis"
#' output: 
#'  html_document:
#'    code_folding: hide
#' ---



#R4.0

#' check the Interaction in Sample_C only
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
suppressMessages(library(ggalluvial))
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
#load(paste0("tmp_data/",TD,"/DEG.cluster.logFCtune.Rdata"),verbose=T)
#load(paste0("tmp_data/",TD,"/DEG.cluster.ov.Rdata"),verbose=T)
#load("tmp_data/gene.meta.Rdata",verbose=T)
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/whole.IT.UMAP.coord.rds")) %>% rows_update(readRDS(paste0("tmp_data/",TD,"/Fibryo.IT.UMAP.coord.rds")) %>% select(cell,EML,cluster_EML),by="cell") 
data.ob.umap <- data.ob.umap %>% mutate(cluster=EML) 
#sel.expG <- rownames(lognormExp.mBN)
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB 
which(CellChatDB.use[["interaction"]]$ligand == "H2-BI") # 1887
CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-1887,]
which(CellChatDB.use[["interaction"]]$ligand == "H2-Ea-ps") #1900
CellChatDB.use[["interaction"]] <- CellChatDB.use[["interaction"]][-1900,]

# temp <- CellChatDB.use$interaction %>% select(interaction_name,ligand,receptor) %>% tbl_df()
# temp.complex <-  CellChatDB.use$complex %>% tibble::rownames_to_column("ID") %>% tbl_df() %>% gather(subunit,gene,-ID) %>% filter(gene!="")
# temp <- temp %>% left_join(temp.complex %>% rename(ligand=ID),by="ligand") %>% left_join(temp.complex %>% rename(receptor=ID),by="receptor") 
# 
# Gene.IT <- temp %>% filter(is.na(gene.x) & is.na(gene.y)) %>% select(interaction_name,ligand,receptor) %>% bind_rows(temp %>% filter(!is.na(gene.x) & is.na(gene.y)) %>% mutate(ligand=gene.x) %>% select(interaction_name,ligand,receptor)) %>% bind_rows(temp %>% filter(is.na(gene.x) & !is.na(gene.y)) %>% mutate(receptor=gene.y) %>% select(interaction_name,ligand,receptor)) %>% bind_rows(temp %>% filter(!is.na(gene.x) & !is.na(gene.y)) %>% mutate(ligand=gene.x,receptor=gene.y) %>% select(interaction_name,ligand,receptor)) %>% unique() %>% inner_join(CellChatDB.use$interaction %>% select(interaction_name,pathway_name) %>% unique(),by="interaction_name")
# Gene.IT.short <- Gene.IT %>% select(interaction_name,pathway_name,ligand) %>% rename(gene=ligand) %>% mutate(type="ligand") %>% unique() %>% bind_rows(Gene.IT %>% select(interaction_name,pathway_name,receptor) %>% rename(gene=receptor) %>% mutate(type="receptor") %>% unique()) %>% unique()

CellChat.ob <-  readRDS(paste0("tmp_data/",TD,"/CellChat.ob.rds"))

sa="Sample_C"
cellchat <- CellChat.ob$Sample_C
netP <- subsetCommunication(cellchat,slot.name="netP",thresh=0.05) %>% tbl_df()%>% mutate(Sample=sa)
groupSize <- as.numeric(table(cellchat@idents))

#+fig.width=16,fig.height=16
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
#+fig.width=16,fig.height=16
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


#ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
#ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
#ht1 + ht2
#+fig.width=7,fig.height=7
netAnalysis_signalingRole_scatter(cellchat)

#' ##check which signals contributing most to outgoing or incoming signaling of certain cell groups.

ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",height=16,font.size=6)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",height=16,font.size=6)
#+fig.width=16,fig.height=10
ht1 
ht2

#' calculate the similaritys
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)


#selectK(cellchat, pattern = "outgoing")

outgoing.npattern <- 3
#+fig.width=15,fig.height=15
plot.new()
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = outgoing.npattern,height=12,font.size=6)
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")
#selectK(cellchat, pattern = "incoming")
incoming.npattern <- 3
plot.new()
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = incoming.npattern,height=12,font.size=6)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")

#+fig.width=9,fig.height=9
for (pathways.show in ( netP %>% group_by(pathway_name) %>% summarise(nPair=n_distinct(source,target)) %>% arrange(desc(nPair)) %>% filter(nPair > 10) %>% pull(pathway_name) %>% head(50))) {
  print(pathways.show)
  levels(cellchat@idents) 
  vertex.receiver = seq(3,11)
  
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 12, height = 4, font.size = 10)
  print(netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show)+ggtitle(pathways.show))
  print(netAnalysis_signalingRole_heatmap(cellchat, signaling = pathways.show)+ggtitle(pathways.show))
  
  
  #netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)  
  netVisual_aggregate(cellchat, signaling = pathways.show,vertex.size = groupSize)[[1]]
  #netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
  #netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
  netAnalysis_contribution(cellchat, signaling = pathways.show)
  netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  
  group.cellType <- c("Dendritic","Endo",rep("Fibro",9),"Kera",rep("Macro",4),"Mast","MelScwann","Neutrophil",rep("NK",2),"SMC","Macro") # grouping cell clusters into fibroblast, DC and TC cells
  names(group.cellType) <- levels(cellchat@idents)
  netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
  netAnalysis_contribution(cellchat, signaling = pathways.show)
  
  #' get the enriched LR pairs
  pairLR<- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
  LR.show <- pairLR[1,] # show one ligand-receptor pair
 
  netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)[[1]]
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")[[1]]
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")[[1]]
  
  
}



