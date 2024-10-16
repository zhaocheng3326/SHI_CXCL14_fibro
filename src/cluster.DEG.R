#' ---
#' title: "DE between Among Sample C, H , M "
#' output: 
#'  html_document:
#'    code_folding: hide
#' ---

#' ## DEG analysis 

#R_test
rm(list=ls())

condaENV <- "/home/chenzh/miniconda3/envs/R_test"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(pheatmap))
suppressMessages(library(Seurat))
suppressMessages(library(topGO))
suppressMessages(library(clusterProfiler))
#suppressMessages(library(STRINGdb))
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


#pt.cutoff <- 0.4
#' loading data
meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds"))
load("tmp_data/gene.meta.Rdata",verbose=T)
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/whole.IT.UMAP.coord.rds")) %>% rows_update(readRDS(paste0("tmp_data/",TD,"/Fibryo.IT.UMAP.coord.rds")) %>% select(cell,EML,cluster_EML),by="cell")
counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))
rownames(lognormExp.mBN) <- gsub("_","-",rownames(lognormExp.mBN))
rownames(counts.filter) <- gsub("_","-",rownames(counts.filter))
sel.expG <- rownames(lognormExp.mBN)

data.ob.umap <- data.ob.umap %>% mutate(cluster=EML) 

#' #### calculate DE between WT vs KO
savefile <- paste0("tmp_data/",TD,"/DEG.cluster.logFCtune.Rdata")
#DEG.cluster.logFCtune.withMT.Rdata

#x,data.ob.umap,data.deg
FunDEGCal <- function(temp.compair,temp.umap,data.deg.raw) {
  temp.merge.out <- list()
  for (x in c("Sample_H","Sample_M")) {
    for (y in c("Sample_C","Sample_H")) {
     
      if (x !=y) {
        c1=temp.umap %>% filter(cluster==temp.compair & sample==x) %>% pull(cell)
        c2=temp.umap %>% filter(cluster==temp.compair & sample==y) %>% pull(cell)
        
        data.deg.temp= subset(data.deg.raw,cell=c(c1,c2))
        G1G2.DEG <- suppressMessages(FindMarkers(data.deg.temp,ident.1=x,ident.2=y,test.use="MAST",assay="RNA",verbose=F,logfc.threshold=0.1,min.pct=0.15,pseudocount.use=0.1) %>% tibble::rownames_to_column(var="gene")) %>% tbl_df() #%>% left_join(Gene.pos %>% dplyr::select(V7,V9)%>% unique()%>% dplyr::rename(gene=V7,discription=V9),by="gene")
        #G1G2.DEG.psd <- suppressMessages(FindMarkers(data.deg.temp,ident.1="WT",ident.2="KO",test.use="MAST",assay="RNA",verbose=F,logfc.threshold=0.1,min.pct=0.15,pseudocount.use=0.1) %>% tibble::rownames_to_column(var="gene")) %>% tbl_df() %>% left_join(Gene.pos %>% dplyr::select(V7,description) %>% dplyr::rename(gene=V7),by="gene")
        # G1G2.DEG.ROC <- FindMarkers(data.deg.temp,ident.1="WT",ident.2="KO",test.use="roc",assay="RNA",verbose=F) %>% tibble::rownames_to_column(var="gene") %>% tbl_df() 
        
        G1.sig.up <- G1G2.DEG %>% filter(avg_log2FC >0.25 & p_val_adj <0.05)  # R4.0
        G2.sig.up <- G1G2.DEG %>% filter(avg_log2FC < -0.25 & p_val_adj <0.05) # R4.0
        #G1.sig.up <- G1G2.DEG %>% filter(avg_logFC >0.25 & p_val_adj <0.05)  # R4.0
        #G2.sig.up <- G1G2.DEG %>% filter(avg_logFC < -0.25 & p_val_adj <0.05) # R4.0
        
        G1G2.DEG.sig <- G1.sig.up %>% bind_rows(G2.sig.up)
        
        temp.up.ID <- G1.sig.up$gene
        temp.down.ID<- G2.sig.up$gene
        
        temp.out <- list()
        
        temp.out[["DEG.all.result"]] <- G1G2.DEG
        temp.out[["DEG.result"]] <- G1G2.DEG.sig
        temp.out[["DEG.result.up"]] <-  G1.sig.up
        temp.out[["DEG.result.down"]] <-  G2.sig.up
        
        if (nrow(G1.sig.up) > 0) {
          temp.out[["DEG.up.GO.result"]] <- topGO_enrichment2(temp.up.ID,ALL_gene,0.05,geneID2GO)
        }else{
          DEG.results[[temp.compair]] [["DEG.up.GO.result"]] <- NA
        }
        if (nrow(G2.sig.up) > 0) {
          temp.out[["DEG.down.GO.result"]] <- topGO_enrichment2(temp.down.ID,ALL_gene,0.05,geneID2GO)
        }else{
          temp.out[["DEG.down.GO.result"]] <- NA
        }
        temp.merge.out[[paste(x,y,sep="_vs_")]] <- temp.out
      }
      
    }
  }
  return( temp.merge.out)
}

  
  

if (file.exists(savefile)){
  load(savefile,verbose = T)
}else{
  
  
  
  #' ### Loading gene annotation
  ALL_gene <- rownames(counts.filter)
  geneID2GO=inverseList(readMappings(file = "big_doc/rats.GO2geneID_ALL.refseq.map"))
  names(geneID2GO)= gsub("_","-",names(geneID2GO))
  
  ## build a dummy seurat object.
  temp.sel.expG <- rownames(lognormExp.mBN )
  temp.M <- meta.filter
  data.deg <-   CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)
  Idents(data.deg) <- factor(data.deg@meta.data$sample)
  data.deg@assays$RNA@data <- as.matrix(lognormExp.mBN[rownames(data.deg@assays$RNA@data),colnames(data.deg@assays$RNA@data)])
  
  
  #' ## Number of cells in different sets
  DEG.stat=matrix(nrow=2,ncol=3*length(unique(data.ob.umap$cluster)))
  rownames(DEG.stat) <- c("up_regulated","down_regulated")
  colnames(DEG.stat) <- paste(rep(unique(data.ob.umap$cluster),each=3),c("Sample_H_vs_Sample_C","Sample_M_vs_Sample_C","Sample_M_vs_Sample_H"),sep=".")
  
  rm(counts.filter)
  DEG.results <- list()
  for (n in unique(data.ob.umap$cluster)) {
    DEG.results[[n]] <- n #DEG.results <- vector("list", ncol(DEG.stat))
  }
  
  # DEG.results <- foreach (x=DEG.results,n=names(DEG.results),.combine=c) %dopar% {
  #   rv=list()
  #   rv[[n]]=FunDEGCal(x,data.ob.umap,data.deg ) ## paralelly running #mast
  #   rv
  # }
  for (n in names(DEG.results)) {
    print(n)
    DEG.results[[n]] <- FunDEGCal(n,data.ob.umap,data.deg )
  }
  # transform the DEG.results type
  DEG.results.trans <- DEG.results
  DEG.results <- list()
  for (n in names(DEG.results.trans)) {
    for (m in names(DEG.results.trans[[n]]) ) {
      DEG.results[[paste(n,m,sep=".")]] <-  DEG.results.trans[[n]][[m]]
    }
  }
  
  for (n in names(DEG.results)) {
      DEG.stat["up_regulated",n] <- nrow(DEG.results[[n]][["DEG.result.up"]])
      DEG.stat["down_regulated",n]  <- nrow(DEG.results[[n]][["DEG.result.down"]])
  }
  
  #' loading path enrichment analysis
  PATH_list <- list()
  for (tp in c("KEGG_2019","WikiPathways_2019","human_Reactome_2016_tm")) {
    temp <- read.delim(paste("~/Genome_new/Mouse/",tp,"_Mouse.mod.tsv",sep=""),head=F,stringsAsFactors=F)%>% mutate(V2=stringr::str_to_title(V2))
    colnames(temp) <- c("TERM","GENE")
    #rownames(temp) <- temp$TERM
    PATH_list[[tp]] <- temp
  }
  genename.bg <-  ALL_gene ### background genes
  
  pt.results <- list()
  for(n in colnames(DEG.stat)) {
    print(n)
    temp.compair <- n
    pt.results[[temp.compair]] <- list()
    
    temp.DEG <- DEG.results[[temp.compair]]$DEG.result
    up.DEG <- DEG.results[[temp.compair]]$DEG.result.up$gene
    down.DEG <- DEG.results[[temp.compair]]$DEG.result.down$gene
    for (tp in c("KEGG_2019","WikiPathways_2019","human_Reactome_2016_tm")) {
      pt.results[[temp.compair]][[tp]]  <- list()
      
      suppressMessages(up.DEG.pt <- enricher(up.DEG, pvalueCutoff = 0.05, pAdjustMethod = "none",genename.bg ,minGSSize = 2,maxGSSize = 1000,qvalueCutoff = 1, TERM2GENE=PATH_list[[tp]]))
      suppressMessages(down.DEG.pt <- enricher(down.DEG, pvalueCutoff = 0.05, pAdjustMethod = "none",genename.bg ,minGSSize = 5,maxGSSize = 1000,qvalueCutoff = 1, TERM2GENE=PATH_list[[tp]]))
      
      pt.results[[temp.compair]][[tp]]$up <- up.DEG.pt
      pt.results[[temp.compair]][[tp]]$down <- down.DEG.pt
    }
  }
  
  # stingDB.results <- list()
  # string_db <- STRINGdb$new( version="11", species=10090,score_threshold=400, input_directory="/home/chenzh/Genome/mouse/")
  # for(n in colnames(DEG.stat)) {
  #   print(n)
  #   temp.compair <- n
  #   stingDB.results [[temp.compair]] <- list()
  # 
  #   temp.DEG <- DEG.results[[temp.compair]]$DEG.result
  #   up.DEG <- DEG.results[[temp.compair]]$DEG.result.up %>% select(gene,avg_log2FC,p_val_adj) %>% as.data.frame()
  #   down.DEG <- DEG.results[[temp.compair]]$DEG.result.down%>% select(gene,avg_log2FC,p_val_adj)%>% as.data.frame()
  #   
  #   st.mapped <- list()
  #   st.mapped$WT <- string_db$map( up.DEG , "gene", removeUnmappedRows = TRUE )
  #   st.mapped$KO <- string_db$map( down.DEG , "gene", removeUnmappedRows = TRUE )
  #   stingDB.results[[temp.compair]] <- st.mapped 
  # }
  #   
  #save(DEG.results,DEG.stat,pt.results, stingDB.results,file=savefile)
  save(DEG.results,DEG.stat,pt.results,file=savefile)
}

#' #### visulization
#' ## Number of DEGs (including non-PC genes)
#+ fig.width=27,fig.height=9
par(mar=c(15,5,5,3))
FunDegNumPlot(DEG.stat,"Number of DEGs")

#' ## Heatmap of DEGs
#+ fig.width=9,fig.height=9

for(n in colnames(DEG.stat)) {
  print(n)
  temp.compair <- n
  temp.cluster <- unlist(strsplit(n,"\\.Sample"))[[1]] 
  temp.s1 <-  unlist(gsub(paste0(temp.cluster,"."),"",n) %>% strsplit("_vs_"))[[1]]
  temp.s2 <-  unlist(gsub(paste0(temp.cluster,"."),"",n) %>% strsplit("_vs_"))[[2]]
  c1 <- data.ob.umap %>% filter(cluster==temp.cluster & sample==temp.s1) %>% FunMaSF(500) %>% pull(cell)
  c2 <-  data.ob.umap %>% filter(cluster==temp.cluster & sample==temp.s2)%>% FunMaSF(500) %>% pull(cell)
  cells.use <- c(c1,c2)
  
  
  up.genes.use <- DEG.results[[temp.compair]]$DEG.result.up$gene %>% head(500)
  down.genes.use <- DEG.results[[temp.compair]]$DEG.result.down$gene%>% head(500)
  
  colorsR <- data.ob.umap %>% filter(cell %in% c(c1)) %>% bind_rows(data.ob.umap %>% filter(cell %in% c(c2))) %>% select(cell,sample,cluster)  %>% tibble::column_to_rownames("cell")
  
  if (nrow(DEG.results[[temp.compair]]$DEG.result) != 0) {
    print(paste("Heatmap of gene expression in",temp.compair))
    
    if (length(up.genes.use) >0) {
      if (length(down.genes.use) >0) {
        pheatmap(FunPreheatmapNoLog(lognormExp.mBN,c(up.genes.use,down.genes.use),c(c1,c2)),cluster_rows=F,cluster_cols=F,border_color = "NA",colorRampPalette(c("royalblue3","white","firebrick4"))(50),gaps_col=rep(length(c1),5),gaps_row = rep(length(up.genes.use),5),show_rownames = F,show_colnames = F,main=temp.compair,annotation_col=colorsR )
      }else{
        pheatmap(FunPreheatmapNoLog(lognormExp.mBN,c(up.genes.use,down.genes.use),c(c1,c2)),cluster_rows=F,cluster_cols=F,border_color = "NA",colorRampPalette(c("royalblue3","white","firebrick4"))(50),gaps_col=rep(length(c1),5),show_rownames = F,show_colnames = F,main=temp.compair,annotation_col=colorsR )
      }
    }else{
      pheatmap(FunPreheatmapNoLog(lognormExp.mBN,c(up.genes.use,down.genes.use),c(c1,c2)),cluster_rows=F,cluster_cols=F,border_color = "NA",colorRampPalette(c("royalblue3","white","firebrick4"))(50),gaps_col=rep(length(c1),5),show_rownames = F,show_colnames = F,main=temp.compair,annotation_col=colorsR ) 
    }
  }else{
    print(paste("No DEGs in",temp.compair))
  }
}


#' #### GO enrichment analysis of DEGs
#+ fig.width=10,fig.height=10
#
par(mar=c(8,5,3,5))
for (n in c(1:length(colnames(DEG.stat)))) {
  temp.compair <- colnames(DEG.stat)[n]
  print(temp.compair)
  if (is.list(DEG.results[[temp.compair]]$DEG.up.GO.result)) {
    print(paste("GO of up regulated genes in",temp.compair))
    FunGOenrichPlot(DEG.results[[temp.compair]]$DEG.up.GO.result$top10,"upRe DEGs")
  }else{
    print(paste("No up regulated genes in",temp.compair))
    print("")
  }
  if (is.list(DEG.results[[temp.compair]]$DEG.down.GO.result)) {
    print(paste("GO of down regulated genes in",temp.compair))
    FunGOenrichPlot(DEG.results[[temp.compair]]$DEG.down.GO.result$top10,"downRe DEGs")
  }else{
    print(paste("No up regulated genes in",temp.compair))
    print("")
  }
}

#' #### Pathway enrichment analysis.

#+ fig.width=15,fig.height=9
for (n in c(1:length(colnames(DEG.stat)))) {
  for (tp in c("KEGG_2019","WikiPathways_2019","human_Reactome_2016_tm")) {
    temp.compair <- colnames(DEG.stat)[n]
    up.DEG.pt <- pt.results[[temp.compair]][[tp]]$up
    down.DEG.pt <- pt.results[[temp.compair]][[tp]]$down
    if (! is.null(up.DEG.pt)) {
      up.DEG.pt.sig=subset(as.data.frame(up.DEG.pt),up.DEG.pt$pvalue <0.05)[,c(1,3,4,5,8,9)]
      print(paste("up-regulated DEGs in",temp.compair,paste("(",tp,")",sep="")))
      print(dotplot(up.DEG.pt,color="pvalue",title=paste("up-regulated DEGs in",temp.compair,paste("(",tp,")",sep=""))))
    }else{
      print(paste("No significant item in up-regulated genes"))
    }
    if (! is.null(down.DEG.pt)) {
      down.DEG.pt.sig=subset(as.data.frame(down.DEG.pt),down.DEG.pt$pvalue <0.05)[,c(1,3,4,5,8,9)]
      print(paste("down-regulated DEGs in",temp.compair,paste("(",tp,")",sep="")))
      print(dotplot(down.DEG.pt,color="pvalue",title=paste("down-regulated DEGs in",temp.compair,paste("(",tp,")",sep=""))))
    }else{
      print(paste("No significant item in down-regulated genes"))
    }
  }
}


#' #### output results
savedir=paste0("tmp_data/",TD,"/DEG/")
for (n in c(1:length(colnames(DEG.stat)))) {
  temp.compair <- colnames(DEG.stat)[n]
  temp.out <- DEG.results[[temp.compair]]
  if (nrow(temp.out$DEG.result.up) > 0) {
    write.table(temp.out$DEG.result.up ,file=paste(savedir,temp.compair,".upDEG",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
    write.table(temp.out[["DEG.up.GO.result"]][["out"]],file=paste(savedir,temp.compair,".upDEG.GO",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
  }else{
    write.table("None",file=paste(savedir,temp.compair,".upDEG",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
    write.table("None",file=paste(savedir,temp.compair,".upDEG.GO",sep=""),col.names = F,row.names = F,sep="\t",quote=F)
  }
  if (nrow(temp.out$DEG.result.down) > 0) {
    write.table(temp.out$DEG.result.down,file=paste(savedir,temp.compair,".downDEG",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
    write.table(temp.out[["DEG.down.GO.result"]][["out"]],file=paste(savedir,temp.compair,".downDEG.GO",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
  }else{
    write.table("None",file=paste(savedir,temp.compair,".downDEG",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
    write.table("None",file=paste(savedir,temp.compair,".downDEG.GO",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
  }
  for (tp in c("KEGG_2019","WikiPathways_2019","human_Reactome_2016_tm")) {
    write.table(pt.results[[temp.compair]][[tp]]$up %>% tbl_df(),file=paste(savedir,temp.compair,".upDEG.",tp,".out",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
    write.table(pt.results[[temp.compair]][[tp]]$down %>% tbl_df(),file=paste(savedir,temp.compair,".downDEG.",tp,".out",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
  }
}


#' #### output results
savedir=paste0("tmp_data/",TD,"/DEG_sel/")
for (temp.compair in c("M0Macro","M1Macro","M2Macro","Spp1_hi_Macro")) {
  
  temp.out <- DEG.results[[paste0(temp.compair,".Sample_M_vs_Sample_C")]]
  temp.out2 <- DEG.results[[paste0(temp.compair,".Sample_H_vs_Sample_C")]]
  
    write.table(temp.out$DEG.result.up %>% mutate(H_vs_C.type=ifelse(gene %in% temp.out2$DEG.result.up$gene,"both","M_unique"))%>% left_join(Gene.pos %>% mutate(gene=V7,annotation=V9) %>% select(gene,annotation),by="gene"),file=paste(savedir,temp.compair,".upDEG",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
    write.table(temp.out[["DEG.up.GO.result"]][["out"]],file=paste(savedir,temp.compair,".upDEG.GO",sep=""),col.names = T,row.names = F,sep="\t",quote=F)

    write.table(temp.out$DEG.result.down %>% mutate(H_vs_C.type=ifelse(gene %in% temp.out2$DEG.result.down$gene,"both","M_unique")) %>% left_join(Gene.pos %>% mutate(gene=V7,annotation=V9) %>% select(gene,annotation),by="gene") ,file=paste(savedir,temp.compair,".downDEG",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
    write.table(temp.out[["DEG.down.GO.result"]][["out"]],file=paste(savedir,temp.compair,".downDEG.GO",sep=""),col.names = T,row.names = F,sep="\t",quote=F)
}




#sessionInfo()