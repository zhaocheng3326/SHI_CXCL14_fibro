#' ---
#' title: "QC and assign cell cycle"
#' output: 
#'  html_document:
#'    code_folding: hide
#' ---




# ### Loading R library

#R3.6
rm(list=ls())
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(pheatmap))
suppressMessages(library(Seurat))
suppressMessages(library(scran))
suppressMessages(library(batchelor))

suppressMessages(library(DoubletDecon))
suppressMessages(library(DoubletFinder))
suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)

# working directory
DIR <- "/home/chenzh/My_project/SHI_Glu"
setwd(DIR)

# Loading R functions
source("~/PC/R_code/functions.R")
source("~/PC/SnkM/SgCell.R")
rename <- dplyr::rename

options(digits = 4)
options(future.globals.maxSize= 3001289600)

qc.nGene.min <- 750
qc.nGene.max <- 6000
qc.mt.perc <- 0.15 #FunMAD(2)
qc.rn45s.perc <- 0.05
TD="Aug_2022"
if (file.exists(paste0("tmp_data/",TD,"/meta.filter.rds"))) {
  load("tmp_data/all.counts.meta.Rdata",verbose=T)
  meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds"))
}else{
  #' loading data
  
  load("tmp_data/all.counts.meta.Rdata",verbose=T)
  load("tmp_data/gene.meta.Rdata",verbose = T)
  
  db_check=TRUE
  #' for cell cycle
  mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran")) %>% lapply(function(x){x %>% tbl_df() %>% inner_join(rats.mouse %>% select(external_gene_name,mmusculus_homolog_ensembl_gene)  %>% rename(first=mmusculus_homolog_ensembl_gene),by="first") %>% mutate(first=external_gene_name)  %>% select(first,second) %>% inner_join(rats.mouse %>% select(external_gene_name,mmusculus_homolog_ensembl_gene)  %>% rename(second=mmusculus_homolog_ensembl_gene),by="second") %>% mutate(second=external_gene_name)  %>% select(first,second) %>% unique() %>% as.data.frame()%>% return()})
  
  meta.filter <- meta.all %>% filter( nGene < qc.nGene.max & mt.perc < qc.mt.perc & nGene > qc.nGene.min) 
  
  #counts.filter <- counts.all[setdiff(rownames(counts.all), mt.gene),meta.filter$cell]
  rownames(counts.all) %>% setdiff(Gene.pos$V7) ## some weird genes
  
  #counts.filter <- counts.all[rownames(counts.all) %>% intersect(Gene.pos$V7),meta.filter$cell] ## stil include the mt genes
  counts.filter <- counts.all[rownames(counts.all) %>% intersect(Gene.pos$V7) %>% setdiff(mt.gene),meta.filter$cell] ## stil include the mt genes
  
  #' get the expressed genes ( expressed in at least 1 datatsets)
  expG.set <- list()
  for (b in unique(meta.filter$sample  %>% unique() %>% as.vector())) { 
    temp.cell <- meta.filter %>% filter(sample==b) %>% pull(cell)
    expG.set[[b]] <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]
  }
  sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()
  
  
  if (file.exists(paste0("tmp_data/",TD,"/meta.filter.dbf.out.rds"))) {
    dbf.out <- readRDS(paste0("tmp_data/",TD,"/meta.filter.dbf.out.rds"))
  }else{
    expG.set <- list()
    for (b in unique(meta.filter$sample  %>% unique() %>% as.vector())) { 
      temp.cell <- meta.filter %>% filter(sample==b) %>% pull(cell)
      expG.set[[b]] <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]
    }
    sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()
    
    #
    dbf.list <- list()
    for (sa in c("Sample_C","Sample_H","Sample_M")) {
      #' check the Doublets ratio ( using hs.mt) to do following analysis
      data.double.stat <- NULL
      temp.sel.expG <- sel.expG
      temp.M <- meta.filter  %>% filter(sample==sa)
      nGene=2000;npcs=25
      data.temp <-   CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)%>% FindVariableFeatures( selection.method = "vst", nfeatures = nGene, verbose = FALSE) %>% ScaleData(verbose=F,vars.to.regress=c("nGene"))%>% RunPCA(verbose=F) %>% RunUMAP( dims = 1:npcs, verbose = FALSE) %>% FindNeighbors( dims = 1:npcs, verbose = FALSE) %>% FindClusters( verbose = FALSE,reso=0.2)
      #nExp <- round(ncol(data.temp) * 0.04)
      #data.temp <- doubletFinder_v3(data.temp, PCs = 1:npcs, pN = 0.25, pK = 0.09, nExp = nExp)
      #dbf.list[[sa]] <- data.temp@meta.data[,ncol(data.temp@meta.data),drop=F]%>% setNames("dbfd") %>% tibble::rownames_to_column("cell")%>% tbl_df()
      ## pK Identification (ground-truth) ------------------------------------------------------------------------------------------
      sweep.res.list_data<- paramSweep_v3(data.temp, PCs = 1:npcs, sct = T,num.cores=20) #long time to run
      sweep.stats_data<- summarizeSweep(sweep.res.list_data, GT = F)
      bcmvn_data<- find.pK(sweep.stats_data)
      
      ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
      homotypic.prop <- modelHomotypic(Idents(data.temp))           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
      nExp_poi <- round(0.05*length(colnames(data.temp)))  ## Assuming 0.05% doublet formation rate - tailor for your dataset
      nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
      
      ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
      opt.pk <- (bcmvn_data %>% arrange(desc(BCmetric)) %>% pull(pK) %>% head(1)) %>% as.vector() %>% as.numeric()
      opt.pn <- sweep.stats_data %>% filter(pK==opt.pk) %>% arrange(desc(BCreal)) %>% pull(pN) %>% head(1) %>% as.vector()%>% as.numeric()
      data.temp <- doubletFinder_v3(data.temp, PCs = 1:npcs, pN = opt.pn, pK = opt.pk, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = T)
      dbf.list[[sa]] <- data.temp@meta.data[,ncol(data.temp@meta.data),drop=F]%>% setNames("dbfd")%>% tibble::rownames_to_column("cell") %>% tbl_df()
    }
    dbf.out <- dbf.list %>% do.call("bind_rows",.) %>% tbl_df()
    saveRDS(dbf.out,paste0("tmp_data/",TD,"/meta.filter.dbf.out.rds"))
  }
  meta.filter <- meta.filter %>% left_join(dbf.out,by="cell") %>% filter(dbfd=="Singlet")
  counts.filter <- counts.filter[,meta.filter$cell] 
  
  expG.set <- list()
  for (b in unique(meta.filter$sample  %>% unique() %>% as.vector())) { 
    temp.cell <- meta.filter %>% filter(sample==b) %>% pull(cell)
    expG.set[[b]] <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]
  }
  sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()
  
  sce.ob <- list()
  sce.ob.phase <- list()
  for (b in unique(meta.filter$sample  %>% unique() %>% as.vector())) { 
    print(b)
    temp.M <- meta.filter %>% filter(sample==b) 
    temp.sce <-  SingleCellExperiment(list(counts=as.matrix(counts.filter[sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors() 
    #temp.assignments <- cyclone(temp.sce, mm.pairs)
    #sce.ob.phase[[b]] <- data.frame(cell=temp.M$cell,phases=temp.assignments$phases) %>% tbl_df() %>% mutate_all(as.vector)
    sce.ob[[b]] <- temp.sce
  }
  #meta.filter <- do.call("bind_rows",sce.ob.phase) %>% full_join(meta.filter,by="cell")
  metta.filter <- meta.filter %>% mutate(phases=NA)
  
  mBN.sce.ob <- multiBatchNorm(sce.ob$Sample_C,sce.ob$Sample_H,sce.ob$Sample_M)
  lognormExp.mBN<- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)
  
  
  saveRDS(counts.filter,file=paste0("tmp_data/",TD,"/counts.filter.rds"))
  saveRDS(meta.filter,file=paste0("tmp_data/",TD,"/meta.filter.rds"))
  saveRDS(lognormExp.mBN,file=paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))
 
  
  
}

#' 

#' #### check the general distribution
meta.all %>% ggplot+geom_histogram(mapping=aes(x=nGene,fill=sample),bins=100)+geom_vline(xintercept=qc.nGene.min)+geom_vline(xintercept=qc.nGene.max)+facet_wrap(.~sample,nrow=2)+ylab("No. of cells")
meta.all %>% ggplot+geom_boxplot(mapping=aes(x=sample,y=nGene,fill=sample))
meta.all %>% ggplot+geom_histogram(mapping=aes(x=mt.perc,fill=sample),bins=100)+geom_vline(xintercept=qc.mt.perc)+facet_wrap(.~sample,nrow=2)+ylab("No. of cells")
meta.all %>% ggplot+geom_point(mapping=aes(x=nGene,y=mt.perc,col=sample))+facet_wrap(.~sample,nrow=2)+ylab("No. of cells")

#FunMAD
#' #### number of cells
meta.all %>% group_by(sample) %>% summarise(BeforeQCnCell=n()) %>% inner_join(meta.filter  %>% group_by(sample) %>% summarise(AfterQCnCell=n()),by="sample") %>% gather(QC,nCell,-sample) %>% ggplot+geom_bar(mapping=aes(x=sample,fill=QC,y=nCell),stat="identity",position="dodge")

#' check the cell number during QC
table(meta.filter$sample)
table(meta.all$sample)