#' ---
#' title: "GSEA analysis of each comparison"
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
  #library(topGO)
  library(clusterProfiler)
  library(ggvenn)
  #library(scran)
  #library(batchelor)
  #library(SeuratWrappers)
  library(DOSE)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(GSVA)
  library(GSEABase)
  library(msigdbr)
})

#' option 
set.seed(1)
DIR=paste0(base_dir,"/My_project/SHI_Glu")
setwd(DIR)
knitr::opts_knit$set(root.dir=DIR)
para_cal <- TRUE
rewrite=FALSE


rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter


TD="Aug_2022"


#pt.cutoff <- 0.4
#' loading data
meta.filter <- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds"))
load("tmp_data/gene.meta.Rdata",verbose=T)
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/whole.IT.UMAP.coord.rds")) %>% rows_update(readRDS(paste0("tmp_data/",TD,"/Fibryo.IT.UMAP.coord.rds")) %>% select(cell,EML,cluster_EML),by="cell")
data.ob.umap <- data.ob.umap %>% mutate(cluster=EML) 
genename.bg <- ALL_gene <- rownames(readRDS(paste0("tmp_data/",TD,"/counts.filter.rds")))


#load DEG results
load(paste0("tmp_data/",TD,"/DEG.cluster.logFCtune.Rdata"),verbose=T)



gsea.results <- list()
pt.results <- list()
if (file.exists(paste0("tmp_data/",TD,"/DEG.cluster.MSigR.pt.results.rds"))) {
  pt.results <- readRDS(paste0("tmp_data/",TD,"/DEG.cluster.MSigR.pt.results.rds"))
  gsea.results <- readRDS(paste0("tmp_data/",TD,"/DEG.cluster.MSigR.gsea.results.rds"))
}else{
  GS_db = readRDS("big_doc/Rats.msigdbr.rds")
  for (m in c("M0Macro","M1Macro","M2Macro","Spp1_hi_Macro")) {
    for (n in c("Sample_M_vs_Sample_C","Sample_H_vs_Sample_C","Sample_M_vs_Sample_H")) {
      temp.compair <- paste(m,n,sep=".")
      temp.fullGene <- DEG.results[[temp.compair]]$DEG.all.result  %>% arrange(desc(avg_log2FC)) ### mast
      temp.fc <-  temp.fullGene$avg_log2FC
      names(temp.fc) <- temp.fullGene$gene
      
      for (tp in c("CP:BIOCARTA", "CP:KEGG", "CP:REACTOME", "CP:WIKIPATHWAYS", "TFT:GTRD")) {
        print(paste(m,n,tp))
        temp.path <- GS_db %>% filter(gs_subcat == tp) %>% mutate(TERM=gs_name,GENE=gene_symbol) %>% select(TERM,GENE) %>% unique() %>% as.data.frame()
        suppressMessages(pt.results[[temp.compair]][[tp]]$up  <- enricher(DEG.results[[temp.compair]]$DEG.result.up$gene, pvalueCutoff = 1, pAdjustMethod = "BH",genename.bg ,minGSSize = 2,maxGSSize = 1000,qvalueCutoff = 1, TERM2GENE=temp.path))
        suppressMessages(pt.results[[temp.compair]][[tp]]$down  <- enricher(DEG.results[[temp.compair]]$DEG.result.down$gene, pvalueCutoff = 1, pAdjustMethod = "BH",genename.bg ,minGSSize = 2,maxGSSize = 1000,qvalueCutoff = 1, TERM2GENE=temp.path))
        if (tp=="CP:WIKIPATHWAYS") {
          suppressMessages(gsea.results[[temp.compair]][[tp]] <-  GSEA(temp.fc,pvalueCutoff = 1, pAdjustMethod = "BH",minGSSize = 2,maxGSSize = 1000, TERM2GENE=temp.path))
        }
      }
    } 
  }
  saveRDS(pt.results,paste0("tmp_data/",TD,"/DEG.cluster.MSigR.pt.results.rds"))
  saveRDS(gsea.results,paste0("tmp_data/",TD,"/DEG.cluster.MSigR.gsea.results.rds"))
  
  temp <- GS_db %>% filter(gs_subcat == "CP:WIKIPATHWAYS") %>% filter(gs_name=="WP_OVERVIEW_OF_PROINFLAMMATORY_AND_PROFIBROTIC_MEDIATORS") %>% select(human_gene_symbol,gene_symbol) %>% unique()
  inflam.list <- list()
  inflam.list$anti <- temp %>% filter(human_gene_symbol %in% read.delim("doc/human.antiinflam.gene.id",head=F)$V1) %>% pull(gene_symbol) %>% intersect(ALL_gene)
  inflam.list$pro <- temp %>% filter(human_gene_symbol %in% read.delim("doc/human.proinflam.gene.id",head=F)$V1) %>% pull(gene_symbol) %>% intersect(ALL_gene)
  saveRDS(inflam.list,file="doc/inflam.geneid.rds")
  
  
}

#" loading overlapped GENE results
load(paste0("tmp_data/",TD,"/DEG.cluster.ov.Rdata"),verbose = T)
gsea.ov.results <- list()
pt.ov.results <- list()
if (file.exists(paste0("tmp_data/",TD,"/DEG.ov.MSigR.pt.results.rds"))) {
  pt.ov.results <- readRDS(paste0("tmp_data/",TD,"/DEG.ov.MSigR.pt.results.rds"))
  gsea.ov.results <- readRDS(paste0("tmp_data/",TD,"/DEG.ov.MSigR.gsea.results.rds"))
}else{
  GS_db = readRDS("big_doc/Rats.msigdbr.rds")
  for (m in c("M0Macro","M1Macro","M2Macro","Spp1_hi_Macro")) {
    pt.ov.results[[m]] <- list()
    gsea.ov.results[[m]] <- list()
    temp.DEG.set <- list()
    temp.DEG.set$MU <- c(DEG.ov[[m]]$out[["MU.upDEG"]],DEG.ov[[m]]$out[["MU.downDEG"]])
    temp.DEG.set$ov <- c(DEG.ov[[m]]$out[["ov.upDEG"]],DEG.ov[[m]]$out[["ov.downDEG"]])
    for (n in c("MU","ov")) {
     
     # temp.DEG <-temp.DEG.set[[n]]
     # temp.fc <-  (DEG.results[[paste0(m,".Sample_M_vs_Sample_C")]]$DEG.result %>% tibble::column_to_rownames("gene"))[temp.DEG,"avg_log2FC"]
     # names(temp.fc) <- temp.DEG
     # temp.fc <- temp.fc[rev(order(temp.fc))]
      
      for (tp in c("CP:BIOCARTA", "CP:KEGG", "CP:REACTOME", "CP:WIKIPATHWAYS", "TFT:GTRD")) {
        print(paste(m,n,tp))
        temp.path <- GS_db %>% filter(gs_subcat == tp) %>% mutate(TERM=gs_name,GENE=gene_symbol) %>% select(TERM,GENE) %>% unique() %>% as.data.frame()
        suppressMessages(pt.ov.results[[m]][[n]][[paste("up",tp)]]  <- enricher(DEG.ov[[m]]$out[[paste0(n,".upDEG")]], pvalueCutoff = 1, pAdjustMethod = "BH",genename.bg ,minGSSize = 2,maxGSSize = 1000,qvalueCutoff = 1, TERM2GENE=temp.path))
        suppressMessages(pt.ov.results[[m]][[n]][[paste("down",tp)]]  <- enricher(DEG.ov[[m]]$out[[paste0(n,".downDEG")]], pvalueCutoff = 1, pAdjustMethod = "BH",genename.bg ,minGSSize = 2,maxGSSize = 1000,qvalueCutoff = 1, TERM2GENE=temp.path))
        #suppressMessages(gsea.ov.results[[m]][[n]] <-  GSEA(temp.fc,pvalueCutoff = 0.5, pAdjustMethod = "BH",minGSSize = 2,maxGSSize = 1000, TERM2GENE=temp.path))
      }
    } 
  }
  saveRDS(pt.ov.results,paste0("tmp_data/",TD,"/DEG.ov.MSigR.pt.results.rds"))
  saveRDS(gsea.ov.results,paste0("tmp_data/",TD,"/DEG.ov.MSigR.gsea.results.rds"))
}



#CP:WIKIPATHWAYS


#pt.results$M1Macro.Sample_M_vs_Sample_C[["CP:WIKIPATHWAYS"]]$down %>% as.data.frame() %>% as_tibble %>% filter(qvalue < 0.05) %>% pull(ID)
#pt.results$M1Macro.Sample_M_vs_Sample_C[["CP:BIOCARTA"]]$down %>% as.data.frame() %>% as_tibble %>% filter(pvalue < 0.05) %>% pull(ID)
#pt.results$M1Macro.Sample_M_vs_Sample_C[["CP:KEGG"]]$down %>% as.data.frame() %>% as_tibble %>% filter(pvalue < 0.05) %>% pull(ID)

pt.results$M1Macro.Sample_M_vs_Sample_C[["CP:BIOCARTA"]]$down %>% as.data.frame() %>% as_tibble %>% filter(pvalue < 0.05) %>% pull(ID)


gsea.results$Spp1_hi_Macro.Sample_M_vs_Sample_C[["CP:WIKIPATHWAYS"]] %>% as.data.frame() %>% View()

temp.gse <- gsea.results$Spp1_hi_Macro.Sample_M_vs_Sample_C[["CP:WIKIPATHWAYS"]]
temp.n <- temp.gse %>% as.data.frame() %>% tbl_df() %>%  tibble::rowid_to_column("n") %>% filter(ID=="WP_TGFBETA_RECEPTOR_SIGNALING") %>% pull(n)
gseaplot(temp.gse, by = "all", title = temp.gse$Description[temp.n], geneSetID = 1)





temp.list <- list()
for (n in c("M0Macro","M1Macro","M2Macro","Spp1_hi_Macro")) {
  temp.list[[n]] <- pt.ov.results[[n]]$ov[["up CP:WIKIPATHWAYS"]] %>% tbl_df() %>% mutate(type="ov",UpDown="UpRe") %>% bind_rows(pt.ov.results[[n]]$MU[["up CP:WIKIPATHWAYS"]] %>% tbl_df() %>% mutate(type="unique",UpDown="UpRe")) %>% bind_rows( pt.ov.results[[n]]$ov[["down CP:WIKIPATHWAYS"]] %>% tbl_df() %>% mutate(type="ov",UpDown="DownRe")) %>% bind_rows( pt.ov.results[[n]]$MU[["down CP:WIKIPATHWAYS"]] %>% tbl_df() %>% mutate(type="unique",UpDown="DownRe"))    
  temp.list[[n]] <- temp.list[[n]] %>% mutate(cluster=n)
}

temp.out <- temp.list %>% do.call("bind_rows",.) %>% mutate(ID=gsub("WP_","",ID))%>% mutate(ID=gsub("_"," ",ID)) %>% mutate(ID=stringr::str_to_title(tolower(ID))) %>% filter(!grepl("Cancer",ID) ) %>% filter(!grepl("Glycolysis",ID))%>% filter(!grepl("Virus",ID))%>% filter(!grepl("Covid19",ID)) %>% filter(!grepl("Disease",ID)) %>% filter(!grepl("Hepatitis",ID)) %>% filter(!grepl("Carcinoma",ID)) %>% filter(!grepl("Vitamin",ID))%>% filter(!grepl("arscov2",ID))%>% filter(!grepl("Prostaglandin",ID)) %>% filter(!grepl("Tumor",ID))%>% filter(!grepl("P53",ID)) %>% filter(!ID %in% c("16p112 Proximal Deletion Syndrome","Acute Viral Myocarditis","Amyotrophic Lateral Sclerosis Als","Apoptosisrelated Network Due To Altered Notch3 In Ovarian Cancer","Cytoplasmic Ribosomal Proteins","Ebola Virus Infection In Host","Electron Transport Chain Oxphos System In Mitochondria","Epo Receptor Signaling","Glutathione Metabolism","Hepatitis B Infection","Hostpathogen Interaction Of Human Coronaviruses Interferon Induction","Immune Response To Tuberculosis","Leptin Signaling Pathway","Mrna Processing","Nonalcoholic Fatty Liver Disease","Pathogenic Escherichia Coli Infection","Pathways Affected In Adenoid Cystic Carcinoma","Photodynamic Therapyinduced Unfolded Protein Response","Prolactin Signaling Pathway","Thyroid Stimulating Hormone Tsh Signaling Pathway","Sarscov2 Innate Immunity Evasion And Cellspecific Immune Response","Translation Factors","Cellular Proteostasis","Oxidative Phosphorylation","Tca Cycle Aka Krebs Or Citric Acid Cycle","Effect Of Progerin On Genes Involved In Hutchinsongilford Progeria Syndrome","Vasopressinregulated Water Reabsorption","Mecp2 And Associated Rett Syndrome","Oncostatin M Signaling Pathway","Mitochondrial Complex I Assembly Model Oxphos System","Copper Homeostasis","Tyrobp Causal Network In Microglia","Microglia Pathogen Phagocytosis Pathway","Tcell Antigen Receptor Tcr Pathway During Staphylococcus Aureus Infection","Modulators Of Tcr Signaling And T Cell Activation","Serotonin Receptor 467 And Nr3c Signaling","Tp53 Network","Lung Fibrosis","Photodynamic Therapyinduced Hif1 Survival Signaling","Cori Cycle","Network Map Of Sarscov2 Signaling Pathway","Nuclear Receptors Metapathway","Prostaglandin Signaling","Photodynamic Therapyinduced Nfkb Survival Signaling","Spinal Cord Injury","Congenital Generalized Lipodystrophy Cgl","Il1 And Megakaryocytes In Obesity","Apoptosis","Selenium Micronutrient Network","Neuroinflammation And Glutamatergic Signaling","Ferroptosis","Tolllike Receptor Signaling Pathway","Antiviral And Antiinflammatory Effects Of Nrf2 On Sarscov2 Pathway","Sudden Infant Death Syndrome Sids Susceptibility Pathways","Thymic Stromal Lymphopoietin Tslp Signaling Pathway","Focal Adhesion Pi3kaktmtorsignaling Pathway","Focal Adhesion","Overview Of Nanoparticle Effects","Pancreatic Adenocarcinoma Pathway","Gastrin Signaling Pathway","Unfolded Protein Response","Mammary Gland Development Pathway Puberty Stage 2 Of 4","Unfolded Protein Response","Adipogenesis","Tgfbeta Receptor Signaling In Skeletal Dysplasias","Nrf2are Regulation","Canonical Nfkb Pathway","Transcriptional Activation By Nrf2 In Response To Phytochemicals","Mirnas Involvement In The Immune Response In Sepsis","Orexin Receptor Pathway","Influence Of Laminopathies On Wnt Signaling","Srebf And Mir33 In Cholesterol And Lipid Homeostasis","Agerage Pathway","Orexin Receptor Pathway","Cholesterol Metabolism With Bloch And Kandutschrussell Pathways","Lncrnamediated Mechanisms Of Therapeutic Resistance","Ltf Danger Signal Response Pathway")) %>% mutate(ID=stringr::str_to_sentence(ID))

temp.sel.GO.ID <- list()
temp.sel.GO.ID$UpRe <- temp.out %>% filter(UpDown=="UpRe") %>% group_by(cluster,UpDown,type) %>% top_n(5,-1*pvalue) %>% tibble::rowid_to_column("nid") %>% ungroup() %>%select(ID) %>% unique()
temp.sel.GO.ID$DownRe <-  temp.out %>% filter(UpDown=="DownRe") %>% group_by(cluster,UpDown,type) %>% top_n(5,-1*pvalue) %>% tibble::rowid_to_column("nid") %>% ungroup() %>%select(ID) %>% unique()

for (m in c("UpRe","DownRe")) {
  temp <- temp.out %>% filter(UpDown==m)  %>% filter(ID %in% temp.sel.GO.ID[[m]]$ID)   %>% mutate(Nlog10Pvalue=-log10(pvalue)) %>% mutate(cluster=paste(cluster,type,sep=":")) %>% select(ID,cluster,Nlog10Pvalue) %>% spread(ID,Nlog10Pvalue) %>% replace(.,is.na(.),0) %>% tibble::column_to_rownames("cluster")
  temp.ph <-temp %>% pheatmap(scale="column",cluster_rows=F,cluster_cols=T)
  temp.sel.ID.od <- colnames(temp)[temp.ph$tree_col$order]
  temp <- temp.out %>% filter(UpDown==m)  %>% filter(ID  %in% temp.sel.ID.od) %>% mutate(ID=factor(ID,temp.sel.ID.od,ordered = T)) %>% arrange(ID)  %>% mutate(Nlog10Pvalue=-log10(pvalue)) 
  
  temp$Nlog10Pvalue[temp$Nlog10Pvalue >10]=10
  temp <- temp %>% filter(pvalue < 0.05)
  plot.results[[paste0("macro",m,"MUO.pt")]] <- temp %>% ggplot()+ geom_point(mapping=aes(y=ID,x=paste(cluster,type,sep=":"),size=Count,col=Nlog10Pvalue))+xlab("")+ylab("WIKIpathway") +theme_classic() +scale_colour_gradient(low="lightgrey",high=UpDown.col[[m]])+theme_bw()+scale_size_continuous(breaks =c(3,10,20,35,50) ,limits = c(2,280),labels = c(3,10,20,35,">50"),range = c(1, 8) )+ theme(axis.text.x=element_text(angle = 90))+ggtitle(m)+theme(plot.title = element_text(hjust=0.5))
}
