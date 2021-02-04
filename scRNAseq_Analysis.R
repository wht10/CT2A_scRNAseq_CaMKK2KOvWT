BiocManager::install("Seurat")
library(tidyverse)
library(Seurat)
library(patchwork)

KO <- Read10X("KOraw_feature_bc_matrix")
WT <- Read10X("WTraw_feature_bc_matrix")

identical(colnames(KO$`Gene Expression`),colnames(KO$`Antibody Capture`))

KOs <- CreateSeuratObject(KO$`Gene Expression`, project = "CaMKK2KO", min.cells = 3, min.features = 200)
WTs <- CreateSeuratObject(WT$`Gene Expression`, project = "WT", min.cells = 3, min.features = 200)

KOHTOs <- KO$`Antibody Capture`[,Cells(KOs)]
WTHTOs <- WT$`Antibody Capture`[,Cells(WTs)]

KOs[["HTO"]] <- CreateAssayObject(KOHTOs)
WTs[["HTO"]] <- CreateAssayObject(WTHTOs)

KOs <- NormalizeData(KOs, assay = "HTO", normalization.method = "CLR")
WTs <- NormalizeData(WTs, assay = "HTO", normalization.method = "CLR")

KOs <- HTODemux(KOs, assay = "HTO", positive.quantile = .99)
WTs <- HTODemux(WTs, assay = "HTO", positive.quantile = .99)

table(KOs$HTO_classification.global)
table(WTs$HTO_classification.global)

Idents(KOs) <- "HTO_maxID"
RidgePlot(KOs, assay = "HTO", features = rownames(KOs[["HTO"]])[1:4],ncol = 4)

Idents(WTs) <- "HTO_maxID"
RidgePlot(WTs, assay = "HTO", features = rownames(WTs[["HTO"]])[1:4],ncol = 4)

FeatureScatter(KOs, feature1 = "hto_M1", feature2 = "hto_M2")

table(KOs@meta.data$HTO_maxID)
table(WTs@meta.data$HTO_maxID)

Idents(KOs) <- "HTO_classification.global"
KOss <- subset(KOs, idents = "Singlet")

Idents(WTs) <- "HTO_classification.global"
WTss <- subset(WTs, idents = "Singlet")

table(KOss@meta.data$HTO_maxID)
table(WTss@meta.data$HTO_maxID)

KOss[["Genotype"]] <- "Camkk2KO"
WTss[["Genotype"]] <- "WT"

#####

m.ct2a <- merge(KOss, y = WTss, add.cell.ids = c("Camkk2KO","WT"), project = "merged.CT2a")
mID <- paste(m.ct2a$Genotype,m.ct2a$HTO_maxID,sep = "_")
m.ct2a[["MouseID"]] <- mID
table(m.ct2a$MouseID)

library(readxl)
MitoCarta <- read_xls("Mouse.MitoCarta3.0.xls", sheet = 2, range = "C1:C894",col_names = T)
i <- MitoCarta$Symbol %in% rownames(m.ct2a)
MitoSymbols <- MitoCarta$Symbol[i]
m.ct2a[["percent.mt"]] <- PercentageFeatureSet(m.ct2a, features = MitoSymbols)
VlnPlot(m.ct2a, features =  c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(m.ct2a, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(m.ct2a, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

library(scater)
m.ct2a.sce <- as.SingleCellExperiment(m.ct2a)

df <- data.frame(NExprs = m.ct2a.sce$nFeature_RNA, LibSize = m.ct2a.sce$nCount_RNA, MitoProp = m.ct2a.sce$percent.mt)
df <- (as.data.frame(df))
df[,1] <- as.numeric(df$NExprs)

qc.lib1 <- isOutlier(log(df$LibSize), type = "lower")
attr(qc.lib1, "threshold")

qc.lib2 <- isOutlier(log(df$LibSize), type = "higher")
attr(qc.lib2, "threshold")

qc.nexp1 <- isOutlier(log(df$NExprs), type = "lower")
attr(qc.nexp1, "threshold")

qc.nexp2 <- isOutlier(log(df$NExprs), type = "higher")
attr(qc.nexp2, "threshold")

qc.mito <- isOutlier(df$MitoProp, type = "higher")
attr(qc.mito, "threshold")


discard <- qc.lib1 | qc.lib2 | qc.nexp1 | qc.nexp2 | qc.mito

i <- Cells(m.ct2a)[!discard]
m.ct2a <- subset(m.ct2a, cells = i)

library(biomaRt)
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

s.genes <- convertHumanGeneList(cc.genes.updated.2019$s.genes) 
g2m.genes <- convertHumanGeneList(cc.genes.updated.2019$g2m.genes)
m.ct2a <- CellCycleScoring(m.ct2a, s.features = s.genes, g2m.features = g2m.genes)

length(m.ct2a@assays$RNA@counts@Dimnames[[1]])

m.ct2a <- SCTransform(m.ct2a, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))

m.ct2a <- RunPCA(m.ct2a)
ElbowPlot(m.ct2a, ndims = 50)
DimHeatmap(m.ct2a, dims = 1:30, cells = 500)

m.ct2a <- FindNeighbors(m.ct2a,dims = 1:30)

m.ct2a <- RunUMAP(m.ct2a, dims = 1:30)

m.ct2a <- FindClusters(m.ct2a, resolution = .5)

m.ct2a <- NormalizeData(m.ct2a, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
all.genes <- rownames(m.ct2a)
m.ct2a <- ScaleData(m.ct2a, features = all.genes, assay = "RNA")

DimPlot(m.ct2a, label = T, repel = T,label.size = 8) + theme(legend.position = "none")
ggsave("m.ct2a_dimplot.pdf",device = "pdf", width = 2, height = 2, units = "in", scale = 6)
DimPlot(m.ct2a, split.by = "Genotype") + theme(legend.position = "none")
DotPlot(m.ct2a, assay = "RNA", features = c("Tcf7","Slamf6","Cd4","Cd8a","Cd3e",
                                            "Ncr1","Klrb1c",
                                            "Cd19",
                                            "Siglech","Tcf4","Fscn1","Ccr7","Itgax","H2-DMb2","Xcr1","Clec9a",
                                            "Flt3","Zbtb46",
                                            "Gata2","Ms4a2",
                                            "S100a9","Mmp9",
                                            "Ace","Adgre4",
                                            "Vcan","Ly6c2",
                                            "P2ry12","Sall1","Crybb1",
                                            "Nos2","Arg1",
                                            "H2-Aa","Cd74",
                                            "Apoe","Mrc1","Mertk","Adgre1")) + RotatedAxis() + theme(axis.text = element_text(size = 20), axis.title = element_blank()) 
ggsave("m.ct2a_featuredotplot.pdf", width = 3, height = 1, units = "in", scale = 6)

DimPlot(m.ct2a, repel = T, group.by = "MouseID")
DimPlot(m.ct2a, repel = T, split.by = "MouseID")
table(m.ct2a$MouseID)


#Cell Annotation
#####

DefaultAssay(m.ct2a) <- "RNA"

markers <- lapply(c(0:(length(levels(Idents(m.ct2a)))-1)), function(i){
  FindConservedMarkers(m.ct2a, ident.1 = i, grouping.var = "Genotype", verbose = T)
})

zerov1 <- FindMarkers(m.ct2a, ident.1 = 0, ident.2 = 1, grouping.var = "Genotype", verbose = T)
onevall <- FindMarkers(m.ct2a, ident.1 = 1, verbose = T)

markers[[1]][1:75,1:10]



#NK
FeaturePlot(m.ct2a,features = "Ncr1")
FeaturePlot(m.ct2a,features = "Klrb1c")
FeaturePlot(m.ct2a,features = "Gzma")
FeaturePlot(m.ct2a,features = "Txk")
FeaturePlot(m.ct2a,features = "Prf1")


#TAM-1
FeaturePlot(m.ct2a,features = "Apoe")
FeaturePlot(m.ct2a,features = "Trem2")
FeaturePlot(m.ct2a,features = "Mrc1")
FeaturePlot(m.ct2a,features = "Spint1")
FeaturePlot(m.ct2a,features = "Spp1")
FeaturePlot(m.ct2a,features = "Itga6")
FeaturePlot(m.ct2a,features = "Cd63")
FeaturePlot(m.ct2a,features = "Mertk")
FeaturePlot(m.ct2a,features = "Adgre1")

#TAM-2
FeaturePlot(m.ct2a,features = "H2-Eb1")
FeaturePlot(m.ct2a,features = "H2-Ab1")
FeaturePlot(m.ct2a,features = "H2-Aa")
FeaturePlot(m.ct2a,features = "Stat1")
FeaturePlot(m.ct2a,features = "Irf7")
FeaturePlot(m.ct2a,features = "Cd74")
FeaturePlot(m.ct2a,features = "Ccl5")
FeaturePlot(m.ct2a,features = "Cxcl9")
FeaturePlot(m.ct2a,features = "Mertk")
FeaturePlot(m.ct2a,features = "Adgre1")
FeaturePlot(m.ct2a,features = "Cxcl10")

#cMono
FeaturePlot(m.ct2a,features = "Vcan")
FeaturePlot(m.ct2a,features = "Ly6c2")
FeaturePlot(m.ct2a,features = "Il1b")
FeaturePlot(m.ct2a,features = "Itgal")
FeaturePlot(m.ct2a,features = "Ccl9")
FeaturePlot(m.ct2a,features = "Lyz2")
FeaturePlot(m.ct2a,features = "Chil3")
FeaturePlot(m.ct2a,features = "Cd44")

#CD8 T Cells
FeaturePlot(m.ct2a,features = "Cd8a")
FeaturePlot(m.ct2a,features = "Cd3e")
FeaturePlot(m.ct2a,features = "Cxcr6")
FeaturePlot(m.ct2a,features = "Ifng")
FeaturePlot(m.ct2a,features = "Tox")
FeaturePlot(m.ct2a,features = "Pdcd1")
FeaturePlot(m.ct2a,features = "Ctla4")
FeaturePlot(m.ct2a,features = "Lag3")
FeaturePlot(m.ct2a,features = "Klrc1")
FeaturePlot(m.ct2a,features = "Havcr2")#Tim3

#microglia
FeaturePlot(m.ct2a,features = "P2ry12")
FeaturePlot(m.ct2a,features = "Sall1")
FeaturePlot(m.ct2a,features = "Crybb1")
FeaturePlot(m.ct2a,features = "Sparc")
FeaturePlot(m.ct2a,features = "Csf1r")
FeaturePlot(m.ct2a,features = "Aif1")
FeaturePlot(m.ct2a,features = "Tmem119")
FeaturePlot(m.ct2a,features = "Cd34")

#cDC2
FeaturePlot(m.ct2a,features = "Flt3")
FeaturePlot(m.ct2a,features = "Itgax")
FeaturePlot(m.ct2a,features = "Ccl17")
FeaturePlot(m.ct2a,features = "Ffar2")
FeaturePlot(m.ct2a,features = "H2-DMb2")
FeaturePlot(m.ct2a,features = "Zbtb46")

#Stem-Like TILs
FeaturePlot(m.ct2a,features = "Cd4")
FeaturePlot(m.ct2a,features = "Cd8a")
FeaturePlot(m.ct2a,features = "Trgc4")
FeaturePlot(m.ct2a,features = "Rora")
FeaturePlot(m.ct2a,features = "Icos")
FeaturePlot(m.ct2a,features = "Tcf7")
FeaturePlot(m.ct2a,features = "Il7r")
FeaturePlot(m.ct2a,features = "Slamf6")

#CD4 TILs
FeaturePlot(m.ct2a,features = "Cd4")
FeaturePlot(m.ct2a,features = "Cxcr3")
FeaturePlot(m.ct2a,features = "Foxp3")
FeaturePlot(m.ct2a,features = "Il10")
FeaturePlot(m.ct2a,features = "Klrg1")
FeaturePlot(m.ct2a,features = "Ctla4")
FeaturePlot(m.ct2a,features = "Icos")
FeaturePlot(m.ct2a,features = "Ccr8")

#B Cells
FeaturePlot(m.ct2a,features = "Cd19")
FeaturePlot(m.ct2a,features = "Cd79a")
FeaturePlot(m.ct2a,features = "Fcmr")
FeaturePlot(m.ct2a,features = "Pax5")
FeaturePlot(m.ct2a,features = "Ebf1")
FeaturePlot(m.ct2a,features = "Iglc2")

#TAM3
FeaturePlot(m.ct2a,features = "Nos2")
FeaturePlot(m.ct2a,features = "Arg1")
FeaturePlot(m.ct2a,features = "Vegfa")
FeaturePlot(m.ct2a,features = "Mmp12")
FeaturePlot(m.ct2a,features = "Egln3")
FeaturePlot(m.ct2a,features = "Mertk")
FeaturePlot(m.ct2a,features = "Adgre1")
FeaturePlot(m.ct2a,features = "Hilpda")
FeaturePlot(m.ct2a,features = "Id2")

#non cMono
FeaturePlot(m.ct2a,features = "Ace")
FeaturePlot(m.ct2a,features = "Adgre4")
FeaturePlot(m.ct2a,features = "Ear2")
FeaturePlot(m.ct2a,features = "Fabp4")
FeaturePlot(m.ct2a,features = "Pparg")
FeaturePlot(m.ct2a,features = "Dnah12")
FeaturePlot(m.ct2a,features = "Cd36")
FeaturePlot(m.ct2a,features = "Ceacam1")
FeaturePlot(m.ct2a,features = "Spn") #Cd43

#DC-3
FeaturePlot(m.ct2a,features = "Batf3")
FeaturePlot(m.ct2a,features = "Ccl22")
FeaturePlot(m.ct2a,features = "Ramp3")
FeaturePlot(m.ct2a,features = "Zbtb46")
FeaturePlot(m.ct2a,features = "Il12b")
FeaturePlot(m.ct2a,features = "Ccr7")
FeaturePlot(m.ct2a,features = "Fscn1")
FeaturePlot(m.ct2a,features = "Sema7a")
FeaturePlot(m.ct2a,features = "Flt3")
FeaturePlot(m.ct2a,features = "Cd70")
FeaturePlot(m.ct2a,features = "Stat4")

#DC1
FeaturePlot(m.ct2a,features = "Xcr1")
FeaturePlot(m.ct2a,features = "Snx22")
FeaturePlot(m.ct2a,features = "Clec9a")
FeaturePlot(m.ct2a,features = "Cd24a")
FeaturePlot(m.ct2a,features = "Tlr3")
FeaturePlot(m.ct2a,features = "Batf3")
FeaturePlot(m.ct2a,features = "Itgae")

#pDC
FeaturePlot(m.ct2a,features = "Klk1")
FeaturePlot(m.ct2a,features = "Ccr9")
FeaturePlot(m.ct2a,features = "Tcf4")
FeaturePlot(m.ct2a,features = "Flt3")
FeaturePlot(m.ct2a,features = "Siglech")
FeaturePlot(m.ct2a,features = "Bcl11a")

#Basophils
FeaturePlot(m.ct2a,features = "Il4")
FeaturePlot(m.ct2a,features = "Gata2")
FeaturePlot(m.ct2a,features = "Ms4a2")
FeaturePlot(m.ct2a,features = "Fcer1a")
FeaturePlot(m.ct2a,features = "Il13")
FeaturePlot(m.ct2a,features = "Slc18a2")
FeaturePlot(m.ct2a,features = "Il6")
FeaturePlot(m.ct2a,features = "Cxcr2")

#Neuts
FeaturePlot(m.ct2a,features = "S100a9")
FeaturePlot(m.ct2a,features = "Mmp9")
FeaturePlot(m.ct2a,features = "Ly6g")
FeaturePlot(m.ct2a,features = "Retnlg")
FeaturePlot(m.ct2a,features = "Stfa2l1")


new.cluster.ids <- c("Apoe+ TAM", "DC-like TAM","NK Cell","cMono","CD8 TIL","Microglia","DC2","Stem-Like TIL","CD4 TIL","B Cell","Nos2+ TAM","non-cMono","DC3","Debris","DC1","pDC","Baso","Neut")
m.ct2a[["orig_clust"]] <- Idents(m.ct2a)
names(new.cluster.ids) <- levels(m.ct2a)
m.ct2a <- RenameIdents(m.ct2a,new.cluster.ids)

debri <- which(Idents(m.ct2a) == "Debris")
i <- Cells(m.ct2a)[debri]
m.ct2a <- m.ct2a[,!colnames(m.ct2a) %in% i]

levels(m.ct2a) <- c("Apoe+ TAM","DC-like TAM","Nos2+ TAM","Microglia",
                    "cMono","non-cMono","Neut","Baso",
                    "DC1","DC2","DC3","pDC",
                    "B Cell","NK Cell","CD8 TIL","CD4 TIL","Stem-Like TIL")

m.ct2a[["Cell.Type"]] <- Idents(m.ct2a)
m.ct2a[["Geno.Ident"]] <- paste(m.ct2a$Genotype,m.ct2a$Cell.Type,sep = "_")

DimPlot(m.ct2a,label = T, repel = T)

#Differential Expression Genotypes
##### 

diffgenolist <- c("Microglia",
                  "cMono","non-cMono","Neut","Baso",
                  "DC1","DC2","DC3","pDC",
                  "B Cell","NK Cell")

Idents(m.ct2a) <- m.ct2a$Geno.Ident
lev <- levels(Idents(m.ct2a))
lev <- lev[order(lev)]
lev[c(2,3,6,8,9,10,11,12,13,14,16)]

g.markers <- lapply(c(2,3,6,8,9,10,11,12,13,14,16), function(i){
  FindMarkers(m.ct2a, ident.1 = lev[i], ident.2 =  lev[i+17], min.pct = .1, test.use = "MAST", verbose = T, assay = "RNA")
})

lev[2]
#B cell
g.markers[[1]][1:100,]
EnhancedVolcano(g.markers[[1]],lab=rownames(g.markers[[1]]),x="avg_logFC",y="p_val_adj", 
                selectLab = (rownames(g.markers[[1]] %>% slice_max(n = 40, order_by = abs(avg_logFC)))), 
                pCutoff = NA,FCcutoff = NA,drawConnectors = T, title = "KO B cell vs. WT B Cell", 
                subtitle = NULL, legendPosition = "none", col = c("black","black","black","red")) 

lev[3]
#Baso - no DEGs
g.markers[[2]][1:100,]
EnhancedVolcano(g.markers[[2]],lab=rownames(g.markers[[2]]),x="avg_logFC",y="p_val_adj", 
                selectLab = (rownames(g.markers[[2]] %>% slice_max(n = 40, order_by = abs(avg_logFC)))), 
                pCutoff = NA,FCcutoff = NA,drawConnectors = T, title = "KO Basophil vs. WT Basophil", 
                subtitle = NULL, legendPosition = "none", col = c("black","black","black","red")) 

lev[6]
#cMono
g.markers[[3]][1:100,]
EnhancedVolcano(g.markers[[3]],lab=rownames(g.markers[[3]]),x="avg_logFC",y="p_val_adj", 
                selectLab = (rownames(g.markers[[3]] %>% slice_max(n = 40, order_by = abs(avg_logFC)))), 
                pCutoff = NA,FCcutoff = NA,drawConnectors = T, title = "KO cMono vs. WT cMono", 
                subtitle = NULL, legendPosition = "none", col = c("black","black","black","red")) 

lev[8]
#DC1
g.markers[[4]][1:100,]
EnhancedVolcano(g.markers[[4]],lab=rownames(g.markers[[4]]),x="avg_logFC",y="p_val_adj", 
                selectLab = (rownames(g.markers[[4]] %>% slice_min(n = 25, order_by = p_val_adj))), 
                pCutoff = NA,FCcutoff = NA,drawConnectors = T, title = "KO DC1 vs. WT DC1", 
                subtitle = NULL, legendPosition = "none", col = c("black","black","black","red"))

lev[9]
#DC2
g.markers[[5]][1:100,]
EnhancedVolcano(g.markers[[5]],lab=rownames(g.markers[[5]]),x="avg_logFC",y="p_val_adj", 
                selectLab = (rownames(g.markers[[5]] %>% slice_max(n = 45, order_by = abs(avg_logFC)))), 
                pCutoff = NA,FCcutoff = NA,drawConnectors = T, title = "KO DC2 vs. WT DC2", 
                subtitle = NULL, legendPosition = "none", col = c("black","black","black","red"))

lev[10]
#DC3
g.markers[[6]][1:100,]
EnhancedVolcano(g.markers[[6]],lab=rownames(g.markers[[6]]),x="avg_logFC",y="p_val_adj", 
                selectLab = (rownames(g.markers[[6]] %>% slice_min(n = 25, order_by = p_val_adj))), 
                pCutoff = NA,FCcutoff = NA,drawConnectors = T, title = "KO DC3 vs. WT DC3", 
                subtitle = NULL, legendPosition = "none", col = c("black","black","black","red"))

lev[11]
#Microglia
g.markers[[7]][1:100,]
EnhancedVolcano(g.markers[[7]],lab=rownames(g.markers[[7]]),x="avg_logFC",y="p_val_adj", 
                selectLab = (c(rownames(g.markers[[7]] %>% slice_min(n = 50, order_by = p_val_adj)),"Spp1")), 
                pCutoff = NA,FCcutoff = NA,drawConnectors = T, title = "KO Microglia vs. WT Microglia", 
                subtitle = NULL, legendPosition = "none", col = c("black","black","black","red"))

lev[12]
#Neutrophil - NO DEGs
g.markers[[8]][1:100,]
EnhancedVolcano(g.markers[[8]],lab=rownames(g.markers[[8]]),x="avg_logFC",y="p_val_adj", 
                selectLab = (c(rownames(g.markers[[8]] %>% slice_min(n = 50, order_by = p_val_adj)))), 
                pCutoff = NA,FCcutoff = NA,drawConnectors = T, title = "KO Neuts vs. WT Neuts", 
                subtitle = NULL, legendPosition = "none", col = c("black","black","black","red"))

lev[13]
#NK
g.markers[[9]][1:100,]
EnhancedVolcano(g.markers[[9]],lab=rownames(g.markers[[9]]),x="avg_logFC",y="p_val_adj", 
                selectLab = (c(rownames(g.markers[[9]] %>% slice_min(n = 50, order_by = p_val_adj)))), 
                pCutoff = NA,FCcutoff = NA,drawConnectors = T, title = "KO NK vs. WT NK", 
                subtitle = NULL, legendPosition = "none", col = c("black","black","black","red"))


lev[14]
#non-cMono
g.markers[[10]][1:100,]
EnhancedVolcano(g.markers[[10]],lab=rownames(g.markers[[10]]),x="avg_logFC",y="p_val_adj", 
                selectLab = (c(rownames(g.markers[[10]] %>% slice_min(n = 50, order_by = p_val_adj)))), 
                pCutoff = NA,FCcutoff = NA,drawConnectors = T, title = "KO non-cMono vs. WT non-cMono", 
                subtitle = NULL, legendPosition = "none", col = c("black","black","black","red"))

lev[16]
#pDC
g.markers[[11]][1:100,]
EnhancedVolcano(g.markers[[11]],lab=rownames(g.markers[[11]]),x="avg_logFC",y="p_val_adj", 
                selectLab = (c(rownames(g.markers[[11]] %>% slice_min(n = 18, order_by = p_val_adj)))), 
                pCutoff = NA,FCcutoff = NA,drawConnectors = T, title = "KO pDC vs. WT pDC", 
                subtitle = NULL, legendPosition = "none", col = c("black","black","black","red"))

#All Cell Type Heatmap
#####

Idents(m.ct2a) <- m.ct2a$Geno.Ident
ct2a.avg <- AverageExpression(m.ct2a,assays = "RNA",add.ident = "hash.ID", return.seurat = T)

DoHeatmap(subset(ct2a.avg, idents = c("Microglia","cMono","non-cMono")), 
          features = c("Apoe","Spint1","Cd74","H2-Aa","Ifi44","Ccl5","Cxcl9","Irf7","Spp1","Stat1","B2m","Chil3"),size = 3,draw.lines = F) + 
  scale_fill_gradientn(colors = c("blue","white","red")) + scale_x_discrete(labels = rep(c("CaMKK2KO","WT")))


#CaMKK2 Counts
##### 

counts <- GetAssayData(m.ct2a, slot = "data", assay = "RNA")

camkk2.counts <- counts["Camkk2",]

camkk2 <- data.frame(counts = camkk2.counts, Genotype = m.ct2a$Genotype, Celltype = Idents(m.ct2a))

VlnPlot(subset(m.ct2a, subset = Genotype == "WT"),features = "Camkk2", assay = "RNA", sort = "increasing", pt.size = 0) + theme(legend.position = "none")


library(kableExtra)
camkk2 %>% group_by(Genotype) %>% summarise(avg = mean(counts),
                                            percent_over_zero = mean(counts != 0)*100) %>% 
  arrange(desc(avg)) %>% 
  kable(col.names = c("Celltype", "Average Expression",
                      "% of Cells with non-Zero Expression"),
        caption = "CaMKK2 Expression - Tumor Mouse") %>% kable_styling() 

#plotting cell frequency
#####

meta <- m.ct2a@meta.data
idents <- m.ct2a@active.ident
meta <- meta[,c("Genotype","MouseID","HTO_maxID")]
meta <- meta %>% add_column(Idents = idents) 
meta[,(c("Genotype","MouseID","HTO_maxID"))] <- lapply(meta[,(c("Genotype","MouseID","HTO_maxID"))], as.factor)

library(stringr)
meta1 <- meta %>% group_by(MouseID,Idents) %>% summarise(n = n()) %>% mutate(freq =  prop.table(n)*100)
geno <- c(rep("Camkk2KO",67),rep("WT",66))
ungroup(meta1) %>% mutate(Geno = geno) %>% ggplot(aes(Idents,freq,fill=Geno)) + geom_boxplot(position = position_dodge())
Cell_freq <- ungroup(meta1) %>% mutate(Geno = geno)
write_excel_csv(Cell_freq,"Cell_freq.csv")

prop.table(table(Idents(m.ct2a),m.ct2a$MouseID),margin = 2)*100


#TAMDEGs and Venn Diagram
#####
DefaultAssay(m.ct2a) <- "RNA"
tam.ct2a <- subset(m.ct2a,idents = c("DC-like TAM","Apoe+ TAM","Nos2+ TAM"))

DimPlot(subset(m.ct2a, idents = c("cMono","Nos2+ TAM","Apoe+ TAM","DC-like TAM","non-cMono", "Microglia")), label = T, repel = T, label.size = 8) + theme(legend.position = "none")
ggsave("TAMDim.pdf",height = 2,width = 2,scale = 3)

tamdegmast <- FindAllMarkers(subset(m.ct2a, idents = c("cMono","Nos2+ TAM","Apoe+ TAM","DC-like TAM","non-cMono", "Microglia")), 
                                    only.pos = T,min.pct = .1,logfc.threshold = .1, test.use = "MAST",verbose = T, assay = "RNA")
tamdegsmast <- tamdegmast %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
#DoHeatmap(tam.ct2a, cells = WhichCells(tam.ct2a,downsample = 150),features = c(tamdegsmast$gene,"Nos2"))
library(viridis)
levels(ct2a.avg) <- c("WT_cMono","Camkk2KO_cMono","WT_Nos2+ TAM","Camkk2KO_Nos2+ TAM","WT_Apoe+ TAM",
                      "Camkk2KO_Apoe+ TAM","WT_DC-like TAM","Camkk2KO_DC-like TAM","WT_non-cMono","Camkk2KO_non-cMono",
                      "WT_Microglia","Camkk2KO_Microglia")
  
DoHeatmap(subset(ct2a.avg, idents = c("WT_cMono","Camkk2KO_cMono","WT_Nos2+ TAM","Camkk2KO_Nos2+ TAM","WT_Apoe+ TAM",
                                      "Camkk2KO_Apoe+ TAM","WT_DC-like TAM","Camkk2KO_DC-like TAM","WT_non-cMono","Camkk2KO_non-cMono",
                                      "WT_Microglia","Camkk2KO_Microglia")), 
          features = c("Apoe","Spint1","Spp1","Cd74","H2-Aa","Ifi44","Cxcl9","Irf7","Spp1","Stat1"),draw.lines = F, size = 3.5,angle = 35, disp.min = 0) + 
  scale_fill_viridis(option = "plasma") 
ggsave("MNPHeatmap.pdf",height = 2,width = 3.5,scale = 3)

#TAM1degs <- FindMarkers(m.ct2a, ident.1 = "DC-like TAM",ident.2 = c("Apoe+ TAM","Nos2+ TAM"), min.pct = .25, test.use = "MAST")
tamdcvapo <- FindMarkers(m.ct2a, ident.1 = "DC-like TAM", ident.2 = "Apoe+ TAM", logfc.threshold = .1, min.pct = .1, test.use = "MAST", verbose = T, assay = "RNA")


library(EnhancedVolcano)

#CD274 is PDL1

EnhancedVolcano(tamdcvapo,lab=rownames(tamdcvapo),x="avg_log2FC",y="p_val_adj", 
                selectLab = c("H2-Aa","Spp1","Stat1","Irf7","Cd74","Ccl5","Cxcl9","Cxcl10","Ccl6","Ccl9", "Vcam1", "Cd83",
                              "C3","Apoe","Cd63","Itga6","Ly6c2","Trem2","Mrc1","Arg1", "C1qa",
                              "Cd9","Ciita","B2m","Tap1","Stat1","Ifi44","Ifi202b","Ifitm3",
                              "Fcrls","Olfml3","Sell","Plac8","Arg2","Irf1","Hif1a","Nos2","Il18bp","Cd47","Cx3cr1","Ccr2","Cd40",
                              "Nfkbia","Egr1","Lpl","Cd63","Ctsl","Cd68","Cadm1","Timp2","Card9"), 
                pCutoff = NA,FCcutoff = NA,drawConnectors = T, title = "DC-like TAM vs. Apoe+ TAM", 
                subtitle = NULL, legendPosition = "none", col = c("red","red","red","red"),colAlpha = .2,ylab = bquote(~-Log[10]~"(adj.p-val)"), 
                xlab =bquote(~Log[2]~"FC(DC-like TAM/Apoe+ TAM)") ,caption = "", labSize = 7) 
ggsave("KOvWTMacVolc.pdf",width = 7,height = 4,scale = 2)

VlnPlot(tam.ct2a,features = c("Trem2","Spp1","Apoe","Cxcl9","Ccl5","H2-Aa","Mrc1","Cx3cr1","Stat1","Irf7","C1qa"), sort = "increasing", pt.size = 0)

#Trem2 Module score
library(Nebulosa)

#genes pulled from FigS5 Upregulated Stage 1 + 2 DAM https://www.cell.com/fulltext/S0092-8674(17)30578-0
damgene2 <- list(c("Ank","Spp1","Axl","Csf1","Cst7","Cd9","Cadm1","Ccl6","Itgax","Cd63",
             "Cd68","Ctsa","Lpl","Serpine2","Cd52","Ctsl","Hif1a"))

damgene <- list(c("Ank","Spp1","Axl","Csf1","Cst7","Cd9","Cadm1","Ccl6","Itgax","Cd63",
              "Cd68","Ctsa","Lpl","Serpine2","Cd52","Ctsl","Hif1a",
             "Apoe","B2m","Cstb","Tyrobp","Timp2","H2-D1","Fth1","Lyz2","Ctsb","Ctsd","Trem2"))

#genes pulled from FigS5 Upregulated DAM4E https://www.cell.com/fulltext/S0092-8674(17)30578-0
consdam <- list(c("Trem2","Axl","Tyrobp","Lpl","Cst7","Ctsb","Apoe","Cd63","Cd9","Fth1","Spp1"))

m.ct2a <- AddModuleScore(m.ct2a,damgene2,name = "DAM_Genes_Adv", assay = "RNA")
m.ct2a <- AddModuleScore(m.ct2a,damgene,name = "DAM_Genes_All", assay = "RNA")
m.ct2a <- AddModuleScore(m.ct2a,consdam,name = "DAM_Genes_Cons", assay = "RNA")

plot_density(subset(m.ct2a, idents = c("cMono","Nos2+ TAM","Apoe+ TAM","DC-like TAM","non-cMono", "Microglia")),
            features = "DAM_Genes_Adv1", pal = "inferno",method = "wkde") + 
  labs(title = "Disease Associated Microglia Signature")

plot_density(subset(m.ct2a, idents = c("cMono","Nos2+ TAM","Apoe+ TAM","DC-like TAM","non-cMono", "Microglia")),
             features = "DAM_Genes_All1", pal = "inferno",method = "wkde") + 
  labs(title = "Disease Associated Microglia Signature")

plot_density(subset(m.ct2a, idents = c("cMono","Nos2+ TAM","Apoe+ TAM","DC-like TAM","non-cMono", "Microglia")),
             features = "DAM_Genes_Cons1", pal = "inferno",method = "wkde") + 
  labs(title = "Disease Associated Microglia Signature")
ggsave("DAMFeatPlot.pdf",width = 3,height = 3,scale = 3)

#TILDegs  and subsetting 
#####
tils <- subset(m.ct2a,idents = c("CD8 TIL","CD4 TIL","Stem-Like TIL"))
DimPlot(tils)

FeaturePlot(tils, features = c("Slamf6","Tcf7"), split.by = "Genotype")
til.count <- GetAssayData(tils, slot = "data", assay = "RNA")
slamf6.count <- til.count["Slamf6",]
tcf7.count <- til.count["Tcf7",]
til.counts <- data.frame(tcf7 = tcf7.count, slamf6 = slamf6.count)
til.counts <- til.counts %>%  filter(tcf7 != 0 & slamf6 != 0)
ggplot(til.counts, aes(tcf7, slamf6)) + geom_point() + geom_smooth(method = "lm")
til.lm <- lm(tcf7 ~ slamf6, data = til.counts)
summary(til.lm)
cor(til.counts,method = "spearman")

Geno_Ident <- paste(tils$Genotype,Idents(tils),sep = "_")
geno_tils <- tils
Idents(geno_tils) <- Geno_Ident
cd8deg <- FindMarkers(geno_tils, ident.1 = "Camkk2KO_CD8 TIL", ident.2 = "WT_CD8 TIL", min.pct = .1, test.use = "MAST")
cd4deg <- FindMarkers(geno_tils, ident.1 = "Camkk2KO_CD4 TIL", ident.2 = "WT_CD4 TIL", min.pct = .1, test.use = "MAST")
stemdeg <- FindMarkers(geno_tils, ident.1 = "Camkk2KO_Stem-Like TIL", ident.2 = "WT_Stem-Like TIL", min.pct = .1, test.use = "MAST")

r.tils <- tils

r.tils <- SCTransform(r.tils, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"))

r.tils <- RunPCA(r.tils)
ElbowPlot(r.tils, ndims = 50)
DimHeatmap(r.tils, dims = 1:30, cells = 500)

r.tils <- FindNeighbors(r.tils,dims = 1:30)

r.tils <- RunUMAP(r.tils, dims = 1:30)

r.tils <- FindClusters(r.tils, resolution = .4)

DimPlot(r.tils, label = T, repel = T,label.size = 8) + theme(legend.position = "none")
ggsave("tilsdimplot.pdf",width = 2,height = 2, scale = 4)
DimPlot(r.tils, split.by = "Genotype") + theme(legend.position = "none")
ggsave("tilsdimplot_genosplit.pdf",width = 4,height = 2, scale = 6)

levels(r.tils) <- c("Gamma-Delta TIL","Stem-Like TIL","CD8 TIL","Effector CD4","Treg")
DotPlot(r.tils, assay = "RNA", features = c("Ikzf2","Foxp3","Icos","Ctla4","Il2ra","Tnfrsf1b","Il10",
                                            "Cd4","Cd40lg","Cd200","Cd82","Zap70","Nfkb1","Il18rap","Ifng","Tnf",
                                            "Cd8a","Cxcr6","Irf8","Pdcd1","Runx3","Lgals3",
                                            "Tcf7","Slamf6","S1pr1","Cd7","Lef1",
                                            "Cd163l1","Trdv4","Rorc","Sox13","Igf1r","Aqp3"
                                           )) + RotatedAxis() + theme(axis.text = element_text(size = 20), axis.title = element_blank()) +
  
ggsave("tilsdotplot.pdf",width = 3, height = 1,scale = 6)

DimPlot(r.tils, repel = T, group.by = "MouseID")
DimPlot(r.tils, repel = T, split.by = "MouseID")

DefaultAssay(r.tils) <- "RNA"

t.markers <- lapply(c(0:(length(levels(Idents(r.tils)))-1)), function(i){
  FindConservedMarkers(r.tils, ident.1 = i, grouping.var = "Genotype", verbose = T)
})


#cluster 0 - CD8 TILs
t.markers[[1]][1:100,c((2:5),(7:10))]
FeaturePlot(r.tils,features = c("Cd8a","Cxcr6","Trgc2","Irf8","Pdcd1","Runx3","Lgals3"))

#cluster 1 - Tregs
t.markers[[2]][1:100,c((2:5),(7:10))]
FeaturePlot(r.tils,features = c("Ikzf2","Foxp3","Icos","Ctla4","Il2ra","Tnfrsf1b","Il10"))

#cluster 2 - Stem-Like Tils
as_tibble(t.markers[[3]][1:100,c((2:5),(7:10))],rownames = NA) %>% filter(WT_avg_logFC > 0) %>% as.data.frame()
FeaturePlot(r.tils,features = c("Tcf7","Slamf6","S1pr1","Cd7","Lef1","Cxcr4"))

#cluster 3 - Th1
VlnPlot(r.tils, features = "Cxcr3", split.by = "Genotype", sort = T)
as_tibble(t.markers[[4]][,c((2:5),(7:10))],rownames = NA) %>% filter(WT_avg_logFC > 0) %>% as.data.frame()
FeaturePlot(r.tils,features = c("Cd4","Cd40lg","Cd200","Cd82","Ccl1","Il1r2","Zap70","Nfkb1","Il18rap","Ifng","Tnf","Il18r1"))

#cluster 4 - gamma delta
as_tibble(t.markers[[5]][1:100,c((2:5),(7:10))],rownames = NA) %>% filter(WT_avg_logFC > 0) %>% as.data.frame()
FeaturePlot(r.tils,features = c("Cd163l1","Trdv4","Trgc1","Trdc","Rorc","Il1r1","Sox13","Il23r","Igf1r","Aqp3"))

#cluster 5 - Cycling TIL
as_tibble(t.markers[[6]][1:100,c((2:5),(7:10))],rownames = NA) %>% filter(WT_avg_logFC > 0) %>% as.data.frame()
FeaturePlot(r.tils,features = c("Cd3e","Nabp1","mt-Atp6","mt-Co2"))

#cluster 6 - TAM that ate TILs
as_tibble(t.markers[[7]][1:100,c((2:5),(7:10))],rownames = NA) %>% filter(WT_avg_logFC > 0) %>% as.data.frame()
FeaturePlot(r.tils,features = c("Trem2","Tmem119","Itgam"))


new.cluster.ids <- c("CD8 TIL", "Treg","Stem-Like TIL","Effector CD4","Gamma-Delta TIL","Cycling TIL","TAM TIL")
r.tils[["orig_clust"]] <- Idents(r.tils)
names(new.cluster.ids) <- levels(r.tils)
r.tils <- RenameIdents(r.tils,new.cluster.ids)

excl <- which(Idents(r.tils) == "Cycling TIL" | Idents(r.tils) == "TAM TIL")
i <- Cells(r.tils)[excl]
r.tils <- r.tils[,!colnames(r.tils) %in% i]

DimPlot(r.tils, label = T, repel = T)
DimPlot(r.tils, label = T, repel = T, group.by = "Genotype")
DimPlot(r.tils, repel = T, split.by = "MouseID")


meta.t <- r.tils@meta.data
idents.t <- r.tils@active.ident
meta.t <- meta.t[,c("Genotype","MouseID","HTO_maxID")]
meta.t <- meta.t %>% add_column(Idents = idents.t) 

library(stringr)
meta.t2 <- meta.t %>% group_by(MouseID,Idents) %>% summarise(n = n()) %>% mutate(freq =  prop.table(n)*100)
meta.t2 <- meta.t2 %>% mutate(geno = gsub('_.*','',MouseID))
ungroup(meta.t2) %>% mutate(Geno = geno) %>% ggplot(aes(Idents,freq,fill=Geno)) + geom_boxplot(position = position_dodge())
til_freq <- ungroup(meta.t2)
write_excel_csv(til_freq,"TIL_freq.csv")

DefaultAssay(r.tils)<- "RNA"
tildegmast <- FindAllMarkers(r.tils, only.pos = T,min.pct = .1,logfc.threshold = .1, test.use = "MAST", assay = "RNA")
tildegsmast <- tildegmast %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)

DoHeatmap(r.tils, features = c(tildegsmast$gene))

VlnPlot(r.tils,
        features = c("Pdcd1","Tox"), 
        idents = c("CD8 TIL"), 
        split.by = "Genotype", split.plot = T)

Geno_Ident <- paste(tils$Genotype,Idents(r.tils),sep = "_")
geno_tils <- r.tils

l <- levels(geno_tils)
l[6] <- l[7]
l[7] <- l6
t.g.markers <- lapply(c(1:(length(levels(Idents(geno_tils)))/2)), function(i){
  FindMarkers(geno_tils, ident.1 = l[i], ident.2 =  l[i+5], min.pct = .1, test.use = "MAST", verbose = T, assay = "RNA")
})

#TIL KO V WT
#CD8
t.g.markers[[1]][1:75,]
 VlnPlot(r.tils,
        features = c("Apoe","Tox","Stat3","Gzmb","Gzma","Ccr2"), 
        idents = "CD8 TIL", 
        split.by = "Genotype", split.plot = F, combine = T, pt.size = 0,cols = c("blue","red"))
ggsave("CD8DEG.pdf",height = 3,width = 2,scale = 2)

EnhancedVolcano(t.g.markers[[1]],lab=rownames(t.g.markers[[1]]),x="avg_logFC",y="p_val_adj", 
                selectLab = c("Apoe","Gzmb","Lgals1","Gzma","Tox","Ccr2","Cd52","Stat3","Ctsw","Ly6c2"), 
                pCutoff = NA,FCcutoff = NA,drawConnectors = T, title = "", labSize = 6,
                subtitle = NULL, legendPosition = "none", col = c("black","black","black","red"), ylab = "-Log[10](adj.p-val)", 
                xlab ="Log[2]FC(CaMKK2KO/WT)" ,caption = "") 
#c(rownames(t.g.markers[[1]] %>% slice_min(n = 35, order_by = p_val_adj)))
# Stem
t.g.markers[[2]][1:75,]
VlnPlot(r.tils,
        features = c("Ly6a","Stat1","Cd4","Irf7"), 
        idents = "Stem-Like TIL", 
        split.by = "Genotype", split.plot = T)

#Treg
t.g.markers[[3]][1:75,]
VlnPlot(r.tils,
        features = c("Ahr","Ccr2","Prf1","Cxcr6","Nfkbid","Ccl5","Gzmb","Klrg1","Foxp3","Cd83"), 
        idents = "Treg", 
        split.by = "Genotype", split.plot = T)

#Gamma-Delta TIL
t.g.markers[[5]][1:75,]
VlnPlot(r.tils,
        features = c("Apoe"), 
        idents = "Gamma-Delta TIL", 
        split.by = "Genotype", split.plot = T)
