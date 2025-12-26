#####scRNA-Seq analysis
library(vctrs)
library(pillar)
library(ggplot2)
library(Seurat)
library(dplyr)
library(tidydr)
library(data.table)

sclist<-list()
##GSE211388--------
count_data <- as.data.frame(fread('~/DATA/blca/GSE211388/GSE211388_count_matrix.txt.gz'))

sclist[[1]]<-CreateSeuratObject(count_data,  min.cells=5, min.features = 500)

cell_ids <- rownames(sclist2[[5]]@meta.data)
prefixes <- sub("^([^_]+_[^_]+).*", "\\1", cell_ids)
unique(prefixes)

sclist2[[5]]$orig.ident <- prefixes
table(sclist2[[5]]@meta.data$orig.ident)

#add sample id
sample_mapping <- c(
  'humanP_155' = 'GSM6468349',
  'humanP_157' = 'GSM6468350', 
  'humanP_158' = 'GSM6468351',
  'humanP_166' = 'GSM6468352',
  'humanP_167' = 'GSM6468353',
  'humanP_168' = 'GSM6468354',
  'humanP_170' = 'GSM6468355',
  'humanN_170' = 'GSM6468356',
  'humanN_171' = 'GSM6468357',
  'humanP_171' = 'GSM6468358',
  'humanN_1a' = 'GSM6468359',
  'humanP_1b' = 'GSM6468360',
  'humanN_2a' = 'GSM6468361',
  'humanP_2b' = 'GSM6468362'
)

# rename
sclist[[1]]$orig.ident <- sample_mapping[sclist[[1]]$orig.ident]
table(sclist[[1]]$orig.ident)



##GSE135337--------
sclis2<-list()

files <- dir(
  path = "~/DATA/blca/GSE135337", 
  recursive = TRUE, 
  pattern = "\\.(txt|xls)\\.gz$"  
)
files

input<-"~/DATA/blca/GSE135337"

# sample id
file_name_gsm <- sub("_.+", "", files)
file_name_gsm


full_files <- file.path(input, files)
full_files

for (i in seq_along(full_files)) {
  count_data <- as.data.frame(fread(full_files[i]))
  count_unique <- count_data[!duplicated(count_data$Symbol ), ] 
  
  rownames(count_unique) <- count_unique$Symbol   
  count_unique <- count_unique[, -c(1, 2)]
  
  seurat_obj <- CreateSeuratObject(
    counts = count_unique,
    min.cells = 5,   
    min.features = 500, 
  )
  
  seurat_obj$orig.ident <- rep(file_name_gsm[i], nrow(seurat_obj@meta.data))
  sclist2[[i]] <- seurat_obj
}
sclist2




##GSE129845--------
setwd('~/DATA/blca/GSE129845/')
dir = c("GSM3723357/","GSM3723358/","GSM3723359/")
names(dir) <-dir

sclist3 <- list()

for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i],gene.column = 2)
  sclist3[[i]] <- CreateSeuratObject(counts,  min.cells=5, min.features = 500)
}

#sample id
gsm_ids <- c('GSM3723357', 'GSM3723358', 'GSM3723359')

for (i in seq_along(sclist3)) {
  sclist3[[i]]$orig.ident <- rep(gsm_ids[i], nrow(sclist3[[i]]@meta.data))
}

combined_list <- c(sclist, sclist2, sclist3)
length(combined_list)

##merge----
scRNA<- merge(x = combined_list[[1]], y = combined_list[2:length(combined_list)],merge.data = TRUE)
scRNA <- JoinLayers(scRNA)

table(scRNA@meta.data$orig.ident)

#QC
scRNA[["mt_percent"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
scRNA$proj<-'BLCA'
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "mt_percent"), ncol = 3,group.by = 'proj')

minGene=500
maxGene=5000
mt=20

scRNA <- subset(scRNA, subset = nFeature_RNA > minGene & nFeature_RNA  < maxGene & mt_percent < mt)
scRNA

##NormalizeData
scRNA <- NormalizeData(scRNA) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% 
  ScaleData() %>% 
  RunPCA(npcs = 20,features = VariableFeatures(object = scRNA))

ElbowPlot(scRNA)

library(harmony)
scRNA <- RunHarmony(scRNA,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")


#RunUMAP
scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:20,reduction.name = "umap")


###DoubletFinder--------
library(DoubletFinder)
sweep.res.list <- paramSweep(scRNA, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
p<-as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK))
p

nExp_poi <- round(0.05*nrow(scRNA@meta.data))## Assuming 5% doublet formation rate - tailor for your dataset

scRNA <- doubletFinder(scRNA, PCs = 1:20, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)
#re-name
colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == "DF.classifications_0.25_5e-04_9189"] <- "doublet_info"
table(scRNA$doublet_info)

DimPlot(scRNA,group.by = 'doublet_info')

scRNA<-subset(scRNA, doublet_info %in% 'Singlet')

##CellCycle---
scRNA<-CellCycleScoring(scRNA,
                        s.features=s.genes,
                        g2m.features=g2m.genes,
                        set.ident=TRUE)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
scRNA<- CellCycleScoring(scRNA, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


##Re-normalize--------
scRNA <- NormalizeData(scRNA) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000) %>% 
  ScaleData(vars.to.regress = c("S.Score","G2M.Score")) %>% 
  RunPCA(npcs = 30,features = VariableFeatures(object = scRNA))
scRNA <-RunPCA(scRNA,npcs = 30,features = VariableFeatures(object = scRNA))

ElbowPlot(scRNA,ndims =30)

scRNA <- RunHarmony(scRNA,reduction = "pca",group.by.vars = "orig.ident",reduction.save = "harmony")


#RunUMAP
scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:30,reduction.name = "umap")

DimPlot(scRNA, reduction = "umap",group.by = "orig.ident")
FeaturePlot(scRNA,features = c('EPCAM','KRT18'))

scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.4)


Marker = list('TNK'=c('CD3D', 'CD3E','NKG7'),
              'B'=c('MS4A1', 'CD79A','CD19'),
              'Plas'=c('IGHG3', 'JCHAIN','IGLC2'),
              'Mac'=c('CD68', 'C1QA','CD163'),
              'Fib'=c('COL1A1','COL3A1','ACTA2'),
              'Epi'=c('EPCAM','KRT18','KRT13'),
              'EC'=c('VWF', 'RAMP2','PTPRB'))

DotPlot(scRNA,features=Marker,cols = c('gray','red')) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

#anno
scRNA$Celltype <-'other'
scRNA$Celltype[scRNA$RNA_snn_res.0.4 %in% c('0','1','5','6','8','10','14','15','16','17','18','19','20')] <- "Epithelial"
scRNA$Celltype[scRNA$RNA_snn_res.0.4 %in% c('2','7')] <- "TNK"
scRNA$Celltype[scRNA$RNA_snn_res.0.4 %in% c('4')] <- "Macrophage"
scRNA$Celltype[scRNA$RNA_snn_res.0.4 %in% c('9')] <- "B"
scRNA$Celltype[scRNA$RNA_snn_res.0.4 %in% c('12')] <- "Plasma"
scRNA$Celltype[scRNA$RNA_snn_res.0.4 %in% c('3','11')] <- "Fibroblast"
scRNA$Celltype[scRNA$RNA_snn_res.0.4 %in% c('13')] <- "Endothelial"

col<-c( "TNK"="#bcbddc","B cell"="#FD6467","Plasma"="#8DD3C7","Macrophage"= "#8968CD" ,"Fibroblast"="#9C964A" ,"Endothelial"="#c2e699",'Epithelial'="steelblue")

DimPlot(scRNA,label=T,cols = col,raster = F,label.box = T,repel = T,label.size = 4)+  theme_classic()+NoLegend()+ theme(
  panel.border = element_rect(
    color = "black",
    size = 1,        
    fill = NA 
  ),
  axis.line = element_blank(),
  legend.text = element_text(size = 12)
)

Marker = c('CD3D', 'CD3E','NKG7', 'CD79A','MS4A1','CD79B', 'JCHAIN','IGLC2','IGHG3', 'CD68','C1QA', 'CD163', 'COL1A1','COL3A1','ACTA2','VWF', 'RAMP2','PECAM1','KRT18','KRT13','EPCAM')

DotPlot(scRNA,features=Marker,cols = c('gray','red')) +RotatedAxis()+ theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"))+scale_colour_gradient2(trans = "reverse")+ylab('Celltype')+
  theme(axis.text.x = element_text(angle = 90,hjust = 1.1))

##rename
sample_mapping <- c(
  # Normal samples
  'GSM3723357' = 'N1',
  'GSM3723358' = 'N2', 
  'GSM3723359' = 'N3',
  'GSM5329919' = 'N4',
  # Tumor samples
  'GSM4006644' = 'T1',
  'GSM4006645' = 'T2',
  'GSM4006646' = 'T3',
  'GSM4006647' = 'T4',
  'GSM4006648' = 'T5',
  'GSM4751267' = 'T6',
  'GSM4751268' = 'T7',
  'GSM6468349' = 'T8',
  'GSM6468350' = 'T9',
  'GSM6468351' = 'T10',
  'GSM6468352' = 'T11',
  'GSM6468353' = 'T12',
  'GSM6468354' = 'T13',
  'GSM6468355' = 'T14',
  'GSM6468356' = 'T15',
  'GSM6468357' = 'T16',
  'GSM6468358' = 'T17',
  'GSM6468359' = 'T18',
  'GSM6468360' = 'T19',
  'GSM6468361' = 'T20',
  'GSM6468362' = 'T21'
)


scRNA$Sample <- sample_mapping[scRNA$orig.ident]

library(scRNAtoolVis)
cellRatioPlot(object = scRNA,
              sample.name = "Sample",
              celltype.name = "Celltype",
              flow.curve = 0.5,fill.col = col)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))


scRNA$Group<-'Tumor'
scRNA$Group[scRNA$Sample %in% c('N1','N2','N3','N4')]<-'Normal'


cellRatioPlot(object = scRNA,
              sample.name = "Group",
              celltype.name = "Celltype",
              flow.curve = 0.5,fill.col = col)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))

##cell ratio across stage
cal_table = function(x,y,prefix ){
  # x = phe$orig.ident
  # y = phe$celltype
  library(sur)
  library(reshape2)
  tbl =  table(x,y)
  pdf(paste0(prefix,'-table.pdf'),width = 10,height = 10)
  gplots::balloonplot( tbl )
  dev.off() 
  df = dcast(as.data.frame(tbl),x~y)
  head(df)
  write.csv(  df ,paste0(prefix,'-table.csv'))
  
  # ptb = round(sur::percent.table(x,y),2)
  ptb = round(100*tbl/rowSums(df[,-1]),2)
  
  pdf(paste0(prefix,'-percent-table.pdf'),width = 10,height = 10)
  gplots::balloonplot( ptb )
  dev.off()
  write.csv(  dcast(as.data.frame(ptb),x~y) ,paste0(prefix,'-percent-table.csv')) 
  
}
phe=scRNA@meta.data

cal_table(phe$orig.ident,phe$Celltype,prefix = 'celltype-vs-orig.ident')
cal_table(phe$Stage,phe$Celltype,prefix = 'celltype-vs-group')
cal_table(phe$Stage,phe$seurat_clusters,prefix = 'seurat_clusters-vs-group')
#cal_table(phe$Stage,phe$RNA_snn_res.0.8,prefix = 'RNA_snn_res.0.8-vs-group')

x='Celltype';y='Stage' 
plot_data <- data.frame(table(phe[, y ],
                              phe[, x ]))
head(plot_data)
plot_data$Total <- apply(plot_data,1,function(x)sum(plot_data[plot_data$Var1 == x[1],3]))
plot_data <- plot_data %>% mutate(Percentage = round(Freq/Total,3) * 100)
colnames(plot_data) <- c("Stage","Celltype","Freq","Total","Percentage") 
head(plot_data)
ggplot(plot_data,aes(Stage,Percentage,group = Celltype,color =Celltype))+geom_point(size=4)+
  geom_line(position = position_dodge(0.1),cex=1)+theme_bw()+theme_test(base_size = 16)  +
  scale_color_manual(values = col)


##Ro/e
library(Startrac)
in.dat<-scRNA@meta.data

R_oe<-calTissueDist(in.dat,
                    byPatient=F,
                    colname.cluster="Celltype",
                    colname.patient="orig.ident",
                    colname.tissue="Stage",
                    method="chisq",
                    min.rowSum=0)
R_oe

at = seq(0, 2, by = 1)

library(ComplexHeatmap)
Celltype =c('TNK','B', 'Plasma','Macrophage','Fibroblast','Endothelial','Epithelial' )
row_anno <- rowAnnotation(
  Celltype = Celltype,
  col = list(Celltype = setNames(col,My_levels)),
  show_legend = TRUE,
  annotation_name_gp = gpar(fontsize = 10)
)

color2 = colorRampPalette(c('#EDEDED', 'lightblue','#27408B'))(50)
Heatmap(as.matrix(R_oe),
        show_heatmap_legend = TRUE, 
        cluster_rows = F, 
        cluster_columns = F,
        row_names_side = 'right', 
        show_column_names = TRUE,
        show_row_names = TRUE,
        col = color2,
        right_annotation = row_anno,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(
          title = "Ro/e"
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          value <- R_oe[i, j]
          
          symbol <- if (value > 1) {
            "+++"
          } else if (value > 0.8 && value <= 1) {
            "++"
          } else if (value > 0.2 && value <= 0.8) {
            "+" 
          } else if (value > 0 && value <= 0.2) {
            "+/−"  
          } else if (value == 0) {
            "−"  
          }
         
          grid.text(symbol, x, y, gp = gpar(fontsize = 12, col = "black"))
        }
)




#####CellChat
scRNA$label<-'NC'
scRNA$label[scRNA$Stage %in% c('Ta','T1')]<-'Ta/T1'
scRNA$label[scRNA$Stage %in% c('T2','T3','T4')]<-'T2/T3/T4'
table(scRNA$label,scRNA$Stage)


library(CellChat)

groups <- list(
  'Low' = 'Ta/T1',
  'High' = 'T2/T3/T4'
)

for (group_name in names(groups)) {
  label_name <- groups[[group_name]]
  cat("Processing", group_name, "group (label:", label_name, ")\n")
  

    cellchat <- createCellChat(
    object = scRNAlist[[label_name]],                           
    meta = scRNAlist[[label_name]]@meta.data,                           
    group.by = "Celltype"
  )
  
  CellChatDB <- CellChatDB.human  
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  cellchat@DB <- CellChatDB.use
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = FALSE) 
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  saveRDS(cellchat, file = paste0('~/DATA/blca/chat/', group_name, '_cellchat.rds'))
  
  cat("Completed", group_name, "analysis\n\n")
}

col<-c( "TNK"="#bcbddc","B"="#FD6467","Plasma"="#8DD3C7","Macrophage"= "#8968CD" ,"Fibroblast"="#9C964A" ,"Endothelial"="#c2e699",'Epithelial'="steelblue")
cols = c("#009ACD","#CD5B45")

cellchat.L<-readRDS('~/DATA/blca/chat/Low_cellchat.rds')
cellchat.H<-readRDS('~/DATA/blca/chat/High_cellchat.rds')
object.list<-list(L=cellchat.L,H=cellchat.H)

cellchat<-mergeCellChat(object.list,add.names=names(object.list))
cellchat

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count",size.text =14,color.use = cols)
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight",size.text =14,color.use = cols)
p <- gg1 + gg2
p

um.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax,color.use = col,dot.alpha = 0.9)
}
patchwork::wrap_plots(plots = gg)

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use = cols)
gg1


netVisual_bubble(cellchat, sources.use = c(7), targets.use =c(1:6) ,  comparison = c(1, 2), angle.x = 45,color.text.use = T,color.text =cols,dot.size.min = 4,dot.size.max = 5)


# Chord diagram
pathway.union <- union(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
pathway.union

pathways.show='ANGPTL'

par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], 
                      signaling = pathways.show, 
                      layout = "chord", 
                      pt.title = 3, 
                      title.space = 0.05,
                      color.use = col, 
                      vertex.label.cex = 0.6, 
                      signaling.name = paste(pathways.show, names(object.list)[i]))
}

#circl
pathways.show <- c("ANGPTL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, point.size =5,layout = "circle",color.use = col,edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}




####Epithelial
col3<-c( "NC"="#B2DF8A", "Ta"="#A6CEE3","T1"="#FDBF6F", "T2"="#CAB2D6","T3"="#1F78B4","T4"="#6A3D9A")

Epi<-subset(scRNA, Celltype %in% 'Epithelial')
DimPlot(Epi,group.by = 'Stage',cols = col3)

Idents(Epi)<-Epi$Stage

Epi<-FindAllMarkers(Epi,only.pos = 0.25,logfc.threshold = 0.25,min.pct = 0.25)

DEG %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top10

averageHeatmap(Epi,
               group.by = 'Stage',
               markerGene = top10$gene,
               gene.order=top10$gene,
               annoCol = T,
               myanCol  = col3,
               fontsize =12)


library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

gid <- bitr(unique(DEG$gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')  
Genelist <- full_join(DEG, gid, by=c('gene' = 'SYMBOL'))

##GO:BP Enrichment
GOP = compareCluster(ENTREZID ~ cluster, 
                     data = Genelist, 
                     fun='enrichGO',
                     OrgDb = 'org.Hs.eg.db',ont='BP')


P<-dotplot(GOP, label_format=40,showCategory=5) + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  ggtitle('GO Enrichment')  +
  scale_fill_distiller(palette = "Blues", direction = -1)
P



####Metabolism
library(scMetabolism)

DimPlot(Epi,group.by = 'Stage',cols = col3)
countexp.Seurat<-sc.metabolism.Seurat(obj = Epi, 
                                      method = "VISION", 
                                      imputation = F, 
                                      ncores = 8, 
                                      metabolism.type = "KEGG")

metabolism <- countexp.Seurat@assays$METABOLISM$score



Epi[["METABOLISM"]] <- CreateAssayObject(counts = metabolism)
DefaultAssay(Epi) <- 'METABOLISM'

FeaturePlot(subMyeloid,features = 'Riboflavin metabolism')+ scale_color_viridis_c(option = "D")

FeaturePlot(subMyeloid,features = 'Glutathione metabolism')+ scale_color_viridis_c(option = "D")

##ANOVA test----
VlnPlot(mal, features = c('Riboflavin'), pt.size = 0, group.by = 'Stage', cols = col3) +
  geom_boxplot(width = .2, col = "black", fill = "white", outlier.shape = NA) +
  NoLegend() +  
  stat_compare_means(method = "anova", size = 6) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    text = element_text(size = 14, color = "black")
  ) +
  ylab('Metabolic activity') + 
  xlab('Stage') +
  ggtitle('Riboflavin')


VlnPlot(mal, features = c('Glutathione'), pt.size = 0, group.by = 'Stage', cols = col3) +
  geom_boxplot(width = .2, col = "black", fill = "white", outlier.shape = NA) +
  NoLegend() +  
  stat_compare_means(method = "anova", size = 6) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16, color = "black"),
    plot.title = element_text(size = 16, color = "black"),
    text = element_text(size = 14, color = "black")
  ) +
  ylab('Metabolic activity') + 
  xlab('Stage') +
  ggtitle('Glutathione')


##HEATMAP
mscore_data = data.frame(t(countexp.Seurat@assays[["METABOLISM"]][["score"]]),sce_Metal_exp$celltype)
avg_sM=aggregate(mscore_data[,1:ncol(mscore_data)-1],list(mscore_data$sce_Metal_exp.celltype),mean)
avg_sM=data.frame(t(avg_sM[,-1]))
avg_sM$KEGG = rownames(sce_Metal_exp@assays[["METABOLISM"]][["score"]])
rownames(avg_sM)=avg_sM$KEGG
avg_sM$KEGG<-NULL

colors = colorRampPalette(c("navy", "white", "firebrick3"))(100)
avg_sM<-avg_sM[1:40,]
pheatmap::pheatmap(avg_sM,show_colnames = T,scale='row',cluster_cols = F,color = colors)


##scFEA
## python
#https://github.com/changwn/scFEA/blob/master/scFEA_tutorial1.ipynb

## R
FeaturePlot(mal,features = 'Fatty Acid −> Acetyl−CoA') + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral")))

FeaturePlot(mal,features = 'Fatty Acid balance') + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral")))


FeaturePlot(mal,features = 'Glutathione −> glutamate') + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral")))
FeaturePlot(mal,features = 'Glutathione balance') + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral")))




library(Seurat)
library(igraph)
library(monocle)


table(Epi$Stage)
Idents(Epi)<-Epi$Stage

allCells=names(Idents(Epi))
allType = levels(Idents(Epi))

set.seed(123)
# downsample 25%
choose_Cells = unlist(lapply(allType, function(x) {
  cgCells = allCells[Idents(Epi) == x]
  numToSelect = floor(length(cgCells) * 0.25)
  cg = sample(cgCells, numToSelect)
  cg
}))


Episub = Epi[, allCells %in% choose_Cells]
Episub
table(Episub$Stage)


data<-GetAssayData(Episub,layer='counts')

#mycds
pd <- new('AnnotatedDataFrame',data = Episub@meta.data)
fData <- data.frame(gene_short_name = row.names(Episub), row.names = row.names(Episub))
fd <- new('AnnotatedDataFrame',data = fData )

mycds<- newCellDataSet(data,
                       phenoData = pd,
                       featureData = fd,
                       expressionFamily = negbinomial.size())
mycds


mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds,cores=8,relative_expr = TRUE)


Episub <- NormalizeData(Episub) %>% 
  FindVariableFeatures(selection.method = "vst",nfeatures = 3000)

Epideg<-FindAllMarkers(Episub,only.pos = T,
                       min.pct = 0.25,
                       logfc.threshold = 0.25)
table(Epideg$cluster)

Epideg %>%
  group_by(cluster) %>%
  top_n(n = 500, wt = avg_log2FC) -> top500

difgenes <-top500$gene


mycds <- setOrderingFilter(mycds,difgenes)
mycds <- reduceDimension(mycds, 
                         max_components = 2,
                         method ='DDRTree')


mycds<-orderCells(mycds)


plot_cell_trajectory(mycds, color_by="Pseudotime",cell_size = 1)+
  scale_color_viridis_c(option = "G")

plot_cell_trajectory(mycds, color_by="State",cell_size = 1)

mycds<-orderCells(mycds,root_state = 2)

plot_cell_trajectory(mycds, color_by="Pseudotime",cell_size = 1)+
  scale_color_viridis_c(option = "G")

plot_cell_trajectory(mycds, color_by="State",cell_size = 1)


plot_cell_trajectory(mycds, color_by="State",cell_size = 1)

plot_complex_cell_trajectory(mycds, x = 1, y = 2,color_by = "State")+ scale_color_manual(values =c('#1965b0','#7aafe3','#b579a2'))


library(ggpubr)
df <- pData(mycds) 
ggplot(df, aes(Pseudotime, colour = Stage, fill=Stage)) +
  geom_density(bw=0.5,size=1,alpha = 0.6)+theme_classic2()+ 
  scale_fill_manual(name= "", values = col3)+
  scale_color_manual(name = "", values = col3)



library(ggpubr)

mycom<-list(c("1", "2"),c("2", "3"),c('1','3'))


VlnPlot(Episub, features = c("Riboflavin"),
        group.by = 'State', 
        pt.size = 0,
        cols = c('#7aafe3','#1965b0','#b579a2')) +
  stat_compare_means(comparisons = mycom, label = "p.signif", size = 7,method = 't.test') 
ylim(0, 4.8) + 
  geom_boxplot(width = 0.1, col = "black", fill = NA, outlier.shape = NA)
theme_minimal() +
  theme(
    text = element_text(size = 18), 
    axis.text = element_text(size = 20),  
    legend.position = "none", 
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, size = 18) 
  )+ylab('Metabolic activity')+xlab('Stage')


Episub$Pseudotime<-mycds@phenoData@data$Pseudotime

cor.test(Episub$Riboflavin,Episub$Pseudotime,method = 'spearman')

df2<-Episub@meta.data
df2$Riboflavin_log <- log(df2$Riboflavin)
df2$Pseudotime_log <- log(df2$Pseudotime)



ggplot(df2, aes(df2$Riboflavin_log,df2$Pseudotime_log)) +       
  xlab('Riboflavin metabolic')+ylab('Pseudotime')+      
  geom_point()+ geom_smooth(method="lm",formula = Pseudotime_log ~ Riboflavin_log) + 
  theme_bw()+      
  stat_cor(method = 'spearman', aes(x =Riboflavin_log, y =Pseudotime_log))


ggplot(df2, aes(x = Riboflavin_log, y = Pseudotime_log)) +
  geom_point(alpha = 0.9, size = 1.5, color = "#7aafe3") +
  geom_smooth(method = "lm", 
              se = TRUE, 
              color = "#E63946", 
              fill = "#F1FAEE", 
              linetype = "solid") +
  stat_cor(
    method = "spearman", 
    label.x = 0.8, 
    label.y = 0.9,  
    color = "#1D3557",
    size = 5,  
    digits = 3  
  ) +
  labs(
    x = "Riboflavin metabolic",
    y = "Pseudotime"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = NA),
    plot.title = element_text(size = 16, face = "bold", color = "#1D3557"),
    plot.subtitle = element_text(size = 12, color = "#457B9D"),
    axis.title.x = element_text(size = 14, color = "#1D3557"),
    axis.title.y = element_text(size = 14, color = "#1D3557"),
    axis.text = element_text(size = 11, color = "#1D3557"),
    panel.grid.minor = element_blank()
  )


ribofl<-c('ACP1','FLAD1','BLVRB','ACP2','ACP5','RFK','ENPP1','ENPP3')

plot_pseudotime_heatmap(mycds[ribofl,], 
                        num_clusters=1,
                        show_rownames=T, 
                        return_heatmap=T, 
                        hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62))


Time_diff <- differentialGeneTest(mycds[difgenes,], cores = 1, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")

Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()

#order
Time_genes <- Time_diff %>%
  arrange(qval) %>% 
  slice_head(n = 300) %>%  
  pull(gene_short_name) %>% 
  as.character() %>%
  unique() 

p<-plot_pseudotime_heatmap(mycds[Time_genes,], 
                           num_clusters=3, 
                           show_rownames=F, 
                           return_heatmap=T, 
                           hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62))
p$tree_row
clusters <- cutree(p$tree_row, k = 3)
clustering <- data.frame(clusters)
clustering[,1] <-as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
clustering$Gene=rownames(clustering)



## Cluster_GO
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

gid <- bitr(unique(clustering$Gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')
Genelist <- full_join(clustering, gid, by=c('Gene' = 'SYMBOL'))

##GO:BP Enrichment
GOP = compareCluster(ENTREZID ~ Gene_Clusters, 
                     data = Genelist, 
                     fun='enrichGO',
                     OrgDb = 'org.Hs.eg.db') 


dotplot(GOP, label_format=40,showCategory=20) + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  ggtitle('GO Enrichment') +  
  scale_color_gradientn(colours = colorRampPalette(colors = brewer.pal(11, name = "RdBu"))(11))




#visium
input_data_dir <- "~/DATA/blca/st"
sample_list <- list.files(input_data_dir)
sample_list
sample_list<-sample_list[5:8]
sample_list

# file path
samples_dir <- sample_list %>% file.path(input_data_dir, .)
samples_dir

##slice id
sample_names <- sample_list
sample_names


sample_objects <- purrr::map(1:length(sample_list), function(x) {
  ## read data
  one_dir <- samples_dir[x]
  sample_id <- sample_list[x]
  slice_id <- sample_names[x]
  sample_object <- Load10X_Spatial(
    data.dir = one_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = slice_id,
    filter.matrix = TRUE
  )
  sample_object@project.name <- sample_id
  sample_object@meta.data$orig.ident <- slice_id
  sample_object <- RenameCells(object = sample_object, add.cell.id = slice_id)
  
  return(sample_object)
})

sample_objects <- lapply(sample_objects, 
                         SCTransform, variable.features.n = 2000,
                         assay = "Spatial", 
                         method = "poisson")

#visium 2
##P5
count<-read.csv('~/DATA/blca/st/P5/GSM7853987_BLCA-B1-count.csv.gz',row.names = 1)

P5<-CreateSeuratObject(counts = t(count),assay = 'Spatial')

position = read.csv('~/DATA/blca/st/P5/GSM7853987_BLCA-B1-Metadata.csv',row.names = 1)
position=position[6:7]
colnames(position)<-c('x','y')
P5@images$P5 =  new(  Class = 'SlideSeq',  assay = "Spatial",  key = "image_",  coordinates = position)

SpatialPlot(P5)

P5<-SCTransform(P5, assay = "Spatial", method = "poisson",variable.features.n = 2000)
SpatialPlot(P5,features = 'EPCAM')

##P6
count<-read.csv('~/DATA/blca/st/P6/GSM7853988_BLCA-B2-count.csv.gz',row.names = 1)

P6<-CreateSeuratObject(counts = t(count),assay = 'Spatial')

position = read.csv('~/DATA/blca/st/P6/GSM7853988_BLCA-B2-Metadata.csv',row.names = 1)
position=position[6:7]
colnames(position)<-c('x','y')
P6@images$P6 =  new(  Class = 'SlideSeq',  assay = "Spatial",  key = "image_",  coordinates = position)
SpatialPlot(P6)

P6<-SCTransform(P6, assay = "Spatial", method = "poisson",variable.features.n = 2000)
SpatialPlot(P6,features = 'EPCAM')


sample_objects[[5]]<-P5
sample_objects[[5]]$orig.ident<-'P5'
sample_objects[[6]]<-P6
sample_objects[[6]]$orig.ident<-'P6'


##merge-----
ST<-merge(sample_objects[[1]], y = sample_objects[2:6])
SpatialPlot(ST)

VariableFeatures(ST) <- unique(unlist(lapply(sample_objects, VariableFeatures)))

ST<-RunPCA(ST, assay = "SCT", verbose = FALSE)
ElbowPlot(ST,ndims = 50)

library(harmony)
table(ST$orig.ident)
ST <- RunHarmony(ST,
                 reduction = "pca",
                 group.by.vars = "orig.ident",r
                 eduction.save = "harmony")

ST <- RunUMAP(ST, reduction="harmony",dims = 1:20)

ST <- FindNeighbors(ST, reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.4)
DimPlot(ST, reduction = "umap",cols = cols)
table(ST$seurat_clusters)

cols <- c(
  '0'="#1f77b4", '1'="#ff7f0e",'2' ="#2ca02c", '3'="#d62728",
  '4'= "#9467bd",'5'="#8c564b", '6'="#e377c2", '7'="#7f7f7f", 
  '8'="#E7298A", '9'= "#B2DF8A", '10'="#55B1B1",'11'="#85d2e1")

plot_gene = function (cluster){
  p1 = SpatialPlot(ST, group.by = cluster,stroke = 0, pt.size.factor =4.2,images = 'P1')& scale_fill_manual(values = cols)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)
  p2 = SpatialPlot(ST, group.by = cluster,stroke = 0,pt.size.factor = 5,images = 'P2')& scale_fill_manual(values = cols)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)
  p3 = SpatialPlot(ST, group.by = cluster,stroke = 0,pt.size.factor = 4,images = 'P3')& scale_fill_manual(values = cols)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)
  p4 = SpatialPlot(ST, group.by = cluster,stroke = 0,pt.size.factor = 4,images = 'P4')& scale_fill_manual(values = cols)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)
  p5 = SpatialPlot(ST, group.by = cluster,stroke = 0,pt.size.factor = 2.5,images = 'P5')& scale_fill_manual(values = cols)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)
  p6 = SpatialPlot(ST, group.by = cluster,stroke = 0, pt.size.factor =2.5,images = 'P6')& scale_fill_manual(values = cols)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)
  ggarrange(p1,p2,p3,p4,p5,p6,ncol  = 6,nrow = 1,common.legend = T,legend = 'bottom')
}

plot_gene(cluster='seurat_clusters')


#RCTD deconvolution
library(spacexr)
DimPlot(scRNA, reduction = "umap",group.by = 'Cellyype',label = T)

#sc-counts
sc_counts<-GetAssayData(scRNA,layer='counts')
Celltype<-as.factor(scRNA$Celltype)
names(Celltype)<-colnames(scRNA)

reference<-Reference(sc_counts,Celltype)


# run RCTD
for(i in 1:6) {
  cat("Processing sample", i, "...\n")
  
  spatial_count <- GetAssayData(sample_objects[[i]], layer = 'counts')
  
  coords <- GetTissueCoordinates(sample_objects[[i]])[, 1:2]
  
  query <- SpatialRNA(coords, spatial_count)
  
  RCTD <- create.RCTD(query, reference, max_cores = 1)
  RCTD <- run.RCTD(RCTD, doublet_mode = "full")
  
  saveRDS(RCTD, file = paste0('~/DATA/blca/st/rctd/RCTDP', i, '.full.rds.gz'))
  
  cat("Completed sample", i, "\n")
}


paths <- list.files("~/DATA/blca/st/rctd")
file_paths = grep('full', paths, value=TRUE)
file_paths = file_paths %>% file.path("~/DATA/blca/st/rctd", .)
file_paths


norm_weights_list <- list()

for (file_path in file_paths) {
  
  myRCTD <- readRDS(file_path)
  
  barcodes <- colnames(myRCTD@spatialRNA@counts)
  weights <- myRCTD@results$weights
  norm_weights <- normalize_weights(weights)
  
  file_name <- basename(file_path)
  norm_weights_list[[file_name]] <- norm_weights
}


merged_df <- data.frame()
sample_name <-c('P1_','P2_','P3_','P4_','P5_','P6_')

for (i in 1:length(file_paths)) {
  file_name <- basename(file_paths[i])
  prefix <- sample_name[i]
  
  df <- as.data.frame(norm_weights_list[[file_name]])
  rownames(df) <- paste0(prefix, rownames(df))
  
  merged_df <- rbind(merged_df, df)
}

merged_df$cb=rownames(merged_df)
merged_df$cb <- ifelse(
  grepl("^.[6]", rownames(merged_df)), 
  paste0(merged_df$cb, "_6"),
  merged_df$cb
)

merged_df$cb <- ifelse(
  grepl("6$", merged_df$cb),          
  sub("^...", "", merged_df$cb), 
  merged_df$cb
)

rownames(merged_df)<-merged_df$cb

ST<-AddMetaData(ST,metadata = merged_df)

SpatialPlot(ST,features = 'Epithelial')


##Cluster-Celltype Correlation
library('tidyr')
library("tidyverse")
library("patchwork")

merged_df$cb=NULL
cell_prop = merged_df %>% 
  as.data.frame() %>%
  rownames_to_column("spot_id") %>%
  pivot_longer(-spot_id)

cell_prop$id<-cell_prop$spot_id

#head(cell_prop)
niche_info = ST@meta.data %>% as.data.frame() %>% 
  rownames_to_column("spot_id") %>%
  select_at(c("spot_id", "orig.ident", "seurat_clusters")) %>%
  dplyr::rename("Niche" = "seurat_clusters") %>%
  mutate(Niche = paste0("Cluster_", Niche))


cellprops_info =  cell_prop  %>%
  left_join(niche_info, c("spot_id" = "spot_id")) %>%
  na.omit()

#head(cellprops_info)
cell_props_summary_CT_pat <- cellprops_info %>%
  group_by(name, Niche) %>%
  summarize(median_CT = median(value))
head(cell_props_summary_CT_pat)

# Check niches that are unique to some patients
cell_prop <-  cell_prop %>%
  left_join(niche_info, c("spot_id" = "spot_id")) %>%
  na.omit() %>%
  dplyr::select(-spot_id) %>%
  group_by(name) %>%
  nest() %>%
  mutate(wres = map(data, function(dat) {
    
    niches <- dat$Niche %>%
      unique() %>%
      set_names()
    
    map(niches, function(g) {
      
      test_data <- dat %>%
        mutate(test_group = ifelse(.data[["Niche"]] == g,
                                   "target", "rest")) %>%
        mutate(test_group = factor(test_group,
                                   levels = c("target", "rest")))
      
      wilcox.test(value ~ test_group, 
                  data = test_data,
                  alternative = "greater") %>%
        broom::tidy()
    }) %>% enframe("Niche") %>%
      unnest()
    
  }))

wilcox_types <- cell_prop %>%
  dplyr::select(wres) %>%
  unnest() %>%
  ungroup() %>%
  dplyr::mutate(adj_pval = p.adjust(p.value)) %>%
  dplyr::mutate(log_adj_pval = -log10(adj_pval)) %>%
  dplyr::mutate(sign = ifelse(adj_pval < 0.005, "*", ""))

#head(wilcox_types)
ct_median_desc <- cell_prop %>%
  dplyr::select(data) %>%
  unnest() %>%
  group_by(name, Niche) %>%
  summarise(median_prop = median(value)) %>%
  mutate(scaled_median_prop = (median_prop - mean(median_prop))/sd(median_prop))
dim(ct_median_desc)

niche_car_df <- left_join(wilcox_types, ct_median_desc) %>% na.omit()
dim(niche_car_df)

ct_median_desc_mat <- niche_car_df %>% dplyr::select(name, Niche, scaled_median_prop) %>%
  pivot_wider(names_from = Niche, values_from = scaled_median_prop) %>%
  column_to_rownames("name") %>%
  as.matrix()

ct_sign_desc_mat <- niche_car_df %>% dplyr::select(name, Niche, sign) %>%
  pivot_wider(names_from = Niche, values_from = sign) %>%
  column_to_rownames("name") %>%
  as.matrix()

ct_median_desc_mat
ct_sign_desc_mat


library("ComplexHeatmap")
library(RColorBrewer)
col3=colorRampPalette(rev(brewer.pal(9, "PRGn")))(62)
Heatmap(ct_median_desc_mat, name = "scaled comp", 
        col = col3,
        rect_gp = gpar(col = "black", lwd = 1),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf(ct_sign_desc_mat[i, j]), x, y, gp = gpar(fontsize = 14))
        })


##cnv----
library(infercnv)

anno<- data.frame(Idents(ST))
counts_matrix = GetAssayData(ST, slot="counts")

##infercnv_obj
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix=counts_matrix,
  annotations_file=anno,
  ref_group_names=c("NormalEpi"),
  gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
  min_max_counts_per_cell = c(100, +Inf),
  delim="\t")

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir='~/DATA/blca/st/inferCNV',  
                             cluster_by_groups=T, 
                             analysis_mode='subclusters', 
                             denoise=F,
                             HMM=F, 
                             num_threads=16,
                             write_expr_matrix=T, 
                             tumor_subcluster_partition_method = "leiden")

expr <- infercnv_obj@expr.data

expr2=expr-1
expr2=expr2 ^ 2
CNV_score=as.data.frame(colMeans(expr2))
colnames(CNV_score)="CNV_score"
CNV_score$CB=rownames(CNV_score)


meta<-C8@meta.data[,c(1,7)]
meta$CB=rownames(meta)
CNV_score=CNV_score%>%inner_join(meta,by="CB")

CNV_score %>% 
  ggplot(aes(seurat_clusters, CNV_score)) +
  geom_violin(aes(fill = seurat_clusters), color = "NA") +
  scale_fill_manual(values = cols) +
  facet_wrap(~ orig.ident, scales = "free_y") +  # Y轴独立
  theme_bw()


ST$CB=rownames(ST@meta.data)
ST@meta.data$CNV_score <- CNV_score$CNV_score[match(ST@meta.data$CB, CNV_score$CB)]

ST$Region<-'nMal'
ST$Region<-as.character(ST$Region)

ST$Region[ST$seurat_clusters %in% c('1','3','7') & ST$orig.ident %in% c('P1')]<-'Mal'
ST$Region[ST$seurat_clusters %in% c('0','3','4','9','10','11') & ST$orig.ident %in% c('P2')]<-'Mal'
ST$Region[ST$orig.ident %in% c('P3')]<-'Mal'
ST$Region[ST$seurat_clusters %in% c('0','3','5','8','11') & ST$orig.ident %in% c('P3')]<-'nMal'

ST$Region[ST$seurat_clusters %in% c('1','3','6','7','8','9','10','11') & ST$orig.ident %in% c('P4')]<-'Mal'
ST$Region[ST$seurat_clusters %in% c('0','3') & ST$orig.ident %in% c('P5')]<-'Mal'
ST$Region[ST$seurat_clusters %in% c('3','5') & ST$orig.ident %in% c('P6')]<-'Mal'

Mal<-subset(ST, ST$Region %in% 'Mal')
plot_gene = function (cluster){
  p1 = SpatialPlot(Mal, features  = cluster,stroke = 0, pt.size.factor =4.2,images = 'P1',max.cutoff = cutoff)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)& scale_fill_viridis_c(option = "A")
  p2 = SpatialPlot(Mal, features = cluster,stroke = 0,pt.size.factor = 5,images = 'P2',max.cutoff = cutoff)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)& scale_fill_viridis_c(option = "A")
  p3 = SpatialPlot(Mal, features = cluster,stroke = 0,pt.size.factor = 4,images = 'P3',max.cutoff = cutoff)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)& scale_fill_viridis_c(option = "A")
  p4 = SpatialPlot(Mal, features = cluster,stroke = 0,pt.size.factor = 4,images = 'P4',max.cutoff = cutoff)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)& scale_fill_viridis_c(option = "A")
  p5 = SpatialPlot(Mal, features = cluster,stroke = 0,pt.size.factor = 2.5,images = 'P5',max.cutoff = cutoff)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)& scale_fill_viridis_c(option = "A")
  p6 = SpatialPlot(Mal, features = cluster,stroke = 0, pt.size.factor =2.5,images = 'P6',max.cutoff = cutoff)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)& scale_fill_viridis_c(option = "A")
  ggarrange(p1,p2,p3,p4,p5,p6,ncol  = 3,nrow = 2,common.legend = T,legend = 'bottom')
}
plot_gene(cluster='CNV_score')


countexp.Seurat <- scmetabolism(obj = Mal2,
                                method = "VISION", 
                                imputation = F, 
                                ncores = 2, 
                                metabolism.type = "KEGG")

metabolism <- as.data.frame(t(countexp.Seurat@assays$METABOLISM$score))
Mal$Riboflavin<-metabolism$`Riboflavin metabolism`

plot_gene = function (cluster){
  p1 = SpatialPlot(Mal, features  = cluster,stroke = 0, pt.size.factor =4.2,images = 'P1',max.cutoff = cutoff)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)& scale_fill_viridis_c(option = "A")
  p2 = SpatialPlot(Mal, features = cluster,stroke = 0,pt.size.factor = 5,images = 'P2',max.cutoff = cutoff)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)& scale_fill_viridis_c(option = "A")
  p3 = SpatialPlot(Mal, features = cluster,stroke = 0,pt.size.factor = 4,images = 'P3',max.cutoff = cutoff)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)& scale_fill_viridis_c(option = "A")
  p4 = SpatialPlot(Mal, features = cluster,stroke = 0,pt.size.factor = 4,images = 'P4',max.cutoff = cutoff)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)& scale_fill_viridis_c(option = "A")
  p5 = SpatialPlot(Mal, features = cluster,stroke = 0,pt.size.factor = 2.5,images = 'P5',max.cutoff = cutoff)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)& scale_fill_viridis_c(option = "A")
  p6 = SpatialPlot(Mal, features = cluster,stroke = 0, pt.size.factor =2.5,images = 'P6',max.cutoff = cutoff)&theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),  aspect.ratio = 1)& scale_fill_viridis_c(option = "A")
  ggarrange(p1,p2,p3,p4,p5,p6,ncol  = 3,nrow = 2,common.legend = T,legend = 'bottom')
}
plot_gene(cluster='Riboflavin')


df=Mal@meta.data

ggplot(df, aes(x = log(Riboflavin), y = log(CNV_score))) +
  geom_point(alpha = 0.9, size = 1.5, color = "#7aafe3") +
  geom_smooth(method = "lm", se = TRUE, color = "#E63946", fill = "#F1FAEE", linetype = "solid") +
  stat_cor(
    method = "spearman", 
    label.x = 0.8, 
    label.y = 0.4,  
    color = "#1D3557",
    size = 5,  
    digits = 3  
  ) +
  labs(
    x = "Riboflavin metabolic",
    y = "CNV score"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.background = element_rect(fill = NA),
    plot.title = element_text(size = 16, face = "bold", color = "#1D3557"),
    plot.subtitle = element_text(size = 12, color = "#457B9D"),
    axis.title.x = element_text(size = 14, color = "#1D3557"),
    axis.title.y = element_text(size = 14, color = "#1D3557"),
    axis.text = element_text(size = 11, color = "#1D3557"),
    panel.grid.minor = element_blank()
  )+
  coord_cartesian(ylim = c(NA, -6)) 


Idents(ST)<-ST$Region
ST<-PrepSCTFindMarkers(ST)

DEG<-FindMarkers(ST,ident.1 = 'Mal',ident.2 = 'nMal',only.pos=F,min.pct = 0.1, logfc.threshold =0.25)

dif=data.frame(
  symbol=rownames(DEG),
  log2FoldChange=DEG$avg_log2FC,
  padj=DEG$p_val_adj
)

VolcanoPlot(dif, padj=0.05, title="Mal vs nMal", label.max = 8)

#KEGG
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

DEG$Gene=rownames(DEG)
gid <- bitr(unique(DEG$Gene), 'SYMBOL', 'ENTREZID', OrgDb= 'org.Hs.eg.db')

KEGG = enrichKEGG(gene  = gid$ENTREZID,pvalueCutoff = 1) 

library(RColorBrewer)
dotplot(KEGG, label_format=40,showCategory=20) + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  ggtitle('KEGG Enrichment') +  
  scale_color_gradientn(colours = colorRampPalette(colors = brewer.pal(11, name = "RdBu"))(11))

selected_pathways <- c("Apoptosis",'Estrogen signaling pathway','Cell adhesion molecules','PI3K-Akt signaling pathway','p53 signaling pathway','ECM-receptor interaction','Cornified envelope formation', "HIF-1 signaling pathway", "Glycolysis / Gluconeogenesis",'Glutathione metabolism','Carbon metabolism','Nucleotide metabolism','Oxidative phosphorylation','Biosynthesis of amino acids','Pyruvate metabolism','Phenylalanine metabolism')



subKEGG<-subset(KEGG@result,Description %in% selected_pathways & KEGG@result$pvalue < 0.05)

write.csv(subKEGG,file = '~/DATA/blca/st/sub_KEGG.csv')