library(Seurat); library(dplyr); library(ggplot2); library(patchwork); library(harmony); library(DoubletFinder)

#################
### 1.Expression/expressions ###
dir()
# [1] "CK-1"                              "CK-1.expression.demo.xls"         
# [3] "CK-2"                              "CK-2.expression.demo.xls"         
# [5] "Mr-DIV1-24h-1"                     "Mr-DIV1-24h-1.expression.demo.xls"
# [7] "Mr-DIV1-24h-2"                     "Mr-DIV1-24h-2.expression.demo.xls"
# [9] "Mr-DIV1-4h-1"                      "Mr-DIV1-4h-1.expression.demo.xls" 
#[11] "Mr-DIV1-4h-2"                      "Mr-DIV1-4h-2.expression.demo.xls" 

#################
dir("CK-1")
#[1] "barcodes.tsv"   "expression.xls" "genes.tsv"      "matrix.mtx"

## doubletfinder
data_name = c("CK-1","CK-2","Mr-DIV1-4h-1","Mr-DIV1-4h-2","Mr-DIV1-24h-1","Mr-DIV1-24h-2")
filter.cells = NULL
for (i in seq(data_name)) {
	mat <- Read10X(data.dir = paste0("rawdata/",data_name[i]), gene.column = 1)
	object.list <- CreateSeuratObject(counts = mat, project = data_name[i], assay = "RNA")
	object <- RenameCells(object.list, add.cell.id = data_name[i])
	object@misc$vars.regress <- NULL
	object <- NormalizeData(object, normalization.method = "LogNormalize",scale.factor = 10000)
	object <- FindVariableFeatures(object, selection.method = "vst",nfeatures = 2000)
	object <- ScaleData(object = object, vars.to.regress = object@misc$vars.regress,features = rownames(object))
	object <- RunPCA(object,assay = "RNA",npcs = 50,features = VariableFeatures(object), verbose = FALSE)
	object <- RunUMAP(object, dims = 1:50, umap.method = "uwot",n.neighbors = 30,n.components = 2L , reduction.name = "umap",reduction = "pca")
	sweep.res.list <- paramSweep_v3(object, PCs = 1:50, sct = FALSE)
	sweep.stats <- summarizeSweep(sweep.res.list)
	bcmvn <- find.pK(sweep.stats)
	pK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)[1]]))
	rate <- 7.6 * 10^-6 * ncol(object) + 5.27 * 10^-4
	nExp_poi <- round(as.numeric(rate) * ncol(object))
	object <- doubletFinder_v3(object, PCs = 1:50, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
	colnames(object@meta.data)[grep('pANN', colnames(object@meta.data))] <- "pANN"
	colnames(object@meta.data)[grep('DF.classifications', colnames(object@meta.data))] <- "classifications"
	doublet <- colnames(object)[(object@meta.data$classifications %in% "Doublet")]
	filter.cells <- c(filter.cells,doublet)
}

## Create Seurat object
data_name = c("CK-1","CK-2","Mr-DIV1-4h-1","Mr-DIV1-4h-2","Mr-DIV1-24h-1","Mr-DIV1-24h-2")
object.list <- list()
for (i in seq(data_name)) {
	mat <- Read10X(data.dir = data_name[i], gene.column = 1)
	object.list[[i]] <- CreateSeuratObject(counts = mat, project = data_name[i], assay = "RNA")
}

## merge Seurat object
object <- merge(x = object.list[[1]], y = object.list[-1],add.cell.ids = data_name)

## filtering
object <- subset(object, subset = nCount_RNA < 11000)
object <- subset(object, subset = nFeature_RNA > 500 & nFeature_RNA < 2600)

cells.use <- setdiff(colnames(object), filter.cells)
object <- object[, cells.use]

## Normalize 
object@misc$vars.regress <- c("nCount_RNA")
object <- NormalizeData(object, normalization.method = "LogNormalize",scale.factor = 10000)
object <- FindVariableFeatures(object, selection.method = "vst",nfeatures = 2000)
object <- ScaleData(object = object, vars.to.regress = object@misc$vars.regress,features = rownames(object))

## integration
object <- RunPCA(object,assay = "RNA",npcs = 50,features = VariableFeatures(object), verbose = FALSE)
object <- RunHarmony(object, group.by.vars = "orig.ident", project.dim = FALSE, assay.use = DefaultAssay(obj))
object <- RunTSNE(object, dims = 1:50,dim.embed = 2L , perplexity = 30, reduction.name = "tsne",reduction="harmony")
object <- RunUMAP(object, dims = 1:50, umap.method = "uwot",n.neighbors = 30,n.components = 2L , reduction.name = "umap",reduction = "harmony")

## Find clusters
object <- FindNeighbors(object, reduction = "harmony", dims = 1:50,force.recalc = TRUE)
object <- FindClusters(object, resolution = 0.5,algorithm = 1, temp.file.location = getwd())

## Draw t-SNE UMAP plot
p1 <- DimPlot(object, reduction = 'tsne', group.by = "orig.ident", cols=rainbow(nlevels(object@meta.data$orig.ident)), label = FALSE)
p2 <- DimPlot(object, reduction = 'tsne', group.by = "seurat_clusters", cols=rainbow(nlevels(object@meta.data$seurat_clusters)), label = TRUE)
ggsave(p1+p2, file = "TSNE.pdf")
p1 <- DimPlot(object, reduction = 'umap', group.by = "orig.ident", cols=rainbow(nlevels(object@meta.data$orig.ident)), label = FALSE)
p2 <- DimPlot(object, reduction = 'umap', group.by = "seurat_clusters", cols=rainbow(nlevels(object@meta.data$seurat_clusters)), label = TRUE)
ggsave(p1+p2, file = "UMAP.pdf")

## Save data object 
DefaultAssay(object) <- "RNA"
save(object, file = "obj.Rda")

## Find positive markers 
Idents(object) <- "seurat_clusters"
obj.markers <- FindAllMarkers(object = object, only.pos = TRUE,
    min.pct = 0.25, logfc.threshold = 0.25,
    return.thresh = 0.01, pseudocount.use = 0,base=exp(1))
save( obj.markers, file = "markers.Rda" )

## group diff
FindGroupMarkers <- function (object, group.data, contrast = NULL,return.thresh = 0.05, use.qvalue = TRUE,
	logfc.threshold = 0.25, test.use = "MAST", min.pct = 0.1, is.bulk = TRUE, ...) {
	
	cells.1 <- rownames(group.data)[group.data[[contrast[[1]]]]]
	cells.2 <- rownames(group.data)[group.data[[contrast[[2]]]]]
	data.1 <- CalAvgExp(object[, cells.1], is.return = TRUE, is.bulk = is.bulk) %>%
		as.data.frame() %>% tibble::rownames_to_column(var = "gene") %>%
		reshape2::melt(id.vars = "gene", variable.name = "cluster", value.name = "mean.1") %>%
		mutate(cluster = as.factor(cluster))
	data.2 <- CalAvgExp(object[, cells.2], is.return = TRUE, is.bulk = is.bulk) %>%
		as.data.frame() %>% tibble::rownames_to_column(var = "gene") %>%
		reshape2::melt(id.vars = "gene", variable.name = "cluster", value.name = "mean.2") %>%
		mutate(cluster = as.factor(cluster))
	pct.1 <- CalPctExp(object[, cells.1], is.return = TRUE, is.bulk = is.bulk) %>% round(digits = 3) %>%
		as.data.frame() %>% tibble::rownames_to_column(var = "gene") %>%
		reshape2::melt(id.vars = "gene", variable.name = "cluster", value.name = "pct.1") %>%
		mutate(cluster = as.factor(cluster))
	pct.2 <- CalPctExp(object[, cells.2], is.return = TRUE, is.bulk = is.bulk) %>% round(digits = 3) %>%
		as.data.frame() %>% tibble::rownames_to_column(var = "gene") %>%
		reshape2::melt(id.vars = "gene", variable.name = "cluster", value.name = "pct.2") %>%
		mutate(cluster = as.factor(cluster))
	c1 <- table(Idents(object[,cells.1]))
	c2 <- table(Idents(object[,cells.2]))
	c1 <- c(bulk = length(cells.1), c1)
	c2 <- c(bulk = length(cells.2), c2)
	data <- full_join(x = data.1, y = data.2, by = c("gene", "cluster")) %>%
		left_join(y = pct.1, by = c("gene", "cluster") ) %>%
		left_join(y = pct.2, by = c("gene", "cluster") ) %>%
		mutate(cells.1 = c1[cluster], cells.2 = c2[cluster]) %>%
		mutate(cluster = factor(cluster, levels = if ( is.bulk ) c("bulk", levels(object)) else levels(object))) %>%
		mutate(avg_logFC = log2(mean.2/mean.1)) %>%
		mutate(p_val = 1, p_val_adj = 1, sig = "nosig")
	data[is.na(data)] <- 0
	idents.all <- c("bulk", levels(droplevels(object)))
	data <- data %>% filter(cluster %in% idents.all)
	features.use <- data %>% filter(mean.1 + mean.2 > 0 & (pct.1 > min.pct | pct.2 > min.pct) & abs(avg_logFC) > logfc.threshold / log(2))
	features.use <- split(features.use, features.use$cluster)
	features.use <- lapply(features.use, function(x) as.character(x$gene))
	genes.de <- list()
	for (i in seq(idents.all)) {
		cells.idents <- if ( idents.all[[i]] == "bulk" ){
			WhichCells(object)
		} else {
				WhichCells(object, idents = idents.all[[i]])
		}
		cells.1.tmp <- intersect(cells.2, cells.idents)
		cells.2.tmp <- intersect(cells.1, cells.idents)
		genes.de[[i]] <- tryCatch(expr = {
			FindMarkers(object = object, ident.1 = cells.1.tmp, ident.2 = cells.2.tmp,
				assay = NULL,features = features.use[[idents.all[[i]]]], logfc.threshold = 0,
				test.use = test.use, slot = "data", min.pct = min.pct,
				min.diff.pct = -Inf, verbose = TRUE,
				only.pos = FALSE, max.cells.per.ident = Inf,
				random.seed = 1, latent.vars = NULL,
				min.cells.feature = 3, min.cells.group = 3,
				pseudocount.use = 0)}, error = function(cond) { return(cond$message) })
	}
	gde.all <- data.frame()
	for (i in seq(idents.all)) {
		gde <- genes.de[[i]]
		gde <- gde[order(gde$p_val_adj, gde$p_val, -gde[, 2]), ]
		gde <- subset(x = gde, subset = p_val_adj < return.thresh)
		genes.de[[i]]$sig <- "nosig"
		genes.de[[i]][rownames(gde)[gde$avg_logFC > 0], "sig"] <- "up"
		genes.de[[i]][rownames(gde)[gde$avg_logFC < 0], "sig"] <- "down"
		genes.de[[i]]$cluster <- idents.all[i]
		genes.de[[i]]$gene <- rownames(x = genes.de[[i]])
		gde.all <- rbind(gde.all, genes.de[[i]])
	}
	gde.all <- gde.all %>% select(cluster, gene, p_val, p_val_adj, sig) %>%
		full_join(y = data, by = c("gene", "cluster")) %>%
		mutate( p_val     = if_else(is.na(p_val.x), p_val.y, p_val.x),
			p_val_adj = if_else(is.na(p_val_adj.x), p_val_adj.y, p_val_adj.x),
			sig       = if_else(is.na(sig.x), sig.y, sig.x)) %>%
		mutate( sig     = factor(sig, levels = c("up", "nosig", "down")),
			cluster = factor(cluster, levels = levels(data$cluster)),
			contrast  = factor(paste0(contrast, collapse = "-vs-")) ) %>%
		select(-p_val.x, -p_val.y, -p_val_adj.x, -p_val_adj.y, -sig.x, -sig.y)
	return(gde.all)
}

sample_col = "orig.ident"
group.data <- object@meta.data %>% select(sample = !! sample_col)
for ( i in unique(group.data$sample)) {
	group.data <- group.data %>% mutate(!! i := sample %in% i )
}
by_sample = list("CK"=c("CK-1", "CK-2"),"Mr-DIV1-24h"=c("Mr-DIV1-24h-1" ,"Mr-DIV1-24h-2"),"Mr-DIV1-4h"=c("Mr-DIV1-4h-1", "Mr-DIV1-4h-2"))
for ( i in names(by_sample) ){
	group.data <- group.data %>% mutate(!! i := sample %in% by_sample[[i]])
}
group.data <- group.data %>% select(- sample)

differ_groups   = list(c("CK", "Mr-DIV1-4h"),c("CK" ,"Mr-DIV1-24h"),c("Mr-DIV1-4h", "Mr-DIV1-24h"))
log2fc = 0.36
logfc.threshold = log2fc / log2(exp(1))
marker <- list()
for ( i in differ_groups ) {
	name <- paste(i, collapse = "-vs-")
	marker[[name]] <- FindGroupMarkers(object, group.data, contrast = i,
		return.thresh = 0.05, use.qvalue = TRUE, is.bulk = TRUE,
		logfc.threshold = logfc.threshold, test.use = "MAST",min.pct = 0.1)
}
# output result
for ( contrast_name in names(marker) ) {
	contrast <- strsplit(contrast_name, "-vs-")[[1]]
	data <- marker[[contrast_name]] %>% 
		select(Cluster = cluster, GeneID = gene,
			!! contrast[[1]] := mean.1, !! contrast[[2]] := mean.2,
			cells.1, cells.2, pct.1, pct.2,
			log2FC = avg_logFC, p_value = p_val, p_val_adjust = p_val_adj, significance = sig)
	write.table(data, file = paste0("DifferMarker.", contrast_name, ".all.xls"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
}


sessionInfo()
#R version 3.5.1 (2018-07-02)
#Platform: x86_64-conda_cos6-linux-gnu (64-bit)
#Running under: CentOS release 6.9 (Final)

#Matrix products: default

# 
#locale:
# [1] LC_CTYPE=en_US.iso885915       LC_NUMERIC=C                  
# [3] LC_TIME=en_US.iso885915        LC_COLLATE=en_US.iso885915    
# [5] LC_MONETARY=en_US.iso885915    LC_MESSAGES=en_US.iso885915   
# [7] LC_PAPER=en_US.iso885915       LC_NAME=C                     
# [9] LC_ADDRESS=C                   LC_TELEPHONE=C                
#[11] LC_MEASUREMENT=en_US.iso885915 LC_IDENTIFICATION=C  
# 
#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     
# 
#other attached packages:
#[1] patchwork_1.0.0.9000 ggplot2_3.2.1        dplyr_0.8.5         
#[4] Seurat_3.1.1        
# 
#loaded via a namespace (and not attached):
# [1] nlme_3.1-147        tsne_0.1-3          bitops_1.0-6       
# [4] RcppAnnoy_0.0.14    RColorBrewer_1.1-2  httr_1.4.0         
# [7] sctransform_0.2.0   tools_3.5.1         R6_2.4.1           
#[10] irlba_2.3.3         KernSmooth_2.23-15  uwot_0.1.5         
#[13] lazyeval_0.2.2      colorspace_1.4-1    npsurv_0.4-0       
#[16] withr_2.1.2         tidyselect_1.0.0    gridExtra_2.3      
#[19] compiler_3.5.1      cli_2.0.2           plotly_4.8.0       
#[22] caTools_1.17.1.2    scales_1.0.0        lmtest_0.9-36      
#[25] ggridges_0.5.1      pbapply_1.4-0       stringr_1.4.0      
#[28] digest_0.6.25       R.utils_2.8.0       pkgconfig_2.0.3    
#[31] htmltools_0.5.1.1   bibtex_0.4.2        htmlwidgets_1.3    
#[34] rlang_1.0.6         zoo_1.8-5           jsonlite_1.6       
#[37] ica_1.0-2           gtools_3.8.1        R.oo_1.22.0        
#[40] magrittr_1.5        Matrix_1.2-17       Rcpp_1.0.4.6       
#[43] munsell_0.5.0       fansi_0.4.1         ape_5.3            
#[46] reticulate_1.11.1   R.methodsS3_1.7.1   stringi_1.4.3      
#[49] gbRd_0.4-11         MASS_7.3-51.3       gplots_3.0.1.1     
#[52] Rtsne_0.15          plyr_1.8.4          grid_3.5.1         
#[55] parallel_3.5.1      gdata_2.18.0        listenv_0.7.0      
#[58] ggrepel_0.8.1       crayon_1.3.4        lattice_0.20-41    
#[61] cowplot_0.9.4       splines_3.5.1       SDMTools_1.1-221   
#[64] pillar_1.4.3        igraph_1.2.4        future.apply_1.3.0 
#[67] reshape2_1.4.3      codetools_0.2-16    leiden_0.3.1       
#[70] glue_1.4.0          lsei_1.2-0          metap_1.1          
#[73] data.table_1.12.2   RcppParallel_4.4.4  png_0.1-7          
#[76] Rdpack_0.10-1       gtable_0.3.0        RANN_2.6.1         
#[79] purrr_0.3.4         tidyr_0.8.3         future_1.15.1      
#[82] assertthat_0.2.1    rsvd_1.0.2          survival_2.44-1.1  
#[85] viridisLite_0.3.0   tibble_2.1.3        cluster_2.0.7-1    
#[88] globals_0.12.4      fitdistrplus_1.0-14 ROCR_1.0-7         

