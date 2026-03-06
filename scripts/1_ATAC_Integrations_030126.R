#### library
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(harmony)
library(ArchR)
library(hexbin)
library(biovizBase)

library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(patchwork)
library(statmod)
library(data.table)
library(foreach)



            #### functions
            
            # removes doublet using DoubletFinder
            Remove_Doublet=function(seurat_objt=EFE001, pseu_rate = 0.04){
              
              DefaultAssay(seurat_objt) <- "RNA"
              
              # standard
              seurat_objt <- NormalizeData(seurat_objt)
              seurat_objt <- FindVariableFeatures(seurat_objt, selection.method = "vst", nfeatures = 2000)
              seurat_objt <- ScaleData(seurat_objt)
              seurat_objt <- RunPCA(seurat_objt) # ElbowPlot(seurat_objt)
              seurat_objt <- RunUMAP(seurat_objt, dims = 1:20)
              
              sweep.res.list_seut <- paramSweep_v3(seurat_objt, PCs = 1:20, sct = TRUE)
              sweep.stats_seut <- summarizeSweep(sweep.res.list_seut, GT = FALSE)
              bcmvn_seut <- find.pK(sweep.stats_seut)
              
              # script for visualizaion of pK parameter selection
              pK=as.numeric(as.character(bcmvn_seut$pK))
              BCmetric=bcmvn_seut$BCmetric
              pK_choose = pK[which(BCmetric %in% max(BCmetric))]
              
              nExp_poi <- round(pseu_rate*nrow(seurat_objt@meta.data))  ## Assuming 4% doublet formation rate - tailor for your dataset
              seurat_objt <- doubletFinder_v3(seurat_objt, pN = 0.25, pK = pK_choose, nExp = nExp_poi, PCs = 1:20) ##update newly
              
              # plot
              DF.name = colnames(seurat_objt@meta.data)[grepl("DF.classification", colnames(seurat_objt@meta.data))]
              plts = cowplot::plot_grid(ncol = 1, DimPlot(seurat_objt, group.by = DF.name) + NoAxes())
              
              
              # remove doublets
              seurat_objt=seurat_objt[, seurat_objt@meta.data[, DF.name] == "Singlet"]
              return(list(objt = seurat_objt, plts = plts))
            }
            
            # removes Mt and Rb genes from the seurat objects
            Remove_MtRb=function(seurat_objt=EFE001){
              
              mito.genes <- grep(pattern = "^MT-", x = rownames(seurat_objt), value = FALSE)
              ribo.genes <- grep(pattern = "^RP[SL]", x = rownames(seurat_objt), value = FALSE)
              retained <- rownames(seurat_objt)[-c(mito.genes, ribo.genes) ]
              seurat_objt[["RNA"]] <- subset(seurat_objt[["RNA"]], features = retained)
              
              return(seurat_objt)
            }
            
            # reads in disease RNA and ATAC data
            Read_Disease=function(counts = counts_efe1, fragpath = fragpath_efe1, projname = "EFE001", annt=annotation){
              dat_BU = CreateSeuratObject(
                counts = counts$`Gene Expression`,
                assay = "RNA",
                project=projname
              )
              
              dat_BU[["ATAC"]] <- CreateChromatinAssay(
                counts = counts$Peaks,
                sep = c(":", "-"),
                fragments = fragpath,
                annotation = annt
              )
              
              DefaultAssay(dat_BU) <- "ATAC"
              dat_BU <- NucleosomeSignal(dat_BU)
              dat_BU <- TSSEnrichment(dat_BU)
              
              DefaultAssay(dat_BU) <- "RNA"
              dat_BU[["percent.mt"]] <- PercentageFeatureSet(dat_BU, pattern = "^MT-")
              dat_BU[['percent.ribo']] <- PercentageFeatureSet(dat_BU, pattern = "^RP[SL]")
              
              return(dat_BU)
            }
            
            # reads in control RNA data
            Read_Control=function(path = y1, projname = "young1"){
              ctrl=Read10X(data.dir = path)
              ctrl_BU <- CreateSeuratObject(counts = ctrl,min.cells = 3, min.features = 200,project = projname)
              
              ctrl_BU[["percent.mt"]] <- PercentageFeatureSet(ctrl_BU, pattern = "^MT-")
              ctrl_BU[['percent.ribo']] <- PercentageFeatureSet(ctrl_BU, pattern = "^RP[SL]")
              
              return(ctrl_BU)
            }
            
            # process of macs2 called peaks for each sample stored in the same folder
            Peak_Macs2=function(macs.path = macs_path, macs.name = name, peaks.macs2 = peaks_macs2){
              
              
              # convert MACS2 narrowPeaks to Granges object
              
              
              for(i in 1:length(macs.name)){
                df <- read.table(file = paste0(macs.path,macs.name[i], "_peaks.narrowPeak"), 
                                 col.names = c("chr","start", "end", "name", "score", "strand", "fold_change", 
                                               "neg_log10pvalue_summit", "neg_log10qvalue_summit","relative_summit_position"))
                gr <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE, 
                                               starts.in.df.are.0based = TRUE)
                
                peaks.macs2[[macs.name[i]]] <- gr
              }
              
              #peaks.macs2
              
              # remove peaks on nonstandard chromosomes
              peaks.macs2 <- lapply(peaks.macs2, function(x){
                x <- keepStandardChromosomes(x, pruning.mode = "coarse")
              })
              
              # encode a genome blacklist regions for hg38 using the reference list from https://github.com/Boyle-Lab/Blacklist/
              blacklist_df <- fread(paste0(macs.path,"hg38-blacklist.v2.bed"))
              blacklist <- makeGRangesFromDataFrame(blacklist_df, ignore.strand = T, seqnames.field = "V1", start.field = "V2", end.field = "V3")
              
              seqlevelsStyle(blacklist) <- 'UCSC'
              
              # remove peaks in genomic blacklist regions
              peaks.macs2 <- lapply(peaks.macs2, function(x){
                x <- subsetByOverlaps(x, ranges = blacklist, invert = TRUE)
              })
              
              
              return(peaks.macs2)
            }
            
            # creating peaks assay for the seurat objects using a combined peak set from processed macs2
            Create_Peaks=function(seurat_objt=EFE001, fragpath = fragpath_efe1, comb.peaks = combined.peaks, annt=annotation){
              
              DefaultAssay(seurat_objt) <- "ATAC"
              
              # quantify counts in each peak
              macs2_counts <- FeatureMatrix(
                fragments = Fragments(seurat_objt),
                features = comb.peaks,
                cells = colnames(seurat_objt)
              )
              
              # create a new assay using the MACS2 peak set and add it to the Seurat object
              seurat_objt[["peaks"]] <- CreateChromatinAssay(
                counts = macs2_counts,
                fragments = fragpath,
                annotation = annt
              )
              
              
              return(seurat_objt)
            }
            
            #### This is preprocessing -----------------------------------------------------------------------------------------------------------------------
            
            #### 1. load the RNA and ATAC data
            
            # get gene annotations for hg38
            annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
            
            ###change the internet to guest, NOT TCH;  very important
            seqlevelsStyle(annotation) <- "UCSC"
            
            #### EFE samples
            counts_efe1 <- Read10X_h5("/Users/yangyu/Desktop/Harvard_ChenLab/EP/EFE_data/EFE001/filtered_feature_bc_matrix.h5")
            fragpath_efe1 <- "/Users/yangyu/Desktop/Harvard_ChenLab/EP/EFE_data/EFE001/atac_fragments.tsv.gz"
            
            counts_efe2 <- Read10X_h5("/Users/yangyu/Desktop/Harvard_ChenLab/EP/EFE_data/EFE002/filtered_feature_bc_matrix.h5")
            fragpath_efe2 <- "/Users/yangyu/Desktop/Harvard_ChenLab/EP/EFE_data/EFE002/atac_fragments.tsv.gz"
            
            counts_efe3 <- Read10X_h5("/Users/yangyu/Desktop/Harvard_ChenLab/EP/EFE_data/EFE003/filtered_feature_bc_matrix.h5")
            fragpath_efe3 <- "/Users/yangyu/Desktop/Harvard_ChenLab/EP/EFE_data/EFE003/atac_fragments.tsv.gz"
            
            #### PuCtrl samples
            counts_puctrl1 <- Read10X_h5("/Users/yangyu/Desktop/Harvard_ChenLab/EP/PuLabCtrl_data/PuCtrl001/filtered_feature_bc_matrix.h5")
            fragpath_puctrl1 <- "/Users/yangyu/Desktop/Harvard_ChenLab/EP/PuLabCtrl_data/PuCtrl001/atac_fragments.tsv.gz"
            
            counts_puctrl2 <- Read10X_h5("/Users/yangyu/Desktop/Harvard_ChenLab/EP/PuLabCtrl_data/PuCtrl002/filtered_feature_bc_matrix.h5")
            fragpath_puctrl2 <- "/Users/yangyu/Desktop/Harvard_ChenLab/EP/PuLabCtrl_data/PuCtrl002/atac_fragments.tsv.gz"
            
            counts_puctrl3 <- Read10X_h5("/Users/yangyu/Desktop/Harvard_ChenLab/EP/PuLabCtrl_data/PuCtrl003/filtered_feature_bc_matrix.h5")
            fragpath_puctrl3 <- "/Users/yangyu/Desktop/Harvard_ChenLab/EP/PuLabCtrl_data/PuCtrl003/atac_fragments.tsv.gz"
            
            # load the RNA and ATAC data
            EFE001 = Read_Disease(counts = counts_efe1, fragpath = fragpath_efe1, projname = "EFE001",annt=annotation)
            EFE002 = Read_Disease(counts = counts_efe2, fragpath = fragpath_efe2, projname = "EFE002",annt=annotation)
            EFE003 = Read_Disease(counts = counts_efe3, fragpath = fragpath_efe3, projname = "EFE003",annt=annotation)
            
            Ctrl001 = Read_Disease(counts = counts_puctrl1, fragpath = fragpath_puctrl1, projname = "Ctrl001",annt=annotation)
            Ctrl002 = Read_Disease(counts = counts_puctrl2, fragpath = fragpath_puctrl2, projname = "Ctrl002",annt=annotation)
            Ctrl003 = Read_Disease(counts = counts_puctrl3, fragpath = fragpath_puctrl3, projname = "Ctrl003",annt=annotation)
            
            # 2.2.1 ......................................................................................................................................... 
            #### removing mt and rb genes 
            EFE001 = Remove_MtRb(seurat_objt = EFE001)
            EFE002 = Remove_MtRb(seurat_objt = EFE002)
            EFE003 = Remove_MtRb(seurat_objt = EFE003)
            Ctrl001 = Remove_MtRb(seurat_objt = Ctrl001)
            Ctrl002 = Remove_MtRb(seurat_objt = Ctrl002)
            Ctrl003 = Remove_MtRb(seurat_objt = Ctrl003)
            
            
            # 2.2.2 .........................................................................................................................................
            # convert MACS2 narrowPeaks to Granges object
            macs_path="/Users/yangyu/Desktop/Harvard_ChenLab/EP/EFE_data/Peak_macs2/"
            name <- c("efe001", "efe002","efe003","ctrl001","ctrl002","ctrl003")
            peaks_macs2 <- list()
            
            peaks_macs2 = Peak_Macs2(macs.path = macs_path, macs.name = name, peaks.macs2 = peaks_macs2)
            
            # Create a unified set of peaks to quantify in each dataset
            combined.peaks <- GenomicRanges::reduce(x = c(peaks_macs2[["efe001"]], peaks_macs2[["efe002"]], peaks_macs2[["efe003"]],peaks_macs2[["ctrl001"]], peaks_macs2[["ctrl002"]], peaks_macs2[["ctrl003"]]))
            
            # Filter out bad peaks based on length
            peakwidths <- width(combined.peaks)
            combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
            combined.peaks
            
            # create a new assay using the MACS2 peak set and add it to the Seurat object
            
            EFE001 = Create_Peaks(seurat_objt=EFE001, fragpath = fragpath_efe1, comb.peaks = combined.peaks, annt=annotation)
            EFE002 = Create_Peaks(seurat_objt=EFE002, fragpath = fragpath_efe2, comb.peaks = combined.peaks, annt=annotation)
            EFE003 = Create_Peaks(seurat_objt=EFE003, fragpath = fragpath_efe3, comb.peaks = combined.peaks, annt=annotation)
            
            Ctrl001 = Create_Peaks(seurat_objt=Ctrl001, fragpath = fragpath_puctrl1, comb.peaks = combined.peaks, annt=annotation)
            Ctrl002 = Create_Peaks(seurat_objt=Ctrl002, fragpath = fragpath_puctrl2, comb.peaks = combined.peaks, annt=annotation)
            Ctrl003 = Create_Peaks(seurat_objt=Ctrl003, fragpath = fragpath_puctrl3, comb.peaks = combined.peaks, annt=annotation)
            
            # gene activities
            DefaultAssay(EFE001) <- DefaultAssay(EFE002) <- DefaultAssay(EFE003) <- DefaultAssay(Ctrl001) <- DefaultAssay(Ctrl002) <- DefaultAssay(Ctrl003) <- "peaks"
            efe1_gene.activities <- GeneActivity(EFE001)
            efe2_gene.activities <- GeneActivity(EFE002)
            efe3_gene.activities <- GeneActivity(EFE003)
            ctrl1_gene.activities <- GeneActivity(Ctrl001)
            ctrl2_gene.activities <- GeneActivity(Ctrl002)
            ctrl3_gene.activities <- GeneActivity(Ctrl003)
            
            # add gene activities as a new assay - unnormalized
            EFE001[["ACTIVITY"]] <- CreateAssayObject(counts = efe1_gene.activities)
            EFE002[["ACTIVITY"]] <- CreateAssayObject(counts = efe2_gene.activities)
            EFE003[["ACTIVITY"]] <- CreateAssayObject(counts = efe3_gene.activities)
            Ctrl001[["ACTIVITY"]] <- CreateAssayObject(counts = ctrl1_gene.activities)
            Ctrl002[["ACTIVITY"]] <- CreateAssayObject(counts = ctrl2_gene.activities)
            Ctrl003[["ACTIVITY"]] <- CreateAssayObject(counts = ctrl3_gene.activities)
            
            #### subset
            EFE001 <- subset(
              x = EFE001,
              subset = nCount_RNA < 7853 &
                nCount_RNA > 527 &
                nFeature_RNA > 300 &
                nCount_ATAC < 6094 &
                nCount_ATAC > 89 &
                percent.mt < 20 &
                percent.ribo < 10 &
                nucleosome_signal < 2 &
                TSS.enrichment > 1
            )
            
            EFE002 <- subset(
              x = EFE002,
              subset = nCount_RNA < 11374 &
                nCount_RNA > 582 &
                nFeature_RNA > 300 &
                nCount_ATAC < 8813 &
                nCount_ATAC > 114 &
                percent.mt < 20 &
                percent.ribo < 10 &
                nucleosome_signal < 2 &
                TSS.enrichment > 1
            )
            
            EFE003 <- subset(
              x = EFE003,
              subset = nCount_RNA < 5567 &
                nCount_RNA > 624 &
                nFeature_RNA > 300 &
                nCount_ATAC < 2456 &
                nCount_ATAC > 153 &
                percent.mt < 20 &
                percent.ribo < 10 &
                nucleosome_signal < 2 &
                TSS.enrichment > 1
            )
            
            Ctrl001 <- subset(
              x = Ctrl001,
              subset = nCount_RNA < 6535 &
                nCount_RNA > 590 &
                nFeature_RNA > 300 &
                nCount_ATAC < 19215 &
                nCount_ATAC > 250 &
                percent.mt < 20 &
                percent.ribo < 10 &
                nucleosome_signal < 4 &
                TSS.enrichment > 1
            )
            
            Ctrl002 <- subset(
              x = Ctrl002,
              subset = nCount_RNA < 4691 &
                nCount_RNA > 647 &
                nFeature_RNA > 300 &
                nCount_ATAC < 10919 &
                nCount_ATAC > 104 &
                percent.mt < 20 &
                percent.ribo < 10 &
                nucleosome_signal < 4 &
                TSS.enrichment > 1
            )
            
            Ctrl003 <- subset(
              x = Ctrl003,
              subset = nCount_RNA < 4786 &
                nCount_RNA > 516 &
                nFeature_RNA > 300 &
                nCount_ATAC < 12764 &
                nCount_ATAC > 102 &
                percent.mt < 20 &
                percent.ribo < 10 &
                nucleosome_signal < 4 &
                TSS.enrichment > 1
            )
            
            # go back to doublet cells
            library(DoubletFinder)
            
            EFE001.sg = Remove_Doublet(seurat_objt = EFE001,pseu_rate = 0.040)
            EFE002.sg = Remove_Doublet(seurat_objt = EFE002,pseu_rate = 0.075)
            EFE003.sg = Remove_Doublet(seurat_objt = EFE003,pseu_rate = 0.100)
            Ctrl001.sg = Remove_Doublet(seurat_objt = Ctrl001,pseu_rate = 0.040)
            Ctrl002.sg = Remove_Doublet(seurat_objt = Ctrl002,pseu_rate = 0.023)
            Ctrl003.sg = Remove_Doublet(seurat_objt = Ctrl003,pseu_rate = 0.016)
            
            #### 3. Integrate samples------------------------------------------------------------------------------------------------------------------------------------------------------------
            # normalize and identify variable features for each dataset independently
            EFE1<-EFE001.sg$objt
            EFE2<-EFE002.sg$objt
            EFE3<-EFE003.sg$objt
            CTRL1<-Ctrl001.sg$objt
            CTRL2<-Ctrl002.sg$objt
            CTRL3<-Ctrl003.sg$objt
            
            EFE001.sg <- EFE002.sg <- EFE003.sg <- EFE001 <- EFE002 <- EFE003 <- Ctrl001.sg <- Ctrl002.sg <- Ctrl003.sg <- Ctrl001 <- Ctrl002 <- Ctrl003 <- NULL
            gc()
            
            # save
            
            #saveRDS(EFE1,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/atac_integrated/atac_only/files/EFE1_adb_raw_100124.rds")
            #saveRDS(EFE2,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/atac_integrated/atac_only/files/EFE2_adb_raw_100124.rds")
            #saveRDS(EFE3,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/atac_integrated/atac_only/files/EFE3_adb_raw_100124.rds")
            #saveRDS(CTRL1,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/atac_integrated/atac_only/files/Ctrl1_adb_raw_100124.rds")
            #saveRDS(CTRL2,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/atac_integrated/atac_only/files/Ctrl2_adb_raw_100124.rds")
            #saveRDS(CTRL3,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/atac_integrated/atac_only/files/Ctrl3_adb_raw_100124.rds")
            
            
            # subset TSS > 3
            
            
            
            # atac compute LSI
            
            # This ordering is essential
            all.list <- list(EFE1, EFE2, EFE3, CTRL1, CTRL2, CTRL3)
            all.list[[1]][["STIM"]] <- "EFE1"
            all.list[[2]][["STIM"]] <- "EFE2"
            all.list[[3]][["STIM"]] <- "EFE3"
            all.list[[4]][["STIM"]] <- "CTRL1"
            all.list[[5]][["STIM"]] <- "CTRL2"
            all.list[[6]][["STIM"]] <- "CTRL3"
            
            for(i in 1:length(all.list)){
              DefaultAssay(all.list[[i]]) <- "peaks"
            }
            
            all.list <- lapply(X = all.list, FUN = function(x) {
              x <- RunTFIDF(x)
              x <- FindTopFeatures(x, min.cutoff = 10)
              x <- RunSVD(x)
            })
            
            # merge
            atac.merged <- merge(all.list[[1]], y = c(all.list[2:6]))
            atac.merged
            
            # process the merged dataset
            atac.merged <- RunTFIDF(atac.merged)
            atac.merged <- FindTopFeatures(atac.merged, min.cutoff = 10)
            atac.merged <- RunSVD(atac.merged)
            atac.merged <- RunUMAP(atac.merged, reduction = "lsi", dims = 2:30)
            p1 <- DimPlot(atac.merged, group.by = "STIM")
            
            # integrate ATAC-seq samples
            # find integration anchors
            integration.anchors <- FindIntegrationAnchors(
              object.list = all.list,
              anchor.features = rownames(all.list[[1]]),
              reduction = "rlsi",
              dims = 2:30
            )
            
            # integrate LSI embeddings
            
            combined.atac <- IntegrateEmbeddings(
              anchorset = integration.anchors,
              reductions = atac.merged[["lsi"]],
              new.reduction.name = "integrated_lsi",
              dims.to.integrate = 1:30
            )
            
            # create a new UMAP using the integrated embeddings
            combined.atac <- RunUMAP(combined.atac, reduction = "integrated_lsi", dims = 2:30)
            p2 <- DimPlot(combined.atac, group.by = "STIM")
            
            (p1 + ggtitle("Merged")) | (p2 + ggtitle("Integrated"))
            
            
            # save
            # saveRDS(combined.atac,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/atac_integrated/atac_only/files/Integrated6samp_atac_raw_092524.rds")
            
            
            
            
            
            # readin rna data This step tries to get annotation information from RNA data, therefore, projecting the query atac onto the rna reference is the way to go
            # "Integrated12samp_combined_refinedannt_121923.rds" is the same as "Integrated12samp_combined_grandannt_010424.rds", with extra annotations
            # matching cells
            rna.data = readRDS("/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/rna_only/inte_12samp/Integrated12samp_combined_refinedannt_121923.rds")
            rna.data = subset(rna.data, subset = SAMPLE %in% c("EFEs","CTRLs"))
            
            # transfer integrated RNA assay to object with integrated ATAC assay
            combined.atac@reductions
            combined.atac@assays
            
            combined.atac[["integrated_RNA"]] <- rna.data[["integrated"]]
            
            
            # combing wnn
            DefaultAssay(combined.atac) <- "integrated_RNA"
            
            # Run the standard workflow for visualization and clustering for RNA
            combined.atac <- ScaleData(combined.atac, verbose = FALSE)
            combined.atac <- RunPCA(combined.atac, npcs = 30, verbose = FALSE)
            combined.atac <- RunUMAP(combined.atac, reduction = "pca", dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
            
            # for atac
            DefaultAssay(combined.atac) <- "peaks"
            combined.atac <- RunTFIDF(combined.atac)
            combined.atac <- FindTopFeatures(combined.atac, min.cutoff = 'q0')
            combined.atac <- RunSVD(combined.atac)
            combined.atac <- RunUMAP(combined.atac, reduction = 'integrated_lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
            
            combined.atac@reductions
            
            combined.atac <- FindMultiModalNeighbors(combined.atac, reduction.list = list("pca", "integrated_lsi"), dims.list = list(1:30, 2:30))
            combined.atac <- RunUMAP(combined.atac, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
            combined.atac <- FindClusters(combined.atac, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
            
            p1 <- DimPlot(combined.atac, reduction = "umap.atac", label = TRUE, label.size = 5, repel = TRUE, group.by = "STIM") + ggtitle("ATAC")
            p2 <- DimPlot(combined.atac, reduction = "wnn.umap", label = TRUE, label.size = 5, repel = TRUE, group.by = "STIM") + ggtitle("WNN")
            p1 + p2 & theme(plot.title = element_text(hjust = 0.5))
            
            
            
            # normalize gene activities
            DefaultAssay(combined.atac) <- "ACTIVITY"
            combined.atac <- NormalizeData(combined.atac)
            combined.atac <- ScaleData(combined.atac, features = rownames(combined.atac))
            
            # Identify anchors
            rna.data <- NormalizeData(rna.data)
            rna.data <- FindVariableFeatures(rna.data, selection.method = "vst", nfeatures = 2000)
            rna.data <- ScaleData(rna.data, verbose = FALSE)
            
            transfer.anchors <- FindTransferAnchors(reference = rna.data, query = combined.atac, features = VariableFeatures(object = rna.data),
                                                    reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")
            
            #transfer.anchors <- FindTransferAnchors(reference = rna.data, query = combined.atac, reduction = "cca")
            
            celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna.data$refined_type,
                                                 weight.reduction = combined.atac[["lsi"]], dims = 2:30)
            
            combined.atac <- AddMetaData(combined.atac, metadata = celltype.predictions)
            
            #combined.atac$annotation_correct <- combined.atac$predicted.id == combined.atac$cell_type
            p1 <- DimPlot(combined.atac, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
            p2 <- DimPlot(combined.atac, group.by = "cell_type", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
            p1 | p2
            
            # suggested coembeddings
            rna.data$tech <- "RNA"
            combined.atac$tech <- "ATAC"
            genes.use = VariableFeatures(rna.data)
            refdata = GetAssayData(rna.data, assay = "RNA", slot = "data")[genes.use, ]
            imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata,
                                       weight.reduction = combined.atac[["lsi"]], dims = 2:30)
            
            # This line adds the imputed data matrix to combined.atac
            combined.atac[["RNA"]] <- imputation
            coembed = merge(x = rna.data, y = combined.atac)
            
            # pca umap
            coembed <- ScaleData(coembed, features = genes.use, do.scale = FALSE)
            coembed <- RunPCA(coembed, features = genes.use, verbose = FALSE)
            coembed <- RunUMAP(coembed, reduction = "pca", dims = 1:30) 
            
            # visualize
            coembed$refined_type <- ifelse(!is.na(coembed$refined_type),coembed$refined_type,coembed$predicted.id)
            
            plot3 <- DimPlot(coembed, group.by = "tech", label = T, repel = T)
            plot4 <- DimPlot(coembed, group.by = "refined_type", label = T, repel = T)
            
            plot3|plot4
            
            #saveRDS(coembed,file="/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/atac_integrated/atac_only/files/Coembed_rna_atac_raw100124.rds")
            
            coembed <- FindNeighbors(coembed)
            coembed <- FindClusters(coembed, resolution = 0.5)
            DimPlot(coembed, reduction = "umap", split.by = "STIM", pt.size = 0.01,label = T) + NoLegend()
            
            DefaultAssay(coembed) <- "ACTIVITY"
            FeaturePlot(coembed, reduction = "umap", features = c("COL1A1","DCN","PECAM1","TTN"))
            DefaultAssay(coembed) <- "RNA"
            FeaturePlot(coembed, reduction = "umap", features = c("COL1A1","DCN","PECAM1","TTN"))
            
            # save processed seurat objects THIS IS THE SAME FILE AS Integrated6samp_atac_raw_092524.rds
            # saveRDS(combined.atac, file = "/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/atac_integrated/atac_only/files/combined_seurat_object_03152022.rds")
            
            
            
            # Assuming your integrated Seurat object is called 'rna.data'
            pca_embeddings <- Embeddings(rna.data, reduction = "pca")
            umap_embeddings <- Embeddings(rna.data, reduction = "umap")
            
            # Assuming your new Seurat object is called 'combined.atac'
            combined.atac[["pcas"]] <- CreateDimReducObject(embeddings = pca_embeddings, key = "PCS_", assay = DefaultAssay(combined.atac))
            combined.atac[["umaps"]] <- CreateDimReducObject(embeddings = umap_embeddings, key = "UMAPS_", assay = DefaultAssay(combined.atac))
            
            # Plot UMAP for the new dataset to verify
            DimPlot(combined.atac, reduction = "umaps", group.by = "seurat_clusters", label = T)
            
            
            
            # Assuming your integrated Seurat object is called 'rna.data'
            # Extract the cell type information from the rna.data dataset
            rna.data <- SetIdent(rna.data, value = "refined_type")
            cell_types <- Idents(rna.data)
            
            # Assuming your new Seurat object is called 'new_dataset'
            # Make sure the cell names match between the two datasets
            common_cells <- intersect(colnames(combined.atac), names(cell_types))
            
            # Transfer cell types to the new dataset
            combined.atac <- AddMetaData(object = combined.atac, metadata = cell_types[common_cells], col.name = "predicted_celltype")
            
            # Optionally, set the identities (Idents) based on the transferred cell types
            Idents(combined.atac) <- combined.atac@meta.data$predicted_celltype
            
            # Visualize UMAP or PCA colored by the transferred cell types
            DimPlot(combined.atac, reduction = "umaps", label = T, group.by = "predicted_celltype")
            
            
            # set idents
            combined.atac$celltype <- Idents(combined.atac)
            
            
            DefaultAssay(combined.atac) <- "ACTIVITY"
            FeaturePlot(combined.atac, reduction = "umaps", features = c("COL1A1","DCN","PECAM1","TTN"))
            
            
            # Annotations
            combined.atac[["SAMPLE"]] <- ""
            combined.atac@meta.data$SAMPLE[combined.atac@meta.data$STIM %in% c("EFE1","EFE2","EFE3")] <- "EFEs"
            combined.atac@meta.data$SAMPLE[combined.atac@meta.data$STIM %in% c("CTRL1","CTRL2","CTRL3")] <- "CTRLs"
            
            combined.atac[["Cell_Type"]] <- ""
            combined.atac@meta.data$Cell_Type[combined.atac@meta.data$predicted_celltype %in% c("Cap.EC","Art.EC","Ven.EC")] <- "EC_1"
            combined.atac@meta.data$Cell_Type[combined.atac@meta.data$predicted_celltype %in% c("Vent.CM","CM")] <- "CM"
            
            combined.atac@meta.data$Cell_Type[combined.atac@meta.data$predicted_celltype %in% c("EFE1.EC")] <- "EC_2"
            combined.atac@meta.data$Cell_Type[combined.atac@meta.data$predicted_celltype %in% c("EFE3.EC")] <- "EC_3"
            combined.atac@meta.data$Cell_Type[combined.atac@meta.data$predicted_celltype %in% c("Pericyte")] <- "PeriC"
            combined.atac@meta.data$Cell_Type[combined.atac@meta.data$predicted_celltype %in% c("FB")] <- "FB"
            combined.atac@meta.data$Cell_Type[combined.atac@meta.data$predicted_celltype %in% c("SMC")] <- "SMC"
            combined.atac@meta.data$Cell_Type[combined.atac@meta.data$predicted_celltype %in% c("Macrophage")] <- "Macrophage"
            combined.atac@meta.data$Cell_Type[combined.atac@meta.data$predicted_celltype %in% c("EndoC")] <- "EndoC"
            combined.atac@meta.data$Cell_Type[combined.atac@meta.data$predicted_celltype %in% c("T-cell")] <- "T-cell"
            combined.atac@meta.data$Cell_Type[combined.atac@meta.data$predicted_celltype %in% c("Neuro")] <- "Neuro"
            combined.atac@meta.data$Cell_Type[combined.atac@meta.data$predicted_celltype %in% c("Mast-cell")] <- "Mast-cell"
            
            # save
            # saveRDS(combined.atac, file = "/Users/yangyu/Desktop/Harvard_ChenLab/EP/Intermediates/atac_integrated/atac_only/files/combined_seurat_object_110624.rds")
            
            
            
            
