
############################ spatial data OPTION 1 ##############################
sObj <- Load10X_Spatial('outs/',
                                  filename = "filtered_feature_bc_matrix.h5",
                                  assay = "Spatial",
                                  slice = "tissue_hires_image.png",
                                  filter.matrix = TRUE,
                                  to.upper = FALSE,
                                  image = NULL)

sObj = SCTransform(sObj, assay = "Spatial", verbose = TRUE, vars.to.regress="nCount_Spatial")
sObj = RunPCA(sObj, npcs=30)
sObj = FindNeighbors(sObj, dims=1:30)
sObj = FindClusters(sObj, resolution=0.2)
sObj = RunUMAP(sObj, dims=1:30, min.dist=0.2)
umap2 = DimPlot(sObj, reduction="umap", label=T)
umap2spatial = SpatialDimPlot(sObj, label=T, pt.size.factor = 1.6, alpha = 1, image.alpha = 0)
umapcomb <- umap2 + umap2spatial

ggsave(filename = "Deconvolution/UMAP_Clusters 0.2_1.pdf", plot = umap2, width=8, height=6)
ggsave(filename = "Deconvolution/UMAP_Clusters 0.2_2.pdf", plot = umap2spatial, width=10, height=10)
ggsave(filename = "Deconvolution/UMAP_Clusters 0.2_combined.pdf", plot = umapcomb, width=20, height=10)
#################################################################################


#########################################################################################################################################################################
# cell cycle scoring Seurat Spatial
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

sObj <- CellCycleScoring(sObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(sObj[[]])

# Visualize the distribution of cell cycle markers across
RidgePlot(sObj, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by phase
sObj <- RunPCA(sObj, features = c(s.genes, g2m.genes))

cell_cycle <- SpatialDimPlot(
  sObj,
  group.by = NULL,
  images = NULL,
  cols = NULL,
  crop = TRUE,
  cells.highlight = NULL,
  #cols.highlight = c("#DE2DFF", "grey10", "grey50"),
  facet.highlight = FALSE,
  label = TRUE,
  label.size = 6,
  label.color = "white",
  repel = TRUE,
  ncol = NULL,
  combine = TRUE,
  pt.size.factor = 1.6,
  alpha = 1, #spot opacity 
  image.alpha = 0, #image opacity
  stroke = 0.25,
  label.box = TRUE,
  interactive = FALSE,
  information = NULL
) & scale_fill_manual(values=c('black', 'blue', 'red'))

ggsave(filename = "Deconvolution/CellCycle.pdf", plot = cell_cycle, width=10, height=10)
#########################################################################################################################################################################


###################### Single cell data option 0 #########################
sObj_original <- sObj
sc_sObj_Original <- readRDS('sObj_aggr.rds')

sObj <- sObj_original
sc_sObj <- sc_sObj_Original
sc_sObj = UpdateSeuratObject(object = sc_sObj)

sc_sObj[[]]
sc_sObj$original <- sc_sObj$seurat_clusters

# check number of cells per subclass
table(sc_sObj$seurat_clusters)

# select 200 cells per subclass, fist set subclass as active.ident
Idents(sc_sObj) <- sc_sObj$seurat_clusters
# sc_sObj <- subset(sc_sObj, cells = WhichCells(sc_sObj, downsample = 200))

# check again number of cells per subclass
table(sc_sObj$seurat_clusters)

sc_sObj@active.assay = "RNA"

markers_sc <- FindAllMarkers(sc_sObj, only.pos = TRUE, max.cells.per.ident = 500, assay = "RNA")
sc_sObj_forLoop <- sc_sObj

# Filter for genes that are also present in the ST data
markers_sc <- markers_sc[markers_sc$gene %in% rownames(sObj), ]

markers_sc_backup <- markers_sc

#"10"
              sObj <- sObj_original
              sc_sObj <- sc_sObj_forLoop
              markers_sc <- markers_sc_backup
              
              # Select top n genes per cluster, select top by first p-value, then absolute diff in pct, then quota of pct.
              markers_sc$pct.diff <- markers_sc$pct.1 - markers_sc$pct.2
              markers_sc$log.pct.diff <- log2((markers_sc$pct.1 * 99 + 1)/(markers_sc$pct.2 * 99 +
                                                                             1))
              markers_sc %>%
                group_by(cluster) %>%
                top_n(-100, p_val) %>%
                top_n(50, pct.diff) %>%
                top_n(20, log.pct.diff) -> top10
              m_feats <- unique(as.character(top10$gene))
              
              eset_SC <- ExpressionSet(assayData = as.matrix(sc_sObj@assays$RNA@counts[m_feats,
              ]), phenoData = AnnotatedDataFrame(sc_sObj@meta.data))
              eset_ST <- ExpressionSet(assayData = as.matrix(sObj@assays$Spatial@counts[m_feats,
              ]), phenoData = AnnotatedDataFrame(sObj@meta.data))
              
              deconvolution_crc <- SCDC::SCDC_prop(bulk.eset = eset_ST, sc.eset = eset_SC, ct.varname = "seurat_clusters",
                                                   ct.sub = as.character(unique(eset_SC$seurat_clusters)))
              head(deconvolution_crc$prop.est.mvw)
              table_deconv <- deconvolution_crc$prop.est.mvw
              write.csv(table_deconv, file='Deconvolution/Deconvolution_matrix_10gene.csv')
              
              sObj@assays[["SCDC"]] <- CreateAssayObject(data = t(deconvolution_crc$prop.est.mvw))
              
              # Seems to be a bug in SeuratData package that the key is not set and any plotting function etc. will throw an error.
              if (length(sObj@assays$SCDC@key) == 0) {
                sObj@assays$SCDC@key = "scdc_"
              }
              
              DefaultAssay(sObj) <- "SCDC"
              
              p1 <- SpatialFeaturePlot(sObj, features = "0", pt.size.factor = 1.6,
                                       crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
              p2 <- SpatialFeaturePlot(sObj, features = "1", pt.size.factor = 1.6,
                                       crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
              p3 <- SpatialFeaturePlot(sObj, features = "2", pt.size.factor = 1.6,
                                       crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
              p4 <- SpatialFeaturePlot(sObj, features = "3", pt.size.factor = 1.6,
                                       crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
              p5 <- SpatialFeaturePlot(sObj, features = "4", pt.size.factor = 1.6,
                                       crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
              p6 <- SpatialFeaturePlot(sObj, features = "5", pt.size.factor = 1.6,
                                       crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
              p7 <- SpatialFeaturePlot(sObj, features = "6", pt.size.factor = 1.6,
                                       crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
              p8 <- SpatialFeaturePlot(sObj, features = "7", pt.size.factor = 1.6,
                                       crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
              p9 <- SpatialFeaturePlot(sObj, features = "8", pt.size.factor = 1.6,
                                       crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
              p10 <- SpatialFeaturePlot(sObj, features = "9", pt.size.factor = 1.6,
                                        crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
              p11 <- SpatialFeaturePlot(sObj, features = "10", pt.size.factor = 1.6,
                                        crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
              p12 <- SpatialFeaturePlot(sObj, features = "11", pt.size.factor = 1.6,
                                        crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
              pcombined = (p1 | p2 | p3) / (p4 | p5 | p6) / (p7 | p8 | p9) / (p10 | p11 | p12)
              #& scale_fill_viridis(limits = c(0.0, 1.0), option = "A")
              ggsave(filename = "Deconvolution/Deconvolution_10gene.pdf", pcombined, width=20, height=20)
              
              sObj <- FindSpatiallyVariableFeatures(sObj, assay = "SCDC", selection.method = "markvariogram",
                                                                      features = rownames(sObj), r.metric = 5, slot = "data")
              top.clusters <- head(SpatiallyVariableFeatures(sObj), 4)
              p_top4 <- SpatialPlot(object = sObj, features = top.clusters, ncol = 2, pt.size.factor = 1.6, crop = TRUE, image.alpha=0) & scale_fill_viridis(option = "A")
              ggsave(filename = "Deconvolution/Deconvolution_10gene_Top4Variable.pdf", p_top4, width=20, height=20)
              pcombined=NULL
              p_top4=NULL
              #eset_SC=NULL
              #eset_ST=NULL
              #deconvolution_crc=NULL
              #m_feats=NULL
              #table_deconv=NULL

#########################################################################################################################################################################
            
            
#########################################################################################################################################################################            
# Cell Chat Spatial
          
            sObj <- readRDS('Deconvolution/sObj.rds')            
            
            sObj[[]]
            sObj <- SetIdent(sObj, value = "seurat_clusters")
            sObj$new <- sObj$old.ident
            sObj$old.ident <- NULL
            
            list1 <- sObj@meta.data[["new"]]
            list1
            list2 <- gsub("5", "C5", list1)
            list3 <- gsub("4", "C4", list2)
            list4 <- gsub("3", "C3", list3)
            list5 <- gsub("2", "C2", list4)
            list6 <- gsub("1", "C1", list5)
            list7 <- gsub("0", "C0", list6)
            list7
            
            sObj <- AddMetaData(sObj, list7, col.name = "old.ident")
            sObj <- SetIdent(sObj, value = "old.ident")

            
            # Prepare input data for CelChat analysis
            data.input = GetAssayData(sObj, slot = "data", assay = "SCT") # normalized data matrix
            meta = data.frame(labels = Idents(sObj), row.names = names(Idents(sObj))) # manually create a dataframe consisting of the cell labels
            unique(meta$labels) # check the cell labels
          
            # Spatial locations of spots from full (NOT high/low) resolution images are required
            spatial.locs = GetTissueCoordinates(sObj, scale = NULL, cols = c("imagerow", "imagecol")) 
            # Scale factors and spot diameters of the full resolution images 
            
            scale.factors = jsonlite::fromJSON(txt = file.path("outs/spatial", 'scalefactors_json.json'))
            scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                                 fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
            )
            # USER can also extract scale factors from a Seurat object, but the `spot` value here is different from the one in Seurat. Thus, USER still needs to get the `spot` value from the json file. 
            
            ###### Applying to different types of spatial imaging data ######
            # `spot.diameter` is dependent on spatial imaging technologies and `spot` is dependent on specific datasets
            
            cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                                       datatype = "spatial", coordinates = spatial.locs, scale.factors = scale.factors)
            
            CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
            
            # use a subset of CellChatDB for cell-cell communication analysis
            # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
            # use all CellChatDB for cell-cell communication analysis
            CellChatDB.use <- CellChatDB # simply use the default CellChatDB
            
            # set the used database in the object
            cellchat@DB <- CellChatDB.use
            
            # subset the expression data of signaling genes for saving computation cost
            cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
            
            cellchat <- identifyOverExpressedGenes(cellchat)
            cellchat <- identifyOverExpressedInteractions(cellchat)
            cellchat <- projectData(cellchat, PPI.human)
            
            #cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
            cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                                          distance.use = TRUE, interaction.length = 200, scale.distance = 0.01)
            
            # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
            #cellchat <- filterCommunication(cellchat, min.cells = 10)
            cellchat <- filterCommunication(cellchat, min.cells = 2)
            
            #df.net <- subsetCommunication(cellchat)
            
            cellchat <- computeCommunProbPathway(cellchat)
            
            cellchat <- aggregateNet(cellchat)
            
            groupSize <- as.numeric(table(cellchat@idents))
            par(mfrow = c(1,2), xpd=TRUE)
            pdf('CellCell Interaction/Number of interactions.pdf', width = 8, height = 10)
            netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
            dev.off()
            
            pdf('CellCell Interaction/Number of interactions, label edge.pdf', width = 8, height = 10)
            netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
            dev.off()
            
            pdf('CellCell Interaction/Interaction weights,strength.pdf', width = 8, height = 10)
            netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
            dev.off()
            
            pathways.show <- c("VCAM")
            # Circle plot
            par(mfrow=c(1,1))
            pdf('CellCell Interaction/VCAM_Circle.pdf', width = 8, height = 10)
            netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
            dev.off()
            
            # Chord diagram
            par(mfrow=c(1,1))
            pdf('CellCell Interaction/VCAM_Chord.pdf', width = 8, height = 10)
            netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
            dev.off()
            
            # Heatmap
            par(mfrow=c(1,1))
            pdf('CellCell Interaction/VCAM_Heatmap.pdf', width = 8, height = 10)
            netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
            dev.off()
            
            # Spatial plot
            par(mfrow=c(1,1))
            pdf('CellCell Interaction/VCAM_Spatial.pdf', width = 8, height = 10)
            netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 3, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)
            dev.off()
            
            # Compute the network centrality scores
            cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
            # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
            par(mfrow=c(1,1))
            netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
            
            # USER can visualize this information on the spatial imaging, e.g., bigger circle indicates larger incoming signaling
            par(mfrow=c(1,1))
            netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, alpha.image = 0.2, vertex.weight = "incoming", vertex.size.max = 3, vertex.label.cex = 3.5)
            
          
            
     
   
            #########################################################################################################################################################################       
            # Cell Chat Spatial with deconvolution clusters
            # run desired gene log cutoff deconvolution and then the script that follows below
            saveRDS(sObj, file = 'Deconvolution/sObj_Deconv.rds')
            
            sObj <- readRDS('Deconvolution/sObj_Deconv.rds')            
            #########################################################################################################################################################################
            ############### things that didn't work
            list1 <- sObj@assays[["SCDC"]]@meta.features[["markvariogram.spatially.variable.rank"]]
            list1
            list1[is.na(list1)] = "0"
            list1
            sObj@assays[["SCDC"]]@meta.features[["markvariogram.spatially.variable.rank"]] <- list1
            
            list1 <- sObj@assays[["SCDC"]]@meta.features[["markvariogram.spatially.variable.rank"]]
            
            
            pbmc_small <- MetaFeature(
              object = sObj,
              features = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"),
              assay = 'SCDC',
              slot = 'data',
              meta.name = "meta.features"
            )
            head(pbmc_small[[]])
            #pbmc_small <- FindSpatiallyVariableFeatures(sObj, assay = \SCDC\, selection.method = \markvariogram\, features = rownames(sObj), r.metric = 5, slot = \data\)
            
            
            
            values <- sObj@assays$SCDC@data["2", ]
            cutoff <- 0.50
            binary_metadata_1 <- ifelse(values > cutoff, "DC2", "No")
            sObj <- AddMetaData(object = sObj, metadata = binary_metadata_1, col.name = "DC2_binary")
            
            values <- sObj@assays$SCDC@data["6", ]
            cutoff <- 0.10
            binary_metadata_2 <- ifelse(values > cutoff, "DC6", "No")
            sObj <- AddMetaData(object = sObj, metadata = binary_metadata, col.name = "DC6_binary")
            
            values <- sObj@assays$SCDC@data["2", ]
            cutoff <- 0.50
            binary_metadata_1 <- ifelse(values > cutoff, 1, 0)
            sObj <- AddMetaData(object = sObj, metadata = binary_metadata, col.name = "DC2_binary")
            
            values <- sObj@assays$SCDC@data["6", ]
            cutoff <- 0.10
            binary_metadata_2 <- ifelse(values > cutoff, 1, 0)
            sObj <- AddMetaData(object = sObj, metadata = binary_metadata, col.name = "DC6_binary")
            
            new_metadata <- binary_metadata_1 | binary_metadata_2
            sObj <- AddMetaData(object = sObj, metadata = new_metadata, col.name = "DC26_binary")
            #########################################################################################################################################################################
            
            
            #################
            #sObj[['merged_string']] <- NULL
            #DefaultAssay(sObj) <- "SCDC"
            #values <- sObj@assays$SCDC@data["4", ]
            #cutoff <- 0.90
            #binary_metadata_1 <- ifelse(values > cutoff, "4", "No")
            
            values_A <- sObj@assays$SCDC@data["4", ]
            values_B <- sObj@assays$SCDC@data["6", ]
            cutoff_A <- 0.80
            cutoff_B <- 0.10
            strr1 <- ifelse(values_A > cutoff_A, "4", "No")
            print(strr1)
            strr2 <- ifelse(values_B > cutoff_B, "6", "No")
            print(strr2)
            merged_string <- character(length(strr1))
            for (i in seq_along(strr1)) {
              if (strr1[i] == "No" && strr2[i] == "No") {
                merged_string[i] <- "No"
              } else if (strr1[i] == "4" && strr2[i] == "No") {
                merged_string[i] <- "4"
              } else if (strr1[i] == "No" && strr2[i] == "6") {
                merged_string[i] <- "6"
              } else {
                merged_string[i] <- paste0(strr1[i], strr2[i])
              }
            }

            print(merged_string)
            #merged_string <- paste0("(", paste(merged_string, collapse = ", "), ")")
            sObj <- AddMetaData(object = sObj, metadata = merged_string, col.name = "merged_string")
            #################
            
            #rm(data.input, meta, string1, string2, values_A, values_B, cutoff_A, cutoff_B)
            sObj <- SetIdent(sObj, value = "merged_string")
            # Prepare input data for CelChat analysis
            data.input = GetAssayData(sObj, slot = "data", assay = "SCT") # normalized data matrix
            meta = data.frame(labels = Idents(sObj), row.names = names(Idents(sObj))) # manually create a dataframe consisting of the cell labels
            unique(meta$labels) # check the cell labels
            #> [1] L2/3 IT Astro   L6b     L5 IT   L6 IT   L6 CT   L4      Oligo  
            #> Levels: Astro L2/3 IT L4 L5 IT L6 IT L6 CT L6b Oligo
            
            # load spatial imaging information
            # Spatial locations of spots from full (NOT high/low) resolution images are required
            spatial.locs = GetTissueCoordinates(sObj, scale = NULL, cols = c("imagerow", "imagecol")) 
            # Scale factors and spot diameters of the full resolution images 
            
            scale.factors = jsonlite::fromJSON(txt = file.path("outs/spatial", 'scalefactors_json.json'))
            scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                                 fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
            )
            # USER can also extract scale factors from a Seurat object, but the `spot` value here is different from the one in Seurat. Thus, USER still needs to get the `spot` value from the json file. 
            
            ###### Applying to different types of spatial imaging data ######
            # `spot.diameter` is dependent on spatial imaging technologies and `spot` is dependent on specific datasets
            
            cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                                       datatype = "spatial", coordinates = spatial.locs, scale.factors = scale.factors)
            
            CellChatDB <- CellChatDB.human # use CellChatDB.human if running on human data
            
            # use a subset of CellChatDB for cell-cell communication analysis
            # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
            # use all CellChatDB for cell-cell communication analysis
            CellChatDB.use <- CellChatDB # simply use the default CellChatDB
            
            # set the used database in the object
            cellchat@DB <- CellChatDB.use
            
            # subset the expression data of signaling genes for saving computation cost
            cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
            
            cellchat <- identifyOverExpressedGenes(cellchat)
            cellchat <- identifyOverExpressedInteractions(cellchat)
            cellchat <- projectData(cellchat, PPI.human)
            
            #cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
            cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                                          distance.use = TRUE, interaction.length = 200, scale.distance = 0.01)
            
            # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
            #cellchat <- filterCommunication(cellchat, min.cells = 10)
            cellchat <- filterCommunication(cellchat, min.cells = 2)
            
            #df.net <- subsetCommunication(cellchat)
            
            cellchat <- computeCommunProbPathway(cellchat)
            
            cellchat <- aggregateNet(cellchat)
            
            groupSize <- as.numeric(table(cellchat@idents))
            par(mfrow = c(1,2), xpd=TRUE)
            pdf('CellCell Interaction/Number of interactions_with deconv clusters trial.pdf', width = 8, height = 10)
            netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
            dev.off()
            
            pdf('CellCell Interaction/Number of interactions, label edge_with deconv clusters trial.pdf', width = 8, height = 10)
            netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
            dev.off()
            
            pdf('CellCell Interaction/Interaction weights,strength_with deconv clusters trial.pdf', width = 8, height = 10)
            netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
            dev.off()
            
            pathways.show <- c("VCAM")
            # Circle plot
            par(mfrow=c(1,1))
            pdf('CellCell Interaction/VCAM_Circle_with deconv clusters trial.pdf', width = 8, height = 10)
            netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
            dev.off()
            
            # Chord diagram
            par(mfrow=c(1,1))
            pdf('CellCell Interaction/VCAM_Chord_with deconv clusters trial.pdf', width = 8, height = 10)
            netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
            dev.off()
            
            # Heatmap
            par(mfrow=c(1,1))
            pdf('CellCell Interaction/VCAM_Heatmap_with deconv clusters trial.pdf', width = 8, height = 10)
            netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
            dev.off()
            
            # Spatial plot
            par(mfrow=c(1,1))
            pdf('CellCell Interaction/VCAM_Spatial_with deconv clusters trial.pdf', width = 8, height = 10)
            netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 3, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)
            dev.off()
            
            # Compute the network centrality scores
            cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
            # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
            par(mfrow=c(1,1))
            netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
            
            # USER can visualize this information on the spatial imaging, e.g., bigger circle indicates larger incoming signaling
            par(mfrow=c(1,1))
            netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, alpha.image = 0.2, vertex.weight = "incoming", vertex.size.max = 3, vertex.label.cex = 3.5)
            
#########################################################################################################################################################################                        
            
            
            
            
