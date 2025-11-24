## File with the function definition used in all the scripts for the single cell analysis 

# This function is created to load data for one patient 
# Input : patient id, data folder path (path for the patient folder ex: P1), days, and name file (from 10X Genomics)
# Output : seurat object list loaded 
load_data_10XG <- function(patient_id, jours, base_dir, name_file) {
  # Construct the path
  paths <- file.path(base_dir, patient_id, paste0("run_count_J", jours), "outs", name_file)
  # Data load
  data_list <- lapply(paths, Read10X_h5)
  names(data_list) <- paste0("CLL", "_D", jours, ".data")
  # Seurat Object creation RNA and ADT assays
  seurat_list <- lapply(seq_along(data_list), function(i) {
    # RNA and ADT extraction
    rna_counts <- data_list[[i]]$`Gene Expression`
    adt_counts <- data_list[[i]]$`Antibody Capture`
    # Align cell barcodes across RNA and ADT
    common_barcodes <- intersect(colnames(rna_counts), colnames(adt_counts))
    rna_counts <- rna_counts[, common_barcodes]
    adt_counts <- adt_counts[, common_barcodes]
    adt_counts <- adt_counts[, match(colnames(rna_counts), colnames(adt_counts))]
    # Seurat Object with RNA
    seurat_obj <- CreateSeuratObject(
      counts = rna_counts,
      project = paste0(patient_id, "_D", jours[i])
    )
    
    # Add ADT in 2nd assay
    seurat_obj[["ADT"]] <- CreateAssayObject(counts = adt_counts)
    return(seurat_obj)
  })
  # Name the data
  names(seurat_list) <- paste0("CLL", "_D", jours)
  return(seurat_list)
}

# This function is created to load data for one patient 
# Input : patient id, data folder path (path for the patient folder ex: P1), days, and name file (from CellBender)
# Output : seurat object list loaded 
load_data_cellbender <- function(patient_id, jours, base_dir, name_file) {
  # Construct the path
  paths <- file.path(base_dir, patient_id, paste0("run_count_J", jours), "outs", name_file)
  # Data load
  data_list <- lapply(paths, Read_CellBender_h5_Mat)
  names(data_list) <- paste0("CLL", "_D", jours, ".data")
  # Seurat Object creation RNA and ADT assays
  seurat_list <- lapply(seq_along(data_list), function(i) {
    # RNA and ADT extraction
    rna_counts <- data_list[[i]]$`Gene Expression`
    adt_counts <- data_list[[i]]$`Antibody Capture`
    # Align cell barcodes across RNA and ADT
    common_barcodes <- intersect(colnames(rna_counts), colnames(adt_counts))
    rna_counts <- rna_counts[, common_barcodes]
    adt_counts <- adt_counts[, common_barcodes]
    adt_counts <- adt_counts[, match(colnames(rna_counts), colnames(adt_counts))]
    # Seurat Object with RNA
    seurat_obj <- CreateSeuratObject(
    counts = rna_counts,
    project = paste0(patient_id, "_D", jours[i])
    )
  
    # Add ADT in 2nd assay
    seurat_obj[["ADT"]] <- CreateAssayObject(counts = adt_counts)
    return(seurat_obj)
  })
# Name the data
  names(seurat_list) <- paste0("CLL", "_D", jours)
  return(seurat_list)
}


# This function takes a single-cell RNA-seq sample (in Seurat object format) as input and performs all the necessary preprocessing steps required 
# by DoubletFinder, including normalization, identification of variable features, scaling, and PCA.
# It then determines the optimal parameters for DoubletFinder (PCs, pN, pK, and nExp), runs the tool, and returns a data frame containing the IDs 
# of all cells classified as doublets in the sample.
complete_DoubletFinder<-function(sample_data, multiplet_rate = NULL){

# pre-processing
sample <- NormalizeData(sample_data)
sample <- FindVariableFeatures(sample)
sample <- ScaleData(sample)
sample <- RunPCA(sample, nfeatures.print = 10)

# Find significant PCs
stdv <- sample[["pca"]]@stdev
percent_stdv <- (stdv/sum(stdv)) * 100
cumulative <- cumsum(percent_stdv)
co1 <- which(cumulative > 90 & percent_stdv < 5)[1] 
co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] - 
                     percent_stdv[2:length(percent_stdv)]) > 0.1), 
            decreasing = T)[1] + 1
min_pc <- min(co1, co2)

# Finish pre-processing with min_pc
sample <- RunUMAP(sample, dims = 1:min_pc)
sample <- FindNeighbors(object = sample, dims = 1:min_pc)              
sample <- FindClusters(object = sample, resolution = 0.1)

# pK identification
sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)   
sweep_stats <- summarizeSweep(sweep_list)
bcmvn <- find.pK(sweep_stats)
optimal.pk <- bcmvn %>% 
  dplyr::filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK)
optimal.pk <- as.numeric(as.character(optimal.pk[[1]])) 

# Homotypic doublet proportion estimate
annotations <- sample@meta.data$seurat_clusters # use the clusters as the user-defined cell types
homotypic.prop <- modelHomotypic(annotations) # get proportions of homotypic doublets

# Get the multiplet rate if not provided
if(is.null(multiplet_rate)){
  print('multiplet_rate not provided....... estimating multiplet rate from cells in dataset')
  # 10X multiplet rates table
  multiplet_rates_10x <- data.frame('Multiplet_rate'= c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
                                    'Loaded_cells' = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
                                    'Recovered_cells' = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000))
  
  #print(multiplet_rates_10x)
  
  multiplet_rate <- multiplet_rates_10x %>% dplyr::filter(Recovered_cells < nrow(sample_data@meta.data)) %>% 
    dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold depending on your number of samples
    dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # get the expected multiplet rate for that number of recovered cells
  
  print(paste('Setting multiplet rate to', multiplet_rate))
}
nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # multiply by number of cells to get the number of expected multiplets
nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of doublets

# run DoubletFinder
sample <- doubletFinder(seu = sample, 
                        PCs = 1:min_pc, 
                        pK = optimal.pk, # the neighborhood size used to compute the number of artificial nearest neighbours
                        nExp = nExp.poi.adj) # number of expected real doublets
# change name of metadata column with Singlet/Doublet information
colnames(sample@meta.data)[grepl('DF.classifications.*', colnames(sample@meta.data))] <- "doublet_finder"
return(sample)
}


# Non-coding genes removal 
# This function accepts a Seurat object at any stage (pre- or post-processing) 
# and returns a version of the object with all non-coding genes filtered out. 
removal_noncoding_gene <- function(seuratobj){
  
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # load the database Ensembl
gene_list <-Features(seuratobj[["RNA"]]) # gene listing of the seurat object
# Query Ensembl for gene biotype
gene_info <- getBM(attributes = c("external_gene_name", "gene_biotype"),
                   filters = "external_gene_name",
                   values = gene_list,
                   mart = mart)
# Filter to keep only protein-coding genes
coding_genes <- gene_info$external_gene_name[gene_info$gene_biotype == "protein_coding"] # list the coding genes
# RNA count matrix
counts_matrix <- GetAssayData(seuratobj, assay = "RNA", layer = "counts")
# Keep only the coding genes from the  RNA counts
counts_matrix_filtered <- counts_matrix[rownames(counts_matrix) %in% coding_genes, ]
# Create the RNA assay with only the coding genes
seuratobj[["RNA"]] <- CreateAssay5Object(counts = counts_matrix_filtered)
return(seuratobj) 
}

# Run graph reduction as UMAP
run_umap <- function(seuratobj, Resolution){
  seuratobj <- NormalizeData(seuratobj)
  seuratobj <- FindVariableFeatures(seuratobj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(seuratobj)
  seuratobj <- ScaleData(seuratobj, features = all.genes)
  seuratobj <- RunPCA(seuratobj, features = VariableFeatures(object = seuratobj))
  seuratobj <- FindNeighbors(seuratobj, dims = 1:10)
  seuratobj <- FindClusters(seuratobj, resolution=Resolution)
  seuratobj <- RunUMAP(seuratobj, dims = 1:10)
  return(seuratobj)
}

# Run UMAP on ADT
run_adt_umap <- function(seuratobj) {
  DefaultAssay(seuratobj) <- "ADT"
  seuratobj <- NormalizeData(seuratobj, normalization.method = "CLR", margin = 2)
  seuratobj <- ScaleData(seuratobj)
  all_feats <- rownames(seuratobj[["ADT"]])
  VariableFeatures(seuratobj) <- all_feats
  seuratobj <- RunPCA(seuratobj, reduction.name = "apca", verbose = FALSE)
  seuratobj <- RunUMAP(seuratobj, reduction = "apca", dims = 1:18,
                       reduction.name = "adt.umap", reduction.key = "ADTUMAP_")
  DefaultAssay(seuratobj) <- "RNA"
  return(seuratobj)
}

# Run wnnUMAP graph 
run_wnnumap <- function(seuratobj, resolution, normalization_type) {
  # Multimodal neighbors & WNN UMAP
  DefaultAssay(seuratobj) <- normalization_type  # "RNA" or "SCT"
  seuratobj <- FindMultiModalNeighbors(
    seuratobj,
    reduction.list = list("pca", "apca"),
    dims.list = list(1:30, 1:18),
    modality.weight.name = "RNA.weight"
  )
  seuratobj <- RunUMAP(seuratobj, nn.name = "weighted.nn",
                       reduction.name = "wnn.umap", reduction.key = "WNNUMAP_")
  seuratobj <- FindClusters(seuratobj, graph.name = "wsnn",
                            algorithm = 3, resolution = resolution, verbose = FALSE)
  return(seuratobj)
}

library(Seurat)


# Function: harmonize_label
harmonize_label <- function(label) {
  label <- as.character(label)
  if (is.na(label) || label == "") return(NA)
  l <- tolower(gsub("[ _-]", "", label))
  
  if (l %in% c("cd4t","cd4tcell","cd4ctl","cd4tcm","cd4tem","cd4naive","cd4proliferating")) return("CD4 T")
  if (l %in% c("cd8t","cd8tcell","cd8tcm","cd8tem","cd8naive","cd8proliferating")) return("CD8 T")
  if (l %in% c("tcell","tcells","t","dnt","gdt","treg","mait")) return("T")
  if (l %in% c("bnaive","bmemory","bintermediate")) return(label)
  if (l %in% c("bcell","bcells","b")) return("B")
  if (grepl("prebcellcd34", l) || grepl("probcellcd34", l)) return("B")
  if (grepl("monocyte", l) || grepl("cd14mono", l) || grepl("cd16mono", l)) return("Monocyte")
  if (grepl("macrophage", l)) return("Macrophage")
  if (grepl("neutrophil", l)) return("Neutrophil")
  if (grepl("^nk", l)) return("NK")
  if (grepl("^dc", l) || grepl("dendriticcell", l) || grepl("asdc", l) || 
      grepl("cdc1", l) || grepl("cdc2", l) || grepl("pdc", l)) return("DC")
  if (grepl("platelet", l)) return("Platelet")
  if (grepl("eryth", l)) return("Erythroblast")
  if (grepl("^myelocyte$", l)) return("Myelocyte")
  if (grepl("^pro-myelo", l)) return("Pro-Myelocyte")
  if (grepl("hspc", l)) return("HSPC")
  if (grepl("ilc", l)) return("ILC")
  if (grepl("plasma", l) || grepl("plasmablast", l)) return("Plasmablast")
  
  return(label)
}

final_consensus <- function(singleR, manual, azimuth) {
  lbls <- list(singleR, manual, azimuth)
  labs_h <- sapply(lbls, harmonize_label)
  not_na <- which(!is.na(labs_h))
  count_labels <- table(labs_h[not_na])
  
  # Define cell families
  tb_group <- c("T", "CD4 T", "CD8 T")
  b_group <- c("B", "B intermediate", "B memory", "B naive")
  
  # If only one valid annotation
  if(length(not_na) == 1) return(labs_h[not_na])
  
  # If two valid annotations
  if(length(not_na) == 2) {
    labs_subset <- labs_h[not_na]
    # If they are identical
    if(labs_subset[1] == labs_subset[2]) return(labs_subset[1])
    # T-family: keep the most specific
    if(all(labs_subset %in% tb_group)) {
      finer_t <- labs_subset[labs_subset %in% c("CD4 T", "CD8 T")]
      if(length(finer_t)) return(finer_t[1]) else return("T")
    }
    # B-family: keep the most specific
    if(all(labs_subset %in% b_group)) {
      finer_b <- labs_subset[labs_subset %in% c("B intermediate", "B memory", "B naive")]
      if(length(finer_b)) return(finer_b[1]) else return("B")
    }
    # No agreement: return manual
    return(harmonize_label(manual))
  }
  
  # If three valid annotations
  if(length(not_na) == 3) {
    labs_all <- labs_h[not_na]
    
    # If exact majority exists (2 or 3 identical)
    if(any(count_labels >= 2)) {
      maj_label <- names(count_labels)[which.max(count_labels)]
      # T-family: keep the most specific
      if(maj_label == "T" && any(labs_all %in% c("CD4 T", "CD8 T"))) {
        finer_t <- labs_all[labs_all %in% c("CD4 T", "CD8 T")]
        return(finer_t[1])
      }
      # B-family: keep the most specific
      if(maj_label == "B" && any(labs_all %in% c("B intermediate", "B memory", "B naive"))) {
        finer_b <- labs_all[labs_all %in% c("B intermediate", "B memory", "B naive")]
        return(finer_b[1])
      }
      return(maj_label)
    }
    
    # If all different, check for majority families
    # Count how many annotations belong to each family
    t_count <- sum(labs_all %in% tb_group)
    b_count <- sum(labs_all %in% b_group)
    
    # If at least 2 annotations are in T-family
    if(t_count >= 2) {
      t_labels <- labs_all[labs_all %in% tb_group]
      finer_t <- t_labels[t_labels %in% c("CD4 T", "CD8 T")]
      if(length(finer_t)) return(finer_t[1]) else return("T")
    }
    
    # If at least 2 annotations are in B-family
    if(b_count >= 2) {
      b_labels <- labs_all[labs_all %in% b_group]
      finer_b <- b_labels[b_labels %in% c("B intermediate", "B memory", "B naive")]
      if(length(finer_b)) return(finer_b[1]) else return("B")
    }
    
    # Otherwise return harmonized manual
    return(harmonize_label(manual))
  }
  
  # If no valid annotation
  return(NA)
}

# Function convert a seurat object v5 like a seurat object v4
ConvertToV4 <- function(seurat_v5, assay = "RNA") {
  # Check input
  if (!inherits(seurat_v5, "Seurat")) {
    stop("Input object is not a Seurat object.")
  }
  
  # Extract data layers from Seurat v5 Assay5
  counts <- LayerData(seurat_v5, assay = assay, layer = "counts")
  data   <- LayerData(seurat_v5, assay = assay, layer = "data")
  scale  <- LayerData(seurat_v5, assay = assay, layer = "scale.data")
  
  # Create a new Seurat object (v4-style) using counts
  seurat_v4 <- Seurat::CreateSeuratObject(
    counts = counts,
    meta.data = seurat_v5@meta.data,
    assay = assay
  )
  
  # Add the other slots to mimic Seurat v4 structure
  seurat_v4[[assay]]@data <- data
  seurat_v4[[assay]]@scale.data <- scale
  
  message("Conversion complete: object is now compatible with v4-based functions")
  return(seurat_v4)
}





