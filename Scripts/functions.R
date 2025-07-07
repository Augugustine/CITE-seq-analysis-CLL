## File with the function definition used in all the scripts for the single cell analysis 

# This function is created to load filtered data for one patient 
# Input : patient number, data path, days 
# Output : seurat object list loaded 
load_data_filtered <- function(patient_id, jours, base_dir = "~/Documents/Singlecell-R/") {
  # Construct the path
  paths <- file.path(base_dir, patient_id, paste0("/run_count_J", jours,"/"), "/outs", "/filtered_feature_bc_matrix.h5")
  # Data load
  data_list <- lapply(paths, Read10X_h5)
  names(data_list) <- paste0("CLL", "_D", jours, "_filtered.data")
  # Seurat Object creation ('Gene Expression')
  seurat_list <- lapply(seq_along(data_list), function(i) {
    CreateSeuratObject(
      counts = data_list[[i]]$`Gene Expression`,
      project = paste0(patient_id, "_D", jours[i], "_filtered"),
      min.cells = 2,
      min.features = 100
    )
  })
  
  # Name the data
  names(seurat_list) <- paste0("CLL", "_D", jours, "_filtered")
  return(seurat_list)
}

# This function is created to load raw data for one patient 
# Input : patient number, data path, days 
# Output : seurat object list loaded 
load_data_raw <- function(patient_id, jours, base_dir = "~/Documents/Singlecell-R/") {
  # Construct the path
  paths <- file.path(base_dir, patient_id, paste0("/run_count_J", jours,"/"), "/outs", "/raw_feature_bc_matrix.h5")
  # Data load
  data_list <- lapply(paths, Read10X_h5)
  names(data_list) <- paste0("CLL_", patient_id, "_D", jours, "_raw.data")
  # Seurat Object creation ('Gene Expression')
  seurat_list <- lapply(seq_along(data_list), function(i) {
    CreateSeuratObject(
      counts = data_list[[i]]$`Gene Expression`,
      project = paste0(patient_id, "_D", jours[i], "_raw"),
      min.cells = 2,
      min.features = 100
    )
  })
  
  # Name the data
  names(seurat_list) <- paste0("CLL_", patient_id, "_D", jours, "_raw")
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

double_finder_res <- sample@meta.data['doublet_finder'] # get the metadata column with singlet, doublet info
double_finder_res <- rownames_to_column(double_finder_res, "row_names") # add the cell IDs as new column to be able to merge correctly
return(double_finder_res)
}



# Non-coding genes removal 
# This function accepts a Seurat object at any stage (pre- or post-processing) 
# and returns a version of the object with all non-coding genes filtered out. 
removal_noncoding_gene <- function(seuratobj){
  
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # load the database Ensembl
gene_list <-Features(seuratobj) # gene listing of the seurat object
# Query Ensembl for gene biotype
gene_info <- getBM(attributes = c("external_gene_name", "gene_biotype"),
                   filters = "external_gene_name",
                   values = gene_list,
                   mart = mart)
# Filter to keep only protein-coding genes
coding_genes <- gene_info$external_gene_name[gene_info$gene_biotype == "protein_coding"] # list the coding genes
seuratobj<-subset(seuratobj, features = coding_genes) # keep only the coding gene of the list, in the seurat object
return(seuratobj) 
}

