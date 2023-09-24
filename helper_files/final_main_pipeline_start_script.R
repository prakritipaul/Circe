# Script loads necessary and commonly used directories, objects and variables.

####### Download Relevant Objects #######
# Note: Seurat and URD objects are downloaded in /tmp directory.

options(timeout = 600)

# Download Seurat Object unless it already exists.
if (file.exists("/tmp/CNS.Seurat.LTB2.lv.anno.new.Robj")) {
  cat("Seurat Object already exists\n")
  
} else {
  cat("Downloading Seurat object\n")
  # Get zip file for Seurat Object.
  neural_seurat_link <- "https://cocktails-r-objects.s3.amazonaws.com/CNS.Seurat.LTB2.lv.anno.new.Robj.zip"
  # Download in /tmp.
  download.file(neural_seurat_link, destfile = "/tmp/CNS.Seurat.LTB2.lv.anno.new.Robj.zip")
  cat("Unzipping\n")
  unzip("/tmp/CNS.Seurat.LTB2.lv.anno.new.Robj.zip", exdir = "/tmp")
}
cat("Loading\n")
load("/tmp/CNS.Seurat.LTB2.lv.anno.new.Robj")

### Get various data (counts, normalized counts, metadata).
CNS_raw_data <- CNS.Seurat.LTB2.lv@raw.data
CNS_data <- CNS.Seurat.LTB2.lv@data
CNS_metadata <- CNS.Seurat.LTB2.lv@meta.data

### Get all CNS CTs.
# 41
CNS_CTs <- unique(CNS_metadata$newcluster.comb.del.anno)
CNS_CTs <- CNS_CTs[-("NA" %in% CNS_CTs)]
CNS_khids <- rownames(CNS_data)

### Make neural_CT_list = [CT_name : c(cells)]
neural_CT_list <- list()
for (CNS_CT in CNS_CTs) {
  # Get dataframe relevant only to that CT. 
  filtered_neural_df <- CNS_metadata %>% filter(newcluster.comb.del.anno == CNS_CT)
  # Get CT's cells.
  CT_cells <- rownames(filtered_neural_df)
  # Update list.
  neural_CT_list[[CNS_CT]] <- CT_cells
}

# Download URD Object unless it already exists.
if (file.exists("/tmp/CNS.URD.Robj")) {
  cat("URD Object already exists\n")
  
} else {
  # Get zip file for URD object.
  neural_URD_link <- "https://cocktails-r-objects.s3.amazonaws.com/CNS.URD.Robj.zip"
  # Download in /tmp.
  cat("Downloading URD Object\n")
  download.file(neural_URD_link, destfile = "/tmp/CNS.URD.Robj.zip")
  cat("Unzipping\n")
  unzip("/tmp/CNS.URD.Robj.zip", exdir = "/tmp")
}
cat("Loading\n")
load("/tmp/CNS.URD.Robj")

URD_data <- CNS.k300s6w4.tree.built.name@logupx.data
URD_raw_data <- CNS.k300s6w4.tree.built.name@count.data

## These are names of all cells in all segments and pseudotimes of all the cells.
cells_in_segments_list <- CNS.k300s6w4.tree.built.name@tree$cells.in.segment
pseudotimes <- CNS.k300s6w4.tree.built.name@tree$pseudotime

###### kai_tf_df ######
# 639
tf_df <- read.table("./helper_files/Ciona_khid_TF.tsv",
                    sep = "\t",
                    header = TRUE,
                    stringsAsFactors = FALSE)
# 14809
URD_khids <- rownames(URD_raw_data)
# 638
tfs <- tf_df$KH.gene.model
# 595
tfs_in_URD <- intersect(URD_khids, tfs)
# 594 x 2
filtered_tf_df <- tf_df %>% dplyr::select(KH.gene.model, Ghost.Name) %>% filter(KH.gene.model %in% tfs_in_URD)



