# Key TFs and Marker Genes.
msxb <- "KH2012:KH.C2.957"
ng <- "KH2012:KH.C6.129"
foxg <- "KH2012:KH.C8.774"
foxc <- "KH2012:KH.L57.25"
prop <- "KH2012:KH.C11.543"
pou4 <- "KH2012:KH.C2.42"
hnf6 <- "KH2012:KH.C6.185"
coe <- "KH2012:KH.L24.10"
pole3 <- "KH2012:KH.L41.26"
alx_arx <- "KH2012:KH.L119.26"
foxHa <- "KH2012:KH.C9.717"
Orphan_bHLH_1 <- "KH2012:KH.C7.269"
islet <- "KH2012:KH.L152.2"
lhx1 <- "KH2012:KH.L107.7"
lhx3 <- "KH2012:KH.S215.4"
orphan_fox2 <- "KH2012:KH.L150.2"
orphan_fox1 <- "KH2012:KH.C14.520"
sp8 <- "KH2012:KH.C13.22"
otx <- "KH2012:KH.C4.84"
mitf <- "KH2012:KH.C10.106"
dmrt <- "KH2012:KH.S544.3"
bsh <- "KH2012:KH.C8.652"
meis <- "KH2012:KH.C10.174"
ptf1a <- "KH2012:KH.L116.39"
ttf1 <- "KH2012:KH.C10.338"
nk4 <- "KH2012:KH.C8.482"
scalloped <- "KH2012:KH.C14.426"
six3_6 <- "KH2012:KH.C10.367"
dmbx <- "KH2012:KH.C1.1212"
xBPd <- "KH2012:KH.C8.877"
# PSC cocktail
dril1_2 <- "KH2012:KH.S406.9"
# dmbx+ GLRA-MG Cocktail
lhx1_5 <- "KH2012:KH.L107.7"
hmg1_2 <- "KH2012:KH.C14.178"

## General Comments ##
# (1) khid is the id for genes in the Ciona genome.
# (2) munged/not munged refers to whether names of TFs are khids (not munged)
#     or Ortholog names (munged)

write_csv_helper <- function(character_vector, filename) {
  # Takes a character vector (e.g. BTN_precursor_pseudostate_DEGs) and writes out
  # a csv file in which each entry of the character vector occupies its own line.
  #
  # Notes: There are no row or column names; There are no commas delimiting the lines.
  # Use case: These contain Candidate Cocktail names and DEG lists.
  write.table(character_vector, filename, row.names = FALSE, col.names = FALSE, 
              sep = ",", quote = FALSE)
}

convert_khids_to_ortholog <- function(DEG_TF_df) {
  # Takes in a DEG_TF_df with rows = khid, cols = logFC, logCPM, F, 
  # pval, and adj_pval, and returns a named_DEG_TF_df, 
  # with cols = KH.gene.model (khid) and Ghost.Name (name of ortholog).
  #
  # Args:
  #   DEG_TF_df: data frame as explained above. It tends to be an edgeR output.
  #
  # Returns:
  #   named_DEG_TF_df: data frame as expalined above.
  # 
  # Example Usage:
  #   named_terminal_BTN_v_siblings_TF_DEGs_df <- convert_khids_to_ortholog(terminal_BTN_v_siblings_TF_DEGs_df)
  #   where terminal_BTN_v_siblings_TF_DEGs_df <- terminal_BTN_v_siblings_TF_edgeR_outputs[[1]].
  #
  # Use case: This makes it easy to know the Ortholog names of the DEG TFs.
  # 
  DEG_TF_df_2 <- rownames_to_column(DEG_TF_df, "KH.gene.model")
  named_DEG_TF_df <- DEG_TF_df_2 %>% inner_join(filtered_tf_df, by = "KH.gene.model") %>%  
    dplyr::select(KH.gene.model, Ghost.Name)
  return(named_DEG_TF_df)
}

get_pseudotimed_cells <- function(cells_in_segments_list, pseudotimes, segment_ID, num_cells, ascending=TRUE) {
  # Returns list of ordered pseudotimed cells from a segment either in ascending
  # or descending order.
  # 
  # Args:
  #   cells_in_segments_list: list with keys = segment_IDs and values = 
  #       cell barcodes of cells that belong to the given segment.
  #   pseudotimes: character vector that contains pseudotimes of all cells from 
  #       all stages (stored in in CNS URD object).
  #   segment_ID: int. ID/number of segment on Virtual Lineage Tree.
  #   num_cells: int; number of pseudotimed cells wanted from segment.
  #   ascending: boolean.
  #     Use ascending for precursor and its sibling(s)
  #       The "first" cells of segment will be chosen.
  #     Use descending for mother and its sibling(s) 
  #       The "last"cells of segment will be chosen.
  #
  # Returns:
  #   ordered_cells: character vector as described above.
  #
  cells <- cells_in_segments_list[[segment_ID]]
  # Named vector of cells and their pseudotimes.
  cells_and_pseudotimes <- pseudotimes[cells]
  if (ascending) {
    # This just gives you the order (int vector).
    order <- order(cells_and_pseudotimes)
  }
  else {
    order <- order(cells_and_pseudotimes, decreasing = TRUE)
  }
  ordered_cells <- names(cells_and_pseudotimes[order])
  final_ordered_cells <- ordered_cells[1:num_cells]
  return(final_ordered_cells)
}

add_segment_state <- function(Cells_list, Cells_list_string) {
  # Adds a segment state (which is any condition that describes a subset of cells
  # present along the segment) as a group column in the group.id data frame in URD Object.
  # For instance, the segment state "BTN_precursor_pseudostate_cells" are all the
  # precursor pseudostate cells present on the BTN precursor segment.
  #
  # Args:
  #   Cells_list: list with key = string that describes cells and
  #     value = cell barcodes of cells.
  #   Cells_list_string: a string that has the same name as Cells_list.
  #
  # Routine:
  #   Iterates through all keys of Cells_list and gets their cell_barcodes.
  #   Each cell_barcode is assigned to its segment state. Both "cell_barcode"
  #   and "segment_state" columns are formed.
  #
  # Returns:
  #   URD_object with modified group ids.
  #   Its group.id dataframe will now have "cell_barcode" and "segment_state" columns.
  #
  # Use case: We use the segment state to visualize customized groups of cells 
  #   along the pseudotimed tree.
  
  neural_URD_location <- paste0(cocktail_data_dir, "CNS.URD.Robj")
  load(neural_URD_location)

  modified_group_ids <- CNS.k300s6w4.tree.built.name@group.ids
  # Make a column with names of cells (barcodes)
  modified_group_ids$cell_barcode <- rownames(modified_group_ids)

  Cells_list_keys <- names(Cells_list)
  # Iterates through all keys of Cells_list and gets their cell_barcodes.
  commands <- parse_exprs(paste0("cell_barcode %in%", Cells_list_string, "$",
                                 Cells_list_keys, "~", "'", Cells_list_keys, "'"))
  # Create segment state column with above commands. 
  modified_group_ids_2 <- modified_group_ids %>% mutate(segment_state = case_when(!!!commands))

  CNS.k300s6w4.tree.built.name@group.ids <- modified_group_ids_2
  return(CNS.k300s6w4.tree.built.name)
}

make_PC_pseudotime_df <- function(pseudotimed_cells, pca_coordinates, num_pcs, pseudotimes) {
  # Creates a PC_pseudotime_df that has PC coordinates and pseudotimes of given 
  # pseudotimed_cells.
  #
  # Args:
  #   psuedotimed_cells: character vector; self-explanatory.
  #   pca_coordinates: data frame with rows = cell barcodes and columns = PCs.
  #     Entries are PC values of cell barcodes. 
  #     It is an Embeddings Seurat object. 
  #     e.g. BTN_PCA_COORDINATES <- Embeddings(BTN_seurat, reduction = "pca") 
  #   num_pcs: int; number of PCs chosen to index pca_coordinates.
  #   pseudotimes: character vector that contains pseudotimes of all cells from 
  #       all stages (stored in in CNS URD object).
  #
  # Returns:
  #   PC_pseudotime_df: data frame with rows = cell barcodes, 
  #     columns = cell_barcode, pseudotime, PC_1... PC_(num_pcs)
  #
  # Use case:
  #   We use the PC_pseudotime_dfs in "get_pseudostates.py" to find Precursor,
  #   Precursor Sibling(s), and Mother pseudostates.
  PC_pseudotime_df <- as.data.frame(pca_coordinates[pseudotimed_cells, 1:num_pcs])
  PC_pseudotime_df$cell_barcode <- rownames(PC_pseudotime_df)
  PC_pseudotime_df$pseudotime <- pseudotimes[pseudotimed_cells]
  PC_pseudotime_df <- PC_pseudotime_df %>% dplyr::select(cell_barcode, pseudotime, everything())
  return(PC_pseudotime_df)
}

get_mother_sibling_pseudostates <- function(mother_pseudostate_cells, pseudotimed_mother_sibling_cells) {
  # Finds which pseudotimed_mother_sibling_cells fall within the upper and lower
  #   pseudotime bounds of mother_pseudostate_cells.
  #
  # Args:
  #   mother_pseudostate_cells: as previously described.
  #   pseudotimed_mother_sibling_cells: as previously described.
  #
  # Routine:
  #   Finds the minimum and maximum pseudotimes in mother_pseudostate_cells and then
  #   the indices of cells in pseudotimed_mother_sibling_cells whose pseudotimes fall
  #   within that range.
  #
  # Notes: pseudotimes are ordered in nonincreasing order.
  # 
  # Returns:
  #   mother_sibling_indices: a list with:
  #       [["first_index"]] <- first index in mother sibling's pseudotimes
  #         that is less than the maximum value of mother's pseudotime.
  #         (because the mother sibling's pseudotimes are ordered in
  #          non-decreasing order)
  #          This corresponds to the larger pseudotime.
  #       [["min_pseudotime"]] <- pseudotime value in mother sibling's pseudotimes
  #         at first_index.
  #       [["last_index"]] <- last index in mother sibling's pseudotimes
  #         that is greater than the minimum value of mother's pseudotime.
  #          This renders the smaller pseudotime.
  #       [["max_pseudotime"]] <- pseudotime value in mother sibling's pseudotimes
  #         at last_index.
  #
  # Example Use Case:
  #   BTN_mother_sibling_63_output <- get_mother_sibling_pseudostates(BTN_mother_pseudostate_cells, pseudotimed_BTN_mother_sibling_63_segment_cells)
  #   mother_sibling_63_first_index <- BTN_mother_sibling_63_output$first_index
  #   mother_sibling_63_last_index <- BTN_mother_sibling_63_output$last_index
  #
  #   (Do also for BTN_mother_sibling_73_cells)
  #
  #   BTN_mother_sibling_pseudostate_cells <- c(pseudotimed_BTN_mother_sibling_63_segment_cells[mother_sibling_63_first_index:mother_sibling_63_last_index],
  #                                             pseudotimed_BTN_mother_sibling_73_segment_cells[mother_sibling_73_first_index:mother_sibling_73_last_index])
  mother_pseudostate_pseudotimes <- pseudotimes[mother_pseudostate_cells]
  mother_sibling_pseudotimes <- pseudotimes[pseudotimed_mother_sibling_cells]
  
  # What are the pseudotime bounds of your estimated mother cells?
  mother_pseudostate_pseudotimes_min <- min(mother_pseudostate_pseudotimes)
  mother_pseudostate_pseudotimes_max <- max(mother_pseudostate_pseudotimes)
  
  # First index of sibling's pseudotime that is less than mother's max pseudotime.
  first_mother_sibling_index <- min(which(mother_sibling_pseudotimes < mother_pseudostate_pseudotimes_max))
  # Last index of sibling's pseudotime that is greater than mother's min pseudotime.
  last_mother_sibling_index <- max(which(mother_sibling_pseudotimes > mother_pseudostate_pseudotimes_min))
  # Cases in which all sibling pseudotimes are greater/less the max/min of mother's pseudotimes.
  if (first_mother_sibling_index == Inf || last_mother_sibling_index == -Inf) {
    first_mother_sibling_index <- 0
    last_mother_sibling_index <- 0
  }
  
  # Now index the sibling's pseudotimes.
  first_mother_sibling_time <- mother_sibling_pseudotimes[first_mother_sibling_index]
  last_mother_sibling_time <- mother_sibling_pseudotimes[last_mother_sibling_index]
  
  # Output.
  mother_sibling_indices <- list()
  mother_sibling_indices[["first_index"]] <- first_mother_sibling_index
  mother_sibling_indices[["min_pseudotime"]] <- first_mother_sibling_time
  mother_sibling_indices[["last_index"]] <- last_mother_sibling_index
  mother_sibling_indices[["max_pseudotime"]] <- last_mother_sibling_time  
  
  return(mother_sibling_indices)
}

edgeR_pipeline <- function(pair_mat, rep_1, rep_2) {
  # Computes differentially expressed genes (DEGs) between two samples.
  # These DEGs are for rep_1 sample with respect to rep_2 sample.
  #
  # Args:
  #   pair_mat: matrix with raw gene expression values of the two samples. 
  #     rows = genes and cols = (cells_1 + cells_2)
  #     e.g. BTN_precursor_v_sibling_pseudostate_mat
  #   rep_1 and rep_2: numbers in each sample.
  #     e.g. rep_1 = BTN precursor and rep_2 = BTN precursor sibling    
  #
  # Returns:
  #   edgeR_res: data frame with logFC, logCPM, F, PValue 
  #     statistics for the DEGs.
  group = factor(c(rep(1, rep_1), rep(2, rep_2)))
  y = DGEList(counts=pair_mat, group=group)
  y = calcNormFactors(y)
  design = model.matrix(~group)
  y = estimateDisp(y, design)
  
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  edgeR_deg_df <- qlf$table
  edgeR_res <- edgeR_deg_df
  
  return(edgeR_res)
}

adj_adjust_p = function(edgeR_res) {
  # Returns edgeR_res data frame with an additional column with adjusted p values. 
  adj_p_vals <- p.adjust(edgeR_res$PValue, method="fdr", nrow(edgeR_res))
  edgeR_res$adj_PValue = adj_p_vals
  return(edgeR_res)
}

filter_edgeR_markers_pval <- function(edgeR_res, adj_p_threshold) {
  # Filters edgeR_res data frame by adjusted p-value < adj_p_threshold.
  filtered_res = edgeR_res[which(edgeR_res$adj_PValue < adj_p_threshold), ]
  return(filtered_res)
}

filter_edgeR_markers <- function(edgeR_res, logFC_threshold) {
  # Filters edgeR_res data frame by logFC < logFC_threshold.
  filtered_res = edgeR_res[which(edgeR_res$logFC < logFC_threshold), ]
  return(filtered_res)
}

complete_edgeR_pipeline <- function(pair_mat, rep_1, rep_2, adj_p_threshold, logFC_threshold, num_degs) {
  # Performs edgeR between 2 samples and returns data frame with all DEGs within adj_p_threshold and
  # above logFC_threshold. Also returns the top num_degs from this data frame. 
  # 
  # Args:
  #   Refer to edgeR_pipeline docstring.
  #   adj_p_threshold is set to 0.05; logFC_threshold is set to 0; num_degs is often 50.
  #
  # Returns:
  #   ans_list: list with 2 values-
  #     1) data frame with all DEGs within adj_p_threshold and above logFC_threshold.
  #     2) character vector of top num_degs from this data frame.
  
  ans_list = list()
  # Do edgeR.
  res <- edgeR_pipeline(pair_mat, rep_1, rep_2)
  # Add adj p_value column.
  res_p <- adj_adjust_p(res)
  # Filter based on adj p_vals, typically < 0.05
  res_p <- filter_edgeR_markers_pval(res_p, adj_p_threshold)
  # Filter and Order by logFC
  res_p <- filter_edgeR_markers(res_p, logFC_threshold)
  res_p_2 <- res_p[order(res_p$logFC), ]
  degs <- rownames(res_p_2)[1:num_degs]
  # We might also be interested in all relevant DEGs.
  ans_list[[1]] <- res_p_2
  ans_list[[2]] <- degs
  return(ans_list)
}

get_av_exp_column <- function(tf_df) {
  # Given a tf_df, returns a modified version of it with a column of average expression of its DEG TFs.
  #
  # Args:
  #   tf_df: data frame with rows = TFs and cols = cells.
  tf_df$mean_exp <- rowMeans(tf_df)
  modified_tf_df <- tf_df %>% dplyr::select(mean_exp) %>% arrange(-mean_exp)
  modified_tf_df <- rownames_to_column(modified_tf_df, "KH.gene.model")
  return(modified_tf_df)
}

get_named_tf_df <- function(modified_tf_df, filtered_tf_df) {
  # Returns the names of khids in modified_tf_df (which has mean expression values
  # of all TFs)
  named_tf_df <- modified_tf_df %>% inner_join(filtered_tf_df, by = "KH.gene.model") %>% dplyr::select(KH.gene.model, Ghost.Name, mean_exp)
  return(named_tf_df)
}

make_named_mean_TF_exp_df <- function(cells) {
  # Returns a data frame with KH.gene.model, Ghost.id and mean expression of all 808 TFs
  # in given cells.
  #
  # Args:
  #   cells: self-explanatory.
  #     Either Precursor or Mother Pseudostate cells.
  #
  # Returns:
  #   named_mean_all_TF_exp_df: as described above.
  # 
  # Use Case: The top 50% of most highly expressed of given DEG TFs in either the
  #   Precursor or Mother will be considered Candidate Cocktail TFs.
  all_TF_exp_df <- as.data.frame(URD_data[tfs_in_URD, cells]) 
  mean_all_TF_exp_df <- get_av_exp_column(all_TF_exp_df)
  named_mean_all_TF_exp_df <- get_named_tf_df(mean_all_TF_exp_df, filtered_tf_df)
  return(named_mean_all_TF_exp_df)
}

munge_combined_df <- function(combined_df) {
  # Cleans up names of columns.
  # Subroutine in "get_cutoff_df".
  munged_combined_df <- combined_df %>% rename_all(~stringr::str_replace_all(., ".x", "_precursor")) %>%
    rename_all(~stringr::str_replace_all(., ".y", "_mother")) %>%
    rename_all(~stringr::str_replace_all(., "mean__precursorp", "mean_exp"))
  return(munged_combined_df)
}

get_cutoff_df <- function(precursor_mean_TF_exp_df, mother_mean_TF_exp_df, khids, cutoff) {
  # Get intersection and union of highly mean expressed TFs in 2 populations at a given 
  #   cutoff. The union is considered the Candidate Cocktail TFs.
  #
  # Args: 
  #   precursor/mother_mean_TF_exp_df: data frame with cols = Kh.gene.model, Ghost.Name and mean_exp
  #     of all TFs in corresponding populations.
  #   khids: khids that will filter mean_TF_exp_dfs.
  #     e.g. BTN_unfiltered_candidate_cocktail_tfs <- unique(c(named_BTN_precursor_v_sibling_pseudostate_TF_DEGs_df$KH.gene.model,
  #                                                            named_BTN_precursor_v_mother_pseudostate_TF_DEGs_df$KH.gene.model,
  #                                                            named_BTN_mother_v_sibling_pseudostate_TF_DEGs_df$KH.gene.model))
  #   cutoff: float; percentage of highly expressed TFs considered.
  #
  # Returns:
  #   cutoff_df_list: [["intersection]] = has intersection of top 50% most highly expressed khids in
  #                       Precursor and Mother.
  #                   [["union]] = has union of top 50% most highly expressed khids in
  #                       Precursor and Mother.
  #                     Note: We will consider the union the Candidate Cocktail Tfs
  
  filtered_precursor_mean_df <- precursor_mean_TF_exp_df %>% filter(KH.gene.model %in% khids)
  filtered_mother_mean_TF_df <- mother_mean_TF_exp_df %>% filter(KH.gene.model %in% khids)
  
  precursor_cutoff_filtered_df <- filtered_precursor_mean_df[1:(cutoff*(length(khids))), ]
  mother_cutoff_filtered_df <- filtered_mother_mean_TF_df[1:(cutoff*(length(khids))), ]
  
  cutoff_df_intersection <- munge_combined_df(inner_join(precursor_cutoff_filtered_df, mother_cutoff_filtered_df, by="KH.gene.model"))
  cutoff_df_union <- munge_combined_df(full_join(precursor_cutoff_filtered_df, mother_cutoff_filtered_df, by="KH.gene.model"))
  
  cutoff_df_list <- list()
  cutoff_df_list[["intersection"]] <- cutoff_df_intersection
  cutoff_df_list[["union"]] <- cutoff_df_union
  
  return(cutoff_df_list)
}

get_TF_num <- function(non0_selected_TF_namenum) {
  # Gets number attached to TF i.e. TF1 -> 1.
  # Subroutine in "get_group_lassoed_tf_list".
  len <- nchar(non0_selected_TF_namenum)
  non0_selected_TF_num <- substring(non0_selected_TF_namenum, 3, len)
  return(non0_selected_TF_num)
}

get_group_lassoed_tf_list <- function(grepreg_cv_fit_object, candidate_cocktail_tfs) {
  # Relevant for Group Lasso
  # 
  # Given a grepreg_cv_fit_object, function iterates through all columns = lambdas,
  # gets the non0 betas and corresponding candidate cocktail TFs ("group lassoed" TFs).
  #
  # Args:
  #   grepreg_cv_fit_object: output of cv.grpreg.
  #   candidate_cocktail_tfs: as previously described.
  #
  # Returns:
  #   group_lassoed_tf_list: list with keys = lambdas
  #     and values = character vectors of TF groups selected at corresponding value of lambdas
  #     ("group lassoed" TFs).
  
  group_lassoed_tf_list <- list()
  # This is a df with dims pk x 100 (num lambdas)
  cv_grepreg_betas <- grepreg_cv_fit_object$fit$beta
  for (i in 1:ncol(cv_grepreg_betas)) {
    
    # betas for given lambda indexed by i.
    lambda_cv_grepreg_betas <- cv_grepreg_betas[, i]
    # without the intercept
    lambda_cv_grepreg_betas <- lambda_cv_grepreg_betas[2:length(lambda_cv_grepreg_betas)]
    non0_beta_indices <- which(lambda_cv_grepreg_betas != 0)
    # actual value of the betas.
    non0_betas <- lambda_cv_grepreg_betas[non0_beta_indices]
    # "TF1_cell_1" "TF4_cell_1"
    non0_beta_names <- names(non0_betas)
    # Each entry: "TF1", "cell", "1"
    non0_beta_split_list <- str_split(non0_beta_names, "_")
    # "TF1" "TF4" (all names for all non0_betas e.g. 553)
    non0_all_selected_namenums <- sapply(non0_beta_split_list,"[[", 3)
    # "TF1" "TF4" (Get the unique ones e.g. 5)
    # These are the groups.
    non0_selected_TF_namenums <- unique(non0_all_selected_namenums)
    non0_selected_TF_nums <- as.numeric(sapply(non0_selected_TF_namenums, get_TF_num))
    # Get the names of the "group lassoed" Candidate Cocktail TFs.
    group_lassoed_candidate_cocktail_tfs <- candidate_cocktail_tfs[non0_selected_TF_nums]
    list_entry_name <- paste0("lambda_", i)
    group_lassoed_tf_list[[list_entry_name]] <- group_lassoed_candidate_cocktail_tfs
  }
  return(group_lassoed_tf_list)
}

get_non_group_lassoed_TFs <- function(group_lassoed_TF_vec, candidate_cocktail_tfs) {
  # Relevant for Group Lasso
  # 
  # Function returns all the TF groups that were not selected at various values of lambda.
  non_group_lassoed_TF_vec <- setdiff(candidate_cocktail_tfs, group_lassoed_TF_vec)
  return(non_group_lassoed_TF_vec)
}

get_TF_importances_list <- function(grepreg_cv_fit_object) {
  # Relevant for Group Lasso
  # 
  # Given grepreg_cv_fit_object, function iterates through all columns = lambdas
  # and gets corresponding betas.
  #
  # Args:
  #   grepreg_cv_fit_object: output of cv.grpreg
  #
  # Returns:
  #   all_TF_importances_lists: list with keys = lambdas and values = list of importances for non0 TFs, 
  #     where keys = TF_nums (e.g. TF1) and values = DEG-specific coefficients.
  all_TF_importances_lists <- list()
  
  cv_grepreg_betas <- grepreg_cv_fit_object$fit$beta
  for (i in 1:ncol(cv_grepreg_betas)) {
    # betas for given lambda indexed by i.
    lambda_cv_grepreg_betas <- cv_grepreg_betas[, i]
    # without the intercept
    lambda_cv_grepreg_betas <- lambda_cv_grepreg_betas[2:length(lambda_cv_grepreg_betas)]
    # Get the non0 values.
    non0_beta_indices <- which(lambda_cv_grepreg_betas != 0)
    #  DEG_1_TF1   DEG_1_TF2
    # 0.001602526 -0.002237238
    non0_betas <- lambda_cv_grepreg_betas[non0_beta_indices]
    # "DEG_1_TF1" "DEG_2_TF4"
    non0_beta_names <- names(non0_betas)
    # Each entry: "DEG" "1"   "TF1"
    non0_beta_split_list <- str_split(non0_beta_names, "_")
    # "TF1" "TF4"
    non0_all_selected_namenums <- sapply(non0_beta_split_list,"[[", 3)
    
    num_betas <- length(non0_betas)
    # cat("length = ", length(num_betas))
    if (num_betas != 0) {
      # cat("num betas = ", num_betas, "\n")
      # "TF1" "TF2" of ONLY the non0 betas.
      names(non0_betas) <- non0_all_selected_namenums
      # TF1         TF2         TF1         
      # 0.10487420  0.10487420 0.20290965
      # -> split -> 
      # $TF1 
      # TF1       TF1 
      # 0.1048742 0.2029097 
      
      # $TF2
      # 0.10487420
      
      # Each TF has k betas. 
      # Thus, irrespective of the number of non0 betas, they will always
      # be grouped in sizes of k! 
      TF_importances_list <- split(non0_betas, names(non0_betas))
      TF_importances_list <- lapply(TF_importances_list, unname)
    }
    else {
      TF_importances_list <- list()
    }
    list_key <- paste0("lambda_", i)
    all_TF_importances_lists[[list_key]] <- TF_importances_list
  }
  return(all_TF_importances_lists)
}

get_total_summed_TF_importance_list <- function(TF_importances_list) {
  # Relevant for Group Lasso
  # 
  # For each TF, function takes absolute value of each DEG-specific beta and sums them.
  # This is the Total Importance Score. 
  #
  # Args:
  #   TF_importances_list: as previously described.
  #
  # Returns:
  #   summed_TF_importances_list: list with keys = TFs and values = Total Importance Scores.
  
  # Some lambdas may not correspond to betas.
  if (length(TF_importances_list) == 0) {
    return(0)
  }
  else {
    summed_TF_importances_list <- lapply(TF_importances_list, abs)
    summed_TF_importances_list <- lapply(summed_TF_importances_list, sum)
    return(summed_TF_importances_list)
  }
  
}
munge_total_TF_importances <- function(total_TF_importances, candidate_cocktail_tfs) {
  # Replaces generic TF_nums with TF Ortholog names and makes a data frame
  # sorted by their Total Importance Scores.
  # Subroutine of "make_sorted_TF_importance_df".
  #
  # Args:
  #   total_TF_importances: named vector of TF and beta value
  #     e.g. BTN_total_TF_importances
  #     TF8         -> isl          
  #     10.12289309 -> 10.12289309
  #   candidate_cocktail_tfs: as previously described.
  #
  # Returns:
  #   total_TF_importances_df: data frame with rownames as TF Ortholog names
  #     and 1 column with non-increasing values of Total Importance Scores.
  #     Rows are sorted by Total Importance Scores.
  total_TF_importances_nums <- as.numeric(sapply(names(total_TF_importances), get_TF_num))
  influential_cocktail_tfs <- candidate_cocktail_tfs[total_TF_importances_nums]
  names(total_TF_importances) <- influential_cocktail_tfs
  # Make the df.
  total_TF_importances_df <- data.frame(as.list(total_TF_importances))
  total_TF_importances_df <- t(data.frame(as.list(total_TF_importances)))
  colnames(total_TF_importances_df) <- "Total Importance"
  return(total_TF_importances_df)
}

make_sorted_TF_importance_df <- function(summed_TF_importances_list, candidate_cocktail_tfs) {
  # Relevant for Group Lasso
  # 
  # Args: 
  #   summed_TF_importances_list: list with keys = lambdas and values as lists with 
  #     keys = TF_nums and values = Total Importance Scores.
  #   candidate_cocktail_tfs: as previously described.
  #
  # Returns:
  #   sorted_TF_importance_df: data frame with rows = TF Ortholog names 
  #     and columns = Total Importance Scores. 
  #     Rows are sorted by Total Importance Score values. 
  #
  # Note: This function is applies to
  summed_total_TF_importances <- sort(unlist(summed_TF_importances_list), decreasing = TRUE)
  sorted_TF_importance_df <- munge_total_TF_importances(summed_total_TF_importances, candidate_cocktail_tfs)
  return(sorted_TF_importance_df)
}

make_coef_df <- function(coef_list, candidate_cocktail_tf_names, precursor_deg_names) {
  # Relevant for GMML
  #
  # Makes a coef_df, which has rows as Candidate Cocktail Precursor DEGs and 
  # columns as Precursor DEGs with entries as regression coefficients present in
  # coef_list.
  # 
  # Args:
  #   coef_list: list that is the output of fitting gaussian multivariate multi-regression
  #     from cv.glmnet at a certain value of lambda.
  #     e.g. BTN_coef_list <- coef(BTN_fit, s = BTN_lambda)
  #       where BTN_fit <- cv.glmnet(glm_BTN_X_mat, glm_BTN_Y_mat...)
  #   candidate_cocktail_tf_names: self-explanatory.
  #   precursor_deg_names: self-explanatory.
  #
  #   Note: Use case of above 2 arguments is to munge and clarify row and column
  #     names.
  #
  # Returns:
  #   coef_df: as described above.
  #
  # Note: Use cases of coef_df are to get average coefficient values and
  #   total importances of Candidate Cocktail TFs.
  coef_df <- do.call(cbind, coef_list)
  coef_df <- as.data.frame(coef_df)
  coef_df <- coef_df[-1, ]
  
  rownames(coef_df) <- candidate_cocktail_tf_names
  colnames(coef_df) <- precursor_deg_names
  return(coef_df)
}

get_average_coefficient_df <- function(coef_df) {
  # Relevant for GMML
  #
  # Given a coef_df, that is a data frame that contains regression coefficients
  # from all Candidate Cocktail TF-DEG pairs, returns average of absolute values of
  # coefficient values.
  #
  # Args:
  #   coef_df: as explained above.
  #
  # Returns:
  #   average_coef_df: Contains average coefficient values for all Candidate
  #     Cocktail TFs. The TFS are sorted in nonincreasing order of coefficient values.
  
  coef_df <- abs(coef_df)
  # Get average coefficients
  average_cocktail_coef <- rowMeans(coef_df)
  average_cocktail_coef <- sort(rowMeans(coef_df), decreasing = TRUE)
  average_cocktail_coef_df <- data.frame(average_cocktail_coef)
  return(average_cocktail_coef_df)
}

get_coefficients_distribution_plot <- function(coef_df) {
  # Relevant for GMML
  #
  # Given a coef_df, that is a data frame that contains regression coefficients
  # from all Candidate Cocktail TF-DEG pairs, visualizes the distribution 
  # of these coefficients in a ggplot.
  # 
  # Args:
  #   coef_df: as explained above.
  #
  # Returns:
  #   coef_distribution_plot: ggplot that visualizes these coefficient values.
  coef_df <- abs(coef_df)
  
  # Make coef_distribution_plot
  gathered_coef_df <- as.data.frame(t(coef_df)) %>% gather()
  coef_distribution_plot <- ggplot(gathered_coef_df, aes(value)) + 
    geom_histogram(aes(y = stat(width*density))) + 
    ylim(0, 0.30) +
    xlim(min(coef_df)-0.1*min(coef_df), max(coef_df)+0.1*max(coef_df)) +
    facet_wrap(~key, scales = "free_x")
  
  return(coef_distribution_plot)
}

get_ddN_coefficients_distribution_plot <- function(coef_df) {
  # Relevant for GMML
  #
  # Given a coef_df, that is a data frame that contains regression coefficients
  # from all Candidate Cocktail TF-DEG pairs, visualizes the distribution 
  # of these coefficients in a ggplot.
  # 
  # Args:
  #   coef_df: as explained above.
  #
  # Returns:
  #   coef_distribution_plot: ggplot that visualizes these coefficient values.
  coef_df <- abs(coef_df)
  
  # Make coef_distribution_plot
  gathered_coef_df <- as.data.frame(t(coef_df)) %>% gather()
  coef_distribution_plot <- ggplot(gathered_coef_df, aes(value)) + 
    geom_histogram(aes(y = stat(width*density))) + 
    ylim(0, 1) +
    xlim(min(coef_df)-0.1*min(coef_df), max(coef_df)+0.1*max(coef_df)) +
    facet_wrap(~key, scales = "free_x")
  
  return(coef_distribution_plot)
}

get_total_GMML_importance_df <- function(coef_df) {
  # Relevant for GMML
  #
  # Given a coef_df, a data frame that contains regression coefficients
  # from all Candidate Cocktail TF-DEG pairs, for a TF, sums the absolute values
  # of its regression coefficients with all DEGs. 
  #
  # Note: coef_df has TFs as rows, khids as columns, and entries as 
  #   coefficient values.
  #
  # Args:
  #   coef_df: as described.
  # 
  # Returns:
  #   total_GMML_importance_df: a data frame sorted in nonincreasing order based
  #   on total
  
  coef_df <- abs(coef_df)
  gathered_coef_df <- as.data.frame(t(coef_df)) %>% gather()
  
  total_GMML_importance_df <- gathered_coef_df %>% group_by(key) %>% summarise(total_GMML_importance = sum(value)) %>% arrange(desc(total_GMML_importance))
  
  return(total_GMML_importance_df)
}

make_grplasso_X_mat <- function(TFs_in_cells_mat, n, p, k) {
  # Makes grplasso_X_mat_2, a matrix of dimensions nk x pk that has Matrix S (TFs_in_cells_mat)
  # repeated k times along its diagonal. All other entries are 0. 
  #
  # Args:
  #   TFs_in_cells_mat: matrix of dimensions n x p (also referred to as Matrix S).
  #     that has expression of p TFs in n cells.
  #     e.g. BTN_TFs_in_precursor_pseudostate_mat <- t(URD_data[BTN_candidate_cocktail_khids, BTN_precursor_pseudostate_cells])
  #   n, p, k = number of cells, TFs, and DEGs.
  #
  # Returns:
  #   grplasso_X_mat_2: as described above.
  S_matrix <- TFs_in_cells_mat
  zero_matrix <- matrix(0, n, p)
  # This is a dummy matrix- we will remove it at the end.
  grplasso_X_mat <- matrix(1, n, p*k)
  
  for (i in 1:k) {
    L_matrix <- do.call(cbind, replicate(i-1, zero_matrix, simplify=FALSE))
    R_matrix <- do.call(cbind, replicate(k-i, zero_matrix, simplify=FALSE))
    grplasso_mat_row <- do.call(cbind, list(L_matrix, S_matrix, R_matrix))
    grplasso_X_mat <- rbind(grplasso_X_mat, grplasso_mat_row)
  }
  grplasso_X_mat_2 <- grplasso_X_mat[-c(1:n), ]
  
  # Rename the rows such that they reflect the cell and DEG.
  grplasso_X_mat_row_cells <- rep(paste0("_cell_", seq(1, n)), k)
  DEG_indices <- rep(c(1:k), each = n)
  new_grplasso_X_mat_2_rownames <- paste0("DEG_", DEG_indices, grplasso_X_mat_row_cells)
  rownames(grplasso_X_mat_2) <- new_grplasso_X_mat_2_rownames
  
  # Rename the columns such that they reflect the TF and DEG.
  grplasso_X_mat_col_TFs <- rep(paste0("_TF", seq(1, p)), k)
  DEG_indices_cols <- rep(c(1:k), each = p)
  new_grplasso_X_mat_2_colnames <- paste0("DEG_", DEG_indices_cols, grplasso_X_mat_col_TFs)
  colnames(grplasso_X_mat_2) <- new_grplasso_X_mat_2_colnames
  
  grplasso_X_mat_2 <- as.matrix(grplasso_X_mat_2)
  return(grplasso_X_mat_2)
}

make_grplasso_Y_mat <- function(DEGs_in_cells_mat, cells, p, k) {
  # Makes Group Lasso Matrix Y with dimensions nk x 1 that contains the expression of all k DEGs 
  # in all n Cells in the following way: The first n rows contain expression 
  # of DEG 1 in all n cells; the next n rows contain expression of DEG 2 in all n cells; 
  # so and so forth for k DEGs.
  #
  # Args:
  #   DEGs_in_cells_mat: matrix of dimensions of n x k that has expression 
  #     of k DEGs in n cells.
  #     BTN_precursor_DEGs_in_precursors_mat <- t(URD_data[cleaned_BTN_precursor_pseudostate_DEGs, BTN_precursor_pseudostate_cells])
  #
  # Returns:
  #   grplasso_Y_mat: as described above.
  n <- length(cells)
  grplasso_Y_vec <- c()
  
  for (deg in 1:k) {
    for (cell in cells) {
      deg_cell_exp <- DEGs_in_cells_mat[cell, deg]
      # cat("exp = ", deg_cell_exp)
      grplasso_Y_vec <- c(grplasso_Y_vec, deg_cell_exp)
    }
  }
  grplasso_Y_mat <- as.matrix(grplasso_Y_vec, ncol=1, nrow=(n*k))
  
  # Change rownames to reflect cell and DEG.
  grplasso_Y_DEG_indices <- rep(c(1:k), each = n)
  grplasso_Y_degs <- paste0("DEG_", grplasso_Y_DEG_indices)
  grplasso_Y_cells <- rep(paste0("_cell_", seq(1, n)), k)
  new_grplasso_Y_rownames <- paste0(grplasso_Y_degs, grplasso_Y_cells)
  
  # Update the rownames.
  rownames(grplasso_Y_mat) <- new_grplasso_Y_rownames
  return(grplasso_Y_mat)
}

make_group_lassoed_TF_string <- function(group_lassoed_TFs) {
  # Given group_lassoed_TFs, a character vector of TFs, 
  # e.g. COE HNF Neurogenin...
  #
  # ... returns "group_lassoed_TF_string", a comma delimited string.
  # e.g. "COE,HNF,Neurogenin". 
  #
  # Subroutine in "make_lambda_to_grouped_tf_df".
  group_lassoed_TF_string <- paste(group_lassoed_TFs, collapse = ",")
  return(group_lassoed_TF_string)
}

make_lambda_to_grouped_tf_df <- function(combined_subsampling_lists) {
  # Given combined_subsampling_lists, e.g. c(BTN_subsampling_pipeline_1, BTN_subsampling_pipeline_1)
  #   where key = run and value is a list of cv_grpreg, group_lasso, and non_group_lasso
  #   and group_lasso/non_group_lasso are lists with key = lambda and value = TFs...
  # (Note: Run refers to a subsample)
  # 
  # ... returns lambda_to_grouped_tf_df: data frame where rownames have run and lamdba 
  #   e.g. "run_13_lambda_2" and column has TF_strings of TF groups selected at
  #         corresponding run and lambda values.
  # 
  # Subroutine of "make_combined_subsampled_tf_df".
  lambda_to_grouped_tf_list <- list()
  for (i in 1:length(combined_subsampling_lists)) {
    # has cv_grpreg, group_lasso, non_group_lasso lists
    subsample_run_output <- combined_subsampling_lists[[i]]
    # goes into group_lasso list- has lambda_*, TFs
    subsample_run_group_lassoed_TF_list <- subsample_run_output$group_lasso
    for (j in 1:length(subsample_run_group_lassoed_TF_list)) {
      lambda_to_grouped_tf_list_key <- paste0("run_", i, "_lambda_", j)
      group_lassoed_TFs <- subsample_run_group_lassoed_TF_list[[j]]
      group_lassoed_TF_string <- make_group_lassoed_TF_string(group_lassoed_TFs)
      lambda_to_grouped_tf_list[[lambda_to_grouped_tf_list_key]] <- group_lassoed_TF_string
    }
  }
  lambda_to_grouped_tf_df <- as.data.frame(do.call(rbind, lambda_to_grouped_tf_list))
  return(lambda_to_grouped_tf_df)
}

make_combined_subsampled_tf_df <- function(subsample_file_dir) {
  # Creates combined_subsampled_tf_df, from directory containing
  # R objects with results of cv.grpreg runs. 
  # 
  # Args:
  #   subsample_file_dir: directory as described. 
  # 
  # Returns:
  #   combined_subsampled_tf_df: data frame with row = run_*_lambda_*, 
  #     where run = the round of subsampling/run of group lasso and 
  #     lambda = one of 100 values of lambda used during cross-validation; 
  #     and column value = string that is a concatenated name of the set of
  #     TF Groups selected at the corresponding value of lambda ("TF-string")
  
  # Read in all pipeline objects. 
  subsample_objects <- list.files(path = subsample_file_dir)
  cat("Getting objects\n")
  # Makes a list of lists (each is from a separate sample/run).
  subsample_objects_list <- lapply(paste0(subsample_file_dir, subsample_objects), function(x) readRDS(x))
  cat("Flattening list\n")
  # Turn into 1 list.
  flattened_subsample_objects_list <- purrr::flatten(subsample_objects_list)
  cat("Making combined subsampled grouped TF DF\n")
  # Make grouped_tf_df (rows = run_n_lambda_p; "TF String")
  combined_subsampled_tf_df <- make_lambda_to_grouped_tf_df(flattened_subsample_objects_list)
  return(combined_subsampled_tf_df)
}

write_out_combined_subsampled_tf_df <- function(combined_subsampled_tf_df, out_dir, out_filename) {
  # Write out combined_subsampled_tf_df as a csv file.
  # This is used for the direct output of make_lambda_to_grouped_tf_df.  
  # 
  # Note: Each line consists of run_*_lambda_* and TF_string.
  #   e.g. "run_1_lambda_2","Neurogenin,COE,HNF6"
  #
  # Args:
  #   combined_subsampled_tf_df: as previously described.
  #   out_dir: directory to write csv file out to.
  #   out_filename: name of csv file.
  # 
  # Returns:
  #   csv file as described.
  #
  # Use case: Will be used as an input to "subsampling_confidence_interval.py"
  write.table(combined_subsampled_tf_df, paste0(out_dir, out_filename), row.names = TRUE, col.names = FALSE, sep = ",", quote = TRUE)
}

write_out_group_lasso_tf_df <- function(group_lasso_tf_df, out_dir, out_filename) {
  # Writes group_lasso_tf_df out as csv file.
  #
  # Args:
  #   group_lasso_tf_df: a data frame with rownames = "run_*_lambda_*" and column is TF_string.
  #   out_dir, out_name: self-explanatory.
  #
  # Use case: 
  #   This is used for the direct output of make_lambda_to_grouped_tf_df and will be
  #   used to get a confidence score for predicted Cocktails. 
  #
  # Note: CT_group_lassoed_tf_list has its own function- 
  #   write_out_munged_CT_group_lassoed_tf_list.
  write.table(group_lasso_tf_df, paste0(out_dir, out_filename), row.names = TRUE, col.names = FALSE, sep = ",", quote = TRUE)
}

write_out_CT_group_lassoed_tf_list <- function(CT_group_lassoed_tf_list, out_dir, out_filename) {
  # Write out a CT_group_lassoed_tf_list as a csv file.
  #
  # Args:
  #   CT_group_lassoed_tf_list: list with keys = lambdas
  #     and values = character vectors of TF groups selected at corresponding value of lambdas
  #     ("group lassoed" TFs).
  #   out_dir: directory to write csv file out to.
  #   out_filename: name of csv file.
  # 
  # Routine:
  #   Makes a data frame with rownames as "lambda_*" and column values as strings
  #     of group lassoed TFs. e.g. "lambda_2","Neurogenin,COE,HNF6"
  #
  #   Then writes df out as a csv file in which every line consists of "lambda_*" and TF_string.
  #
  # Returns:
  #   csv file as described.
  #
  # Use case: Will be used as an input to "subsampling_confidence_interval.py"
  munged_CT_group_lassoed_tf_list <- lapply(CT_group_lassoed_tf_list, make_group_lassoed_TF_string)
  munged_CT_group_lassoed_tf_df <- as.data.frame(do.call(rbind, munged_CT_group_lassoed_tf_list))
  write_out_group_lasso_tf_df(munged_CT_group_lassoed_tf_df, out_dir, out_filename)
}

get_dmbx_expressing_cells <- function(target_cells, pseudotimes) {
  # Finds dmbx-expressing cells within target_cells.
  #
  # Args:
  #   target_cells: character vector of cell barcodes.
  #   pseudotimes: self-explanatory.
  #
  # Returns:
  #   ordered_dmbx_exp_cells: character vector of cell barcodes
  #     of cells that express dmbx+ and ordered by nondecreasing
  #     pseudotime values.
  #
  # Use case: GLRA segment or GLRA Precursor Pseudostate cells 
  #   are inputted to find ddN Precursor Pseudostate.
  dmbx_exp_df <- as.data.frame(URD_data[dmbx, ])
  dmbx_exp_df <- add_rownames(dmbx_exp_df, var = "rowname")
  colnames(dmbx_exp_df) <- c("cell", "dmbx_exp")
  filtered_dmbx_exp_df <- dmbx_exp_df %>% filter(dmbx_exp != 0, cell %in% target_cells)
  unordered_dmbx_exp_cells <- filtered_dmbx_exp_df$cell
  # Add a column with pseudotime values
  cell_pseudotimes <- pseudotimes[unordered_dmbx_exp_cells]
  filtered_dmbx_exp_df$cell_pseudotimes <- cell_pseudotimes
  # Order by pseudotimes
  ordered_filtered_dmbx_exp_df <- filtered_dmbx_exp_df %>% arrange(cell_pseudotimes)
  ordered_dmbx_exp_cells <- ordered_filtered_dmbx_exp_df$cell
  return(ordered_dmbx_exp_cells)
}

make_precursor_mats_list <- function(subsampled_precursor_cells, TFs_in_precursor_mat, precursor_degs_in_precursors_mat) {
  # Relevant for ddN subsampling and runs of Group Lasso with cross-validation.
  # Subroutine of "grepreg_routine".
  #
  # Given subsampled_precursor_cells (and candidate_cocktail_tfs and cleaned_precursor_DEGs), 
  # function subsets TFs_in_precursor_mat and precursor_degs_in_precursors_mat using those cells.
  # e.g. ddN_TFs_in_precursor_pseudostate_mat and ddN_precursor_DEGs_in_precursors_mat will be subsetted
  #       by some number of ddN_precursor_pseudostate_cells.
  #
  # Returns:
  #   precursor_mats_list: list with keys = TFs and prec_degs.
  #     Values are subsets of TFs_in_precursor_mat and precursor_degs_in_precursors_mat.
  grplasso_TFs_in_precursor_mat <- as.matrix(TFs_in_precursor_mat[subsampled_precursor_cells, ])
  grplasso_precursor_degs_in_precursors_mat <- as.matrix(precursor_degs_in_precursors_mat[subsampled_precursor_cells, ])
  
  precursor_mats_list <- list()
  precursor_mats_list[["TFs"]] <- grplasso_TFs_in_precursor_mat
  precursor_mats_list[["prec_degs"]] <- grplasso_precursor_degs_in_precursors_mat
  
  return(precursor_mats_list)
}

seed_do_cv_grpreg <- function(grplasso_X_mat_2, grplasso_Y_mat, p, k, grpreg_seed) {
  # Function performs Group Lasso with cross-validation with the given arguments, 
  # namely a grepreg seed (int), which allows for reproducibility of results.
  # Subroutine of "grepreg_routine".
  grplasso_groups <- rep(seq(1, p), k)
  group_id <- grplasso_groups
  grplasso_cv_grepreg_fit <- cv.grpreg(grplasso_X_mat_2, grplasso_Y_mat, group_id, penalty="grLasso", seed=grpreg_seed)
  return(grplasso_cv_grepreg_fit)
}
 
grepreg_routine <- function(subsampled_precursor_cells, 
                            candidate_cocktail_tfs, 
                            cleaned_precursor_DEGs,
                            TFs_in_precursor_mat, 
                            prec_degs_in_precursors_mat,
                            grpreg_seed) {
  # Function performs Group Lasso with cross-validation for subsampled_precursor_cells
  # with given arguments, namely a grepreg_seed (as previously described).
  #
  # Returns:
  #   output: list with 3 keys = "cv_grepreg", "group_lasso", and "non_group_lasso".
  #     "cv_grpreg has value" cv.grpreg object. "group_lasso" and "non_group_lasso" are
  #     lists with keys = lambda and values = TF groups selected/not selected respectively.
  
  n_samples <- length(subsampled_precursor_cells)
  p <- length(candidate_cocktail_tfs)
  k <- length(cleaned_precursor_DEGs)
  
  output <- list()
  # Step 1
  precursor_mats_list <- make_precursor_mats_list(subsampled_precursor_cells,
                                                  TFs_in_precursor_mat,
                                                  prec_degs_in_precursors_mat)
  
  subsampled_grplasso_TFs_in_precursor_mat <- precursor_mats_list$TFs
  subsampled_grplasso_prec_degs_in_precursors_mat <- precursor_mats_list$prec_degs
  
  # Step 2
  cat("Making grplasso_X_mat_2\n")
  subsampled_grplasso_X_mat_2 <- make_grplasso_X_mat(subsampled_grplasso_TFs_in_precursor_mat, 
                                                     n_samples, p, k)
  # Step 3
  cat("Making grplasso_Y_vec\n")
  subsampled_grplasso_Y_vec <- make_grplasso_Y_mat(prec_degs_in_precursors_mat, 
                                                   subsampled_precursor_cells, p, k) 
  # Step 4
  cat("Doing cv_grepreg\n")
  subsampled_grplasso_cv_grepreg <- seed_do_cv_grpreg(subsampled_grplasso_X_mat_2, 
                                                      subsampled_grplasso_Y_vec, p, k, grpreg_seed)
  output[["cv_grepreg"]] <- subsampled_grplasso_cv_grepreg
  
  # Step 5
  cat("Munging all results\n")
  subsampled_group_lassoed_tf_list <- get_group_lassoed_tf_list(subsampled_grplasso_cv_grepreg, 
                                                                candidate_cocktail_tfs)
  output[["group_lasso"]] <- subsampled_group_lassoed_tf_list
  
  subsampled_non_group_lassoed_TF_list <- lapply(subsampled_group_lassoed_tf_list, 
                                                 get_non_group_lassoed_TFs,
                                                 candidate_cocktail_tfs)
  output[["non_group_lasso"]] <- subsampled_non_group_lassoed_TF_list
  return(output)
}

