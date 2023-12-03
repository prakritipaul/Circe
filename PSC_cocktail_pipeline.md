PSC Cocktail Inference Pipeline
================
This RMarkdown includes all analyses resulting from the computational
methods developed in Chapter 3 of Prakriti Paul Chacha’s PhD Thesis
‘Computational Approaches to Reprogram Neuronal Cell Identities in Ciona
intestinalis’.

- <a href="#brief-introduction" id="toc-brief-introduction">Brief
  Introduction</a>
- <a href="#import-relevant-libraries"
  id="toc-import-relevant-libraries">Import Relevant Libraries</a>
- <a
  href="#1-identify-precursor-mother-and-their-siblings-segment-cell-populations-on-virtual-neural-lineage-tree"
  id="toc-1-identify-precursor-mother-and-their-siblings-segment-cell-populations-on-virtual-neural-lineage-tree">(1)
  Identify Precursor, Mother, and their Sibling(s) Segment Cell
  Populations on Virtual Neural Lineage Tree</a>
  - <a href="#precursor-segment" id="toc-precursor-segment">Precursor
    Segment</a>
    - <a href="#get-ids-of-psc-precursor-segment-and-its-sibling-segment"
      id="toc-get-ids-of-psc-precursor-segment-and-its-sibling-segment">Get
      IDs of PSC Precursor Segment and its Sibling Segment</a>
    - <a href="#get-psc-precursor-and-sibling-segment-cells"
      id="toc-get-psc-precursor-and-sibling-segment-cells">Get PSC Precursor
      and Sibling Segment Cells</a>
    - <a href="#order-cells-by-pseudotime"
      id="toc-order-cells-by-pseudotime">Order Cells by Pseudotime</a>
  - <a href="#mother-segment" id="toc-mother-segment">Mother Segment</a>
    - <a href="#get-ids-of-psc-mother-segment-and-its-siblings"
      id="toc-get-ids-of-psc-mother-segment-and-its-siblings">Get IDs of PSC
      Mother Segment and its Siblings</a>
    - <a href="#get-psc-mother-and-sibling-segment-cells"
      id="toc-get-psc-mother-and-sibling-segment-cells">Get PSC Mother and
      Sibling Segment Cells</a>
    - <a href="#order-cells-by-pseudotime-1"
      id="toc-order-cells-by-pseudotime-1">Order Cells by Pseudotime</a>
  - <a href="#visualize-cells-on-virtual-neural-lineage-tree"
    id="toc-visualize-cells-on-virtual-neural-lineage-tree">Visualize Cells
    on Virtual Neural Lineage Tree</a>
- <a
  href="#2-get-pc-coordinates-and-corresponding-pseudotimes-for-psc-cells"
  id="toc-2-get-pc-coordinates-and-corresponding-pseudotimes-for-psc-cells">(2)
  Get PC Coordinates and Corresponding Pseudotimes for PSC Cells</a>
  - <a href="#make-seurat-object" id="toc-make-seurat-object">Make Seurat
    Object</a>
  - <a href="#use-elbow-plot-to-choose-the-number-of-pcs-12"
    id="toc-use-elbow-plot-to-choose-the-number-of-pcs-12">Use Elbow Plot to
    choose the number of PCs (12)</a>
  - <a href="#key-variables" id="toc-key-variables">Key Variables</a>
- <a href="#3-get-pseudostates" id="toc-3-get-pseudostates">(3) Get
  Pseudostates</a>
  - <a href="#make-pc_pseudotime_dfs" id="toc-make-pc_pseudotime_dfs">Make
    PC_pseudotime_dfs</a>
  - <a href="#export-pc_pseudotime_dfs"
    id="toc-export-pc_pseudotime_dfs">Export PC_pseudotime_dfs</a>
  - <a href="#load-get_pseudostatespy" id="toc-load-get_pseudostatespy">Load
    “get_pseudostates.py”</a>
    - <a href="#create-a-new-reticulate-environment"
      id="toc-create-a-new-reticulate-environment">Create a new reticulate
      environment</a>
    - <a
      href="#install-relevant-packages-for-get_pseudostatespy-into-given-r-reticulate-environment"
      id="toc-install-relevant-packages-for-get_pseudostatespy-into-given-r-reticulate-environment">Install
      Relevant Packages for “get_pseudostates.py” into given r-reticulate
      environment.</a>
  - <a href="#key-variables-1" id="toc-key-variables-1">Key Variables</a>
  - <a href="#psc-precursor" id="toc-psc-precursor">PSC Precursor</a>
  - <a href="#psc-precursor-sibling" id="toc-psc-precursor-sibling">PSC
    Precursor Sibling</a>
  - <a href="#psc-mother" id="toc-psc-mother">PSC Mother</a>
  - <a href="#mothers-siblings" id="toc-mothers-siblings">Mother’s
    Siblings</a>
  - <a href="#visualize-pseudostates-on-virtual-neural-lineage-tree"
    id="toc-visualize-pseudostates-on-virtual-neural-lineage-tree">Visualize
    Pseudostates on Virtual Neural Lineage Tree</a>
- <a href="#4-get-differentially-expressed-genes"
  id="toc-4-get-differentially-expressed-genes">(4) Get Differentially
  Expressed Genes</a>
  - <a href="#precursor-identity-genes"
    id="toc-precursor-identity-genes">Precursor Identity Genes</a>
  - <a href="#precursor-deg-tfs" id="toc-precursor-deg-tfs">Precursor DEG
    TFs</a>
  - <a href="#mother-deg-tfs" id="toc-mother-deg-tfs">Mother DEG TFs</a>
  - <a href="#raw-matrices" id="toc-raw-matrices">Raw Matrices</a>
  - <a
    href="#precursor-vs-precursor-sibling--deg-tfs-and-precursor-identity-genes"
    id="toc-precursor-vs-precursor-sibling--deg-tfs-and-precursor-identity-genes">Precursor
    vs. Precursor Sibling- DEG TFs and Precursor Identity Genes</a>
    - <a href="#precursor-identity-genes-1"
      id="toc-precursor-identity-genes-1">Precursor Identity Genes</a>
      - <a href="#write-out-precursor-degs-precursor-identity-genes"
        id="toc-write-out-precursor-degs-precursor-identity-genes">Write out
        Precursor DEGs (Precursor Identity Genes)</a>
    - <a href="#precursor-deg-tfs-1" id="toc-precursor-deg-tfs-1">Precursor
      DEG TFs</a>
      - <a href="#write-out-tf-degs" id="toc-write-out-tf-degs">Write out TF
        DEGs</a>
  - <a href="#precursor-vs-mother--deg-tfs"
    id="toc-precursor-vs-mother--deg-tfs">Precursor vs. Mother- DEG TFs</a>
    - <a href="#write-out-tf-degs-1" id="toc-write-out-tf-degs-1">Write out TF
      DEGs</a>
  - <a href="#mother-vs-mothers-sibling--deg-tfs"
    id="toc-mother-vs-mothers-sibling--deg-tfs">Mother vs. Mother’s Sibling-
    DEG TFs</a>
    - <a href="#write-out-tf-degs-2" id="toc-write-out-tf-degs-2">Write out TF
      DEGs</a>
- <a href="#5-get-candidate-cocktail-tfs"
  id="toc-5-get-candidate-cocktail-tfs">(5) Get Candidate Cocktail TFs</a>
  - <a href="#get-mean-expression-of-all-deg-tfs"
    id="toc-get-mean-expression-of-all-deg-tfs">Get Mean Expression of all
    DEG TFs</a>
  - <a
    href="#take-the-union-of-top-50-most-highly-expressed-unfiltered-candidate-cocktail-tfs-in-precursor-and-mother-these-are-the-candidate-cocktail-tfs"
    id="toc-take-the-union-of-top-50-most-highly-expressed-unfiltered-candidate-cocktail-tfs-in-precursor-and-mother-these-are-the-candidate-cocktail-tfs">Take
    the union of top 50% most highly expressed Unfiltered Candidate Cocktail
    TFs in Precursor and Mother. These are the Candidate Cocktail TFs</a>
    - <a href="#write-out-the-candidate-cocktail-tfs"
      id="toc-write-out-the-candidate-cocktail-tfs">Write out the Candidate
      Cocktail TFs</a>
- <a href="#6-do-cocktail-tf-inference-using-group-lasso-model"
  id="toc-6-do-cocktail-tf-inference-using-group-lasso-model">(6) Do
  Cocktail TF Inference using Group Lasso Model</a>
  - <a href="#get-a-set-of-precursor-degs-excluding-deg-tfs"
    id="toc-get-a-set-of-precursor-degs-excluding-deg-tfs">Get a set of
    Precursor DEGs excluding DEG TFs</a>
  - <a
    href="#get-numbers-of-precursor-cells-n-candidate-cocktail-tfs-p-and-precusor-degs-k"
    id="toc-get-numbers-of-precursor-cells-n-candidate-cocktail-tfs-p-and-precusor-degs-k">Get
    numbers of Precursor cells (n), Candidate Cocktail TFs (p), and Precusor
    DEGs (k)</a>
  - <a href="#make-group-lasso-matrix-x"
    id="toc-make-group-lasso-matrix-x">Make Group Lasso Matrix X</a>
    - <a href="#make-matrix-s" id="toc-make-matrix-s">Make Matrix S</a>
    - <a href="#concatenate-k-s-matrices-to-make-group-lasso-matrix-x"
      id="toc-concatenate-k-s-matrices-to-make-group-lasso-matrix-x">Concatenate
      k S Matrices to make Group Lasso Matrix X</a>
    - <a href="#rename-row-and-column-names-for-convenience"
      id="toc-rename-row-and-column-names-for-convenience">Rename Row and
      Column names for Convenience</a>
      - <a href="#rename-row-names-to-reflect-cell-number-and-its-deg-status"
        id="toc-rename-row-names-to-reflect-cell-number-and-its-deg-status">Rename
        Row Names to reflect cell number and its DEG status</a>
      - <a href="#rename-column-names-to-reflect-tf-number-and-its-deg-status"
        id="toc-rename-column-names-to-reflect-tf-number-and-its-deg-status">Rename
        Column Names to reflect TF number and its DEG status</a>
    - <a href="#final-group-lasso-matrix-x"
      id="toc-final-group-lasso-matrix-x">Final Group Lasso Matrix X</a>
  - <a href="#make-group-lasso-matrix-y"
    id="toc-make-group-lasso-matrix-y">Make Group Lasso Matrix Y</a>
    - <a
      href="#make-a-matrix-that-has-the-expression-of-all-precursor-degs-in-precursor-cells-of-dimensions-k-x-n"
      id="toc-make-a-matrix-that-has-the-expression-of-all-precursor-degs-in-precursor-cells-of-dimensions-k-x-n">Make
      a matrix that has the expression of all Precursor DEGs in Precursor
      Cells of dimensions k x n</a>
    - <a href="#index-above-matrix-so-as-so-make-group-lasso-matrix-y"
      id="toc-index-above-matrix-so-as-so-make-group-lasso-matrix-y">Index
      above matrix so as so make Group Lasso Matrix Y</a>
    - <a href="#rename-row-names-to-reflect-cell-and-deg-status"
      id="toc-rename-row-names-to-reflect-cell-and-deg-status">Rename Row
      names to reflect cell and DEG status</a>
    - <a href="#final-matrix-y" id="toc-final-matrix-y">Final Matrix Y</a>
  - <a href="#run-group-lasso" id="toc-run-group-lasso">Run Group Lasso</a>
    - <a href="#variables-group-lasso-matrices-x-and-y-and-group-ids"
      id="toc-variables-group-lasso-matrices-x-and-y-and-group-ids">Variables:
      Group Lasso Matrices X and Y and Group IDs</a>
    - <a href="#do-cross-validation" id="toc-do-cross-validation">Do
      Cross-Validation</a>
    - <a href="#get-sets-of-selected-tf-groups-at-various-values-of-lambda"
      id="toc-get-sets-of-selected-tf-groups-at-various-values-of-lambda">Get
      Sets of Selected TF Groups at Various Values of Lambda</a>
      - <a
        href="#let-us-visualize-the-number-of-groups-selected-at-various-values-of-lambda"
        id="toc-let-us-visualize-the-number-of-groups-selected-at-various-values-of-lambda">Let
        us Visualize the Number of Groups Selected at Various Values of
        Lambda</a>
    - <a href="#get-sorted-total-tf-importance-lists"
      id="toc-get-sorted-total-tf-importance-lists">Get Sorted Total TF
      Importance Lists</a>
- <a href="#final-result-inferred-psc-cocktail"
  id="toc-final-result-inferred-psc-cocktail">Final Result: Inferred PSC
  Cocktail</a>
  - <a
    href="#upon-inspection-of-psc_group_lassoed_tf_list-we-see-that-the-set-of-3-consists-of-xbpd-dril12-and-foxg"
    id="toc-upon-inspection-of-psc_group_lassoed_tf_list-we-see-that-the-set-of-3-consists-of-xbpd-dril12-and-foxg">Upon
    Inspection of “PSC_group_lassoed_tf_list”, we see that the set of 3
    consists of xBPd, DRIL1/2, and Foxg.</a>
  - <a href="#tf-set-of-size-4-consists-of-these-3-tfs-along-with-isl"
    id="toc-tf-set-of-size-4-consists-of-these-3-tfs-along-with-isl">TF Set
    of size 4 consists of these 3 TFs along with Isl.</a>
- <a href="#extended-analyses" id="toc-extended-analyses">Extended
  Analyses</a>
  - <a
    href="#1-comparison-to-gaussian-multivariate-multi-response-linear-regression-gmml"
    id="toc-1-comparison-to-gaussian-multivariate-multi-response-linear-regression-gmml">(1)
    Comparison to Gaussian Multivariate Multi-Response Linear Regression
    (GMML)</a>
    - <a href="#matrices-x-and-y" id="toc-matrices-x-and-y">Matrices X and
      Y</a>
    - <a href="#perform-regression" id="toc-perform-regression">Perform
      Regression</a>
      - <a href="#without-penalty" id="toc-without-penalty">Without penalty</a>
      - <a href="#with-penalty-selecting-for-3-tf-variables"
        id="toc-with-penalty-selecting-for-3-tf-variables">With penalty,
        selecting for 3 TF Variables</a>
    - <a href="#strategy-1" id="toc-strategy-1">Strategy 1</a>
      - <a href="#get-regression-coefficients-for-all-candidate-cocktail-tfs"
        id="toc-get-regression-coefficients-for-all-candidate-cocktail-tfs">Get
        Regression Coefficients for all Candidate Cocktail TFs</a>
      - <a href="#visualize-distribution-of-coefficients-and-get-their-averages"
        id="toc-visualize-distribution-of-coefficients-and-get-their-averages">Visualize
        Distribution of Coefficients and get their Averages</a>
      - <a href="#get-total-importance-scores-for-candidate-cocktail-tfs"
        id="toc-get-total-importance-scores-for-candidate-cocktail-tfs">Get
        Total Importance Scores for Candidate Cocktail TFs</a>
    - <a href="#strategy-2" id="toc-strategy-2">Strategy 2</a>
      - <a
        href="#first-get-the-ranks-of-the-absolute-values-of-the-regression-coefficients"
        id="toc-first-get-the-ranks-of-the-absolute-values-of-the-regression-coefficients">First
        Get the ranks of the absolute values of the regression coefficients</a>
      - <a
        href="#visualize-distribution-of-coefficients-and-get-their-averages-1"
        id="toc-visualize-distribution-of-coefficients-and-get-their-averages-1">Visualize
        Distribution of Coefficients and get their Averages</a>
      - <a href="#get-total-importance-scores-for-candidate-cocktail-tfs-1"
        id="toc-get-total-importance-scores-for-candidate-cocktail-tfs-1">Get
        Total Importance Scores for Candidate Cocktail TFs</a>
    - <a href="#notes-on-final-result" id="toc-notes-on-final-result">Notes on
      Final Result</a>
  - <a href="#2-comparison-with-grnboost2-inferences"
    id="toc-2-comparison-with-grnboost2-inferences">(2) Comparison with
    GRNBoost2 Inferences</a>
    - <a href="#psc-precursor-deg-counts-matrix"
      id="toc-psc-precursor-deg-counts-matrix">PSC Precursor DEG Counts
      Matrix</a>
    - <a href="#inputs-to-grnboost2-pipeline"
      id="toc-inputs-to-grnboost2-pipeline">Inputs to GRNBoost2 Pipeline</a>
    - <a href="#do-grnboost2-50-times" id="toc-do-grnboost2-50-times">Do
      GRNBoost2 50 times.</a>
    - <a
      href="#install-relevant-packages-for-grnboost2_pipelinepy-into-given-r-reticulate-environment"
      id="toc-install-relevant-packages-for-grnboost2_pipelinepy-into-given-r-reticulate-environment">Install
      Relevant Packages for “GRNBoost2_pipeline.py” into given r-reticulate
      environment.</a>
      - <a href="#sorted_total_importance_tfs"
        id="toc-sorted_total_importance_tfs">sorted_total_importance_tfs</a>
      - <a href="#sorted_percent_top_list"
        id="toc-sorted_percent_top_list">sorted_percent_top_list</a>
      - <a href="#inspect-output-contents"
        id="toc-inspect-output-contents">Inspect Output contents.</a>
    - <a href="#notes-on-final-result-1"
      id="toc-notes-on-final-result-1">Notes on Final Result</a>
  - <a href="#3-comparison-with-method-applied-on-terminal-pscs"
    id="toc-3-comparison-with-method-applied-on-terminal-pscs">(3)
    Comparison with Method Applied on Terminal PSCs</a>
    - <a href="#identify-terminal-pscs-and-their-siblings"
      id="toc-identify-terminal-pscs-and-their-siblings">Identify Terminal
      PSCs and their Siblings</a>
    - <a href="#raw-mats" id="toc-raw-mats">Raw Mats</a>
    - <a href="#calculate-terminal-degs-with-respect-to-siblings"
      id="toc-calculate-terminal-degs-with-respect-to-siblings">Calculate
      Terminal DEGs with respect to Siblings</a>
    - <a href="#get-terminal-degs" id="toc-get-terminal-degs">Get Terminal
      DEGs</a>
      - <a href="#write-them-out" id="toc-write-them-out">Write them out.</a>
    - <a href="#get-terminal-deg-tfs" id="toc-get-terminal-deg-tfs">Get
      Terminal DEG TFs</a>
      - <a href="#write-them-out-1" id="toc-write-them-out-1">Write them
        out.</a>
    - <a href="#get-terminal-candidate-cocktail-tfs"
      id="toc-get-terminal-candidate-cocktail-tfs">Get Terminal Candidate
      Cocktail TFs</a>
      - <a href="#write-out-the-candidate-cocktail-tfs-1"
        id="toc-write-out-the-candidate-cocktail-tfs-1">Write out the Candidate
        Cocktail TFs</a>
      - <a href="#brief-aside" id="toc-brief-aside">Brief Aside</a>
    - <a href="#a-do-cocktail-tf-inference-using-group-lasso-model"
      id="toc-a-do-cocktail-tf-inference-using-group-lasso-model">(a) Do
      Cocktail TF Inference using Group Lasso Model</a>
      - <a href="#variables" id="toc-variables">Variables</a>
      - <a href="#make-group-lasso-x-matrix"
        id="toc-make-group-lasso-x-matrix">Make Group Lasso X Matrix</a>
      - <a href="#make-group-lasso-y-matrix"
        id="toc-make-group-lasso-y-matrix">Make Group Lasso Y Matrix</a>
      - <a href="#variables-group-lasso-matrices-x-and-y-and-group-ids-1"
        id="toc-variables-group-lasso-matrices-x-and-y-and-group-ids-1">Variables:
        Group Lasso Matrices X and Y and Group IDs</a>
      - <a href="#run-group-lasso-1" id="toc-run-group-lasso-1">Run Group
        Lasso</a>
    - <a href="#get-sets-of-selected-tf-groups-at-various-values-of-lambda-1"
      id="toc-get-sets-of-selected-tf-groups-at-various-values-of-lambda-1">Get
      Sets of Selected TF Groups at Various Values of Lambda</a>
      - <a
        href="#let-us-visualize-the-number-of-groups-selected-at-various-values-of-lambda-1"
        id="toc-let-us-visualize-the-number-of-groups-selected-at-various-values-of-lambda-1">Let
        us Visualize the Number of Groups Selected at Various Values of
        Lambda</a>
    - <a href="#get-sorted-total-tf-importance-lists-1"
      id="toc-get-sorted-total-tf-importance-lists-1">Get Sorted Total TF
      Importance Lists</a>
      - <a href="#result-we-are-unable-to-uncover-the-psc-cocktail"
        id="toc-result-we-are-unable-to-uncover-the-psc-cocktail">Result: We are
        unable to uncover the PSC Cocktail.</a>
    - <a href="#b-do-gmml" id="toc-b-do-gmml">(b) Do GMML</a>
    - <a href="#perform-regression-1" id="toc-perform-regression-1">Perform
      Regression</a>
      - <a href="#without-penalty-1" id="toc-without-penalty-1">Without
        penalty</a>
      - <a href="#with-penalty-selecting-for-3-tf-variables-1"
        id="toc-with-penalty-selecting-for-3-tf-variables-1">With penalty,
        selecting for 3 TF Variables</a>
    - <a href="#strategy-1-1" id="toc-strategy-1-1">Strategy 1</a>
      - <a href="#get-regression-coefficients-for-all-candidate-cocktail-tfs-1"
        id="toc-get-regression-coefficients-for-all-candidate-cocktail-tfs-1">Get
        Regression Coefficients for all Candidate Cocktail TFs</a>
      - <a
        href="#visualize-distribution-of-coefficients-and-get-their-averages-2"
        id="toc-visualize-distribution-of-coefficients-and-get-their-averages-2">Visualize
        Distribution of Coefficients and get their Averages</a>
      - <a href="#get-total-importance-scores-for-candidate-cocktail-tfs-2"
        id="toc-get-total-importance-scores-for-candidate-cocktail-tfs-2">Get
        Total Importance Scores for Candidate Cocktail TFs</a>
    - <a href="#strategy-2-1" id="toc-strategy-2-1">Strategy 2</a>
      - <a
        href="#first-get-the-ranks-of-the-absolute-values-of-the-regression-coefficients-1"
        id="toc-first-get-the-ranks-of-the-absolute-values-of-the-regression-coefficients-1">First
        Get the ranks of the absolute values of the regression coefficients</a>
      - <a
        href="#visualize-distribution-of-coefficients-and-get-their-averages-3"
        id="toc-visualize-distribution-of-coefficients-and-get-their-averages-3">Visualize
        Distribution of Coefficients and get their Averages</a>
      - <a href="#get-total-importance-scores-for-candidate-cocktail-tfs-3"
        id="toc-get-total-importance-scores-for-candidate-cocktail-tfs-3">Get
        Total Importance Scores for Candidate Cocktail TFs</a>
    - <a href="#c-comparison-with-grnboost2-inferences"
      id="toc-c-comparison-with-grnboost2-inferences">(c) Comparison with
      GRNBoost2 Inferences</a>
      - <a href="#terminal-psc-deg-counts-matrix"
        id="toc-terminal-psc-deg-counts-matrix">Terminal PSC DEG Counts
        Matrix</a>
      - <a href="#inputs-to-grnboost2-pipeline-1"
        id="toc-inputs-to-grnboost2-pipeline-1">Inputs to GRNBoost2 Pipeline</a>
      - <a href="#do-grnboost2-50-times-1" id="toc-do-grnboost2-50-times-1">Do
        GRNBoost2 50 times.</a>
        - <a href="#sorted_total_importance_tfs-1"
          id="toc-sorted_total_importance_tfs-1">sorted_total_importance_tfs</a>
        - <a href="#sorted_percent_top_list-1"
          id="toc-sorted_percent_top_list-1">sorted_percent_top_list</a>
        - <a href="#inspect-output-contents-1"
          id="toc-inspect-output-contents-1">Inspect Output contents.</a>
      - <a href="#note-on-final-result" id="toc-note-on-final-result">Note on
        Final Result</a>
  - <a href="#4-subsampling" id="toc-4-subsampling">(4) Subsampling</a>
    - <a href="#get-50-runs-of-group-lasso-model-on-subsampled-data"
      id="toc-get-50-runs-of-group-lasso-model-on-subsampled-data">Get 50 runs
      of Group Lasso Model on Subsampled Data.</a>
      - <a href="#access-subsamples" id="toc-access-subsamples">Access
        subsamples.</a>
      - <a href="#make-combined_subsampled_tf_df"
        id="toc-make-combined_subsampled_tf_df">Make
        combined_subsampled_tf_df.</a>
      - <a href="#write-out-combined_subsampled_tf_df"
        id="toc-write-out-combined_subsampled_tf_df">Write out
        combined_subsampled_tf_df.</a>
      - <a
        href="#also-write-out-pscs-group_lassoed_tf_list-has-predicted-cocktail"
        id="toc-also-write-out-pscs-group_lassoed_tf_list-has-predicted-cocktail">Also
        write out PSC’s group_lassoed_tf_list (Has Predicted Cocktail).</a>
    - <a href="#run-subsampling-prediction-confidence-pipeline"
      id="toc-run-subsampling-prediction-confidence-pipeline">Run Subsampling
      Prediction Confidence Pipeline.</a>
    - <a href="#note-on-final-result-1" id="toc-note-on-final-result-1">Note
      on Final Result</a>
- <a href="#rsession-info" id="toc-rsession-info">RSession Info</a>

# Brief Introduction

The analysis pipeline present here demonstrates **Circe**, a suite of
computational methods that predict reprogramming Cocktails for neuronal
cell types in Ciona intestinalis, as **applied to the Palp Sensory Cells
(PSCs)**. A **reprogramming Cocktail is a combination of transcription
factors** that when overexpressed, can induce a particular cellular
identity in the cell type in which it was introduced.

**There are 3 Main Computational Methods:**  
**Method 1: Identification of Pseudostates, detailed in:**  
(1) Identify Precursor, Mother, and their Sibling(s) Segment Cell
Populations on Virtual Neural Lineage Tree  
(2) Get PC Coordinates and Corresponding Pseudotimes for PSC Cells  
(3) Get Pseudostates

**Method 2: Calculation of Differentially Expressed Genes, detailed
in:**  
(4) Get Differentially Expressed Genes  
(5) Get Candidate Cocktail TFs  

**3. Group Lasso, detailed in:**  
(6) Do Cocktail TF Inference using Group Lasso Model

**Extended Analyses**  
Finally, a series of extended analyses were performed to assess the
**strength of both Circe’s biological logic** (studying gene expression
levels at the time of specification vs. at the terminally differentiated
state) **and computational results**.

- Extended Analyses 1 and 2 compare Cocktail Inference results from
  another Linear Model, Gaussian Multivariate Multi-Response Linear
  Regression (GMML), and a Nonlinear Model based on GRNBoost2
  respectively.
- Extended Analysis 3 compares Cocktail Inference results with Terminal
  BTNs (fully differentiated BTNs).
- Extended Analysis 4 performs subsampling of BTN Precursor data to
  build confidence for Group Lasso results.

**The final results of TF Sets of sizes 3 and 4 are also presented
here.**  

**Experimental Validation of xBPd, Foxg, and Isl (in TF Set of size 4)
result can be found** in Figure 3-S1 in Chapter 3 Supplementary
Appendix.

# Import Relevant Libraries

``` r
library(URD)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(stats)
library(Matrix)
library(stringr)
library(edgeR)
library(rlang)
library(grpreg)
library(glmnet)
library(dplyr)
library(purrr)
library(stringr)
library(reticulate)
source("./helper_files/final_main_pipeline_start_script.R")
source("./helper_files/final_helper_functions.R")
```

# (1) Identify Precursor, Mother, and their Sibling(s) Segment Cell Populations on Virtual Neural Lineage Tree

Segment cells refers to the cells along the tree segment/branch that
will contain Precursor, Precursor Sibling, and Mother pseudostates.

Note: “CNS.k300s6w4.tree.built.name” is the URD object generated from
time-course single-cell RNA-seq data in Chen et. al\* for the Central
Nervous System (CNS) of Ciona intestinalis. It contains cell names,
their pseudotimes, the Virtual Lineage Tree, and all URD object
attributes.

\*Cao, Chen, et al. “Comprehensive single-cell transcriptome lineages of
a proto-vertebrate.” Nature 571.7765 (2019): 349-354.

## Precursor Segment

### Get IDs of PSC Precursor Segment and its Sibling Segment

``` r
# Get ID from URD Object.
PSC_ID <- "10"
# 56
PSC_sibling_ID <- segSiblings(CNS.k300s6w4.tree.built.name, PSC_ID, include.self = F)
```

### Get PSC Precursor and Sibling Segment Cells

``` r
# 551
PSC_precursor_segment_cells <- cellsInCluster(CNS.k300s6w4.tree.built.name, clustering="segment", PSC_ID)
# 394
PSC_sibling_segment_cells <- cellsInCluster(CNS.k300s6w4.tree.built.name, clustering="segment", PSC_sibling_ID)
```

### Order Cells by Pseudotime

``` r
pseudotimed_PSC_precursor_segment_cells <- get_pseudotimed_cells(cells_in_segments_list, pseudotimes, PSC_ID, length(PSC_precursor_segment_cells), ascending=TRUE)

pseudotimed_PSC_precursor_sibling_segment_cells <- get_pseudotimed_cells(cells_in_segments_list, pseudotimes, PSC_sibling_ID, length(PSC_sibling_segment_cells), ascending=TRUE)
```

## Mother Segment

### Get IDs of PSC Mother Segment and its Siblings

``` r
PSC_mother_ID <- "63"
PSC_mother_sibling_IDs <- segSiblings(CNS.k300s6w4.tree.built.name, PSC_mother_ID, include.self = F)

PSC_mother_sibling_70_ID <- PSC_mother_sibling_IDs[1]
PSC_mother_sibling_73_ID <- PSC_mother_sibling_IDs[2]
```

### Get PSC Mother and Sibling Segment Cells

``` r
# 853
PSC_mother_segment_cells <- cellsInCluster(CNS.k300s6w4.tree.built.name, clustering="segment", PSC_mother_ID)
# 463
PSC_mother_sibling_70_segment_cells <- cellsInCluster(CNS.k300s6w4.tree.built.name, clustering="segment", PSC_mother_sibling_70_ID)
# 25
PSC_mother_sibling_73_segment_cells <- cellsInCluster(CNS.k300s6w4.tree.built.name, clustering="segment", PSC_mother_sibling_73_ID)
```

### Order Cells by Pseudotime

``` r
pseudotimed_PSC_mother_segment_cells <- get_pseudotimed_cells(cells_in_segments_list, pseudotimes, PSC_mother_ID, length(PSC_mother_segment_cells), ascending=FALSE)

pseudotimed_PSC_mother_sibling_70_segment_cells <- get_pseudotimed_cells(cells_in_segments_list, pseudotimes, PSC_mother_sibling_70_ID, length(PSC_mother_sibling_70_segment_cells), ascending=FALSE)

pseudotimed_PSC_mother_sibling_73_segment_cells <- get_pseudotimed_cells(cells_in_segments_list, pseudotimes, PSC_mother_sibling_73_ID, length(PSC_mother_sibling_73_segment_cells), ascending=FALSE)
```

## Visualize Cells on Virtual Neural Lineage Tree

``` r
PSC_segment_cells_list <- list("pseudotimed_PSC_precursor_segment_cells" = pseudotimed_PSC_precursor_segment_cells,
                               "pseudotimed_PSC_precursor_segment_cells" = pseudotimed_PSC_precursor_segment_cells,
                               "pseudotimed_PSC_mother_segment_cells" = pseudotimed_PSC_mother_segment_cells,
                               "pseudotimed_PSC_mother_sibling_70_segment_cells" = pseudotimed_PSC_mother_sibling_70_segment_cells,
                               "pseudotimed_PSC_mother_sibling_73_segment_cells" = pseudotimed_PSC_mother_sibling_73_segment_cells)

modified_URD <- add_segment_state(PSC_segment_cells_list, "PSC_segment_cells_list")
plotTree(modified_URD, label.segments = F, label = "segment_state", label.type = "group")
```

# (2) Get PC Coordinates and Corresponding Pseudotimes for PSC Cells

The PC coordinates and corresponding pseudotimes for PSC Cells will help
determine where to cut the segments and find pseudostates for PSC
Precursors, their Sibling(s) and Mothers. Pseudostates are
approximations of true Precursor, Precursor Sibling(s) and Mother
states. They will be used to calculate Differentially Expressed Genes in
the next part of the pipeline.

## Make Seurat Object

``` r
PSC_cells <- c(pseudotimed_PSC_precursor_segment_cells, pseudotimed_PSC_precursor_sibling_segment_cells, 
               pseudotimed_PSC_mother_segment_cells, pseudotimed_PSC_mother_sibling_70_segment_cells, 
               pseudotimed_PSC_mother_sibling_73_segment_cells)
  
PSC_raw_mat <- URD_raw_data[, PSC_cells]
PSC_seurat <- CreateSeuratObject(counts = PSC_raw_mat,
                                min.cells = 3,
                                min.features = 200)
# Normalize Data.
PSC_seurat <- NormalizeData(PSC_seurat)
# Find Variable Feature.
PSC_seurat <- FindVariableFeatures(PSC_seurat, selection.method = "vst", nfeatures = 2000)
# Scale data.
all_PSC_seurat_genes <- rownames(PSC_seurat)
# use.umi argument to regress on umi count data.
PSC_seurat <- ScaleData(PSC_seurat, features = all_PSC_seurat_genes, use.umi = TRUE)
# Do PCA.
PSC_seurat <- RunPCA(PSC_seurat, features = VariableFeatures(object = PSC_seurat))
```

## Use Elbow Plot to choose the number of PCs (12)

``` r
ElbowPlot(PSC_seurat)
```

## Key Variables

``` r
PSC_PCA_COORDINATES <- Embeddings(PSC_seurat, reduction = "pca")
PSC_NUM_PCS <- 12
```

# (3) Get Pseudostates

As described, pseudostates are approximations of true Precursor,
Precursor Sibling, and Mother populations. We will first make
“PC_pseudotime_dfs”, which are data frames that contain PC coordinates
and pseudotimes for the PSC cells.

The data frames of Precursor Precursor Sibling, and Mother segments will
exported as csv files, which will then be used as inputs to modules in
“get_pseudostates.py”.

PC_pseudotime_dfs of Mother’s Siblings will be subjected to a separate
function (get_est_mother_sibling_cells) to get their pseudostates.

## Make PC_pseudotime_dfs

``` r
PSC_precursor_segment_PC_pseudotime_df <- make_PC_pseudotime_df(pseudotimed_PSC_precursor_segment_cells, PSC_PCA_COORDINATES, PSC_NUM_PCS, pseudotimes)

PSC_sibling_segment_PC_pseudotime_df <- make_PC_pseudotime_df(pseudotimed_PSC_precursor_sibling_segment_cells, PSC_PCA_COORDINATES, PSC_NUM_PCS, pseudotimes)

PSC_mother_segment_PC_pseudotime_df <- make_PC_pseudotime_df(pseudotimed_PSC_mother_segment_cells, PSC_PCA_COORDINATES, PSC_NUM_PCS, pseudotimes)
```

## Export PC_pseudotime_dfs

``` r
PC_pseudotime_dfs_out_dir <- "./PC_pseudotime_dfs/"

write.table(PSC_precursor_segment_PC_pseudotime_df, file = paste0(PC_pseudotime_dfs_out_dir, "PSC_precursor_segment_PC_pseudotime_df.csv"), sep = "\t", row.names = FALSE)

write.table(PSC_sibling_segment_PC_pseudotime_df, file = paste0(PC_pseudotime_dfs_out_dir, "PSC_sibling_segment_PC_pseudotime_df.csv"), sep = "\t", row.names = FALSE)

write.table(PSC_mother_segment_PC_pseudotime_df, file = paste0(PC_pseudotime_dfs_out_dir, "PSC_mother_segment_PC_pseudotime_df.csv"), sep = "\t", row.names = FALSE)
```

## Load “get_pseudostates.py”

“get_pseudostates.py” consists of modules that, given a
PC_pseudotime_dict of cells on a segment, finds a cutoff that demarcates
the pseudostate on that segment. That is, all the cells on a segment
until that cutoff are considered a part of the pseudostate.

### Create a new reticulate environment

This is where we will install all necessary python packages. use

``` r
# create a new environment 
conda_create("r-reticulate")
use_condaenv("r-reticulate")
```

### Install Relevant Packages for “get_pseudostates.py” into given r-reticulate environment.

``` r
py_install("python-igraph", pip = T)
conda_install("r-reticulate", "scipy")
```

``` r
source_python("get_pseudostates.py")
```

## Key Variables

> n_leiden_iterations = number of iterations to be used by the Leiden
> algorithm.

> n_routine_iterations = number of iterations for the entire routine.

> verbose = refers to whether cutoffs will be printed as routine runs or
> not.

Please refer to Code documentation for further details.

``` r
n_leiden_iterations <- 5L
n_routine_iterations <- 50L
verbose <- FALSE
```

## PSC Precursor

``` r
set.seed(0)
# Took ~11 mins to run.
PSC_precursor_segment_df_file <- "./PC_pseudotime_dfs/PSC_precursor_segment_PC_pseudotime_df.csv"

PSC_precursor_segment_cutoff_output = get_pseudostates_pipeline(PSC_precursor_segment_df_file, n_leiden_iterations, n_routine_iterations, verbose)
```

``` r
set.seed(0)
PSC_precursor_pseudostate_cutoff <- PSC_precursor_segment_cutoff_output[[2]]
# 339
PSC_precursor_pseudostate_cells <- pseudotimed_PSC_precursor_segment_cells[1:PSC_precursor_pseudostate_cutoff]
```

## PSC Precursor Sibling

``` r
set.seed(0)
# Took ~4 mins to run.
PSC_precursor_sibing_segment_df_file <- "./PC_pseudotime_dfs/PSC_sibling_segment_PC_pseudotime_df.csv"
PSC_precursor_segment_sibling_cutoff_output <- get_pseudostates_pipeline(PSC_precursor_sibing_segment_df_file, n_leiden_iterations, n_routine_iterations, verbose)
```

``` r
set.seed(0)
PSC_precursor_sibling_pseudostate_cutoff <- PSC_precursor_segment_sibling_cutoff_output[[2]]
# 180
PSC_precursor_sibling_pseudostate_cells <- pseudotimed_PSC_precursor_sibling_segment_cells[1:PSC_precursor_sibling_pseudostate_cutoff]
```

## PSC Mother

``` r
set.seed(0)
# Took ~30 mins to run.
PSC_mother_segment_df_file <- "./PC_pseudotime_dfs/PSC_mother_segment_PC_pseudotime_df.csv"
PSC_mother_segment_cutoff_output <- get_pseudostates_pipeline(PSC_mother_segment_df_file, n_leiden_iterations, n_routine_iterations, verbose)
```

``` r
set.seed(0)
PSC_mother_pseudostate_cutoff <- PSC_mother_segment_cutoff_output[[2]]
# 353
PSC_mother_pseudostate_cells <- pseudotimed_PSC_mother_segment_cells[1:PSC_mother_pseudostate_cutoff]
```

## Mother’s Siblings

The Mother’s Siblings must be comparable to it in pseudotime. Thus, we
use the “get_mother_sibling_pseudostates” function to find their
pseudostates. We combine them into one pseudostate.

``` r
# 41
PSC_mother_sibling_70_output <- get_mother_sibling_pseudostates(PSC_mother_pseudostate_cells, pseudotimed_PSC_mother_sibling_70_segment_cells)

mother_sibling_70_first_index <- PSC_mother_sibling_70_output$first_index
mother_sibling_70_last_index <- PSC_mother_sibling_70_output$last_index

# 0
PSC_mother_sibling_73_output <- get_mother_sibling_pseudostates(PSC_mother_pseudostate_cells, pseudotimed_PSC_mother_sibling_73_segment_cells)

mother_sibling_73_first_index <- PSC_mother_sibling_73_output$first_index
mother_sibling_73_last_index <- PSC_mother_sibling_73_output$last_index

# 41
PSC_mother_sibling_pseudostate_cells <- c(pseudotimed_PSC_mother_sibling_70_segment_cells[mother_sibling_70_first_index:mother_sibling_70_last_index],
                                          pseudotimed_PSC_mother_sibling_73_segment_cells[mother_sibling_73_first_index:mother_sibling_73_last_index])
```

## Visualize Pseudostates on Virtual Neural Lineage Tree

``` r
PSC_pseudostates_list <- list("PSC_precursor_pseudostate_cells" = PSC_precursor_pseudostate_cells,
                              "PSC_precursor_sibling_pseudostate_cells" = PSC_precursor_sibling_pseudostate_cells,
                              "PSC_mother_pseudostate_cells" = PSC_mother_pseudostate_cells,
                              "PSC_mother_sibling_pseudostate_cells" = PSC_mother_sibling_pseudostate_cells)

modified_URD <- add_segment_state(PSC_pseudostates_list, "PSC_pseudostates_list")
plotTree(modified_URD, label.segments = F, label = "segment_state", label.type = "group")
```

# (4) Get Differentially Expressed Genes

After identifying the Precursor, Mother, and their respective Sibling(s)
populations, we obtain the following gene lists by calculating
Differentially Expressed Genes between the specified populations.

### Precursor Identity Genes

> These are genes that define the overall state of the Precursor. They
> will be activated by the Cocktail TFs.

Comparison 1: Precursor vs. Sibling

### Precursor DEG TFs

> This set contains TFs relevant for specification that are present in
> the Precursor.

Comparison 2: Precursor vs. Sibling(s)

Comparison 3: Precursor vs. Mother

> The DEG TFs that result from the comparison between the Precursor and
> its Mother are those that get activated in the Precursor from the
> previous time point.

### Mother DEG TFs

> This set contains TFs relevant for specification that are present in
> the Mother. Comparison 4: Mother vs. Sibling(s)

Finally, we will take the union of DEG TFs from Comparisons 2, 3, and 4
and refer to it as “Unfiltered Candidate Cocktail TFs”.

## Raw Matrices

``` r
# Raw Matrices with expression levels of only TFs. 
# 339
PSC_precursor_pseudostate_TF_raw_mat <- URD_raw_data[tfs_in_URD, PSC_precursor_pseudostate_cells]
# 180
PSC_precursor_sibling_pseudostate_TF_raw_mat <- URD_raw_data[tfs_in_URD, PSC_precursor_sibling_pseudostate_cells]
# 353
PSC_mother_pseudostate_TF_raw_mat <- URD_raw_data[tfs_in_URD, PSC_mother_pseudostate_cells]
# 41
PSC_mother_sibling_pseudostate_TF_raw_mat <- URD_raw_data[tfs_in_URD, PSC_mother_sibling_pseudostate_cells]

# Raw Matrices with expression levels of all genes.
PSC_precursor_pseudostate_raw_mat <- URD_raw_data[, PSC_precursor_pseudostate_cells]

PSC_precursor_sibling_pseudostate_raw_mat <- URD_raw_data[, PSC_precursor_sibling_pseudostate_cells]
```

All DEG lists will be written out to this dir.

``` r
DEG_lists_out_dir <- "./DEG_lists/"
```

## Precursor vs. Precursor Sibling- DEG TFs and Precursor Identity Genes

### Precursor Identity Genes

Perform edgeR to get DEGs between Precursor and Precursor Sibling
Pseudostates.

``` r
PSC_precursor_v_sibling_pseudostate_mat <- as.matrix(cbind(PSC_precursor_pseudostate_raw_mat, PSC_precursor_sibling_pseudostate_raw_mat))
```

``` r
PSC_precursor_v_sibling_pseudostate_edgeR_outputs <- complete_edgeR_pipeline(PSC_precursor_v_sibling_pseudostate_mat, ncol(PSC_precursor_pseudostate_raw_mat), ncol(PSC_precursor_sibling_pseudostate_raw_mat), 0.05, -1, 50)
```

Get the DEGs.

``` r
# 84
PSC_precursor_pseudostate_DEGs_df <- PSC_precursor_v_sibling_pseudostate_edgeR_outputs[[1]]
PSC_precursor_pseudostate_DEGs <- rownames(PSC_precursor_pseudostate_DEGs_df)
```

#### Write out Precursor DEGs (Precursor Identity Genes)

``` r
write_csv_helper(PSC_precursor_pseudostate_DEGs, paste0(DEG_lists_out_dir, "PSC_precursor_pseudostate_DEGs.csv"))
```

### Precursor DEG TFs

``` r
PSC_precursor_v_sibling_pseudostate_TF_mat <- as.matrix(cbind(PSC_precursor_pseudostate_TF_raw_mat, PSC_precursor_sibling_pseudostate_TF_raw_mat))

PSC_precursor_v_sibling_pseudostate_TF_edgeR_outputs <- complete_edgeR_pipeline(PSC_precursor_v_sibling_pseudostate_TF_mat, ncol(PSC_precursor_pseudostate_TF_raw_mat), ncol(PSC_precursor_sibling_pseudostate_TF_raw_mat), 0.05, -1, 50)
```

``` r
# 7
PSC_precursor_v_sibling_pseudostate_TF_DEGs_df <- PSC_precursor_v_sibling_pseudostate_TF_edgeR_outputs[[1]]
# Get TF orthologs.
named_PSC_precursor_v_sibling_pseudostate_TF_DEGs_df <- convert_khids_to_ortholog(PSC_precursor_v_sibling_pseudostate_TF_DEGs_df)

PSC_precursor_v_sibling_pseudostate_TF_DEGs <- named_PSC_precursor_v_sibling_pseudostate_TF_DEGs_df$Ghost.Name
```

#### Write out TF DEGs

``` r
write_csv_helper(PSC_precursor_v_sibling_pseudostate_TF_DEGs, paste0(DEG_lists_out_dir, "PSC_precursor_v_sibling_pseudostate_TF_DEGs.csv"))
```

## Precursor vs. Mother- DEG TFs

``` r
PSC_precursor_v_mother_pseudostate_TF_mat <- as.matrix(cbind(PSC_precursor_pseudostate_TF_raw_mat, PSC_mother_pseudostate_TF_raw_mat))

PSC_precursor_v_mother_pseudostate_TF_edgeR_outputs <- complete_edgeR_pipeline(PSC_precursor_v_mother_pseudostate_TF_mat, ncol(PSC_precursor_pseudostate_TF_raw_mat), ncol(PSC_mother_pseudostate_TF_raw_mat), 0.05, -1, 50)
```

``` r
# 13
PSC_precursor_v_mother_pseudostate_TF_DEGs_df <- PSC_precursor_v_mother_pseudostate_TF_edgeR_outputs[[1]]
# Get TF orthologs.
named_PSC_precursor_v_mother_pseudostate_TF_DEGs_df <- convert_khids_to_ortholog(PSC_precursor_v_mother_pseudostate_TF_DEGs_df)

PSC_precursor_v_mother_pseudostate_TF_DEGs <- named_PSC_precursor_v_mother_pseudostate_TF_DEGs_df$Ghost.Name
```

#### Write out TF DEGs

``` r
write_csv_helper(PSC_precursor_v_mother_pseudostate_TF_DEGs, paste0(DEG_lists_out_dir, "PSC_precursor_v_mother_pseudostate_TF_DEGs.csv"))
```

## Mother vs. Mother’s Sibling- DEG TFs

``` r
PSC_mother_v_sibling_pseudostate_TF_mat <- as.matrix(cbind(PSC_mother_pseudostate_TF_raw_mat, PSC_mother_sibling_pseudostate_TF_raw_mat))

PSC_mother_v_sibling_pseudostate_TF_edgeR_outputs <- complete_edgeR_pipeline(PSC_mother_v_sibling_pseudostate_TF_mat, ncol(PSC_mother_pseudostate_TF_raw_mat), ncol(PSC_mother_sibling_pseudostate_TF_raw_mat), 0.05, -1, 50)
```

``` r
# 9
PSC_mother_v_sibling_pseudostate_TF_DEGs_df <- PSC_mother_v_sibling_pseudostate_TF_edgeR_outputs[[1]]
# Get TF orthologs.
named_PSC_mother_v_sibling_pseudostate_TF_DEGs_df <- convert_khids_to_ortholog(PSC_mother_v_sibling_pseudostate_TF_DEGs_df)

PSC_mother_v_sibling_pseudostate_TF_DEGs <- named_PSC_mother_v_sibling_pseudostate_TF_DEGs_df$Ghost.Name
```

#### Write out TF DEGs

``` r
write_csv_helper(PSC_mother_v_sibling_pseudostate_TF_DEGs, paste0(DEG_lists_out_dir, "PSC_mother_v_sibling_pseudostate_TF_DEGs.csv"))
```

# (5) Get Candidate Cocktail TFs

We assume that the most relevant Unfiltered Candidate Cocktail TFs in
either the Precursor or Mother will be those that are most highly
expressed. Furthermore, based on the idea that TFs that are highly
expressed in the Mother are likely to continue to be highly expressed in
the Precursor, we will take the union of the top 50% most highly
expressed Unfiltered Candidate Cocktail TFs in either the Mother or
Precursor, and consider these the Candidate Cocktail TFs.

These Candidate Cocktail TFs, along with the Precursor Identity Genes
will be used as inputs to the Group Lasso Model.

## Get Mean Expression of all DEG TFs

``` r
PSC_precursor_pseudostate_mean_TF_exp_df <- make_named_mean_TF_exp_df(PSC_precursor_pseudostate_cells)

PSC_mother_pseudostate_mean_TF_exp_df <- make_named_mean_TF_exp_df(PSC_mother_pseudostate_cells)
```

######## REF SUPP FIGURE PG 17

## Take the union of top 50% most highly expressed Unfiltered Candidate Cocktail TFs in Precursor and Mother. These are the Candidate Cocktail TFs

``` r
# 19
PSC_unfiltered_candidate_cocktail_tfs <- unique(c(named_PSC_precursor_v_sibling_pseudostate_TF_DEGs_df$KH.gene.model,
                                                  named_PSC_precursor_v_mother_pseudostate_TF_DEGs_df$KH.gene.model,
                                                  named_PSC_mother_v_sibling_pseudostate_TF_DEGs_df$KH.gene.model))

# Take 50% cutoff.
# This returns both a union and intersection of top 50% most highly expressed unfiltered candidate cocktail TFs.
PSC_cutoff_df_list <- get_cutoff_df(PSC_precursor_pseudostate_mean_TF_exp_df, PSC_mother_pseudostate_mean_TF_exp_df, PSC_unfiltered_candidate_cocktail_tfs, 0.5)

# Let us take the union.
# 13
PSC_cutoff_union <- PSC_cutoff_df_list$union
# These are the Candidate Cocktail TFs.
PSC_candidate_cocktail_tfs <- unique(c(PSC_cutoff_union$Ghost.Name_precursor, PSC_cutoff_union$Ghost.Name_mother))
# Remove any NA's!
PSC_candidate_cocktail_tfs <- PSC_candidate_cocktail_tfs[!is.na(PSC_candidate_cocktail_tfs)]

# We also want their khids.
PSC_candidate_cocktail_khids <- PSC_cutoff_union$KH.gene.model
```

### Write out the Candidate Cocktail TFs

``` r
write_csv_helper(PSC_candidate_cocktail_tfs, paste0(DEG_lists_out_dir, "PSC_candidate_cocktail_tfs.csv"))
```

# (6) Do Cocktail TF Inference using Group Lasso Model

The final component of the computational method is a Group Lasso Model
in which individual TFs are groups and subsets of these are selected at
various values of the regularization parameter lambda to “explain”
expression of Precursor DEGs. The model captures the idea that
expression of Precursor DEGs will be regulated by a common set of TFs.
We consider the first set of 3 TFs as the Cocktail.

Please refer to the Methods of paper for a full description of the
model.

Briefly, Group Lasso Matrix X contains expression of Candidate Cocktail
TFs in all Precursor cells; Matrix Y contains expression of Precursor
DEGs in all Precursor cells; The Group Lasso Model learns a coefficient
vector that has “DEG-specific” TF coefficients, that is, it learns a
unique set of TF weights for each Precursor DEG as each DEG will be
affected by the same TF to a different extent.

The Group Lasso Model is implemented using the grpreg package.

## Get a set of Precursor DEGs excluding DEG TFs

Their expression will be explained by the model.

``` r
# 6
common_PSC_DEG_TFs <- intersect(PSC_candidate_cocktail_khids, PSC_precursor_pseudostate_DEGs)
# 78
cleaned_PSC_precursor_pseudostate_DEGs <- setdiff(PSC_precursor_pseudostate_DEGs, common_PSC_DEG_TFs)
```

## Get numbers of Precursor cells (n), Candidate Cocktail TFs (p), and Precusor DEGs (k)

``` r
# 339, 13, 78
PSC_n <- length(PSC_precursor_pseudostate_cells)
PSC_p <- length(PSC_candidate_cocktail_tfs)
PSC_k <- length(cleaned_PSC_precursor_pseudostate_DEGs)
```

## Make Group Lasso Matrix X

Matrix S is a matrix that has the expression of all Candidate Cocktail
TFs in all Precursor Cells of dimensions n x p (below). It will be
repeated k times along the diagonal of Group Lasso Matrix X. Thus, Group
Lasso Matrix X will have dimensions nk x pk. All other entries in Group
Lasso Matrix X will be 0.

### Make Matrix S

Matrix S is the same as PSC_TFs_in_precursor_pseudostate_mat below.

``` r
# 339 x 13
PSC_TFs_in_precursor_pseudostate_mat <- t(URD_data[PSC_candidate_cocktail_khids, PSC_precursor_pseudostate_cells])
```

### Concatenate k S Matrices to make Group Lasso Matrix X

``` r
PSC_S_matrix <- PSC_TFs_in_precursor_pseudostate_mat
PSC_zero_matrix <- matrix(0, PSC_n, PSC_p)
# This is a dummy matrix- we will remove it at the end.
PSC_grplasso_X_mat <- matrix(1, PSC_n, PSC_p*PSC_k)

for (i in 1:PSC_k) {
  PSC_L_matrix <- do.call(cbind, replicate(i-1, PSC_zero_matrix, simplify=FALSE))
  PSC_R_matrix <- do.call(cbind, replicate(PSC_k-i, PSC_zero_matrix, simplify=FALSE))
  # At each iteration, we create a row in which Matrix S is sandwiched between 
  # Left "L" Zero Matrix(ces) and Right "R" Zero Matrix(ces) given its position on
  # the diagonal of Group Lasso Matrix X.
  PSC_grplasso_mat_row <- do.call(cbind, list(PSC_L_matrix, PSC_S_matrix, PSC_R_matrix))
  PSC_grplasso_X_mat <- rbind(PSC_grplasso_X_mat, PSC_grplasso_mat_row)
}

# 26442 x 1014
PSC_grplasso_X_mat_2 <- PSC_grplasso_X_mat[-c(1:PSC_n), ]
```

### Rename Row and Column names for Convenience

#### Rename Row Names to reflect cell number and its DEG status

The Row Names have the following pattern: the first n rows are
DEG_1\_cell_1… DEG_1\_cell_n. The following n rows are DEG_2\_cell_1…
DEG_2\_cell_n for all k DEGs.

``` r
PSC_grplasso_X_mat_row_cells <- rep(paste0("_cell_", seq(1, PSC_n)), PSC_k)
PSC_DEG_statuses <- rep(c(1:PSC_k), each = PSC_n)
PSC_new_grplasso_X_mat_2_rownames <- paste0("DEG_", PSC_DEG_statuses, PSC_grplasso_X_mat_row_cells)

rownames(PSC_grplasso_X_mat_2) <- PSC_new_grplasso_X_mat_2_rownames
```

#### Rename Column Names to reflect TF number and its DEG status

Similarly, the Column Names have the following pattern: the first p rows
are DEG_1\_TF1… DEG_1\_TFp. The following p rows are DEG_2\_TF1…
DEG_2\_TFp for all p TFs.

``` r
PSC_grplasso_X_mat_col_TFs <- rep(paste0("_TF", seq(1, PSC_p)), PSC_k)
PSC_DEG_statuses_cols <- rep(c(1:PSC_k), each = PSC_p)
PSC_new_grplasso_X_mat_2_colnames <- paste0("DEG_", PSC_DEG_statuses_cols, PSC_grplasso_X_mat_col_TFs)

colnames(PSC_grplasso_X_mat_2) <- PSC_new_grplasso_X_mat_2_colnames
```

### Final Group Lasso Matrix X

``` r
# 26442 x 1014
final_grplasso_X_mat <- as.matrix(PSC_grplasso_X_mat_2)
```

## Make Group Lasso Matrix Y

Group Lasso Matrix Y contains the expression of all Precursor DEGs in
all Precursor Cells in the following way: The first n rows contain
expression of DEG 1 in all n cells; the next n rows contain expression
of DEG 2 in all n cells; so and so forth for k DEGs. It has dimensions
nk x 1.

### Make a matrix that has the expression of all Precursor DEGs in Precursor Cells of dimensions k x n

``` r
# 339 x 78
PSC_precursor_DEGs_in_precursors_mat <- t(URD_data[cleaned_PSC_precursor_pseudostate_DEGs, PSC_precursor_pseudostate_cells])
```

### Index above matrix so as so make Group Lasso Matrix Y

``` r
PSC_grplasso_Y_vec <- c()

for (DEG in 1:PSC_k) {
  for (cell in 1:PSC_n) {
    PSC_DEG_cell_exp <- PSC_precursor_DEGs_in_precursors_mat[cell, DEG]
    PSC_grplasso_Y_vec <- c(PSC_grplasso_Y_vec, PSC_DEG_cell_exp)
  }
}

# 26442 x 1
PSC_grplasso_Y_mat <- as.matrix(PSC_grplasso_Y_vec, ncol=1, nrow=(PSC_n*PSC_k))
```

### Rename Row names to reflect cell and DEG status

The Row Names follow the same pattern as the Row Names of Group Lasso
Matrix X.

``` r
PSC_grplasso_Y_DEG_statuses <- rep(c(1:PSC_k), each = PSC_n)
PSC_grplasso_Y_DEGs <- paste0("DEG_", PSC_grplasso_Y_DEG_statuses)
PSC_grplasso_Y_cells <- rep(paste0("_cell_", seq(1, PSC_n)), PSC_k)
PSC_new_grplasso_Y_rownames <- paste0(PSC_grplasso_Y_DEGs, PSC_grplasso_Y_cells)

rownames(PSC_grplasso_Y_mat) <- PSC_new_grplasso_Y_rownames
```

### Final Matrix Y

``` r
# 26442 x 1
final_PSC_grplasso_Y_mat <- PSC_grplasso_Y_mat
```

## Run Group Lasso

### Variables: Group Lasso Matrices X and Y and Group IDs

Each TF is a group. Thus, there are p groups each of size k. Given the
way Group Lasso Matrix X is organized, the Group IDs follow the pattern
of 1, 2, … p; 1, 2, … p; k times.

``` r
grplasso_PSC_groups <- rep(seq(1, PSC_p), PSC_k)

PSC_X <- final_grplasso_X_mat
PSC_Y <- final_PSC_grplasso_Y_mat
PSC_group_id <- grplasso_PSC_groups
```

### Do Cross-Validation

``` r
# Took ~6 mins to run.
grplasso_PSC_cv_grepreg_fit <- cv.grpreg(PSC_X, PSC_Y, PSC_group_id, penalty="grLasso")
```

### Get Sets of Selected TF Groups at Various Values of Lambda

``` r
# Key: lambda; Values: Set of Selected TF Groups
PSC_group_lassoed_tf_list <- get_group_lassoed_tf_list(grplasso_PSC_cv_grepreg_fit, PSC_candidate_cocktail_tfs)

# Key: lambda; Values: Set of Non-selected TF Groups
PSC_non_group_lassoed_TF_list <- lapply(PSC_group_lassoed_tf_list, get_non_group_lassoed_TFs, PSC_candidate_cocktail_tfs)
```

#### Let us Visualize the Number of Groups Selected at Various Values of Lambda

``` r
plot(grplasso_PSC_cv_grepreg_fit)
```

### Get Sorted Total TF Importance Lists

The Total Importance Score of a Cocktail TF is the sum of the absolute
value of its DEG-specific coefficients. It thus gives a quantitative
measure of a Cocktail TF’s strength over the Gene Regulatory Network (as
defined by the Precursor DEGs) underlying the Precursor’s state.

``` r
# (For all lambdas) key = lambda; value is a list with key = TF and value = non0 betas.
PSC_TF_importances_lists <- get_TF_importances_list(grplasso_PSC_cv_grepreg_fit)

# Sum the absolute values of these betas for all TFs for all lambdas.
PSC_sum_total_TF_importance_lists <- PSC_TF_importances_lists %>% map(get_total_summed_TF_importance_list)

# Make sorted dfs given these Total Importance Scores.
PSC_sorted_sum_total_TF_importance_lists <- PSC_sum_total_TF_importance_lists %>% map(make_sorted_TF_importance_df, PSC_candidate_cocktail_tfs)
```

# Final Result: Inferred PSC Cocktail

#### Upon Inspection of “PSC_group_lassoed_tf_list”, we see that the set of 3 consists of xBPd, DRIL1/2, and Foxg.

#### TF Set of size 4 consists of these 3 TFs along with Isl.

# Extended Analyses

> Extended Analyses 1 and 2 compare Cocktail Inference results from
> another Linear Model, Gaussian Multivariate Multi-Response Linear
> Regression (GMML), and a Nonlinear Model based on GRNBoost2
> respectively.

> Extended Analysis 3 compares Cocktail Inference results with Terminal
> BTNs (fully differentiated BTNs).

> Extended Analysis 4 performs subsampling of BTN Precursor data to
> build confidence for Group Lasso results.

## (1) Comparison to Gaussian Multivariate Multi-Response Linear Regression (GMML)

DEG expression levels can be modeled as linear combinations of Candidate
Cocktail TF expression levels. In our case, we have multiple predictors
(TFs) and multiple response variables (expression of DEGs) and both the
predictors and response variables are correlated. Thus, an appropriate
way to model these interactions is via a gaussian multivariate
multi-response linear regression (GMML).

To infer the Cocktail, we employ 2 strategies:

Strategy 1: We calculate the average coefficient values for each
Candidate Cocktail TF, get a Total Importance Score (described
afterwards), and consider the top 3 based on the Total Importance Score
as the Cocktail.

Strategy 2: We repeat the same process, but on the **ranks** of the
coefficient values for each Candidate Cocktail TF. This is because DEGs
with very high expression levels will bias coefficient values to be
higher. We prefer Strategy 2.

Refer to Methods of paper for a full description of the model.

### Matrices X and Y

Matrix X has expression of Cocktail Candidate TFs in Precursor cells and
has dimensions n x p.  Matrix Y has expression of all Precursor DEGs in
Precursor cells and has dimensions n x k.

``` r
glm_PSC_X_mat <- as.matrix(PSC_TFs_in_precursor_pseudostate_mat)
glm_PSC_Y_mat <- as.matrix(PSC_precursor_DEGs_in_precursors_mat)
```

### Perform Regression

#### Without penalty

We will refer to this as the complete GMML model.

``` r
# ~Took 5s to run.
PSC_fit <- cv.glmnet(glm_PSC_X_mat, glm_PSC_Y_mat, intercept = FALSE, family = "mgaussian")
PSC_lambda <- PSC_fit$lambda.min
# PSC_coef is a list with keys = khids of Precursor DEGs and values = series of attributes
# one of which (@ x) contains weights for all 14 input TF variables.
PSC_coef_list <- coef(PSC_fit, s = PSC_lambda)
```

#### With penalty, selecting for 3 TF Variables

Here, we set the dfmax parameter to 2, which applies a penalty and
selects for 3 TF variables.

One can observe the results given different values of this parameter,
such as 3 or 4, which respectively returns 4 or 5 variables, but we
study 3 variables in order to compare this result to our inferred
cocktail, which contains 3 TFs.

``` r
# Also took ~5s to run.
PSC_fit_3TFs <- cv.glmnet(glm_PSC_X_mat, glm_PSC_Y_mat, intercept = FALSE, family = "mgaussian", dfmax = 2)
PSC_lambda_3TFs <- PSC_fit_3TFs$lambda.min
PSC_coef_3TFs_list <- coef(PSC_fit_3TFs, s = PSC_lambda_3TFs)
```

### Strategy 1

#### Get Regression Coefficients for all Candidate Cocktail TFs

Make a coef_df, which is a data frame with rows = Candidate Cocktail
Precursor DEGs, columns = Precursor DEGs, and entries as regression
coefficients present in coef_list.

``` r
PSC_coef_df <- make_coef_df(PSC_coef_list,  PSC_candidate_cocktail_tfs, cleaned_PSC_precursor_pseudostate_DEGs)
PSC_coef_df_3TFs <- make_coef_df(PSC_coef_3TFs_list,  PSC_candidate_cocktail_tfs, cleaned_PSC_precursor_pseudostate_DEGs)
```

#### Visualize Distribution of Coefficients and get their Averages

We notice that Islet, Serum Response Factor, and xBPd have the highest
average coefficient values in the complete GMML model.

######## REF SUPP FIGURE PG 18

``` r
PSC_average_coef_df <- get_average_coefficient_df(PSC_coef_df)
PSC_coef_distribution_plot <- get_coefficients_distribution_plot(PSC_coef_df)
PSC_coef_distribution_plot

PSC_average_coef_df_3TFs <- get_average_coefficient_df(PSC_coef_df_3TFs)
```

#### Get Total Importance Scores for Candidate Cocktail TFs

> Sum the absolute values of the regression coefficients

######## REF SUPP FIGURE PG 19

``` r
PSC_total_GMML_importance_df <- get_total_GMML_importance_df(PSC_coef_df)
PSC_total_GMML_importance_df_3TFs <- get_total_GMML_importance_df(PSC_coef_df_3TFs)
```

### Strategy 2

#### First Get the ranks of the absolute values of the regression coefficients

``` r
ranked_PSC_coef_df <- PSC_coef_df %>% mutate(across(everything(), abs), across(everything(), min_rank))
ranked_PSC_coef_df_3TFs <- PSC_coef_df_3TFs %>% mutate(across(everything(), abs), across(everything(), min_rank))
```

#### Visualize Distribution of Coefficients and get their Averages

We notice that Islet, Serum Response Factor, and xBPd have the highest
ranked average coefficient values in the complete GMML model.

######## REF SUPP FIGURE PG 19, 20

``` r
ranked_PSC_average_coef_df <- get_average_coefficient_df(ranked_PSC_coef_df)
ranked_PSC_coef_distribution_plot <- get_coefficients_distribution_plot(ranked_PSC_coef_df)
ranked_PSC_coef_distribution_plot

ranked_PSC_average_coef_df_3TFs <- get_average_coefficient_df(ranked_PSC_coef_df_3TFs)
```

#### Get Total Importance Scores for Candidate Cocktail TFs

######## REF SUPP FIGURE PG 21

``` r
ranked_PSC_total_GMML_importance_df <- get_total_GMML_importance_df(ranked_PSC_coef_df)
ranked_PSC_total_GMML_importance_df_3TFs <- get_total_GMML_importance_df(ranked_PSC_coef_df_3TFs)
```

### Notes on Final Result

> As Foxg is missing, the PSC Cocktail is not correctly inferred from
> the complete GMML model. This demonstrates the need for the Group
> Lasso..

## (2) Comparison with GRNBoost2 Inferences

The purpose of this part of the analysis is to compare Cocktail
predictions derived from nonlinear models (GRNBoost2) with those from
the linear Group Lasso model, as TF-gene relationships may be better
modeled as nonlinear ones.

This Pipeline returns 2 lists: sorted_total_importance_tfs and
sorted_percent_top_list, which are used in tandem to infer the Cocktail.
Please refer to documentation in “GRNBoost2_pipeline.py” to see the
methods behind deriving the lists.

``` r
PSC_GRN_dir <- "GRN/PSC/"
```

#### PSC Precursor DEG Counts Matrix

``` r
# Matrix must have counts from both DEGs and Candidate Cocktail TFs.
PSC_precursor_pseudostate_DEGs_plus <- unique(c(PSC_precursor_pseudostate_DEGs, PSC_candidate_cocktail_khids))
PSC_precursor_pseudostate_GRN_mat <- URD_raw_data[PSC_precursor_pseudostate_DEGs_plus, PSC_precursor_pseudostate_cells]
# 339 x 91
PSC_precursor_pseudostate_GRN_df <- as.data.frame(t(PSC_precursor_pseudostate_GRN_mat))
```

#### Inputs to GRNBoost2 Pipeline

``` r
write_tsv(PSC_precursor_pseudostate_GRN_df, paste0("GRN/", "PSC_precursor_pseudostate_GRN_df.tsv"))
```

``` r
write_csv_helper(PSC_candidate_cocktail_khids, paste0("GRN/", "PSC_candidate_cocktail_khids.csv"))
```

#### Do GRNBoost2 50 times.

### Install Relevant Packages for “GRNBoost2_pipeline.py” into given r-reticulate environment.

``` r
conda_install("r-reticulate", "pandas")
conda_install("r-reticulate", "numpy")
py_install("arboreto", pip = T)
```

``` r
source_python("GRNBoost2_pipeline.py")
```

``` r
set.seed(16)

expression_count_matrix_tsv <- paste0("GRN/", "PSC_precursor_pseudostate_GRN_df.tsv")
tf_khid_csv <- paste0("GRN/", "PSC_candidate_cocktail_khids.csv")
outname <- "GRN/PSC/PSC"
min_i <- 0L
max_i <- 50L

# Took ~6 mins to run.
do_GRNBoost2(expression_count_matrix_tsv, tf_khid_csv, outname, min_i, max_i)
```

##### sorted_total_importance_tfs

Each element of the list is a tuple of TF and list with Mean and
Standard Deviation of Total Importance Scores. The TFs are sorted in
nonincreasing order based on the mean.

##### sorted_percent_top_list

Each element of the list is a TF and its Percentage Covered as a Top 3
Regulator for all DEGs.

######## REF SUPP FIGURE PG 21, 22

``` r
khid_tf_tsv_file <- "./helper_files/Ciona_khid_TF.tsv"
PSC_GRN_path <- "/Users/prakritipaul/Git/Cocktails/FINAL_code_and_results/GRN/PSC/*.tsv"
n <- 3L

# Returns sorted_total_importance_tfs and sorted_percent_top_list
PSC_GRNBoost2_outputs <- GRNBoost2_pipeline(PSC_GRN_path, khid_tf_tsv_file, n)
```

##### Inspect Output contents.

``` r
PSC_sorted_total_importance_tfs <- PSC_GRNBoost2_outputs[1]
PSC_sorted_percent_top_list <- PSC_GRNBoost2_outputs[2]
```

### Notes on Final Result

> The Cocktail appears within the first 5 TFs, which supports the
> results obtained from the Group Lasso Model.

## (3) Comparison with Method Applied on Terminal PSCs

A necessary comparison is with Terminal PSCs, or fully differentiated
PSCs, as a key claim in this work is that the Cocktail cannot be
inferred without studying cell states at the time of specification.

Potentially, we can find the Cocktail either in the top 50% mostly
highly expressed DEG TFs in the Terminal DEGs (these are Terminal
Candidate Cocktail TFs) or in the output of running the Group Lasso
Model using Terminal DEGs and Terminal Candidate Cocktail Factors. We
will also use the GMML.

From the following analyses, we demonstrate that the Cocktail cannot be
inferred from either way and build evidence for the power of our
computational method and the strength of the assumptions we made.

### Identify Terminal PSCs and their Siblings

Siblings were identified from observing the Virtual Lineage Tree.

``` r
terminal_PSCs <- neural_CT_list$PSCs
# 578
terminal_PSC_siblings <- c(neural_CT_list$`glia cells`, neural_CT_list$pATENs, neural_CT_list$`PSCs related`, neural_CT_list$CESNs)
```

### Raw Mats

``` r
terminal_PSC_raw_mat <- URD_raw_data[, terminal_PSCs]
terminal_PSC_TF_raw_mat <- URD_raw_data[tfs_in_URD, terminal_PSCs]

terminal_PSC_sibling_raw_mat <- URD_raw_data[, terminal_PSC_siblings]
terminal_PSC_sibling_TF_raw_mat <- URD_raw_data[tfs_in_URD, terminal_PSC_siblings]
```

### Calculate Terminal DEGs with respect to Siblings

Terminal DEGs are analogous to Precursor DEGs.

We calculate Terminal DEG TFs with respect to Siblings because
comparison to all other neural cell types did not render meaningful TFs.
Instead, TFs known to play roles in specification appeared as
differentially expressed when Terminal PSCs were compared to its
siblings. This made these results more competitive with the results
obtained from applying the computational method to Precursors at the
time of specification.

### Get Terminal DEGs

``` r
terminal_PSC_v_siblings_mat <- as.matrix(cbind(terminal_PSC_raw_mat, terminal_PSC_sibling_raw_mat))

terminal_PSC_v_siblings_edgeR_outputs <- complete_edgeR_pipeline(terminal_PSC_v_siblings_mat, ncol(terminal_PSC_raw_mat), ncol(terminal_PSC_sibling_raw_mat), 0.05, -1, 50)
```

``` r
# 494
terminal_PSC_v_siblings_DEGs_df <- terminal_PSC_v_siblings_edgeR_outputs[[1]]
terminal_PSC_v_siblings_DEGs <- rownames(terminal_PSC_v_siblings_DEGs_df)
```

#### Write them out.

``` r
write_csv_helper(terminal_PSC_v_siblings_DEGs, paste0(DEG_lists_out_dir, "terminal_PSC_v_siblings_DEGs.csv"))
```

### Get Terminal DEG TFs

We will filter these downstream to get Terminal Candidate Cocktail TFs.

``` r
terminal_PSC_v_siblings_TF_mat <- as.matrix(cbind(terminal_PSC_TF_raw_mat, terminal_PSC_sibling_TF_raw_mat))

terminal_PSC_v_siblings_TF_edgeR_outputs <- complete_edgeR_pipeline(terminal_PSC_v_siblings_TF_mat, ncol(terminal_PSC_TF_raw_mat), ncol(terminal_PSC_sibling_TF_raw_mat), 0.05, -1, 50)
```

``` r
# 20
terminal_PSC_v_siblings_TF_DEGs_df <- terminal_PSC_v_siblings_TF_edgeR_outputs[[1]]
# Get TF orthologs.
named_terminal_PSC_v_siblings_TF_DEGs_df <- convert_khids_to_ortholog(terminal_PSC_v_siblings_TF_DEGs_df)
terminal_PSC_v_sibling_TF_DEGs <- named_terminal_PSC_v_siblings_TF_DEGs_df$Ghost.Name
```

#### Write them out.

``` r
write_csv_helper(terminal_PSC_v_sibling_TF_DEGs, paste0(DEG_lists_out_dir, "terminal_PSC_v_sibling_TF_DEGs.csv"))
```

### Get Terminal Candidate Cocktail TFs

As before, we will take the top 50% most highly expressed of the
Terminal DEG TFs and consider them as the Terminal Candidate Cocktail
TFs.

######## REF SUPP FIGURE PG 23

``` r
# Get mean Expression of all DEG TFs
terminal_PSC_mean_TF_exp_df <- make_named_mean_TF_exp_df(terminal_PSCs)

# Get khids of Terminal DEG TFs
terminal_PSC_v_siblings_TF_DEG_khids <- rownames(terminal_PSC_v_siblings_TF_DEGs_df)

# Get Mean expression of Terminal DEG TFs in terminal cells.
terminal_PSC_cocktail_TF_mean_exp_df <- terminal_PSC_mean_TF_exp_df %>% filter(KH.gene.model %in% terminal_PSC_v_siblings_TF_DEG_khids)

# Get top 50% most highly expressed of above DEG TFS
terminal_cutoff <- 0.5

terminal_PSC_cutoff_df <- terminal_PSC_cocktail_TF_mean_exp_df[1:(terminal_cutoff*(length(terminal_PSC_v_siblings_TF_DEG_khids))), ]

# Get the TFs and khids
# 10 
terminal_PSC_candidate_cocktail_khids <- terminal_PSC_cutoff_df$KH.gene.model
terminal_PSC_candidate_cocktail_tfs <- terminal_PSC_cutoff_df$Ghost.Name
```

#### Write out the Candidate Cocktail TFs

``` r
write_csv_helper(terminal_PSC_candidate_cocktail_tfs, paste0(DEG_lists_out_dir, "terminal_PSC_candidate_cocktail_tfs.csv"))
```

#### Brief Aside

An interesting observation is that the common Candidate Cocktail TFs
between the Precursor and Terminal PSC are Islet and Foxg (Key
Determinant)

``` r
common_PSC_cocktail_TFs <- intersect(terminal_PSC_candidate_cocktail_tfs, PSC_candidate_cocktail_tfs)
common_PSC_cocktail_TFs
```

### (a) Do Cocktail TF Inference using Group Lasso Model

#### Variables

``` r
# Get rid of TFs in DEGs.
# 10
common_terminal_PSC_DEG_TFs <- intersect(terminal_PSC_candidate_cocktail_khids, terminal_PSC_v_siblings_DEGs)
# 484
cleaned_terminal_PSC_DEGs <- setdiff(terminal_PSC_v_siblings_DEGs, common_terminal_PSC_DEG_TFs)
```

``` r
# 138, 10, 484
terminal_PSC_n <- length(terminal_PSCs)
terminal_PSC_p <- length(terminal_PSC_candidate_cocktail_tfs)
terminal_PSC_k <- length(cleaned_terminal_PSC_DEGs)
```

#### Make Group Lasso X Matrix

``` r
# 138 x 10
terminal_PSC_TFs_in_terminal_mat <- t(URD_data[terminal_PSC_candidate_cocktail_khids, terminal_PSCs])
```

``` r
# 66792 x 4840
terminal_PSC_grplasso_X_mat <- make_grplasso_X_mat(terminal_PSC_TFs_in_terminal_mat, terminal_PSC_n, terminal_PSC_p, terminal_PSC_k)
```

#### Make Group Lasso Y Matrix

``` r
# 138 x 484
terminal_PSC_DEGs_in_terminal_mat <- t(URD_data[cleaned_terminal_PSC_DEGs, terminal_PSCs])
```

``` r
# 66792 x 1
terminal_PSC_grplasso_Y_mat <- make_grplasso_Y_mat(terminal_PSC_DEGs_in_terminal_mat, terminal_PSCs, terminal_PSC_p, terminal_PSC_k)
```

#### Variables: Group Lasso Matrices X and Y and Group IDs

``` r
terminal_grplasso_PSC_groups <- rep(seq(1, terminal_PSC_p), terminal_PSC_k)

terminal_PSC_X <- terminal_PSC_grplasso_X_mat
terminal_PSC_Y <- terminal_PSC_grplasso_Y_mat
terminal_PSC_group_id <- terminal_grplasso_PSC_groups
```

#### Run Group Lasso

``` r
# 6h 
terminal_grplasso_PSC_cv_grepreg_fit <- cv.grpreg(terminal_PSC_X, terminal_PSC_Y, terminal_PSC_group_id, penalty="grLasso")
```

### Get Sets of Selected TF Groups at Various Values of Lambda

Only Isl appears in TF Sets up to size 6. Thus, Cocktail cannot be
uncovered.

######## REF SUPP FIGURE PG 24

``` r
# Key: lambda; Values: Set of Selected TF Groups
terminal_PSC_group_lassoed_tf_list <- get_group_lassoed_tf_list(terminal_grplasso_PSC_cv_grepreg_fit, terminal_PSC_candidate_cocktail_tfs)

# Key: lambda; Values: Set of Non-selected TF Groups
terminal_PSC_non_group_lassoed_TF_list <- lapply(terminal_PSC_group_lassoed_tf_list, get_non_group_lassoed_TFs, terminal_PSC_candidate_cocktail_tfs)
```

#### Let us Visualize the Number of Groups Selected at Various Values of Lambda

``` r
plot(terminal_grplasso_PSC_cv_grepreg_fit)
```

### Get Sorted Total TF Importance Lists

The Total Importance Scores within these lists give us a quantitative
measure of a Cocktail TF’s strength over the Gene Regulatory Network (as
defined by the Terminal DEGs) underlying the Terminally differentiated
cell type’s state.

``` r
# (For all lambdas) key = lambda; value is a list with key = TF and value = non0 betas.
terminal_PSC_TF_importances_lists <- get_TF_importances_list(terminal_grplasso_PSC_cv_grepreg_fit)

# Sum the absolute values of these betas for all TFs for all lambdas.
terminal_PSC_sum_total_TF_importance_lists <- terminal_PSC_TF_importances_lists %>% map(get_total_summed_TF_importance_list)

# Make sorted dfs given these Total Importance Scores.
terminal_PSC_sorted_sum_total_TF_importance_lists <- terminal_PSC_sum_total_TF_importance_lists %>% map(make_sorted_TF_importance_df, terminal_PSC_candidate_cocktail_tfs)
```

#### Result: We are unable to uncover the PSC Cocktail.

### (b) Do GMML

``` r
glm_terminal_PSC_X_mat <- as.matrix(terminal_PSC_TFs_in_terminal_mat)
glm_terminal_PSC_Y_mat <- as.matrix(terminal_PSC_DEGs_in_terminal_mat)
```

### Perform Regression

#### Without penalty

We will refer to this as the complete GMML model.

``` r
# ~Took 30 to run.
terminal_PSC_fit <- cv.glmnet(glm_terminal_PSC_X_mat, glm_terminal_PSC_Y_mat, intercept = FALSE, family = "mgaussian")
terminal_PSC_lambda <- terminal_PSC_fit$lambda.min
# terminal_PSC_coef is a list with keys = khids of Precursor DEGs and values = series of attributes
# one of which (@ x) contains weights for all 14 input TF variables.
terminal_PSC_coef_list <- coef(terminal_PSC_fit, s = terminal_PSC_lambda)
```

#### With penalty, selecting for 3 TF Variables

Here, we set the dfmax parameter to 2, which applies a penalty and
selects for 3 TF variables.

One can observe the results given different values of this parameter,
such as 3 or 4, which respectively returns 4 or 5 variables, but we
study 3 variables in order to compare this result to our inferred
cocktail, which contains 3 TFs.

######## REF SUPP FIGURE PG 25, 26

``` r
# Took ~20s to run.
terminal_PSC_fit_3TFs <- cv.glmnet(glm_terminal_PSC_X_mat, glm_terminal_PSC_Y_mat, intercept = FALSE, family = "mgaussian", dfmax = 2)
terminal_PSC_lambda_3TFs <- terminal_PSC_fit_3TFs$lambda.min
terminal_PSC_coef_3TFs_list <- coef(terminal_PSC_fit_3TFs, s = terminal_PSC_lambda_3TFs)
```

### Strategy 1

#### Get Regression Coefficients for all Candidate Cocktail TFs

Make a coef_df, which is a data frame with rows = Candidate Cocktail
Precursor DEGs, columns = Precursor DEGs, and entries as regression
coefficients present in coef_list.

``` r
terminal_PSC_coef_df <- make_coef_df(terminal_PSC_coef_list,  terminal_PSC_candidate_cocktail_tfs, cleaned_terminal_PSC_DEGs)
terminal_PSC_coef_df_3TFs <- make_coef_df(terminal_PSC_coef_3TFs_list,  terminal_PSC_candidate_cocktail_tfs, cleaned_terminal_PSC_DEGs)
```

#### Visualize Distribution of Coefficients and get their Averages

We notice that Isl, Ci-ZF220 and Tbx2/3 have the highest average
coefficient values in the complete GMML model.

######## REF SUPP FIGURE PG 25, 26

``` r
terminal_PSC_average_coef_df <- get_average_coefficient_df(terminal_PSC_coef_df)
terminal_PSC_coef_distribution_plot <- get_coefficients_distribution_plot(terminal_PSC_coef_df)
terminal_PSC_coef_distribution_plot

terminal_PSC_average_coef_df_3TFs <- get_average_coefficient_df(terminal_PSC_coef_df_3TFs)
```

#### Get Total Importance Scores for Candidate Cocktail TFs

######## REF SUPP FIGURE PG 27

``` r
terminal_terminal_PSC_total_GMML_importance_df <- get_total_GMML_importance_df(terminal_PSC_coef_df)
terminal_PSC_total_GMML_importance_df_3TFs <- get_total_GMML_importance_df(terminal_PSC_coef_df_3TFs)
```

### Strategy 2

#### First Get the ranks of the absolute values of the regression coefficients

``` r
ranked_terminal_PSC_coef_df <- terminal_PSC_coef_df %>% mutate(across(everything(), abs), across(everything(), min_rank))
ranked_terminal_PSC_coef_df_3TFs <- terminal_PSC_coef_df_3TFs %>% mutate(across(everything(), abs), across(everything(), min_rank))
```

#### Visualize Distribution of Coefficients and get their Averages

We notice that Isl, Ci-ZF220 and Tbx2/3 have the highest ranked average
coefficient values in the complete GMML model.

######## REF SUPP FIGURE PG 25, 26

``` r
ranked_terminal_PSC_average_coef_df <- get_average_coefficient_df(ranked_terminal_PSC_coef_df)
ranked_terminal_PSC_coef_distribution_plot <- get_coefficients_distribution_plot(ranked_terminal_PSC_coef_df)
ranked_terminal_PSC_coef_distribution_plot

ranked_terminal_PSC_average_coef_df_3TFs <- get_average_coefficient_df(ranked_terminal_PSC_coef_df_3TFs)
```

#### Get Total Importance Scores for Candidate Cocktail TFs

######## REF SUPP FIGURE PG 27

``` r
ranked_terminal_PSC_total_GMML_importance_df <- get_total_GMML_importance_df(ranked_terminal_PSC_coef_df)
ranked_terminal_PSC_total_GMML_importance_df_3TFs <- get_total_GMML_importance_df(ranked_terminal_PSC_coef_df_3TFs)
```

### (c) Comparison with GRNBoost2 Inferences

The purpose of this part of the analysis is to compare Cocktail
predictions derived from nonlinear models (GRNBoost2) with those from
the linear Group Lasso model, as TF-gene relationships may be better
modeled as nonlinear ones.

This Pipeline returns 2 lists: sorted_total_importance_tfs and
sorted_percent_top_list, which are used in tandem to infer the Cocktail.
Please refer to documentation in “GRNBoost2_pipeline.py” to see the
methods behind deriving the lists.

``` r
terminal_PSC_GRN_dir <- "GRN/terminal_PSC/"
```

##### Terminal PSC DEG Counts Matrix

``` r
# Matrix must have counts from both DEGs and Candidate Cocktail TFs.
terminal_PSC_DEGs_plus <- unique(c(terminal_PSC_v_siblings_DEGs, terminal_PSC_candidate_cocktail_khids))
terminal_PSC_GRN_mat <- URD_raw_data[terminal_PSC_DEGs_plus, terminal_PSCs]

# 65 x 713
terminal_PSC_GRN_df <- as.data.frame(t(terminal_PSC_GRN_mat))
```

#### Inputs to GRNBoost2 Pipeline

``` r
write_tsv(terminal_PSC_GRN_df, paste0("GRN/", "terminal_PSC_GRN_df.tsv"))
```

``` r
write_csv_helper(terminal_PSC_candidate_cocktail_khids, paste0("GRN/", "terminal_PSC_candidate_cocktail_khids.csv"))
```

#### Do GRNBoost2 50 times.

``` r
set.seed(16)

expression_count_matrix_tsv <- paste0("GRN/", "terminal_PSC_GRN_df.tsv")
tf_khid_csv <- paste0("GRN/", "terminal_PSC_candidate_cocktail_khids.csv")
outname <- "GRN/terminal_PSC/terminal_PSC"
min_i <- 0L
max_i <- 50L

# Took ~23 mins to run.
do_GRNBoost2(expression_count_matrix_tsv, tf_khid_csv, outname, min_i, max_i)
```

##### sorted_total_importance_tfs

Each element of the list is a tuple of TF and list with Mean and
Standard Deviation of Total Importance Scores. The TFs are sorted in
nonincreasing order based on the mean.

##### sorted_percent_top_list

Each element of the list is a TF and its Percentage Covered as a Top 3
Regulator for all DEGs.

######## REF SUPP FIGURE PG 28

``` r
khid_tf_tsv_file <- "./helper_files/Ciona_khid_TF.tsv"
terminal_PSC_GRN_path <- "/Users/prakritipaul/Git/Cocktails/FINAL_code_and_results/GRN/terminal_PSC/*.tsv"
n <- 3L

# Returns sorted_total_importance_tfs and sorted_percent_top_list
terminal_PSC_GRNBoost2_outputs <- GRNBoost2_pipeline(terminal_PSC_GRN_path, khid_tf_tsv_file, n)
```

##### Inspect Output contents.

``` r
terminal_PSC_sorted_total_importance_tfs <- terminal_PSC_GRNBoost2_outputs[1]
terminal_PSC_sorted_percent_top_list <- terminal_PSC_GRNBoost2_outputs[2]
```

#### Note on Final Result

As our Cocktail TFs are not among the top TFs given these two measures,
we are once again able to uncover the Cocktail from studying the
terminal state. This further supports the need to study gene regulatory
changes occurring at the time of specification.

## (4) Subsampling

> Finally, in order to build evidence for our predictions, we randomly
> subsample 80% of the data, run the Group Lasso Model, and keep track
> of its inferred Cocktail. We count the number of times that a TF in
> our predicted Cocktail appears in a subsample’s inferred Cocktail and
> thus get a percentage/confidence for the overall prediction.

### Get 50 runs of Group Lasso Model on Subsampled Data.

> The results of this subsampling can be found in the given directory-
> each object has 5 runs of cv.grpreg, resulting in 50 samples/runs. The
> following code makes a “combined_subsampled_tf_df”, a data frame with
> row = run\_**lambda**, where run = the round of subsampling/run of
> group lasso and lambda = one of 100 values of lambda used during
> cross-validation; and column value = string that is a concatenated
> name of the set of TF Groups selected at the corresponding value of
> lambda (“TF-string”). These “combined_subsampled_tf_df”s are
> subsequently processed in “subsampling_prediction_confidence.py”.

##### Access subsamples.

``` r
if (!(file.exists("/tmp/PSC_subsamples"))) {
  cat("Downloading PSC_subsamples directory\n")
  options(timeout = 600)
  PSC_subsampling_link <- "https://cocktails-r-objects.s3.amazonaws.com/PSC_subsamples.zip"
  # Download in /tmp.
  download.file(PSC_subsampling_link, destfile = "/tmp/PSC_subsamples.zip")
  cat("Unzipping\n")
  unzip("/tmp/PSC_subsamples.zip", exdir = "/tmp")
  
} else {
  cat("PSC_subsamples directory is present.\n")
}
```

##### Make combined_subsampled_tf_df.

``` r
PSC_combined_subsampled_tf_df <- make_combined_subsampled_tf_df("/tmp/PSC_subsamples/")
```

``` r
subsampling_out_dir <- "./subsampling/"
```

##### Write out combined_subsampled_tf_df.

``` r
write_out_combined_subsampled_tf_df(PSC_combined_subsampled_tf_df, subsampling_out_dir, "PSC_combined_subsampled_tf_df.csv")
```

##### Also write out PSC’s group_lassoed_tf_list (Has Predicted Cocktail).

``` r
write_out_CT_group_lassoed_tf_list(PSC_group_lassoed_tf_list, subsampling_out_dir, "PSC_group_lassoed_tf_df.csv")
```

### Run Subsampling Prediction Confidence Pipeline.

``` r
source_python("subsampling_prediction_confidence.py")
```

``` r
PSC_group_lassoed_tf_df_csv <- paste0(subsampling_out_dir, "PSC_group_lassoed_tf_df.csv")
PSC_subsampled_tf_df_csv <- paste0(subsampling_out_dir, "PSC_combined_subsampled_tf_df.csv")
PSC_tf_set_sizes <- list(3, 4, 5, 6)
PSC_ideal_set_size <- 4

PSC_cocktail_subsampling_confidence <- subsampling_prediction_confidence(PSC_group_lassoed_tf_df_csv, PSC_subsampled_tf_df_csv, PSC_tf_set_sizes, PSC_ideal_set_size)
```

### Note on Final Result

Foxg, Islet, DRIL1/2, and xBPd each appear 100% time in TF Sets of size
4.

# RSession Info

> Refer to Supplement for complete information.

``` r
sessionInfo()
```

    ## R version 4.3.0 (2023-04-21)
    ## Platform: x86_64-apple-darwin20 (64-bit)
    ## Running under: macOS Ventura 13.4.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.3.0  fastmap_1.1.1   cli_3.6.1       tools_4.3.0    
    ##  [5] htmltools_0.5.5 rstudioapi_0.14 yaml_2.3.7      rmarkdown_2.21 
    ##  [9] knitr_1.42      xfun_0.39       digest_0.6.31   rlang_1.1.1    
    ## [13] evaluate_0.21
