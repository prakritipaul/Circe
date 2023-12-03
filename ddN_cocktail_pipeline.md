ddN Cocktail Inference Pipeline
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
    - <a href="#get-ids-of-glra-precursor-segment-and-its-sibling-segment"
      id="toc-get-ids-of-glra-precursor-segment-and-its-sibling-segment">Get
      IDs of GLRA Precursor Segment and its Sibling Segment</a>
    - <a href="#get-glra-precursor-and-sibling-segment-cells"
      id="toc-get-glra-precursor-and-sibling-segment-cells">Get GLRA Precursor
      and Sibling Segment Cells</a>
    - <a href="#order-cells-by-pseudotime"
      id="toc-order-cells-by-pseudotime">Order Cells by Pseudotime</a>
  - <a href="#mother-segment" id="toc-mother-segment">Mother Segment</a>
    - <a href="#get-ids-of-glra-mother-segment-and-its-siblings"
      id="toc-get-ids-of-glra-mother-segment-and-its-siblings">Get IDs of GLRA
      Mother Segment and its Siblings</a>
    - <a href="#get-glra-mother-and-sibling-segment-cells"
      id="toc-get-glra-mother-and-sibling-segment-cells">Get GLRA Mother and
      Sibling Segment Cells</a>
    - <a href="#order-cells-by-pseudotime-1"
      id="toc-order-cells-by-pseudotime-1">Order Cells by Pseudotime</a>
  - <a href="#visualize-cells-on-virtual-neural-lineage-tree"
    id="toc-visualize-cells-on-virtual-neural-lineage-tree">Visualize Cells
    on Virtual Neural Lineage Tree</a>
- <a
  href="#2-get-pc-coordinates-and-corresponding-pseudotimes-for-glra-cells"
  id="toc-2-get-pc-coordinates-and-corresponding-pseudotimes-for-glra-cells">(2)
  Get PC Coordinates and Corresponding Pseudotimes for GLRA Cells</a>
  - <a href="#make-seurat-object" id="toc-make-seurat-object">Make Seurat
    Object</a>
  - <a href="#use-elbow-plot-to-choose-the-number-of-pcs-15"
    id="toc-use-elbow-plot-to-choose-the-number-of-pcs-15">Use Elbow Plot to
    choose the number of PCs (15)</a>
  - <a href="#key-variables" id="toc-key-variables">Key Variables</a>
- <a href="#3-get-pseudostates" id="toc-3-get-pseudostates">(3) Get
  Pseudostates</a>
  - <a href="#make-pc_pseudotime_dfs" id="toc-make-pc_pseudotime_dfs">Make
    PC_pseudotime_dfs</a>
  - <a href="#export-pc_pseudotime_dfs"
    id="toc-export-pc_pseudotime_dfs">Export PC_pseudotime_dfs</a>
  - <a href="#load-get_pseudostatespy" id="toc-load-get_pseudostatespy">Load
    “get_pseudostates.py”</a>
  - <a href="#key-variables-1" id="toc-key-variables-1">Key Variables</a>
  - <a href="#ddn-precursor" id="toc-ddn-precursor">ddN Precursor</a>
  - <a href="#glra-precursor-sibling-31"
    id="toc-glra-precursor-sibling-31">GLRA Precursor Sibling 31</a>
  - <a href="#glra-precursor-sibling-45"
    id="toc-glra-precursor-sibling-45">GLRA Precursor Sibling 45</a>
  - <a href="#glra-mother" id="toc-glra-mother">GLRA Mother</a>
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
- <a href="#final-result-inferred-ddn-cocktail"
  id="toc-final-result-inferred-ddn-cocktail">Final Result: Inferred ddN
  Cocktail</a>
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
      - <a href="#without-lasso-penalty" id="toc-without-lasso-penalty">Without
        Lasso penalty</a>
      - <a href="#with-lasso-penalty-selecting-for-3-tf-variables"
        id="toc-with-lasso-penalty-selecting-for-3-tf-variables">With Lasso
        penalty, selecting for 3 TF Variables</a>
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
    - <a href="#ddn-precursor-deg-counts-matrix"
      id="toc-ddn-precursor-deg-counts-matrix">ddN Precursor DEG Counts
      Matrix</a>
    - <a href="#inputs-to-grnboost2-pipeline"
      id="toc-inputs-to-grnboost2-pipeline">Inputs to GRNBoost2 Pipeline</a>
    - <a href="#do-grnboost2-50-times" id="toc-do-grnboost2-50-times">Do
      GRNBoost2 50 times.</a>
      - <a href="#sorted_total_importance_tfs"
        id="toc-sorted_total_importance_tfs">sorted_total_importance_tfs</a>
      - <a href="#sorted_percent_top_list"
        id="toc-sorted_percent_top_list">sorted_percent_top_list</a>
      - <a href="#inspect-output-contents"
        id="toc-inspect-output-contents">Inspect Output contents.</a>
    - <a href="#notes-on-final-result-1"
      id="toc-notes-on-final-result-1">Notes on Final Result</a>
  - <a href="#3-comparison-with-method-applied-on-terminal-ddns"
    id="toc-3-comparison-with-method-applied-on-terminal-ddns">(3)
    Comparison with Method Applied on Terminal ddNs</a>
    - <a href="#identify-terminal-ddns-and-their-siblings"
      id="toc-identify-terminal-ddns-and-their-siblings">Identify Terminal
      ddNs and their Siblings</a>
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
      - <a href="#error-group-lasso-does-not-solve"
        id="toc-error-group-lasso-does-not-solve">Error: Group Lasso does not
        solve.</a>
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
      - <a href="#result" id="toc-result">Result</a>
    - <a href="#b-do-gmml" id="toc-b-do-gmml">(b) Do GMML</a>
    - <a href="#perform-regression-1" id="toc-perform-regression-1">Perform
      Regression</a>
      - <a href="#without-lasso-penalty-1"
        id="toc-without-lasso-penalty-1">Without Lasso penalty</a>
      - <a href="#with-lasso-penalty-selecting-for-3-tf-variables-1"
        id="toc-with-lasso-penalty-selecting-for-3-tf-variables-1">With Lasso
        penalty, selecting for 3 TF Variables</a>
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
      - <a href="#terminal-ddn-deg-counts-matrix"
        id="toc-terminal-ddn-deg-counts-matrix">Terminal ddN DEG Counts
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
        - <a
          href="#we-are-once-again-able-to-uncover-the-cocktail-from-studying-the-terminal-state-this-further-supports-the-need-to-study-gene-regulatory-changes-occurring-at-the-time-of-specification"
          id="toc-we-are-once-again-able-to-uncover-the-cocktail-from-studying-the-terminal-state-this-further-supports-the-need-to-study-gene-regulatory-changes-occurring-at-the-time-of-specification">We
          are once again able to uncover the Cocktail from studying the terminal
          state. This further supports the need to study gene regulatory changes
          occurring at the time of specification.</a>
  - <a href="#4-subsampling" id="toc-4-subsampling">(4) Subsampling</a>
    - <a href="#subsample-8-cells-9-times"
      id="toc-subsample-8-cells-9-times">Subsample 8 cells 9 times.</a>
- <a href="#do-subsampling-by-leaving-one-cell-out"
  id="toc-do-subsampling-by-leaving-one-cell-out">Do subsampling by
  leaving one cell out.</a>
- <a href="#lets-make-the-subsampled-group-tf-df"
  id="toc-lets-make-the-subsampled-group-tf-df">Let’s make the subsampled
  group tf df.</a>
- <a href="#write-it-out" id="toc-write-it-out">Write it out.</a>
- <a href="#also-write-out-ct_group_lassoed_tf_list"
  id="toc-also-write-out-ct_group_lassoed_tf_list">Also write out
  CT_group_lassoed_tf_list.</a>
  - <a href="#run-subsampling-prediction-confidence-pipeline"
    id="toc-run-subsampling-prediction-confidence-pipeline">Run Subsampling
    Prediction Confidence Pipeline</a>
- <a href="#rsession-info" id="toc-rsession-info">RSession Info</a>

# Brief Introduction

The analysis pipeline present here demonstrates **Circe**, a suite of
computational methods that predict reprogramming Cocktails for neuronal
cell types in Ciona intestinalis, as **applied to the Decussating
Neurons (ddNs)**. A **reprogramming Cocktail is a combination of
transcription factors** that when overexpressed, can induce a particular
cellular identity in the cell type in which it was introduced.

**There are 3 Main Computational Methods:**  
**Method 1: Identification of Pseudostates, detailed in:**  
(1) Identify Precursor, Mother, and their Sibling(s) Segment Cell
Populations on Virtual Neural Lineage Tree  
(2) Get PC Coordinates and Corresponding Pseudotimes for BTN Cells   (3)
Get Pseudostates  

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

**The final results (“Collier BTN Cocktail” and “POU IV” BTN Cocktail”)
also presented here.**  
**Experimental Validations of both Cocktails can be found** in “Circe
Predicts Neurogenin, HNF6, and Collier as the BTN Cocktail” and
“Neurogenin, HNF6, and POU IV are Sufficient to Reprogram Cell Types
into BTNs” subsections.

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

> Note: GLRA refers to the cell segment that contains both GLRA and
> non-GLRA cells. GLRAs express Dmbx+, while the latter does not. Thus,
> throughout this pipeline, GLRA ID’s and siblings will be used, but
> later subsetted to get GLRAs.

## Precursor Segment

### Get IDs of GLRA Precursor Segment and its Sibling Segment

``` r
# Get ID from URD Object.
GLRA_ID <- "38"
# 31, 45
GLRA_sibling_IDs <- segSiblings(CNS.k300s6w4.tree.built.name, GLRA_ID, include.self = F)
GLRA_sibling_31_ID <- GLRA_sibling_IDs[1]
GLRA_sibling_45_ID <- GLRA_sibling_IDs[2]
```

### Get GLRA Precursor and Sibling Segment Cells

> Key Note: GLRA_segment_cells contain both ddN (Dmbx+ cells) and
> non-ddN cells. We will need to differentiate between the two below in
> order to correctly identify the dd Precursor Pseudostate.

``` r
# 82
GLRA_segment_cells <- cellsInCluster(CNS.k300s6w4.tree.built.name, clustering="segment", GLRA_ID)
# 139
GLRA_sibling_31_segment_cells <- cellsInCluster(CNS.k300s6w4.tree.built.name, clustering="segment", GLRA_sibling_31_ID)
# 27 (Will need to take everything.) 
GLRA_sibling_45_segment_cells <- cellsInCluster(CNS.k300s6w4.tree.built.name, clustering="segment", GLRA_sibling_45_ID)
```

### Order Cells by Pseudotime

``` r
pseudotimed_GLRA_segment_cells <- get_pseudotimed_cells(cells_in_segments_list, pseudotimes, GLRA_ID, length(GLRA_segment_cells), ascending=TRUE)

# Get GLRA cells that express Dmbx- these are the GLRAs.
# 16/82
pseudotimed_ddN_precursor_cells <- get_dmbx_expressing_cells(pseudotimed_GLRA_segment_cells, pseudotimes)

pseudotimed_GLRA_sibling_31_segment_cells <- get_pseudotimed_cells(cells_in_segments_list, pseudotimes, GLRA_sibling_31_ID, length(GLRA_sibling_31_segment_cells), ascending=TRUE)

pseudotimed_GLRA_sibling_45_segment_cells <- get_pseudotimed_cells(cells_in_segments_list, pseudotimes, GLRA_sibling_45_ID, length(GLRA_sibling_45_segment_cells), ascending=TRUE)
```

## Mother Segment

### Get IDs of GLRA Mother Segment and its Siblings

``` r
GLRA_mother_ID <- "59"
# 64, 15
GLRA_mother_sibling_IDs <- segSiblings(CNS.k300s6w4.tree.built.name, GLRA_mother_ID, include.self = F)
GLRA_mother_sibling_64_ID <- GLRA_mother_sibling_IDs[1]
GLRA_mother_sibling_15_ID <- GLRA_mother_sibling_IDs[2]
```

### Get GLRA Mother and Sibling Segment Cells

``` r
# 395
GLRA_mother_segment_cells <- cellsInCluster(CNS.k300s6w4.tree.built.name, clustering="segment", GLRA_mother_ID)
# 147
GLRA_mother_sibling_64_segment_cells <- cellsInCluster(CNS.k300s6w4.tree.built.name, clustering="segment", GLRA_mother_sibling_64_ID)
# 89
GLRA_mother_sibling_15_segment_cells <- cellsInCluster(CNS.k300s6w4.tree.built.name, clustering="segment", GLRA_mother_sibling_15_ID)
```

### Order Cells by Pseudotime

``` r
pseudotimed_GLRA_mother_segment_cells <- get_pseudotimed_cells(cells_in_segments_list, pseudotimes, GLRA_mother_ID, length(GLRA_mother_segment_cells), ascending=FALSE)

pseudotimed_GLRA_mother_sibling_64_segment_cells <- get_pseudotimed_cells(cells_in_segments_list, pseudotimes, GLRA_mother_sibling_64_ID, length(GLRA_mother_sibling_64_segment_cells), ascending=FALSE)

pseudotimed_GLRA_mother_sibling_15_segment_cells <- get_pseudotimed_cells(cells_in_segments_list, pseudotimes, GLRA_mother_sibling_15_ID, length(GLRA_mother_sibling_15_segment_cells), ascending=FALSE)
```

## Visualize Cells on Virtual Neural Lineage Tree

``` r
GLRA_segment_cells_list <- list("pseudotimed_ddN_precursor_cells" = pseudotimed_ddN_precursor_cells,
                               "pseudotimed_GLRA_sibling_31_segment_cells" = pseudotimed_GLRA_sibling_31_segment_cells,
                               "pseudotimed_GLRA_sibling_45_segment_cells" = pseudotimed_GLRA_sibling_45_segment_cells,
                               "pseudotimed_GLRA_mother_segment_cells" = pseudotimed_GLRA_mother_segment_cells,
                               "pseudotimed_GLRA_mother_sibling_64_segment_cells" = pseudotimed_GLRA_mother_sibling_64_segment_cells,
                               "pseudotimed_GLRA_mother_sibling_15_segment_cells" = pseudotimed_GLRA_mother_sibling_15_segment_cells)

modified_URD <- add_segment_state(GLRA_segment_cells_list, "GLRA_segment_cells_list")
plotTree(modified_URD, label.segments = F, label = "segment_state", label.type = "group")
```

# (2) Get PC Coordinates and Corresponding Pseudotimes for GLRA Cells

The PC coordinates and corresponding pseudotimes for GLRA Cells will
help determine where to cut the segments and find pseudostates for ddN
Precursors, GLRA Sibling(s) and Mothers. Pseudostates are approximations
of true Precursor, Precursor Sibling(s) and Mother states. They will be
used to calculate Differentially Expressed Genes in the next part of the
pipeline.

## Make Seurat Object

``` r
GLRA_cells <- c(pseudotimed_ddN_precursor_cells, 
                pseudotimed_GLRA_sibling_31_segment_cells,
                pseudotimed_GLRA_sibling_45_segment_cells,
                pseudotimed_GLRA_mother_segment_cells, 
                pseudotimed_GLRA_mother_sibling_64_segment_cells, 
                pseudotimed_GLRA_mother_sibling_15_segment_cells)
  
GLRA_raw_mat <- URD_raw_data[, GLRA_cells]
GLRA_seurat <- CreateSeuratObject(counts = GLRA_raw_mat,
                                min.cells = 3,
                                min.features = 200)
# Normalize Data.
GLRA_seurat <- NormalizeData(GLRA_seurat)
# Find Variable Feature.
GLRA_seurat <- FindVariableFeatures(GLRA_seurat, selection.method = "vst", nfeatures = 2000)
# Scale data.
all_GLRA_seurat_genes <- rownames(GLRA_seurat)
# use.umi argument to regress on umi count data.
GLRA_seurat <- ScaleData(GLRA_seurat, features = all_GLRA_seurat_genes, use.umi = TRUE)
# Do PCA.
GLRA_seurat <- RunPCA(GLRA_seurat, features = VariableFeatures(object = GLRA_seurat))
```

## Use Elbow Plot to choose the number of PCs (15)

``` r
ElbowPlot(GLRA_seurat)
```

## Key Variables

``` r
GLRA_PCA_COORDINATES <- Embeddings(GLRA_seurat, reduction = "pca")
GLRA_NUM_PCS <- 15
```

# (3) Get Pseudostates

As described, pseudostates are approximations of true Precursor,
Precursor Sibling, and Mother populations. We will first make
“PC_pseudotime_dfs”, which are data frames that contain PC coordinates
and pseudotimes for the GLRA cells.

The data frames of Precursor Precursor Sibling, and Mother segments will
exported as csv files, which will then be used as inputs to modules in
“get_pseudostates.py”.

PC_pseudotime_dfs of Mother’s Siblings will be subjected to a separate
function (get_est_mother_sibling_cells) to get their pseudostates.

## Make PC_pseudotime_dfs

``` r
ddN_PC_pseudotime_df <- make_PC_pseudotime_df(pseudotimed_ddN_precursor_cells, GLRA_PCA_COORDINATES, GLRA_NUM_PCS, pseudotimes)

GLRA_sibling_31_segment_PC_pseudotime_df <- make_PC_pseudotime_df(pseudotimed_GLRA_sibling_31_segment_cells, GLRA_PCA_COORDINATES, GLRA_NUM_PCS, pseudotimes)

GLRA_sibling_45_segment_PC_pseudotime_df <- make_PC_pseudotime_df(pseudotimed_GLRA_sibling_45_segment_cells, GLRA_PCA_COORDINATES, GLRA_NUM_PCS, pseudotimes)

GLRA_mother_segment_PC_pseudotime_df <- make_PC_pseudotime_df(pseudotimed_GLRA_mother_segment_cells, GLRA_PCA_COORDINATES, GLRA_NUM_PCS, pseudotimes)
```

## Export PC_pseudotime_dfs

``` r
PC_pseudotime_dfs_out_dir <- "./PC_pseudotime_dfs/"

write.table(ddN_PC_pseudotime_df, file = paste0(PC_pseudotime_dfs_out_dir, "ddN_PC_pseudotime_df.csv"), sep = "\t", row.names = FALSE)

write.table(GLRA_sibling_31_segment_PC_pseudotime_df, file = paste0(PC_pseudotime_dfs_out_dir, "GLRA_sibling_31_segment_PC_pseudotime_df.csv"), sep = "\t", row.names = FALSE)

write.table(GLRA_sibling_45_segment_PC_pseudotime_df, file = paste0(PC_pseudotime_dfs_out_dir, "GLRA_sibling_45_segment_PC_pseudotime_df.csv"), sep = "\t", row.names = FALSE)

write.table(GLRA_mother_segment_PC_pseudotime_df, file = paste0(PC_pseudotime_dfs_out_dir, "GLRA_mother_segment_PC_pseudotime_df.csv"), sep = "\t", row.names = FALSE)
```

## Load “get_pseudostates.py”

“get_pseudostates.py” consists of modules that, given a
PC_pseudotime_dict of cells on a segment, finds a cutoff that demarcates
the pseudostate on that segment. That is, all the cells on a segment
until that cutoff are considered a part of the pseudostate.

Refer to Code documentation and Methods of paper for details on
implementation details.

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

## ddN Precursor

``` r
set.seed(0)
# Took ~3s to run.
ddN_df_file <- "./PC_pseudotime_dfs/ddN_PC_pseudotime_df.csv"

ddN_precursor_segment_cutoff_output = get_pseudostates_pipeline(ddN_df_file, n_leiden_iterations, n_routine_iterations, verbose)
```

``` r
set.seed(0)
ddN_precursor_pseudostate_cutoff <- ddN_precursor_segment_cutoff_output[[2]]
# 9
ddN_precursor_pseudostate_cells <- pseudotimed_ddN_precursor_cells[1:ddN_precursor_pseudostate_cutoff]
```

## GLRA Precursor Sibling 31

``` r
set.seed(0)
# Took ~1.5 mins to run.
GLRA_precursor_sibling_31_segment_df_file <- "./PC_pseudotime_dfs/GLRA_sibling_31_segment_PC_pseudotime_df.csv"
GLRA_precursor_sibling_31_segment_cutoff_output <- get_pseudostates_pipeline(GLRA_precursor_sibling_31_segment_df_file, n_leiden_iterations, n_routine_iterations, verbose)
```

``` r
set.seed(0)
GLRA_precursor_sibling_31_pseudostate_cutoff <- GLRA_precursor_sibling_31_segment_cutoff_output[[2]]
# 80
GLRA_precursor_sibling_31_pseudostate_cells <- pseudotimed_GLRA_sibling_31_segment_cells[1:GLRA_precursor_sibling_31_pseudostate_cutoff]
```

## GLRA Precursor Sibling 45

``` r
set.seed(0)
# Took ~5s to run.
GLRA_precursor_sibling_45_segment_df_file <- "./PC_pseudotime_dfs/GLRA_sibling_45_segment_PC_pseudotime_df.csv"
GLRA_precursor_sibling_45_segment_cutoff_output <- get_pseudostates_pipeline(GLRA_precursor_sibling_45_segment_df_file, n_leiden_iterations, n_routine_iterations, verbose)
```

``` r
set.seed(0)
GLRA_precursor_sibling_45_pseudostate_cutoff <- GLRA_precursor_sibling_45_segment_cutoff_output[[2]]
# 17
GLRA_precursor_sibling_45_pseudostate_cells <- pseudotimed_GLRA_sibling_45_segment_cells[1:GLRA_precursor_sibling_45_pseudostate_cutoff]
```

``` r
# 80/139 + 17/27 = 97/166
GLRA_precursor_sibling_pseudostate_cells <- c(GLRA_precursor_sibling_31_pseudostate_cells, GLRA_precursor_sibling_45_pseudostate_cells)
```

## GLRA Mother

``` r
set.seed(0)
# Took ~4.5 mins to run.
GLRA_mother_segment_df_file <- "./PC_pseudotime_dfs/GLRA_mother_segment_PC_pseudotime_df.csv"
GLRA_mother_segment_cutoff_output <- get_pseudostates_pipeline(GLRA_mother_segment_df_file, n_leiden_iterations, n_routine_iterations, verbose)
```

``` r
set.seed(0)
GLRA_mother_pseudostate_cutoff <- GLRA_mother_segment_cutoff_output[[2]]
# 171
GLRA_mother_pseudostate_cells <- pseudotimed_GLRA_mother_segment_cells[1:GLRA_mother_pseudostate_cutoff]
```

## Mother’s Siblings

The Mother’s Siblings must be comparable to it in pseudotime. Thus, we
use the “get_mother_sibling_pseudostates” function to find their
pseudostates. We combine them into one pseudostate.

``` r
# 12
GLRA_mother_sibling_64_output <- get_mother_sibling_pseudostates(GLRA_mother_pseudostate_cells, pseudotimed_GLRA_mother_sibling_64_segment_cells)

mother_sibling_64_first_index <- GLRA_mother_sibling_64_output$first_index
mother_sibling_64_last_index <- GLRA_mother_sibling_64_output$last_index

# 15
GLRA_mother_sibling_15_output <- get_mother_sibling_pseudostates(GLRA_mother_pseudostate_cells, pseudotimed_GLRA_mother_sibling_15_segment_cells)

mother_sibling_15_first_index <- GLRA_mother_sibling_15_output$first_index
mother_sibling_15_last_index <- GLRA_mother_sibling_15_output$last_index

# 27
GLRA_mother_sibling_pseudostate_cells <- c(pseudotimed_GLRA_mother_sibling_64_segment_cells[mother_sibling_64_first_index:mother_sibling_64_last_index],              pseudotimed_GLRA_mother_sibling_15_segment_cells[mother_sibling_15_first_index:mother_sibling_15_last_index])
```

## Visualize Pseudostates on Virtual Neural Lineage Tree

``` r
GLRA_pseudostates_list <- list("ddN_precursor_pseudostate_cells" = ddN_precursor_pseudostate_cells,
                              "GLRA_precursor_sibling_pseudostate_cells" = GLRA_precursor_sibling_pseudostate_cells,
                              "GLRA_mother_pseudostate_cells" = GLRA_mother_pseudostate_cells,
                              "GLRA_mother_sibling_pseudostate_cells" = GLRA_mother_sibling_pseudostate_cells)

modified_URD <- add_segment_state(GLRA_pseudostates_list, "GLRA_pseudostates_list")
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
# 595 x 9
ddN_precursor_pseudostate_TF_raw_mat <- URD_raw_data[tfs_in_URD, ddN_precursor_pseudostate_cells]
# 97
GLRA_precursor_sibling_pseudostate_TF_raw_mat <- URD_raw_data[tfs_in_URD, GLRA_precursor_sibling_pseudostate_cells]
# 171
GLRA_mother_pseudostate_TF_raw_mat <- URD_raw_data[tfs_in_URD, GLRA_mother_pseudostate_cells]
# 27
GLRA_mother_sibling_pseudostate_TF_raw_mat <- URD_raw_data[tfs_in_URD, GLRA_mother_sibling_pseudostate_cells]

# Raw Matrices with expression levels of all genes.
# 14809 x 9
ddN_precursor_pseudostate_raw_mat <- URD_raw_data[, ddN_precursor_pseudostate_cells]
# 97
GLRA_precursor_sibling_pseudostate_raw_mat <- URD_raw_data[, GLRA_precursor_sibling_pseudostate_cells]
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
ddN_precursor_v_sibling_pseudostate_mat <- as.matrix(cbind(ddN_precursor_pseudostate_raw_mat, GLRA_precursor_sibling_pseudostate_raw_mat))
```

``` r
ddN_precursor_v_sibling_pseudostate_edgeR_outputs <- complete_edgeR_pipeline(ddN_precursor_v_sibling_pseudostate_mat, ncol(ddN_precursor_pseudostate_raw_mat), ncol(GLRA_precursor_sibling_pseudostate_raw_mat), 0.05, -1, 50)
```

Get the DEGs.

``` r
# 196
ddN_precursor_pseudostate_DEGs_df <- ddN_precursor_v_sibling_pseudostate_edgeR_outputs[[1]]
ddN_precursor_pseudostate_DEGs <- rownames(ddN_precursor_pseudostate_DEGs_df)
```

#### Write out Precursor DEGs (Precursor Identity Genes)

``` r
write_csv_helper(ddN_precursor_pseudostate_DEGs, paste0(DEG_lists_out_dir, "ddN_precursor_pseudostate_DEGs.csv"))
```

### Precursor DEG TFs

``` r
ddN_precursor_v_sibling_pseudostate_TF_mat <- as.matrix(cbind(ddN_precursor_pseudostate_TF_raw_mat, GLRA_precursor_sibling_pseudostate_TF_raw_mat))

ddN_precursor_v_sibling_pseudostate_TF_edgeR_outputs <- complete_edgeR_pipeline(ddN_precursor_v_sibling_pseudostate_TF_mat, ncol(ddN_precursor_pseudostate_TF_raw_mat), ncol(GLRA_precursor_sibling_pseudostate_TF_raw_mat), 0.05, -1, 50)
```

``` r
# 16
ddN_precursor_v_sibling_pseudostate_TF_DEGs_df <- ddN_precursor_v_sibling_pseudostate_TF_edgeR_outputs[[1]]
# Get TF orthologs.
named_ddN_precursor_v_sibling_pseudostate_TF_DEGs_df <- convert_khids_to_ortholog(ddN_precursor_v_sibling_pseudostate_TF_DEGs_df)

ddN_precursor_v_sibling_pseudostate_TF_DEGs <- named_ddN_precursor_v_sibling_pseudostate_TF_DEGs_df$Ghost.Name
```

#### Write out TF DEGs

``` r
write_csv_helper(ddN_precursor_v_sibling_pseudostate_TF_DEGs, paste0(DEG_lists_out_dir, "ddN_precursor_v_sibling_pseudostate_TF_DEGs.csv"))
```

## Precursor vs. Mother- DEG TFs

``` r
ddN_precursor_v_mother_pseudostate_TF_mat <- as.matrix(cbind(ddN_precursor_pseudostate_TF_raw_mat, GLRA_mother_pseudostate_TF_raw_mat))

ddN_precursor_v_mother_pseudostate_TF_edgeR_outputs <- complete_edgeR_pipeline(ddN_precursor_v_mother_pseudostate_TF_mat, ncol(ddN_precursor_pseudostate_TF_raw_mat), ncol(GLRA_mother_pseudostate_TF_raw_mat), 0.05, -1, 50)
```

``` r
# 12
ddN_precursor_v_mother_pseudostate_TF_DEGs_df <- ddN_precursor_v_mother_pseudostate_TF_edgeR_outputs[[1]]
# Get TF orthologs.
named_ddN_precursor_v_mother_pseudostate_TF_DEGs_df <- convert_khids_to_ortholog(ddN_precursor_v_mother_pseudostate_TF_DEGs_df)

ddN_precursor_v_mother_pseudostate_TF_DEGs <- named_ddN_precursor_v_mother_pseudostate_TF_DEGs_df$Ghost.Name
```

#### Write out TF DEGs

``` r
write_csv_helper(ddN_precursor_v_mother_pseudostate_TF_DEGs, paste0(DEG_lists_out_dir, "ddN_precursor_v_mother_pseudostate_TF_DEGs.csv"))
```

## Mother vs. Mother’s Sibling- DEG TFs

``` r
GLRA_mother_v_sibling_pseudostate_TF_mat <- as.matrix(cbind(GLRA_mother_pseudostate_TF_raw_mat, GLRA_mother_sibling_pseudostate_TF_raw_mat))

GLRA_mother_v_sibling_pseudostate_TF_edgeR_outputs <- complete_edgeR_pipeline(GLRA_mother_v_sibling_pseudostate_TF_mat, ncol(GLRA_mother_pseudostate_TF_raw_mat), ncol(GLRA_mother_sibling_pseudostate_TF_raw_mat), 0.05, -1, 50)
```

``` r
# 17
GLRA_mother_v_sibling_pseudostate_TF_DEGs_df <- GLRA_mother_v_sibling_pseudostate_TF_edgeR_outputs[[1]]
# Get TF orthologs.
named_GLRA_mother_v_sibling_pseudostate_TF_DEGs_df <- convert_khids_to_ortholog(GLRA_mother_v_sibling_pseudostate_TF_DEGs_df)

GLRA_mother_v_sibling_pseudostate_TF_DEGs <- named_GLRA_mother_v_sibling_pseudostate_TF_DEGs_df$Ghost.Name
```

#### Write out TF DEGs

``` r
write_csv_helper(GLRA_mother_v_sibling_pseudostate_TF_DEGs, paste0(DEG_lists_out_dir, "GLRA_mother_v_sibling_pseudostate_TF_DEGs.csv"))
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
ddN_precursor_pseudostate_mean_TF_exp_df <- make_named_mean_TF_exp_df(ddN_precursor_pseudostate_cells)

GLRA_mother_pseudostate_mean_TF_exp_df <- make_named_mean_TF_exp_df(GLRA_mother_pseudostate_cells)
```

## Take the union of top 50% most highly expressed Unfiltered Candidate Cocktail TFs in Precursor and Mother. These are the Candidate Cocktail TFs

######## REF SUPP FIGURE PG 30

``` r
# 33
ddN_unfiltered_candidate_cocktail_tfs <- unique(c(named_ddN_precursor_v_sibling_pseudostate_TF_DEGs_df$KH.gene.model,
                                                  named_ddN_precursor_v_mother_pseudostate_TF_DEGs_df$KH.gene.model,
                                                  named_GLRA_mother_v_sibling_pseudostate_TF_DEGs_df$KH.gene.model))

# Take 50% cutoff.
# This returns both a union and intersection of top 50% most highly expressed unfiltered candidate cocktail TFs.
ddN_cutoff_df_list <- get_cutoff_df(ddN_precursor_pseudostate_mean_TF_exp_df, GLRA_mother_pseudostate_mean_TF_exp_df, ddN_unfiltered_candidate_cocktail_tfs, 0.5)

# Let us take the union.
# 25
ddN_cutoff_union <- ddN_cutoff_df_list$union
# These are the Candidate Cocktail TFs.
ddN_candidate_cocktail_tfs <- unique(c(ddN_cutoff_union$Ghost.Name_precursor, ddN_cutoff_union$Ghost.Name_mother))
# Remove any NA's!
ddN_candidate_cocktail_tfs <- ddN_candidate_cocktail_tfs[!is.na(ddN_candidate_cocktail_tfs)]

# We also want their khids.
ddN_candidate_cocktail_khids <- ddN_cutoff_union$KH.gene.model
```

### Write out the Candidate Cocktail TFs

``` r
write_csv_helper(ddN_candidate_cocktail_tfs, paste0(DEG_lists_out_dir, "ddN_candidate_cocktail_tfs.csv"))
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
# 11
common_ddN_DEG_TFs <- intersect(ddN_candidate_cocktail_khids, ddN_precursor_pseudostate_DEGs)
# 185
cleaned_ddN_precursor_pseudostate_DEGs <- setdiff(ddN_precursor_pseudostate_DEGs, common_ddN_DEG_TFs)
```

## Get numbers of Precursor cells (n), Candidate Cocktail TFs (p), and Precusor DEGs (k)

``` r
# 9, 25, 185
ddN_n <- length(ddN_precursor_pseudostate_cells)
ddN_p <- length(ddN_candidate_cocktail_tfs)
ddN_k <- length(cleaned_ddN_precursor_pseudostate_DEGs)
```

### Make Matrix S

Matrix S is the same as ddN_TFs_in_precursor_pseudostate_mat below.

``` r
# 9 x 25
ddN_TFs_in_precursor_pseudostate_mat <- t(URD_data[ddN_candidate_cocktail_khids, ddN_precursor_pseudostate_cells])
```

### Concatenate k S Matrices to make Group Lasso Matrix X

``` r
ddN_S_matrix <- ddN_TFs_in_precursor_pseudostate_mat
ddN_zero_matrix <- matrix(0, ddN_n, ddN_p)
# This is a dummy matrix- we will remove it at the end.
ddN_grplasso_X_mat <- matrix(1, ddN_n, ddN_p*ddN_k)

for (i in 1:ddN_k) {
  ddN_L_matrix <- do.call(cbind, replicate(i-1, ddN_zero_matrix, simplify=FALSE))
  ddN_R_matrix <- do.call(cbind, replicate(ddN_k-i, ddN_zero_matrix, simplify=FALSE))
  # At each iteration, we create a row in which Matrix S is sandwiched between 
  # Left "L" Zero Matrix(ces) and Right "R" Zero Matrix(ces) given its position on
  # the diagonal of Group Lasso Matrix X.
  ddN_grplasso_mat_row <- do.call(cbind, list(ddN_L_matrix, ddN_S_matrix, ddN_R_matrix))
  ddN_grplasso_X_mat <- rbind(ddN_grplasso_X_mat, ddN_grplasso_mat_row)
}

# 1665 x 4625
ddN_grplasso_X_mat_2 <- ddN_grplasso_X_mat[-c(1:ddN_n), ]
```

### Rename Row and Column names for Convenience

#### Rename Row Names to reflect cell number and its DEG status

The Row Names have the following pattern: the first n rows are
DEG_1\_cell_1… DEG_1\_cell_n. The following n rows are DEG_2\_cell_1…
DEG_2\_cell_n for all k DEGs.

``` r
ddN_grplasso_X_mat_row_cells <- rep(paste0("_cell_", seq(1, ddN_n)), ddN_k)
ddN_DEG_statuses <- rep(c(1:ddN_k), each = ddN_n)
ddN_new_grplasso_X_mat_2_rownames <- paste0("DEG_", ddN_DEG_statuses, ddN_grplasso_X_mat_row_cells)

rownames(ddN_grplasso_X_mat_2) <- ddN_new_grplasso_X_mat_2_rownames
```

#### Rename Column Names to reflect TF number and its DEG status

Similarly, the Column Names have the following pattern: the first p rows
are DEG_1\_TF1… DEG_1\_TFp. The following p rows are DEG_2\_TF1…
DEG_2\_TFp for all p TFs.

``` r
ddN_grplasso_X_mat_col_TFs <- rep(paste0("_TF", seq(1, ddN_p)), ddN_k)
ddN_DEG_statuses_cols <- rep(c(1:ddN_k), each = ddN_p)
ddN_new_grplasso_X_mat_2_colnames <- paste0("DEG_", ddN_DEG_statuses_cols, ddN_grplasso_X_mat_col_TFs)

colnames(ddN_grplasso_X_mat_2) <- ddN_new_grplasso_X_mat_2_colnames
```

### Final Group Lasso Matrix X

``` r
# 1665 x 4625
final_grplasso_X_mat <- as.matrix(ddN_grplasso_X_mat_2)
```

## Make Group Lasso Matrix Y

Group Lasso Matrix Y contains the expression of all Precursor DEGs in
all Precursor Cells in the following way: The first n rows contain
expression of DEG 1 in all n cells; the next n rows contain expression
of DEG 2 in all n cells; so and so forth for k DEGs. It has dimensions
nk x 1.

### Make a matrix that has the expression of all Precursor DEGs in Precursor Cells of dimensions k x n

``` r
# 9 x 185
ddN_precursor_DEGs_in_precursors_mat <- t(URD_data[cleaned_ddN_precursor_pseudostate_DEGs, ddN_precursor_pseudostate_cells])
```

### Index above matrix so as so make Group Lasso Matrix Y

``` r
ddN_grplasso_Y_vec <- c()

for (DEG in 1:ddN_k) {
  for (cell in 1:ddN_n) {
    ddN_DEG_cell_exp <- ddN_precursor_DEGs_in_precursors_mat[cell, DEG]
    ddN_grplasso_Y_vec <- c(ddN_grplasso_Y_vec, ddN_DEG_cell_exp)
  }
}

# 1665 x 1
ddN_grplasso_Y_mat <- as.matrix(ddN_grplasso_Y_vec, ncol=1, nrow=(ddN_n*ddN_k))
```

### Rename Row names to reflect cell and DEG status

The Row Names follow the same pattern as the Row Names of Group Lasso
Matrix X.

``` r
ddN_grplasso_Y_DEG_statuses <- rep(c(1:ddN_k), each = ddN_n)
ddN_grplasso_Y_DEGs <- paste0("DEG_", ddN_grplasso_Y_DEG_statuses)
ddN_grplasso_Y_cells <- rep(paste0("_cell_", seq(1, ddN_n)), ddN_k)
ddN_new_grplasso_Y_rownames <- paste0(ddN_grplasso_Y_DEGs, ddN_grplasso_Y_cells)

rownames(ddN_grplasso_Y_mat) <- ddN_new_grplasso_Y_rownames
```

### Final Matrix Y

``` r
# 1665 x 1
final_ddN_grplasso_Y_mat <- ddN_grplasso_Y_mat
```

## Run Group Lasso

### Variables: Group Lasso Matrices X and Y and Group IDs

Each TF is a group. Thus, there are p groups each of size k. Given the
way Group Lasso Matrix X is organized, the Group IDs follow the pattern
of 1, 2, … p; 1, 2, … p; k times.

``` r
grplasso_ddN_groups <- rep(seq(1, ddN_p), ddN_k)

ddN_X <- final_grplasso_X_mat
ddN_Y <- final_ddN_grplasso_Y_mat
ddN_group_id <- grplasso_ddN_groups
```

### Do Cross-Validation

``` r
# ~3 mins to Run 
grplasso_ddN_cv_grepreg_fit <- cv.grpreg(ddN_X, ddN_Y, ddN_group_id, penalty="grLasso")
```

### Get Sets of Selected TF Groups at Various Values of Lambda

######## REF SUPP FIGURE PG 31

``` r
# Key: lambda; Values: Set of Selected TF Groups
ddN_group_lassoed_tf_list <- get_group_lassoed_tf_list(grplasso_ddN_cv_grepreg_fit, ddN_candidate_cocktail_tfs)

# Key: lambda; Values: Set of Non-selected TF Groups
ddN_non_group_lassoed_TF_list <- lapply(ddN_group_lassoed_tf_list, get_non_group_lassoed_TFs, ddN_candidate_cocktail_tfs)
```

#### Let us Visualize the Number of Groups Selected at Various Values of Lambda

``` r
plot(grplasso_GLRA_cv_grepreg_fit)
```

### Get Sorted Total TF Importance Lists

The Total Importance Score of a Cocktail TF is the sum of the absolute
value of its DEG-specific coefficients. It thus gives a quantitative
measure of a Cocktail TF’s strength over the Gene Regulatory Network (as
defined by the Precursor DEGs) underlying the Precursor’s state.

``` r
# (For all lambdas) key = lambda; value is a list with key = TF and value = non0 betas.
ddN_TF_importances_lists <- get_TF_importances_list(grplasso_ddN_cv_grepreg_fit)

# Sum the absolute values of these betas for all TFs for all lambdas.
ddN_sum_total_TF_importance_lists <- ddN_TF_importances_lists %>% map(get_total_summed_TF_importance_list)

# Make sorted dfs given these Total Importance Scores.
ddN_sorted_sum_total_TF_importance_lists <- ddN_sum_total_TF_importance_lists %>% map(make_sorted_TF_importance_df, ddN_candidate_cocktail_tfs)
```

# Final Result: Inferred ddN Cocktail

Upon Inspection of “ddN_group_lassoed_tf_list”, we see that the set of
3, the Inferred Cocktail, consists of Dmbx, Lhx1, and HMG1/2.

# Extended Analyses

> Extended Analyses 1 and 2 compare Cocktail Inference results from
> another Linear Model, Gaussian Multivariate Multi-Response Linear
> Regression (GMML), and a Nonlinear Model based on GRNBoost2
> respectively to both demonstrate the superiority of the Group Lasso
> and validate its results.

> Extended Analysis 3 compares Cocktail Inference results with Terminal
> ddNs (fully differentiated ddNs) to demonstrate the need to study
> regulatory changes at the time of specification in order to correctly
> infer Cocktails.

> Extended Analysis 4 performs subsampling of ddN Precursor data to
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
glm_ddN_X_mat <- as.matrix(ddN_TFs_in_precursor_pseudostate_mat)
glm_ddN_Y_mat <- as.matrix(ddN_precursor_DEGs_in_precursors_mat)
```

### Perform Regression

#### Without Lasso penalty

We will refer to this as the complete GMML model.

``` r
# ~Took 5s to run.
ddN_fit <- cv.glmnet(glm_ddN_X_mat, glm_ddN_Y_mat, intercept = FALSE, family = "mgaussian")
ddN_lambda <- ddN_fit$lambda.min
# ddN_coef is a list with keys = khids of Precursor DEGs and values = series of attributes
# one of which (@ x) contains weights for all 14 input TF variables.
ddN_coef_list <- coef(ddN_fit, s = ddN_lambda)
```

#### With Lasso penalty, selecting for 3 TF Variables

Here, we set the dfmax parameter to 2, which applies a Lasso penalty and
selects for 3 TF variables.

One can observe the results given different values of this parameter,
such as 3 or 4, which respectively returns 4 or 5 variables, but we
study 3 variables in order to compare this result to our inferred
cocktail, which contains 3 TFs.

``` r
# Took ~15s to run.
ddN_fit_3TFs <- cv.glmnet(glm_ddN_X_mat, glm_ddN_Y_mat, intercept = FALSE, family = "mgaussian", dfmax = 2)
ddN_lambda_3TFs <- ddN_fit_3TFs$lambda.min
ddN_coef_3TFs_list <- coef(ddN_fit_3TFs, s = ddN_lambda_3TFs)
```

### Strategy 1

#### Get Regression Coefficients for all Candidate Cocktail TFs

Make a coef_df, which is a data frame with rows = Candidate Cocktail
Precursor DEGs, columns = Precursor DEGs, and entries as regression
coefficients present in coef_list.

``` r
ddN_coef_df <- make_coef_df(ddN_coef_list,  ddN_candidate_cocktail_tfs, cleaned_ddN_precursor_pseudostate_DEGs)

ddN_coef_df_3TFs <- make_coef_df(ddN_coef_3TFs_list,  ddN_candidate_cocktail_tfs, cleaned_ddN_precursor_pseudostate_DEGs)
```

### Visualize Distribution of Coefficients and get their Averages

We notice that Dmbx, Lhx1, and HMG1/2 have the highest 1st, 2nd, and 4th
average coefficient values in the complete GMML model.

######## REF SUPP FIGURE PG 32, 33

``` r
ddN_average_coef_df <- get_average_coefficient_df(ddN_coef_df)
ddN_coef_distribution_plot <- get_coefficients_distribution_plot(ddN_coef_df)
ddN_coef_distribution_plot

ddN_average_coef_df_3TFs <- get_average_coefficient_df(ddN_coef_df_3TFs)
```

#### Get Total Importance Scores for Candidate Cocktail TFs

> Sum the absolute values of the regression coefficients

######## REF SUPP FIGURE PG 34

``` r
ddN_total_GMML_importance_df <- get_total_GMML_importance_df(ddN_coef_df)
ddN_total_GMML_importance_df_3TFs <- get_total_GMML_importance_df(ddN_coef_df_3TFs)
```

### Strategy 2

#### First Get the ranks of the absolute values of the regression coefficients

``` r
ranked_ddN_coef_df <- ddN_coef_df %>% mutate(across(everything(), abs), across(everything(), min_rank))
ranked_ddN_coef_df_3TFs <- ddN_coef_df_3TFs %>% mutate(across(everything(), abs), across(everything(), min_rank))
```

#### Visualize Distribution of Coefficients and get their Averages

We notice that Dmbx, Lhx1, and HMG1/2 have the highest 1st, 2nd, and 4th
average ranked coefficient values in the complete GMML model.

######## REF SUPP FIGURE PG 32, 34

``` r
ranked_ddN_average_coef_df <- get_average_coefficient_df(ranked_ddN_coef_df)
ranked_ddN_coef_distribution_plot <- get_coefficients_distribution_plot(ranked_ddN_coef_df)
ranked_ddN_coef_distribution_plot

ranked_ddN_average_coef_df_3TFs <- get_average_coefficient_df(ranked_ddN_coef_df_3TFs)
```

#### Get Total Importance Scores for Candidate Cocktail TFs

######## REF SUPP FIGURE PG 35

``` r
ranked_ddN_total_GMML_importance_df <- get_total_GMML_importance_df(ranked_ddN_coef_df)

ranked_ddN_total_GMML_importance_df_3TFs <- get_total_GMML_importance_df(ranked_ddN_coef_df_3TFs)
```

### Notes on Final Result

> In this case, given the above results, the ddN Cocktail is likely
> inferred from the GMML model.

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
ddN_GRN_dir <- "GRN/ddN/"
```

#### ddN Precursor DEG Counts Matrix

``` r
# Matrix must have counts from both DEGs and Candidate Cocktail TFs.
ddN_precursor_pseudostate_DEGs_plus <- unique(c(ddN_precursor_pseudostate_DEGs, ddN_candidate_cocktail_khids))
ddN_precursor_pseudostate_GRN_mat <- URD_raw_data[ddN_precursor_pseudostate_DEGs_plus, ddN_precursor_pseudostate_cells]
# 9 x 210
ddN_precursor_pseudostate_GRN_df <- as.data.frame(t(ddN_precursor_pseudostate_GRN_mat))
```

#### Inputs to GRNBoost2 Pipeline

``` r
write_tsv(ddN_precursor_pseudostate_GRN_df, paste0("GRN/", "ddN_precursor_pseudostate_GRN_df.tsv"))
```

``` r
write_csv_helper(ddN_candidate_cocktail_khids, paste0("GRN/", "ddN_candidate_cocktail_khids.csv"))
```

#### Do GRNBoost2 50 times.

``` r
source_python("GRNBoost2_pipeline.py")
```

``` r
set.seed(16)

expression_count_matrix_tsv <- paste0("GRN/", "ddN_precursor_pseudostate_GRN_df.tsv")
tf_khid_csv <- paste0("GRN/", "ddN_candidate_cocktail_khids.csv")
outname <- "GRN/ddN/ddN"
min_i <- 0L
max_i <- 50L

# ~Took 15 mins to run.
do_GRNBoost2(expression_count_matrix_tsv, tf_khid_csv, outname, min_i, max_i)
```

##### sorted_total_importance_tfs

Each element of the list is a tuple of TF and list with Mean and
Standard Deviation of Total Importance Scores. The TFs are sorted in
nonincreasing order based on the mean.

##### sorted_percent_top_list

Each element of the list is a TF and its Percentage Covered as a Top 3
Regulator for all DEGs.

``` r
khid_tf_tsv_file <- "/Users/prakritipaul/Git/Cocktails/useful_lists/kai_tf.tsv"
ddN_GRN_path <- "/Users/prakritipaul/Git/Cocktails/FINAL_code_and_results/GRN/ddN/*.tsv"
n <- 3L

# Returns sorted_total_importance_tfs and sorted_percent_top_list
ddN_GRNBoost2_outputs <- GRNBoost2_pipeline(ddN_GRN_path, khid_tf_tsv_file, n)
```

##### Inspect Output contents.

``` r
ddN_sorted_total_importance_tfs <- ddN_GRNBoost2_outputs[1]
ddN_sorted_percent_top_list <- ddN_GRNBoost2_outputs[2]
```

### Notes on Final Result

> Upon inspection of the sorted_total_importance_tfs and
> sorted_percent_top-list, only Dmbx and Lhx1 appear in the Top 3. In
> this case, this method is not sufficient. However, this result does
> further support the results obtained from the Group Lasso Model.

## (3) Comparison with Method Applied on Terminal ddNs

A necessary comparison is with Terminal ddNs, or fully differentiated
ddNs, as a key claim in this work is that the Cocktail cannot be
inferred without studying cell states at the time of specification.

Potentially, we can find the Cocktail either in the top 50% mostly
highly expressed DEG TFs in the Terminal DEGs (these are Terminal
Candidate Cocktail TFs) or in the output of running the Group Lasso
Model using Terminal DEGs and Terminal Candidate Cocktail Factors. We
will also use the GMML.

From the following analyses, we demonstrate that the Cocktail cannot be
inferred from either way and build evidence for the power of our
computational method and the strength of the assumptions we made.

### Identify Terminal ddNs and their Siblings

Siblings were identified from observing the Virtual Lineage Tree.

``` r
# 64
terminal_GLRAs <- neural_CT_list$`GLRA1.2.3+ motor ganglion`
# 18
terminal_ddNs <- get_dmbx_expressing_cells(terminal_GLRAs, pseudotimes)
# 117
terminal_ddN_siblings <- c(neural_CT_list$`KCNB1+ motor ganglion`, neural_CT_list$`AMD+ motor ganglion` )
```

### Raw Mats

``` r
terminal_ddN_raw_mat <- URD_raw_data[, terminal_ddNs]
terminal_ddN_TF_raw_mat <- URD_raw_data[tfs_in_URD, terminal_ddNs]

terminal_ddN_sibling_raw_mat <- URD_raw_data[, terminal_ddN_siblings]
terminal_ddN_sibling_TF_raw_mat <- URD_raw_data[tfs_in_URD, terminal_ddN_siblings]
```

### Calculate Terminal DEGs with respect to Siblings

Terminal DEGs are analogous to Precursor DEGs.

We calculate Terminal DEG TFs with respect to Siblings because
comparison to all other neural cell types did not render meaningful TFs.
Instead, TFs known to play roles in specification appeared as
differentially expressed when Terminal ddNs were compared to its
siblings. This made these results more competitive with the results
obtained from applying the computational method to Precursors at the
time of specification.

### Get Terminal DEGs

``` r
terminal_ddN_v_siblings_mat <- as.matrix(cbind(terminal_ddN_raw_mat, terminal_ddN_sibling_raw_mat))

terminal_ddN_v_siblings_edgeR_outputs <- complete_edgeR_pipeline(terminal_ddN_v_siblings_mat, ncol(terminal_ddN_raw_mat), ncol(terminal_ddN_sibling_raw_mat), 0.05, -1, 50)
```

``` r
# 297
terminal_ddN_v_siblings_DEGs_df <- terminal_ddN_v_siblings_edgeR_outputs[[1]]
terminal_ddN_v_siblings_DEGs <- rownames(terminal_ddN_v_siblings_DEGs_df)
```

#### Write them out.

``` r
write_csv_helper(terminal_ddN_v_siblings_DEGs, paste0(DEG_lists_out_dir, "terminal_ddN_v_siblings_DEGs.csv"))
```

### Get Terminal DEG TFs

We will filter these downstream to get Terminal Candidate Cocktail TFs.

``` r
terminal_ddN_v_siblings_TF_mat <- as.matrix(cbind(terminal_ddN_TF_raw_mat, terminal_ddN_sibling_TF_raw_mat))

terminal_ddN_v_siblings_TF_edgeR_outputs <- complete_edgeR_pipeline(terminal_ddN_v_siblings_TF_mat, ncol(terminal_ddN_TF_raw_mat), ncol(terminal_ddN_sibling_TF_raw_mat), 0.05, -1, 50)
```

``` r
# 14
terminal_ddN_v_siblings_TF_DEGs_df <- terminal_ddN_v_siblings_TF_edgeR_outputs[[1]]
# Get TF orthologs.
named_terminal_ddN_v_siblings_TF_DEGs_df <- convert_khids_to_ortholog(terminal_ddN_v_siblings_TF_DEGs_df)
terminal_ddN_v_sibling_TF_DEGs <- named_terminal_ddN_v_siblings_TF_DEGs_df$Ghost.Name
```

#### Write them out.

``` r
write_csv_helper(terminal_ddN_v_sibling_TF_DEGs, paste0(DEG_lists_out_dir, "terminal_ddN_v_sibling_TF_DEGs.csv"))
```

### Get Terminal Candidate Cocktail TFs

As before, we will take the top 50% most highly expressed of the
Terminal DEG TFs and consider them as the Terminal Candidate Cocktail
TFs.

######## REF SUPP FIGURE PG 35

``` r
# Get mean Expression of all DEG TFs
terminal_ddN_mean_TF_exp_df <- make_named_mean_TF_exp_df(terminal_ddNs)

# Get khids of Terminal DEG TFs
terminal_ddN_v_siblings_TF_DEG_khids <- rownames(terminal_ddN_v_siblings_TF_DEGs_df)

# Get Mean expression of Terminal DEG TFs in terminal cells.
terminal_ddN_cocktail_TF_mean_exp_df <- terminal_ddN_mean_TF_exp_df %>% filter(KH.gene.model %in% terminal_ddN_v_siblings_TF_DEG_khids)

# Get top 50% most highly expressed of above DEG TFS
terminal_cutoff <- 0.5

terminal_ddN_cutoff_df <- terminal_ddN_cocktail_TF_mean_exp_df[1:(terminal_cutoff*(length(terminal_ddN_v_siblings_TF_DEG_khids))), ]

# Get the TFs and khids
# 7
terminal_ddN_candidate_cocktail_khids <- terminal_ddN_cutoff_df$KH.gene.model
terminal_ddN_candidate_cocktail_tfs <- terminal_ddN_cutoff_df$Ghost.Name
```

#### Write out the Candidate Cocktail TFs

``` r
write_csv_helper(terminal_ddN_candidate_cocktail_tfs, paste0(DEG_lists_out_dir, "terminal_ddN_candidate_cocktail_tfs.csv"))
```

#### Brief Aside

An interesting observation is that the common Candidate Cocktail TFs
between the Precursor and Terminal ddN are the Cocktail itself except
for HMG1/2 (Pioneer factor).

``` r
common_ddN_cocktail_TFs <- intersect(terminal_ddN_candidate_cocktail_tfs, ddN_candidate_cocktail_tfs)
common_ddN_cocktail_TFs
```

### (a) Do Cocktail TF Inference using Group Lasso Model

#### Variables

``` r
# Get rid of TFs in DEGs.
# 7
common_terminal_ddN_DEG_TFs <- intersect(terminal_ddN_candidate_cocktail_khids, terminal_ddN_v_siblings_DEGs)
# 290
cleaned_terminal_ddN_DEGs <- setdiff(terminal_ddN_v_siblings_DEGs, common_terminal_ddN_DEG_TFs)
```

``` r
# 18, 7, 290
terminal_ddN_n <- length(terminal_ddNs)
terminal_ddN_p <- length(terminal_ddN_candidate_cocktail_tfs)
terminal_ddN_k <- length(cleaned_terminal_ddN_DEGs)
```

#### Make Group Lasso X Matrix

``` r
# 18 x 7
terminal_ddN_TFs_in_terminal_mat <- t(URD_data[terminal_ddN_candidate_cocktail_khids, terminal_ddNs])
```

``` r
# 5220 x 2030
terminal_ddN_grplasso_X_mat <- make_grplasso_X_mat(terminal_ddN_TFs_in_terminal_mat, terminal_ddN_n, terminal_ddN_p, terminal_ddN_k)
```

#### Make Group Lasso Y Matrix

``` r
# 18 x 290
terminal_ddN_DEGs_in_terminal_mat <- t(URD_data[cleaned_terminal_ddN_DEGs, terminal_ddNs])
```

``` r
# 5220    1
terminal_ddN_grplasso_Y_mat <- make_grplasso_Y_mat(terminal_ddN_DEGs_in_terminal_mat, terminal_ddNs, terminal_ddN_p, terminal_ddN_k)
```

#### Variables: Group Lasso Matrices X and Y and Group IDs

``` r
terminal_grplasso_ddN_groups <- rep(seq(1, terminal_ddN_p), terminal_ddN_k)

terminal_ddN_X <- terminal_ddN_grplasso_X_mat
terminal_ddN_Y <- terminal_ddN_grplasso_Y_mat
terminal_ddN_group_id <- terminal_grplasso_ddN_groups
```

#### Error: Group Lasso does not solve.

``` r
# 
terminal_grplasso_ddN_cv_grepreg_fit <- cv.grpreg(terminal_ddN_X, terminal_ddN_Y, terminal_ddN_group_id, penalty="grLasso")
```

### Get Sets of Selected TF Groups at Various Values of Lambda

``` r
# Key: lambda; Values: Set of Selected TF Groups
terminal_ddN_group_lassoed_tf_list <- get_group_lassoed_tf_list(terminal_grplasso_ddN_cv_grepreg_fit, terminal_ddN_candidate_cocktail_tfs)

# Key: lambda; Values: Set of Non-selected TF Groups
terminal_ddN_non_group_lassoed_TF_list <- lapply(terminal_ddN_group_lassoed_tf_list, get_non_group_lassoed_TFs, terminal_ddN_candidate_cocktail_tfs)
```

#### Let us Visualize the Number of Groups Selected at Various Values of Lambda

``` r
plot(terminal_grplasso_ddN_cv_grepreg_fit)
```

### Get Sorted Total TF Importance Lists

The Total Importance Scores within these lists give us a quantitative
measure of a Cocktail TF’s strength over the Gene Regulatory Network (as
defined by the Precursor DEGs) underlying the Precursor’s state.

``` r
# (For all lambdas) key = lambda; value is a list with key = TF_num and value = betas.
terminal_ddN_TF_importances_lists <- get_TF_importances_list(terminal_grplasso_ddN_cv_grepreg_fit)

# Sum the absolute values of these betas for all TFs for all lambdas.
terminal_ddN_sum_total_TF_importance_lists <- terminal_ddN_TF_importances_lists %>% map(get_total_summed_TF_importance_list)

# Make sorted dfs given these Total Importance Scores.
terminal_ddN_sorted_sum_total_TF_importance_lists <- terminal_ddN_sum_total_TF_importance_lists %>% map(make_sorted_TF_importance_df, terminal_ddN_candidate_cocktail_tfs)
```

#### Result

The model does not converge. We cannot know if the Group Lasso was
successful.

### (b) Do GMML

``` r
glm_terminal_ddN_X_mat <- as.matrix(terminal_ddN_TFs_in_terminal_mat)
glm_terminal_ddN_Y_mat <- as.matrix(terminal_ddN_DEGs_in_terminal_mat)
```

### Perform Regression

#### Without Lasso penalty

We will refer to this as the complete GMML model.

``` r
# ~Took 40s to run.
terminal_ddN_fit <- cv.glmnet(glm_terminal_ddN_X_mat, glm_terminal_ddN_Y_mat, intercept = FALSE, family = "mgaussian")
terminal_ddN_lambda <- terminal_ddN_fit$lambda.min
# terminal_ddN_coef is a list with keys = khids of Precursor DEGs and values = series of attributes
# one of which (@ x) contains weights for all 14 input TF variables.
terminal_ddN_coef_list <- coef(terminal_ddN_fit, s = terminal_ddN_lambda)
```

#### With Lasso penalty, selecting for 3 TF Variables

Here, we set the dfmax parameter to 2, which applies a Lasso penalty and
selects for 3 TF variables.

One can observe the results given different values of this parameter,
such as 3 or 4, which respectively returns 4 or 5 variables, but we
study 3 variables in order to compare this result to our inferred
cocktail, which contains 3 TFs.

``` r
# Took ~30s to run.
terminal_ddN_fit_3TFs <- cv.glmnet(glm_terminal_ddN_X_mat, glm_terminal_ddN_Y_mat, intercept = FALSE, family = "mgaussian", dfmax = 2)
terminal_ddN_lambda_3TFs <- terminal_ddN_fit_3TFs$lambda.min
terminal_ddN_coef_3TFs_list <- coef(terminal_ddN_fit_3TFs, s = terminal_ddN_lambda_3TFs)
```

### Strategy 1

#### Get Regression Coefficients for all Candidate Cocktail TFs

Make a coef_df, which is a data frame with rows = Candidate Cocktail
Precursor DEGs, columns = Precursor DEGs, and entries as regression
coefficients present in coef_list.

``` r
terminal_ddN_coef_df <- make_coef_df(terminal_ddN_coef_list,  terminal_ddN_candidate_cocktail_tfs, cleaned_terminal_ddN_DEGs)

terminal_ddN_coef_df_3TFs <- make_coef_df(terminal_ddN_coef_3TFs_list,  terminal_ddN_candidate_cocktail_tfs, cleaned_terminal_ddN_DEGs)
```

#### Visualize Distribution of Coefficients and get their Averages

Only Dmbx has a solved coefficient.

######## REF SUPP FIGURE PG 36, 37

``` r
terminal_ddN_average_coef_df <- get_average_coefficient_df(terminal_ddN_coef_df)
terminal_ddN_coef_distribution_plot <- get_coefficients_distribution_plot(terminal_ddN_coef_df)
terminal_ddN_coef_distribution_plot

terminal_ddN_average_coef_df_3TFs <- get_average_coefficient_df(terminal_ddN_coef_df_3TFs)
```

#### Get Total Importance Scores for Candidate Cocktail TFs

######## REF SUPP FIGURE PG 38

``` r
terminal_terminal_ddN_total_GMML_importance_df <- get_total_GMML_importance_df(terminal_ddN_coef_df)

terminal_ddN_total_GMML_importance_df_3TFs <- get_total_GMML_importance_df(terminal_ddN_coef_df_3TFs)
```

### Strategy 2

#### First Get the ranks of the absolute values of the regression coefficients

``` r
ranked_terminal_ddN_coef_df <- terminal_ddN_coef_df %>% mutate(across(everything(), abs), across(everything(), min_rank))

ranked_terminal_ddN_coef_df_3TFs <- terminal_ddN_coef_df_3TFs %>% mutate(across(everything(), abs), across(everything(), min_rank))
```

#### Visualize Distribution of Coefficients and get their Averages

Only Dmbx has a solved coefficient.

######## REF SUPP FIGURE PG 37

``` r
ranked_terminal_ddN_average_coef_df <- get_average_coefficient_df(ranked_terminal_ddN_coef_df)

ranked_terminal_ddN_coef_distribution_plot <- get_ddN_coefficients_distribution_plot(ranked_terminal_ddN_coef_df)
ranked_terminal_ddN_coef_distribution_plot

ranked_terminal_ddN_average_coef_df_3TFs <- get_average_coefficient_df(ranked_terminal_ddN_coef_df_3TFs)
```

#### Get Total Importance Scores for Candidate Cocktail TFs

######## REF SUPP FIGURE PG 38

``` r
ranked_terminal_ddN_total_GMML_importance_df <- get_total_GMML_importance_df(ranked_terminal_ddN_coef_df)

ranked_terminal_ddN_total_GMML_importance_df_3TFs <- get_total_GMML_importance_df(ranked_terminal_ddN_coef_df_3TFs)
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
terminal_ddN_GRN_dir <- "./GRN/terminal_ddN/"
```

##### Terminal ddN DEG Counts Matrix

``` r
# Matrix must have counts from both DEGs and Candidate Cocktail TFs.
terminal_ddN_DEGs_plus <- unique(c(terminal_ddN_v_siblings_DEGs, terminal_ddN_candidate_cocktail_khids))
terminal_ddN_GRN_mat <- URD_raw_data[terminal_ddN_DEGs_plus, terminal_ddNs]

# 18 x 297
terminal_ddN_GRN_df <- as.data.frame(t(terminal_ddN_GRN_mat))
```

#### Inputs to GRNBoost2 Pipeline

``` r
write_tsv(terminal_ddN_GRN_df, paste0("GRN/", "terminal_ddN_GRN_df.tsv"))
```

``` r
write_csv_helper(terminal_ddN_candidate_cocktail_khids, paste0("GRN/", "terminal_ddN_candidate_cocktail_khids.csv"))
```

#### Do GRNBoost2 50 times.

``` r
set.seed(16)

expression_count_matrix_tsv <- paste0("GRN/", "terminal_ddN_GRN_df.tsv")
tf_khid_csv <- paste0("GRN/", "terminal_ddN_candidate_cocktail_khids.csv")
outname <- "GRN/terminal_ddN/terminal_ddN"
min_i <- 0L
max_i <- 50L

# Took ~9 mins to run.
do_GRNBoost2(expression_count_matrix_tsv, tf_khid_csv, outname, min_i, max_i)
```

##### sorted_total_importance_tfs

Each element of the list is a tuple of TF and list with Mean and
Standard Deviation of Total Importance Scores. The TFs are sorted in
nonincreasing order based on the mean.

##### sorted_percent_top_list

Each element of the list is a TF and its Percentage Covered as a Top 3
Regulator for all DEGs.

``` r
khid_tf_tsv_file <- "./helper_files/Ciona_khid_TF.tsv"
terminal_ddN_GRN_path <- "/Users/prakritipaul/Git/Cocktails/FINAL_code_and_results/GRN/terminal_ddN/*.tsv"
n <- 3L

# Returns sorted_total_importance_tfs and sorted_percent_top_list
terminal_ddN_GRNBoost2_outputs <- GRNBoost2_pipeline(terminal_ddN_GRN_path, khid_tf_tsv_file, n)
```

##### Inspect Output contents.

``` r
terminal_ddN_sorted_total_importance_tfs <- terminal_ddN_GRNBoost2_outputs[1]
terminal_ddN_sorted_percent_top_list <- terminal_ddN_GRNBoost2_outputs[2]
```

#### Note on Final Result

##### We are once again able to uncover the Cocktail from studying the terminal state. This further supports the need to study gene regulatory changes occurring at the time of specification.

## (4) Subsampling

> Finally, in order to build evidence for our predictions, we sample 8
> cells 9 times, run the Group Lasso Model, and keep track of its
> inferred Cocktail. We count the number of times that a TF in our
> predicted Cocktail appears in a subsample’s inferred Cocktail and thus
> get a percentage/confidence for the overall prediction.

> In our analyses, we obtain a 100% confidence for each of the TFs Dmbx,
> Lhx1, and HMG1/2 in the ddN Cocktail.

### Subsample 8 cells 9 times.

``` r
subsampling_out_dir <- "/Users/prakritipaul/Git/Cocktails/FINAL_code_and_results/subsampling/"
```

# Do subsampling by leaving one cell out.

``` r
ddN_seed <- 8
# This list contains output of grepreg_routine for all 9 rounds of subsampling.
ddN_subsampling_output_list <- list()

i <- 1
for (ddN_cell in ddN_precursor_pseudostate_cells) {
  cat("Doing for i = ", i, "cell = ", ddN_cell, "\n")
  ddN_cells <- ddN_precursor_pseudostate_cells[!ddN_precursor_pseudostate_cells %in% ddN_cell]
  
  # Sometimes the model doesn't converge.
  tryCatch({ddN_subsampling_output <- grepreg_routine(ddN_cells,
                                  ddN_candidate_cocktail_tfs,
                                  cleaned_ddN_precursor_pseudostate_DEGs,
                                  ddN_TFs_in_precursor_pseudostate_mat,
                                  ddN_precursor_DEGs_in_precursors_mat,
                                  ddN_seed)
  ddN_subsampling_output_list[[i]] <- ddN_subsampling_output
  
  }, error = function(error) {
    cat(paste("An error occurred for ddN cell", ddN_cell, "and i =", i, "\n"))
    ddN_subsampling_output_list[[i]] <- ddN_cell
  })
  i <- i+1
}
```

# Let’s make the subsampled group tf df.

``` r
ddN_combined_subsampled_grouped_tf_df <- make_lambda_to_grouped_tf_df(ddN_subsampling_output_list)
```

# Write it out.

``` r
write_out_group_lasso_tf_df(ddN_combined_subsampled_grouped_tf_df, subsampling_out_dir, "ddN_combined_subsampled_tf_df.csv")
```

# Also write out CT_group_lassoed_tf_list.

``` r
write_out_CT_group_lassoed_tf_list(ddN_group_lassoed_tf_list, subsampling_out_dir, "ddN_group_lassoed_tf_df.csv")
```

### Run Subsampling Prediction Confidence Pipeline

``` r
source_python("subsampling_prediction_confidence.py")
```

``` r
ddN_group_lassoed_tf_df_csv <- paste0(subsampling_out_dir, "ddN_group_lassoed_tf_df.csv")
ddN_subsampled_tf_df_csv <- paste0(subsampling_out_dir, "ddN_combined_subsampled_tf_df.csv")
ddN_tf_set_sizes <- list(3, 4, 5, 6)
ddN_ideal_set_size <- 3

ddN_cocktail_subsampling_confidence <- subsampling_prediction_confidence(ddN_group_lassoed_tf_df_csv, ddN_subsampled_tf_df_csv, ddN_tf_set_sizes, ddN_ideal_set_size)

ddN_cocktail_subsampling_confidence
```

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
