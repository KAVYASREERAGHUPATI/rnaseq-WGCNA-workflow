############################################################
# COMPLETE GENERAL WGCNA PIPELINE FOR TWO-GROUP RNA-seq DATA
#
# Purpose:
#   Perform a complete Weighted Gene Co-expression Network
#   Analysis (WGCNA) using a normalized expression matrix and
#   a sample-trait file.
#
# REQUIRED INPUT FILES — ONLY TWO
#   1. Normalized expression matrix:
#        - Rows: genes
#        - Columns: samples
#        - First column: gene IDs
#
#   2. Trait file:
#        - First column: sample
#        - Second column: condition
#        - Example:
#            sample,condition
#            Sample_01,Control
#            Sample_02,Treatment
#
# IMPORTANT:
#   Use normalized/transformed expression values such as VST,
#   rlog, log2(normalized_count + 1), or another suitable
#   approximately homoscedastic expression matrix.
#   Do not use raw integer read counts directly.
#
# MAIN OUTPUTS
#   - Input/QC tables
#   - Sample clustering
#   - Soft-threshold analysis
#   - Gene dendrogram and module colors
#   - Module eigengene clustering/correlation
#   - Module-trait correlation and P-values
#   - Gene significance and module membership
#   - Intramodular connectivity
#   - All genes from every module
#   - Top 50 hub genes from every module
#   - GS-versus-MM plots for significant modules
#   - Correlation, adjacency, and TOM matrices
#   - TOM heatmap
#   - TOM-based Cytoscape edge/node files
#   - Reusable R objects and run summary
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

# ==========================================================
# 0. USER SETTINGS — EDIT THIS SECTION ONLY
# ==========================================================

# Input 1: normalized expression matrix
expression_file <- "C:/path/to/normalized_expression_matrix.csv"

# Input 2: sample trait file
trait_file <- "C:/path/to/sample_traits.csv"

# Output directory
output_dir <- file.path(path.expand("~"), "Desktop", "General_WGCNA_Results")

# The condition name that should be coded as 0.
# Every other condition will be coded as 1.
reference_condition <- "Control"

# Maximum number of most-variable genes retained for WGCNA.
# Increase only when sufficient RAM is available.
max_genes_by_variance <- 10000

# Network settings
network_type <- "signed"
tom_type <- "signed"
minimum_module_size <- 30
merge_cut_height <- 0.25
minimum_scale_free_R2 <- 0.80

# Significant module definition for GS-versus-MM plots
module_trait_p_cutoff <- 0.05
module_trait_cor_cutoff <- 0.30

# Number of hub genes exported from each module
top_hub_genes_per_module <- 50

# Matrix export safeguard.
# Full correlation, adjacency, and TOM matrices contain genes²
# values and can become extremely large.
matrix_export_gene_limit <- 8000

# Number of highly connected genes displayed in the TOM heatmap
tom_heatmap_gene_limit <- 500

# Cytoscape settings
cytoscape_top_genes_per_module <- 50
cytoscape_edge_threshold <- 0.10

# ==========================================================
# 1. INSTALL AND LOAD REQUIRED PACKAGES
# ==========================================================

cran_packages <- c(
  "WGCNA",
  "dynamicTreeCut",
  "dplyr"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

allowWGCNAThreads()

# ==========================================================
# 2. CREATE OUTPUT FOLDERS
# ==========================================================

folders <- c(
  "01_QC",
  "02_SoftThreshold",
  "03_Network",
  "04_ModuleTrait",
  "05_HubGenes",
  "06_Cytoscape",
  "07_Matrices",
  "08_Tables",
  "09_Objects"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (folder in folders) {
  dir.create(
    file.path(output_dir, folder),
    recursive = TRUE,
    showWarnings = FALSE
  )
}

cat("\nResults will be saved in:\n", output_dir, "\n\n")

# ==========================================================
# 3. CHECK THE TWO REQUIRED INPUT FILES
# ==========================================================

if (!file.exists(expression_file)) {
  stop("Expression file was not found: ", expression_file)
}

if (!file.exists(trait_file)) {
  stop("Trait file was not found: ", trait_file)
}

# ==========================================================
# 4. READ AND PREPARE THE EXPRESSION MATRIX
# ==========================================================

# Expected input:
#   First column = gene ID
#   Remaining columns = samples

expression_raw <- read.csv(
  expression_file,
  row.names = 1,
  check.names = FALSE
)

expression_raw <- as.data.frame(expression_raw)

if (nrow(expression_raw) < 2 || ncol(expression_raw) < 6) {
  stop(
    "The expression matrix appears too small. ",
    "It must contain genes in rows and samples in columns."
  )
}

if (anyDuplicated(rownames(expression_raw)) > 0) {
  stop("Duplicate gene IDs were detected in the expression matrix.")
}

if (anyDuplicated(colnames(expression_raw)) > 0) {
  stop("Duplicate sample names were detected in the expression matrix.")
}

# Convert all expression columns to numeric.
# Non-numeric entries become NA and are handled below.
expression_raw[] <- lapply(
  expression_raw,
  function(x) suppressWarnings(as.numeric(as.character(x)))
)

# Remove genes containing missing/non-finite values.
finite_gene <- apply(
  expression_raw,
  1,
  function(x) all(is.finite(x))
)

expression_clean <- expression_raw[finite_gene, , drop = FALSE]

# Remove genes with zero variance because they cannot contribute
# to a correlation-based co-expression network.
gene_variance <- apply(expression_clean, 1, var, na.rm = TRUE)
expression_clean <- expression_clean[
  is.finite(gene_variance) & gene_variance > 0,
  ,
  drop = FALSE
]

if (nrow(expression_clean) < minimum_module_size * 2) {
  stop("Too few usable genes remain after expression filtering.")
}

# Retain the most-variable genes if the dataset is very large.
gene_variance <- apply(expression_clean, 1, var, na.rm = TRUE)
number_to_keep <- min(max_genes_by_variance, nrow(expression_clean))

selected_genes <- names(
  sort(gene_variance, decreasing = TRUE)
)[seq_len(number_to_keep)]

expression_filtered <- expression_clean[
  selected_genes,
  ,
  drop = FALSE
]

# WGCNA requires samples in rows and genes in columns.
datExpr <- as.data.frame(t(expression_filtered))

write.csv(
  data.frame(
    GeneID = rownames(expression_filtered),
    expression_filtered,
    check.names = FALSE
  ),
  file.path(
    output_dir,
    "01_QC",
    "Filtered_Normalized_Expression_Matrix_Used.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(GeneID = colnames(datExpr)),
  file.path(
    output_dir,
    "01_QC",
    "Genes_Retained_For_WGCNA.csv"
  ),
  row.names = FALSE
)

# ==========================================================
# 5. READ, VALIDATE, AND MATCH THE TRAIT FILE
# ==========================================================

# Expected input:
#   First column  = sample
#   Second column = condition

traits <- read.csv(
  trait_file,
  check.names = FALSE
)

if (ncol(traits) < 2) {
  stop(
    "The trait file must contain at least two columns: ",
    "sample and condition."
  )
}

colnames(traits)[1:2] <- c("sample", "condition")

traits$sample <- trimws(as.character(traits$sample))
traits$condition <- trimws(as.character(traits$condition))

if (anyDuplicated(traits$sample) > 0) {
  stop("Duplicate sample IDs were detected in the trait file.")
}

common_samples <- intersect(
  rownames(datExpr),
  traits$sample
)

if (length(common_samples) < 6) {
  stop(
    "Fewer than six sample names match between the expression ",
    "matrix and trait file. Check spelling and capitalization."
  )
}

# Preserve expression-matrix sample order.
datExpr <- datExpr[common_samples, , drop = FALSE]

traits <- traits[
  match(common_samples, traits$sample),
  ,
  drop = FALSE
]

rownames(traits) <- traits$sample

if (!(tolower(reference_condition) %in%
      tolower(unique(traits$condition)))) {
  stop(
    "The reference condition '", reference_condition,
    "' was not found in the trait file."
  )
}

# Binary trait coding:
#   Reference condition = 0
#   Every other condition = 1
traits$Trait_numeric <- ifelse(
  tolower(traits$condition) == tolower(reference_condition),
  0,
  1
)

if (length(unique(traits$Trait_numeric)) != 2) {
  stop(
    "The trait file must contain the reference condition and ",
    "at least one non-reference condition."
  )
}

write.csv(
  traits,
  file.path(
    output_dir,
    "01_QC",
    "Matched_Samples_and_Traits.csv"
  ),
  row.names = FALSE
)

writeLines(
  c(
    paste0("Reference condition coded as 0: ", reference_condition),
    "All non-reference conditions are coded as 1.",
    paste0(
      "Positive module-trait correlation: module eigengene is ",
      "higher in the non-reference group."
    ),
    paste0(
      "Negative module-trait correlation: module eigengene is ",
      "higher in the reference group."
    )
  ),
  file.path(
    output_dir,
    "01_QC",
    "Trait_Encoding_Key.txt"
  )
)

# ==========================================================
# 6. QUALITY CONTROL OF SAMPLES AND GENES
# ==========================================================

gsg <- goodSamplesGenes(
  datExpr,
  verbose = 3
)

if (!gsg$allOK) {
  removed_samples <- rownames(datExpr)[!gsg$goodSamples]
  removed_genes <- colnames(datExpr)[!gsg$goodGenes]

  datExpr <- datExpr[
    gsg$goodSamples,
    gsg$goodGenes,
    drop = FALSE
  ]

  traits <- traits[
    rownames(datExpr),
    ,
    drop = FALSE
  ]

  write.csv(
    data.frame(Sample = removed_samples),
    file.path(output_dir, "01_QC", "Removed_Samples.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(GeneID = removed_genes),
    file.path(output_dir, "01_QC", "Removed_Genes.csv"),
    row.names = FALSE
  )
}

# Sample hierarchical clustering is used to inspect possible
# outlier samples before network construction.
sample_tree <- hclust(
  dist(datExpr),
  method = "average"
)

png(
  file.path(output_dir, "01_QC", "Sample_Clustering.png"),
  width = 2200,
  height = 1400,
  res = 200
)

par(mar = c(5, 4, 4, 2))
plot(
  sample_tree,
  main = "Sample clustering to detect possible outliers",
  sub = "",
  xlab = "",
  cex.lab = 1.2,
  cex.axis = 0.8
)

dev.off()

# ==========================================================
# 7. SELECT THE SOFT-THRESHOLDING POWER
# ==========================================================

candidate_powers <- c(1:20, seq(22, 30, by = 2))

sft <- pickSoftThreshold(
  datExpr,
  powerVector = candidate_powers,
  networkType = network_type,
  verbose = 5
)

fit_indices <- sft$fitIndices

signed_R2 <- -sign(fit_indices[, 3]) * fit_indices[, 2]

eligible_power_index <- which(
  is.finite(signed_R2) &
  signed_R2 >= minimum_scale_free_R2
)

if (length(eligible_power_index) > 0) {
  soft_power <- fit_indices$Power[
    eligible_power_index[1]
  ]
} else {
  # Fallback: select the power with the highest valid signed R².
  best_index <- which.max(
    ifelse(is.finite(signed_R2), signed_R2, -Inf)
  )
  soft_power <- fit_indices$Power[best_index]

  warning(
    "No tested power reached the requested scale-free R². ",
    "The power with the highest observed signed R² was selected."
  )
}

write.csv(
  fit_indices,
  file.path(
    output_dir,
    "02_SoftThreshold",
    "Soft_Threshold_Fit_Indices.csv"
  ),
  row.names = FALSE
)

writeLines(
  paste0("Selected soft-thresholding power: ", soft_power),
  file.path(
    output_dir,
    "02_SoftThreshold",
    "Selected_Soft_Power.txt"
  )
)

png(
  file.path(
    output_dir,
    "02_SoftThreshold",
    "Soft_Threshold_Plots.png"
  ),
  width = 2600,
  height = 1300,
  res = 200
)

par(mfrow = c(1, 2))

plot(
  fit_indices[, 1],
  signed_R2,
  xlab = "Soft-thresholding power",
  ylab = "Scale-free topology fit, signed R²",
  type = "n",
  main = "Scale independence"
)

text(
  fit_indices[, 1],
  signed_R2,
  labels = fit_indices[, 1],
  cex = 0.8
)

abline(
  h = minimum_scale_free_R2,
  lty = 2
)

plot(
  fit_indices[, 1],
  fit_indices[, 5],
  xlab = "Soft-thresholding power",
  ylab = "Mean connectivity",
  type = "n",
  main = "Mean connectivity"
)

text(
  fit_indices[, 1],
  fit_indices[, 5],
  labels = fit_indices[, 1],
  cex = 0.8
)

dev.off()

# ==========================================================
# 8. CONSTRUCT THE NETWORK AND IDENTIFY MODULES
# ==========================================================

net <- blockwiseModules(
  datExpr,
  power = soft_power,
  networkType = network_type,
  TOMType = tom_type,
  minModuleSize = minimum_module_size,
  mergeCutHeight = merge_cut_height,
  reassignThreshold = 0,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  maxBlockSize = ncol(datExpr),
  saveTOMs = TRUE,
  saveTOMFileBase = file.path(
    output_dir,
    "09_Objects",
    "Blockwise_TOM"
  ),
  verbose = 3
)

module_colors <- labels2colors(net$colors)

# Recalculate module eigengenes using color names.
MEs_initial <- moduleEigengenes(
  datExpr,
  colors = module_colors
)$eigengenes

MEs <- orderMEs(MEs_initial)

gene_module_table <- data.frame(
  GeneID = colnames(datExpr),
  Module = module_colors,
  stringsAsFactors = FALSE
)

write.csv(
  gene_module_table,
  file.path(
    output_dir,
    "08_Tables",
    "Gene_Module_Assignment.csv"
  ),
  row.names = FALSE
)

module_sizes <- as.data.frame(
  table(module_colors),
  stringsAsFactors = FALSE
)

colnames(module_sizes) <- c(
  "Module",
  "Number_of_Genes"
)

module_sizes <- module_sizes[
  order(-module_sizes$Number_of_Genes),
  ,
  drop = FALSE
]

write.csv(
  module_sizes,
  file.path(
    output_dir,
    "08_Tables",
    "Module_Sizes.csv"
  ),
  row.names = FALSE
)

# Gene dendrogram with module colors.
png(
  file.path(
    output_dir,
    "03_Network",
    "Gene_Dendrogram_and_Module_Colors.png"
  ),
  width = 2800,
  height = 1800,
  res = 200
)

plotDendroAndColors(
  net$dendrograms[[1]],
  module_colors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and detected co-expression modules"
)

dev.off()

# ==========================================================
# 9. MODULE EIGENGENE ANALYSIS
# ==========================================================

ME_correlation <- cor(
  MEs,
  use = "pairwise.complete.obs"
)

write.csv(
  data.frame(
    ModuleEigengene = rownames(ME_correlation),
    ME_correlation,
    check.names = FALSE
  ),
  file.path(
    output_dir,
    "03_Network",
    "Module_Eigengene_Correlation.csv"
  ),
  row.names = FALSE
)

# Module eigengene clustering.
ME_tree <- hclust(
  as.dist(1 - ME_correlation),
  method = "average"
)

png(
  file.path(
    output_dir,
    "03_Network",
    "Module_Eigengene_Clustering.png"
  ),
  width = 2000,
  height = 1400,
  res = 200
)

plot(
  ME_tree,
  main = "Clustering of module eigengenes",
  xlab = "",
  sub = ""
)

abline(
  h = merge_cut_height,
  lty = 2
)

dev.off()

# Module eigengene correlation heatmap.
png(
  file.path(
    output_dir,
    "03_Network",
    "Module_Eigengene_Correlation_Heatmap.png"
  ),
  width = 1900,
  height = 1700,
  res = 200
)

labeledHeatmap(
  Matrix = ME_correlation,
  xLabels = colnames(ME_correlation),
  yLabels = rownames(ME_correlation),
  ySymbols = rownames(ME_correlation),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = signif(ME_correlation, 2),
  setStdMargins = FALSE,
  cex.text = 0.7,
  zlim = c(-1, 1),
  main = "Module eigengene correlation"
)

dev.off()

# ==========================================================
# 10. MODULE–TRAIT RELATIONSHIPS
# ==========================================================

trait_numeric <- as.data.frame(
  traits$Trait_numeric
)

colnames(trait_numeric) <- "Condition"

module_trait_cor <- cor(
  MEs,
  trait_numeric,
  use = "pairwise.complete.obs"
)

module_trait_p <- corPvalueStudent(
  module_trait_cor,
  nSamples = nrow(datExpr)
)

write.csv(
  data.frame(
    ModuleEigengene = rownames(module_trait_cor),
    module_trait_cor,
    check.names = FALSE
  ),
  file.path(
    output_dir,
    "04_ModuleTrait",
    "Module_Trait_Correlations.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    ModuleEigengene = rownames(module_trait_p),
    module_trait_p,
    check.names = FALSE
  ),
  file.path(
    output_dir,
    "04_ModuleTrait",
    "Module_Trait_Pvalues.csv"
  ),
  row.names = FALSE
)

module_trait_summary <- data.frame(
  ModuleEigengene = rownames(module_trait_cor),
  Module = sub("^ME", "", rownames(module_trait_cor)),
  Trait = "Condition",
  Correlation = as.numeric(module_trait_cor[, 1]),
  P_value = as.numeric(module_trait_p[, 1]),
  stringsAsFactors = FALSE
)

module_trait_summary$Significant <- (
  module_trait_summary$P_value < module_trait_p_cutoff &
  abs(module_trait_summary$Correlation) >=
    module_trait_cor_cutoff
)

module_trait_summary$Interpretation <- ifelse(
  module_trait_summary$Correlation > 0,
  paste0(
    "Positive association: module eigengene is higher in ",
    "the non-reference condition."
  ),
  paste0(
    "Negative association: module eigengene is higher in ",
    "the reference condition."
  )
)

module_trait_summary <- module_trait_summary[
  order(module_trait_summary$P_value),
  ,
  drop = FALSE
]

write.csv(
  module_trait_summary,
  file.path(
    output_dir,
    "08_Tables",
    "Module_Trait_Summary.csv"
  ),
  row.names = FALSE
)

heatmap_text <- paste0(
  signif(module_trait_cor, 2),
  "\n(",
  signif(module_trait_p, 2),
  ")"
)

dim(heatmap_text) <- dim(module_trait_cor)

png(
  file.path(
    output_dir,
    "04_ModuleTrait",
    "Module_Trait_Relationship_Heatmap.png"
  ),
  width = 1700,
  height = 2200,
  res = 200
)

par(mar = c(6, 9, 3, 3))

labeledHeatmap(
  Matrix = module_trait_cor,
  xLabels = "Condition",
  yLabels = rownames(module_trait_cor),
  ySymbols = rownames(module_trait_cor),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = heatmap_text,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1, 1),
  main = "Module–trait relationships"
)

dev.off()

# ==========================================================
# 11. GENE SIGNIFICANCE AND MODULE MEMBERSHIP
# ==========================================================

# Gene significance:
# correlation between each gene and the numeric trait.
GS <- as.numeric(
  cor(
    datExpr,
    traits$Trait_numeric,
    use = "pairwise.complete.obs"
  )
)

GS_pvalue <- as.numeric(
  corPvalueStudent(
    GS,
    nSamples = nrow(datExpr)
  )
)

# Module membership:
# correlation between each gene and each module eigengene.
MM <- as.data.frame(
  cor(
    datExpr,
    MEs,
    use = "pairwise.complete.obs"
  )
)

MM_pvalue <- as.data.frame(
  corPvalueStudent(
    as.matrix(MM),
    nSamples = nrow(datExpr)
  )
)

colnames(MM) <- paste0(
  "MM_",
  sub("^ME", "", colnames(MM))
)

colnames(MM_pvalue) <- paste0(
  "MM_pvalue_",
  sub("^ME", "", colnames(MM_pvalue))
)

# Intramodular connectivity:
# kWithin = connectivity to genes in the same module.
# kOut    = connectivity to genes outside the module.
# kTotal  = total network connectivity.
# kDiff   = kWithin - kOut.
connectivity <- intramodularConnectivity.fromExpr(
  datExpr,
  colors = module_colors,
  power = soft_power,
  networkType = network_type
)

all_gene_statistics <- data.frame(
  GeneID = colnames(datExpr),
  Module = module_colors,
  Gene_Significance = GS,
  GS_Pvalue = GS_pvalue,
  kTotal = connectivity$kTotal,
  kWithin = connectivity$kWithin,
  kOut = connectivity$kOut,
  kDiff = connectivity$kDiff,
  stringsAsFactors = FALSE
)

all_gene_statistics <- cbind(
  all_gene_statistics,
  MM,
  MM_pvalue
)

# Extract each gene's module-membership value and P-value
# for its own assigned module.
all_gene_statistics$MM_Own_Module <- NA_real_
all_gene_statistics$MM_Own_Module_Pvalue <- NA_real_

for (i in seq_len(nrow(all_gene_statistics))) {
  module_name <- all_gene_statistics$Module[i]
  mm_column <- paste0("MM_", module_name)
  p_column <- paste0("MM_pvalue_", module_name)

  if (mm_column %in% colnames(all_gene_statistics)) {
    all_gene_statistics$MM_Own_Module[i] <-
      all_gene_statistics[i, mm_column]
  }

  if (p_column %in% colnames(all_gene_statistics)) {
    all_gene_statistics$MM_Own_Module_Pvalue[i] <-
      all_gene_statistics[i, p_column]
  }
}

write.csv(
  all_gene_statistics,
  file.path(
    output_dir,
    "08_Tables",
    "All_Genes_ModuleMembership_GeneSignificance_Connectivity.csv"
  ),
  row.names = FALSE
)

# ==========================================================
# 12. EXPORT ALL GENES AND TOP 50 HUB GENES PER MODULE
# ==========================================================

module_names <- sort(unique(module_colors))
top_hub_combined <- list()

for (module_name in module_names) {

  module_data <- all_gene_statistics[
    all_gene_statistics$Module == module_name,
    ,
    drop = FALSE
  ]

  # Rank genes using:
  #   1. absolute own-module membership
  #   2. intramodular connectivity
  #   3. absolute gene significance
  module_data <- module_data[
    order(
      -abs(module_data$MM_Own_Module),
      -module_data$kWithin,
      -abs(module_data$Gene_Significance)
    ),
    ,
    drop = FALSE
  ]

  write.csv(
    module_data,
    file.path(
      output_dir,
      "05_HubGenes",
      paste0(
        "All_Genes_Module_",
        module_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )

  # Grey genes are exported for completeness but are not
  # interpreted as a biologically coherent module.
  if (module_name != "grey") {
    number_of_hubs <- min(
      top_hub_genes_per_module,
      nrow(module_data)
    )

    top_hubs <- module_data[
      seq_len(number_of_hubs),
      ,
      drop = FALSE
    ]

    top_hubs$Hub_Rank <- seq_len(nrow(top_hubs))

    write.csv(
      top_hubs,
      file.path(
        output_dir,
        "05_HubGenes",
        paste0(
          "Top",
          top_hub_genes_per_module,
          "_Hub_Genes_Module_",
          module_name,
          ".csv"
        )
      ),
      row.names = FALSE
    )

    top_hub_combined[[module_name]] <- top_hubs
  }
}

if (length(top_hub_combined) > 0) {
  combined_hubs <- dplyr::bind_rows(
    top_hub_combined
  )

  write.csv(
    combined_hubs,
    file.path(
      output_dir,
      "05_HubGenes",
      paste0(
        "Top",
        top_hub_genes_per_module,
        "_Hub_Genes_All_Modules.csv"
      )
    ),
    row.names = FALSE
  )
}

# ==========================================================
# 13. GS-VERSUS-MM PLOTS FOR SIGNIFICANT MODULES
# ==========================================================

significant_modules <- module_trait_summary$Module[
  module_trait_summary$Significant &
  module_trait_summary$Module != "grey"
]

for (module_name in unique(significant_modules)) {

  genes_in_module <- module_colors == module_name
  mm_column <- paste0("MM_", module_name)

  if (
    sum(genes_in_module) >= 3 &&
    mm_column %in% colnames(all_gene_statistics)
  ) {

    png(
      file.path(
        output_dir,
        "04_ModuleTrait",
        paste0(
          "GS_vs_MM_Module_",
          module_name,
          ".png"
        )
      ),
      width = 1700,
      height = 1500,
      res = 200
    )

    verboseScatterplot(
      abs(
        all_gene_statistics[
          genes_in_module,
          mm_column
        ]
      ),
      abs(
        all_gene_statistics$Gene_Significance[
          genes_in_module
        ]
      ),
      xlab = paste0(
        "Absolute module membership in ",
        module_name,
        " module"
      ),
      ylab = "Absolute gene significance",
      main = paste0(
        "Gene significance versus module membership: ",
        module_name
      ),
      col = module_name
    )

    dev.off()
  }
}

# ==========================================================
# 14. CORRELATION, ADJACENCY, TOM, AND TOM HEATMAP
# ==========================================================

if (ncol(datExpr) <= matrix_export_gene_limit) {

  gene_correlation <- cor(
    datExpr,
    use = "pairwise.complete.obs"
  )

  adjacency_matrix <- adjacency(
    datExpr,
    power = soft_power,
    type = network_type
  )

  TOM_matrix <- TOMsimilarity(
    adjacency_matrix,
    TOMType = tom_type
  )

  rownames(gene_correlation) <- colnames(datExpr)
  colnames(gene_correlation) <- colnames(datExpr)

  rownames(adjacency_matrix) <- colnames(datExpr)
  colnames(adjacency_matrix) <- colnames(datExpr)

  rownames(TOM_matrix) <- colnames(datExpr)
  colnames(TOM_matrix) <- colnames(datExpr)

  saveRDS(
    gene_correlation,
    file.path(
      output_dir,
      "09_Objects",
      "Gene_Correlation_Matrix.rds"
    )
  )

  saveRDS(
    adjacency_matrix,
    file.path(
      output_dir,
      "09_Objects",
      "Adjacency_Matrix.rds"
    )
  )

  saveRDS(
    TOM_matrix,
    file.path(
      output_dir,
      "09_Objects",
      "TOM_Matrix.rds"
    )
  )

  # Compressed CSV files reduce disk space.
  write.csv(
    gene_correlation,
    gzfile(
      file.path(
        output_dir,
        "07_Matrices",
        "Gene_Correlation_Matrix.csv.gz"
      )
    )
  )

  write.csv(
    adjacency_matrix,
    gzfile(
      file.path(
        output_dir,
        "07_Matrices",
        "Adjacency_Matrix.csv.gz"
      )
    )
  )

  write.csv(
    TOM_matrix,
    gzfile(
      file.path(
        output_dir,
        "07_Matrices",
        "TOM_Matrix.csv.gz"
      )
    )
  )

  # Display only the most highly connected genes in the TOM heatmap.
  mean_connectivity <- rowMeans(
    TOM_matrix,
    na.rm = TRUE
  )

  heatmap_gene_count <- min(
    tom_heatmap_gene_limit,
    length(mean_connectivity)
  )

  heatmap_genes <- names(
    sort(
      mean_connectivity,
      decreasing = TRUE
    )
  )[seq_len(heatmap_gene_count)]

  TOM_subset <- TOM_matrix[
    heatmap_genes,
    heatmap_genes,
    drop = FALSE
  ]

  TOM_dissimilarity <- 1 - TOM_subset
  TOM_tree <- hclust(
    as.dist(TOM_dissimilarity),
    method = "average"
  )

  plot_order <- TOM_tree$order

  png(
    file.path(
      output_dir,
      "03_Network",
      "TOM_Heatmap_TopConnectedGenes.png"
    ),
    width = 2200,
    height = 2000,
    res = 200
  )

  TOMplot(
    TOM_dissimilarity[
      plot_order,
      plot_order
    ],
    TOM_tree,
    module_colors[
      match(
        heatmap_genes[plot_order],
        colnames(datExpr)
      )
    ],
    main = "TOM heatmap of highly connected genes"
  )

  dev.off()

} else {

  warning(
    "Full matrices were not exported because the analysis ",
    "contains more than ", matrix_export_gene_limit,
    " genes. Blockwise TOM files and network objects are ",
    "still saved."
  )

  writeLines(
    c(
      paste0(
        "Full correlation, adjacency, and TOM matrices were ",
        "not exported because the dataset contained ",
        ncol(datExpr),
        " genes."
      ),
      paste0(
        "Current matrix export limit: ",
        matrix_export_gene_limit
      ),
      paste0(
        "Increase matrix_export_gene_limit only when sufficient ",
        "RAM and disk space are available."
      )
    ),
    file.path(
      output_dir,
      "07_Matrices",
      "Matrix_Export_Skipped.txt"
    )
  )

  TOM_matrix <- NULL
}

# ==========================================================
# 15. TOM-BASED CYTOSCAPE NETWORKS
# ==========================================================

modules_for_cytoscape <- unique(
  module_trait_summary$Module[
    module_trait_summary$Significant &
    module_trait_summary$Module != "grey"
  ]
)

# If no module reaches the selected significance thresholds,
# export the strongest non-grey module so the workflow still
# produces a Cytoscape example.
if (length(modules_for_cytoscape) == 0) {
  non_grey_summary <- module_trait_summary[
    module_trait_summary$Module != "grey",
    ,
    drop = FALSE
  ]

  if (nrow(non_grey_summary) > 0) {
    modules_for_cytoscape <- non_grey_summary$Module[1]
  }
}

for (module_name in modules_for_cytoscape) {

  module_table <- all_gene_statistics[
    all_gene_statistics$Module == module_name,
    ,
    drop = FALSE
  ]

  module_table <- module_table[
    order(
      -abs(module_table$MM_Own_Module),
      -module_table$kWithin,
      -abs(module_table$Gene_Significance)
    ),
    ,
    drop = FALSE
  ]

  selected_gene_count <- min(
    cytoscape_top_genes_per_module,
    nrow(module_table)
  )

  selected_genes <- module_table$GeneID[
    seq_len(selected_gene_count)
  ]

  selected_genes <- intersect(
    selected_genes,
    colnames(datExpr)
  )

  if (length(selected_genes) < 3) {
    next
  }

  # Calculate TOM similarity specifically for the selected
  # genes in this module.
  module_TOM <- TOMsimilarityFromExpr(
    datExpr[, selected_genes, drop = FALSE],
    power = soft_power,
    networkType = network_type,
    TOMType = tom_type
  )

  rownames(module_TOM) <- selected_genes
  colnames(module_TOM) <- selected_genes

  # Remove self-connections.
  diag(module_TOM) <- 0

  edge_index <- which(
    module_TOM >= cytoscape_edge_threshold,
    arr.ind = TRUE
  )

  edge_index <- edge_index[
    edge_index[, 1] < edge_index[, 2],
    ,
    drop = FALSE
  ]

  if (nrow(edge_index) == 0) {
    warning(
      "No Cytoscape edges passed the threshold for module ",
      module_name,
      ". Consider lowering cytoscape_edge_threshold."
    )
    next
  }

  cytoscape_edges <- data.frame(
    fromNode = rownames(module_TOM)[edge_index[, 1]],
    toNode = colnames(module_TOM)[edge_index[, 2]],
    TOM_similarity = module_TOM[edge_index],
    stringsAsFactors = FALSE
  )

  cytoscape_edges <- cytoscape_edges[
    order(
      -cytoscape_edges$TOM_similarity
    ),
    ,
    drop = FALSE
  ]

  connected_genes <- unique(
    c(
      cytoscape_edges$fromNode,
      cytoscape_edges$toNode
    )
  )

  cytoscape_nodes <- module_table[
    module_table$GeneID %in% connected_genes,
    c(
      "GeneID",
      "Module",
      "MM_Own_Module",
      "MM_Own_Module_Pvalue",
      "Gene_Significance",
      "GS_Pvalue",
      "kWithin",
      "kTotal",
      "kOut",
      "kDiff"
    ),
    drop = FALSE
  ]

  write.csv(
    cytoscape_edges,
    file.path(
      output_dir,
      "06_Cytoscape",
      paste0(
        "Cytoscape_Edges_Module_",
        module_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )

  write.csv(
    cytoscape_nodes,
    file.path(
      output_dir,
      "06_Cytoscape",
      paste0(
        "Cytoscape_Nodes_Module_",
        module_name,
        ".csv"
      )
    ),
    row.names = FALSE
  )
}

# ==========================================================
# 16. SAVE REUSABLE R OBJECTS
# ==========================================================

saveRDS(
  datExpr,
  file.path(
    output_dir,
    "09_Objects",
    "Expression_Data_Used.rds"
  )
)

saveRDS(
  traits,
  file.path(
    output_dir,
    "09_Objects",
    "Traits_Used.rds"
  )
)

saveRDS(
  net,
  file.path(
    output_dir,
    "09_Objects",
    "WGCNA_Network_Object.rds"
  )
)

saveRDS(
  MEs,
  file.path(
    output_dir,
    "09_Objects",
    "Module_Eigengenes.rds"
  )
)

saveRDS(
  all_gene_statistics,
  file.path(
    output_dir,
    "09_Objects",
    "All_Gene_Statistics.rds"
  )
)

# ==========================================================
# 17. AUTOMATIC ANALYSIS SUMMARY
# ==========================================================

strongest_module <- module_trait_summary[
  module_trait_summary$Module != "grey",
  ,
  drop = FALSE
]

if (nrow(strongest_module) > 0) {
  strongest_module <- strongest_module[1, , drop = FALSE]
}

summary_lines <- c(
  "GENERAL WGCNA ANALYSIS SUMMARY",
  "========================================",
  "",
  paste0("Expression input: ", expression_file),
  paste0("Trait input: ", trait_file),
  paste0("Reference condition: ", reference_condition),
  paste0("Samples analyzed: ", nrow(datExpr)),
  paste0("Genes analyzed: ", ncol(datExpr)),
  paste0("Selected soft-thresholding power: ", soft_power),
  paste0(
    "Non-grey modules detected: ",
    length(setdiff(unique(module_colors), "grey"))
  ),
  ""
)

if (nrow(strongest_module) > 0) {
  summary_lines <- c(
    summary_lines,
    paste0(
      "Strongest non-grey module-trait association: ",
      strongest_module$Module
    ),
    paste0(
      "Correlation: ",
      round(strongest_module$Correlation, 4)
    ),
    paste0(
      "P-value: ",
      signif(strongest_module$P_value, 4)
    ),
    strongest_module$Interpretation,
    ""
  )
}

summary_lines <- c(
  summary_lines,
  paste0(
    "Hub genes were ranked using absolute module membership, ",
    "intramodular connectivity, and absolute gene significance."
  ),
  paste0(
    "Cytoscape networks were generated using TOM similarity."
  ),
  paste0(
    "Grey-module genes were exported but should not normally ",
    "be interpreted as a coherent biological module."
  )
)

writeLines(
  summary_lines,
  file.path(
    output_dir,
    "Analysis_Summary.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    output_dir,
    "Session_Info.txt"
  )
)

cat("\n====================================================\n")
cat("WGCNA ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("Results folder:\n", output_dir, "\n")
cat("====================================================\n")
