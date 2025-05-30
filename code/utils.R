

#' Wrapper function to load data from STOmicsDB as a `SpatialExperiment` object
#'
#' This function downloads an H5AD file, and reads it as a `SpatalExperiment`
#' object
#'
#' @param path Character with file path.
#' 
#' @return A `SpatialExperiment` object.
#' 
#' @importFrom zellkonverter readH5AD
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom scater computeLibraryFactors logNormCounts
#' 
stomics2spe <- function(path, remote = FALSE) {
    
    outfile <- path
    if(remote) {
        # Download file
        outfile <- tempfile(fileext = ".h5ad")
        download.file(path, destfile = outfile)
    }
    
    # Read H5AD file as a SingleCellExperiment object
    sce <- zellkonverter::readH5AD(outfile)
    
    # Convert SingleCellExperiment to SpatialExperiment
    ## Get spatial coords
    coords <- as.matrix(SingleCellExperiment::reducedDim(sce, "spatial"))
    colnames(coords) <- c("x_coord", "y_coord")

    ## Check if assay contains counts or logcounts
    samples <- as.numeric(SummarizedExperiment::assay(sce)[1:3, 1:3])
    check_counts <- all.equal(samples, as.integer(samples))
    
    ## Create SpatialExperiment object
    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = assay(sce)),
        colData = SummarizedExperiment::colData(sce),
        rowData = SummarizedExperiment::rowData(sce),
        spatialCoords = coords
    )
    
    ## If assay contains logcounts, add another assay named `logcounts`
    if(check_counts != TRUE) {
        SummarizedExperiment::assay(spe, "logcounts") <- SummarizedExperiment::assay(sce)
    } else {
        spe <- scater::computeLibraryFactors(spe)
        spe <- scater::logNormCounts(spe)
    }
    
    return(spe)
}


#' Wrapper function to filter a `SpatialExperiment` object in a standard way
#'
#' @param spe A `SpatialExperiment` object
#' @param min_exp Numeric vector of length 2 with the minimum expression and
#' the minimum percentage of samples with the minimum expression, respectively.
#' Default: \code{c(3, 0.5)}, indicating that genes must have a minimum 
#' expression of 5 in at least 0.5% of the spots. 
#' @param domain_col Character indicating the name of the column in
#' the colData slot with information on spatial domains or cell types.
#' Default: "cell_type".
#' @param remove_missing Logical indicating whether to remove samples
#' that have missing values in \strong{domain_col}. Default: TRUE.
#'
#' @return A `SpatialExperiment` object.
#'
#' @details
#' This function filters a `SpatialExperiment` object to:
#' (i). remove genes with zero counts across spots; (ii) remove spots
#' with zero counts across all genes; and (iii) remove lowly expressed
#' genes.
#'
#' @importFrom nnSVG filter_genes
#' @importFrom SummarizedExperiment assay colData
#'
#'
process_spe <- function(
    spe, min_exp = c(3, 0.5), domain_col = "cell_type", remove_missing = TRUE
) {
    
    # Remove genes with zero counts across all spots
    zeroc_genes <- rowSums(SummarizedExperiment::assay(spe)) == 0
    if(sum(zeroc_genes) > 0) { spe <- spe[!zeroc_genes, ] }
    
    # Remove spots with zero counts across all genes
    zeroc_spots <- colSums(SummarizedExperiment::assay(spe)) == 0
    if(sum(zeroc_spots) > 0) { spe <- spe[, !zeroc_spots] }
    
    # Keep only genes with x counts in at least y% of the samples
    spe <- nnSVG::filter_genes(spe, 
        filter_genes_ncounts = min_exp[1],
        filter_genes_pcspots = min_exp[2],
        filter_mito = FALSE
    )
    
    # Remove samples with missing values in domain/cell type annotation
    if(remove_missing) {
        if(domain_col %in% colnames(SummarizedExperiment::colData(spe))) {
            spe <- spe <- spe[, !is.na(SummarizedExperiment::colData(spe)[[domain_col]])]
        } else {
            message("Column '", domain_col, "' was not found in colData.")
        }
    }
    
    return(spe)
}


compare <- function(data, form, ref = NULL) {
    
    # Wilcoxon test - greater and less alternatives
    wilcoxtest_greater <- tibble::as_tibble(data) %>%
        rstatix::wilcox_test(
            formula(form), p.adjust.method = "BH", ref.group = ref,
            alternative = "greater"
        )
    pg <- ifelse("p.adj" %in% names(wilcoxtest_greater), "p.adj", "p")
    wilcoxtest_greater <- wilcoxtest_greater %>% dplyr::select(
        group1, group2, n1, n2, padj_greater = all_of(pg)
    )
    
    wilcoxtest_less <- tibble::as_tibble(data) %>%
        rstatix::wilcox_test(
            formula(form), p.adjust.method = "BH", ref.group = ref,
            alternative = "less"
        )
    pl <- ifelse("p.adj" %in% names(wilcoxtest_less), "p.adj", "p")
    wilcoxtest_less <- wilcoxtest_less %>% dplyr::select(
        group1, group2, n1, n2, padj_less = all_of(pl)
    )
    
    wilcox_summary <- dplyr::inner_join(wilcoxtest_greater, wilcoxtest_less) %>%
        dplyr::mutate(padj_interpretation = dplyr::case_when(
            padj_less < 0.05 ~ "less",
            padj_greater < 0.05 ~ "greater",
            TRUE ~ "ns"
        ))
    
    # Effect sizes for Wilcoxon test - greater and less alternatives
    
    effsize <- tibble::as_tibble(data) %>%
        rstatix::wilcox_effsize(
            formula(form), ref.group = ref,
        ) %>%
        dplyr::select(
            group1, group2, effsize, magnitude
        )
    
    
    result <- as.data.frame(inner_join(wilcox_summary, effsize))
    
    return(result)
}


filter_comparison <- function(compare_output, effsize_cutoff = 0) {
    
    filtered_df <- compare_output |>
        dplyr::filter(padj_interpretation != "ns") |>
        dplyr::filter(effsize >= effsize_cutoff) |>
        mutate(padj = case_when(
            padj_interpretation == "greater" ~ padj_greater,
            padj_interpretation == "less" ~ padj_less
        )) |>
        dplyr::select(
            group1, group2, n1, n2, `p.adj` = padj, 
            interpretation = padj_interpretation,
            effsize, magnitude
        )
    
    return(filtered_df)
}


#' Wrapper function to create violin plot of correlation distros by mode
#' 
#' @param cor_df A data frame with correlations for each gene pair.
#' @param wilcox_df A data frame with Wilcoxon test statistics (P-value
#' and effect sizes) as returned by \code{compare} and \code{filter_comparison}.
#' All comparisons in this data frame will be plotted as significance bars,
#' so any kind of filtering should be done beforehand.
#' 
plot_cor_violin <- function(
        cor_df, wilcox_df, x = "type", y = "cor", ...
    ) {
    
    # Get a table of contrasts and P-values to show
    pval_df <- wilcox_df |>
        rstatix::add_significance() |>
        rstatix::add_xy_position(x = "type")
    
    # Get violin plot
    bg <- grid::linearGradient(colorRampPalette(c("gray90", "white"))(100))
    p_violin <- cor_df |>
        filter(!is.na(.data[[y]])) |>
        ggplot(aes(y = .data[[y]], x = .data[[x]])) +
        geom_violin(aes(fill = .data[[x]])) +
        ggsci::scale_fill_jama() +
        geom_boxplot(width = 0.3) +
        theme_bw() +
        theme(
            panel.grid = element_blank(),
            panel.background = element_rect(fill = bg)
        ) +
        ggpubr::stat_pvalue_manual(
            pval_df, label = "p.adj.signif",
            hide_ns = TRUE,
            ...
        )
    
    return(p_violin)
}


#' Wrapper to compare distributions using Kruskal-Wallis + Dunn
#'
#' @param df A data frame with values and variables.
#' @param var Character, name of column indicating variable levels.
#' @param value Numeric, values to be used in comparisons.
#'
#' @return A data frame with CLD for each level.
#' @noRd
#' 
cld_kw_dunn <- function(df, var = "type", value = "csum") {
    
    # Perform Kruskal-Wallis test followed by a post-hoc Dunn's test
    dunn <- FSA::dunnTest(
        df[[value]] ~ df[[var]], method = "bh"
    ) |>
        purrr::pluck("res") |>
        arrange(Z)
    
    # Get compact letter display (CLD)
    cld <- rcompanion::cldList(
        comparison = dunn$Comparison,
        p.value = dunn$P.adj,
        threshold = 0.05
    )
    
    return(cld)
}



