

#' Aggregate spots into metaspots
#'
#' @param spe A `SpatialExperiment` object.
#' @param spatial_domain Character indicating the name of the colData variable
#' containing information on spatial domains (or cell types).
#' @param nspots Numeric indicating how many spots each meta-spot should
#' have (tentatively). Default: 7.
#' @param assay Character indicating the name of the assay to use for
#' aggregation. Default: "counts".
#'
#' @return A `SpatialExperiment` object with aggregated counts.
#' @details
#' For each spatial domain, this function creates 'meta-spots' by using k-means
#' clustering with k = total number of spots / `nspots`. Then, counts are
#' aggregated bu summing the counts in each meta-spot.
#' @noRd
#'
aggregate_metaspots <- function(
        spe, spatial_domain, nspots = 7, assay = "counts", ...
) {
    
    # For each spatial domain, get meta spots containing `nspots` spots
    domains <- as.character(unique(spe[[spatial_domain]]))
    metaspots_coldata <- Reduce(rbind, lapply(domains, function(x) {
        
        fspe <- spe[, spe[[spatial_domain]] == x]
        
        # Get metaspot IDs with k-means clustering
        coords <- spatialCoords(fspe)
        km <- setNames(rep(1, nrow(coords)), rownames(coords))
        if(nrow(coords) > nspots) {
            km <- kmeans(
                coords, centers = round(nrow(coords) / nspots), ...
            )$cluster
        }
        fspe$metaspot_id <- paste(fspe[[spatial_domain]], km, sep = "_")
        
        return(colData(fspe))
    }))
    
    colData(spe) <- metaspots_coldata
    
    # Aggregate counts to metaspots
    metaspots <- unique(spe$metaspot_id)
    counts_metaspots <- Reduce(cbind, lapply(metaspots, function(x) {
        
        fspe <- spe[, spe$metaspot_id == x]
        ag_counts <- Matrix::Matrix(
            rowSums(assay(fspe, assay)),
            dimnames = list(rownames(fspe), x)
        )
        
        return(ag_counts)
    }))
    
    # Build the `SpatialExperiment` object again
    coldata <- metaspots_coldata[!duplicated(metaspots_coldata$metaspot_id), ]
    rownames(coldata) <- coldata$metaspot_id
    
    final_spe <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts_metaspots),
        rowData = rowData(spe),
        colData = coldata
    )

    return(final_spe)
}


#' Aggregate spots by spatial domain (a.k.a. pseudobulk)
#'
#' @param spe A `SpatialExperiment` object.
#' @param spatial_domain Character indicating the name of the colData variable
#' containing information on spatial domains (or cell types).
#' @param assay Character indicating the name of the assay to use for
#' aggregation. Default: "counts".
#'
#' @return A `SpatialExperiment` object with aggregated counts.
#' @noRd
#'
aggregate_pseudobulk <- function(
        spe, spatial_domain, assay = "counts"
) {
    
    # Aggregate counts using pseudobulk
    domains <- as.character(unique(spe[[spatial_domain]]))
    counts_pseudobulk <- Reduce(cbind, lapply(domains, function(x) {

        ## Keep only spots for domain {{x}}
        fspe <- spe        
        fspe[[spatial_domain]] <- as.character(fspe[[spatial_domain]])
        fspe <- fspe[, fspe[[spatial_domain]] == x]
        
        ## Aggregate counts
        ag_counts <- Matrix::Matrix(
            rowSums(assay(fspe, assay)),
            dimnames = list(rownames(fspe), x)
        )
        
        return(ag_counts)
    }))
    colnames(counts_pseudobulk) <- paste0("sd_", colnames(counts_pseudobulk))
    
    # Rebuild the `SpatialExperiment` object
    coldata <- data.frame(
        row.names = colnames(counts_pseudobulk), 
        domain = colnames(counts_pseudobulk)
    )
    final_spe <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts_pseudobulk),
        rowData = rowData(spe),
        colData = coldata
    )

    return(final_spe)
}



#' Aggregate counts using metaspots or pseudobulk
#' 
#' @param spe A `SpatialExperiment` object.
#' @param spatial_domain Character indicating the name of the colData variable
#' containing information on spatial domains (or cell types).
#' @param method Character indicating the aggregation method to use.
#' One of 'metaspot' or 'pseudobulk'. Default: 'metaspot'.
#' @param nspots Numeric indicating how many spots each meta-spot should
#' have (tentatively). Default: 7.
#' @param assay Character indicating the name of the assay to use for
#' aggregation. Default: "counts".
#' @param add_logcounts Logical indicating whether to add a `logcounts` assay
#' by computing size factors from aggregated counts and normalizing by
#' library size with `scater::logNormCounts()`. Default: TRUE.
#' @param keep_raw Logical indicating whether to keep the assay with raw 
#' aggregated counts. Default: FALSE.
#'
#' @return A `SpatialExperiment` object with aggregated counts.
#' 
#' @details
#' For each spatial domain, this function creates 'meta-spots' by using k-means
#' clustering with k = total number of spots / `nspots`. Then, counts are
#' aggregated bu summing the counts in each meta-spot.
#'
#'
aggregate_counts <- function(
        spe, spatial_domain, method = "metaspot",
        assay = "counts", nspots = 7,  
        add_logcounts = TRUE, keep_raw = FALSE, ...
) {
    
    # Aggregate counts
    if(method == "metaspot") {
        final_spe <- aggregate_metaspots(spe, spatial_domain, nspots, assay, ...)
    } else {
        final_spe <- aggregate_pseudobulk(spe, spatial_domain, assay)
    }

    if(add_logcounts) {
        final_spe <- scater::computeLibraryFactors(final_spe)
        final_spe <- scater::logNormCounts(final_spe)
    }
    
    if(keep_raw == FALSE) {
        assay(final_spe, assay) <- NULL
    }
    
    return(final_spe)
}
