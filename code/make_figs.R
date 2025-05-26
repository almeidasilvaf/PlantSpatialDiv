
set.seed(123)

# Load required packages ----
library(here)
library(tidyverse)
library(patchwork)
library(AngioWGD)
library(ggtree)

# Load all plots ----
files <- list.files(here("products", "plots"), full.names = TRUE)
plots <- lapply(files, readRDS)
names(plots) <- gsub("\\.rds", "", basename(files))

# Figure 1: gene-level stats ----
## Create plot
fig1 <- wrap_plots(
    plots$mean_expression_by_duplication_mode,
    plots$expression_breadth_by_duplication_mode + labs(
        subtitle = "CLD from Kruskal-Wallis + post-hoc Dunn's test; P <0.05"
    ),
    plots$heatmap_svgs_by_dupmode,
    nrow = 3
) +
    plot_annotation(tag_levels = "a")

## Save to file
ggsave(
    fig1, width = 11, height = 10,
    file = here("products", "figs", "Figure_01.pdf")
)


# Figure 2: impact of duplication mode on gene-gene correlations----
## Get a species tree with WGDs highlighted
data(tree)
data(species_metadata)
data(wgd_dates)

### Subset tree and metadata
keep <- c(
    "Arabidopsis_thaliana",
    "Brassica_rapa",
    "Glycine_max",
    "Medicago_truncatula",
    "Hordeum_vulgare",
    "Phalaenopsis_equestris",
    "Zea_mays",
    "Amborella_trichopoda"
    #"Pharus_latifolius"
)
included <- c(
    "Arabidopsis_thaliana", 
    "Glycine_max",
    "Hordeum_vulgare",
    "Zea_mays",
    "Phalaenopsis_equestris"
)
ftree <- tidytree::keep.tip(tree, keep)
fmeta <- species_metadata |>
    inner_join(data.frame(
        latin_name = keep, 
        highlight = ifelse(keep %in% included, TRUE, FALSE)
    )) |>
    mutate(species_name = str_replace_all(
        species_name, "Phalaenopsis equestris", "Phalaenopsis aphrodite"
    ))

### Plot tree
p_tree <- revts(ggtree(ftree)) %<+% fmeta +
    geom_tiplab(
        aes(label = species_name, color = highlight),
        fontface = "italic"
    ) +
    scale_color_manual(values = c("gray10", "brown3")) +
    deeptime::coord_geo(
        neg = TRUE,
        abbrv = TRUE,
        alpha = 0.6,
        expand = TRUE,
        xlim = c(-round(max(ape::node.depth.edgelength(ftree@phylo)) + 20, -1), 200),
        size = 3.5, height = unit(0.8, "cm")
    ) +
    geom_range(range = 'length_95_HPD', color = "cornflowerblue", alpha = 0.6) +
    ggtree::theme_tree2() +
    scale_y_continuous(guide = NULL) +
    scale_x_continuous(
        breaks = -c(300, 201.4, 145.0, 66.0, 23.03),
        labels = format(c(300, 201.4, 145.0, 66.0, 23.03), drop0trailing = FALSE)
    ) +
    theme(
        legend.position = "none"
    ) +
    labs(
        title = "Timetree with WGD events",
        subtitle = "Red: species analyzed in this study"
    )

p_tree

### Add WGD rectangles
place_wgds <- function(tree, wgd_dates) {
    root_id <- length(tree@phylo$tip.label) + 1
    wgds <- wgd_dates[!duplicated(wgd_dates$wgd_id), ]
    rect <- Reduce(rbind, lapply(seq_len(nrow(wgds)), function(x) {
        
        # Get WGD age estimates
        mu <- wgds$consensus_mean[x]
        hcr_min <- as.numeric(gsub("-.*", "", wgds$x90_percent_hcr[x]))
        hcr_max <- as.numeric(gsub(".*-", "", wgds$x90_percent_hcr[x]))
        
        # Get species affected by WGD in the tree
        sps <- unique(unlist(strsplit(wgds$full_species[x], ", ")))
        fsps <- sps[sps %in% tree@phylo$tip.label]
        
        # Find node
        if(length(fsps) == 1) {
            node <- which(tree@phylo$tip.label == fsps)
        } else if(length(fsps) >1) {
            node <- ape::getMRCA(tree@phylo, fsps)
        } else {
            node <- NA
        }
        
        df <- NULL
        if(!is.na(node)) {
            df <- data.frame(
                node = node, 
                wgd_id = wgds$wgd_id[x], 
                xmin = -hcr_min, xmax = -hcr_max
            )
        }
        
        return(df)
    }))
    
    return(rect)
}

wgd_ids <- c(
    "BRCE",
    "BRAS a", "BRAS b", 
    "GLYC", "PAPI",
    "ZEAM", "POAC",
    "ORCH"
)
fwgd <- wgd_dates |> filter(wgd_id %in% wgd_ids) |> distinct(wgd_id, .keep_all = TRUE)

wgd_pos <- place_wgds(ftree, fwgd) |>
    inner_join(p_tree$data) |>
    mutate(ymin = y - 0.20, ymax = y + 0.20) |>
    as.data.frame()

p_tree_final <- p_tree
for(i in seq_len(nrow(wgd_pos))) {
    bg <- AngioWGD:::make_gradient_fill(wgd_pos$xmax[i], wgd_pos$xmin[i])
    p_tree_final <- p_tree_final + geom_rect(
        data = wgd_pos[i, ], inherit.aes = FALSE, 
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = bg, color = NA
    )
}
    
## Create plot
fig2 <- wrap_plots(
    wrap_plots(
        p_tree_final, 
        plots$mean_zrho_comp, 
        widths = c(0.35, 0.65), ncol = 2
    ),
    plots$IQR_rho_per_gene_family, 
    nrow = 2
) +
    plot_annotation(tag_levels = "a")

fig2

## Save to file
ggsave(
    fig2, width = 12, height = 8,
    file = here("products", "figs", "Figure_02.pdf")
)


# Figure 3: network-based expression divergence ----
## Create plot
fig3 <- wrap_plots(
    plots$network_based_divergence + labs(y = NULL),
    plots$GCN_hubs_by_duplication_mode + labs(y = NULL),
    nrow = 2
) +
    plot_annotation(tag_levels = "a")

## Save to file
ggsave(
    fig3, width = 8, height = 6,
    file = here("products", "figs", "Figure_03.pdf")
)

# Figure 4: REB-based divergence classes ----
fig4 <- wrap_plots(
    plots$smoothed_densities_relative_expression_breadth,
    plots$ORA_dupmode_and_divergence_class +
        theme(
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            legend.position = "right"
        ),
    nrow = 2,
    heights = c(0.85, 0.15)
) +
    plot_annotation(tag_levels = "a")

## Save to file
ggsave(
    fig4, width = 10, height = 9,
    file = here("products", "figs", "Figure_04.pdf")
)

# Fig. S1: summary stats ----
sf1 <- plots$p_spe_summary_stats

sf1 <- wrap_plots(
    plots$p_spots_one_slide_per_species,
    plots$p_spe_summary_stats +
        ggh4x::facetted_pos_scales(
            x = list(
                scale_x_continuous(
                    limits = c(0, 30e3),
                    breaks = c(0, 10e3, 20e3, 30e3),
                    labels = c("0", "10K", "20K", "30K")
                ),
                scale_x_continuous(
                    limits = c(0, 1e4),
                    breaks = c(0, 2.5e3, 5e3, 7.5e3, 10e3),
                    labels = c("0", "2.5K", "5K", "7.5K", "10K")
                ),
                scale_x_continuous(limits =  c(0, 20))
            )
        )
        ,
    nrow = 2, heights = c(0.8, 0.2)
) +
    plot_annotation(
        tag_levels = list(c("a", rep("", 4), "b"))
    )

## Save to file
ggsave(
    sf1, width = 6, height = 7,
    file = here("products", "figs", "Sup_Figure_01.pdf")
)

# Fig. S2: gene-gene correlations for increasingly large metaspots ----
## Create figure
sf2 <- plots$simulation_rho_metaspots

## Save to file
ggsave(
    sf2, width = 9, height = 5,
    file = here("products", "figs", "Sup_Figure_02.pdf")
)

# Fig. S3: dissimilarity between module eigengenes ----
## Create figure
sf3 <- plots$ME_similarities_diverged_pairs

## Save to file
ggsave(
    sf3, width = 8, height = 4,
    file = here("products", "figs", "Sup_Figure_03.pdf")
)



# PDF to PNG ---------------------------------------------
pdf_files <- list.files(
    here("products", "figs"), pattern = ".pdf", full.names = TRUE
)
pdf_figs <- lapply(pdf_files, magick::image_read_pdf)
write <- lapply(seq_along(pdf_figs), function(x) {
    fname <- gsub("\\.pdf", "\\.png", pdf_files[x])
    w <- magick::image_write(pdf_figs[[x]], path = fname, density = 300)
    return(w)
})



