# Prototyping script for gene specific diversity analysis
##############
### setup
rm(list = setdiff(ls(), c("master_script.data", "commandArgs_custom")))

# Set global options
set.seed(1)
options(scipen = 6)
bg_color <- "grey97"

# Install all necessary packages
install <- F
if(install) {
    install.packages('tidyverse')
    install.packages('data.table')
    install.packages('janitor')
    install.packages('ggExtra')
    install.packages('openxlsx')
    install.packages('broom')
    install.packages('scales')
    install.packages('ggdist')
    install.packages('ggtext')
    install.packages('glue')
    install.packages('patchwork')
    install.packages('plotly')
    install.packages('flexdashboard')
    install.packages('crosstalk')
    install.packages('ComplexUpset')
}; rm(install)

# load packages, silently
suppressPackageStartupMessages({
    library(tidyverse)
    library(data.table)
    library(janitor)
    library(ggExtra)
    library(openxlsx)
    library(broom)
    library(scales)
    library(ggdist)
    library(ggtext)
    library(glue)
    library(patchwork)
    library(ComplexUpset)
    library(ggpattern)
    library(ggpubr)
    library(gridExtra)
    library(corrplot)
})

# Load custom functions
source("/Users/thomasdebruijn/Documents/PhD/R_PhD/srr_library_load_func.R")

# Set dataset names
dataset.vec <- c("bodega_all", "glacial_all",
                 "santosviromes_all", "landuse_all")

# Run data loading loop
first_iteration <- TRUE
for (dataset_i in 1:length(dataset.vec)){
    # Set data config
    #dataset_i <- 5
    dataset_id <- dataset.vec[dataset_i]
    print(dataset_id)
    commandArgs_custom <- list(
        instrain_dir = paste("instrain", dataset_id, sep = "_"),
        lifestyle_phatyp_file = paste(str_split_i(dataset_id, "_", 1), "phatyp_prediction.tsv", sep = "_"),
        lifestyle_phastyle_file = paste(str_split_i(dataset_id, "_", 1), "phastyle_prediction.tsv", sep = "_"),
        checkv_file = paste(str_split_i(dataset_id, "_", 1), "quality_summary.tsv", sep = "_"),
        coverm_file = paste("coverm_relabun", dataset_id, "vOTUs.tsv", sep = "_"),
        genomad_file = paste(str_split_i(dataset_id, "_", 1), "vOTUs_virus_summary.tsv", sep = "_"),
        wd = "/Users/thomasdebruijn/Documents/PhD/R_PhD",
        data_wd = "/Users/thomasdebruijn/Documents/PhD/DATASETS",
        install = F,
        bioproject_id = ""
    )
    sys.args <- commandArgs_custom
    
    # Set files and params to load
    wd <- sys.args$wd
    setwd(wd)
    data_wd <- sys.args$data_wd
    input.file.instrain.dir <- sys.args$instrain_dir
    input.file.phatyp.lifestyle.name <- sys.args$lifestyle_phatyp_file
    input.file.phastyle.lifestyle.name <- sys.args$lifestyle_phastyle_file
    input.file.coverm.name <- sys.args$coverm_file
    input.file.checkv.name <- sys.args$checkv_file
    input.file.genomad.name <- sys.args$genomad_file
    bioproject_ID <- sys.args$bioproject_id
    
    #####
    ### Data loading
    input.files.instrain.list <- list.files(path = paste(data_wd, input.file.instrain.dir, sep = "/"),
                                            full.names = F)
    input.files.SRR.vector <- sapply(unlist(input.files.instrain.list), SRR_extract)
    test.input.files.instrain.vec <- sapply(input.files.instrain.list, function(x){
        paste(data_wd, input.file.instrain.dir, x[1], "output", x[1], sep = "/")
    })
    
    # Load inStrain gene diversity data and clean
    raw.instrain.gene.data <- read.csv(file = paste0(test.input.files.instrain.vec[1], "_gene_info.tsv"), 
                                       sep = "\t", header = T)
    
    clean.instrain.gene.data <- raw.instrain.gene.data |>
        dplyr::mutate(Accession = str_split_i(scaffold, "_length_", i = 1), .keep = "unused", .before = 1) |>
        dplyr::group_by(Accession) # |>
    # dplyr::summarise(gene_nucl_diversity = median(nucl_diversity, na.rm = T), 
    #                  gene_dNdS_substitutions = median(dNdS_substitutions, na.rm = T))
    
    rm("raw.instrain.gene.data")
    # Load geNomad predictions and clean
    raw.genomad.data <- read.csv(file = paste(data_wd, input.file.genomad.name, sep = "/"), 
                                 sep = "\t", header = T)
    
    clean.genomad.data <- raw.genomad.data |>
        dplyr::mutate(contig_id = str_split_i(seq_name, "_length_", i = 1), .keep = "unused") |>
        dplyr::select(contig_id, taxonomy)
    
    rm("raw.genomad.data")
    # PROTOTYPING - Load raw instrain data and clean
    tmp.raw.instrain.data <- read.csv(file = paste0(test.input.files.instrain.vec[1], "_scaffold_info.tsv"), 
                                      sep = "\t", header = T)
    
    tmp.clean.instrain.data <- tmp.raw.instrain.data |>
        dplyr::mutate(Accession = str_split_i(scaffold, "_length_", i = 1), .after = 1) |>
        dplyr::mutate(Accession = str_split_i(Accession, "\\|\\|", i = 1)) |>
        dplyr::group_by(Accession) |>
        dplyr::slice_max(coverage) |>
        dplyr::ungroup() #|>
    #dplyr::select(Accession, conANI_reference, popANI_reference)
    
    rm("tmp.raw.instrain.data")
    # Load CheckV predictions and clean
    raw.checkv.data <- read.csv(file = paste(data_wd, input.file.checkv.name, sep = "/"), 
                                sep = "\t", header = T)
    
    clean.checkv.data <- raw.checkv.data |>
        dplyr::mutate(contig_id = str_split_i(contig_id, "_length_", i = 1)) |>
        dplyr::group_by(contig_id) |>
        dplyr::slice_max(completeness) |>
        dplyr::ungroup()
    rm("raw.checkv.data")
    
    # Load CoverM relative abundance data
    tmp.raw.data.relabund <- read.csv(file = paste(data_wd, input.file.coverm.name, sep = "/"),
                                      sep = "\t", header = T)
    
    tmp.clean.data.relabund <- tmp.raw.data.relabund |>
        dplyr::rename(genome_ID = 1, rel_abundance = 2) |>
        dplyr::slice(-1) |>
        dplyr::mutate(genome_ID = ifelse(grepl("single_", genome_ID), as.character(gsub("single_","", genome_ID)), NA)) |>
        dplyr::mutate(Accession = str_split_i(genome_ID, "_length_", i = 1), .keep = "unused", .before = 1) |>
        dplyr::group_by(Accession) |>
        dplyr::slice_max(rel_abundance) |>
        dplyr::ungroup()
    rm("tmp.raw.data.relabund")
    
    # Load DIAMOND blastp integrase results and clean
    raw.blast.data <- read.csv(file = paste(data_wd, paste0(str_split_i(dataset_id, "_", 1), "_against_integrases.tsv"), sep = "/"),
                               sep = "\t", header = F)
    
    clean.blast.data <- raw.blast.data |>
        purrr::set_names(c("qseqid","sseqid","qcovhsp","pident","length","mismatch",
                           "gapopen","qstart","qend","sstart","send","evalue","bitscore")) |>
        dplyr::filter(evalue < 1e-20 & qcovhsp > 80) |>
        dplyr::arrange(length, evalue) |>
        dplyr::distinct(qseqid, .keep_all = T) |>
        dplyr::mutate(qseqid = as.character(qseqid)) |>
        dplyr::mutate(Accession = str_remove(qseqid, pattern = "_\\d{1,5}$"), .before = 1) |>
        dplyr::mutate(Accession = str_split_i(Accession, "_length_", i = 1)) |>
        dplyr::distinct(Accession, .keep_all = T) |>
        dplyr::select(Accession, evalue)
    rm("raw.blast.data")
    
    # Load Phatyp predictions and clean
    raw.phatyp.lifestyle.data <- read.csv(file = paste(data_wd, input.file.phatyp.lifestyle.name, sep = "/"), 
                                          sep = "\t", header = T)
    
    clean.phatyp.lifestyle.data <- raw.phatyp.lifestyle.data |>
        dplyr::mutate(PhaTYPType = as.factor(TYPE), .keep = "unused") |>
        dplyr::mutate(Accession = str_split_i(Accession, "_length_", i = 1)) |>
        dplyr::distinct(Accession, .keep_all = T)
    rm("raw.phatyp.lifestyle.data")
    clean.phatyp.virulent.set <- dplyr::filter(clean.phatyp.lifestyle.data, PhaTYPType == "virulent")$Accession
    clean.phatyp.temperate.set <- dplyr::filter(clean.phatyp.lifestyle.data, PhaTYPType == "temperate")$Accession
    
    # Load PhaStyle predictions and clean
    raw.phastyle.lifestyle.data <- read.csv(file = paste(data_wd, input.file.phastyle.lifestyle.name, sep = "/"), 
                                            sep = "\t", header = T)
    
    clean.phastyle.lifestyle.data <- raw.phastyle.lifestyle.data |>
        dplyr::mutate(PhaStyleType = as.factor(predicted_label), .keep = "unused") |>
        dplyr::rename(Accession = fasta_id) |>
        dplyr::mutate(Accession = str_split_i(Accession, "_length_", i = 1)) |>
        dplyr::distinct(Accession, .keep_all = T) |>
        dplyr::mutate(PhaStyleScore = ifelse(PhaStyleType == "temperate", p_temperate,
                                             ifelse(PhaStyleType == "virulent", p_virulent, NA))) |>
        dplyr::select(Accession, PhaStyleType, PhaStyleScore)
    rm("raw.phastyle.lifestyle.data")
    clean.phastyle.virulent.set <- dplyr::filter(clean.phastyle.lifestyle.data, PhaStyleType == "virulent")$Accession
    clean.phastyle.temperate.set <- dplyr::filter(clean.phastyle.lifestyle.data, PhaStyleType == "temperate")$Accession
    
    # Combine lifestyle predictions
    # clean.lifestyle.data.all <- clean.phatyp.lifestyle.data |>
    #     dplyr::inner_join(clean.phastyle.lifestyle.data, by = join_by(Accession == Accession), 
    #                       keep = F, relationship = "one-to-one", multiple = "any") |>
    #     dplyr::inner_join(clean.checkv.data, by = join_by(Accession == contig_id),
    #                       keep = F, relationship = "one-to-one", multiple = "any") |>
    #     dplyr::left_join(clean.blast.data, by = join_by(Accession == Accession), keep = F,
    #                      relationship = "one-to-one") |>
    #     dplyr::left_join(tmp.clean.data.relabund, by = join_by(Accession), keep = F,
    #                      relationship = "one-to-one", multiple = "any") |>
    #     dplyr::left_join(clean.instrain.gene.data, by = join_by(Accession), keep = F,
    #                      relationship = "one-to-one", multiple = "any") |>
    #     dplyr::mutate(blasthit = ifelse(!is.na(evalue), T, F)) |>
    #     dplyr::mutate(phatyp_lytic = ifelse(PhaTYPType == "virulent", TRUE, FALSE),
    #                   phatyp_temperate = ifelse(PhaTYPType == "temperate", TRUE, FALSE),
    #                   phastyle_lytic = ifelse(PhaStyleType == "virulent", TRUE, FALSE),
    #                   phastyle_temperate = ifelse(PhaStyleType == "temperate", TRUE, FALSE)) |>
    #     dplyr::filter((phatyp_temperate == T & phastyle_temperate == T) | (phastyle_lytic == T & phatyp_lytic == T)) |>
    #     dplyr::filter((phatyp_temperate == T & blasthit == T) | (phatyp_lytic == T & blasthit == F)) |>
    #     tidyr::drop_na(completeness) |>
    #     dplyr::mutate(completeness_group = ifelse(completeness > 80, ">80%",
    #                                               ifelse(completeness >= 10, "10-80%", "<10%"))) |>
    #     dplyr::mutate(completeness_group = factor(completeness_group, levels = c(">80%","10-80%","<10%"),
    #                                               ordered = T))
    
    clean.lifestyle.data.all <- clean.instrain.gene.data |>
        dplyr::inner_join(clean.phastyle.lifestyle.data, by = join_by(Accession == Accession), 
                          keep = F, relationship = "many-to-one") |>
        dplyr::inner_join(clean.checkv.data, by = join_by(Accession == contig_id),
                          keep = F, relationship = "many-to-one") |>
        dplyr::left_join(clean.blast.data, by = join_by(Accession == Accession), keep = F,
                         relationship = "many-to-one") |>
        dplyr::left_join(tmp.clean.data.relabund, by = join_by(Accession), keep = F,
                         relationship = "many-to-one") |>
        dplyr::left_join(clean.phatyp.lifestyle.data, by = join_by(Accession), keep = F,
                         relationship = "many-to-one") |>
        dplyr::left_join(tmp.clean.instrain.data, by = join_by(Accession), keep = F,
                         relationship = "many-to-one", suffix = c("","_extra")) |>
        dplyr::left_join(clean.genomad.data, by = join_by(Accession == contig_id),
                         keep = F, relationship = "many-to-one", multiple = "any") |>
        dplyr::mutate(tailed_phage = if_else(grepl("Caudoviricetes", taxonomy), TRUE, FALSE)) |>
        dplyr::filter(tailed_phage == TRUE) |>
        dplyr::mutate(blasthit = ifelse(!is.na(evalue), T, F)) |>
        dplyr::mutate(phatyp_lytic = ifelse(PhaTYPType == "virulent", TRUE, FALSE),
                      phatyp_temperate = ifelse(PhaTYPType == "temperate", TRUE, FALSE),
                      phastyle_lytic = ifelse(PhaStyleType == "virulent", TRUE, FALSE),
                      phastyle_temperate = ifelse(PhaStyleType == "temperate", TRUE, FALSE)) |>
        dplyr::filter((phatyp_temperate == T & phastyle_temperate == T) | (phastyle_lytic == T & phatyp_lytic == T)) |>
        dplyr::filter((phatyp_temperate == T & blasthit == T) | (phatyp_lytic == T & blasthit == F)) |>
        tidyr::drop_na(completeness) |>
        dplyr::mutate(completeness_group = ifelse(completeness > 80, ">80%",
                                                  ifelse(completeness >= 10, "10-80%", "<10%"))) |>
        dplyr::mutate(completeness_group = factor(completeness_group, levels = c(">80%","10-80%","<10%"),
                                                  ordered = T))
    rm("tmp.clean.instrain.data")
    # Select PhaStyle as we pick the overlap anyway
    clean.lifestyle.data <- clean.lifestyle.data.all |>
        dplyr::select(Accession, Length, PhaStyleScore, PhaStyleType) |>
        dplyr::rename(PhaTYPScore = PhaStyleScore, Type = PhaStyleType)
    
    # Load instrain data per dataset
    tmp.instrain.data.list <- load_srr_libraries(srr_file_vector = test.input.files.instrain.vec,
                                                 lifestyle_df = clean.lifestyle.data,
                                                 checkv_df = clean.checkv.data,
                                                 thresholds = list(
                                                     breadth = 0.8,
                                                     coverage = 10,
                                                     completeness = 80
                                                 ))
    
    tmp.clean.data <- tmp.instrain.data.list[[1]]
    tmp.clean.stat.data <- tmp.instrain.data.list[[2]]
    tmp.list.clean.data <- tmp.instrain.data.list[[3]]
    tmp.raw.stat.data.all <- tmp.instrain.data.list[[4]]
    
    if(first_iteration == T){
        main.clean.data <- tmp.clean.data
        main.clean.stat.data <- tmp.clean.stat.data
        seco.clean.lifestyle.data <- clean.lifestyle.data.all
        
        first_iteration <- F
    } else {
        main.clean.data <- rbind(main.clean.data, tmp.clean.data)
        seco.clean.lifestyle.data <- rbind(seco.clean.lifestyle.data, clean.lifestyle.data.all)
        
        if (paste(colnames(main.clean.stat.data), collapse = "") == paste(colnames(tmp.clean.stat.data), collapse = "")){
            main.clean.stat.data <- rbind(main.clean.stat.data, tmp.clean.stat.data)    
        }
    }
    print(paste(paste0("PhaStyle hits: ", length(unique(sort(clean.phastyle.lifestyle.data$Accession)))),
                paste0("Blast hits: ", length(unique(sort(clean.blast.data$Accession)))),
                paste0("Tmp main data: ", length(unique(sort(tmp.clean.data$genome_ID)))),
                paste0("Total main data: ", length(unique(sort(main.clean.data$genome_ID)))),
                paste0("CoverM hits: ", length(unique(sort(tmp.clean.data.relabund$Accession)))),
                paste0("Lifestyle df: ", length(unique(sort(clean.lifestyle.data.all$Accession)))),
                sep = "   "))
    rm("tmp.clean.data", "tmp.clean.stat.data", "tmp.list.clean.data", "tmp.raw.stat.data.all",
       "clean.lifestyle.data", "clean.lifestyle.data.all","clean.checkv.data", "tmp.clean.data.relabund",
       "tmp.instrain.data.list","clean.instrain.gene.data",
       "clean.phastyle.lifestyle.data", "clean.phatyp.lifestyle.data")
}

# Conversion to one big dataframe
main.clean.data.plot <- main.clean.data |>
    #dplyr::filter(variable == "nucl_diversity") |>
    dplyr::mutate(SRR_ID_lifestyle = paste(SRR_ID, lifestyle)) |>
    dplyr::group_by(SRR_ID)

# Visualize the nucl_diversity of genes specifically by lifestyle
# SRR_TO_PLOT.data <- main.clean.data.plot |>
#     dplyr::ungroup() |>
#     tidyr::pivot_wider(names_from = variable, values_from = value) |>
#     dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, 
#                                    Accession, gene_nucl_diversity, gene_dNdS_substitutions),
#                      by = join_by(genome_ID == Accession),
#                      relationship = "one-to-one", keep = F) |>
#     tidyr::pivot_longer(cols = c("nucl_diversity", "length", 
#                                  "breadth", "coverage",
#                                  "gene_nucl_diversity","gene_dNdS_substitutions"),
#                         names_to = "variable", values_to = "value") #|>
# #dplyr::filter(variable == "gene_dNdS_substitutions") |>
# #tidyr::drop_na(value)
# SRR_TO_PLOT <- unique(sort(SRR_TO_PLOT.data$SRR_ID))
# 
# SRR_TO_PLOT.data |>
#     tidyr::pivot_wider(names_from = variable, values_from = value) |>
#     ggplot( aes(x = gene_nucl_diversity, 
#                 y = nucl_diversity,
#                 color = lifestyle)) +
#     geom_point(alpha = 0.5) +
#     facet_wrap(~lifestyle) +
#     geom_abline()
# 
# SRR_TO_PLOT.data |>
#     tidyr::pivot_wider(names_from = variable, values_from = value) |>
#     tidyr::drop_na(gene_nucl_diversity, gene_dNdS_substitutions) |>
#     ggplot( aes(x = gene_nucl_diversity, 
#                 y = gene_dNdS_substitutions,
#                 color = lifestyle)) +
#     geom_point() +
#     facet_wrap(~lifestyle)

SRR_TO_PLOT.data.wide <- main.clean.data.plot |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(seco.clean.lifestyle.data,
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-many", keep = F, suffix = c("_genome", "_gene")) |>
    dplyr::mutate(gene_length = if_else(is.na(gene_length), end - start + 1, gene_length)) |>
    dplyr::filter(coverage_gene > 10 & breadth_gene > 0.8)
cor(SRR_TO_PLOT.data.wide$nucl_diversity_genome, SRR_TO_PLOT.data.wide$nucl_diversity_gene,
    method = "spearman", use = "complete")
wilcox.test(x = dplyr::filter(SRR_TO_PLOT.data.wide, lifestyle == "temperate")$nucl_diversity_gene, 
            y = dplyr::filter(SRR_TO_PLOT.data.wide, lifestyle == "virulent")$nucl_diversity_gene,
            conf.int = T)

wilcox.test(x = dplyr::filter(SRR_TO_PLOT.data.wide, lifestyle == "temperate" & !is.na(dNdS_substitutions))$dNdS_substitutions, 
            y = dplyr::filter(SRR_TO_PLOT.data.wide, lifestyle == "virulent" & !is.na(dNdS_substitutions))$dNdS_substitutions,
            conf.int = T)

SRR_TO_PLOT.data.wide |>
    ggplot( aes(x = gene_length, colour = lifestyle)) +
    geom_density() +
    scale_x_log10()

p.bar.dnds.all <- SRR_TO_PLOT.data.wide |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID) |>
    ggplot( aes(x = study_ID, y = dNdS_substitutions, fill = lifestyle)) +
    geom_jitter(aes(colour = lifestyle), position = position_jitterdodge(), alpha = 0.7, stroke = F) +
    geom_boxplot(outliers = F) +
    geom_text(data = group_by(ungroup(SRR_TO_PLOT.data.wide), SRR_ID) |>
                  dplyr::filter(!is.na(dNdS_substitutions)) |>
                  dplyr::summarise(sum = n()) |>
                  dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                                   if_else(SRR_ID == "SRR000004", "AgriSoil",
                                                           if_else(SRR_ID == "SRR000002","BodegaBay",
                                                                   if_else(SRR_ID == "SRR000001", "Land-use",
                                                                           if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                                   if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                                .after = SRR_ID),
              aes(x = study_ID, y = 2.05, label = sum), inherit.aes = F) +
    scale_fill_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                 "below_thresh" = "lightgrey", "unknown" = "grey"),
                      breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                      labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                      name = "vOTU lifestyle") +
    scale_color_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                  "below_thresh" = "lightgrey", "unknown" = "grey"),
                       breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                       labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                       name = "vOTU lifestyle") +
    scale_y_continuous(expand = expansion(c(0.01,0.01))) +
    geom_hline(yintercept = 1, colour = "red") +
    labs(title = "(x). Distribution of S/N mutations",
         y = "dN/dS ratio",
         x = "Dataset ID") +
    theme_classic() +
    theme(legend.position = "right",
          legend.background = element_rect(fill = "white"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"),
          axis.title.x = element_blank())

# ggsave(plot = p.bar.dnds.all,
#        filename = paste(sys.args$data_wd, "genes_dnds_all.png", sep = "/"),
#        width = 700, height = 500,
#        units = "px", dpi = 400, scale = 2.5)

##### First collection plot, dNds stuff ####
pnps.stat.data.clean <- SRR_TO_PLOT.data.wide |>
    base::split(~SRR_ID, drop = T) |>
    purrr::map(~ wilcox.test(pNpS_variants ~ lifestyle, data = ., conf.int = T)) |>
    purrr::map_dfr( ~ broom::tidy(.), .id = "SRR_ID") |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID) |>
    dplyr::rename(wilcox_pvalue = p.value, wilcox_estimate = estimate) |>
    dplyr::mutate(markup = ifelse(wilcox_pvalue > 0.05, "",
                                  ifelse(wilcox_pvalue > 0.01, "*",
                                         ifelse(wilcox_pvalue > 0.001, "**",
                                                ifelse(wilcox_pvalue > 0.0001, "***",
                                                       ifelse(wilcox_pvalue < 0.0001, "****", NA)))))) |>
    dplyr::ungroup() |>
    dplyr::mutate(markup_estimate = ifelse(!markup %in% c("",NA), wilcox_estimate, ""))

p.bar.pnps.all <- SRR_TO_PLOT.data.wide |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID) |>
    dplyr::group_by(study_ID, SRR_ID, lifestyle, genome_ID) |>
    dplyr::summarise(pNpS_variants_median = median(pNpS_variants, na.rm = T), .groups = "keep") |>
    ggplot( aes(x = study_ID, y = pNpS_variants_median, fill = lifestyle)) +
    geom_jitter(aes(colour = lifestyle), position = position_jitterdodge(), alpha = 0.7, stroke = F) +
    geom_boxplot(outliers = F) +
    geom_text(data = SRR_TO_PLOT.data.wide |>
                  dplyr::group_by(SRR_ID, lifestyle, genome_ID) |>
                  dplyr::summarise(pNpS_variants_median = median(pNpS_variants, na.rm = T), .groups = "keep") |>
                  dplyr::filter(!is.na(pNpS_variants_median)) |>
                  dplyr::group_by(SRR_ID, lifestyle) |>
                  dplyr::summarise(sum = n(), .groups = "keep") |>
                  dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                                   if_else(SRR_ID == "SRR000004", "AgriSoil",
                                                           if_else(SRR_ID == "SRR000002","BodegaBay",
                                                                   if_else(SRR_ID == "SRR000001", "Land-use",
                                                                           if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                                   if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                                .after = SRR_ID),
              aes(x = study_ID, y = 40, label = sum, group = lifestyle, hjust = 0.5, angle = 30), 
              inherit.aes = F, position = position_dodge(width = 0.9), color = "black") +
    geom_text(inherit.aes = F, data = pnps.stat.data.clean,
              aes(x = study_ID, y = 11, label = markup), hjust = 0.5, nudge_x = 0) +
    scale_fill_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                 "below_thresh" = "lightgrey", "unknown" = "grey"),
                      breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                      labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                      name = "vOTU lifestyle") +
    scale_color_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                  "below_thresh" = "lightgrey", "unknown" = "grey"),
                       breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                       labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                       name = "vOTU lifestyle") +
    scale_y_log10() +
    geom_hline(yintercept = 1, colour = "red") +
    labs(title = "Distribution of pN/pS mutations (median value)",
         y = "pN/pS ratio",
         x = "Dataset ID") +
    theme_classic() +
    theme(legend.position = "right",
          legend.background = element_rect(fill = "white"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"),
          axis.title.x = element_blank())
tmp.genes.legend <- get_legend(p.bar.pnps.all)
p.bar.pnps.all <- p.bar.pnps.all +
    theme(legend.position = "none")

# ggsave(plot = p.bar.pnps.all,
#        filename = paste(sys.args$data_wd, "genes_pnps_all_median.png", sep = "/"),
#        width = 700, height = 500,
#        units = "px", dpi = 400, scale = 2.5)

p.bar.pnps.all.mean <- SRR_TO_PLOT.data.wide |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID) |>
    dplyr::group_by(study_ID, SRR_ID, lifestyle, genome_ID) |>
    dplyr::summarise(pNpS_variants_median = mean(pNpS_variants, na.rm = T), .groups = "keep") |>
    ggplot( aes(x = study_ID, y = pNpS_variants_median, fill = lifestyle)) +
    geom_jitter(aes(colour = lifestyle), position = position_jitterdodge(), alpha = 0.7, stroke = F) +
    geom_boxplot(outliers = F) +
    geom_text(data = SRR_TO_PLOT.data.wide |>
                  dplyr::group_by(SRR_ID, lifestyle, genome_ID) |>
                  dplyr::summarise(pNpS_variants_median = mean(pNpS_variants, na.rm = T), .groups = "keep") |>
                  dplyr::filter(!is.na(pNpS_variants_median)) |>
                  dplyr::group_by(SRR_ID, lifestyle) |>
                  dplyr::summarise(sum = n(), .groups = "keep") |>
                  dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                                   if_else(SRR_ID == "SRR000004", "AgriSoil",
                                                           if_else(SRR_ID == "SRR000002","BodegaBay",
                                                                   if_else(SRR_ID == "SRR000001", "Land-use",
                                                                           if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                                   if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                                .after = SRR_ID),
              aes(x = study_ID, y = 40, label = sum, group = lifestyle, hjust = 0.5, angle = 30), 
              inherit.aes = F, position = position_dodge(width = 0.9), color = "black") +
    scale_fill_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                 "below_thresh" = "lightgrey", "unknown" = "grey"),
                      breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                      labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                      name = "vOTU lifestyle") +
    scale_color_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                  "below_thresh" = "lightgrey", "unknown" = "grey"),
                       breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                       labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                       name = "vOTU lifestyle") +
    scale_y_log10() +
    geom_hline(yintercept = 1, colour = "red") +
    labs(title = "(x). Distribution of pN/pS mutations (average value)",
         y = "pN/pS ratio",
         x = "Dataset ID") +
    theme_classic() +
    theme(legend.position = "right",
          legend.background = element_rect(fill = "white"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"),
          axis.title.x = element_blank())

# ggsave(plot = p.bar.pnps.all.mean,
#        filename = paste(sys.args$data_wd, "genes_pnps_all_mean.png", sep = "/"),
#        width = 700, height = 500,
#        units = "px", dpi = 400, scale = 2.5)

p.bar.dnds.all <- SRR_TO_PLOT.data.wide |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID) |>
    ggplot( aes(x = study_ID, y = dNdS_substitutions, fill = lifestyle)) +
    geom_jitter(aes(colour = lifestyle), position = position_jitterdodge(), alpha = 0.7, stroke = F) +
    geom_boxplot(outliers = F) +
    geom_text(data = group_by(ungroup(SRR_TO_PLOT.data.wide), SRR_ID, lifestyle) |>
                  dplyr::filter(!is.na(dNdS_substitutions)) |>
                  dplyr::summarise(sum = n(), .groups = "keep") |>
                  dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                                   if_else(SRR_ID == "SRR000004", "AgriSoil",
                                                           if_else(SRR_ID == "SRR000002","BodegaBay",
                                                                   if_else(SRR_ID == "SRR000001", "Land-use",
                                                                           if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                                   if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                                .after = SRR_ID),
              aes(x = study_ID, y = 1.75, label = sum, group = lifestyle, hjust = 0.5, angle = 0), 
              inherit.aes = F, position = position_dodge(width = 0.9), color = "black") +
    scale_fill_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                 "below_thresh" = "lightgrey", "unknown" = "grey"),
                      breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                      labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                      name = "vOTU lifestyle") +
    scale_color_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                  "below_thresh" = "lightgrey", "unknown" = "grey"),
                       breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                       labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                       name = "vOTU lifestyle") +
    #scale_y_log10() +
    geom_hline(yintercept = 1, colour = "red") +
    labs(title = "(b). Distribution of dN/dS mutations",
         y = "dN/dS ratio") +
    theme_classic() +
    theme(legend.position = "none",
          legend.background = element_rect(fill = "white"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"),
          axis.title.x = element_blank())

p.box.gene.gene_count.all <- SRR_TO_PLOT.data.wide |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID) |>
    dplyr::group_by(study_ID, lifestyle, genome_ID) |>
    dplyr::distinct(genome_ID, .keep_all = T) |>
    ggplot( aes(x = study_ID, y = gene_count, fill = lifestyle)) +
    geom_jitter(position = position_jitterdodge(), alpha = 0.3, aes(colour = lifestyle), stroke = F) +
    geom_boxplot(outliers = F) +
    geom_text(data = group_by(ungroup(SRR_TO_PLOT.data.wide), SRR_ID, lifestyle) |>
                  dplyr::distinct(genome_ID, .keep_all = T) |>
                  dplyr::summarise(sum = n(), .groups = "keep") |>
                  dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                                   if_else(SRR_ID == "SRR000004", "AgriSoil",
                                                           if_else(SRR_ID == "SRR000002","BodegaBay",
                                                                   if_else(SRR_ID == "SRR000001", "Land-use",
                                                                           if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                                   if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                                .after = SRR_ID),
              aes(x = study_ID, y = 750, label = sum, group = lifestyle, hjust = 0.5, angle = 30), 
              inherit.aes = F, position = position_dodge(width = 0.9), color = "black") +
    scale_fill_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                 "below_thresh" = "lightgrey", "unknown" = "grey"),
                      breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                      labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                      name = "vOTU lifestyle") +
    scale_color_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                  "below_thresh" = "lightgrey", "unknown" = "grey"),
                       breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                       labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                       name = "vOTU lifestyle") +
    labs(title = "(c). Number of genes on contig",
         y = "Number of genes",
         x = "Dataset ID") +
    theme_classic() +
    theme(legend.position = "none",
          legend.background = element_rect(fill = "white"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"))

p.box.gene.host_gene_count.all <- SRR_TO_PLOT.data.wide |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID) |>
    dplyr::group_by(study_ID, lifestyle, genome_ID) |>
    dplyr::distinct(genome_ID, .keep_all = T) |>
    ggplot( aes(x = study_ID, y = host_genes, fill = lifestyle)) +
    geom_jitter(position = position_jitterdodge(), alpha = 0.3, aes(color = lifestyle), stroke = F) +
    geom_boxplot(outliers = F) +
    scale_fill_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                 "below_thresh" = "lightgrey", "unknown" = "grey"),
                      breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                      labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                      name = "vOTU lifestyle") +
    scale_color_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                  "below_thresh" = "lightgrey", "unknown" = "grey"),
                       breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                       labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                       name = "vOTU lifestyle") +
    labs(title = "(d). Number of host genes detected",
         y = "Number of host genes",
         x = "Dataset ID") +
    theme_classic() +
    theme(legend.position = "none",
          legend.background = element_rect(fill = "white"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"))

p_collected_genes <- gridExtra::grid.arrange(p.bar.pnps.all,
                                             p.box.gene.gene_count.all,
                                             p.bar.dnds.all,
                                             p.box.gene.host_gene_count.all,
                                             tmp.genes.legend, 
                                             top = NULL,
                                             ncol = 3, nrow = 2,
                                             layout_matrix = cbind(c(1,2),c(3,4),c(5,5)),
                                             widths = c(2.5,2.5,.7))

# ggsave(plot = p_collected_genes,
#        filename = paste(sys.args$data_wd, "dNdS_combination_plot_genes_all.png", sep = "/"),
#        width = 2000, height = 1500,
#        units = "px", dpi = 400, scale = 1.5)

#####

SRR_TO_PLOT.data.wide |>
    base::split(~SRR_ID, drop = T) |>
    purrr::map(~ wilcox.test(gene_count ~ lifestyle, data = ., conf.int = T)) |>
    purrr::map_dfr( ~ broom::tidy(.), .id = "SRR_ID") |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID)

SRR_TO_PLOT.data.wide |>
    base::split(~SRR_ID, drop = T) |>
    purrr::map(~ wilcox.test(nucl_diversity_gene ~ lifestyle, data = ., conf.int = T)) |>
    purrr::map_dfr( ~ broom::tidy(.), .id = "SRR_ID") |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID)

SRR_TO_PLOT.data.wide |>
    #dplyr::filter(SRR_ID %in% c("SRR000001","SRR000006")) |>
    base::split(~SRR_ID, drop = T) |>
    purrr::map(~ wilcox.test(pNpS_variants ~ lifestyle, data = ., conf.int = T)) |>
    purrr::map_dfr( ~ broom::tidy(.), .id = "SRR_ID") |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID)

SRR_TO_PLOT.data.wide |>
    #dplyr::filter(SRR_ID %in% c("SRR000001","SRR000006")) |>
    dplyr::group_by(SRR_ID, lifestyle, genome_ID) |>
    dplyr::summarise(pNpS_variants_median = median(pNpS_variants, na.rm = T), .groups = "keep") |>
    base::split(~SRR_ID, drop = T) |>
    purrr::map(~ wilcox.test(pNpS_variants_median ~ lifestyle, data = ., conf.int = T)) |>
    purrr::map_dfr( ~ broom::tidy(.), .id = "SRR_ID") |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID)

cor.data.virulent <- cor(dplyr::filter(SRR_TO_PLOT.data.wide, lifestyle == "virulent" & SRR_ID == "SRR000006") |>
                             dplyr::distinct(genome_ID, .keep_all = T) |>
                             dplyr::select(length_genome, coverage_genome, breadth_genome, completeness, nucl_diversity_genome,
                                           conANI_reference, popANI_reference, divergent_site_count_extra, SNV_count_extra, SNS_count_extra,
                                           consensus_divergent_sites, population_divergent_sites),
                         method = "spearman")
cor.data.temperate <- cor(dplyr::filter(SRR_TO_PLOT.data.wide, lifestyle == "temperate" & SRR_ID == "SRR000006") |>
                              dplyr::distinct(genome_ID, .keep_all = T) |>
                              dplyr::select(length_genome, coverage_genome, breadth_genome, completeness, nucl_diversity_genome,
                                            conANI_reference, popANI_reference, divergent_site_count_extra, SNV_count_extra, SNS_count_extra,
                                            consensus_divergent_sites, population_divergent_sites),
                          method = "spearman")

cor.data.all <- cor(dplyr::distinct(SRR_TO_PLOT.data.wide, genome_ID, .keep_all = T) |>
                        dplyr::select(length_genome, coverage_genome, breadth_genome, completeness, nucl_diversity_genome,
                                      conANI_reference, popANI_reference, divergent_site_count_extra, SNV_count_extra, SNS_count_extra,
                                      consensus_divergent_sites, population_divergent_sites),
                    method = "spearman")

corrplot(cor.data.temperate, method = "number")
corrplot(cor.data.virulent, method = "number")
corrplot(cor.data.all, method = "number", type = "upper", tl.srt = 45)
#p.corr.all.v.all <- recordPlot()

# ggsave(plot = replayPlot(p.corr.all.v.all),
#        filename = paste(sys.args$data_wd, "cor_all_v_all.png", sep = "/"),
#        width = 500, height = 500,
#        units = "px")

##### Second collection plot, general statistics ####
p.bar.lifestyle_total.all.per_study.per_lifestyle <- SRR_TO_PLOT.data.wide |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID) |>
    dplyr::filter(completeness >= 80 & coverage_genome >= 10 & breadth_genome > 0.8) |>
    dplyr::group_by(study_ID, lifestyle) |>
    ggplot( aes(x = study_ID, fill = lifestyle)) +
    geom_bar(position = "fill", alpha = 1, width = 0.3) +
    geom_text(data = dplyr::mutate(SRR_TO_PLOT.data.wide, study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                                                             if_else(SRR_ID == "SRR000004", "AgriSoil",
                                                                                     if_else(SRR_ID == "SRR000002","BodegaBay",
                                                                                             if_else(SRR_ID == "SRR000001", "Land-use",
                                                                                                     if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                                                             if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                                   .after = SRR_ID) |>
                  dplyr::filter(completeness >= 80 & coverage_genome >= 10 & breadth_genome > 0.8) |>
                  #ungroup() |>
                  dplyr::group_by(SRR_ID, study_ID) |>
                  dplyr::summarise(sum = n()),
              aes(x = study_ID, y = 1.09, label = sum), inherit.aes = F, vjust = 1) +
    scale_fill_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                 "below_thresh" = "lightgrey", "unknown" = "grey",
                                 "used" = "black", "discarded" = "grey"),
                      breaks = c("temperate", "virulent", "below_thresh", "unknown", "used", "discarded"),
                      labels = c("temperate", "virulent", "below quality\nthreshold", 
                                 "unknown lifestyle", "used", "discarded"),
                      name = "vOTU\nlifestyle") +
    scale_y_continuous(labels = scales::label_percent(), expand = expansion(c(0,0.1), 0), breaks = c(0, 0.25, 0.50, 0.75, 1)) +
    labs(title = "(a). Genes per dataset",
         y = "Percentage of genes",
         x = "Dataset ID") +
    theme_classic() +
    theme(legend.position = "right",
          legend.background = element_rect(fill = "white"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"),
          axis.title.x = element_blank())
tmp.lifestyle.legend <- get_legend(p.bar.lifestyle_total.all.per_study.per_lifestyle)
p.bar.lifestyle_total.all.per_study.per_lifestyle <- p.bar.lifestyle_total.all.per_study.per_lifestyle +
    theme(legend.position = "none")

# Boxplot of length per lifestyle per study 
p.box.length.all.per_study.per_lifestyle <- SRR_TO_PLOT.data.wide |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID) |>
    dplyr::filter(completeness >= 80 & coverage_genome >= 10 & breadth_genome > 0.8) |>
    dplyr::group_by(study_ID, lifestyle, genome_ID) |>
    dplyr::summarise(median_gene_length = median(gene_length), .groups = "keep") |>
    ggplot( aes(x = study_ID, y = median_gene_length, fill = lifestyle)) +
    geom_jitter(position = position_jitterdodge(), alpha = 0.3, aes(color = lifestyle), stroke = F) +
    geom_boxplot(outliers = F) +
    scale_fill_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                 "below_thresh" = "lightgrey", "unknown" = "grey"),
                      breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                      labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                      name = "vOTU lifestyle") +
    scale_color_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                  "below_thresh" = "lightgrey", "unknown" = "grey"),
                       breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                       labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                       name = "vOTU lifestyle") +
    scale_y_continuous(expand = expansion(c(0, 0.01), 0), limits = c(0, 1200)) +
    labs(title = "(b). Median gene length",
         y = "Median gene length \nper vOTU (nucleotides)",
         x = "Dataset ID") +
    theme_classic() +
    theme(legend.position = "none",
          legend.background = element_rect(fill = "grey93"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"),
          axis.title.x = element_blank())

# Boxplot of coverage per lifestyle per study 
p.box.coverage.all.per_study.per_lifestyle <- SRR_TO_PLOT.data.wide |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID) |>
    dplyr::filter(completeness >= 80 & coverage_genome >= 10 & breadth_genome > 0.8) |>
    dplyr::group_by(study_ID, lifestyle, genome_ID) |>
    dplyr::summarise(median_coverage_gene = median(coverage_gene), .groups = "keep") |>
    ggplot( aes(x = study_ID, y = median_coverage_gene, fill = lifestyle)) +
    geom_jitter(position = position_jitterdodge(), alpha = 0.3, aes(color = lifestyle), stroke = F) +
    geom_boxplot(outliers = F) +
    scale_fill_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                 "below_thresh" = "lightgrey", "unknown" = "grey"),
                      breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                      labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                      name = "vOTU lifestyle") +
    scale_color_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                  "below_thresh" = "lightgrey", "unknown" = "grey"),
                       breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                       labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                       name = "vOTU lifestyle") +
    scale_y_continuous(limits = c(0,2500), expand = expansion(c(0, 0))) +
    labs(title = "(d). Median gene coverage",
         y = "\nCoverage of gene",
         x = "Dataset ID") +
    #geom_hline(yintercept = 10, color = "red") +
    theme_classic() +
    theme(legend.position = "none",
          legend.background = element_rect(fill = "grey93"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"))

# Boxplot of completeness per lifestyle per study 
p.box.completeness.all.per_study.per_lifestyle <- SRR_TO_PLOT.data.wide |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID) |>
    dplyr::filter(completeness >= 80 & coverage_genome >= 10 & breadth_genome > 0.8) |>
    dplyr::group_by(study_ID, lifestyle, genome_ID) |>
    dplyr::summarise(median_breadth_gene = median(breadth_minCov), .groups = "keep") |>
    dplyr::mutate(median_breadth_gene = median_breadth_gene * 100) |>
    ggplot( aes(x = study_ID, y = median_breadth_gene, fill = lifestyle)) +
    geom_jitter(position = position_jitterdodge(), alpha = 0.3, aes(color = lifestyle), stroke = F) +
    geom_boxplot(outliers = F, position = "dodge") +
    geom_text(data = dplyr::mutate(SRR_TO_PLOT.data.wide, study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                                                             if_else(SRR_ID == "SRR000004", "AgriSoil",
                                                                                     if_else(SRR_ID == "SRR000002","BodegaBay",
                                                                                             if_else(SRR_ID == "SRR000001", "Land-use",
                                                                                                     if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                                                             if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                                   .after = SRR_ID) |>
                  dplyr::filter(completeness >= 80 & coverage_genome >= 10 & breadth_genome > 0.8) |>
                  #ungroup() |>
                  dplyr::group_by(study_ID, lifestyle) |>
                  dplyr::summarise(sum = n()),
              aes(x = study_ID, y = 90, label = sum), inherit.aes = F, vjust = 1, position = position_dodge2(0.9), angle = 45) +
    scale_fill_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                 "below_thresh" = "lightgrey", "unknown" = "grey"),
                      breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                      labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                      name = "vOTU lifestyle") +
    scale_color_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                  "below_thresh" = "lightgrey", "unknown" = "grey"),
                       breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                       labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                       name = "vOTU lifestyle") +
    labs(title = "(c). Median breadth of gene per vOTU",
         y = "Gene breadth",
         x = "Dataset ID") +
    scale_y_continuous(labels = scales::label_percent(scale = 1), limits = c(0, 100),
                       expand = expansion(0, c(0,0))) +
    theme_classic() +
    theme(legend.position = "none",
          legend.background = element_rect(fill = "grey93"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"))

# Assemble above graphs into a single block using gridExtra
p_collected_population <- gridExtra::grid.arrange(p.bar.lifestyle_total.all.per_study.per_lifestyle,
                                                  p.box.completeness.all.per_study.per_lifestyle,
                                                  p.box.length.all.per_study.per_lifestyle,
                                                  p.box.coverage.all.per_study.per_lifestyle,
                                                  tmp.lifestyle.legend, 
                                                  top = NULL,
                                                  ncol = 3, nrow = 2,
                                                  layout_matrix = cbind(c(1,2),c(3,4),c(5,5)),
                                                  widths = c(2.5,2.5,.7))

# ggsave(plot = p_collected_population,
#        filename = paste(sys.args$data_wd, "lifestyle_plot_genes_all.png", sep = "/"),
#        width = 2000, height = 1500,
#        units = "px", dpi = 400, scale = 1.5)

##### Making line graphs of relationships between genome metrics and gene metrics ####
p.point.microdiversity.gene_vs_genome.all <- SRR_TO_PLOT.data.wide |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID) |>
    dplyr::filter(completeness >= 80 & coverage_genome >= 10 & breadth_genome > 0.8) |>
    dplyr::group_by(study_ID, lifestyle, genome_ID) |>
    tidyr::drop_na(nucl_diversity_gene) |>
    ggplot( aes(x = nucl_diversity_genome, y = nucl_diversity_gene, colour = lifestyle)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method=lm , color="red", se=FALSE) +
    geom_abline() +
    scale_x_log10(expand = expansion(c(0,0)), limits = c(0.0001,0.1), n.breaks = 5) +
    scale_y_log10(expand = expansion(c(0,0)), limits = c(0.0001,0.1)) +
    labs(title = "(S.figure.X). Relation gene vs genome microdiversity",
         y = "Gene microdiversity (π)",
         x = "Genome microdiversity (π)") +
    theme_classic() +
    theme(legend.position = "right",
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"))

p.point.coverage.gene_vs_genome.all <- SRR_TO_PLOT.data.wide |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID) |>
    dplyr::filter(completeness >= 80 & coverage_genome >= 10 & breadth_genome > 0.8) |>
    dplyr::group_by(study_ID, lifestyle, genome_ID) |>
    tidyr::drop_na(nucl_diversity_gene) |>
    ggplot( aes(x = coverage_genome, y = coverage_gene, colour = lifestyle)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method=lm , color="red", se=FALSE) +
    geom_abline() +
    scale_x_log10(expand = expansion(c(0,0))) +
    scale_y_log10(expand = expansion(c(0,0))) +
    labs(title = "(S.figure.X). Relation gene vs genome coverage",
         y = "Gene coverage",
         x = "Genome coverage") +
    theme_classic() +
    theme(legend.position = "right",
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"))

SRR_TO_PLOT.data.wide |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID) |>
    dplyr::filter(completeness >= 80 & coverage_genome >= 10 & breadth_genome > 0.8) |>
    dplyr::group_by(study_ID, lifestyle, genome_ID) |>
    tidyr::drop_na(nucl_diversity_gene) |>
    ggplot( aes(x = breadth_genome, y = breadth_gene, colour = lifestyle)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method=lm , color="red", se=FALSE) +
    geom_abline() +
    scale_x_log10(expand = expansion(c(0,0))) +
    scale_y_log10(expand = expansion(c(0,0))) +
    labs(title = "(S.figure.X). Relation gene vs genome coverage",
         y = "Gene coverage",
         x = "Genome coverage") +
    theme_classic() +
    theme(legend.position = "right",
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"))

