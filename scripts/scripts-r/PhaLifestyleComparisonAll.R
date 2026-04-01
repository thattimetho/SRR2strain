# Analyse lifestyle predictions of all-at-once mapping
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
    library(grid)
})

# Load custom functions
source("/Users/thomasdebruijn/Documents/PhD/R_PhD/srr_library_load_func.R")

# Set dataset names
dataset.vec <- c("glacial_all","santosviromes_all", "landuse_all", 
                 "wildfire_all", "hopland_all", "wetup_all", "ncalhabitat_all", 
                 "bodega_all", "permathaw_all", "newmexico_all",
                 "medigrass_all","conifer_all")
# "spruce_all","rhizo_all","ukagri_all"

# Run data loading loop
first_iteration <- TRUE
for (dataset_i in 1:length(dataset.vec)){
    # Set data config
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
    
    # Load geNomad predictions and clean
    raw.genomad.data <- read.csv(file = paste(data_wd, input.file.genomad.name, sep = "/"), 
                                 sep = "\t", header = T)
    
    clean.genomad.data <- raw.genomad.data |>
        dplyr::mutate(contig_id = str_split_i(seq_name, "_length_", i = 1), .keep = "unused") |>
        dplyr::select(contig_id, taxonomy)
    
    rm("raw.genomad.data")
    # Load CheckV predictions and clean
    raw.checkv.data <- read.csv(file = paste(data_wd, input.file.checkv.name, sep = "/"), 
                                sep = "\t", header = T)
    
    clean.checkv.data <- raw.checkv.data |>
        dplyr::mutate(contig_id = str_split_i(contig_id, "_length_", i = 1))
    rm("raw.checkv.data")
    
    # Load CoverM relative abundance data
    tmp.raw.data.relabund <- read.csv(file = paste(data_wd, input.file.coverm.name, sep = "/"),
                                      sep = "\t", header = T)
    
    tmp.clean.data.relabund <- tmp.raw.data.relabund |>
        dplyr::rename(genome_ID = 1, rel_abundance = 2) |>
        #dplyr::slice(-1) |>
        dplyr::mutate(genome_ID = ifelse(grepl("single_", genome_ID), as.character(gsub("single_","", genome_ID)), genome_ID)) |>
        dplyr::mutate(Accession = str_split_i(genome_ID, "_length_", i = 1), .keep = "unused")
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
    clean.lifestyle.data.all <- clean.phatyp.lifestyle.data |>
        dplyr::inner_join(clean.phastyle.lifestyle.data, by = join_by(Accession == Accession), 
                          keep = F, relationship = "one-to-one", multiple = "any") |>
        dplyr::inner_join(clean.checkv.data, by = join_by(Accession == contig_id),
                          keep = F, relationship = "one-to-one", multiple = "any") |>
        dplyr::left_join(clean.blast.data, by = join_by(Accession == Accession), keep = F,
                         relationship = "one-to-one") |>
        dplyr::left_join(tmp.clean.data.relabund, by = join_by(Accession), keep = F,
                         relationship = "one-to-one", multiple = "any") |>
        dplyr::left_join(clean.genomad.data, by = join_by(Accession == contig_id),
                         keep = F, relationship = "one-to-one", multiple = "any") |>
        dplyr::mutate(tailed_phage = if_else(grepl("Caudoviricetes", taxonomy), TRUE, FALSE)) |>
        dplyr::mutate(blasthit = ifelse(!is.na(evalue), T, F)) |>
        dplyr::mutate(phatyp_lytic = ifelse(PhaTYPType == "virulent", TRUE, FALSE),
                      phatyp_temperate = ifelse(PhaTYPType == "temperate", TRUE, FALSE),
                      phastyle_lytic = ifelse(PhaStyleType == "virulent", TRUE, FALSE),
                      phastyle_temperate = ifelse(PhaStyleType == "temperate", TRUE, FALSE)) |>
        dplyr::filter((phatyp_temperate == T & phastyle_temperate == T) | (phastyle_lytic == T & phatyp_lytic == T)) |>
        dplyr::filter((phatyp_temperate == T & blasthit == T) | (phatyp_lytic == T & blasthit == F)) |>
        dplyr::filter(tailed_phage == TRUE) |>
        tidyr::drop_na(completeness) |>
        dplyr::mutate(completeness_group = ifelse(completeness > 80, ">80%",
                                                  ifelse(completeness >= 10, "10-80%", "<10%"))) |>
        dplyr::mutate(completeness_group = factor(completeness_group, levels = c(">80%","10-80%","<10%"),
                                                  ordered = T))
    
    clean.lifestyle.data.all.inverse <- clean.phastyle.lifestyle.data |>
        dplyr::inner_join(clean.phatyp.lifestyle.data, by = join_by(Accession == Accession), 
                          keep = F, relationship = "one-to-one", multiple = "any") |>
        dplyr::inner_join(clean.checkv.data, by = join_by(Accession == contig_id),
                          keep = F, relationship = "one-to-one", multiple = "any") |>
        dplyr::left_join(clean.blast.data, by = join_by(Accession == Accession), keep = F,
                         relationship = "one-to-one") |>
        dplyr::left_join(tmp.clean.data.relabund, by = join_by(Accession), keep = F,
                         relationship = "one-to-one", multiple = "any") |>
        dplyr::left_join(clean.genomad.data, by = join_by(Accession == contig_id),
                         keep = F, relationship = "one-to-one", multiple = "any") |>
        dplyr::mutate(tailed_phage = if_else(grepl("Caudoviricetes", taxonomy), TRUE, FALSE)) |>
        dplyr::mutate(blasthit = ifelse(!is.na(evalue), T, F)) |>
        dplyr::mutate(phatyp_lytic = ifelse(PhaTYPType == "virulent", TRUE, FALSE),
                      phatyp_temperate = ifelse(PhaTYPType == "temperate", TRUE, FALSE),
                      phastyle_lytic = ifelse(PhaStyleType == "virulent", TRUE, FALSE),
                      phastyle_temperate = ifelse(PhaStyleType == "temperate", TRUE, FALSE)) |>
        dplyr::filter((as.character(PhaStyleType) != as.character(PhaTYPType)) | 
                          (PhaStyleType == "temperate" & PhaTYPType == "virulent") |
                          (PhaStyleType == "virulent" & PhaTYPType == "temperate")) |>
        tidyr::drop_na(completeness) |>
        dplyr::mutate(completeness_group = ifelse(completeness > 80, ">80%",
                                                  ifelse(completeness >= 10, "10-80%", "<10%"))) |>
        dplyr::mutate(completeness_group = factor(completeness_group, levels = c(">80%","10-80%","<10%"),
                                                  ordered = T))
    
    # Select PhaStyle as we pick the overlap anyway
    clean.lifestyle.data <- clean.lifestyle.data.all |>
        dplyr::select(Accession, Length, PhaStyleScore, PhaStyleType) |>
        dplyr::rename(PhaTYPScore = PhaStyleScore, Type = PhaStyleType)
    
    clean.lifestyle.data.inverse <- clean.lifestyle.data.all.inverse |>
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
    tmp.instrain.data.list.below_thresh <- load_srr_libraries(srr_file_vector = test.input.files.instrain.vec,
                                                              lifestyle_df = rbind(clean.lifestyle.data,
                                                                                   clean.lifestyle.data.inverse),
                                                              checkv_df = clean.checkv.data,
                                                              thresholds = list(
                                                                  breadth = 0,
                                                                  coverage = 0,
                                                                  completeness = 0
                                                              ))
    tmp.instrain.data.list.inverse <- load_srr_libraries(srr_file_vector = test.input.files.instrain.vec,
                                                         lifestyle_df = clean.lifestyle.data.inverse,
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
        seco.clean.lifestyle.data.inverse <- clean.lifestyle.data.all.inverse
        main.clean.data.inverse <- tmp.instrain.data.list.inverse[[1]]
        main.clean.data.below_thresh <- tmp.instrain.data.list.below_thresh[[1]]
        
        first_iteration <- F
    } else {
        main.clean.data <- rbind(main.clean.data, tmp.clean.data)
        seco.clean.lifestyle.data <- rbind(seco.clean.lifestyle.data, clean.lifestyle.data.all)
        seco.clean.lifestyle.data.inverse <- rbind(seco.clean.lifestyle.data.inverse, clean.lifestyle.data.all.inverse)
        main.clean.data.inverse <- rbind(main.clean.data.inverse, tmp.instrain.data.list.inverse[[1]])
        main.clean.data.below_thresh <- rbind(main.clean.data.below_thresh, tmp.instrain.data.list.below_thresh[[1]])
        
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
                paste0("Ratio in main data: ", length(unique(sort(tmp.clean.data$genome_ID)))/length(unique(sort(tmp.clean.data.relabund$Accession)))),
                sep = "   "))
    rm("tmp.clean.data", "tmp.clean.stat.data", "tmp.list.clean.data", "tmp.raw.stat.data.all",
       "clean.lifestyle.data", "clean.lifestyle.data.all", "clean.lifestyle.data.all.inverse",
       "clean.lifestyle.data.inverse", "clean.checkv.data", "tmp.clean.data.relabund",
       "tmp.instrain.data.list",
       "tmp.instrain.data.list.inverse", 
       "clean.phastyle.lifestyle.data", "clean.phatyp.lifestyle.data")
}

# Conversion to one big dataframe
main.clean.data.plot <- main.clean.data |>
    dplyr::filter(variable == "nucl_diversity") |>
    dplyr::mutate(SRR_ID_lifestyle = paste(SRR_ID, lifestyle)) |>
    dplyr::group_by(SRR_ID)

#######################
# Statistics markup
main.clean.stat.data.new <- main.clean.stat.data |>
    dplyr::mutate(markup = ifelse(wilcox_pvalue > 0.05, "",
                                  ifelse(wilcox_pvalue > 0.01, "*",
                                         ifelse(wilcox_pvalue > 0.001, "**",
                                                ifelse(wilcox_pvalue > 0.0001, "***",
                                                       ifelse(wilcox_pvalue < 0.0001, "****", NA)))))) |>
    dplyr::ungroup() |>
    dplyr::mutate(markup_estimate = ifelse(!markup %in% c("",NA), round(wilcox_estimate, digits = 4), "")) |>
    dplyr::distinct(SRR_ID, .keep_all = T) |>
    dplyr::mutate(SRR_ID = as.character(SRR_ID)) |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID)

SRR_TO_PLOT <- unique(sort(main.clean.data.plot$SRR_ID))
SRR_TO_PLOT.data <- dplyr::filter(main.clean.data.plot, value > 0 &
                                      SRR_ID %in% SRR_TO_PLOT) |>
    dplyr::ungroup() |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID) |>
    #dplyr::mutate(study_ID = as.factor(study_ID)) |>
    dplyr::mutate(study_ID = fct_relevel(study_ID, c("Wildfire","WetUp","PermaThaw","NewMexico","NCalHabitat","MediGrass","Land-use",
                                                     "Intertidal","Hopland","Glacial","Conifer","BodegaBay","AgriSoil")))

##### Make interval plots, by lifestyle and by study ####
p.interval.microdiversity_bysample_subset <- SRR_TO_PLOT.data |>
    ggplot( aes(x = study_ID, y = value, group = lifestyle)) +
    stat_halfeye(data = dplyr::filter(SRR_TO_PLOT.data, lifestyle == "virulent"), fill = "#6b200c",
                 position = position_nudge(x = .08), alpha = 0.4, scale = 0.40) +
    stat_halfeye(data = dplyr::filter(SRR_TO_PLOT.data, lifestyle == "temperate"), side = "bottom",
                 fill = "#133e7e", position = position_nudge(x = -0.08), alpha = 0.4, scale = 0.40) +
    stat_interval(position = position_dodgejust(width = 0.4), linewidth = 2, width = 1, show.legend = NA,
                  colour = c(rep(head(MetBrewer::met.brewer("OKeeffe1", direction = -1), 3), length(SRR_TO_PLOT)),
                             rep(head(MetBrewer::met.brewer("OKeeffe1", direction = 1), 3), length(SRR_TO_PLOT)))) +
    stat_summary(geom = "point", fun = median, position = position_dodgejust(width = 0.4), colour = "white", size = 0.8, fill = NA) +
    geom_text(inherit.aes = F, data = dplyr::filter(main.clean.stat.data.new, SRR_ID %in% SRR_TO_PLOT),
              aes(x = study_ID, y = 0.0001, label = markup), hjust = 0.5, nudge_x = 0) +
    geom_text(inherit.aes = F, data = dplyr::filter(main.clean.stat.data.new, SRR_ID %in% SRR_TO_PLOT),
              aes(x = study_ID, y = 0.00012, label = markup_estimate), hjust = 0) +
    coord_flip(clip = "on", ylim = c(0.0001,0.035)) +
    scale_x_discrete(expand = expansion(c(0,0))) +
    scale_y_log10(guide = "axis_logticks") +
    theme_classic() +
    labs(title = "**B** Microdiversity per dataset, split by lifestyle",
         #subtitle = "P values adjusted using BH method. Top bar = virulent, bottom bar = temperate\n *: p<0.05, **: p<0.01, ***: p<0.001, ****: p<0.0001",
         y = "Microdiversity (\u03c0)",
         x = "Dataset ID") +
    theme(#plot.background = element_rect(color = NA, fill = bg_color),
        panel.grid.major.y = element_line(linewidth = 1, linetype = 2),
        axis.text.y = element_text(hjust = 1),
        plot.title = element_markdown(),
        plot.title.position = "plot",
        legend.position = "right")
p.interval.microdiversity_bysample_subset

p.interval.microdiversity_all <- SRR_TO_PLOT.data |>
    ggplot( aes(x = study_ID, y = value)) +
    stat_halfeye(data = dplyr::filter(SRR_TO_PLOT.data), fill = "black",
                 position = position_nudge(x = 0), alpha = 0.4, scale = 0.45) +
    stat_interval(position = position_dodgejust(width = 0.4), linewidth = 3, width = 1, show.legend = NA) +
    stat_summary(geom = "point", fun = median, colour = "black", size = 0.8, fill = NA) +
    coord_flip(clip = "on", ylim = c(0.0001,0.035)) +
    scale_x_discrete(expand = expansion(c(0.095,0.095))) +
    scale_y_log10(guide = "axis_logticks") +
    scale_color_viridis_d(breaks = c("0.5", "0.8", "0.95"),
                          labels = c("50%", "80%", "95%"),
                          name = "Proportion of\nmicrodiversity") +
    theme_classic() +
    labs(title = "**A** Microdiversity per dataset",
         #subtitle = "P values adjusted using BH method. Top bar = virulent, bottom bar = temperate\n *: p<0.05, **: p<0.01, ***: p<0.001, ****: p<0.0001",
         y = "Microdiversity (\u03c0)",
         x = "Dataset ID") +
    theme(#plot.background = element_rect(color = NA, fill = bg_color),
        panel.grid.major.y = element_line(linewidth = 1, linetype = 2),
        axis.text.y = element_text(hjust = 1),
        plot.title = element_markdown(),
        plot.title.position = "plot",
        legend.box.just = "center",
        legend.position = "right",
        legend.margin = margin())

tmp.interval.legend <- get_legend(p.interval.microdiversity_all)
tmp.interval.legend.2 <- get_legend(p.interval.microdiversity_all + scale_color_manual(values = c("0.5" = "#133e7e", "0.8" = "#225bb2", "0.95" = "#447fdd"),
                                                                                       breaks = c("0.5", "0.8", "0.95"),
                                                                                       labels = c("50%", "80%", "95%"),
                                                                                       name = "temperate"))
tmp.interval.legend.3 <- get_legend(p.interval.microdiversity_all + scale_color_manual(values = c("0.5" = "#6b200c", "0.8" = "#973d21", "0.95" = "#da6c42"),
                                                                                       breaks = c("0.5", "0.8", "0.95"),
                                                                                       labels = c("50%", "80%", "95%"),
                                                                                       name = "virulent") +
                                        theme(legend.box.just = "center"))
p.interval.microdiversity_all <- p.interval.microdiversity_all +
    theme(legend.position = "none")

p_collected_interval <- gridExtra::grid.arrange(p.interval.microdiversity_all,
                                                p.interval.microdiversity_bysample_subset,
                                                tmp.interval.legend,
                                                tmp.interval.legend.2,
                                                tmp.interval.legend.3,
                                                rectGrob(gp=gpar(col=NA)),
                                                rectGrob(gp=gpar(col=NA)),
                                                top = NULL,
                                                ncol = 4, nrow = 4,
                                                layout_matrix = cbind(c(1,1,1,1),c(2,2,2,2),c(7,3,4,6),c(7,3,5,6)),
                                                widths = c(2.4,2.4,.45,.45))

# ggsave(plot = p_collected_interval,
#        filename = paste(sys.args$data_wd, "lifestyle_microdiversity_plot_all.png", sep = "/"),
#        width = 3335, height = 1400,
#        units = "px", dpi = 600, scale = 2.0)
# 
# ggsave(plot = p_collected_interval,
#        filename = paste(sys.args$data_wd, "lifestyle_microdiversity_plot_all.tiff", sep = "/"),
#        width = 1700, height = 600,
#        units = "px", dpi = 300, scale = 1.8)

# Summary statistics
main.data.summary.stats <- main.clean.data |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::ungroup() |>
    dplyr::group_by(SRR_ID, lifestyle) |>
    dplyr::distinct(genome_ID, .keep_all = T) |>
    dplyr::summarise(count = n(),
                     avg_length = mean(length),
                     med_length = median(length),
                     avg_nucl_diversity = mean(nucl_diversity),
                     med_nucl_diversity = median(nucl_diversity),
                     avg_breadth = mean(breadth),
                     med_breadth = median(breadth),
                     avg_coverage = mean(coverage),
                     med_coverage = median(coverage),
                     quant95 = quantile(nucl_diversity, probs = 0.95),
                     .groups = "keep") |>
    dplyr::select(SRR_ID, lifestyle, count, med_length, med_nucl_diversity, med_breadth, med_coverage) |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID)

#####PROTOTYPING####

# Plot bar graphs with relative abundance correction
# p.bar.lifestyle_total.all.per_study.per_lifestyle.corr_rel_abundance <- main.clean.data |>
#     tidyr::pivot_wider(names_from = variable, values_from = value) |>
#     dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness),
#                      by = join_by(genome_ID == Accession),
#                      relationship = "one-to-one", keep = F) |>
#     dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
#                                      if_else(SRR_ID == "SRR000004", "AgriSoil",
#                                              if_else(SRR_ID == "SRR000002","BodegaBay",
#                                                      if_else(SRR_ID == "SRR000001", "Land-use",
#                                                              if_else(SRR_ID == "SRR000003", "Intertidal",
#                                                                      if_else(SRR_ID == "SRR000007", "Wildfire",
#                                                                              if_else(SRR_ID == "SRR000008", "Antarctic", NA))))))),
#                   .after = SRR_ID) |>
#     dplyr::filter(completeness >= 80) |>
#     dplyr::group_by(study_ID, lifestyle) |>
#     ggplot( aes(x = study_ID, fill = lifestyle)) +
#     geom_bar(data = tidyr::pivot_wider(main.clean.data, names_from = variable, values_from = value) |>
#                  dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness),
#                                   by = join_by(genome_ID == Accession),
#                                   relationship = "one-to-one", keep = F) |>
#                  dplyr::bind_rows(tmp.inverse.clean.data, tmp.below_thresh.clean.data) |>
#                  dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
#                                                   if_else(SRR_ID == "SRR000004", "AgriSoil",
#                                                           if_else(SRR_ID == "SRR000002","BodegaBay",
#                                                                   if_else(SRR_ID == "SRR000001", "Land-use",
#                                                                           if_else(SRR_ID == "SRR000003", "Intertidal",
#                                                                                   if_else(SRR_ID == "SRR000007", "Wildfire",
#                                                                                           if_else(SRR_ID == "SRR000008", "Antarctic", NA))))))),
#                                .after = SRR_ID) |>
#                  dplyr::mutate(used = ifelse(lifestyle %in% c("virulent","temperate"), "used","discarded")) |>
#                  dplyr::group_by(study_ID, used),
#              aes(x = study_ID, fill = used, weight = rel_abundance), position = "fill", just = 1, width = 0.3) +
#     geom_bar(position = "fill", alpha = 1, just = 0, width = 0.3, aes(weight = rel_abundance)) +
#     geom_text(data = summarise(group_by(ungroup(main.clean.stat.data), SRR_ID), sum = sum(count)) |>
#                   dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
#                                                    if_else(SRR_ID == "SRR000004", "AgriSoil",
#                                                            if_else(SRR_ID == "SRR000002","BodegaBay",
#                                                                    if_else(SRR_ID == "SRR000001", "Land-use",
#                                                                            if_else(SRR_ID == "SRR000003", "Intertidal",
#                                                                                    if_else(SRR_ID == "SRR000007", "Wildfire",
#                                                                                            if_else(SRR_ID == "SRR000008", "Antarctic", NA))))))),
#                                 .after = SRR_ID),
#               aes(x = study_ID, y = 1.09, label = sum), inherit.aes = F, vjust = 1, nudge_x = 0.18, angle = 30) +
#     geom_text(data = dplyr::bind_rows(tmp.inverse.clean.data, tmp.below_thresh.clean.data) |>
#                   dplyr::distinct(genome_ID, .keep_all = T) |>
#                   dplyr::filter(!genome_ID %in% SRR_TO_PLOT.data$genome_ID) |>
#                   dplyr::group_by(SRR_ID) |>
#                   dplyr::summarise(sum = n()) |>
#                   dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
#                                                    if_else(SRR_ID == "SRR000004", "AgriSoil",
#                                                            if_else(SRR_ID == "SRR000002","BodegaBay",
#                                                                    if_else(SRR_ID == "SRR000001", "Land-use",
#                                                                            if_else(SRR_ID == "SRR000003", "Intertidal",
#                                                                                    if_else(SRR_ID == "SRR000007", "Wildfire",
#                                                                                            if_else(SRR_ID == "SRR000008", "Antarctic", NA))))))),
#                                 .after = SRR_ID),
#               aes(x = study_ID, y = 1.09, label = sum), inherit.aes = F, vjust = 1, nudge_x = -0.18, angle = 30) +
#     scale_fill_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e",
#                                  "below_thresh" = "lightgrey", "unknown" = "grey",
#                                  "used" = "black", "discarded" = "grey"),
#                       breaks = c("temperate", "virulent", "below_thresh", "unknown", "used", "discarded"),
#                       labels = c("temperate", "virulent", "below quality\nthreshold",
#                                  "unknown lifestyle", "used", "discarded"),
#                       name = "vOTU lifestyle") +
#     scale_y_continuous(labels = scales::label_percent(), expand = expansion(c(0,0.1), 0), breaks = c(0, 0.25, 0.50, 0.75, 1)) +
#     labs(title = "Ratio of vOTU lifestyle",
#          subtitle = "corrected by relative abundance",
#          y = "Percentage of vOTUs",
#          x = "Dataset ID") +
#     theme_classic() +
#     theme(legend.position = "right",
#           legend.background = element_rect(fill = "white"),
#           panel.grid.minor.y = element_line(colour = "grey95"),
#           panel.grid.major.y = element_line(colour = "grey90"),
#           plot.title.position = "plot",
#           plot.background = element_rect(fill = "white"),
#           axis.title.x = element_blank())
# 
# ggsave(plot = p.bar.lifestyle_total.all.per_study.per_lifestyle.corr_rel_abundance,
#        filename = paste(sys.args$data_wd, "lifestyle_plot_ratios_rel_abundance.png", sep = "/"),
#        width = 700, height = 500,
#        units = "px", dpi = 400, scale = 2.5)

# Plot coverage (or relative abundance) x microdiversity relationship
library(ggpointdensity)
library(corrplot)

main.clean.data |>
    dplyr::filter(variable == "coverage") |>
    ggplot( aes(x=value, group=SRR_ID, color=SRR_ID)) +
    geom_density() +
    labs(x = "coverage") +
    scale_x_log10()

main.clean.data |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness), 
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-one", keep = F) |>
    dplyr::group_by(lifestyle) |>
    dplyr::filter(completeness >= 80) |> 
    dplyr::slice_sample(n = 2000) |>
    ggplot( aes(x = rel_abundance, y = nucl_diversity, z = lifestyle)) +
    geom_pointdensity(adjust = 1) +
    scale_color_viridis_c() +
    scale_y_log10() +
    scale_x_log10() +
    facet_wrap(~lifestyle, nrow = 2, ncol = 1)

tmp.corr.matrix <- main.clean.data |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness), 
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-one", keep = F) |>
    dplyr::select(4:9)
testRes <- cor.mtest(tmp.corr.matrix, conf.level = 0.95, method = "sp", exact = F)
# png(filename = paste(sys.args$data_wd, "corr_plot_small_all_datasets.png", sep = "/"),
#        width = 1200, height = 1200,
#        units = "px", pointsize = 40)
# corrplot(cor(tmp.corr.matrix, use = "pairwise.complete.obs", method = "sp"), p.mat = testRes$p, 
#          method = 'color', diag = FALSE, type = 'upper', addCoef.col ='grey50', 
#          sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, col = rev(COL2("RdBu")),
#          insig = 'blank', pch.col = 'grey20', order = 'AOE')
# dev.off()

main.clean.data |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness), 
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-one", keep = F) |>
    dplyr::filter(completeness >= 80 & SRR_ID == "SRR000001") |>
    summary()

# Plot numbers / bar charts for overview plots
tmp.inverse.clean.data <- main.clean.data.inverse |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(dplyr::select(seco.clean.lifestyle.data.inverse, Accession, rel_abundance, completeness), 
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-one", keep = F) |>
    dplyr::mutate(lifestyle = "unknown")

tmp.below_thresh.clean.data <- main.clean.data.below_thresh |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::filter(!genome_ID %in% main.clean.data$genome_ID) |>
    dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness), 
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-one", keep = F) |>
    dplyr::mutate(lifestyle = "below_thresh")

##### Combination plot of general vOTU statistics ####
# Percentage fill of ratio lifestyle per study 
p.bar.lifestyle_total.all.per_study.per_lifestyle <- main.clean.data |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness), 
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-one", keep = F) |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID) |>
    dplyr::filter(completeness >= 80) |>
    dplyr::group_by(study_ID, lifestyle) |>
    ggplot( aes(x = study_ID, fill = lifestyle)) +
    geom_bar(data = tidyr::pivot_wider(main.clean.data, names_from = variable, values_from = value) |>
                 dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness), 
                                  by = join_by(genome_ID == Accession),
                                  relationship = "one-to-one", keep = F) |>
                 dplyr::bind_rows(dplyr::bind_rows(tmp.inverse.clean.data, tmp.below_thresh.clean.data) |>
                                      #dplyr::distinct(genome_ID, .keep_all = T) |>
                                      dplyr::filter(!genome_ID %in% main.clean.data$genome_ID)) |>
                 dplyr::mutate(study_ID = SRR2name(SRR_ID),
                               .after = SRR_ID) |>
                 dplyr::mutate(used = ifelse(lifestyle %in% c("virulent","temperate"), "used","discarded")) |>
                 dplyr::group_by(study_ID, used), 
             aes(x = study_ID, fill = used), position = "fill", just = 1, width = 0.3) +
    geom_bar(position = "fill", alpha = 1, just = 0, width = 0.3) +
    geom_text(data = summarise(group_by(ungroup(main.clean.stat.data), SRR_ID), sum = sum(count)) |>
                  dplyr::mutate(study_ID = SRR2name(SRR_ID),
                                .after = SRR_ID),
              aes(x = study_ID, y = 1.09, label = sum), inherit.aes = F, vjust = 1, nudge_x = 0.1, angle = 25) +
    # geom_text(data = dplyr::bind_rows(tmp.inverse.clean.data, tmp.below_thresh.clean.data) |>
    #               dplyr::distinct(genome_ID, .keep_all = T) |>
    #               dplyr::filter(!genome_ID %in% SRR_TO_PLOT.data$genome_ID) |>
    #               dplyr::group_by(SRR_ID) |>
    #               dplyr::summarise(sum = n()) |>
    #               dplyr::mutate(study_ID = SRR2name(SRR_ID),
    #                             .after = SRR_ID),
    #           aes(x = study_ID, y = 1.11, label = sum), colour = "grey40", inherit.aes = F, vjust = 1, nudge_x = -0.18, angle = 30) +
    scale_fill_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                 "below_thresh" = "lightgrey", "unknown" = "grey",
                                 "used" = "black", "discarded" = "grey"),
                      breaks = c("temperate", "virulent", "below_thresh", "unknown", "used", "discarded"),
                      labels = c("temperate", "virulent", "below quality\nthreshold", 
                                 "unknown lifestyle", "used", "discarded"),
                      name = "vOTU lifestyle") +
    scale_y_continuous(labels = scales::label_percent(), expand = expansion(c(0,0.1), 0), breaks = c(0, 0.25, 0.50, 0.75, 1)) +
    labs(title = "**A** Ratio of vOTU lifestyle",
         y = "Percentage of vOTUs",
         x = "Dataset ID") +
    theme_classic() +
    theme(legend.position = "right",
          legend.background = element_rect(fill = "white"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.title = element_markdown(),
          plot.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
tmp.lifestyle.legend <- get_legend(p.bar.lifestyle_total.all.per_study.per_lifestyle)
p.bar.lifestyle_total.all.per_study.per_lifestyle <- p.bar.lifestyle_total.all.per_study.per_lifestyle +
    theme(legend.position = "none")

# ggsave(plot = p.bar.lifestyle_total.all.per_study.per_lifestyle,
#        filename = paste(sys.args$data_wd, "lifestyle_plot_ratios.png", sep = "/"),
#        width = 1000, height = 750,
#        units = "px", dpi = 400, scale = 2)

# Boxplot of length per lifestyle per study 
p.box.length.all.per_study.per_lifestyle <- main.clean.data |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness), 
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-one", keep = F) |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID) |>
    dplyr::filter(completeness >= 80) |>
    dplyr::group_by(study_ID, lifestyle) |>
    ggplot( aes(x = study_ID, y = length, fill = lifestyle)) +
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
    #scale_y_continuous(transform = "log10", n.breaks = 8, minor_breaks = c(15000, 40000, 75000, 150000, 400000)) +
    scale_y_continuous(expand = expansion(c(0, 0.01), 0), limits = c(0, 700000), labels = c("0","200kbp","400kbp","600kbp"),
                       breaks = c(0,200000,400000,600000)) +
    labs(title = "**B** vOTU length",
         y = "Length of vOTU",
         x = "Dataset ID") +
    theme_classic() +
    theme(legend.position = "none",
          legend.background = element_rect(fill = "grey93"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title = element_markdown(),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

# Boxplot of coverage per lifestyle per study 
p.box.coverage.all.per_study.per_lifestyle <- main.clean.data |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness), 
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-one", keep = F) |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID) |>
    dplyr::filter(completeness >= 80) |>
    dplyr::group_by(study_ID, lifestyle) |>
    ggplot( aes(x = study_ID, y = coverage, fill = lifestyle)) +
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
    #scale_y_log10(guide = "axis_logticks", breaks = c(10, 50, 100, 500, 1000)) +
    labs(title = "**D** vOTU coverage",
         y = "Coverage of vOTU",
         x = "Dataset ID") +
    #geom_hline(yintercept = 10, color = "red") +
    theme_classic() +
    theme(legend.position = "none",
          legend.background = element_rect(fill = "grey93"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title = element_markdown(),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"),
          axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

# ggsave(plot = p.box.coverage.all.per_study.per_lifestyle,
#        filename = paste(sys.args$data_wd, "lifestyle_plot_coverage.png", sep = "/"),
#        width = 1000, height = 750,
#        units = "px", dpi = 400, scale = 2)

# Boxplot of completeness per lifestyle per study 
p.box.completeness.all.per_study.per_lifestyle <- main.clean.data |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness), 
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-one", keep = F) |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID) |>
    dplyr::filter(completeness >= 80) |>
    dplyr::group_by(study_ID, lifestyle) |>
    ggplot( aes(x = study_ID, y = completeness, fill = lifestyle)) +
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
    labs(title = "**C** vOTU completeness",
         y = "Completeness of vOTU",
         x = "Dataset ID") +
    scale_y_continuous(labels = scales::label_percent(scale = 1), limits = c(0, 100),
                       expand = expansion(0, c(0,2))) +
    geom_hline(yintercept = 80, color = "red") +
    theme_classic() +
    theme(legend.position = "none",
          legend.background = element_rect(fill = "grey93"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title = element_markdown(),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"),
          axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))

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
#        filename = paste(sys.args$data_wd, "lifestyle_plot_all.png", sep = "/"),
#        width = 3000, height = 2000,
#        units = "px", dpi = 400, scale = 1.5)
# 
# ggsave(plot = p_collected_population,
#        filename = paste(sys.args$data_wd, "lifestyle_plot_all.tiff", sep = "/"),
#        width = 3000, height = 2000,
#        units = "px", dpi = 400, scale = 1.5)

# Bar chart with lifestyles and unknown and below thresh data
main.clean.data |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness), 
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-one", keep = F) |>
    dplyr::bind_rows(tmp.inverse.clean.data) |>
    dplyr::bind_rows(tmp.below_thresh.clean.data) |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID) |>
    dplyr::group_by(study_ID, lifestyle) |>
    ggplot( aes(x = study_ID, 
                fill = factor(lifestyle, levels = c("temperate", "virulent", 
                                                    "below_thresh", "unknown"), ordered = T))) +
    #geom_bar(position = "fill", aes(weight = rel_abundance), width = 0.3) +
    geom_bar(position = "dodge", width = 0.5) +
    scale_fill_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                 "below_thresh" = "lightgrey", "unknown" = "grey"),
                      breaks = c("temperate", "virulent", "below_thresh", "unknown"),
                      labels = c("temperate", "virulent", "below quality threshold", "unknown"),,
                      name = "vOTU lifestyle") +
    labs(title = "Total number of vOTUs per dataset, split by lifestyle",
         y = "Number of vOTUs",
         x = "Dataset ID") +
    theme_classic() +
    scale_y_continuous(transform = "log10", limits = c(1, 50000), n.breaks = 6, minor_breaks = c(5, 50, 500, 5000, 50000)) +
    theme(legend.position = "inside",
          legend.position.inside = c(0.2, 0.90),
          legend.background = element_rect(fill = "grey93"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot")

stop("Main script done, stopping...")
####Wilcox test for length####
tmp.clean.data <- main.clean.data %>%
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    # dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness), 
    #                  by = join_by(genome_ID == Accession),
    #                  relationship = "one-to-one", keep = F) |>
    #dplyr::mutate(nucl_diversity_cor = nucl_diversity * rel_abundance) |>
    #dplyr::filter(variable == "length") |>
    #dplyr::group_by(SRR_ID, lifestyle) %>%
    #dplyr::select(-variable) %>%
    # dplyr::filter(SRR_ID %in% as.character({dplyr::summarise(., count = n(), .groups = "keep") |>
    #         dplyr::group_by(SRR_ID) |>
    #         dplyr::filter(any(lifestyle == "virulent") & any(lifestyle == "temperate")) |>
    #         dplyr::filter(all(count > 10))}$SRR_ID)) |>
    #dplyr::ungroup() |>
    base::split(~SRR_ID, drop = T) |>
    purrr::map(~ wilcox.test(nucl_diversity ~ lifestyle, data = ., conf.int = T)) |>
    purrr::map_dfr( ~ broom::tidy(.), .id = "SRR_ID") |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID) |>
    dplyr::mutate(wilcox_pvalue = p.value, wilcox_estimate = estimate) |>
    dplyr::mutate(markup = ifelse(wilcox_pvalue > 0.05, "",
                                  ifelse(wilcox_pvalue > 0.01, "*",
                                         ifelse(wilcox_pvalue > 0.001, "**",
                                                ifelse(wilcox_pvalue > 0.0001, "***",
                                                       ifelse(wilcox_pvalue < 0.0001, "****", NA)))))) |>
    dplyr::mutate(markup_estimate = ifelse(!markup %in% c("",NA), wilcox_estimate, ""))

tmp.clean.data.corr_rel_abundance <- main.clean.data %>%
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness), 
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-one", keep = F) |>
    dplyr::mutate(nucl_diversity_cor = nucl_diversity * rel_abundance)

tmp.clean.data.corr_rel_abundance |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID) |>
    ggplot( aes(x = study_ID, y = nucl_diversity, fill = lifestyle)) +
    geom_point(position = position_jitterdodge(), aes(colour = lifestyle)) +
    geom_boxplot(outliers = F, alpha = 0.8) +
    scale_y_log10()

main.clean.stat.data |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID) |>
    View()

###### Upset plot for tool mismatch visualisation ####
p.upset.lifestyle.mismatch <- seco.clean.lifestyle.data |>
    dplyr::union(seco.clean.lifestyle.data.inverse) |>
    dplyr::filter(Accession %in% c(main.clean.data$genome_ID, main.clean.data.inverse$genome_ID)) |>
    dplyr::filter(completeness > 80) |>
    ComplexUpset::upset(intersect = c("phatyp_lytic","phatyp_temperate",
                                      "phastyle_lytic","phastyle_temperate"),
                        name = element_blank(),
                        base_annotations=list(
                            'Intersection size'=intersection_size(
                                counts=FALSE,
                                mapping=aes(fill=blasthit))
                            +theme(axis.title.y = element_text(angle = 70, vjust = 0.5),
                                   legend.position = "top", plot.title.position = "plot")
                            +scale_fill_discrete(name = "BLAST hit w/ integrase")
                            +scale_y_continuous(position = "left")
                            +labs(title = "Comparison between lifestyle prediction tools \n(PhaTYP and PhaStyle)")
                        ),
                        annotations = list(),
                        set_sizes=(
                            upset_set_size(
                                geom=geom_bar(aes(fill = blasthit), width=0.4))
                            + geom_label(aes(label=after_stat(count)), size = 2.5, nudge_y = 5000,
                                         hjust = 0.5, stat='count', label.size = unit(0, "lines"))
                            + theme(axis.text.x=element_text(angle=45, vjust = 0.5),
                                    axis.ticks.x=element_line(),
                                    axis.title.x = element_blank(),
                                    legend.position = "none")
                            #+ ylim(20000, 0)
                        ),
                        sort_sets=F) & theme(plot.background = element_rect(fill = "white", colour = NA))

ggsave(plot = p.upset.lifestyle.mismatch,
       filename = paste(sys.args$data_wd, "lifestyle_upset_plot_mismatch.png", sep = "/"),
       width = 3000, height = 2000,
       units = "px", dpi = 400, scale = 1.0)

##### Bar plots using taxonomy data ####
ggplotly(main.clean.data |>
             tidyr::pivot_wider(names_from = variable, values_from = value) |>
             dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness, taxonomy), 
                              by = join_by(genome_ID == Accession),
                              relationship = "one-to-one", keep = F) |> 
             dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                              if_else(SRR_ID == "SRR000004", "AgriSoil",
                                                      if_else(SRR_ID == "SRR000002","BodegaBay",
                                                              if_else(SRR_ID == "SRR000001", "Land-use",
                                                                      if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                              if_else(SRR_ID == "SRR000007", "Wildfire", 
                                                                                      if_else(SRR_ID == "SRR000008", "Antarctic", 
                                                                                              if_else(SRR_ID == "SRR000009", "Hopland", NA)))))))),
                           .after = SRR_ID) |>
             dplyr::group_by(study_ID, lifestyle) |>
             dplyr::mutate(tailed_phage = if_else(grepl("Caudoviricetes", taxonomy), TRUE, FALSE)) |>
             dplyr::mutate(taxonomy = as.factor(taxonomy)) |>
             ggplot( aes(x = study_ID, fill = taxonomy)) +
             geom_bar(stat = "count") +
             facet_wrap(~tailed_phage, scales = "free_y")
)

ggplotly(seco.clean.lifestyle.data |>
             bind_rows(seco.clean.lifestyle.data.inverse) |>
             dplyr::mutate(tailed_phage = if_else(grepl("Caudoviricetes", taxonomy), TRUE, FALSE)) |>
             dplyr::mutate(taxonomy = as.factor(taxonomy)) |>
             ggplot( aes(x = taxonomy, fill = taxonomy)) +
             geom_bar(stat = "count") +
             facet_wrap(~tailed_phage, scales = "free_y") +
             theme(axis.text.x = element_blank())
)

SRR_TO_PLOT.data.wide <- main.clean.data.plot |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(seco.clean.lifestyle.data,
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-many", keep = F, suffix = c("_genome", "_gene")) |>
    dplyr::mutate(tailed_phage = if_else(grepl("Caudoviricetes", taxonomy), TRUE, FALSE)) |>
    dplyr::filter(tailed_phage == TRUE) |>
    dplyr::mutate(taxonomy = as.factor(taxonomy))

main.clean.stat.data.new <- SRR_TO_PLOT.data.wide |>
    base::split(~SRR_ID, drop = T) |>
    purrr::map(~ wilcox.test(nucl_diversity ~ lifestyle, data = ., conf.int = T)) |>
    purrr::map_dfr( ~ broom::tidy(.), .id = "SRR_ID") |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID) |>
    dplyr::rename(wilcox_pvalue = p.value, wilcox_estimate = estimate) |>
    dplyr::mutate(markup = ifelse(wilcox_pvalue > 0.05, "",
                                  ifelse(wilcox_pvalue > 0.01, "*",
                                         ifelse(wilcox_pvalue > 0.001, "**",
                                                ifelse(wilcox_pvalue > 0.0001, "***",
                                                       ifelse(wilcox_pvalue < 0.0001, "****", NA)))))) |>
    dplyr::ungroup() |>
    dplyr::mutate(markup_estimate = ifelse(!markup %in% c("",NA), round(wilcox_estimate, digits = 4), "")) |>
    dplyr::distinct(SRR_ID, .keep_all = T) |>
    dplyr::mutate(SRR_ID = as.character(SRR_ID))

SRR_TO_PLOT <- unique(sort(main.clean.data.plot$SRR_ID))
SRR_TO_PLOT.data <- dplyr::filter(main.clean.data.plot, value > 0 &
                                      SRR_ID %in% SRR_TO_PLOT & genome_ID %in% SRR_TO_PLOT.data.wide$genome_ID)

p.interval.microdiversity_bysample_subset_tailed <- SRR_TO_PLOT.data |>
    ggplot( aes(x = SRR_ID, y = value, group = lifestyle)) +
    stat_halfeye(data = dplyr::filter(SRR_TO_PLOT.data, lifestyle == "virulent"), fill = "#6b200c",
                 position = position_nudge(x = .06), alpha = 0.4, scale = 0.5) +
    stat_halfeye(data = dplyr::filter(SRR_TO_PLOT.data, lifestyle == "temperate"), side = "bottom",
                 fill = "#133e7e", position = position_nudge(x = -0.06), alpha = 0.4, scale = 0.5) +
    stat_interval(position = position_dodgejust(width = 0.3), linewidth = 3, width = 1, show.legend = NA,
                  colour = c(rep(head(MetBrewer::met.brewer("OKeeffe1", direction = -1), 3), length(SRR_TO_PLOT)),
                             rep(head(MetBrewer::met.brewer("OKeeffe1", direction = 1), 3), length(SRR_TO_PLOT)))) +
    stat_summary(geom = "point", fun = median, position = position_dodgejust(width = 0.3), colour = "white") +
    geom_text(inherit.aes = F, data = dplyr::filter(main.clean.stat.data.new, SRR_ID %in% SRR_TO_PLOT),
              aes(x = SRR_ID, y = 0.0001, label = markup), hjust = 0.5, nudge_x = 0) +
    geom_text(inherit.aes = F, data = dplyr::filter(main.clean.stat.data.new, SRR_ID %in% SRR_TO_PLOT),
              aes(x = SRR_ID, y = 0.00012, label = markup_estimate), hjust = 0) +
    coord_flip(clip = "on", ylim = c(0.0001,0.02)) +
    scale_x_discrete(breaks = waiver(), labels = c("SRR000006" = "Glacial",
                                                   "SRR000004" = "AgriSoil",
                                                   "SRR000002" = "BodegaBay",
                                                   "SRR000001" = "Land-use",
                                                   "SRR000003" = "Intertidal",
                                                   "SRR000007" = "Wildfire",
                                                   "SRR000008" = "Antarctic"),
                     expand = expansion(c(0,0))) +
    scale_y_log10(guide = "axis_logticks") +
    theme_classic() +
    labs(title = "(b). Microdiversity levels per dataset, tailed only",
         #subtitle = "P values adjusted using BH method. Top bar = virulent, bottom bar = temperate\n *: p<0.05, **: p<0.01, ***: p<0.001, ****: p<0.0001",
         y = "Microdiversity (log10(\u03c0))",
         x = "Dataset ID") +
    theme(#plot.background = element_rect(color = NA, fill = bg_color),
        panel.grid.major.y = element_line(linewidth = 1, linetype = 2),
        axis.text.y = element_text(hjust = 1),
        plot.title.position = "plot",
        legend.position = "right")

p_collected_interval_tailed <- gridExtra::grid.arrange(p.interval.microdiversity_bysample_subset,
                                                       p.interval.microdiversity_bysample_subset_tailed,
                                                       tmp.interval.legend,
                                                       tmp.interval.legend.2,
                                                       tmp.interval.legend.3,
                                                       rectGrob(gp=gpar(col=NA)),
                                                       top = NULL,
                                                       ncol = 4, nrow = 3,
                                                       layout_matrix = cbind(c(1,1,1),c(2,2,2),c(3,4,6),c(3,5,6)),
                                                       widths = c(2.5,2.5,.3,.3))

ggsave(plot = p_collected_interval_tailed,
       filename = paste(sys.args$data_wd, "lifestyle_microdiversity_plot_tailed_vs_all.png", sep = "/"),
       width = 2300, height = 750,
       units = "px", dpi = 400, scale = 2)

p.bar.lifestyle_total.all.per_study.per_lifestyle.tailed <- main.clean.data |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness, taxonomy), 
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-one", keep = F) |>
    dplyr::mutate(tailed_phage = if_else(grepl("Caudoviricetes", taxonomy), TRUE, FALSE)) |>
    dplyr::filter(tailed_phage == TRUE) |>
    dplyr::mutate(taxonomy = as.factor(taxonomy)) |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID) |>
    dplyr::group_by(study_ID, lifestyle) |>
    ggplot( aes(x = study_ID, fill = lifestyle)) +
    geom_bar(data = tidyr::pivot_wider(main.clean.data, names_from = variable, values_from = value) |>
                 dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness, taxonomy), 
                                  by = join_by(genome_ID == Accession),
                                  relationship = "one-to-one", keep = F) |>
                 dplyr::mutate(tailed_phage = if_else(grepl("Caudoviricetes", taxonomy), TRUE, FALSE)) |>
                 dplyr::mutate(taxonomy = as.factor(taxonomy)) |>
                 dplyr::mutate(study_ID = SRR2name(SRR_ID),
                               .after = SRR_ID) |>
                 dplyr::group_by(study_ID, lifestyle) |>
                 dplyr::mutate(tailed = ifelse(tailed_phage, "tailed","non-tailed")) |>
                 dplyr::group_by(study_ID, tailed), 
             aes(x = study_ID, fill = tailed), position = "fill", just = 1, width = 0.3) +
    geom_bar(position = "fill", alpha = 1, just = 0, width = 0.3) +
    geom_text(data = ungroup(SRR_TO_PLOT.data.wide) |>
                  group_by(SRR_ID) |>
                  summarise(sum = n()) |>
                  dplyr::mutate(study_ID = SRR2name(SRR_ID),
                                .after = SRR_ID),
              aes(x = study_ID, y = 1.09, label = sum), inherit.aes = F, vjust = 1, nudge_x = 0.18, angle = 30) +
    geom_text(data = ungroup(main.clean.stat.data) |>
                  group_by(SRR_ID) |>
                  summarise(sum = sum(count)) |>
                  dplyr::mutate(study_ID = SRR2name(SRR_ID),
                                .after = SRR_ID),
              aes(x = study_ID, y = 1.09, label = sum), inherit.aes = F, vjust = 1, nudge_x = -0.18, angle = 30) +
    scale_fill_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                 "below_thresh" = "lightgrey", "unknown" = "grey",
                                 "tailed" = "black", "non-tailed" = "grey"),
                      breaks = c("temperate", "virulent", "below_thresh", "unknown", "used", "discarded","tailed","non-tailed"),
                      labels = c("temperate", "virulent", "below quality\nthreshold", 
                                 "unknown lifestyle", "used", "discarded","tailed","non-tailed"),
                      name = "vOTU lifestyle") +
    scale_y_continuous(labels = scales::label_percent(), expand = expansion(c(0,0.1), 0), breaks = c(0, 0.25, 0.50, 0.75, 1)) +
    labs(title = "Ratio of vOTU lifestyle and tailed phages",
         y = "Percentage of vOTUs",
         x = "Dataset ID") +
    theme_classic() +
    theme(legend.position = "right",
          legend.background = element_rect(fill = "white"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"),
          axis.title.x = element_blank())

ggsave(plot = p.bar.lifestyle_total.all.per_study.per_lifestyle.tailed,
       filename = paste(sys.args$data_wd, "lifestyle_plot_ratios_tailed_vs_nontailed.png", sep = "/"),
       width = 1000, height = 750,
       units = "px", dpi = 400, scale = 2)

p.bar.lifestyle_total.all.per_study.per_used.tailed <- main.clean.data |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness, taxonomy), 
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-one", keep = F) |>
    dplyr::mutate(tailed_phage = if_else(grepl("Caudoviricetes", taxonomy), TRUE, FALSE)) |>
    dplyr::filter(tailed_phage == TRUE) |>
    dplyr::mutate(taxonomy = as.factor(taxonomy)) |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID) |>
    dplyr::group_by(study_ID, lifestyle) |>
    ggplot( aes(x = study_ID, fill = lifestyle)) +
    geom_bar(data = tidyr::pivot_wider(main.clean.data, names_from = variable, values_from = value) |>
                 dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness, taxonomy), 
                                  by = join_by(genome_ID == Accession),
                                  relationship = "one-to-one", keep = F) |>
                 dplyr::bind_rows(tmp.inverse.clean.data, tmp.below_thresh.clean.data) |>
                 dplyr::mutate(tailed_phage = if_else(grepl("Caudoviricetes", taxonomy), TRUE, FALSE)) |>
                 dplyr::mutate(taxonomy = as.factor(taxonomy)) |>
                 dplyr::mutate(study_ID = SRR2name(SRR_ID),
                               .after = SRR_ID) |>
                 dplyr::mutate(used = ifelse(lifestyle %in% c("virulent","temperate") & tailed_phage, "used, tail","discarded")) |>
                 dplyr::group_by(study_ID, used), 
             aes(x = study_ID, fill = used), position = "fill", just = -0.1, width = 0.3) +
    geom_bar(data = tidyr::pivot_wider(main.clean.data, names_from = variable, values_from = value) |>
                 dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness, taxonomy), 
                                  by = join_by(genome_ID == Accession),
                                  relationship = "one-to-one", keep = F) |>
                 dplyr::bind_rows(tmp.inverse.clean.data, tmp.below_thresh.clean.data) |>
                 dplyr::mutate(tailed_phage = if_else(grepl("Caudoviricetes", taxonomy), TRUE, FALSE)) |>
                 dplyr::mutate(taxonomy = as.factor(taxonomy)) |>
                 dplyr::mutate(study_ID = SRR2name(SRR_ID),
                               .after = SRR_ID) |>
                 dplyr::mutate(used = ifelse(lifestyle %in% c("virulent","temperate"), "used","discarded")) |>
                 dplyr::group_by(study_ID, used), 
             aes(x = study_ID, fill = used), position = "fill", just = 1.1, width = 0.3) +
    scale_fill_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                 "below_thresh" = "lightgrey", "unknown" = "grey",
                                 "used" = "black", "discarded" = "grey", "used, tail" = "#2480e6"),
                      breaks = c("temperate", "virulent", "below_thresh", "unknown", "discarded", "used", "tailed","non-tailed","used, tail"),
                      labels = c("temperate", "virulent", "below quality\nthreshold", 
                                 "unknown lifestyle", "discarded", "used", "tailed","non-tailed","used, tail"),
                      name = "vOTU lifestyle") +
    scale_y_continuous(labels = scales::label_percent(), expand = expansion(c(0,0.1), 0), breaks = c(0, 0.25, 0.50, 0.75, 1)) +
    labs(title = "Ratio of vOTU lifestyle and tailed phages",
         y = "Percentage of vOTUs",
         x = "Dataset ID") +
    theme_classic() +
    theme(legend.position = "right",
          legend.background = element_rect(fill = "white"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"),
          axis.title.x = element_blank())

SRR_TO_PLOT.data.wide |>
    dplyr::group_by(SRR_ID, tailed_phage) |>
    dplyr::summarise(count = n()) |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID) |>
    dplyr::ungroup() |>
    dplyr::select(-SRR_ID) |>
    knitr::kable()

#### Bootstrapping/randomizing vOTU selection ####

# Setup data structures
bootstrap.clean.data <- data.frame(matrix(ncol = 16, nrow = 0))
bootstrap.clean.data.all <- data.frame(matrix(ncol = 8, nrow = 0))
bootstrap_iterations <- 50

# Go through SRR_IDs
for (bootstrap_i in 1:bootstrap_iterations){
    for (SRR_i in 1:length(SRR_TO_PLOT)){
        skip_to_next <<- F
        #bootstrap_i <- 1
        #SRR_i <- 1
        
        # Set names and select df
        tmp.SRR_ID <- SRR_TO_PLOT[SRR_i]
        tmp.SRR_ID.clean.data <- main.clean.data |>
            dplyr::filter(SRR_ID == tmp.SRR_ID) |>
            dplyr::mutate(bootstrap_iteration = bootstrap_i, .after = SRR_ID)
        
        # Create sub-dfs
        tmp.SRR_ID.clean.data.virulent <- dplyr::filter(tmp.SRR_ID.clean.data, lifestyle == "virulent") |>
            tidyr::pivot_wider(names_from = variable, values_from = value)
        tmp.SRR_ID.clean.data.temperate <- dplyr::filter(tmp.SRR_ID.clean.data, lifestyle == "temperate") |>
            tidyr::pivot_wider(names_from = variable, values_from = value)
        
        # Set sample size
        if (nrow(tmp.SRR_ID.clean.data.virulent) < nrow(tmp.SRR_ID.clean.data.temperate)){
            sample_size_length <- nrow(tmp.SRR_ID.clean.data.virulent)
        } else {
            sample_size_length <- nrow(tmp.SRR_ID.clean.data.temperate)
        }
        sample_size <- as.integer(sample_size_length * 0.5)
        
        # Sample x out of each
        tryCatch({
            tmp.SRR_ID.clean.data.virulent <- dplyr::slice_sample(tmp.SRR_ID.clean.data.virulent, n = sample_size)
            tmp.SRR_ID.clean.data.temperate <- dplyr::slice_sample(tmp.SRR_ID.clean.data.temperate, n = sample_size)
        }, error = function(e) {skip_to_next <<- TRUE})
        if(skip_to_next){
            rm(skip_to_next)
            next
        }
        
        # Test 
        tmp.wilcox.test <- wilcox.test(tmp.SRR_ID.clean.data.virulent$nucl_diversity, 
                                       tmp.SRR_ID.clean.data.temperate$nucl_diversity, conf.int = T)
        
        # Calculate summary stats on lifestyle data
        tmp.summary.stats <- tmp.SRR_ID.clean.data.virulent %>%
            rbind(tmp.SRR_ID.clean.data.temperate) %>%
            dplyr::group_by(SRR_ID, lifestyle) %>%
            dplyr::summarise(count = n(),
                             avg_length = mean(length),
                             med_length = median(length),
                             avg_nucl_diversity = mean(nucl_diversity),
                             med_nucl_diversity = median(nucl_diversity),
                             avg_breadth = mean(breadth),
                             med_breadth = median(breadth),
                             avg_coverage = mean(coverage),
                             med_coverage = median(coverage),
                             quant95 = quantile(nucl_diversity, probs = 0.95),
                             .groups = "keep") %>%
            dplyr::mutate(bootstrap_iteration = bootstrap_i, .after = SRR_ID) %>%
            dplyr::mutate(wilcox_pvalue = tmp.wilcox.test$p.value) %>%
            dplyr::mutate(wilcox_estimate = round(tmp.wilcox.test$estimate, digits = 7)) %>%
            dplyr::mutate(lifestyle_ratio = length(tmp.SRR_ID.clean.data.virulent$nucl_diversity)/length(tmp.SRR_ID.clean.data.temperate), .after = count)
        
        # Save to big df
        bootstrap.clean.data <- rbind(bootstrap.clean.data, tmp.summary.stats)
        bootstrap.clean.data.all <- bootstrap.clean.data.all %>%
            rbind(tmp.SRR_ID.clean.data.temperate) %>%
            rbind(tmp.SRR_ID.clean.data.virulent)
    }
}

# setup df for plotting
colnames(bootstrap.clean.data) <- colnames(tmp.summary.stats)

#bootstrap.clean.data$p.value.adj <- p.adjust(bootstrap.clean.data$wilcox_pvalue, method = "BH")
bootstrap.clean.data.adj <- bootstrap.clean.data %>%
    dplyr::group_by(SRR_ID) %>%
    dplyr::mutate(p.value.adj = p.adjust(wilcox_pvalue, method='BH')) |>
    #dplyr::mutate(p.value.adj = wilcox_pvalue) |>
    dplyr::mutate(markup = ifelse(p.value.adj > 0.05, "",
                                  ifelse(p.value.adj > 0.01, "*",
                                         ifelse(p.value.adj > 0.001, "**",
                                                ifelse(p.value.adj > 0.0001, "***",
                                                       ifelse(p.value.adj < 0.0001, "****", NA)))))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(estimate = round(wilcox_estimate, digits = 7)) %>%
    dplyr::mutate(markup_estimate = ifelse(!markup %in% c("","NS",NA), estimate, "")) %>%
    #dplyr::distinct(SRR_ID, .keep_all = T) %>%
    dplyr::mutate(SRR_ID = as.character(SRR_ID))

# Plot
p.heatmap.microdiversity_bysample.bootstrapped <- bootstrap.clean.data.adj |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID) |>
    ggplot( aes(y = study_ID, x = bootstrap_iteration, fill = estimate)) +
    geom_tile() +
    geom_text(inherit.aes = F, 
              aes(y = study_ID, x = bootstrap_iteration, label = markup), color = "black", hjust = 0.5, nudge_x = 0, angle = 30) +
    scale_fill_gradientn(colours = MetBrewer::met.brewer("OKeeffe1", type = "continuous", direction = -1), limits = c(- 0.005, 0.005), na.value = "grey80") +
    theme_classic() +
    labs(title = "Estimated difference in microdiversity between samples (bootstrapped)",
         #subtitle = "P values adjusted using BH method. \n*: p<0.05, **: p<0.01, ***: p<0.001, ****: p<0.0001",
         x = "Bootstrap iteration",
         y = "Dataset ID",
         fill = "Estimated difference \nbetween lifestyle",
         caption = paste(paste0("Data = ", "all"),
                         paste0("Bootstraps = ", bootstrap_iterations), 
                         "Adjusted = Yes (BH)",
                         sep = "\n")) +
    scale_x_continuous(expand = expansion(c(0,0))) +
    scale_y_discrete(expand = expansion(c(0,0))) +
    theme(legend.position = "right",
          plot.caption.position = "plot",
          plot.title.position = "plot")
p.heatmap.microdiversity_bysample.bootstrapped

ggsave(plot = p.heatmap.microdiversity_bysample.bootstrapped,
       filename = paste(sys.args$data_wd, "microdiversity_bootstrap_plot.png", sep = "/"),
       width = 2100, height = 1500,
       units = "px", dpi = 400, scale = 1.5)

##### Statistical prototypes ####
tmp.data <- main.clean.data |> 
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(seco.clean.lifestyle.data,
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-many", keep = F, suffix = c("_genome", "_gene")) |>
    dplyr::mutate(tailed_phage = if_else(grepl("Caudoviricetes", taxonomy), TRUE, FALSE)) |>
    dplyr::filter(tailed_phage == TRUE) |>
    dplyr::mutate(taxonomy = as.factor(taxonomy))

tmp.data |>
    base::split(~SRR_ID, drop = T) |>
    purrr::map(~ cor(.$nucl_diversity, .$coverage, method = "sp")) |>
    purrr::map_dfr( ~ broom::tidy(.), .id = "SRR_ID") |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID) |>
    rbind(tmp.data |>
              base::split(~SRR_ID, drop = T) |>
              purrr::map(~ cor(data = ., formula = nucl_diversity ~ breadth, method = "sp")) |>
              purrr::map_dfr( ~ broom::tidy(.), .id = "SRR_ID") |>
              dplyr::mutate(study_ID = SRR2name(SRR_ID),
                            .after = SRR_ID) |>
              dplyr::filter(term != "(Intercept)")) |>
    rbind(tmp.data |>
              base::split(~SRR_ID, drop = T) |>
              purrr::map(~ lm(data = ., formula = nucl_diversity ~ breadth)) |>
              purrr::map_dfr( ~ broom::tidy(.), .id = "SRR_ID") |>
              dplyr::mutate(study_ID = SRR2name(SRR_ID),
                            .after = SRR_ID) |>
              dplyr::filter(term != "(Intercept)"))
rm(tmp.data)

SRR_TO_PLOT.data.wide <- main.clean.data |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(seco.clean.lifestyle.data,
                     by = join_by(genome_ID == Accession),
                     relationship = "many-to-one", keep = F, multiple = "any") |>
    dplyr::left_join(main.clean.data |>
                         dplyr::ungroup() |>
                         tidyr::pivot_wider(names_from = variable, values_from = value) |>
                         dplyr::group_by(genome_ID) |>
                         dplyr::summarise(occurences = n()),
                     by = join_by(genome_ID))

SRR_TO_PLOT.data.wide |>
    dplyr::group_by(SRR_ID, lifestyle, genome_ID) |>
    base::split(~SRR_ID, drop = T) |>
    purrr::map(~ wilcox.test(length ~ lifestyle, data = ., conf.int = T)) |>
    purrr::map_dfr( ~ broom::tidy(.), .id = "SRR_ID") |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID)

SRR_TO_PLOT.data.wide |>
    dplyr::group_by(SRR_ID, lifestyle, genome_ID) |>
    base::split(~SRR_ID, drop = T) |>
    purrr::map(~ wilcox.test(completeness ~ lifestyle, data = ., conf.int = T)) |>
    purrr::map_dfr( ~ broom::tidy(.), .id = "SRR_ID") |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID)

SRR_TO_PLOT.data.wide |>
    dplyr::group_by(SRR_ID, lifestyle, genome_ID) |>
    base::split(~SRR_ID, drop = T) |>
    purrr::map(~ wilcox.test(breadth ~ lifestyle, data = ., conf.int = T)) |>
    purrr::map_dfr( ~ broom::tidy(.), .id = "SRR_ID") |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID)

SRR_TO_PLOT.data.wide |>
    dplyr::group_by(SRR_ID, lifestyle, genome_ID) |>
    base::split(~SRR_ID, drop = T) |>
    purrr::map(~ wilcox.test(coverage ~ lifestyle, data = ., conf.int = T)) |>
    purrr::map_dfr( ~ broom::tidy(.), .id = "SRR_ID") |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID)

SRR_TO_PLOT.data.wide.inverse <- main.clean.data.inverse |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(seco.clean.lifestyle.data.inverse,
                     by = join_by(genome_ID == Accession),
                     relationship = "many-to-one", keep = F, multiple = "any") |>
    dplyr::left_join(main.clean.data |>
                         dplyr::ungroup() |>
                         tidyr::pivot_wider(names_from = variable, values_from = value) |>
                         dplyr::group_by(genome_ID) |>
                         dplyr::summarise(occurences = n()),
                     by = join_by(genome_ID))

SRR_TO_PLOT.data.wide.below_thresh <- main.clean.data.below_thresh |>
    dplyr::ungroup() |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(rbind(seco.clean.lifestyle.data.inverse, seco.clean.lifestyle.data),
                     by = join_by(genome_ID == Accession),
                     relationship = "many-to-one", keep = F, multiple = "any") |>
    dplyr::filter(SRR_ID %in% c("SRR000006","SRR000004","SRR000002","SRR000001")) |>
    # dplyr::filter(PhaStyleType %in% c("virulent","temperate") & PhaTYPType %in% c("virulent","temperate")) |>
    # dplyr::filter((phatyp_temperate == T & phastyle_temperate == T) | (phastyle_lytic == T & phatyp_lytic == T)) |>
    # dplyr::filter((phatyp_temperate == T & blasthit == T) | (phatyp_lytic == T & blasthit == F)) |>
    dplyr::filter(!is.na(taxonomy)) |>
    dplyr::filter(tailed_phage == TRUE) |>
    dplyr::filter(coverage >= 10 & completeness >= 80 & breadth >= 0.8) |>
    group_by(SRR_ID) |>
    summarise(count = n(), avg_cov = mean(coverage), avg_length = mean(length)) |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID) |>
    dplyr::mutate(percent = (count / 30762) * 100)

main.clean.data |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, evalue, completeness, taxonomy, tailed_phage), 
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-one", keep = F) |>
    dplyr::group_by(SRR_ID) |>
    summarise(count = n(), avg_length = mean(length), avg_cov = mean(coverage), 
              avg_breadth = mean(breadth), avg_completeness = mean(completeness)) |>
    dplyr::mutate(study_ID = SRR2name(SRR_ID),
                  .after = SRR_ID)
    View()
    
    
