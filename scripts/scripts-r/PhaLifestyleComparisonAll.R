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
dataset.vec <- c("bodega_all", "glacial_all",
                 "santosviromes_all", "landuse_all")

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
    bioproject_ID <- sys.args$bioproject_id
    
    #####
    ### Data loading
    input.files.instrain.list <- list.files(path = paste(data_wd, input.file.instrain.dir, sep = "/"),
                                            full.names = F)
    input.files.SRR.vector <- sapply(unlist(input.files.instrain.list), SRR_extract)
    test.input.files.instrain.vec <- sapply(input.files.instrain.list, function(x){
        paste(data_wd, input.file.instrain.dir, x[1], "output", x[1], sep = "/")
    })
    
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
        dplyr::slice(-1) |>
        dplyr::mutate(genome_ID = ifelse(grepl("single_", genome_ID), as.character(gsub("single_","", genome_ID)), NA)) |>
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
    
    clean.lifestyle.data.all.inverse <- clean.phastyle.lifestyle.data |>
        dplyr::inner_join(clean.phatyp.lifestyle.data, by = join_by(Accession == Accession), 
                          keep = F, relationship = "one-to-one", multiple = "any") |>
        dplyr::inner_join(clean.checkv.data, by = join_by(Accession == contig_id),
                          keep = F, relationship = "one-to-one", multiple = "any") |>
        dplyr::left_join(clean.blast.data, by = join_by(Accession == Accession), keep = F,
                         relationship = "one-to-one") |>
        dplyr::left_join(tmp.clean.data.relabund, by = join_by(Accession), keep = F,
                         relationship = "one-to-one", multiple = "any") |>
        dplyr::mutate(blasthit = ifelse(!is.na(evalue), T, F)) |>
        dplyr::mutate(phatyp_lytic = ifelse(PhaTYPType == "virulent", TRUE, FALSE),
                      phatyp_temperate = ifelse(PhaTYPType == "temperate", TRUE, FALSE),
                      phastyle_lytic = ifelse(PhaStyleType == "virulent", TRUE, FALSE),
                      phastyle_temperate = ifelse(PhaStyleType == "temperate", TRUE, FALSE)) |>
        dplyr::filter((as.character(PhaStyleType) != as.character(PhaTYPType)) | 
                          (PhaStyleType == "temperate" & PhaTYPType == "temperate" & blasthit == F) |
                          (PhaStyleType == "virulent" & PhaTYPType == "virulent" & blasthit == T)) |>
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
                                                              lifestyle_df = clean.lifestyle.data,
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
                sep = "   "))
    rm("tmp.clean.data", "tmp.clean.stat.data", "tmp.list.clean.data", "tmp.raw.stat.data.all",
       "clean.lifestyle.data", "clean.lifestyle.data.all", "clean.lifestyle.data.all.inverse",
       "clean.lifestyle.data.inverse", "clean.checkv.data", "tmp.clean.data.relabund",
       "tmp.instrain.data.list",
       "tmp.instrain.data.list.inverse", "tmp.inverse.clean.data", 
       "clean.phastyle.lifestyle.data", "clean.phatyp.lifestyle.data")
}

# Conversion to one big dataframe
main.clean.data.plot <- main.clean.data |>
    dplyr::filter(variable == "nucl_diversity") |>
    dplyr::mutate(SRR_ID_lifestyle = paste(SRR_ID, lifestyle)) |>
    dplyr::group_by(SRR_ID)

#######################
### World map #########
library(giscoR)
library(sf)

tmp.raw.data.locations <- read.csv(paste0(sys.args$data_wd, "/", "SRR_all_locations.tsv"),
                                   sep = "\t", header = T)

tmp.clean.data.locations <- tmp.raw.data.locations |>
    tidyr::separate_wider_delim(Lat_Long, delim = " N ", 
                                names = c("lat", "long")) |>
    dplyr::mutate(lat = as.numeric(lat)) |>
    dplyr::mutate(long = ifelse(grepl("E", long), as.double(gsub(" E","",long)),
                                ifelse(grepl("W", long), (as.double(gsub(" W","",long)) * -1),
                                       NA))) |>
    dplyr::filter(Dataset_ID != "Wildfire") |>
    st_as_sf(coords = c("long", "lat"), crs = "EPSG:4326") |>
    dplyr::mutate(SRR_ID = if_else(Dataset_ID == "Glacial", "SRR000006",
                                   if_else(Dataset_ID == "AgriSoil", "SRR000004",
                                           if_else(Dataset_ID == "BodegaBay","SRR000002",
                                                   if_else(Dataset_ID == "Landuse", "SRR000001",
                                                           if_else(Dataset_ID == "Intertidal", "SRR000003", NA))))),
                  .after = Dataset_ID)

tmp.clean.data.locations.distinct <- tmp.clean.data.locations |>
    dplyr::distinct(Dataset_ID, .keep_all = T) |>
    dplyr::select(-Run_SRR_ID) |>
    dplyr::left_join(main.clean.stat.data, by = join_by(SRR_ID)) |>
    dplyr::ungroup() |>
    dplyr::group_by(Dataset_ID)

p.bar.lifestyle.percent.vOTU.all <- tmp.clean.data.locations.distinct |>
    dplyr::filter(!Dataset_ID %in% c("Watershed","Intertidal")) |>
    ggplot( aes(x = factor(Dataset_ID, levels=c("Glacial", "Intertidal", "Landuse","BodegaBay","AgriSoil","Watershed")), 
                y = count, fill = lifestyle)) +
    geom_bar(position="fill", stat="identity") +
    facet_wrap(~ factor(Dataset_ID, levels=c("Glacial", "Intertidal", "Landuse","BodegaBay","AgriSoil","Watershed")), 
               scales = "free_x", nrow = 1) +
    labs(y = "Ratio of lifestyle") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "bottom") +
    scale_y_continuous(labels = scales::percent)

tmp.world.map <- giscoR::gisco_get_countries(year = "2024", epsg = "4326", resolution = "60")

p.map.vOTU_position.all <- ggplot() +
    geom_sf(data = st_wrap_dateline(st_transform(tmp.world.map, "+proj=lonlat +lon_0=160")), fill = "grey", alpha = 0.3) +
    geom_sf(data = summarize(tmp.clean.data.locations.distinct, count = sum(count)), 
            aes(geometry = geometry, color = Dataset_ID, size = count), alpha = 0.7) +
    coord_sf(ylim = c(-20,75), xlim = c(-80, 120)) +
    scale_size_continuous(range = c(4, 15), name = "Number of vOTUs\nin dataset") +
    scale_color_discrete(name = "Dataset") +
    theme_classic() +
    theme(axis.line.y = element_blank(),
          legend.position = "top")
p.map.vOTU_position.all

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
    dplyr::mutate(SRR_ID = as.character(SRR_ID))

SRR_TO_PLOT <- unique(sort(main.clean.data.plot$SRR_ID))
SRR_TO_PLOT.data <- dplyr::filter(main.clean.data.plot, value > 0 &
                                      SRR_ID %in% SRR_TO_PLOT)

##### Make interval plots, by lifestyle and by study ####
p.interval.microdiversity_bysample_subset <- SRR_TO_PLOT.data |>
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
                                                   "SRR000007" = "Wildfire"),
                     expand = expansion(c(0,0))) +
    scale_y_log10(guide = "axis_logticks") +
    theme_classic() +
    labs(title = "(b). Microdiversity levels per dataset, split by lifestyle",
         #subtitle = "P values adjusted using BH method. Top bar = virulent, bottom bar = temperate\n *: p<0.05, **: p<0.01, ***: p<0.001, ****: p<0.0001",
         y = "Microdiversity (log10(\u03c0))",
         x = "Dataset ID") +
    theme(#plot.background = element_rect(color = NA, fill = bg_color),
          panel.grid.major.y = element_line(linewidth = 1, linetype = 2),
          axis.text.y = element_text(hjust = 1),
          plot.title.position = "plot",
          legend.position = "right")
p.interval.microdiversity_bysample_subset

p.interval.microdiversity_all <- SRR_TO_PLOT.data |>
    ggplot( aes(x = SRR_ID, y = value)) +
    stat_halfeye(data = dplyr::filter(SRR_TO_PLOT.data), fill = "black",
                 position = position_nudge(x = 0), alpha = 0.4, scale = 0.5) +
    stat_interval(position = position_dodgejust(width = 0.4), linewidth = 3, width = 1, show.legend = NA) +
    stat_summary(geom = "point", fun = median, colour = "black") +
    coord_flip(clip = "on", ylim = c(0.0001,0.02)) +
    scale_x_discrete(breaks = waiver(), labels = c("SRR000006" = "Glacial",
                                                   "SRR000004" = "AgriSoil",
                                                   "SRR000002" = "BodegaBay",
                                                   "SRR000001" = "Land-use",
                                                   "SRR000003" = "Intertidal",
                                                   "SRR000007" = "Wildfire"),
                     expand = expansion(c(0.4,0.05))) +
    scale_y_log10(guide = "axis_logticks") +
    scale_color_viridis_d(breaks = c("0.5", "0.8", "0.95"),
                          labels = c("50%", "80%", "95%"),
                          name = "Proportion of\nmicrodiversity") +
    theme_classic() +
    labs(title = "(a). Microdiversity levels per dataset",
         #subtitle = "P values adjusted using BH method. Top bar = virulent, bottom bar = temperate\n *: p<0.05, **: p<0.01, ***: p<0.001, ****: p<0.0001",
         y = "Microdiversity (log10(\u03c0))",
         x = "Dataset ID") +
    theme(#plot.background = element_rect(color = NA, fill = bg_color),
          panel.grid.major.y = element_line(linewidth = 1, linetype = 2),
          axis.text.y = element_text(hjust = 1),
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
                                                                                       name = "/virulent") +
                                        theme(legend.box.just = "center"))
p.interval.microdiversity_all <- p.interval.microdiversity_all +
    theme(legend.position = "none")

p_collected_interval <- gridExtra::grid.arrange(p.interval.microdiversity_all,
                                                p.interval.microdiversity_bysample_subset,
                                                tmp.interval.legend,
                                                tmp.interval.legend.2,
                                                tmp.interval.legend.3,
                                                rectGrob(gp=gpar(col=NA)),
                                                  top = NULL,
                                                  ncol = 4, nrow = 3,
                                                  layout_matrix = cbind(c(1,1,1),c(2,2,2),c(3,4,6),c(3,5,6)),
                                                  widths = c(2.5,2.5,.3,.3))

ggsave(plot = p_collected_interval,
       filename = paste(sys.args$data_wd, "lifestyle_microdiversity_plot_all.png", sep = "/"),
       width = 2300, height = 750,
       units = "px", dpi = 400, scale = 2)

# ggsave(plot = p.interval.microdiversity_bysample_subset, 
#        filename = paste(sys.args$data_wd, "microdiversity_plot_all.png", sep = "/"),
#        width = 3000, height = 2000,
#        units = "px", dpi = 400)

#####

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
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID)

#####PROTOTYPING####
####################
# Plot coverage (or relative abundance) x microdiversity relationship
library(ggpointdensity)
main.clean.data |>
    dplyr::filter(variable == "coverage") |>
    ggplot( aes(x=value, group=SRR_ID, color=SRR_ID)) +
    geom_density() +
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
                     relationship = "one-to-one", keep = F)
stats::cor(tmp.corr.matrix$rel_abundance, y= tmp.corr.matrix$nucl_diversity,
           method = "spearman")

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
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID) |>
    dplyr::filter(completeness >= 80) |>
    dplyr::group_by(study_ID, lifestyle) |>
    ggplot( aes(x = study_ID, fill = lifestyle)) +
    geom_bar(data = tidyr::pivot_wider(main.clean.data, names_from = variable, values_from = value) |>
                 dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness), 
                                  by = join_by(genome_ID == Accession),
                                  relationship = "one-to-one", keep = F) |>
                 dplyr::bind_rows(tmp.inverse.clean.data, tmp.below_thresh.clean.data) |>
                 dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                                  if_else(SRR_ID == "SRR000004", "AgriSoil",
                                                          if_else(SRR_ID == "SRR000002","BodegaBay",
                                                                  if_else(SRR_ID == "SRR000001", "Land-use",
                                                                          if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                                  if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                               .after = SRR_ID) |>
                 dplyr::mutate(used = ifelse(lifestyle == "virulent", "used","discarded")) |>
                 dplyr::group_by(study_ID, used), 
             aes(x = study_ID, fill = used), position = "fill", just = 1, width = 0.3) +
    geom_bar(position = "fill", alpha = 1, just = 0, width = 0.3) +
    geom_text(data = summarise(group_by(ungroup(main.clean.stat.data), SRR_ID), sum = sum(count)) |>
                  dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                                   if_else(SRR_ID == "SRR000004", "AgriSoil",
                                                           if_else(SRR_ID == "SRR000002","BodegaBay",
                                                                   if_else(SRR_ID == "SRR000001", "Land-use",
                                                                           if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                                   if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                                .after = SRR_ID),
              aes(x = study_ID, y = 1.09, label = sum), inherit.aes = F, vjust = 1, nudge_x = 0.18, angle = 30) +
    geom_text(data = tidyr::pivot_wider(main.clean.data, names_from = variable, values_from = value) |>
                  dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness), 
                                   by = join_by(genome_ID == Accession),
                                   relationship = "one-to-one", keep = F) |>
                  dplyr::bind_rows(tmp.inverse.clean.data, tmp.below_thresh.clean.data) |>
                  dplyr::group_by(SRR_ID) |>
                  dplyr::summarise(sum = n()) |>
                  dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                                   if_else(SRR_ID == "SRR000004", "AgriSoil",
                                                           if_else(SRR_ID == "SRR000002","BodegaBay",
                                                                   if_else(SRR_ID == "SRR000001", "Land-use",
                                                                           if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                                   if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                                .after = SRR_ID),
              aes(x = study_ID, y = 1.09, label = sum), inherit.aes = F, vjust = 1, nudge_x = -0.18, angle = 30) +
    scale_fill_manual(values = c("temperate" = "#2480e6", "virulent" = "#ec612e", 
                                 "below_thresh" = "lightgrey", "unknown" = "grey",
                                 "used" = "black", "discarded" = "grey"),
                      breaks = c("temperate", "virulent", "below_thresh", "unknown", "used", "discarded"),
                      labels = c("temperate", "virulent", "below quality\nthreshold", 
                                 "unknown lifestyle", "used", "discarded"),
                      name = "vOTU lifestyle") +
    scale_y_continuous(labels = scales::label_percent(), expand = expansion(c(0,0.1), 0), breaks = c(0, 0.25, 0.50, 0.75, 1)) +
    labs(title = "(a). Ratio of vOTU lifestyle",
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
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
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
    labs(title = "(b). vOTU length",
         y = "Length of vOTU",
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
p.box.coverage.all.per_study.per_lifestyle <- main.clean.data |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness), 
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-one", keep = F) |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
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
    labs(title = "(d). vOTU coverage",
         y = "Coverage of vOTU (log10)",
         x = "Dataset ID") +
    #geom_hline(yintercept = 10, color = "red") +
    theme_classic() +
    theme(legend.position = "none",
          legend.background = element_rect(fill = "grey93"),
          panel.grid.minor.y = element_line(colour = "grey95"),
          panel.grid.major.y = element_line(colour = "grey90"),
          plot.title.position = "plot",
          plot.background = element_rect(fill = "white"))

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
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
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
    labs(title = "(c). vOTU completeness",
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

ggsave(plot = p_collected_population,
       filename = paste(sys.args$data_wd, "lifestyle_plot_all.png", sep = "/"),
       width = 2000, height = 1500,
       units = "px", dpi = 400, scale = 1.5)

# Bar chart with lifestyles and unknown and below thresh data
main.clean.data |>
    tidyr::pivot_wider(names_from = variable, values_from = value) |>
    dplyr::left_join(dplyr::select(seco.clean.lifestyle.data, Accession, rel_abundance, completeness), 
                     by = join_by(genome_ID == Accession),
                     relationship = "one-to-one", keep = F) |>
    dplyr::bind_rows(tmp.inverse.clean.data) |>
    dplyr::bind_rows(tmp.below_thresh.clean.data) |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", NA))))),
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

stop()
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
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
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
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID) |>
    ggplot( aes(x = study_ID, y = nucl_diversity, fill = lifestyle)) +
    geom_point(position = position_jitterdodge(), aes(colour = lifestyle)) +
    geom_boxplot(outliers = F, alpha = 0.8) +
    scale_y_log10()

main.clean.stat.data |>
    dplyr::mutate(study_ID = if_else(SRR_ID == "SRR000006", "Glacial",
                                     if_else(SRR_ID == "SRR000004", "AgriSoil",
                                             if_else(SRR_ID == "SRR000002","BodegaBay",
                                                     if_else(SRR_ID == "SRR000001", "Land-use",
                                                             if_else(SRR_ID == "SRR000003", "Intertidal", 
                                                                     if_else(SRR_ID == "SRR000007", "Wildfire", NA)))))),
                  .after = SRR_ID) |>
    View()

###### Upset plot for tool mismatch visualisation ####
seco.clean.lifestyle.data |>
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
                            + geom_label(aes(label=after_stat(count)), size = 2.5, nudge_y = 2000,
                                         hjust = 0.5, stat='count', label.size = unit(0, "lines"))
                            + theme(axis.text.x=element_text(angle=45, vjust = 0.5),
                                    axis.ticks.x=element_line(),
                                    axis.title.x = element_blank(),
                                    legend.position = "none")
                            #+ ylim(20000, 0)
                        ),
                        sort_sets=F) & theme(plot.background = element_rect(fill = "grey97", colour = NA))








