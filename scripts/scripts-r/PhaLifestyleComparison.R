# Analyse lifestyle predictions
##############
### setup
rm(list = setdiff(ls(), c("master_script.data", "commandArgs_custom")))

dataset_id <- "glacial"
commandArgs_custom <- list(
    instrain_dir = paste("instrain", dataset_id, sep = "_"),
    lifestyle_phatyp_file = paste(dataset_id, "phatyp_prediction.tsv", sep = "_"),
    lifestyle_phastyle_file = paste(dataset_id, "phastyle_prediction.tsv", sep = "_"),
    checkv_file = paste(dataset_id, "quality_summary.tsv", sep = "_"),
    wd = "/Users/thomasdebruijn/Documents/PhD/R_PhD",
    data_wd = "/Users/thomasdebruijn/Documents/PhD/DATASETS",
    install = F,
    bioproject_id = "PRJNA882246"
)

# Load command arguments
sys.args <- commandArgs_custom

# Set global options
set.seed(1)
options(scipen = 6)
bg_color <- "grey97"

# Install all necessary packages
install <- sys.args$install
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
    install.packages('UpSetR')
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
    library(plotly)
    library(flexdashboard)
    library(crosstalk)
    library(UpSetR)
    library(VennDiagram)
    library(ComplexUpset)
})

# Load custom functions
source("/Users/thomasdebruijn/Documents/PhD/R_PhD/srr_library_load_func.R")

# Set files and params to load
wd <- sys.args$wd
setwd(wd)
data_wd <- sys.args$data_wd
input.file.instrain.dir <- sys.args$instrain_dir
input.file.phatyp.lifestyle.name <- sys.args$lifestyle_phatyp_file
input.file.phastyle.lifestyle.name <- sys.args$lifestyle_phastyle_file
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

clean.checkv.data <- raw.checkv.data %>%
    dplyr::mutate(contig_id = str_split_i(contig_id, "_length_", i = 1))
rm("raw.checkv.data")

# Load DIAMOND blastp integrase results and clean
raw.blast.data <- read.csv(file = paste(data_wd, paste0(dataset_id, "_against_integrases.tsv"), sep = "/"),
                           sep = "\t", header = F)

clean.blast.data <- raw.blast.data %>%
    purrr::set_names(c("qseqid","sseqid","qcovhsp","pident","length","mismatch",
                       "gapopen","qstart","qend","sstart","send","evalue","bitscore")) %>%
    dplyr::filter(evalue < 1e-20 & qcovhsp > 80) %>%
    dplyr::arrange(length, evalue) %>%
    dplyr::distinct(qseqid, .keep_all = T) %>%
    dplyr::mutate(qseqid = as.character(qseqid)) %>%
    dplyr::mutate(Accession = str_remove(qseqid, pattern = "_\\d{1,5}$"), .before = 1) %>%
    dplyr::mutate(Accession = str_split_i(Accession, "_length_", i = 1)) %>%
    dplyr::distinct(Accession, .keep_all = T) %>%
    dplyr::select(Accession, evalue)
rm("raw.blast.data")

# Load Phatyp predictions and clean
raw.phatyp.lifestyle.data <- read.csv(file = paste(data_wd, input.file.phatyp.lifestyle.name, sep = "/"), 
                                      sep = "\t", header = T)

clean.phatyp.lifestyle.data <- raw.phatyp.lifestyle.data %>%
    dplyr::mutate(PhaTYPType = as.factor(TYPE), .keep = "unused") %>%
    dplyr::mutate(Accession = str_split_i(Accession, "_length_", i = 1)) %>%
    dplyr::distinct(Accession, .keep_all = T)
rm("raw.phatyp.lifestyle.data")
clean.phatyp.virulent.set <- dplyr::filter(clean.phatyp.lifestyle.data, PhaTYPType == "virulent")$Accession
clean.phatyp.temperate.set <- dplyr::filter(clean.phatyp.lifestyle.data, PhaTYPType == "temperate")$Accession

# Load PhaStyle predictions and clean
raw.phastyle.lifestyle.data <- read.csv(file = paste(data_wd, input.file.phastyle.lifestyle.name, sep = "/"), 
                                        sep = "\t", header = T)

clean.phastyle.lifestyle.data <- raw.phastyle.lifestyle.data %>%
    dplyr::mutate(PhaStyleType = as.factor(predicted_label), .keep = "unused") %>%
    dplyr::rename(Accession = fasta_id) %>%
    dplyr::mutate(Accession = str_split_i(Accession, "_length_", i = 1)) %>%
    dplyr::distinct(Accession, .keep_all = T) %>%
    dplyr::mutate(PhaStyleScore = ifelse(PhaStyleType == "temperate", p_temperate,
                                         ifelse(PhaStyleType == "virulent", p_virulent, NA))) %>%
    dplyr::select(Accession, PhaStyleType, PhaStyleScore)
rm("raw.phastyle.lifestyle.data")
clean.phastyle.virulent.set <- dplyr::filter(clean.phastyle.lifestyle.data, PhaStyleType == "virulent")$Accession
clean.phastyle.temperate.set <- dplyr::filter(clean.phastyle.lifestyle.data, PhaStyleType == "temperate")$Accession

# Combine lifestyle predictions
clean.lifestyle.data.all <- clean.phatyp.lifestyle.data %>%
    dplyr::inner_join(clean.phastyle.lifestyle.data, by = join_by(Accession == Accession), 
                      keep = F, relationship = "one-to-one", multiple = "any") %>%
    dplyr::inner_join(clean.checkv.data, by = join_by(Accession == contig_id),
                      keep = F, relationship = "one-to-one", multiple = "any") %>%
    # dplyr::left_join(clean.blast.data, by = join_by(Accession == Accession), keep = F,
    #                   relationship = "one-to-one") %>%
    dplyr::mutate(phatyp_lytic = ifelse(PhaTYPType == "virulent", TRUE, FALSE),
                  phatyp_temperate = ifelse(PhaTYPType == "temperate", TRUE, FALSE),
                  phastyle_lytic = ifelse(PhaStyleType == "virulent", TRUE, FALSE),
                  phastyle_temperate = ifelse(PhaStyleType == "temperate", TRUE, FALSE)) %>%
    tidyr::drop_na(completeness) %>%
    dplyr::mutate(completeness_group = ifelse(completeness > 80, ">80%",
                                              ifelse(completeness >= 10, "10-80%", "<10%"))) %>%
    dplyr::mutate(completeness_group = factor(completeness_group, levels = c(">80%","10-80%","<10%"),
                                              ordered = T)) #%>%
#tidyr::drop_na(evalue)

# Visualize comparison of lifestyle predictions
# venn.data <- venn.diagram(x = list("phatyp_lytic" = clean.phatyp.virulent.set, 
#                                        "phatyp_temperate" = clean.phatyp.temperate.set,
#                                        "phastyle_lytic" = clean.phastyle.virulent.set, 
#                                        "phastyle_temperate" = clean.phastyle.temperate.set),
#                               filename = NULL)
# UpSetR::upset(fromList(list("phatyp_lytic" = clean.phatyp.virulent.set, 
#                             "phatyp_temperate" = clean.phatyp.temperate.set,
#                             "phastyle_lytic" = clean.phastyle.virulent.set, 
#                             "phastyle_temperate" = clean.phastyle.temperate.set)),
#               order.by = "freq")

ComplexUpset::upset(clean.lifestyle.data.all, intersect = c("phatyp_lytic","phatyp_temperate",
                                                            "phastyle_lytic","phastyle_temperate"),
                    name = element_blank(),
                    base_annotations=list(
                        'Intersection size'=intersection_size(
                            counts=FALSE,
                            mapping=aes(fill=completeness_group))
                        +theme(#axis.title.y = element_text(angle = 45, vjust = 0.5),
                            legend.position = "none")
                        +scale_y_continuous(position = "left")
                    ),
                    annotations = list(
                        "Completeness proportion"=list(
                            aes=aes(x=intersection, fill = completeness_group),
                            geom=list(
                                geom_bar(stat = "count", position = "fill"),
                                theme(#axis.title.y = element_text(angle = 45, vjust = 0.5),
                                    legend.position = "top"),
                                scale_y_continuous(labels=scales::percent_format(),
                                                   position = "left"),
                                labs(fill = "Estimated completeness (CheckV)") 
                            )
                        )
                    ),
                    set_sizes=(
                        upset_set_size(
                            geom=geom_bar(width=0.4))
                        + geom_label(aes(label=after_stat(count)), size = 3, hjust=-1, stat='count')
                        + theme(axis.text.x=element_text(angle=45, vjust = 0.5),
                                axis.ticks.x=element_line(),
                                axis.title.x = element_blank())
                    ),
                    sort_sets=F)

p.point.completeness_score_byPhaStyletype <- clean.lifestyle.data.all %>%
    ggplot( aes(x = completeness, y = PhaStyleScore, color = PhaStyleType)) +
    geom_point(alpha = 0.7) +
    labs(title = "Certainty of PhaStyle vs vOTU completeness")
p.point.completeness_score_byPhaStyletype.mp <- ggMarginal(p.point.completeness_score_byPhaStyletype, 
                                                           type = "violin", groupColour = T)

p.point.completeness_score_byPhaTYPtype <- clean.lifestyle.data.all %>%
    dplyr::filter(PhaTYPScore >= 0.5) %>%
    ggplot( aes(x = completeness, y = PhaTYPScore, color = PhaTYPType)) +
    geom_point(alpha = 0.7) +
    labs(title = "Certainty of PhaTYP vs vOTU completeness")
p.point.completeness_score_byPhaTYPtype.mp <- ggMarginal(p.point.completeness_score_byPhaTYPtype, 
                                                         type = "violin", groupColour = T)

ggMarginal(p.point.completeness_score_byPhaTYPtype / p.point.completeness_score_byPhaStyletype)

clean.lifestyle.data <- clean.lifestyle.data.all %>%
    dplyr::select(Accession, Length, PhaStyleScore, PhaStyleType) %>%
    dplyr::rename(PhaTYPScore = PhaStyleScore, Type = PhaStyleType)

# clean.lifestyle.data <- clean.lifestyle.data.all %>%
#     dplyr::select(Accession, Length, PhaTYPScore, PhaTYPType) %>%
#     dplyr::rename(Type = PhaTYPType)

#####

test.output.data <- load_srr_libraries(srr_file_vector = test.input.files.instrain.vec,
                                       lifestyle_df = clean.lifestyle.data,
                                       checkv_df = clean.checkv.data)

main.clean.data <- test.output.data[[1]]
main.clean.stat.data <- test.output.data[[2]]
main.list.clean.data <- test.output.data[[3]]
main.raw.stat.data.all <- test.output.data[[4]]

#####
# # Setup main data structures
# main.list.clean.data <- list()
# main.clean.data <- data.frame(matrix(ncol = 4, nrow = 0))
# main.clean.stat.data <- data.frame(matrix(ncol = 15, nrow = 0))
# main.raw.stat.data.all <- data.frame(matrix(ncol = 13, nrow = 0))
# colnames(main.clean.data) <- c("SRR_ID","genome_ID","variable","value")
# 
# #####
# # Loading data in data loop
# for (SRR_i in 1:length(input.files.SRR.vector)){
#     skip_to_next <<- F
#     
#     # Set keys
#     temp.SRR_ID <- input.files.SRR.vector[SRR_i]
#     temp_instrain_dir <- input.files.instrain.list[SRR_i]
#     
#     # Load instrain data
#     tryCatch({
#         raw.instrain.data <- read.csv(file = paste(data_wd, input.file.instrain.dir, temp_instrain_dir, "output", 
#                                                    paste0(temp_instrain_dir, "_scaffold_info.tsv"), sep = "/"), 
#                                       sep = "\t", header = T)
#     }, error = function(e) {skip_to_next <<- TRUE})
#     if(skip_to_next){
#         rm(skip_to_next)
#         next
#     }
#     
#     # Clean instrain data
#     clean.instrain.data <- raw.instrain.data %>%
#         dplyr::mutate(short_genome_ID = str_split_i(scaffold, "_length_", i = 1), .after = 1) %>%
#         dplyr::mutate(short_genome_ID = str_split_i(short_genome_ID, "\\|\\|", i = 1))
#     
#     tmp.clean.data.all <- clean.instrain.data %>%
#         dplyr::distinct(short_genome_ID, .keep_all = T) %>%
#         dplyr::inner_join(clean.lifestyle.data, by = join_by(short_genome_ID == Accession), 
#                           keep = F, relationship = "one-to-one", multiple = "any") %>%
#         dplyr::inner_join(clean.checkv.data, by = join_by(short_genome_ID == contig_id),
#                           keep = F, relationship = "one-to-one", multiple = "any") %>%
#         dplyr::mutate(Type = as.factor(Type))
#     
#     tmp.clean.data <- tmp.clean.data.all %>%
#         dplyr::filter(Type %in% c("virulent", "temperate")) %>%
#         tidyr::drop_na(nucl_diversity) %>%
#         dplyr::filter(breadth > 0.8 & coverage > 5 & completeness >= 80)
#     
#     rm(list = c("raw.instrain.data", "clean.instrain.data"))
#     #####
#     tmp.p.point.coverage_breadth_byType <- tmp.clean.data %>%
#         ggplot( aes(x = coverage, y = breadth, color = Type)) +
#         geom_point(alpha = 0.7)
#     
#     tmp.p.box.length_byType <- tmp.clean.data %>%
#         ggplot( aes(x = Type, y = Length, color = Type)) +
#         geom_boxplot()
#     
#     tmp.p.point.phatypscore_length_byType <- tmp.clean.data %>%
#         ggplot( aes(x = PhaTYPScore, y = Length, color = Type)) +
#         geom_point(alpha = 0.7)
#     
#     #####
#     # Statistics
#     tmp.virulent.wt <- dplyr::filter(tmp.clean.data, Type == "virulent")$nucl_diversity
#     tmp.temperate.wt <- dplyr::filter(tmp.clean.data, Type == "temperate")$nucl_diversity
#     
#     if (!length(tmp.temperate.wt) > 10 | !length(tmp.virulent.wt) > 10){
#         rm(list = c("tmp.clean.data", "tmp.temperate.wt",
#                     "tmp.virulent.wt",
#                     "temp.SRR_ID", "temp_instrain_dir"))
#         next
#     }
#     
#     tmp.wilcox.test <- wilcox.test(tmp.virulent.wt, tmp.temperate.wt, conf.int = T)
#     
#     # Transpose and construct df for addition to main df
#     tmp.clean.data <- tmp.clean.data %>%
#         dplyr::select(short_genome_ID, nucl_diversity, length, breadth, coverage, Type) %>%
#         dplyr::rename(genome_ID = short_genome_ID,
#                       lifestyle = Type,
#                       nucl_diversity = nucl_diversity) %>%
#         dplyr::mutate(SRR_ID = temp.SRR_ID, .before = 1)
#     
#     tmp.clean.data.long <- tmp.clean.data %>%
#         tidyr::pivot_longer(cols = c("nucl_diversity", "length", 
#                                      "breadth", "coverage"),
#                             names_to = "variable", values_to = "value")
#     
#     # Calculate summary stats on lifestyle data
#     tmp.summary.stats <- tmp.clean.data %>%
#         dplyr::group_by(SRR_ID, lifestyle) %>%
#         dplyr::summarise(count = n(),
#                          avg_length = mean(length),
#                          med_length = median(length),
#                          avg_nucl_diversity = mean(nucl_diversity),
#                          med_nucl_diversity = median(nucl_diversity),
#                          avg_breadth = mean(breadth),
#                          med_breadth = median(breadth),
#                          avg_coverage = mean(coverage),
#                          med_coverage = median(coverage),
#                          quant95 = quantile(nucl_diversity, probs = 0.95),
#                          .groups = "keep") %>%
#         dplyr::mutate(wilcox_pvalue = tmp.wilcox.test$p.value) %>%
#         dplyr::mutate(wilcox_estimate = round(tmp.wilcox.test$estimate, digits = 7)) %>%
#         dplyr::mutate(lifestyle_ratio = length(tmp.temperate.wt)/length(tmp.virulent.wt), .after = count)
#     
#     # Calculate summary stats on raw data
#     tmp.summary.stats.raw <- tmp.clean.data.all %>%
#         dplyr::rename(genome_ID = short_genome_ID,
#                       lifestyle = Type,
#                       nucl_diversity = nucl_diversity) %>%
#         dplyr::mutate(SRR_ID = temp.SRR_ID, .before = 1) %>%
#         dplyr::group_by(SRR_ID, lifestyle) %>%
#         dplyr::summarise(count = n(),
#                          avg_length = mean(length),
#                          med_length = median(length),
#                          avg_nucl_diversity = mean(nucl_diversity, na.rm = T),
#                          med_nucl_diversity = median(nucl_diversity, na.rm = T),
#                          avg_breadth = mean(breadth),
#                          med_breadth = median(breadth),
#                          avg_coverage = mean(coverage),
#                          med_coverage = median(coverage),
#                          quant95 = quantile(nucl_diversity, probs = 0.95, na.rm = T),
#                          .groups = "keep")
#     
#     # Save data in main variables
#     main.list.clean.data[[temp.SRR_ID]][["raw_data"]] <- tmp.clean.data.all
#     main.list.clean.data[[temp.SRR_ID]][["filtered_data"]] <- tmp.clean.data
#     main.list.clean.data[[temp.SRR_ID]][["plots"]] <- list(p.point.coverage_breadth_byType = tmp.p.point.coverage_breadth_byType,
#                                                            p.box.length_byType = tmp.p.box.length_byType,
#                                                            p.point.phatypscore_length_byType = tmp.p.point.phatypscore_length_byType)
#     main.clean.data <- rbind(main.clean.data, tmp.clean.data.long)
#     main.clean.stat.data <- rbind(main.clean.stat.data, tmp.summary.stats)
#     main.raw.stat.data.all <- rbind(main.raw.stat.data.all, tmp.summary.stats.raw)
#     
#     rm(list = c("tmp.clean.data", "tmp.temperate.wt", "tmp.clean.data.all",
#                 "tmp.virulent.wt", "tmp.clean.data.long", "tmp.summary.stats.raw",
#                 "temp.SRR_ID", "temp_instrain_dir",
#                 "tmp.summary.stats", "tmp.wilcox.test",
#                 "tmp.p.point.coverage_breadth_byType","tmp.p.box.length_byType",
#                 "tmp.p.point.phatypscore_length_byType"))
# };rm(SRR_i)

#####

# Conversion to one big dataframe
main.clean.data.plot <- main.clean.data %>%
    dplyr::filter(variable == "nucl_diversity") %>%
    dplyr::mutate(SRR_ID_lifestyle = paste(SRR_ID, lifestyle)) %>%
    dplyr::group_by(SRR_ID)

# Single summary
main.data.summary.stats <- main.clean.data %>%
    tidyr::pivot_wider(names_from = variable, values_from = value) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(lifestyle) %>%
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
    dplyr::select(lifestyle, count, med_length, med_nucl_diversity, med_breadth, med_coverage) %>%
    tidyr::pivot_longer(cols = c("count","med_length","med_nucl_diversity","med_breadth","med_coverage"), 
                        names_to = "variable", values_to = "value")

# Statistics markup
main.clean.stat.data$wilcox_pvalue <- p.adjust(main.clean.stat.data$wilcox_pvalue, method = "BH")
main.clean.stat.data.new <- main.clean.stat.data %>%
    dplyr::mutate(markup = ifelse(wilcox_pvalue > 0.05, "",
                                  ifelse(wilcox_pvalue > 0.01, "*",
                                         ifelse(wilcox_pvalue > 0.001, "**",
                                                ifelse(wilcox_pvalue > 0.0001, "***",
                                                       ifelse(wilcox_pvalue < 0.0001, "****", NA)))))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(markup_estimate = ifelse(!markup %in% c("",NA), wilcox_estimate, "")) %>%
    dplyr::distinct(SRR_ID, .keep_all = T) %>%
    dplyr::mutate(SRR_ID = as.character(SRR_ID))

# Draw plot for general statistics
main.clean.stat.data.shared <- SharedData$new(data = main.clean.stat.data, key = ~SRR_ID)
p.bar.breadth_bysample <- main.clean.stat.data.shared %>%
    ggplot( aes(x = SRR_ID, y = avg_breadth, fill = lifestyle)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = "")

p.bar.coverage_bysample <- main.clean.stat.data.shared %>%
    ggplot( aes(x = SRR_ID, y = avg_coverage, fill = lifestyle)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = "bottom")

p.bar.length_bysample <- main.clean.stat.data.shared %>%
    ggplot( aes(x = SRR_ID, y = avg_length, fill = lifestyle)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.position = "bottom")

# Draw plot for microdiversity
p.interval.microdiversity_bysample <- main.clean.data.plot %>%
    dplyr::filter(value > 0) %>%
    ggplot( aes(x = SRR_ID, y = value, group = lifestyle)) +
    stat_interval(position = position_dodgejust(width = 0.8), linewidth = 3, width = 1, show.legend = NA,
                  colour = c(rep(head(MetBrewer::met.brewer("OKeeffe1", direction = -1), 3), length(main.clean.stat.data.new$SRR_ID)), 
                             rep(head(MetBrewer::met.brewer("OKeeffe1", direction = 1), 3), length(main.clean.stat.data.new$SRR_ID)))) +
    stat_summary(geom = "point", fun = median, position = position_dodgejust(width = 0.8), colour = "white") +
    geom_text(inherit.aes = F, data = main.clean.stat.data.new,
              aes(x = SRR_ID, y = 0.0001, label = markup), hjust = 0.5, nudge_x = 0) +
    geom_text(inherit.aes = F, data = main.clean.stat.data.new,
              aes(x = SRR_ID, y = 0.00015, label = markup_estimate), hjust = 0) +
    coord_flip(clip = "on", ylim = c(0.0001,0.1)) +
    scale_y_log10(guide = "axis_logticks") +
    theme_classic() +
    labs(title = "Microdiversity levels based on per sample instrain data",
         subtitle = "P values adjusted using BH method. \n *: p<0.05, **: p<0.01, ***: p<0.001, ****: p<0.0001",
         y = "Microdiversity (\u03c0)",
         x = "SRR IDs",
         caption = paste(paste0("Data = ", input.file.instrain.dir),
                         paste0("Bioproject = ", bioproject_ID), 
                         sep = "\n")) +
    theme(plot.background = element_rect(color = NA, fill = bg_color),
          panel.grid.major.y = element_line(linewidth = 1, linetype = 2),
          axis.text.y = element_text(hjust = 0),
          plot.title.position = "plot")

#####

# Lifestyle prediction check with thresholds
tmp.clean.data <- main.clean.data %>%
    tidyr::pivot_wider(names_from = variable, values_from = value) %>%
    dplyr::left_join(clean.lifestyle.data.all, relationship = "many-to-one",
                     by = join_by(genome_ID == Accession)) %>%
    dplyr::left_join(clean.blast.data, by = join_by(genome_ID == Accession), keep = F,
                     relationship = "many-to-one") %>%
    tidyr::drop_na(nucl_diversity) %>%
    dplyr::mutate(blasthit = ifelse(!is.na(evalue), T, F)) %>%
    dplyr::filter(breadth > 0.8 & coverage > 5 & completeness > 80) %>%
    dplyr::filter((phatyp_temperate == T & phastyle_temperate == T) | (phastyle_lytic == T & phatyp_lytic == T)) %>%
    dplyr::filter((phatyp_temperate == T & blasthit == T) | (phatyp_lytic == T & blasthit == F)) %>%
    dplyr::select(SRR_ID, genome_ID, lifestyle, nucl_diversity, 
                  length, breadth, coverage) %>%
    tidyr::pivot_longer(cols = c("nucl_diversity", "length", 
                                 "breadth", "coverage"),
                        names_to = "variable", values_to = "value") %>%
    dplyr::filter(variable == "nucl_diversity") %>%
    dplyr::select(-variable) %>%
    dplyr::group_by(SRR_ID, lifestyle)

tmp.clean.stat.data <- tmp.clean.data %>%
    dplyr::filter(SRR_ID %in% as.character({dplyr::summarise(., count = n(), .groups = "keep") %>% dplyr::filter(count > 10)}$SRR_ID)) %>%
    base::split(~SRR_ID, drop = T) %>%
    purrr::map(~ wilcox.test(value ~ lifestyle, data = ., conf.int = T)) %>%
    purrr::map_dfr( ~ broom::tidy(.), .id = "SRR_ID")

tmp.clean.stat.data$p.value.adj <- p.adjust(tmp.clean.stat.data$p.value, method = "BH")
#tmp.clean.data$SRR_ID <- unique(sort(main.clean.data$SRR_ID))
tmp.clean.stat.data <- tmp.clean.stat.data %>%
    dplyr::mutate(markup = ifelse(p.value.adj > 0.05, "",
                                  ifelse(p.value.adj > 0.01, "*",
                                         ifelse(p.value.adj > 0.001, "**",
                                                ifelse(p.value.adj > 0.0001, "***",
                                                       ifelse(p.value.adj < 0.0001, "****", NA)))))) %>%
    dplyr::ungroup() %>%
    #dplyr::mutate(estimate = round(estimate, 2)) %>%
    dplyr::mutate(markup_estimate = ifelse(!markup %in% c("",NA), round(estimate, digits = 5), "")) %>%
    #dplyr::mutate(markup_estimate = round(as.numeric(markup_estimate), digits = 5)) %>%
    dplyr::distinct(SRR_ID, .keep_all = T) %>%
    dplyr::mutate(SRR_ID = as.character(SRR_ID))

p.interval.microdiversity_bysample <- tmp.clean.data %>%
    dplyr::filter(value > 0 & SRR_ID %in% tmp.clean.stat.data$SRR_ID) %>%
    ggplot( aes(x = SRR_ID, y = value, group = lifestyle)) +
    stat_interval(position = position_dodgejust(width = 0.8), linewidth = 3, width = 1, show.legend = NA,
                  colour = c(rep(head(MetBrewer::met.brewer("OKeeffe1", direction = -1), 3), length(tmp.clean.stat.data$SRR_ID)), 
                             rep(head(MetBrewer::met.brewer("OKeeffe1", direction = 1), 3), length(tmp.clean.stat.data$SRR_ID)))) +
    stat_summary(geom = "point", fun = median, position = position_dodgejust(width = 0.8), colour = "white") +
    geom_text(inherit.aes = F, data = tmp.clean.stat.data,
              aes(x = SRR_ID, y = 0.00008, label = markup), hjust = 0.5, nudge_x = 0) +
    geom_text(inherit.aes = F, data = tmp.clean.stat.data,
              aes(x = SRR_ID, y = 0.0001, label = markup_estimate), hjust = 0) +
    coord_flip(clip = "on", ylim = c(0.0001,0.02)) +
    scale_y_log10(guide = "axis_logticks", expand = expansion(c(0.1,0))) +
    theme_classic() +
    labs(title = paste0("Microdiversity levels based on per sample instrain (",dataset_id,") data"),
         y = "Microdiversity (log10(\u03c0))",
         x = "SRR IDs",
         caption = paste(paste0("Data = ", input.file.instrain.dir),
                         paste0("Bioproject = ", bioproject_ID, "\nP values adjusted using BH method. *: p<0.05, **: p<0.01, ***: p<0.001, ****: p<0.0001"), 
                         sep = " | ")) +
    theme(plot.background = element_rect(color = NA),
          panel.grid.major.y = element_line(linewidth = 1, linetype = 2),
          axis.text.y = element_text(hjust = 0),
          plot.title.position = "plot")
p.interval.microdiversity_bysample

ggsave(plot = p.interval.microdiversity_bysample,
       filename = paste(sys.args$data_wd, paste(dataset_id,"microdiversity_interval_plot.png", sep = "_"), sep = "/"),
       width = 1600, height = 1700,
       units = "px", dpi = 400, scale = 1.5)

#####
ComplexUpset::upset(tmp.clean.data, intersect = c("phatyp_lytic","phatyp_temperate",
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
                    annotations = list(
                        # "Coverage dist."=list(
                        #     aes=aes(x=intersection, y = coverage),
                        #     geom=list(
                        #         geom_boxplot(),
                        #         coord_cartesian(ylim = c(0, 150)),
                        #         theme(axis.title.y = element_text(angle = 70, vjust = 0.5))
                        #     )
                        # ),
                        # "Completeness dist."=list(
                        #     aes=aes(x=intersection, y = completeness),
                        #     geom=list(
                        #         geom_violin(),
                        #         coord_cartesian(ylim = c(80,100)),
                        #         theme(axis.title.y = element_text(angle = 70, vjust = 0.5))
                        #     )
                        # )
                        # ,"Breadth dist."=list(
                        #     aes=aes(x=intersection, y = breadth),
                        #     geom=list(
                        #         geom_violin(),
                        #         coord_cartesian(ylim = c(0.8, 1)),
                        #         theme(axis.title.y = element_text(angle = 70, vjust = 0.5))
                        #     )
                        # )
                    ),
                    set_sizes=(
                        upset_set_size(
                            geom=geom_bar(aes(fill = blasthit), width=0.4))
                        + geom_label(aes(label=after_stat(count)), size = 2.5, nudge_y = -1000,
                                     hjust = 0.5, stat='count', label.size = unit(0, "lines"))
                        + theme(axis.text.x=element_text(angle=45, vjust = 0.5),
                                axis.ticks.x=element_line(),
                                axis.title.x = element_blank(),
                                legend.position = "none")
                        + ylim(20000, 0)
                    ),
                    sort_sets=F) & theme(plot.background = element_rect(fill = "grey97", colour = NA))

ggsave(filename = paste(wd, paste0(input.file.instrain.dir, "_lifestyle_upset_plot.png"), sep = "/"),
       width = 3000, height = 1500,
       units = "px", dpi = 400)

#####
# Check for unique vOTUs in each SRR_ID
main.data.summary.all <- main.clean.data %>%
    dplyr::ungroup() %>%
    dplyr::filter(variable == "length") %>%
    dplyr::mutate(SRR_ID = as.character(SRR_ID)) %>%
    dplyr::group_by(genome_ID) %>%
    dplyr::summarise(unique_count = n()) %>%
    dplyr::filter(unique_count == 1)

tmp.unique.data <- main.clean.data %>%
    dplyr::mutate(unique = ifelse(genome_ID %in% main.data.summary.all$genome_ID, T, F)) %>%
    dplyr::mutate(lifestyle = as.factor(lifestyle)) %>%
    dplyr::group_by(SRR_ID, lifestyle) %>%
    dplyr::filter(unique == T) %>%
    dplyr::mutate(SRR_ID = as.factor(as.character(SRR_ID))) %>%
    tidyr::pivot_wider(names_from = variable, values_from = value) %>%
    dplyr::ungroup() %>%
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
                     .groups = "keep")

tmp.clean.data <- main.clean.data %>%
    dplyr::mutate(unique = ifelse(genome_ID %in% main.data.summary.all$genome_ID, T, F)) %>%
    dplyr::filter(variable == "nucl_diversity", unique == T) %>%
    dplyr::mutate(SRR_ID = as.character(SRR_ID)) %>%
    dplyr::group_by(SRR_ID, lifestyle) %>%
    dplyr::select(-variable) %>%
    dplyr::filter(SRR_ID %in% as.character({dplyr::summarise(., count = n(), .groups = "keep") %>% 
            dplyr::group_by(SRR_ID) %>%
            dplyr::filter(any(lifestyle == "virulent") & any(lifestyle == "temperate")) %>%
            dplyr::filter(all(count > 10))}$SRR_ID)) %>%
    dplyr::ungroup() %>%
    base::split(~SRR_ID, drop = T) %>%
    purrr::map(~ wilcox.test(value ~ lifestyle, data = ., conf.int = T)) %>%
    purrr::map_dfr( ~ broom::tidy(.), .id = "SRR_ID")

tmp.clean.data$p.value.adj <- p.adjust(tmp.clean.data$p.value, method = "BH")
#tmp.clean.data$SRR_ID <- unique(sort(main.clean.data$SRR_ID))
tmp.clean.data <- tmp.clean.data %>%
    dplyr::mutate(markup = ifelse(p.value.adj > 0.05, "NS",
                                  ifelse(p.value.adj > 0.01, "*",
                                         ifelse(p.value.adj > 0.001, "**",
                                                ifelse(p.value.adj > 0.0001, "***",
                                                       ifelse(p.value.adj < 0.0001, "****", NA)))))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(estimate = round(estimate, digits = 7)) %>%
    dplyr::mutate(markup_estimate = ifelse(!markup %in% c("","NS",NA), estimate, "")) %>%
    dplyr::distinct(SRR_ID, .keep_all = T) %>%
    dplyr::mutate(SRR_ID = as.character(SRR_ID))

p.interval.microdiversity_bysample_unique <- main.clean.data %>%
    dplyr::mutate(unique = ifelse(genome_ID %in% main.data.summary.all$genome_ID, T, F)) %>%
    dplyr::group_by(SRR_ID, lifestyle) %>%
    dplyr::filter(unique == T) %>%
    dplyr::filter(variable == "nucl_diversity") %>%
    dplyr::filter(value > 0) %>%
    dplyr::filter(SRR_ID %in% tmp.clean.data$SRR_ID) %>%
    ggplot( aes(x = SRR_ID, y = value, group = lifestyle)) +
    stat_interval(position = position_dodgejust(width = 0.8), linewidth = 3, width = 1, show.legend = NA,
                  colour = c(rep(head(MetBrewer::met.brewer("OKeeffe1", direction = -1), 3), length(tmp.clean.data$SRR_ID)), 
                             rep(head(MetBrewer::met.brewer("OKeeffe1", direction = 1), 3), length(tmp.clean.data$SRR_ID)))) +
    stat_summary(geom = "point", fun = median, position = position_dodgejust(width = 0.8), colour = "white") +
    geom_text(inherit.aes = F, data = tmp.clean.data,
              aes(x = SRR_ID, y = 0.0001, label = markup), hjust = 0.5, nudge_x = 0) +
    geom_text(inherit.aes = F, data = tmp.clean.data,
              aes(x = SRR_ID, y = 0.00015, label = markup_estimate), hjust = 0) +
    coord_flip(clip = "on", ylim = c(0.0001,0.05)) +
    scale_y_log10(guide = "axis_logticks") +
    theme_classic() +
    labs(title = "Microdiversity levels based on per sample instrain data",
         subtitle = "P values adjusted using BH method. Top bar = virulent, bottom bar = temperate\n *: p<0.05, **: p<0.01, ***: p<0.001, ****: p<0.0001",
         y = "Microdiversity (\u03c0)",
         x = "SRR IDs",
         caption = paste(paste0("Data = ", input.file.instrain.dir),
                         paste0("Bioproject = ", bioproject_ID), 
                         sep = "\n")) +
    theme(plot.background = element_rect(color = NA, fill = bg_color),
          panel.grid.major.y = element_line(linewidth = 1, linetype = 2),
          axis.text.y = element_text(hjust = 0),
          plot.title.position = "plot")
p.interval.microdiversity_bysample_unique
p.interval.microdiversity_bysample

#####
# Bootstrapping/randomizing vOTU selection

# Setup data structures
bootstrap.clean.data <- data.frame(matrix(ncol = 16, nrow = 0))
bootstrap.clean.data.all <- data.frame(matrix(ncol = 8, nrow = 0))
sample_size <- 50
bootstrap_iterations <- 50

# Go through SRR_IDs
for (bootstrap_i in 1:bootstrap_iterations){
    for (SRR_i in 1:length(main.list.clean.data)){
        skip_to_next <<- F
        
        # Set names and select df
        tmp.SRR_ID <- names(main.list.clean.data)[SRR_i]
        tmp.SRR_ID.clean.data <- main.list.clean.data[SRR_i][[tmp.SRR_ID]][["filtered_data"]] %>%
            dplyr::mutate(bootstrap_iteration = bootstrap_i, .after = SRR_ID)
        
        # Create sub-dfs
        tmp.SRR_ID.clean.data.virulent <- dplyr::filter(tmp.SRR_ID.clean.data, lifestyle == "virulent")
        tmp.SRR_ID.clean.data.temperate <- dplyr::filter(tmp.SRR_ID.clean.data, lifestyle == "temperate")
        
        # Set sample size
        if (nrow(tmp.SRR_ID.clean.data.virulent) < nrow(tmp.SRR_ID.clean.data.temperate)){
            sample_size_length <- nrow(tmp.SRR_ID.clean.data.virulent)
        } else {
            sample_size_length <- nrow(tmp.SRR_ID.clean.data.temperate)
        }
        sample_size <- as.integer(sample_size_length * 0.25)
        
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

bootstrap.clean.data$p.value.adj <- p.adjust(bootstrap.clean.data$wilcox_pvalue, method = "BH")
bootstrap.clean.data <- bootstrap.clean.data %>%
    dplyr::group_by(bootstrap_iteration) %>%
    dplyr::mutate(p.value.adj = p.adjust(wilcox_pvalue, method='BH')) %>%
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
p.heatmap.microdiversity_bysample.bootstrapped <- bootstrap.clean.data %>%
    ggplot( aes(y = SRR_ID, x = bootstrap_iteration, fill = estimate)) +
    geom_tile() +
    geom_text(inherit.aes = F, data = bootstrap.clean.data,
              aes(y = SRR_ID, x = bootstrap_iteration, label = markup), color = "black", hjust = 0.5, nudge_x = 0, angle = 30) +
    scale_fill_gradientn(colours = MetBrewer::met.brewer("OKeeffe1", type = "continuous", direction = -1), limits = c(- 0.005, 0.005), na.value = "grey80") +
    theme_classic() +
    labs(title = "Estimated difference in microdiversity between samples (bootstrapped)",
         subtitle = "P values adjusted using BH method. \n*: p<0.05, **: p<0.01, ***: p<0.001, ****: p<0.0001",
         x = "Bootstrap iteration",
         y = "SRR IDs",
         fill = "Estimated difference \nbetween lifestyle",
         caption = paste(paste0("Data = ", input.file.instrain.dir),
                         paste0("Bioproject = ", bioproject_ID), 
                         sep = "\n")) +
    scale_x_continuous(expand = expansion(c(0,0))) +
    theme(legend.position = "right",
          plot.caption.position = "plot",
          plot.title.position = "plot")
p.heatmap.microdiversity_bysample.bootstrapped


