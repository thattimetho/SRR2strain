# Analyse SRRun2strain data per study/publication
##############
### setup
rm(list = setdiff(ls(), c("master_script.data", "commandArgs_custom")))

# commandArgs_custom <- list(
#   instrain_dir = "instrain_glacial",
#   lifestyle_file = "glacial_phatyp_prediction.tsv",
#   wd = "/Users/thomasdebruijn/Documents/PhD/R_PhD",
#   data_wd = "/Users/thomasdebruijn/Documents/PhD/DATASETS",
#   install = F
# )

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
})

# Custom SRR ID function
SRR_extract <- function(in.string){
  SRR.pattern <- "^SRR\\d{6,10}$"
  in.string = str_split_1(in.string, "_|-")
  out.string = in.string[str_detect(in.string, SRR.pattern)]
  out.string = as.character(out.string)
  return(out.string)
}

# Set files and params to load
wd <- sys.args$wd
setwd(wd)
data_wd <- sys.args$data_wd
input.file.instrain.dir <- sys.args$instrain_dir
input.file.lifestyle.name <- sys.args$lifestyle_file
bioproject_ID <- sys.args$bioproject_id

#####
### Data loading
input.files.instrain.list <- list.files(path = paste(data_wd, input.file.instrain.dir, sep = "/"),
                                        full.names = F)
input.files.SRR.vector <- sapply(unlist(input.files.instrain.list), SRR_extract)

# Load Phatyp predictions and clean
raw.lifestyle.data <- read.csv(file = paste(data_wd, input.file.lifestyle.name, sep = "/"), 
                               sep = "\t", header = T)

clean.lifestyle.data <- raw.lifestyle.data %>%
  dplyr::mutate(Type = as.factor(TYPE), .keep = "unused") %>%
  dplyr::mutate(Accession = str_split_i(Accession, "_length_", i = 1))
rm("raw.lifestyle.data")

# Setup main data structures
main.list.clean.data <- list()
main.clean.data <- data.frame(matrix(ncol = 4, nrow = 0))
main.clean.stat.data <- data.frame(matrix(ncol = 15, nrow = 0))
main.raw.stat.data.all <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(main.clean.data) <- c("SRR_ID","genome_ID","variable","value")

#####
# Loading data in data loop
for (SRR_i in 1:length(input.files.SRR.vector)){
  skip_to_next <<- F
  
  # Set keys
  temp.SRR_ID <- input.files.SRR.vector[SRR_i]
  temp_instrain_dir <- input.files.instrain.list[SRR_i]
  
  # Load instrain data
  tryCatch({
    raw.instrain.data <- read.csv(file = paste(data_wd, input.file.instrain.dir, temp_instrain_dir, "output", 
                                               paste0(temp_instrain_dir, "_scaffold_info.tsv"), sep = "/"), 
                                  sep = "\t", header = T)
  }, error = function(e) {skip_to_next <<- TRUE})
  if(skip_to_next){
    rm(skip_to_next)
    next
  }
  
  # Clean instrain data
  clean.instrain.data <- raw.instrain.data %>%
    dplyr::mutate(short_genome_ID = str_split_i(scaffold, "_length_", i = 1), .after = 1) %>%
    dplyr::mutate(short_genome_ID = str_split_i(short_genome_ID, "\\|\\|", i = 1))
  
  tmp.clean.data.all <- clean.instrain.data %>%
    dplyr::distinct(short_genome_ID, .keep_all = T) %>%
    dplyr::inner_join(clean.lifestyle.data, by = join_by(short_genome_ID == Accession), 
                      keep = F, relationship = "one-to-one", multiple = "any") %>%
    dplyr::mutate(Type = as.factor(Type))
  
  tmp.clean.data <- tmp.clean.data.all %>%
    dplyr::filter(Type %in% c("virulent", "temperate")) %>%
    dplyr::filter(breadth > 0.8 & coverage > 10 & PhaTYPScore > 0.8)
  
  rm(list = c("raw.instrain.data", "clean.instrain.data"))
  #####
  tmp.p.point.coverage_breadth_byType <- tmp.clean.data %>%
      ggplot( aes(x = coverage, y = breadth, color = Type)) +
      geom_point(alpha = 0.7)

  tmp.p.box.length_byType <- tmp.clean.data %>%
      ggplot( aes(x = Type, y = Length, color = Type)) +
      geom_boxplot()

  tmp.p.point.phatypscore_length_byType <- tmp.clean.data %>%
      ggplot( aes(x = PhaTYPScore, y = Length, color = Type)) +
      geom_point(alpha = 0.7)
  
  #####
  # Statistics
  tmp.virulent.wt <- dplyr::filter(tmp.clean.data, Type == "virulent")$nucl_diversity
  tmp.temperate.wt <- dplyr::filter(tmp.clean.data, Type == "temperate")$nucl_diversity
  
  if (!length(tmp.temperate.wt) > 10 | !length(tmp.virulent.wt) > 10){
    rm(list = c("tmp.clean.data", "tmp.temperate.wt",
                "tmp.virulent.wt",
                "temp.SRR_ID", "temp_instrain_dir"))
    next
  }
  
  tmp.wilcox.test <- wilcox.test(tmp.virulent.wt, tmp.temperate.wt, conf.int = T)
  
  # Transpose and construct df for addition to main df
  tmp.clean.data <- tmp.clean.data %>%
    dplyr::select(short_genome_ID, nucl_diversity, length, breadth, coverage, Type) %>%
    dplyr::rename(genome_ID = short_genome_ID,
                  lifestyle = Type,
                  nucl_diversity = nucl_diversity) %>%
    dplyr::mutate(SRR_ID = temp.SRR_ID, .before = 1)
  
  tmp.clean.data.long <- tmp.clean.data %>%
    tidyr::pivot_longer(cols = c("nucl_diversity", "length", 
                                 "breadth", "coverage"),
                        names_to = "variable", values_to = "value")
  
  # Calculate summary stats on lifestyle data
  tmp.summary.stats <- tmp.clean.data %>%
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
    dplyr::mutate(wilcox_pvalue = tmp.wilcox.test$p.value) %>%
    dplyr::mutate(wilcox_estimate = round(tmp.wilcox.test$estimate, digits = 5)) %>%
    dplyr::mutate(lifestyle_ratio = length(tmp.temperate.wt)/length(tmp.virulent.wt), .after = count)
  
  # Calculate summary stats on raw data
  tmp.summary.stats.raw <- tmp.clean.data.all %>%
    dplyr::rename(genome_ID = short_genome_ID,
                  lifestyle = Type,
                  nucl_diversity = nucl_diversity) %>%
    dplyr::mutate(SRR_ID = temp.SRR_ID, .before = 1) %>%
    dplyr::group_by(SRR_ID, lifestyle) %>%
    dplyr::summarise(count = n(),
                     avg_length = mean(length),
                     med_length = median(length),
                     avg_nucl_diversity = mean(nucl_diversity, na.rm = T),
                     med_nucl_diversity = median(nucl_diversity, na.rm = T),
                     avg_breadth = mean(breadth),
                     med_breadth = median(breadth),
                     avg_coverage = mean(coverage),
                     med_coverage = median(coverage),
                     quant95 = quantile(nucl_diversity, probs = 0.95, na.rm = T),
                     .groups = "keep")
  
  # Save data in main variables
  main.list.clean.data[[temp.SRR_ID]][["raw_data"]] <- tmp.clean.data.all
  main.list.clean.data[[temp.SRR_ID]][["filtered_data"]] <- tmp.clean.data
  main.list.clean.data[[temp.SRR_ID]][["plots"]] <- list(p.point.coverage_breadth_byType = tmp.p.point.coverage_breadth_byType,
                                                         p.box.length_byType = tmp.p.box.length_byType,
                                                         p.point.phatypscore_length_byType = tmp.p.point.phatypscore_length_byType)
  main.clean.data <- rbind(main.clean.data, tmp.clean.data.long)
  main.clean.stat.data <- rbind(main.clean.stat.data, tmp.summary.stats)
  main.raw.stat.data.all <- rbind(main.raw.stat.data.all, tmp.summary.stats.raw)
  
  rm(list = c("tmp.clean.data", "tmp.temperate.wt", "tmp.clean.data.all",
              "tmp.virulent.wt", "tmp.clean.data.long", "tmp.summary.stats.raw",
              "temp.SRR_ID", "temp_instrain_dir",
              "tmp.summary.stats", "tmp.wilcox.test",
              "tmp.p.point.coverage_breadth_byType","tmp.p.box.length_byType",
              "tmp.p.point.phatypscore_length_byType"))
};rm(SRR_i)

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
main.clean.stat.data$wilcox_pvalue <- p.adjust(main.clean.stat.data$wilcox_pvalue, method = "none")
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
                colour = c(rep(MetBrewer::met.brewer("OKeeffe2", n = 3, direction = -1), length(main.clean.stat.data.new$SRR_ID)), 
                           rep(MetBrewer::met.brewer("Peru2", n = 3), length(main.clean.stat.data.new$SRR_ID)))) +
  stat_summary(geom = "point", fun = median, position = position_dodgejust(width = 0.8), colour = "white") +
  #stat_summary(geom = "point", fun = mean, position = position_dodgejust(width = 0.8), colour = "white", shape = 17) +
  geom_text(inherit.aes = F, data = main.clean.stat.data.new,
            aes(x = SRR_ID, y = 0.0001, label = markup), hjust = 0.5, nudge_x = -0.1) +
  geom_text(inherit.aes = F, data = main.clean.stat.data.new,
            aes(x = SRR_ID, y = 0.00015, label = markup_estimate), hjust = 0) +
  coord_flip(clip = "on", ylim = c(0.0001,0.1)) +
  scale_y_log10(guide = "axis_logticks") +
  theme_classic() +
  labs(title = "Microdiversity levels based on per sample instrain data",
       subtitle = "P values not(!) adjusted using BH method. \n *: p<0.05, **: p<0.01, ***: p<0.001, ****: p<0.0001",
       y = "Microdiversity (\u03c0)",
       x = "SRR IDs",
       caption = paste(paste0("Data = ", input.file.instrain.dir),
                       paste0("Bioproject = ", bioproject_ID), 
                       sep = "\n")) +
  theme(plot.background = element_rect(color = NA, fill = bg_color),
        panel.grid.major.y = element_line(linewidth = 1, linetype = 2),
        axis.text.y = element_text(hjust = 0),
        plot.title.position = "plot")
#p.interval.microdiversity_bysample
# 
# legend.plot.data <- main.clean.data.plot %>%
#   dplyr::filter(value > 0 & SRR_ID == "SRR14194044") %>%
#   ggplot( aes(x = reorder(SRR_ID, value), y = value, group = lifestyle)) +
#   stat_interval(position = position_dodgejust(width = 0.3), linewidth = 3, width = 1, show.legend = NA,
#                 colour = c(rep(MetBrewer::met.brewer("OKeeffe2", n = 3, direction = -1), 1),
#                            rep(MetBrewer::met.brewer("Peru2", n = 3), 1))) +
#   stat_summary(geom = "point", fun = median, position = position_dodgejust(width = 0.3), colour = "white") +
#   #stat_summary(geom = "point", fun = mean, position = position_dodgejust(width = 0.3), colour = "white", shape = 17) +
#   geom_text(inherit.aes = F, data = dplyr::filter(main.clean.stat.data.new, SRR_ID == "SRR14194044"),
#             aes(x = SRR_ID, y = 0.000015, label = markup), hjust = 0.5) +
#   geom_text(inherit.aes = F, data = dplyr::filter(main.clean.stat.data.new, SRR_ID == "SRR14194044"),
#             aes(x = SRR_ID, y = 0.00002, label = markup_estimate), hjust = 0) +
#   # Legend of intervals
#   annotate(
#     "richtext",
#     x = c(1.8, 0.7, 1.95, 0.7),
#     y = c(0.01, 0.0045, 0.0018, 0.001),
#     label = c("95% of vOTUs","80% of vOTUs","50% of vOTUs<br>fall within this range", "median"),
#     fill = NA, label.size = 0, size = 3, vjust = 1,
#   ) +
#   geom_curve(
#     inherit.aes = F,
#     data = data.frame(
#       xend = c(1.6, 1.6, 0.65, 0.65, 1.6, 1.6, 0.7),
#       x = c(1.15, 0.95, 1.025, 0.85, 1.15, 0.95, 1.05),
#       yend = c(0.009, 0.009, 0.0045, 0.0045, 0.002, 0.002, 0.001),
#       y = c(0.012, 0.011, 0.0055, 0.0055, 0.0035, 0.003, 0.0017)),
#     aes(x = x, xend = xend, y = y, yend = yend),
#     stat = "unique", curvature = 0.2, size = 1, color = "grey12",
#     arrow = arrow(angle = 20, length = unit(1, "mm"), ends = "first")
#   ) +
#   # Legend of lifestyle
#   annotate(
#     "richtext",
#     x = c(1.5, 0.6),
#     y = c(0.04, 0.03),
#     label = c("Lytic","Lysogenic"),
#     fill = NA, label.size = 0, size = 3, vjust = 1,
#   ) +
#   geom_curve(
#     inherit.aes = F,
#     data = data.frame(
#       x = c(1.4, 0.6),
#       xend = c(1.12, 0.88),
#       y = c(0.03, 0.03),
#       yend = c(0.016, 0.016)),
#     aes(x = x, xend = xend, y = y, yend = yend),
#     stat = "unique", curvature = 0.2, size = 1, color = "grey12",
#     arrow = arrow(angle = 20, length = unit(1, "mm"))
#   ) +
#   coord_flip(clip = "on", ylim = c(0.0005,0.05)) +
#   scale_y_log10(guide = "axis_logticks") +
#   theme_void() +
#   labs(title = "Legend") +
#   theme(plot.background = element_rect(color = "grey30", linewidth = 0.3, fill = bg_color),
#         panel.grid.major.y = element_line(linewidth = 0.5, linetype = 2),
#         axis.text.y = element_blank(),
#         axis.ticks.x = element_blank(),
#         plot.title = element_text(hjust = 0.01))
# 
# #legend.plot.data
# p.interval.microdiversity_bysample + inset_element(legend.plot.data, l = 0.80, r = 0.99, t = 0.25, b = 0.10, clip = F, align_to = "full")
# 
# ggsave(filename = paste(wd, paste0(input.file.instrain.dir, "_hist_plot.png"), sep = "/"),
#        plot = p.interval.microdiversity_bysample + inset_element(legend.plot.data, l = 0.80, r = 0.999, t = 0.25, b = 0.10, clip = F, align_to = "full"),
#        width = 2000, height = 2200,
#        units = "px",
#        scale = 1.7)

p.box.length_bylifestyle <- main.clean.data %>%
  dplyr::filter(variable == "length") %>%
  ggplot( aes(x = lifestyle, y = value, colour = lifestyle)) +
  geom_boxplot(alpha = 0.8) +
  coord_flip() +
  scale_y_continuous(n.breaks = 10, ) +
  labs(title = "Distribution of genome (vOTU) lengths per lifestyle",
       subtitle = "data of all samples combined",
       x = "Lifestyle", y = "Genome length (vOTU)") +
  theme(plot.title.position = "plot",
        legend.position = "none",
        axis.text.x = element_text(angle = -40))

p.box.coverage_bylifestyle <- main.clean.data %>%
  dplyr::filter(variable == "coverage") %>%
  ggplot( aes(x = lifestyle, y = value, colour = lifestyle)) +
  geom_boxplot(alpha = 0.8) +
  coord_flip() +
  scale_y_continuous(n.breaks = 10, ) +
  labs(title = "Distribution of genome (vOTU) coverage per lifestyle",
       subtitle = "data of all samples combined",
       x = "Lifestyle", y = "Genome coverage (x times)") +
  theme(plot.title.position = "plot",
        legend.position = "none",
        axis.text.x = element_text(angle = -40))

p.box.breadth_bylifestyle <- main.clean.data %>%
  dplyr::filter(variable == "breadth") %>%
  ggplot( aes(x = lifestyle, y = value, colour = lifestyle)) +
  geom_boxplot(alpha = 0.8) +
  coord_flip() +
  scale_y_continuous(n.breaks = 10, ) +
  labs(title = "Distribution of genome (vOTU) coverage breadth per lifestyle",
       subtitle = "data of all samples combined",
       x = "Lifestyle", y = "Mapping/coverage breadth") +
  theme(plot.title.position = "plot",
        legend.position = "none",
        axis.text.x = element_text(angle = -40))
