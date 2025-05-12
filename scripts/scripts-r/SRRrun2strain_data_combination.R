# Analyse PhaBOX output (lifestyle)
##############
### setup
rm(list = ls())
wd <- getwd()
wd <- c("/Users/thomasdebruijn/Documents/PhD/R_PhD")
setwd(wd)
set.seed(1)
options(scipen = 6)
bg_color <- "grey97"
font_family <- "Fira Sans"

# Install all necessary packages
install <- FALSE
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
})

# Custom SRR ID function
SRR_extract <- function(in.string){
    SRR.pattern <- "^SRR\\d{6,10}$"
    in.string = str_split_1(in.string, "_")
    out.string = in.string[str_detect(in.string, SRR.pattern)]
    out.string = as.character(out.string)
    return(out.string)
}

# Load command arguments
#sys.args <- commandArgs_custom

# Set files to load
#wd <- sys.args$wd
wd <- "/Users/thomasdebruijn/Documents/PhD/DATASETS"
#input.file.instrain.dir <- sys.args$instrain_dir
input.file.instrain.dir <- "instrain_test_dir"
input.files.instrain.list <- list.files(path = paste(wd, input.file.instrain.dir, sep = "/"),
                                        full.names = F)
input.files.SRR.vector <- sapply(unlist(input.files.instrain.list), SRR_extract)

# Load Phatyp predictions
raw.lifestyle.landuse.data <- read.csv(file = "../DATASETS/landuse_phatyp_prediction.tsv", sep = "\t", header = T)
raw.lifestyle.santosv.data <- read.csv(file = "../DATASETS/santosviromes_phatyp_prediction.tsv", sep = "\t", header = T)
clean.lifestyle.data <- raw.lifestyle.landuse.data %>%
    rbind(raw.lifestyle.santosv.data) %>%
    dplyr::mutate(Type = as.factor(TYPE), .keep = "unused") %>%
    dplyr::mutate(Accession = str_split_i(Accession, "_length_", i = 1))
rm("raw.lifestyle.landuse.data","raw.lifestyle.santosv.data")

# Setup main data structures
main.list.clean.data <- list()
main.clean.data <- data.frame(matrix(ncol = 4, nrow = 0))
main.clean.stat.data <- data.frame(matrix(ncol = 13, nrow = 0))
colnames(main.clean.data) <- c("SRR_ID","genome_ID","variable","value")

#####
# Loading data in data loop
for (SRR_i in 1:length(input.files.SRR.vector)){
    # Set keys
    temp.SRR_ID <- input.files.SRR.vector[SRR_i]
    temp_instrain_dir <- input.files.instrain.list[SRR_i]
    
    # Load instrain data
    raw.instrain.data <- read.csv(file = paste(wd, input.file.instrain.dir, temp_instrain_dir, "output", 
                                               paste0(temp_instrain_dir, "_scaffold_info.tsv"), sep = "/"), 
                                  sep = "\t", header = T)
    
    # Clean instrain data
    clean.instrain.data <- raw.instrain.data %>%
        dplyr::mutate(short_genome_ID = str_split_i(scaffold, "_length_", i = 1), .after = 1)
    
    clean.data <- clean.instrain.data %>%
        dplyr::distinct(short_genome_ID, .keep_all = T) %>%
        dplyr::inner_join(clean.lifestyle.data, by = join_by(short_genome_ID == Accession), 
                          keep = F, relationship = "one-to-one", multiple = "any") %>%
        dplyr::filter(nucl_diversity_rarefied > 0) %>%
        dplyr::filter(Type %in% c("virulent", "temperate")) %>%
        dplyr::mutate(Type = as.factor(Type)) %>%
        dplyr::filter(breadth > 0.8 & coverage > 10 & PhaTYPScore > 0.8)
    
    rm(list = c("raw.instrain.data", "clean.instrain.data"))
    ###
    # clean.data %>%
    #     ggplot( aes(x = coverage, y = breadth, color = Type)) +
    #     geom_point() +
    #     coord_cartesian(ylim = c(0, 1))
    # 
    # clean.data %>%
    #     ggplot( aes(x = coverage, y = breadth, color = Type)) +
    #     geom_point() +
    #     coord_cartesian(ylim = c(0.8, 1))
    # 
    # clean.data %>%
    #     ggplot( aes(x = Type, y = Length, color = Type)) +
    #     geom_boxplot() +
    #     theme(legend.position = "none")
    # 
    # clean.data %>%
    #     ggplot( aes(x = PhaTYPScore, y = Length, color = Type)) +
    #     geom_point() +
    #     theme(legend.position = "none")
    # 
    # clean.data %>%
    #     ggplot( aes(x = coverage, y = breadth, color = Type)) +
    #     geom_point(alpha = 0.5) +
    #     coord_cartesian(xlim = c(0, 100))
    # 
    # clean.data %>%
    #     ggplot( aes(x = Type, y = nucl_diversity, color = Type)) +
    #     geom_boxplot()
    # 
    # clean.data %>%
    #     ggplot( aes(x = Type, y = nucl_diversity, fill = Type)) +
    #     geom_violin() +
    #     geom_boxplot(alpha = 0.5) +
    #     scale_y_continuous(trans=log2_trans(),
    #                        breaks = trans_breaks("log2", function(x) 2^x),
    #                        labels = trans_format("log2", math_format(2^.x))) +
    #     annotation_logticks(sides = "l")
    # 
    # clean.data %>%
    #     ggplot( aes(x = Type, y = nucl_diversity_rarefied, fill = Type)) +
    #     geom_violin() +
    #     geom_boxplot(alpha = 0.5) +
    #     scale_y_continuous(transform = log2_trans(),
    #                        breaks = trans_breaks("log2", function(x) 2^x),
    #                        labels = trans_format("log2", math_format(2^.x))) +
    #     annotation_logticks(sides = "l")
    
    #var.test(nucl_diversity_rarefied ~ Type, data = clean.data)
    
    # clean.data %>%
    #     ggplot( aes(nucl_diversity, fill = Type)) +
    #     geom_density(alpha = 0.6) +
    #     scale_x_continuous(transform = log2_trans(),
    #                        breaks = trans_breaks("log2", function(x) 2^x),
    #                        labels = trans_format("log2", math_format(2^.x)))
    # 
    # clean.data %>%
    #     ggplot( aes(x = coverage, y = breadth, color = Type)) +
    #     coord_cartesian(xlim = c(0, 2000)) +
    #     geom_point(alpha = 0.5)
    # 
    # clean.data %>%
    #     dplyr::mutate(relative_divergent_site_count = divergent_site_count / length) %>%
    #     ggplot( aes(x = relative_divergent_site_count, fill = Type)) +
    #     geom_density(alpha = 0.5)
    
    # Statistics
    tmp.virulent.wt <- dplyr::filter(clean.data, Type == "virulent")$nucl_diversity_rarefied
    tmp.temperate.wt <- dplyr::filter(clean.data, Type == "temperate")$nucl_diversity_rarefied
    
    if (!length(tmp.temperate.wt) > 50 | !length(tmp.virulent.wt) > 50){
        next
    }
    
    # cat(c(temp.SRR_ID, "\n"))
    # cat(paste(paste("Sample size virulent:", length(tmp.virulent.wt), sep = " "),
    #           paste("Median virulent nucl. diversity:", median(tmp.virulent.wt), "\n", sep = " "),
    #           sep = ", "))
    # cat(paste(paste("Sample size temperate:", length(tmp.temperate.wt), sep = " "),
    #           paste("Median temperate nucl. diversity:", median(tmp.temperate.wt), "\n", sep = " "),
    #           sep = ", "))
    tmp.wilcox.test <- wilcox.test(tmp.virulent.wt, tmp.temperate.wt, conf.int = T)
    
    # Transpose and construct df for addition to main df
    clean.data <- clean.data %>%
        dplyr::select(short_genome_ID, nucl_diversity_rarefied, length, breadth, coverage, Type) %>%
        dplyr::rename(genome_ID = short_genome_ID,
                      lifestyle = Type,
                      nucl_diversity = nucl_diversity_rarefied) %>%
        dplyr::mutate(SRR_ID = temp.SRR_ID, .before = 1)
    
    tmp.clean.data <- clean.data %>%
        tidyr::pivot_longer(cols = c("nucl_diversity", "length", 
                                     "breadth", "coverage"),
                            names_to = "variable", values_to = "value")
        
    
    tmp.summary.stats <- clean.data %>%
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
                         quant95 = quantile(nucl_diversity, probs = 0.98),
                         .groups = "keep") %>%
        dplyr::mutate(wilcox_pvalue = tmp.wilcox.test$p.value)
    
    # Save data in main variables
    main.list.clean.data[[temp.SRR_ID]] <- clean.data
    main.clean.data <- rbind(main.clean.data, tmp.clean.data)
    main.clean.stat.data <- rbind(main.clean.stat.data, tmp.summary.stats)
    
    rm(list = c("clean.data", "tmp.temperate.wt",
                "tmp.virulent.wt", "tmp.clean.data",
                "temp.SRR_ID", "temp_instrain_dir",
                "tmp.summary.stats", "tmp.wilcox.test"))
};rm(SRR_i)

#####

# Conversion to one big dataframe
main.clean.data.plot <- main.clean.data %>%
    dplyr::filter(variable == "nucl_diversity") %>%
    dplyr::mutate(SRR_ID_lifestyle = paste(SRR_ID, lifestyle)) %>%
    dplyr::group_by(SRR_ID)

# Statistics markup
main.clean.stat.data.new <- main.clean.stat.data %>%
    dplyr::mutate(markup = ifelse(wilcox_pvalue > 0.05, "NS",
                                  ifelse(wilcox_pvalue > 0.01, "*",
                                         ifelse(wilcox_pvalue > 0.001, "**",
                                                ifelse(wilcox_pvalue > 0.0001, "***",
                                                       ifelse(wilcox_pvalue < 0.0001, "****", NA)))))) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(SRR_ID, .keep_all = T) %>%
    dplyr::mutate(SRR_ID = as.character(SRR_ID))

# Draw plot
plot.data <- main.clean.data.plot %>%
    ggplot( aes(x = SRR_ID, y = value, group = lifestyle)) +
    stat_interval(position = position_dodgejust(width = 0.8), linewidth = 3, width = 1, show.legend = NA,
                  colour = c(rep(MetBrewer::met.brewer("OKeeffe2", n = 3, direction = -1), length(main.clean.stat.data.new$SRR_ID)), 
                             rep(MetBrewer::met.brewer("Peru2", n = 3), length(main.clean.stat.data.new$SRR_ID)))) +
    stat_summary(geom = "point", fun = median, position = position_dodgejust(width = 0.8), colour = "white") +
    stat_summary(geom = "point", fun = mean, position = position_dodgejust(width = 0.8), colour = "white", shape = 17) +
    geom_text(inherit.aes = F, data = main.clean.stat.data.new, 
              aes(x = SRR_ID, y = max(quant95)*0.9, label = markup)) +
    coord_flip(clip = "on") +
    scale_y_continuous(transform = log2_trans(),
                       breaks = trans_breaks("log2", function(x) 2^x),
                       labels = trans_format("log2", math_format(2^.x))) +
    theme_classic() +
    theme(plot.background = element_rect(color = NA, fill = bg_color),
          panel.grid.major.y = element_line(linewidth = 1, linetype = 2),
          axis.text.y = element_text(hjust = 0)) +
    scale_y_continuous(n.breaks = 10) +
    labs(y = "Microdiversity (\u03c0)",
         x = "SRR IDs")
plot.data





