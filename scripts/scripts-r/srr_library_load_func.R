# Custom SRR ID function
SRR_extract <- function(in.string){
  SRR.pattern <- "^SRR\\d{6,10}$"
  in.string = str_split_1(in.string, "_|-")
  out.string = in.string[str_detect(in.string, SRR.pattern)]
  out.string = as.character(out.string)
  
  ifelse(length(out.string) > 1,
                out.string <- out.string[1],
                out.string)
  
  return(out.string)
}

# Main SRR function for loading instrain scaffold output by SRR ID
load_srr_libraries <- function(srr_file_vector, lifestyle_df, checkv_df,
                               thresholds = list(
                                 breadth = 0.8,
                                 coverage = 5,
                                 completeness = 80), 
                               minimum_votus = 10){
  
  #' Title: load_srr_libraries
  #'
  #' Description:
  #' Loads a specified list of instrain output files by SRR library
  #'
  #' Usage:
  #' load_srr_libraries(srr_file_vector = list(), lifestyle_df = df())
  #'
  #' Arguments:
  #' @param srr_file_vector character vector, each string an instrain file
  #' @param lifestyle_df data frame containing lifestyle details per vOTU
  #' @param checkv_df data frame containing checkv completeness estimation
  #'
  #' Returns:
  #' Combined filtered instrain output files in several data frames
  #'
  #'
  #' Author: Thomas de Bruijn
  #' Date: 2025-07-09
  #' Version: 0.1
  #'
  #' Dependencies:
  #' tidyverse
  #'
  #' @export
  
  # Setup main data structures
  main.list.clean.data <- list()
  main.clean.data <- data.frame(matrix(ncol = 4, nrow = 0)) |>
    purrr::set_names(c("SRR_ID","genome_ID","variable","value"))
  main.clean.stat.data <- data.frame(matrix(ncol = 15, nrow = 0))
  main.raw.stat.data.all <- data.frame(matrix(ncol = 13, nrow = 0))
  
  #####
  # Loading data in data loop
  for (SRR_i in 1:length(srr_file_vector)){
    skip_to_next <<- F
    
    # Set keys
    temp.SRR_ID <- SRR_extract(srr_file_vector[SRR_i])
    temp.instrain.dir <- srr_file_vector[SRR_i]
    # Load instrain data
    tryCatch({
      raw.instrain.data <- read.csv(file = paste0(temp.instrain.dir, "_scaffold_info.tsv"), 
                                    sep = "\t", header = T)
    }, error = function(e) {skip_to_next <<- TRUE})
    if(skip_to_next){
      rm(skip_to_next)
      next
    }
    
    # Clean instrain data
    clean.instrain.data <- raw.instrain.data |>
      dplyr::mutate(short_genome_ID = str_split_i(scaffold, "_length_", i = 1), .after = 1) |>
      dplyr::mutate(short_genome_ID = str_split_i(short_genome_ID, "\\|\\|", i = 1))
    
    tmp.clean.data.all <- clean.instrain.data |>
      dplyr::distinct(short_genome_ID, .keep_all = T) |>
      dplyr::inner_join(lifestyle_df, by = join_by(short_genome_ID == Accession), 
                        keep = F, relationship = "one-to-one", multiple = "any") |>
      dplyr::inner_join(checkv_df, by = join_by(short_genome_ID == contig_id),
                        keep = F, relationship = "one-to-one", multiple = "any") |>
      dplyr::mutate(Type = as.factor(Type))
    
    tmp.clean.data <- tmp.clean.data.all |>
      dplyr::filter(Type %in% c("virulent", "temperate")) |>
      tidyr::drop_na(nucl_diversity) |>
      dplyr::filter(breadth > thresholds$breadth & 
                      coverage > thresholds$coverage & 
                      completeness >= thresholds$completeness)
    
    rm("raw.instrain.data", "clean.instrain.data")
    #####
    tmp.p.point.coverage_breadth_byType <- tmp.clean.data |>
      ggplot( aes(x = coverage, y = breadth, color = Type)) +
      geom_point(alpha = 0.7)
    
    tmp.p.box.length_byType <- tmp.clean.data |>
      ggplot( aes(x = Type, y = Length, color = Type)) +
      geom_boxplot()
    
    tmp.p.point.phatypscore_length_byType <- tmp.clean.data |>
      ggplot( aes(x = PhaTYPScore, y = Length, color = Type)) +
      geom_point(alpha = 0.7)
    
    #####
    # Statistics
    tmp.virulent.wt <- dplyr::filter(tmp.clean.data, Type == "virulent")$nucl_diversity
    tmp.temperate.wt <- dplyr::filter(tmp.clean.data, Type == "temperate")$nucl_diversity
    
    if (!length(tmp.temperate.wt) > minimum_votus | !length(tmp.virulent.wt) > minimum_votus){
      rm(list = c("tmp.clean.data", "tmp.temperate.wt", "tmp.clean.data.all",
                  "tmp.virulent.wt", "clean.instrain.data", "raw.instrain.data",
                  "temp.SRR_ID", "temp.instrain.dir"))
      next
    }
    
    tmp.wilcox.test <- wilcox.test(tmp.virulent.wt, tmp.temperate.wt, conf.int = T)
    
    # Transpose and construct df for addition to main df
    tmp.clean.data <- tmp.clean.data |>
      dplyr::select(short_genome_ID, nucl_diversity, length, breadth, coverage, Type) |>
      dplyr::rename(genome_ID = short_genome_ID,
                    lifestyle = Type,
                    nucl_diversity = nucl_diversity) |>
      dplyr::mutate(SRR_ID = temp.SRR_ID, .before = 1)
    
    tmp.clean.data.long <- tmp.clean.data |>
      tidyr::pivot_longer(cols = c("nucl_diversity", "length", 
                                   "breadth", "coverage"),
                          names_to = "variable", values_to = "value")
    
    # Calculate summary stats on lifestyle data
    tmp.summary.stats <- tmp.clean.data |>
      dplyr::group_by(SRR_ID, lifestyle) |>
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
      dplyr::mutate(wilcox_pvalue = tmp.wilcox.test$p.value) |>
      dplyr::mutate(wilcox_estimate = round(tmp.wilcox.test$estimate, digits = 7)) |>
      dplyr::mutate(lifestyle_ratio = length(tmp.temperate.wt)/length(tmp.virulent.wt), .after = count)
    
    # Calculate summary stats on raw data
    tmp.summary.stats.raw <- tmp.clean.data.all |>
      dplyr::rename(genome_ID = short_genome_ID,
                    lifestyle = Type,
                    nucl_diversity = nucl_diversity) |>
      dplyr::mutate(SRR_ID = temp.SRR_ID, .before = 1) |>
      dplyr::group_by(SRR_ID, lifestyle) |>
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
                "temp.SRR_ID", "temp.instrain.dir",
                "tmp.summary.stats", "tmp.wilcox.test",
                "tmp.p.point.coverage_breadth_byType","tmp.p.box.length_byType",
                "tmp.p.point.phatypscore_length_byType"))
  };rm(SRR_i)
  
  return(list(main.clean.data, main.clean.stat.data, 
              main.list.clean.data, main.raw.stat.data.all))
}
