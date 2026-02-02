# proteomics_utils.R
#
# Utility functions for reproducible proteomics and phosphoproteomics analysis
#
# Author: Fabian Kliem
#
# Description:
# A collection of helper functions developed during doctoral research for
# processing, quality control, statistical analysis and visualization of
# large-scale proteomics and phosphoproteomics data (MaxQuant and DIA-NN).
#
# The focus of this code is on transparent and reusable data processing
# pipelines that can be applied consistently across multiple projects.

options(stringsAsFactors = FALSE)

# =============================================================================
# Wrapper function for quality control and exploratory plots
# =============================================================================

#' Run a standard QC report for a MaxQuant proteome run
#'
#' Reads proteinGroups.txt and evidence.txt from a directory and
#' generates a multi-page PDF with common quality control plots.
#'
#' @param dir Directory containing MaxQuant output files
#' @return The input directory (invisibly)
QC_proteome <- function(dir = choose.dir()){
  #specify path to save PDF to
  destination = paste0(dir, "\\QC_" , as.character(basename(dir)), ".pdf")
  cat(paste0("Destination:\n",destination,"\n\nimporting data"))
  proteinGroups = read_Perseus_file(paste0(dir, "\\proteinGroups.txt"))
  evidence <- read_Perseus_file(paste0(dir, "\\evidence.txt"))
  cat("Data imported.")

  #open PDF
  pdf(file=destination, onefile = T)

  #save plots to PDF
  cat("Plotting quantifications.")
  print(plot_quant_proteinGroups(quant_proteinGroups(proteinGroups)))
  cat("Plotting missed cleavages.")
  print(plot_missed_cleavages(evidence))
  cat("Plotting retention time.")
  print(Med_Retention_time(evidence, plot = T))
  cat("Plotting ranked Intensities.")
  print(plot_rankedIntensities_V2(proteinGroups))
  # cat("Plotting m/z range.")
  # print(plotMZrange(evidence))  does not work for single samples
  # cat("Plotting overlapping proteins.")
  # print(plot_Overlapping_proteins(proteinGroups)) does not work for single samples
  # cat("Plotting PCA")
  # print(plot_QC_PCA(proteinGroups)) does not work for single samples

  #turn off PDF plotting
  dev.off()
  # dev.off()

  print("Finished.")
  return(dir)
}


# =============================================================================
# Quality control and exploratory plots
# =============================================================================

plot_quants <- function(quants,
                        label.names = T,
                        regex.rm = T,
                        scale.lim.y = NA){

  if(regex.rm){
    require(stringr)
    regex <- detect_prefix_regex(quants$Experiment)
    quants$Experiment <- str_remove(quants$Experiment, regex)
  }

  if("prot.quants" %in% colnames(quants)){
    ggplot(quants, aes(x=Experiment,
                       y=prot.quants))+
      geom_col()+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      {if(!label.names)theme(axis.text.x=element_blank(),
                             axis.ticks.x=element_blank())}+
      {if(!is.na(scale.lim.y))scale_y_continuous(expand = c(0, 0), limits = c(0, scale.lim.y))}+
      labs(x = "sample", y="protein quantifications", title = "Protein Quantifications")+
      {if(regex.rm)labs(subtitle = regex)}+
      geom_text(aes(label = prot.quants), vjust = -0.5)
  }else if("phos.quants" %in% colnames(quants)){
    ggplot(quants, aes(x=Experiment,
                       y=phos.quants))+
      geom_col()+
      theme_classic()+
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
      {if(!label.names)theme(axis.text.x=element_blank(),
                             axis.ticks.x=element_blank())}+
      {if(!is.na(scale.lim.y))scale_y_continuous(expand = c(0, 0), limits = c(0, scale.lim.y))}+
      labs(x = "sample", y="phosphopeptide quantifications", title = "Phosphopeptide Quantifications")+
      {if(regex.rm)labs(subtitle = regex)}+
      geom_text(aes(label = phos.quants), vjust = -0.5)
  }
}


#' #' Plot time-resolved protein or phosphosite intensities
#'
#' Creates time-course plots from MaxQuant-style protein or phosphosite
#' tables using an external sample annotation table.
#'
#' The function supports optional time rescaling to a 24-hour cycle and
#' different smoothing models for trend estimation.
#'
#' @param df A MaxQuant-style data frame containing quantitative values.
#' @param IdOfInterest Optional character vector of identifiers to subset
#'   (proteome: \code{id}, phosphoproteome: \code{UID}).
#' @param annotation Data frame containing sample annotation with at least
#'   the columns \code{Name}, \code{time}, and \code{group}.
#' @param scaleTo24h Logical. If \code{TRUE}, values larger than 24 are
#'   repeatedly reduced by 24 until all time points fall within a 24-hour
#'   cycle. This is useful when combining experiments with different
#'   sampling schemes.
#' @param fit Character string specifying the smoothing method to use.
#'   Supported values are \code{"harmonic"} and \code{"loess"}.
#' @param facet Character string controlling faceting of the plot.
#'   Use \code{"id"} to facet by identifier, or \code{"id+group"} to facet
#'   by identifier and experimental group.
#' @param datatype Character string specifying the data type.
#'   Either \code{"prot"} or \code{"phospho"} (requires columns
#'   \code{id} or \code{UID}, respectively).
#'
#' @return A \code{ggplot} object.
cyclic_plotting <- function(df, IdOfInterest=NULL, annotation, scaleTo24h = F, zScore = T, fit = "", facet = "", datatype){
  require(tidyverse)
  require(svDialogs)
  require(janitor)

  sc_cols <- make.names(annotation$Name)



  if(typeof(IdOfInterest)=="character"){
    IDs <- unlist(str_split(IdOfInterest, "\n")) # make an iterable value of GOIs
    # isolate ID rows from original table:
    if(grepl("UID", IDs[1], fixed = TRUE)){
      POI <- df %>%
        filter(Unique.identifier %in% IDs)
    }else{
      cat("Make sure your protein data contains a column called 'id'.")
      POI <-df %>%
        filter(id %in% IDs)
    }

    # pivot the table from wide to long so that ggplot can work with it
    if(datatype=="prot"){
      POI_long <- POI %>%
        select('Gene.names', "id", all_of(sc_cols)) %>%
        pivot_longer(!c('Gene.names',"id"), names_to = "Sample.Nr", values_to = "value")
    }else if(datatype == "phospho"){
      POI_long <- POI %>%
        select('Gene.names', "Unique.identifier", starts_with("Intensity.")) %>%
        pivot_longer(!c('Gene.names',"Unique.identifier"), names_to = "Sample.Nr", values_to = "value")
    }
  }else if (typeof(IdOfInterest) %in% c("numerical","double")){
    print("Please format your IdOfInterest as a character.")
    }else if (datatype == "phospho"){
    POI_long <- df %>%
      select('Gene.names', "Unique.identifier", starts_with("Intensity.")) %>%
      pivot_longer(!c('Gene.names',"Unique.identifier"), names_to = "Sample.Nr", values_to = "value")
  }else if (datatype == "prot"){
    cat("If you get the error #`cols` must select at least one column# it was assumed that we are dealing with proteomes but there is no columns starting with LFQ.intensity.")
    POI_long <- df %>%
      select('Gene.names', "id", starts_with("LFQ.intensity.")) %>%
      pivot_longer(!c('Gene.names',"id"), names_to = "Sample.Nr", values_to = "value")
  }

  # make syntactically valid names in case there were e.g. spaces
  annotation$Name <- make.names(annotation$Name)


  # add the time & group measure from the annotation file to POI_long
  if(!all(c("time","Name") %in% colnames(annotation))){
    print("Make sure that your df contains the column names 'Name' and 'time'.")
  }
  POI_long$time <- as.numeric(annotation$time[match(POI_long$Sample.Nr,annotation$Name)])
  if("group" %in% colnames(annotation)){
    POI_long$group <- annotation$group[match(POI_long$Sample.Nr,annotation$Name)]
  }

  # subtract 24 from times >24 until none remain >24
  if(scaleTo24h){
    while(max(POI_long$time)>24){
      POI_long <- POI_long %>% mutate(time = case_when(time>24 ~ time-24,
                                                       TRUE ~ time))
    }
  }

  if("Unique.identifier" %in% colnames(df)){
    POI_long$id <- POI_long$Unique.identifier
  }

  # for group wise z-scoring combine the UID and group column (if there is one)
  if(zScore){
    if("group" %in% colnames(annotation)){
      POI_long <- POI_long %>% mutate(id_group = paste(id,group, sep = "_"))
      POI_long$value <- ave(POI_long$value, POI_long$id_group, FUN=scale)
    }else{
      POI_long$value <- ave(POI_long$value, POI_long$id, FUN=scale)
      }
    }

  # plot all LFQ intensities
  if("group" %in% colnames(annotation)){
    p <- ggplot(POI_long, aes(x=time, y=value, color = group)) +
      geom_point() +
      {if(facet == "id+group")facet_grid(group~Gene.names+id)}+
      {if(facet == "id")facet_grid(.~Gene.names+id)}+
      {if(facet != "id" & facet != "id+group")facet_grid(. ~ .data[[facet]])}+
      {if(fit == "harmonic")geom_smooth(method = "lm", formula = y ~ cos(2*pi/24*x)+sin(2*pi/24*x), se=F)}+
      {if(fit == "loess")geom_smooth(se=F, span = .6)}+
      stat_summary(fun.data = mean_se,
                   geom = "errorbar"); p
  }else{
    p <- ggplot(POI_long, aes(x=time, y=value)) +
      geom_point() +
      {if(facet == "id")facet_grid(.~Gene.names+id)}+
      {if(fit == "harmonic")geom_smooth(method = "lm", formula = y ~ cos(2*pi/24*x)+sin(2*pi/24*x), se=F)}+
      {if(fit == "loess")geom_smooth(se=F, span = .6)}+
      stat_summary(fun.data = mean_se,
                   geom = "errorbar"); p
  }


  return(p)
}

enrichAnalPlot <- function(df = read_Perseus_file(),
                           Categorical_Column = "",
                           p_cutoff = 1,
                           FDR_cutoff = 1,
                           FC_min = 1,
                           FC_max = Inf,
                           intersection_cutoff = 1,
                           text_size = 7,
                           size_range = c(1,3),
                           top_n = Inf,
                           facet_wrap = T,
                           rm_label = F){
  # this function plots the results of a Perseus Fisher's Exact Test exported by "Generic matrix export"
  # Categorical_Column: per default it uses all Categories in the Category.column column; to change this use e.g. Categorical_Column = "Keywords"
  # p_cutoff: per default no p-value cutoff, you can set one
  # output: a list with the element "image" being the plot and "df" being the filtered dataframe used for the plot
  # example use: enrichAnalPlot()$image

  require(tidyverse)
  require(tidytext)
  
  

  if(all(Categorical_Column != "")){
    df <- df %>% filter(Category.column %in% Categorical_Column)
  }
  
  # Update "Category value" column
  df <- df %>%
    mutate(Category.value = if_else(Category.value == "+", Category.column, Category.value),
           Selection.value = if_else(Selection.value == "+", Selection.column, Selection.value))

  df <- df %>% 
    filter(Benj..Hoch..FDR<=FDR_cutoff) %>% 
    filter(P.value<=p_cutoff) %>% 
    filter(Intersection.size>=intersection_cutoff) %>% 
    filter(Enrichment.factor>FC_min & Enrichment.factor<FC_max)
  
  # Add the top_n parameter to filter the top n entries by Enrichment.factor
  df <- df %>%
    group_by(Selection.value) %>%
    arrange(desc(Enrichment.factor)) %>%
    slice_head(n = top_n) %>%
    ungroup()

  # Basic barplot
  p1<-ggplot(data = df,
             aes(x=reorder_within(Category.value, Enrichment.factor, Selection.value),
                 y=Enrichment.factor,
                 color=P.value,
                 size=Intersection.size)) +
    geom_point()+
    scale_size(limits = c(1, max(df$Intersection.size)),
                    breaks = round(seq(1, max(df$Intersection.size), length.out = 5)),
                    range = size_range)+
    scale_x_reordered()+
    scale_color_gradient(low="blue", high="red")+
    labs(x= Categorical_Column, y="Enrichment", color="p-value", size="count")+
    theme_bw()+
    theme(text = element_text(size = text_size),
          legend.key.size = unit(text_size, "pt"))+
    # geom_hline(aes(yintercept=1),
    #            linetype="dashed",size=1, color="red")+
    ylim(ifelse(min(df$Enrichment.factor) > 1, 1, min(df$Enrichment.factor)), max(df$Enrichment.factor)) +
    {if(facet_wrap)facet_wrap(~Selection.value, scales = "free")}+
    {if(rm_label)theme(
      strip.background = element_blank(),
      strip.text.x = element_blank())}+
    coord_flip();p1

  return(list("image" = p1, "df" = df))
}

plot_missed_cleavages <- function(evidence = read_Perseus_file(),
                                  rm.regex=T){
  require(stringr)
  require(dplyr)
  require(ggplot2)
  
  regex <- detect_prefix_regex(evidence$Experiment)
  
  if(rm.regex){
    evidence$Experiment <- str_remove(evidence$Experiment, regex)
  }

  df <- evidence %>% 
    group_by(Experiment, Missed.cleavages) %>%
    dplyr::summarise(n = n(), .groups = "drop") %>%
    group_by(Experiment) %>% 
    mutate(freq = n / sum(n))
  ggplot(df, 
         aes(x=Experiment, 
             y=freq, 
             color=as.factor(Missed.cleavages), 
             group=as.factor(Missed.cleavages)))+
    geom_line(size=1) + 
    theme(axis.text.x = element_text(angle = 90))+
    labs(color = "Missed cleavage", 
         subtitle = regex, 
         title = "Missed cleavages")
}

Med_Retention_time <- 
  function(evidence = read_Perseus_file(),
                               min_Retention = 20, max_Retention = 90,
                               plot = T){
  # calculate the median retention time of everthing eluting in a certain Retention time window
  # from Maxquant output use the evidence file, from DIA-NN output use the general report file
  
  # Load required packages
  require(tidyverse)
  require(progress)
  
  # Read evidence file if a character string is provided
  if(all(class(evidence) == "character")){
    evidence <- read_Perseus_file(evidence)
  }
    
  # Handle missing Retention.time column, which is a column in MaxQuant output and not present in DIA-NN output
  if(!"Retention.time" %in% colnames(evidence)){
    evidence <- evidence %>% mutate(Retention.time = RT,
                                    Retention.length = RT.Stop-RT.Start,
                                    Experiment = Run)
  }
  
  # Filter evidence within the specified retention time window
    evidence_filtered <- evidence %>% filter(Retention.time>min_Retention & Retention.time<max_Retention) %>% mutate(Retention.length.min = Retention.length*60)
  df <- evidence_filtered %>%
    group_by(Experiment) %>%
    dplyr::summarize(median_second = median(Retention.length.min, na.rm = T)) %>% 
    left_join(evidence %>%
                group_by(Experiment) %>%
                dplyr::summarise(measured_RT_range_min = paste0(round(min(Retention.time),1), "-", round(max(Retention.time),1))),
              by = "Experiment")
  
  cat(paste0("mRT calculated within the range of ",min_Retention,"-",max_Retention," min"))
  
  if(!plot){return(df)} else{
    require(stringr)
    regex <- detect_prefix_regex(df$Experiment)
    df$Experiment <- str_remove(df$Experiment, regex)
    evidence_filtered$Experiment <- str_remove(evidence_filtered$Experiment, regex)

    p <- ggplot(df, aes(x=Experiment, y=median_second, group=1))+
      geom_line() +
      theme(axis.text.x = element_text(angle = 90)) +
      labs(y = "median Retention time [s]",
           title = paste0("Median Retention time (", as.character(min_Retention), "-", as.character(max_Retention), "min)"),
           subtitle = regex) +
      expand_limits(y = 0)
    p2 <- ggplot(evidence_filtered, aes(x=Retention.length.min))+
      geom_histogram() +
      facet_wrap(vars(Experiment)) +
      theme(axis.text.x = element_text(angle = 90)) +
      labs(title = paste0("Median Retention times (", as.character(min_Retention), "-", as.character(max_Retention), "min)"),
           subtitle = regex)
    plot(p2)
    plot(p)
    cat("Two figures have been plotted, retention time histograms and medRT profile. Use the arrows in the plot pane to navigate them.\n")
    return(df)
  }
}

plot_rankedIntensities_V1 <- function(proteinGroups = read_Perseus_file(),
                                   rm.regex = T){
  # this function will mix all your samples together to an overall ranking
  # it will also label the top and low 5 gene names
  # the mixing probably only makes sense with replicates (if even)
  df <- select_intensity_columns(proteinGroups, addCols = c("Gene.names")) %>%
    tibble::rownames_to_column("number") %>%
    pivot_longer(!c('Gene.names','number'), names_to = "Experiment", values_to = "Intensity") %>%
    filter(Intensity != 0) %>%
    # mutate(labelthis = ifelse(Intensity==max(Intensity)|Intensity==min(Intensity), Gene.names, NA)) %>%
    tibble::rownames_to_column("id")
  df <- df %>%
    left_join(df %>% arrange(desc(Intensity)) %>% distinct(Gene.names,.keep_all = TRUE) %>%
                           head(5) %>%
                           mutate(labelthis = Gene.names)) %>%
    left_join(df %>% arrange(Intensity) %>% distinct(Gene.names,.keep_all = TRUE) %>%
                head(5) %>%
                mutate(labelthis = Gene.names),by = c("id", "number", "Gene.names", "Experiment", "Intensity")) %>%
    mutate(labelthis = ifelse(is.na(labelthis.x),labelthis.y,labelthis.x))
  require(stringr)
  regex <- detect_prefix_regex(df$Experiment)
  if(rm.regex){
    df$Experiment <- str_remove(df$Experiment, regex)
  }
  require(tidytext)
  ggplot(df, aes(x=reorder_within(number, Intensity, Experiment), y=Intensity, color=Experiment, group=as.factor(Experiment)))+
    geom_point()+
    scale_y_continuous(trans='log10')+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(x = "rank", title = "ranked Intensities", subtitle = regex)+
    {if(nlevels(as.factor("Experiment")))theme(legend.position="none")}+
    ggrepel::geom_label_repel(aes(label = labelthis), min.segment.length = 0, show.legend = FALSE)
}

plot_rankedIntensities_V2 <- function(proteinGroups = read_Perseus_file(),
                                      rm.regex = T,
                                      subgroup = NA,
                                      GOI = NA,
                                      alpha = 1){
  # this function will calculate a rank percentile for all proteins within the Experiments
  # the ranked intensity curves will be colored by Experiment
  require(tidyverse)
  try(proteinGroups <- proteinGroups %>% rename(Genes = Gene.name), silent = T)
  if(!is.na(subgroup)){
    df <- select_intensity_columns(proteinGroups, addCols = c("Genes",subgroup)) %>%
      pivot_longer(-c('Genes', all_of(subgroup)), names_to = "Experiment", values_to = "Intensity") %>%
      filter(Intensity != 0 & !is.nan(Intensity)) %>%
      arrange(Experiment, Intensity) %>%
      group_by(Experiment) %>%
      mutate(rank = rank(Intensity),
             percrank=rank(Intensity)/length(Intensity))
  }else{
    df <- select_intensity_columns(proteinGroups, addCols = c("Genes")) %>%
      pivot_longer(-'Genes', names_to = "Experiment", values_to = "Intensity") %>%
      filter(Intensity != 0 & !is.nan(Intensity)) %>%
      arrange(Experiment, Intensity) %>%
      group_by(Experiment) %>%
      mutate(rank = rank(Intensity),
             percrank=rank(Intensity)/length(Intensity))
  }
  

  require(stringr)
  regex <- detect_prefix_regex(df$Experiment)
  if(rm.regex){
    df$Experiment <- str_remove(df$Experiment, regex)
  }
  require(tidytext)

  ggplot(df,
         aes(x=rank, y=Intensity, color = Experiment))+
    geom_point(alpha = alpha)+
    scale_y_continuous(trans='log10')+
    theme_classic()+
    {if(!is.na(subgroup))geom_point(data = df %>% filter(df[subgroup] == "+"), color = "blue", alpha = alpha)}+
    labs(x = "intensity rank", title = "ranked Intensities", subtitle = regex)
  
  ggplot(df,
         aes(x=rank, y=Intensity))+
    geom_point(alpha = alpha)+
    facet_wrap(vars(Experiment))+
    scale_y_continuous(trans='log10')+
    theme_classic()+
    {if(!is.na(GOI))geom_point(data = df %>% filter(grepl(GOI, Genes)), color = "blue", alpha = 1)}+
    {if(!is.na(GOI))ggrepel::geom_text_repel(data = df %>% filter(grepl(GOI, Genes)),
                                              aes(label = Genes), min.segment.length = 0)}+
    labs(x = "intensity rank", title = "ranked Intensities", subtitle = regex)
}

plot_QC_PCA <- function(proteinGroups = read_Perseus_file(),
                        groups = NA){
  # this script is supposed to take the proteinGroups.txt output from MaxQuant and potentially an annotation file
  # it will log2 transform and only use rows with 100% valid values

  int <- select_intensity_columns(proteinGroups)
  log2_df <- log2(int)
  is.na(log2_df) <- sapply(log2_df, is.infinite) # the log2 mutation makes 0 values to "-inf", I want them as NA

  df_NoNa <- log2_df[complete.cases(log2_df), ]
  pca_ma <- prcomp(t(df_NoNa), scale. = FALSE)

  varExp <- round(pca_ma$sdev^2/sum(pca_ma$sdev^2)*100, 1)
  df_pca <- data.frame(PC1 = pca_ma$x[,1], PC2 = pca_ma$x[,2]) %>%
    tibble::rownames_to_column(var = "Name")

  require(stringr)
  regex <- detect_prefix_regex(df_pca$Name)
  df_pca$Name <- str_remove(df_pca$Name, regex)

  if(!is.na(groups)){
    groups$Name <- make.names(groups$Name)
    df_pca <- df_pca %>%
      tibble::rownames_to_column("Name") %>%
      left_join(groups)
  }
  ggplot(df_pca, aes(PC1, PC2)) + ggrepel::geom_label_repel(aes(label=Name), size=3) +
    geom_point() +
    theme_bw(base_size=10) +   scale_color_brewer(name = "Group", palette="Dark2") +
    theme(aspect.ratio=1) + xlab(paste0("PC1, VarExp: ", varExp[1], "%")) +
    ylab(paste0("PC2, VarExp: ", varExp[2], "%")) +
    labs(title = "Principal component analysis", subtitle = regex)
}

plot_circular_histogram <- function(phases, color = "none", color_12 = "#B788BD", color_0 = "#F16130", bar_line_width = 1){
  require(cowplot)
  # takes a list of values and plots them in a circular histogram
  phases<-phases[!is.na(phases)]
  # create a df with the amount of entries with phases adjacent to each full hour of the day
  DF <- data.frame(
    variable=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23),
    value=c(sum(phases>23.5 & phases<0.5), sum(phases>0.5 & phases<1.5), sum(phases>1.5 & phases<2.5),sum(phases>2.5 & phases<3.5),
            sum(phases>3.5 & phases<4.5),sum(phases>4.5 & phases<5.5),sum(phases>5.5 & phases<6.5),sum(phases>6.5 & phases<7.5),
            sum(phases>7.5 & phases<8.5),sum(phases>8.5 & phases<9.5),sum(phases>9.5 & phases<10.5),sum(phases>10.5 & phases<11.5),
            sum(phases>11.5 & phases<12.5),sum(phases>12.5 & phases<13.5),sum(phases>13.5 & phases<14.5),sum(phases>14.5 & phases<15.5),
            sum(phases>15.5 & phases<16.5),sum(phases>16.5 & phases<17.5),sum(phases>17.5 & phases<18.5), sum(phases>18.5 & phases<19.5),
            sum(phases>19.5 & phases<20.5),sum(phases>20.5 & phases<21.5),sum(phases>21.5&phases<22.5),sum(phases>22.5&phases<23.5)))

  DF <- DF %>% mutate(day_color = case_when(color != "none" ~ color,
                                            color == "none" & variable<18 & variable >6  ~ color_12,
                                            color == "none" & (variable>18 | variable <6) ~ color_0)) # Phase 6 & 18 get grey now, what should they have?

  # make a bar plot to be able to extract the y-axis limit (important for making the rose plot more pretty)
  bar <- ggplot(DF, aes(x=as.factor(variable), y=value))  +
    geom_bar(width = 1, fill = "blue", colour ="black", stat="identity", alpha=0.7)

  # extract the y-axis limit
  roseLimits <- layer_scales(bar)$y$range$range

  # make the same plot without theme then plot the line grid in a standardized and pretty way, then the bar plot again on top and make it circular
  p <- ggplot(DF, aes(x=as.factor(variable), y=value,fill = day_color)) + theme_bw() +
    geom_bar(width = 1, stat="identity", size = bar_line_width) +
    geom_hline(yintercept = seq(0, roseLimits[2], by = roseLimits[2]/5), colour = "grey60", size = 0.2) +
    geom_vline(xintercept = seq(0, 24, by = 1), colour = "grey60", size = 0.2) +
    geom_bar(width = 1, colour ="black",fill = DF$day_color, stat="identity", size = bar_line_width) +
    #scale_y_continuous(breaks = c(0,50,100,150,200,250,300)) +
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(color="black", size = 10),
          axis.ticks = element_line(color = "black"),
          axis.text.y = element_text(color="black", size = 10)) +
    coord_polar(start = - pi / 24) +
    labs(x = "", y = "")

  ggdraw(p) +
    draw_label(paste0("n = ", as.character(length(phases))), x = 0.95, y = 0.1, hjust = 1, vjust = 0, size = 10)

}

plotMZrange <- function(evidence = read_Perseus_file()){
  quantiles <- quantile(evidence$m.z, probs = c(0.025,0.975))
  p <- ggplot(evidence, aes(x=m.z))+
    geom_histogram(fill = "blue", bins = 100)+
    theme_classic();p
  y.Range <- ggplot_build(p)$layout$panel_params[[1]]$y.range
  p <- ggplot(evidence, aes(x=m.z))+
    annotate("rect", xmin = quantiles[1], xmax = quantiles[2],
             ymin = 0, ymax = y.Range[2],
             fill="light blue") +
    geom_histogram(fill = "blue", bins = 100)+
    geom_text(y = y.Range[2], vjust=1, x = quantiles[2], hjust=0,
              label=paste0("95% interval \n", round(quantiles[1]), "-", round(quantiles[2]), " m.z"))+
    theme_classic()+
    labs(title = "Global m/z distribution");p

  return(p)
}

plot_Overlapping_proteins <- function(proteinGroups = read_Perseus_file()){
  require(UpSetR)
  int <- select_intensity_columns(proteinGroups)
  int[int != 0] <- 1
  maxSetSize <- max(colSums(int))
  upset(int, order.by = "freq", set_size.show = T, set_size.scale_max = maxSetSize*1.2)
}















# =============================================================================
# Data preprocessing and transformation
# =============================================================================

' Read Perseus-style tables (including annotation rows)
#'
#' Supports standard Perseus text files and parquet files
#' (e.g. from DIA-NN exports). Annotation rows starting with "#!"
#' are automatically removed.
#'
#' @param path Path to the file
#' @param changeWD Set working directory to the file location
#' @return A data.frame
read_Perseus_file <- function(path = NULL, changeWD = TRUE) {
  cat("The read_Perseus_file function was updated on 12.6.24 to run faster and support .parquet files which are much faster to load than the report.tsv from DIA-NN. Tell Fabian if you encounter errors.")
  require(tidyverse)
  require(stringr)
  require(arrow)
  
  if (is.null(path)) {
    path <- file.choose()
  }
  
  if (changeWD) {
    setwd(dirname(path))
  }
  
  print(path)
  
  if(grepl("\\.parquet$", path)){
    df <- read_parquet(file = path, 
                       col_select = NULL, as_data_frame = TRUE, props = ParquetArrowReaderProperties$create(), mmap = TRUE) # this creates a tbl, not da df; if that ever becomes a problem I could use df <- as.data.frame(df)
  }else{
    # Read all lines of the file
    all_lines <- readLines(path)
    
    # Identify annotation rows 
    annot_rows <- grep("^#!", all_lines)
    
    # Check if annotation rows exist
    if (length(annot_rows) == 0) {
      # If no annotation rows found, read the entire file
      df <- read.delim(path, header = TRUE)
    } else {
      # Extract non-annotation rows
      data_lines <- all_lines[-annot_rows]
      
      # Create a connection with the remaining lines
      con <- textConnection(data_lines)
      
      # Read the data into a data frame
      df <- read.delim(con, header = TRUE)
      
      # Close the connection
      close(con)}
  }
  
  return(df)
}

read_Perseus_file_archived_V1 <- function(path = NULL, changeWD = T){
  require(tidyverse)
  require(stringr)
  if(is.null(path)){
    path <- file.choose()
  }
  if(changeWD==T){
    setwd(dirname(path))
  }
  print(path)
  df <- read.delim(path) # this contains annotation rows were entries start with #! in the first column and its entries make all columns parse as char
  # count annotation rows
  annotRowNr <- str_count(paste(df[,1], collapse = ""), "#!")
  header <- colnames(df)
  # read file again, skipping annotation rows
  df <- read.delim(path, header = T,skip = annotRowNr)
  colnames(df) <- header
  return(df)
}

extract_UniMods <- function(pr_matrix = read_Perseus_file()) {
  pr_matrix <- pr_matrix %>%
    mutate(
      # Extracting numbers following "UniMod:" and splitting them into a list
      UniMod_numbers = str_extract_all(Modified.Sequence, "(?<=UniMod:)[0-9]+")
    )
  return(pr_matrix)
}

filter_data <- function(df,dataType){
  # remove potential contaminants, reverse & only identified by site
  df <- df %>%
    filter(Potential.contaminant != "+") %>%
    filter(Reverse != "+")
  if(dataType == "prot"){
    df <- df %>%
      filter(Only.identified.by.site != "+")
  }
  return(df)
}

expand_site_table <- function(data) {
  expanded_data <- data %>%
    # gather all Intensity___ columns into one sample column
    gather(sample, Intensity, contains("___"),-starts_with("Intensity___")) %>%
    # separate this new sample column, to get a multiplicity column
    separate(sample, into = c("expanded", "multiplicity"), sep="___") %>%
    # spread the "expanded" column again into indiv sample columns, now having 3x the entries after isolating the multiplicity into its own column
    spread(expanded,Intensity, sep="_")
  return(expanded_data)
}

remove_RegEx_ColName <- function(df, regularExpression){
  # remove expressions from column names
  for(expression in regularExpression){
    colnames(df)<-gsub(expression,"",colnames(df))
  }
  return(df)
}

select_intensity_columns <- function(df, addCols = c()){
  require(tidyverse)
  if(any(grepl("LFQ.intensity.", colnames(df)))){
    df <- df %>% dplyr::select(contains("LFQ.intensity."), all_of(addCols))
  }else if(any(grepl("Intensity.", colnames(df)))){
    df <- df %>% dplyr::select(contains("Intensity."), all_of(addCols))
  }else{
    df <- df %>% dplyr::select(contains(".raw"), all_of(addCols))
  }
  return(df)
}

detect_prefix_regex <- function(str_list){
  # this function takes the first string from a list of strings and truncates its end until it finds a leading regex that matches all strings in list
  require(stringr)

  first_str <- str_list[1]

  for (i in 1:nchar(first_str)) {
    truncated_str <- str_sub(first_str, end=-i)
    if(all(str_detect(str_list, truncated_str))==T) break
  }
  return(truncated_str)
}

profilePlot <- function(df, annot, profile.id){
  # df is a Perseus matrix where replicates of samples at one circadian time have been averaged
  # it also contains a column "Selection" where a "+" indicates correlating profiles to be colored
  # annot is a data frame with the columns column.names and time
  # column.names is a list of the names of the columns of the averaged replicates in df
  # time is a list of the circadian time in the same order as the averaged replicates in df
  # profile is the id of the entry that was used for the profile correlation

  # plots a line plot of the profiles


  df_long <- df %>% pivot_longer(cols = all_of(annot$column.names), names_to = "name", values_to = "value") %>% right_join(annot, by=c("name"="column.names"))

  ggplot(df_long %>% filter(Selection == ""), aes(time,value,group=id))+
    geom_line(color = "lightgrey")+
    geom_line(data = df_long %>% filter(Selection == "+"), color = "darkorchid1")+
    geom_line(data = df_long %>% filter(id == profile.id), color = "darkorchid4", size = 3)+
    theme_classic()

}


pearson_dist <- function(data, profile.id, plot.cutoff) {

  # transpose df
  df_transpose <- df %>%
    select(all_of(annot$column.names), id) %>%
    gather(variable, value, -id) %>%
    spread(id, value) %>%
    select(-variable)

  # loop through columns, apply `cor` vs 'd' column
  cor <- colnames(df_transpose) %>%
    set_names() %>%
    map(~ cor(df_transpose[, .x], df_transpose[, '1792'], use = "complete.obs")) %>%
    map_dfr(., broom::tidy, .id = "id")
  cor$id <- as.integer(cor$id)
  cor <- cor %>% mutate(pearson_dist = 1-x,
                        x = NULL)

  df <- df %>% left_join(cor, by=c("id"="id"))

  df_long <- df %>% pivot_longer(cols = all_of(annot$column.names), names_to = "name", values_to = "value") %>% right_join(annot, by=c("name"="column.names"))


  print(
    ggplot(df_long %>% filter(Selection == ""), aes(time,value,group=id))+
      geom_line(color = "lightgrey")+
      geom_line(data = df_long %>% filter(pearson_dist<plot.cutoff), color = "darkorchid1")+
      geom_line(data = df_long %>% filter(id == profile.id), color = "darkorchid4", size = 3)+
      theme_classic()
  )


  # return the resulting distance vector
  return(df)
}


category_counting <- function (df, column_name, distinction_col = NA){
  if (!is.na(distinction_col)){df <- df %>% distinct(!!as.name(distinction_col), .keep_all = TRUE)}
  count <- data.frame(name = str_split(paste(df[,column_name], collapse = ";"), # first collapse a column from the df
                                       ";", simplify = T) %>% t) %>% # then split it along ";" and make a df out of it
    add_count(name) %>%
    distinct %>%
    arrange(desc(n))
  return(count)
}

Perseus_left_join_indexing <- function(df1,df2,
                                       join_style = "inner",
                                       unnest_column1 = "Protein.Group",
                                       unnest_column2 = "Protein.Group"){
  # writing a function to imitate the Perseus matching
  # not sure if this works
  
  require(dplyr)
  
  # Sample data
  df1 <- df1 %>% tibble::rowid_to_column("ID_1")
  df2 <- df2 %>% tibble::rowid_to_column("ID_2")
  
  # Function to split and unnest the DataFrame
  split_and_unnest <- function(df, col) {
    df %>%
      mutate(!!col := strsplit(as.character(!!sym(col)), ";")) %>%
      unnest(!!sym(col))
  }
  
  # Split and unnest both dataframes
  df1_unnested <- split_and_unnest(df1, unnest_column1)
  df2_unnested <- split_and_unnest(df2, unnest_column2)
  
  if(join_style == "left.index"){
    # Perform a left join style indexing
    result <- df1_unnested %>% mutate(match = case_when(Protein.Group %in% df2_unnested$Protein.Group ~ "+",
                                                        T ~ ""))
  }else if(join_style == "outer"){
    # Perform the inner join
    result <- outer_join(df1_unnested, df2_unnested, by = c("Protein.Group" = "Protein.Group"), keep = T)
  }else{
    cat("unsupported joining style")
    break
  }
  
  # Aggregate the results back to the original format
  aggregated_result <- result %>%
    group_by(ID_1) %>% 
    mutate(Protein.Group = paste0(Protein.Group, collapse = ";"),
           match = if(any(match == "+")) "+" else match) %>% 
    ungroup %>% 
    distinct(ID_1, .keep_all = TRUE)
  
  return(aggregated_result)
}

Perseus_matching <- function(df1,df2,
                             join_style = "left",
                             join_columns = c("Genes","Sequence.window")){
  # writing a function to imitate the Perseus matching, by splitting matching columns at ";" separators
  # be aware that this is not a perfect solution for matching Phospho datasets, neither is Perseus
  
  require(tidyverse)
  
  # Sample data
  df1 <- df1 %>% tibble::rowid_to_column("ID_1")
  df2 <- df2 %>% tibble::rowid_to_column("ID_2")
  
  # Function to split and unnest multiple columns
  split_and_unnest <- function(df, col_list) {
    # Iterate over the list of column names
    for (col in col_list) {
      df <- df %>%
        mutate(!!sym(col) := strsplit(as.character(.data[[col]]), ";")) %>%
        unnest(!!sym(col))
    }
    return(df)
  }
  
  
  # Split and unnest both dataframes
  df1_unnested <- split_and_unnest(df1, join_columns)
  df2_unnested <- split_and_unnest(df2, join_columns)
  
  if(join_style == "left"){
    # Perform a left join style indexing
    result <- df1_unnested %>% left_join(df2_unnested,
                                         by = join_columns)
  }else{
    cat("unsupported joining style")
    return()
  }
  
  # Aggregate the results back to the original format
  aggregated_result <- result %>%
    group_by(ID_1) %>% 
    mutate(across(all_of(join_columns), ~ paste(unique(.), collapse = ";"))) %>% 
    ungroup %>% 
    distinct(ID_1, .keep_all = TRUE) %>% 
    select(-ID_1, -ID_2)
  
  return(aggregated_result)
}

# =============================================================================
# Imputation
# =============================================================================

set_floor_min_if_all_group_nan <- function(report, metadata, group_col) {
  
  # Ensure that the metadata has both 'Name' and the grouping column (e.g., background.condition)
  if (!("Name" %in% colnames(metadata)) || !(group_col %in% colnames(metadata))) {
    stop("The 'metadata' must contain both 'Name' and the specified grouping column.")
  }
  
  # Find the intensity columns in report based on metadata$Name
  intensity_cols <- metadata$Name
  
  # Select only the columns from report that are in the intensity_cols
  report_selected <- report[, intensity_cols, drop = FALSE]
  
  # Group the intensity columns by the group_col in metadata
  grouped_cols <- split(metadata$Name, metadata[[group_col]])
  
  # Calculate the minimum value across all columns in the dataframe, ignoring NaNs
  min_value <- min(as.matrix(report_selected), na.rm = TRUE)
  floor_min_value <- floor(min_value)  # Take the floor of the minimum value
  
  # Iterate over each group
  for (group in names(grouped_cols)) {
    group_columns <- grouped_cols[[group]]  # Get the columns in the current group
    
    # Check if all values in these columns are NaN for each row
    na_rows <- apply(report_selected[, group_columns, drop = FALSE], 1, function(row) all(is.nan(row)))
    
    # Set the NaN values to the floor of the minimum value for the rows where all values are NaN
    report_selected[na_rows, group_columns] <- floor_min_value
  }
  
  # Update the original report with the modified values
  report[, intensity_cols] <- report_selected
  
  # Return the modified report
  return(report)
}

impute_perseus_normal <- function(df, width = 0.3, downshift = 1.8) {
  # Iterate over each column in the dataframe
  for (col in colnames(df)) {
    # Identify missing values (e.g., NaN or NA)
    missing_indices <- is.na(df[[col]])
    
    # If there are missing values, perform imputation
    if (any(missing_indices)) {
      # Calculate the mean and standard deviation of non-missing values
      non_missing_values <- df[[col]][!missing_indices]
      mean_val <- mean(non_missing_values, na.rm = TRUE)
      sd_val <- sd(non_missing_values, na.rm = TRUE)
      
      # Calculate the mean and standard deviation for the normal distribution
      impute_mean <- mean_val - downshift * sd_val
      impute_sd <- width * sd_val
      
      # Impute missing values with random draws from the normal distribution
      imputed_values <- rnorm(sum(missing_indices), mean = impute_mean, sd = impute_sd)
      
      # Replace missing values with the imputed values
      df[[col]][missing_indices] <- imputed_values
    }
  }
  
  # Return the dataframe with imputed values
  return(df)
}

row_average_imputation <- function(data_matrix, groups_matrix, group){
  # groups_matrix has to contain the columns "Name" and another for the groups you want to average
  # the column name for the groups you want to average you can define with the group variable
  
  # Tidy up names
  groups_matrix$Name <- make.names(groups_matrix$Name)
  
  # Create a list of lists where each list contains the names from "Name" column for each group
  numerical_column_groups <- split(groups_matrix$Name, groups_matrix[[group]])
  
  
  # Function to replace missing values with row averages for specified columns
  replace_missing_with_row_average <- function(data, columns) {
    imputed_data <- data %>%
      mutate(across(all_of(columns), ~ifelse(is.na(.), rowMeans(across(all_of(columns)), na.rm = TRUE), .)))
    return(imputed_data)
  }
  
  # Iterate over each list of columns and apply the transformation
  imputed_matrix <- data_matrix
  for (columns_list in numerical_column_groups) {
    imputed_matrix <- replace_missing_with_row_average(imputed_matrix, columns_list)
  }
  
  # Replace any NaN values with NA
  imputed_matrix <- imputed_matrix %>% mutate_all(~ifelse(is.nan(.), NA, .))
  return(imputed_matrix)
}

# =============================================================================
# Statistical testing
# =============================================================================

# conditions for t_test() are extracted by dplyr::pull() from a design table and contain names of columns from the same experimental group, e.g.:
# WT_0 <- metadata_phos %>% filter(condition == "ZT0" & background == "Camk2b (WT)") %>% pull(Name)

#' Row-wise two-sample t-tests for proteomics-style intensity tables
#'
#' Performs row-wise t-tests between two sets of sample columns and
#' returns effect sizes, p-values and Benjamini–Hochberg adjusted p-values.
#'
#' @param report Data frame containing quantitative intensities
#' @param condition_1 Character vector of column names
#' @param condition_2 Character vector of column names
#' @param min_valid Minimum number of non-missing values per group
#' @param apply_imputation Logical, apply Perseus-style imputation before testing
#' @return Original table with additional result columns
t_test <- function(report, condition_1, condition_2, min_valid = 3, apply_imputation = FALSE){
  # internal function for individual t-test
  t_test_protein <- function(row, condition_1, condition_2) {
    group1 <- as.numeric(row[condition_1])
    group2 <- as.numeric(row[condition_2])
    
    # Check if both groups have at least 2 non-NA values for the t-test
    if (sum(!is.na(group1)) >= 2 & sum(!is.na(group2)) >= 2) {
      
      # Check if the groups have identical values
      if (all(group1 == group1[1], na.rm = TRUE) & all(group2 == group2[1], na.rm = TRUE)) {
        # If both groups are constant (or identical), return NA
        return(c(NA, NA, NA))
      }
      
      # Perform the t-test if the values are not constant
      test <- t.test(group1, group2, var.equal = TRUE)
      return(c(test$estimate[1], test$estimate[2], test$p.value))
    } else {
      # Return NaN for mean values and p-value if not enough observations
      return(c(NA, NA, NA))
    }
  }
  
  # Apply the t-test across all proteins
  # Adding columns for the mean of each group and p-value
  
  # Capture the names of condition_1 and condition_2 as strings
  condition_1_name <- deparse(substitute(condition_1))
  condition_2_name <- deparse(substitute(condition_2))
  
  # Filter for rows where at least one of the groups has `min_valid` non-NA values
  filtered_report <- report %>%
    select(-matches("T_test_Difference_"), -matches("T_test_Pval_"), -matches("neg_log10_T_test_Pval_"), -matches("Adj_Pval_")) %>% 
    rowwise() %>%
    filter(
      sum(!is.na(c_across(all_of(condition_1)))) >= min_valid |
        sum(!is.na(c_across(all_of(condition_2)))) >= min_valid
    ) %>%
    ungroup()
  cat("Number of rows after filtering: ", as.character(nrow(filtered_report)), "\n")
  
  # Apply imputation if specified
  if (apply_imputation) {
    filtered_report <- impute_perseus_normal(filtered_report)
  }
  
  # Compute the t-test for rows that pass the filter
  results <- filtered_report %>%
    rowwise() %>%
    mutate(
      T_test_Difference = mean(c_across(all_of(condition_1)), na.rm = TRUE) - mean(c_across(all_of(condition_2)), na.rm = TRUE),
      T_test_Pval = t_test_protein(cur_data(),condition_1,condition_2)[3],
      neg_log10_T_test_Pval = -log(T_test_Pval, 10)
    ) %>%
    ungroup()
  
  
  # 5. Adjust p-values for multiple testing using Benjamini-Hochberg method (FDR)
  results <- results %>%
    mutate(Adj_Pval = p.adjust(T_test_Pval, method = "BH")) %>%
    # Rename columns using dynamic names
    rename(
      !!paste0("T_test_Difference_", condition_1_name, "_vs_", condition_2_name) := T_test_Difference,
      !!paste0("T_test_Pval_", condition_1_name, "_vs_", condition_2_name) := T_test_Pval,
      !!paste0("neg_log10_T_test_Pval_", condition_1_name, "_vs_", condition_2_name) := neg_log10_T_test_Pval,
      !!paste0("Adj_Pval_", condition_1_name, "_vs_", condition_2_name) := Adj_Pval
    )
  
  # Merge the original report_pr with results to retain all rows, filling NAs for rows that didn't pass the filter
  final_results <- left_join(report, results %>%
                               select(ID, matches("T_test_Difference_"), matches("T_test_Pval_"), matches("neg_log10_T_test_Pval_"), matches("Adj_Pval_")), 
                             by = "ID")
  
  return(final_results)
}

check_exclusive_presence <- function(report, condition_1, condition_2, min_valid = 3) {
  # Capture the names of condition_1 and condition_2 as strings
  condition_1_name <- deparse(substitute(condition_1))
  condition_2_name <- deparse(substitute(condition_2))
  
  # Check for exclusive presence and create dynamically named columns
  results <- report %>%
    rowwise() %>%
    mutate(
      # Check for exclusive presence in condition_1 (all values are NA in condition_2, but not in condition_1)
      !!paste0("exclusively_present_in_", condition_1_name, "_vs_", condition_2_name) := 
        (sum(!is.na(c_across(all_of(condition_2)))) == 0 & sum(!is.na(c_across(all_of(condition_1)))) >= min_valid),
      
      # Check for exclusive presence in condition_2 (all values are NA in condition_1, but not in condition_2)
      !!paste0("exclusively_present_in_", condition_2_name, "_vs_", condition_1_name) := 
        (sum(!is.na(c_across(all_of(condition_1)))) == 0 & sum(!is.na(c_across(all_of(condition_2)))) >= min_valid)
    ) %>%
    ungroup()
  
  return(results)
}


# Function to plot T-test results (only works on dataframes produced by t_test() using the same conditions)
plot_t_test_results <- function(results, condition_1, condition_2,
                                show_cutoff = T, show_n = T, label = "Genes") {
  # Capture the names of condition_1 and condition_2 as strings
  condition_1_name <- deparse(substitute(condition_1))
  condition_2_name <- deparse(substitute(condition_2))
  
  # Create dynamic column names based on conditions
  T_test_Difference_col <- paste0("T_test_Difference_", condition_1_name, "_vs_", condition_2_name)
  T_test_Pval_col <- paste0("T_test_Pval_", condition_1_name, "_vs_", condition_2_name)
  neg_log10_T_test_Pval_col <- paste0("neg_log10_T_test_Pval_", condition_1_name, "_vs_", condition_2_name)
  Adj_Pval_col <- paste0("Adj_Pval_", condition_1_name, "_vs_", condition_2_name)
  
  lowest_neglog_pValue_at_qCutoff <- results %>% filter(!!sym(Adj_Pval_col)<0.05) %>% select(!!sym(neg_log10_T_test_Pval_col)) %>% min(na.rm = T)
  cat("Minimum p-value at adjusted p-value < 0.05: ", as.character(round(10^(-lowest_neglog_pValue_at_qCutoff),4)), " for ", condition_1_name, " vs ", condition_2_name, "\n")
  cutoff <- data.frame(yintercept=lowest_neglog_pValue_at_qCutoff, Lines='adj. p-value < 0.05')
  
  # Print dynamic column names for debugging
  cat("T_test_Difference column:", T_test_Difference_col, "\n")
  cat("neg_log10_T_test_Pval column:", neg_log10_T_test_Pval_col, "\n")
  
  # Count the number of non-NA points plotted
  number_of_datapoints <- results %>%
    filter(!is.na(!!sym(T_test_Difference_col)) & !is.na(!!sym(neg_log10_T_test_Pval_col))) %>%
    nrow()
  
  # Generate the plot
  p <- ggplot(results,
              aes(
                x = !!sym(T_test_Difference_col), 
                y = !!sym(neg_log10_T_test_Pval_col), 
                text = !!sym(label),
                label = !!sym(label)
              )) +
    geom_point(color = "lightcyan3") +
    {if(show_cutoff)geom_hline(aes(yintercept=yintercept, linetype=Lines), cutoff)}+
    {if(show_n)annotate(geom = 'text', label = paste0("n = ", number_of_datapoints), x = Inf, y = -Inf, hjust = 1, vjust = -1)}+
    labs(y="-log10 p-value", x=paste0("Difference ", condition_1_name, "-", condition_2_name))+
    theme_minimal()
  
  return(p)
}

plot_heatmap <- function(dataframe, samples, row_annotation_cols = NA, show_colnames = F) {
  # Capture the original variable names of the lists
  sample_names <- sapply(substitute(samples)[-1], deparse)
  
  # Flatten the list of lists into a single vector of column names
  selected_columns <- unlist(samples)
  
  # Filter the dataframe to keep only the selected columns and the Genes column
  filtered_data <- dataframe %>%
    select(all_of(selected_columns), Genes)
  
  # Set 'Genes' column as row names
  filtered_data <- filtered_data %>% 
    column_to_rownames(var = "Genes")
  
  # Replace NA with the floor of the lowest value in the matrix -1
  if (any(is.na(filtered_data))) {
    impute_value <- floor(min(filtered_data, na.rm = TRUE)) - 1  # Find the lowest value, floor it, and subtract 1
    filtered_data[is.na(filtered_data)] <- impute_value  # Replace NA
    cat("NA values detected; replacing with ", as.character(impute_value), "\n")
  }
  
  # Convert to a matrix (required for pheatmap)
  intensity_matrix <- as.matrix(filtered_data)
  
  # Check if matrix is empty
  if (nrow(intensity_matrix) == 0) {
    stop("No rows with sufficient variability to cluster. Check your input data.")
  }
  
  # Create a mask for imputed values
  impute_mask <- intensity_matrix == impute_value
  # Perform row scaling (ignoring the mask during scaling)
  scaled_matrix <- t(apply(intensity_matrix, 1, scale, center = TRUE, scale = TRUE))
  # Replace imputed values to 0 for proper clustering
  scaled_matrix[impute_mask] <- 0 # Temporarily set them to 0 or some neutral value for clustering
  # Perform clustering manually (Calculate Euclidean distance)
  distance_matrix <- dist(scaled_matrix, method = "euclidean") # Calculate Euclidean distance
  # Perform hierarchical clustering
  clustering <- hclust(distance_matrix, method = "complete") # Perform hierarchical clustering
  # Replace imputed values back to NA after clustering
  scaled_matrix[impute_mask] <- NA
  
  # Create a column annotation dataframe
  annotation_col <- data.frame(Group = rep(sample_names, lengths(samples)))
  rownames(annotation_col) <- selected_columns
  
  # Define colors for annotation_col
  unique_samples <- unique(annotation_col$Group)
  annotation_colors <- list(Group = setNames(
    rainbow(length(unique_samples)), # Generate colors for each group
    unique_samples
  ))
  
  # Generate row annotations
  annotation_row <- NA
  if (any(!is.na(row_annotation_cols))) {
    # Filter for annotation columns and create row annotation
    row_annotation_data <- dataframe %>%
      select(all_of(row_annotation_cols), Genes) %>%
      column_to_rownames(var = "Genes")
    
    # Convert "+" entries to TRUE/FALSE
    annotation_row <- as.data.frame(apply(row_annotation_data, 2, function(x) x == "+"))
    
    # Add colors for row annotations to the annotation_colors list
    row_annotation_colors <- lapply(
      colnames(annotation_row),
      function(col) setNames(c("grey", "blue"), c(FALSE, TRUE))
    )
    names(row_annotation_colors) <- colnames(annotation_row)
    annotation_colors <- c(annotation_colors, row_annotation_colors)
  } else {
    annotation_row <- NA # Explicitly set to NA
  }
  
  # Generate the heatmap
  pheatmap(scaled_matrix,
           scale = "none", # No additional row scaling for visualization
           cluster_rows = clustering, # Use custom clustering
           cluster_cols = F,
           clustering_method = "complete",
           color = colorRampPalette(c("navy", "white", "firebrick"))(50), # Color gradient
           annotation_col = annotation_col,
           annotation_colors = annotation_colors, # Combine both column and row annotation colors
           annotation_row = annotation_row, # Add row annotations (NA if not provided)
           na_col = "grey",
           show_colnames = show_colnames) # Toggle column names display
}


# =============================================================================
# Quantification helpers
# =============================================================================

quant_phos_STY <- function(phos_STY = read_Perseus_file(), software = "MaxQuant"){
  #determine the amount of quantified phosphopeptides in a MaxQuant output Phospho(STY).txt
  if(class(phos_STY) == "character"){
    phos_STY <- read_Perseus_file(phos_STY)
  }
  require(tidyverse)
  
  if(software == "MaxQuant"){
    tryCatch({df <- expand_site_table(phos_STY)}, error = function(e) {message("If you quantified with DIA-NN please specify 'software = \"DIA-NN\"'")})
  df <- df %>%
    select(contains("expanded_Intensity."))}
  if(software == "DIA-NN"){
    df <- phos_STY %>% 
      filter(grepl("UniMod:21", Modified.Sequence)) %>%
      select(contains(".raw"))}
  
  df[df != 0] <- 1
  df <- df %>%
    summarise_all(sum, na.rm=T) %>%
    pivot_longer(everything(), names_to = "Experiment", values_to = "phos.quants")
  df$Experiment<-gsub("expanded_Intensity.","",df$Experiment)
  return(df)
}

quant_proteinGroups <- function(proteinGroups = read_Perseus_file(),
                                regularExpression = c(),
                                allCols = F){
  #determine the amount of quantified proteinGroups in a MaxQuant output proteinGroups.txt
  require(tidyverse)
  if(class(proteinGroups) == "character"){
    proteinGroups <- read_Perseus_file(proteinGroups)
  }
  if(!allCols){
    if(any(grep("LFQ", colnames(proteinGroups), fixed = TRUE))){ # will be True when there is column names including "LFQ" (MQ proteome)
      df <- proteinGroups %>%
        select(contains("LFQ.intensity."))
      }else if(any(grep("Intensity.", colnames(proteinGroups), fixed = TRUE))){ # will be True when there is column names including "Intensity." (MQ output)
        df <- proteinGroups %>%
          select(contains("Intensity."))
      }else if(any(grep(".raw", colnames(proteinGroups), fixed = TRUE))){ # will be True when there is column names including ".raw" (DIA-NN output)
        df <- proteinGroups %>%
          select(contains(".raw"))
      }
  }else{
      df <- proteinGroups
    }

  df[df != 0] <- 1
  df <- df %>%
    summarise_all(sum, na.rm = T) %>%
    pivot_longer(everything(), names_to = "Experiment", values_to = "prot.quants")

  # remove regex from Quants
  for(expression in regularExpression){
    df$Experiment<-gsub(expression,"",df$Experiment)
  }
  return(df)
}

phos_ratio <- function(modificationSpecificPeptides = read_Perseus_file()){
  # determine the ratio of phosphorylated peptides/all in a MaxQuant output modificationSpecificPeptides.txt
  phos_count_df <- modificationSpecificPeptides %>%
    filter(str_detect(Modifications, "Phospho")) %>%
    select(matches("Experiment")) %>%
    replace(is.na(.), 0) %>%
    summarise_all(sum) %>%
    pivot_longer(everything(), names_to = "Experiment", values_to = "count.phos")
  all_count_df <- modificationSpecificPeptides %>%
    select(matches("Experiment")) %>%
    replace(is.na(.), 0) %>%
    summarise_all(sum) %>%
    pivot_longer(everything(), names_to = "Experiment", values_to = "count.all")
  phos_df <- phos_count_df %>%
    left_join(all_count_df, by = "Experiment") %>%
    mutate(phos.ratio = count.phos/count.all) %>%
    select(Experiment,phos.ratio)
  phos_df$Experiment<-gsub("Experiment.","",phos_df$Experiment)
  return(phos_df)
}




# =============================================================================
# Miscellaneous
# =============================================================================

convertGeneList <- function(genes, to){
  # from Human/Mouse to the other
  require("biomaRt")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

  genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol",
                    values = genes , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)

  genes_df <- data.frame(matrix(unlist(genes), nrow=length(genes), byrow=TRUE))

  if(to == "mouse"){
    genes_converted <- genes_df %>% left_join(genesV2, by=c("matrix.unlist.genes...nrow...length.genes...byrow...TRUE."="HGNC.symbol"))
    genes_converted <- genes_converted %>% mutate(x = case_when(is.na(MGI.symbol) ~ matrix.unlist.genes...nrow...length.genes...byrow...TRUE.,
                                                  TRUE ~ MGI.symbol))
    return(genes_converted$x)
  }else if(to == "human"){
    genes_converted <- genes_df %>% left_join(genesV2, by=c("matrix.unlist.genes...nrow...length.genes...byrow...TRUE."="MGI.symbol"))
    genes_converted <- genes_converted %>% mutate(x = case_when(is.na(HGNC.symbol) ~ matrix.unlist.genes...nrow...length.genes...byrow...TRUE.,
                                                                TRUE ~ HGNC.symbol))
    return(genes_converted$x)
  }
}

get.ppiNCBI <- function(g.n) {
  # this is a function to parse the web page of NCBI Entrez gene, which compiles interactions from BioGRID and HPRD
  # for a given gene list it creates a protein-protein-interaction network that can be imported to cytoscape
  # courtesy David Ruau: https://www.r-bloggers.com/2012/06/obtaining-a-protein-protein-interaction-network-for-a-gene-list-in-r/
  require(XML)
  ppi <- data.frame()
  for(i in 1:length(g.n)){
    o <- htmlParse(paste("http://www.ncbi.nlm.nih.gov/gene/", g.n[i], sep=''))
    # check if interaction table exists
    exist <- length(getNodeSet(o, "//table//th[@id='inter-prod']"))>0
    if(exist){
      p <- getNodeSet(o, "//table")
      ## need to know which table is the good one
      for(j in 1:length(p)){
        int <- readHTMLTable(p[[j]])
        if(colnames(int)[2]=="Interactant"){break}
      }
      ppi <- rbind(ppi, data.frame(egID=g.n[i], intSymbol=int$`Other Gene`))
    }
    # play nice! and avoid being kicked out from NCBI servers
    Sys.sleep(1)
  }
  if(dim(ppi)[1]>0){
    ppi <- unique(ppi)
    print(paste(dim(ppi)[1], "interactions found"))
    return(ppi)
  } else{
    print("No interaction found")
  }
}



plot_volcano <- function(df_volcano = read_Perseus_file(),
                         p_value_treshold = 1,
                         pos_difference_treshold = 1.5,
                         neg_difference_treshold = -1.5,
                         plot_title = "Give your plot a name please",
                         condition_name = "condition",
                         control_name = "control",
                         color_control = "#14B1BB",
                         color_condition = "#EE611B",
                         label_genes = T,
                         gene_list = c("PDCD6"),
                         gene_list_color = "#dde329",
                         save_image = T,
                         image_name = "volcano_X",
                         treshold_lines = T,
                         highligh_significants = T,
                         force = 1,
                         text_size = 5,
                         general_point_size = 1,
                         highlighted_point_size = 1,
                         line_width = .2,
                         x_breaks = 2,
                         custom_x_limits = F,
                         x_limits = c(-10,10),
                         point_transparency_significants = 0.7,
                         point_transparency_labeled_genes = 0.7){                    
  
  #this functions needs 1 dataframe, the plain matrix export from the volcano (containing the columns "Difference" & "X.Log.P.value.")
  

  require(ggplot2)
  require(tidyverse)
  require(dplyr)
  require(ggpubr)
  require(ggrepel)
  require(ggtext)
  require(ggiraph)
  
  if("Genes" %in% colnames(df_volcano)){
    gene_column = "Genes"
  } else if("Gene.names" %in% colnames(df_volcano)){
    gene_column = "Gene.names"
  }
  
  gene_vector <- as.vector(gene_list)
  
  beautiful_volcano <- ggplot(data = df_volcano,
                              aes(x = Difference,
                                  y = X.Log.P.value.)) +
    geom_point(color = if(highligh_significants){
      ifelse(df_volcano$Difference > pos_difference_treshold & df_volcano$X.Log.P.value > p_value_treshold, color_condition,
             ifelse(df_volcano$Difference < neg_difference_treshold & df_volcano$X.Log.P.value > p_value_treshold, color_control,"gray69"))}
      else{"gray69"},
      size = general_point_size,
      alpha = ifelse(df_volcano$Difference > pos_difference_treshold & df_volcano$X.Log.P.value. > p_value_treshold,point_transparency_significants,
                     ifelse(df_volcano$Difference < neg_difference_treshold & df_volcano$X.Log.P.value. > p_value_treshold,point_transparency_significants, 0.2))) +
    {if (label_genes)
      geom_point(data = df_volcano %>% filter(.data[[gene_column]] %in% gene_list),
                 color = gene_list_color,
                 size = highlighted_point_size,
                 alpha = point_transparency_labeled_genes)}+
    ylab("-log (p-value)")+
    {if (custom_x_limits)
      scale_x_continuous(breaks = seq(min(x_limits),min(x_limits),x_breaks),
                         limits = x_limits)}+
    theme_classic() +
    labs(title = plot_title,
         x = paste0("log2 fold change (<span style = 'color:", color_condition, ";'>**", condition_name, "**</span> - ",
                    "<span style = 'color:",color_control,";'>**",control_name, "**</span>)"))+
    {if (label_genes)
      geom_label_repel(aes(label = ifelse(.data[[gene_column]] %in% gene_list, .data[[gene_column]], "")),
                       segment.color = "black",
                       max.overlaps = 10000000,
                       segment.size = line_width,
                       segment.alpha = 1,
                       nudge_x = -1,
                       size=text_size/4)}+ # the text here has a wierd size input that produces ~ 1/3 of a pt
    {if (treshold_lines)
      geom_vline(xintercept = pos_difference_treshold, alpha = 0.4, linetype = "dashed", linewidth = line_width)}+
    {if (treshold_lines)
      geom_vline(xintercept = neg_difference_treshold, alpha = 0.4, linetype = "dashed", linewidth = line_width)}+
    {if (treshold_lines)
      geom_hline(yintercept = p_value_treshold, alpha = 0.4, linetype = "dashed", linewidth = line_width)}+
    theme(legend.position = "none",
          plot.title = element_text(size = text_size,
                                    face = "bold.italic"),
          axis.title.x = element_markdown(),
          axis.text=element_text(size=text_size),
          axis.title=element_text(size=text_size),
          axis.line = element_line(colour = 'black', linewidth = line_width),
          axis.ticks = element_line(colour = "black", linewidth = line_width))+
    {if (custom_x_limits)
      scale_x_continuous(breaks = seq(min(x_limits),max(x_limits),x_breaks),
                         limits = x_limits) else {scale_x_continuous(breaks = seq(floor(min(df_volcano$Difference)), ceiling(max(df_volcano$Difference)), x_breaks))}} 

  if(save_image){
    save_wd <- choose.dir(default = "",
                          caption = "Choose in what folder you want to save the image")
    ggsave(filename = paste0(image_name, ".pdf"),
           path = save_wd,
           height = 5,
           width = 5,
           dpi = 1000,
           units = "in")
  }
  return(beautiful_volcano)
}


PTM_localization <- function(path = NULL,
                             changeWD = T,
                             organism = "mouse"){
  #this function needs to have the latest fasta & additional fasta files in the WR that you choose
  #it is important to specify organism in the arguments for the correct fasta files to be selected

  if(is.null(path)){
    path <- choose.dir()
  }
  if(changeWD==T){
    setwd(dirname(path))
  }
  #import the data
  df_ptm <- read_Perseus_file()

  #import fasta files
  require(phylotools)

  if(organism == "mouse"){
    fasta <- read.fasta(file = "UP000000589_10090_2022_mouse.fasta",
                        clean_name = F)
    fasta_add <- read.fasta(file = "UP000000589_10090_2022_mouse_additional.fasta",
                            clean_name = F)} else {
                              fasta <- read.fasta(file = "UP000005640_9606_human.fasta",
                                                  clean_name = F)
                              fasta_add <- read.fasta(file = "UP000005640_9606_additional_human.fasta",
                                                      clean_name = F)}
  fasta$seq.name <- sub(".*sp\\|", "", fasta$seq.name)

  fasta$seq.name <- sub(".*tr\\|", "", fasta$seq.name)

  fasta$seq.name <- sub("\\|.*", "", fasta$seq.name)

  colnames(fasta)[1] = "Protein.ID"

  fasta_add$seq.name <- sub(".*sp\\|", "", fasta_add$seq.name)

  fasta_add$seq.name <- sub(".*tr\\|", "", fasta_add$seq.name)

  fasta_add$seq.name <- sub("\\|.*", "", fasta_add$seq.name)

  colnames(fasta_add)[1] = "Protein.ID"

  fasta_full <- rbind(fasta,
                      fasta_add)

  df_ptm$Protein.Group <- gsub("\\;.*", "", df_ptm$Protein.Group)

  df_ptm <- left_join(df_ptm,
                      fasta_full,
                      by = c("Protein.Group" = "Protein.ID"))

  df_ptm$Modified.Sequence <- gsub("(.)(?=\\()", "\\L\\1", df_ptm$Modified.Sequence , perl = TRUE)

  df_ptm$Modified.Sequence <- gsub("[^[:alpha:]]", "", df_ptm$Modified.Sequence)

  df_ptm$Modified.Sequence <- gsub("UniMod", "", df_ptm$Modified.Sequence)

  pattern <- paste0("(?i)", df_ptm$Modified.Sequence)

  df_ptm$seq.text <- mapply(gsub,
                            pattern,
                            df_ptm$Modified.Sequence,
                            df_ptm$seq.text)

  a <- gregexpr("[a-z]", df_ptm$seq.text)

  matches_a <- regmatches(df_ptm$seq.text, a)

  l1 <- stringi::stri_locate_all(df_ptm$seq.text, regex = "[a-z]")

  df_ptm <- df_ptm %>%
    mutate(Modified_position = gregexpr("[a-z]", seq.text))

  names(df_ptm)[names(df_ptm) == 'seq.text'] <- 'Protein_sequence'

  df_ptm$Modified_position <- sub(".*c\\(", "", df_ptm$Modified_position)
  df_ptm$Modified_position <- sub("\\).*", "", df_ptm$Modified_position)
  df_ptm$Modified_position <- sub(":", ",", df_ptm$Modified_position)
  df_ptm$Modified_position <- sub(", ", ",", df_ptm$Modified_position)

  df_ptm <- df_ptm %>%
    separate_rows(Modified_position, sep = ",")

  df_ptm$Modified_position[df_ptm$Modified_position == -1] <- NaN
  
  
  df_ptm <- df_ptm %>%
    mutate(Peptide_inprotein_position = sapply(seq_len(nrow(df_ptm)), function(i) {
      peptide <- df_ptm$Modified.Sequence[i]
      protein_seq <- df_ptm$Protein_sequence[i]
      match_position <- regexpr(peptide, protein_seq)
      if (match_position == -1) {
        return(NA)  # Return NA if the peptide is not found
      } else {
        return(match_position)
      }
    }))
  
  return(df_ptm)
}

PTM_localization_V2 <- function(df_ptm = read_Perseus_file(), 
                                fasta = read.fasta(file.choose(), clean_name = F),
                                fasta_add = read.fasta(file.choose(), clean_name = F)) {
  
  # file selector not working, I dont know why (seems to fail to load df_ptm, if only that gets supplied as dataframe it works)
  # please manually load and supply files like this:
  
  # library(phylotools)
  # fasta = read.fasta("path.fasta", clean_name = F)
  # fasta_add = read.fasta("path.fasta", clean_name = F)
  # 
  # df_ptm <- read_Perseus_file("path.pr_matrix.tsv")
  # df_ptm_annot <- PTM_localization_V2(df_ptm, fasta, fasta_add)
  # 
  # write.table(df_ptm_annot, file = "path.tsv", row.names = F, col.names = T, sep = "\t", quote = FALSE, na = "")
  
  
  
  
  require(phylotools)
  require(stringr)  # For easier regex operations
  require(dplyr)    # For data manipulation
  
  # Preprocess fasta and fasta_add
  fasta$seq.name <- sub(".*sp\\|", "", fasta$seq.name)
  fasta$seq.name <- sub(".*tr\\|", "", fasta$seq.name)
  fasta$seq.name <- sub("\\|.*", "", fasta$seq.name)
  
  colnames(fasta)[1] = "Protein.ID"
  
  fasta_add$seq.name <- sub(".*sp\\|", "", fasta_add$seq.name)
  fasta_add$seq.name <- sub(".*tr\\|", "", fasta_add$seq.name)
  fasta_add$seq.name <- sub("\\|.*", "", fasta_add$seq.name)
  
  colnames(fasta_add)[1] = "Protein.ID"
  
  # Combine fasta and fasta_add
  fasta_full <- rbind(fasta, fasta_add)
  
  # Clean Protein.Group
  df_ptm$Protein.Group <- gsub("\\;.*", "", df_ptm$Protein.Group)
  
  # Join df_ptm with fasta sequences
  df_ptm <- left_join(df_ptm, fasta_full, by = c("Protein.Group" = "Protein.ID"))
  # Check if sequences were matched
  if (all(is.na(df_ptm$seq.text))) { print("No sequence could be matched. Have you used the correct Fasta files?"); return() }
  
  # Extract all UniMod numbers from the "Modified.Sequence"
  df_ptm$unimods <- str_extract_all(df_ptm$Modified.Sequence, "\\(UniMod:\\d+\\)")
  
  # Extract unique UniMod numbers
  unique_unimods <- unique(unlist(lapply(df_ptm$unimods, function(x) as.numeric(str_extract(x, "\\d+")))))
  
  # Initialize empty columns for each unique UniMod number
  for (unimod in unique_unimods) {
    df_ptm[[paste0("Position_UniMod_", unimod)]] <- NA
    df_ptm[[paste0("SequenceWindow_UniMod_", unimod)]] <- NA  # New column for sequence window
  }
  
  # Function to calculate UniMod positions relative to seq.text
  calculate_unimod_positions <- function(mod_seq, seq_text, unimod_number) {
    # Regex to match amino acid followed by the UniMod modification
    pattern <- paste0("[A-Z]\\(UniMod:", unimod_number, "\\)")
    
    # Find positions of modifications in the modified sequence
    mod_positions <- gregexpr(pattern, mod_seq)[[1]]
    if (mod_positions[1] == -1) return(NA)  # No matches found
    
    # Get clean modified sequence (remove UniMod annotations)
    clean_mod_seq <- gsub("\\(UniMod:\\d+\\)", "", mod_seq)
    
    # Align modified sequence with seq.text
    aligned_positions <- stringr::str_locate_all(seq_text, clean_mod_seq)[[1]]
    if (is.null(aligned_positions) || nrow(aligned_positions) == 0) return(NA)  # No match
    
    # Adjust mod_positions to calculate accurate positions in seq.text
    seq_text_pos <- aligned_positions[1, 1]  # The starting position of the sequence alignment in seq.text
    mod_pos_in_seq_text <- c()
    
    # Determine positions of each modification
    for (mod_pos in mod_positions) {
      if (mod_pos > 0) {
        clean_pos <- mod_pos - sum(nchar(unlist(str_extract_all(substring(mod_seq, 1, mod_pos - 1), "\\(UniMod:\\d+\\)"))))  # Adjust the position by subtracting the length of UniMod annotations up to this point
        corresponding_pos_in_seq_text <- seq_text_pos + clean_pos - 1  # Adjust relative to seq.text starting point
        mod_pos_in_seq_text <- c(mod_pos_in_seq_text, corresponding_pos_in_seq_text)
      }
    }
    
    return(paste(mod_pos_in_seq_text, collapse = ";"))
  }
  
  # Function to extract sequence window around the UniMod position
  extract_sequence_window <- function(seq_text, position) {
    position <- as.numeric(position)
    if (is.na(position)) return(NA)
    
    start_pos <- max(1, position - 15)
    end_pos <- min(nchar(seq_text), position + 15)
    
    window <- substring(seq_text, start_pos, end_pos)
    
    # Add padding with underscores if sequence is shorter than 31 characters
    padding_left <- ifelse(start_pos == 1, 15 - (position - 1), 0)
    padding_right <- ifelse(end_pos == nchar(seq_text), 15 - (nchar(seq_text) - position), 0)
    
    window <- paste0(strrep("_", padding_left), window, strrep("_", padding_right))
    
    return(window)
  }
  
  # Iterate over each row and calculate positions and sequence windows for each UniMod column
  df_ptm <- df_ptm %>%
    rowwise() %>%
    mutate(across(starts_with("Position_UniMod_"), 
                  ~ {
                    unimod_number <- as.numeric(str_replace(cur_column(), "Position_UniMod_", ""))
                    calculate_unimod_positions(Modified.Sequence, seq.text, unimod_number)
                  },
                  .names = "{col}")) %>%
    mutate(across(starts_with("Position_UniMod_"), 
                  ~ {
                    if (!is.na(.)) {
                      positions <- unlist(strsplit(as.character(.), ";"))
                      windows <- sapply(positions, function(pos) extract_sequence_window(seq.text, pos))
                      paste(windows, collapse = ";")
                    } else {
                      NA
                    }
                  },
                  .names = "{gsub('Position_', 'SequenceWindow_', col)}"))  # Ensure correct column names
  
  
  # remove "NA" strings that were created for some reason...
  df_ptm <- df_ptm %>%
    mutate(
      # Replace exact "NA", "NA;NA", and "NA;NA;NA" with NA in columns containing "SequenceWindow"
      across(contains("SequenceWindow"), ~ ifelse(. %in% c("NA", "NA;NA", "NA;NA;NA"), NA, .)),
      
      # Replace any value containing "NA" as a substring with NA in columns containing "Position"
      across(contains("Position"), ~ ifelse(grepl("NA", .), NA, .))
    )
  
  df_ptm <- df_ptm %>% select(-unimods) %>% tibble::rowid_to_column("ID")
  
  return(df_ptm)
}

SequenceWindow_size_reduction <- function(df_ptm = read_Perseus_file()){
  # this function needs a dataframe with a SequenceWindow_UniMod_21 column containing 31 amino acids (renamed from MaxQuant phospho output or my function PTM_localization_V2)
  # it will create a new column SequenceWindow_UniMod_21_size_15, that only has +-7aa around the phosphosite
  # with this you can match predicted target sequences from An atlas of substrate specificities for the human serine/threonine kinome (https://doi.org/10.1038/s41586-022-05575-3)
  
  # safe the output like this and select the input file in the pop up window:
  # write.table(SequenceWindow_size_reduction(), file = "Your path.tsv", row.names = F, col.names = T, sep = "\t", quote = FALSE, na = "")
  
  # or write a full analysis sequence (change the "Your path" before running):
  # write.table(
  #   SequenceWindow_size_reduction(
  #     PTM_localization_V2(read_Perseus_file("Your path\\report.pr_matrix.tsv"),
  #                         fasta = read.fasta("Your path.fasta", clean_name = F),
  #                         fasta_add = read.fasta("Your path.fasta", clean_name = F))), 
  #   file = "Your output path.tsv", row.names = F, col.names = T, sep = "\t", quote = FALSE, na = "")
  
  require(tidyverse)
  
  df_ptm2 <- df_ptm %>%
    mutate(SequenceWindow_UniMod_21 = strsplit(as.character(SequenceWindow_UniMod_21), ";")) %>%
    unnest(SequenceWindow_UniMod_21) %>%
    mutate(SequenceWindow_UniMod_21_size_15 = substring(SequenceWindow_UniMod_21, 9, 23)) %>%
    group_by(ID) %>%  # Group by ID before collapsing
    mutate(
      SequenceWindow_UniMod_21_size_15 = paste(SequenceWindow_UniMod_21_size_15, collapse = ";"),
      SequenceWindow_UniMod_21 = paste(SequenceWindow_UniMod_21, collapse = ";")
    ) %>%
    distinct(ID, .keep_all = TRUE) %>%  # Keep one row per group
    ungroup() %>%
    mutate(
      SequenceWindow_UniMod_21 = if_else(is.na(Position_UniMod_21), NA_character_, SequenceWindow_UniMod_21),
      SequenceWindow_UniMod_21_size_15 = if_else(is.na(Position_UniMod_21), NA_character_, SequenceWindow_UniMod_21_size_15)
    ) %>% 
    select(-ends_with(".raw"), ends_with(".raw")) # Move columns ending with ".raw" to the end
  
  return(df_ptm2)
}

calculate_coverage <- function(protein_sequence, short_sequences) {
  # this function calculates how much of a big sequence is covered by a list of small sequences
  # example use 1 is getting a protein sequence from uniprot and calculating the coverage of self specified short sequences:
  # result <- calculate_coverage(protein_sequence = "MATITCTRFTEEYQLFEELGKGAFSVVRRCVKVLAGQEYAAKI", short_sequences = c("EYQLFEELGKGAFSV","GKGAFSVVRRCVKVLAG"))
  # example use 2 is getting a protein sequence from uniprot and calculating the coverage of all quantified precursors in your DIA-NN pr.matrix:
  # result <- calculate_coverage(protein_sequence = "MATITCTRFTEEYQLFEELGKGAFSVVRRCVKVLAGQEYAAKI", short_sequences = (pr.matrix %>% filter(gene_column=="Camk2b"))$Stripped.Sequence)
  
  # DIRECTLY USE THIS (copy the following line without the "#" into the Console and then change the protein sequence to your protein sequence of interest (from Uniprot -> go to sequence -> copy sequence), change Genes to target gene name, press enter and in the window that opens select your DIA-NN pr.matrix file):
  # result <- calculate_coverage(protein_sequence = "MNIRNARPEDLMNMQHCNLLCLPENYQMKYYFYHGLSWPQLSYIAEDENGKIVGYVLAKMEEDPDDVPHGHITSLAVKRSHRRLGLAQKLMDQASRAMIENFNAKYVSLHVRKSNRAALHLYSNTLNFQISEVEPKYYADGEDAYAMKRDLTQMADELRRHLELKEKGRHVVLGAIENKVESKGNSPPSSGEACREEKGLAAEDSGGDSKDLSEVSETTESTDVKDSSEASDSAS", short_sequences = (read_Perseus_file() %>% filter(Genes=="NAA10"))$Stripped.Sequence)
  
  # Troubleshooting:
  # if you get 0 coverage you may not have quantified the specified protein, you may have a typo in Gene name or sequence, you may have the wrong organism (Sequence & Gene name)
  
  
  require(tidyverse)
  
  # Convert protein sequence to upper case
  protein_sequence <- toupper(gsub("[^A-Z]", "", protein_sequence))
  
  # Initialize a vector to store coverage information
  coverage_vector <- rep(0, nchar(protein_sequence))
  
  # Loop through each short sequence
  for (short_seq in short_sequences) {
    # Convert short sequence to upper case
    short_seq <- toupper(gsub("[^A-Z]", "", short_seq))
    
    # Check if short sequence is present in the protein sequence
    if (grepl(short_seq, protein_sequence)) {
      
      # Find positions of the short sequence in the protein sequence
      positions <- regexpr(short_seq, protein_sequence)
      
      # Update coverage vector based on the positions
      coverage_vector[positions:(positions + nchar(short_seq) - 1)] <- 1
    }
  }
  
  # Define ANSI escape codes for console text color
  black <- "\033[30m"
  grey <- "\033[37m"
  reset <- "\033[0m"
  
  # Function to print coverage vector with color formatting
  print_colored_coverage <- function(coverage_vector) {
    for (i in 1:length(coverage_vector)) {
      if (coverage_vector[i] == 1) {
        cat(black, substring(protein_sequence, i, i), reset)
      } else {
        cat(grey, substring(protein_sequence, i, i), reset)
      }
    }
    cat("\n")
  }
  
  # Calculate coverage percentage
  coverage_percentage <- sum(coverage_vector) / nchar(protein_sequence) * 100
  
  # Function to find the ranges of 1s in a vector
  find_ranges_of_ones <- function(coverage_vector) {
    # Find indices where the vector has 1s
    one_indices <- which(coverage_vector == 1)
    
    # Initialize variables
    ranges <- character(0)
    start_index <- NA
    end_index <- NA
    
    # Loop through the indices
    for (i in 1:length(one_indices)) {
      if (is.na(start_index)) {
        # Set start index if not already set
        start_index <- one_indices[i]
      } else if (i == length(one_indices) || one_indices[i] != one_indices[i + 1] - 1) {
        # Set end index if last element or if consecutive sequence breaks
        end_index <- one_indices[i]
        ranges <- c(ranges, paste(start_index, end_index, sep = "-"))
        start_index <- NA
      }
    }
    
    # Return ranges as a string
    return(paste(ranges, collapse = ", "))
  }
  
  coverage_ranges <- find_ranges_of_ones(coverage_vector)
  
  # Print coverage information
  cat("Coverage Vector:\n")
  print_colored_coverage(coverage_vector)
  cat("Coverage Ranges:\n", coverage_ranges, "\n")
  cat("Coverage Percentage:", coverage_percentage, "%\n")
  
  # Return coverage vector, ranges, and percentage
  return(list(coverage_vector = coverage_vector, coverage_ranges = coverage_ranges, coverage_percentage = coverage_percentage))
}


plot_PCA <- function(projections = read_Perseus_file(changeWD = F), ColourVariable, ShapeVariable = NA, PercentageC1, PercentageC2, plotHulls = T){
  
  require(ggalt)
  
  # both ColourVariable and ShapeVariable have to be strings
  
  
  # set colour coding variable to factor to get a discrete and not continuous colour scale
  projections[,ColourVariable] <- as.factor(projections[,ColourVariable])
  if(!is.na(ShapeVariable)){
    projections[,ShapeVariable] <- as.factor(projections[,ShapeVariable])
  } 
  
  
  # plot projections
  p_projections <- ggplot(projections, aes_string(x="Component.1",y="Component.2", col=ColourVariable, group = ColourVariable)) +
    {if(!is.na(ShapeVariable))geom_point(aes_string(x="Component.1",y="Component.2", col=ColourVariable, shape = ShapeVariable))}+
    {if(is.na(ShapeVariable))geom_point(aes_string(x="Component.1",y="Component.2", col=ColourVariable))}+
    # theme(aspect.ratio = 1/1.66)+
    theme_bw()+
    {if(plotHulls)geom_encircle(aes_string(fill = ColourVariable), s_shape = 1, expand = 0,
                                            alpha = 0.2, color = "black", show.legend = FALSE)}+
    labs(x= paste0("Component 1 [",PercentageC1,"%]"),
         y=paste0("Component 2 [",PercentageC2,"%]"))
    
  
  
  return(p_projections)
}



