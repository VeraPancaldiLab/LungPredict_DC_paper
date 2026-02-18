compute.anova.test = function(data, coldata, trait, pval = 0.05){
  
  sig = c()
  coldata[,trait] = as.factor(coldata[,trait])
  
  for (j in 1:ncol(data)) {
    df = data.frame("Value" = data[,j], "Trait" = coldata[,trait])
    
    # skip if feature is constant
    if (length(unique(df$Value)) < 2) next
    
    # try ANOVA safely
    res.aov <- tryCatch({
      df %>% rstatix::anova_test(Value ~ Trait)
    }, error = function(e) NULL)
    
    if(!is.null(res.aov)){
      if(round(res.aov$p, 5) <= pval){
        sig = c(sig, j)
      }
    }
  }
  
  features_sig = data[,sig, drop = FALSE]
  
  if(length(sig)==0){
    message("No significant features (pvalue < ", pval, ") after Anova test")
    return(NULL)
  } else {
    return(features_sig)
  }
}


# Function to extract markers for the group
extract_markers <- function(res, lfc_thresh = 1, padj_thresh = 0.05, direction = c("up", "down", "both")) {
  direction <- match.arg(direction)
  
  res <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene")
  
  if (direction == "up") {
    genes <- subset(res, log2FoldChange >  lfc_thresh & padj < padj_thresh)$gene
  } else if (direction == "down") {
    genes <- subset(res, log2FoldChange < -lfc_thresh & padj < padj_thresh)$gene
  } else { # both
    genes <- subset(res, abs(log2FoldChange) > lfc_thresh & padj < padj_thresh)$gene
  }
  
  return(genes)
}

plot_features_across_clusters <- function(data_groups, features_df, file_title, filter_significant = FALSE, pval_threshold = 0.05) {
  
  cluster_group <- data_groups %>%
    dplyr::mutate(Sample = rownames(.)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(MaxCluster = names(.)[which.max(c_across(-Sample))]) %>%
    dplyr::select(Sample, MaxCluster)
  
  features_long <- features_df %>%
    data.frame() %>%
    dplyr::mutate(Sample = rownames(.)) %>%
    tidyr::pivot_longer(cols = -Sample, names_to = "Feature", values_to = "Value") %>%
    dplyr::left_join(cluster_group, by = "Sample")
  
  if (filter_significant) {
    feature_pvals <- features_long %>%
      group_by(Feature) %>%
      summarise(pvalue = summary(aov(Value ~ MaxCluster))[[1]][["Pr(>F)"]][1], .groups = "drop")
    
    sig_features <- feature_pvals %>%
      filter(pvalue < pval_threshold) %>%
      pull(Feature)
    
    features_long <- features_long %>%
      filter(Feature %in% sig_features)
  }
  
  features_summary <- features_long %>%
    dplyr::group_by(Feature, MaxCluster) %>%
    dplyr::summarise(mean = mean(Value), .groups = "drop")
  
  p <- ggplot(features_summary, aes(x = MaxCluster, y = mean, group = Feature, color = Feature)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    theme_bw() +
    labs(
      x = "Cluster",
      y = "Mean feature value",
      title = file_title
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

plot_features_weighted_across_clusters <- function(data_groups, features_df, file_title, filter_significant = FALSE, pval_threshold = 0.05) {
  
  # Ensure sample names exist
  data_groups <- data_groups %>%
    as.data.frame() %>%
    dplyr::mutate(Sample = rownames(.))
  
  features_df <- as.data.frame(scale(features_df, center = TRUE, scale = TRUE))
  
  features_df <- features_df %>%
    as.data.frame() %>%
    dplyr::mutate(Sample = rownames(.))
  
  # Convert both to long format
  memberships_long <- data_groups %>%
    tidyr::pivot_longer(cols = -Sample, names_to = "Cluster", values_to = "Weight")
  
  features_long <- features_df %>%
    tidyr::pivot_longer(cols = -Sample, names_to = "Feature", values_to = "Value")
  
  # Merge to associate each sample's feature value with its cluster weights
  combined_long <- dplyr::left_join(memberships_long, features_long, by = "Sample")
  
  # Compute weighted mean per feature per cluster
  features_summary <- combined_long %>%
    dplyr::group_by(Feature, Cluster) %>%
    dplyr::summarise(
      weighted_mean = sum(Value * Weight, na.rm = TRUE) / sum(Weight, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Optionally filter features by weighted ANOVA significance
  if (filter_significant) {
    feature_pvals <- combined_long %>%
      group_by(Feature) %>%
      summarise(
        pvalue = {
          # approximate significance using weighted linear model
          tryCatch({
            fit <- lm(Value ~ Cluster, weights = Weight)
            anova(fit)$`Pr(>F)`[1]
          }, error = function(e) NA)
        },
        .groups = "drop"
      )
    
    sig_features <- feature_pvals %>%
      filter(!is.na(pvalue), pvalue < pval_threshold) %>%
      pull(Feature)
    
    features_summary <- features_summary %>%
      filter(Feature %in% sig_features)
  }
  
  # Plot weighted means
  p <- ggplot(features_summary, aes(x = Cluster, y = weighted_mean, group = Feature, color = Feature)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    theme_bw() +
    labs(
      x = "Cluster",
      y = "Weighted mean feature value",
      title = file_title
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
}

# Volcano plot with labels
plot_volcano <- function(res, group_name, lfc_thresh = 2, padj_thresh = 0.05) {
  
  res <- as.data.frame(res)
  
  res$Significance <- "Not Significant"
  res$Significance[res$log2FoldChange >  lfc_thresh & res$padj < padj_thresh] <- "Up in group"
  res$Significance[res$log2FoldChange < -lfc_thresh & res$padj < padj_thresh] <- "Up in rest"
  
  ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up in group" = "red", 
                                  "Up in rest" = "blue", 
                                  "Not Significant" = "grey")) +
    geom_vline(xintercept = c(-lfc_thresh, lfc_thresh), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", color = "black") +
    geom_text_repel(
      data = subset(res, padj < padj_thresh & abs(log2FoldChange) > lfc_thresh) %>%
        tibble::rownames_to_column("gene"),
      aes(label = gene),
      size = 3,
      max.overlaps = 15
    ) +
    labs(title = paste0("Volcano Plot: ", group_name, " vs Rest"),
         x = "log2 Fold Change",
         y = "-log10 adjusted p-value") +
    theme_minimal()
}

minMax <- function(x) {
  #columns: features
  x = data.matrix(x)
  for(i in 1:ncol(x)){
    x[,i] = (x[,i] - min(x[,i], na.rm = T)) / (max(x[,i], na.rm = T) - min(x[,i], na.rm = T))
  }
  
  return(x)
  
}

diff_analysis = function(RNA.counts.normalized, feature, test, ref, pval = 0.05, TF.collection = "CollecTRI", aracne_path = NULL, file.name){
  
  if(TF.collection == "CollecTRI"){
    net = decoupleR::get_collectri(organism = 'human', split_complexes = F) #Get universe

    collectri2viper_regulons <- function(df) {
      regulon_list <- split(df, df$source)
      regulons <- lapply(regulon_list, function(regulon) {
        tfmode <- stats::setNames(regulon$mor, regulon$target)
        list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
      })

      return(regulons)
    }

    net_regulons = collectri2viper_regulons(net)
  }else if(TF.collection == "ARACNE"){
    #network_file = read.delim2("~/Documents/ARACNe-AP/outputFolder/network.txt")
    #write.table(network_file, file = "~/Documents/ARACNe-AP/outputFolder/network.txt", sep = "\t", row.names = F, quote = F, col.names = FALSE)
    net_regulons <- viper::aracne2regulon(aracne_path, as.matrix(RNA.counts.normalized), format = "3col")
  }

  ####TFs differential active 
  # Generating test and ref data
  test_i <- which(feature == test)
  ref_i <- which(feature == ref)
  
  mat_test <- as.matrix(RNA.counts.normalized[,test_i, drop = F])
  mat_ref <- as.matrix(RNA.counts.normalized[,ref_i, drop = F])
  
  # Generating NULL model (test, reference)
  dnull <- ttestNull(mat_test, mat_ref, per=1000, verbose = T)
  
  # Generating signature: Observed DEG scores
  signature <- rowTtest(mat_test, mat_ref) #Computes t-statistic for each gene
  signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) * sign(signature$statistic)) #Transform pvalues to z scores to standarize strength and direction
  signature <- na.omit(signature) #Remove NAs
  signature <- signature[,1] #Take scores as a numeric vector
  
  # Running msVIPER
  mra <- msviper(signature, net_regulons, dnull, verbose = F, minsize = 5)
  
  tfs_da = mra$es$p.value
  tfs_da = tfs_da[tfs_da < pval]
  tfs_names = names(tfs_da)
  
  nes = mra$es$nes
  nes = nes[tfs_names] #only significant
  up = names(nes[nes > 0]) #up active TFs
  down = names(nes[nes < 0]) #down active TFs
  
  # Plot DiffActive TFs
  svg(paste0("Results/Differential_TFs_", file.name, ".svg"),
      width = 20, height = 15)
  print(plot(mra, mrs = up, cex = 4, sep = 1.5, include = c("activity")))
  dev.off()
  
  return(list(TFs_UP = up, TFs_DOWN = down, mra))
  
}

compute.survival.analysis = function(survival.data, PFS, PFS_event, file_name = NULL, 
                                     features = NULL, p.value = 0.05, thres = 0.5, group_column = NULL) {
  require(survival)
  
  # Prepare the data for survival analysis
  data_for_model <- data.frame("time" = survival.data[,PFS],
                               "status" = survival.data[,PFS_event]) 
  
  ###### CASE 1: Use predefined groups (like patient clusters or risk groups)
  if (!is.null(group_column)) {
    
    # Use an existing column in survival.data to define patient groups
    data_for_model$coxHL <- paste0("Group_", survival.data[[group_column]])
    
    # Fit Kaplan–Meier survival curves across groups
    km_fit <- survival::survfit(survival::Surv(time, status) ~ coxHL, data = data_for_model)
    
    #median PFS 
    medians <- survminer::surv_median(km_fit)
    
    median_text <- paste0(
      paste0(medians$strata, ": mPFS = ", round(medians$median, 2), 
             " (", round(medians$lower,2), "-", round(medians$upper, 2), ")"),
      collapse = "\n"
    )
    
    #median follow-up (reverse KM) 
    rev_km <- survfit(Surv(time, 1 - status) ~ 1, data = data_for_model)
    median_follow_up <- summary(rev_km)$table["median"]
    
    followup_text <- paste0("Median follow-up = ", round(median_follow_up, 2))
    
    # Log-rank test for survival difference among groups
    pval <- survminer::surv_pvalue(km_fit, data = data_for_model)$pval 
    
    # Extract group names for plot legend
    strata_names <- gsub("coxHL=", "", names(km_fit$strata)) 
    legend_labels <- strata_names
    legend.labs = legend_labels
    
    p = survminer::ggsurvplot(km_fit,
                                data = data_for_model,
                                size = 1,
                                conf.int.style = "step",
                                pval = TRUE,
                                pval.size = 8,
                                risk.table = TRUE,
                                risk.table.col = "strata",
                                fontsize = 9,
                                legend.labs = legend.labs,
                                risk.table.height = 0.3,
                                ggtheme = theme_grey(),
                                title = paste0("Cox PH for Surv(time, status) ~ ", group_column),
                                xlab = "Time to death/recurrence/progression"
    )
    
    p$table <- p$table + theme(legend.position = "none")
    
    combined <- gridExtra::arrangeGrob(p$plot, p$table, ncol = 1, heights = c(2/3, 1/3))
    ggsave(
      filename = paste0("Results/SurvPlot_", group_column, "_", file_name, ".svg"),
      plot = combined,
      width = 10,
      height = 6,
      units = "in",
      device = "svg"
    )

    return(list(
      km_plot = p,
      median_PFS = medians,
      p_value = pval,
      median_follow_up = median_follow_up
    ))
    
  } 
  
  ###### CASE 2: Use features to define high/low groups automatically
  else if (!is.null(features)) {
    significant_combinations <- list() # To store significant feature combinations
    contador = 1
    
    n_features = ncol(features)
    
    for (n in 1:n_features) {
      
      # Create a formula 
      formula <- as.formula(paste("Surv(time, status) ~", colnames(features)[n]))
      
      # Prepare the dataset containing survival time, event status, and feature values
      data_for_model <- data.frame("time" = survival.data[,PFS],
                                   "status" = survival.data[,PFS_event]) 
      
      # Compute the quantile cutoff used to split the feature into high vs. low groups
      quantiles <- quantile(features[,n], thres) 
      
      # Divide patients into high and low groups based on the quantile threshold of the feature
      data_for_model$coxHL <- factor(ifelse(features[,n] >= quantiles, paste0('High_', colnames(features)[n]), paste0("Low_", colnames(features)[n]))) 
      
      # Fit Kaplan–Meier survival curves for the high and low feature groups
      km_fit <- survival::survfit(survival::Surv(time, status) ~ coxHL, data = data_for_model)
      
      #  #Performs log-rank test to see whether both survival curves are significantly different between the groups
      pval <- survminer::surv_pvalue(km_fit, data = data_for_model)$pval
      
      if (!is.na(pval) && pval < p.value) {
        
        significant_combinations[[contador]] <- formula
        names(significant_combinations)[contador] = colnames(features)[n]
        
        ### Set colors of each curve with the correct order 
        strata_names <- gsub("coxHL=", "", names(km_fit$strata)) 
        legend_labels <- strata_names
        legend.labs = legend_labels
        
        p = survminer::ggsurvplot(km_fit,
                                    data = data_for_model,
                                    size = 1,
                                    palette = c("#E7B800", "#2E9FDF"),
                                    conf.int.style = "step",
                                    pval = TRUE,
                                    pval.size = 8,
                                    risk.table = TRUE,
                                    risk.table.col = "strata",
                                    legend.labs = legend.labs,
                                    risk.table.height = 0.3,
                                    fontsize = 9,
                                    ggtheme = ggplot2::theme_grey(),
                                    title = paste0("Cox PH for Surv(time, status) ~", colnames(features)[n]),
                                    xlab = "Time to death/recurrence/progression")
        
        p$table <- p$table + theme(legend.position = "none")
        combined <- gridExtra::arrangeGrob(p$plot, p$table, ncol = 1, heights = c(2/3, 1/3))
        ggsave(
          filename = paste0("Results/SurvPlot_", names(significant_combinations)[contador], "_", file_name, ".svg"),
          plot = combined,
          width = 10,
          height = 6,
          units = "in",
          device = "svg"
        )

        contador = contador + 1
      }
    }
    
    if (length(significant_combinations) == 0) {
      print("No significant combinations found.")
    } else {
      return(significant_combinations)
    }
    
  } else {
    stop("Please provide either 'group_column' for predefined groups or 'features' for feature-based analysis.")
  }
}

compute.survival.models = function(features, survival.data, PFS, PFS_event, file_name = NULL, p.value = 0.05, thres = 0.5, max_factors = Inf) {
  n_features <- ncol(features)
  significant_combinations <- list() # To store significant feature combinations
  
  # Generate all possible combinations of the features
  contador = 1
  for (n in 1:min(n_features, max_factors)) {
    combinations <- combn(1:n_features, n, simplify = FALSE)
    
    for (comb in combinations) {
      # Create a formula dynamically based on the combination
      formula <- as.formula(paste("Surv(time, status) ~", paste(colnames(features)[comb], collapse = " + ")))
      
      # Prepare the data for survival analysis
      data_for_model <- data.frame("time" = survival.data[,PFS],
                                   "status" = survival.data[,PFS_event]) 
      
      data_for_model = cbind(data_for_model, features[,comb, drop=F])
      
      # Fit the Cox PH model with the combination of features (cox PH take into account covariates and measure the impact of each variable in the survival time)
      cox <- rms::cph(formula, data = data_for_model)
      data_for_model$CoxPredictors <- cox$linear.predictors #linear predictors is the risk score for each individual in the dataset
      
      # Compute the quantile cutoff used to split between high vs. low risk groups
      quantiles <- quantile(data_for_model$CoxPredictors, thres) 
      
      # Binarize the Cox model output to draw two KM lines (linear predictors are used to stratify between high-risk and low-risk groups)
      data_for_model$coxHL <- ifelse(cox$linear.predictors >= quantiles, 'High', "Low") 
      
      # Perform Kaplan-Meier based on coxHL
      km_fit <- survival::survfit(Surv(time, status) ~ coxHL, data = data_for_model)
      
      pval <- survminer::surv_pvalue(km_fit, data = data_for_model)$pval #Performs log-rank test to see whether both survival curves are significantly different
      
      if (!is.na(pval) && pval < p.value) {
        significant_combinations[[contador]] <- formula
        names(significant_combinations)[contador] = paste(colnames(features)[comb], collapse = " + ")
        
        pdf(paste0("Results/SurvPlot_", names(significant_combinations)[contador], "_", file_name), width = 10, height = 5, onefile = FALSE)
        print(survminer::ggsurvplot(km_fit,
                                    data = data_for_model,
                                    size = 1,
                                    palette = c("#E7B800", "#2E9FDF"),
                                    conf.int.style = "step",
                                    pval = TRUE,
                                    risk.table = TRUE,
                                    risk.table.col = "strata",
                                    legend.labs = c("High", "Low"),
                                    risk.table.height = 0.3,
                                    ggtheme = theme_grey(),
                                    title = paste0("Cox PH for ", names(significant_combinations)[contador]),
                                    xlab = "Time to death/recurrence/progression")
        )
        dev.off()
        contador = contador + 1
      }
    }
  }
  
  if (length(significant_combinations) == 0) {
    print("No significant combinations found.")
  } else {
    return(significant_combinations)
  }
  
}

cell.groups.wilcox.test <- function(cell.groups, coldata, trait, pval = 0.05) {
  sig = c()
  coldata[, trait] = as.factor(coldata[, trait])
  
  if (length(unique(coldata[, trait])) != 2) {
    stop("Wilcoxon test requires a binary trait (exactly two groups).")
  }
  
  for (j in 1:ncol(cell.groups[[1]])) {
    data = data.frame(Value = cell.groups[[1]][, j], Trait = coldata[, trait])
    
    # Perform Wilcoxon rank-sum test
    res.test <- stats::wilcox.test(Value ~ Trait, data = data, exact = FALSE)
    
    if (round(res.test$p.value, 5) <= pval) {
      cat("Significant p-value after Wilcoxon test for", colnames(cell.groups[[1]])[j], "\n")
      
      # Save boxplot
      pdf(paste0("Results/Wilcoxon_", trait, "_", colnames(cell.groups[[1]])[j], ".pdf"),
          width = 12, height = 9)
      print(
        ggplot(data, aes(x = Trait, y = Value, fill = Trait)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "gray30") +  # hide boxplot outliers
          geom_jitter(width = 0.15, size = 2.5, alpha = 0.8, color = "black") +  # show all sample points
          stat_summary(fun = median, geom = "point", shape = 23,
                       size = 3, fill = "white") +
          scale_fill_brewer(palette = "Set2") +
          labs(
            title = paste0("Wilcoxon Test - ", colnames(cell.groups[[1]])[j]),
            subtitle = paste0("P-value: ", round(res.test$p.value, 5)),
            x = paste0("Clinical trait: ", trait),
            y = "Cell group score"
          ) +
          theme_minimal(base_size = 16) +
          theme(
            axis.text.x = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 12),
            axis.title.y = element_text(size = 20, angle = 90),
            plot.title = element_text(size = 20),
            plot.subtitle = element_text(size = 15, face = "bold")
          ) +
          scale_fill_discrete(name = trait)
      )
      dev.off()
      
      sig = c(sig, j)
    }
  }
  
  # Collect significant cell groups
  if (length(sig) == 0) {
    message("No significant cell groups (p-value < ", pval, ") after Wilcoxon test.")
    return(invisible(NULL))
  } else {
    cell.groups.sig = list()
    cell.groups.sig[[1]] = cell.groups[[1]][, sig, drop = FALSE]
    cell.groups.sig[[2]] = cell.groups[[2]][sig]
    cell.groups.sig[[3]] = cell.groups[[3]][sig]
    return(cell.groups.sig)
  }
}

cell.groups.anova.test = function(cell.groups, coldata, trait, pval = 0.05){
  sig = c()
  coldata[, trait] = as.factor(coldata[, trait])
  
  for (j in 1:ncol(cell.groups[[1]])) {
    data = data.frame(Value = cell.groups[[1]][, j], Trait = coldata[, trait])
    res.aov <- rstatix::anova_test(data = data, dv = Value, between = Trait)
    
    if (res.aov$p < pval) {
      cat("Significant p-value after ANOVA for", colnames(cell.groups[[1]])[j], "\n")
      
      # Tukey post-hoc test + position for significance annotation
      pwc <- data %>% 
        rstatix::tukey_hsd(Value ~ Trait) %>% 
        rstatix::add_xy_position(x = "Trait")
      
      # Generate annotated boxplot (same style as first function)
      pdf(paste0("Results/ANOVA_", trait, "_", colnames(cell.groups[[1]])[j], ".pdf"), 
          width = 12, height = 9)
      print(
        ggpubr::ggboxplot(data, x = "Trait", y = "Value", fill = "Trait", add = "jitter") +
          ggpubr::stat_pvalue_manual(pwc, hide.ns = TRUE) +
          labs(
            title = "One-way ANOVA with Tukey HSD",
            subtitle = rstatix::get_test_label(res.aov, detailed = TRUE),
            caption = rstatix::get_pwc_label(pwc),
            x = paste0("Clinical trait: ", trait),
            y = paste0("Values for ", colnames(cell.groups[[1]])[j])
          ) +
          theme(
            axis.text.x = element_text(size = 20, angle = 0),
            axis.title.x = element_text(size = 20),
            legend.position = "none",
            axis.title.y = element_text(size = 20, angle = 90),
            plot.title = element_text(size = 20),
            plot.subtitle = element_text(size = 15, face = "bold")
          )
      )
      dev.off()
      
      sig = c(sig, j)
    }
  }
  
  cell.groups.sig = list()
  if (length(sig) == 0) {
    message("No significant cell groups (p-value < ", pval, ") after ANOVA test")
    return(NULL)
  } else {
    cell.groups.sig[[1]] = cell.groups[[1]][, sig, drop = FALSE]
    cell.groups.sig[[2]] = cell.groups[[2]][sig]
    cell.groups.sig[[3]] = cell.groups[[3]][sig]
    return(cell.groups.sig)
  }
}


