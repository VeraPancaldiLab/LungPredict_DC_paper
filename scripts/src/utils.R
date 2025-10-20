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

diff_analysis = function(RNA.counts.normalized, feature, test, ref, pval = 0.05, file.name){
  
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

  # Plot DiffActive TFs
  pdf(paste0("Results/Differential_TFs_", file.name, ".pdf"))
  print(plot(mra, mrs=15, cex=1, include = c("expression","activity")))
  dev.off()
  
  tfs_da = mra$es$p.value
  tfs_da = tfs_da[tfs_da < pval]
  tfs_names = names(tfs_da)
  
  nes = mra$es$nes
  nes = nes[tfs_names] #only significant
  up = names(nes[nes > 0]) #up active TFs
  down = names(nes[nes < 0]) #down active TFs
  
  return(list(TFs_UP = up, TFs_DOWN = down, mra))
  
}

compute.survival.by.groups = function(survival.data, PFS, PFS_event, time_unit, file_name) {
  library(survival)
  
  # Prepare the data for survival analysis
  data_for_model <- data.frame("time" = survival.data[,PFS],
                               "status" = survival.data[,PFS_event]) 
  
  data_for_model$coxHL <- paste0("Group_", survival.data$Groups)
  
  # Perform Kaplan-Meier based on coxHL
  km_fit <- survival::survfit(survival::Surv(time, status) ~ coxHL, data = data_for_model)
  
  pval <- survminer::surv_pvalue(km_fit, data = data_for_model)$pval #Performs log-rank test to see whether both survival curves are significantly different
  
  ### Set colors of each curve with the correct order 
  strata_names <- gsub("coxHL=", "", names(km_fit$strata)) 
  legend_labels <- strata_names
  legend.labs = legend_labels
  
  pdf(paste0("Results/SurvPlot_", file_name), width = 10, height = 5, onefile = FALSE)
  print(survminer::ggsurvplot(km_fit,
                              data = data_for_model,
                              size = 1,
                              conf.int.style = "step",
                              pval = TRUE,
                              risk.table = TRUE,
                              risk.table.col = "strata",
                              legend.labs = legend.labs,
                              risk.table.height = 0.3,
                              ggtheme = theme_grey(),
                              title = paste0("Cox PH for Surv(time, status) ~ patient groups"),
                              xlab = paste0("Time to death/recurrence/progression (", time_unit, ")")
  ))
  dev.off()
}

compute.survival.by.features = function(features, survival.data, PFS, PFS_event, time_unit, p.value = 0.05, thres = 0.5, max_factors = Inf) {
  library(survival)
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
      
      quantiles <- quantile(features[,comb], thres) 
      
      # Binarize the Cox model output to draw two KM lines (linear predictors are used to stratify between high-risk and low-risk groups)
      data_for_model$coxHL <- factor(ifelse(features[,comb] >= quantiles, paste0('High_', colnames(features)[comb]), paste0("Low_", colnames(features)[comb]))) 
      
      feature_name <- colnames(features)[comb]
      
      # Perform Kaplan-Meier based on coxHL
      km_fit <- survival::survfit(survival::Surv(time, status) ~ coxHL, data = data_for_model)
      
      pval <- survminer::surv_pvalue(km_fit, data = data_for_model)$pval #Performs log-rank test to see whether both survival curves are significantly different
      
      if (!is.na(pval) && pval < p.value) {
        
        significant_combinations[[contador]] <- formula
        names(significant_combinations)[contador] = paste0("Formula_", contador)
        
        ### Set colors of each curve with the correct order 
        strata_names <- gsub("coxHL=", "", names(km_fit$strata)) 
        legend_labels <- strata_names
        legend.labs = legend_labels
        
        pdf(paste0("Results/SurvPlot_", names(significant_combinations)[contador]), width = 10, height = 5, onefile = FALSE)
        print(survminer::ggsurvplot(km_fit,
                                    data = data_for_model,
                                    size = 1,
                                    palette = c("#E7B800", "#2E9FDF"),
                                    conf.int.style = "step",
                                    pval = TRUE,
                                    risk.table = TRUE,
                                    risk.table.col = "strata",
                                    legend.labs = legend.labs,
                                    risk.table.height = 0.3,
                                    ggtheme = theme_grey(),
                                    title = paste0("Cox PH for ", paste("Surv(time, status) ~", paste(colnames(features)[comb], collapse = " + "))),
                                    subtitle = direction,
                                    xlab = paste0("Time to death/recurrence/progression (", time_unit, ")")
        ))
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
