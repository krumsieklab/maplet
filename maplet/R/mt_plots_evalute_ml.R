#' Plot Evaluation Measures for Repeated Machine Learning Predictions
#'
#' This function creates plots of evaluation measures for repeated machine learning predictions. The types of plots
#' created vary depending on the cross-validation / sampling method used. Below is a list of plot types generated
#' for each method:
#' \itemize{
#'    \item CV
#'    \itemize{
#'       \item Confusion matrix - all indices
#'       \item Confusion matrices per fold (if per_fold_plots = T)
#'       \item
#'    }
#' }
#'
#'
#' @param D \code{SummarizedExperiment} input.
#' @param ml_name Name under which the machine learning results are stored; must be unique to all other machine
#'    learning results.
#' @param plot_measures A vector of measures to include; can be one or more of the following: spec, sens, f1,
#'    acc, and/or ppv. Default: c("spec", "sens", "ppv").
#' @param cutoff Threshold cutoff for prediction probability class assignment; function will select threshold
#'    closest to provided value. Default: 0.5.
#' @param per_fold_plots Whether to generate per fold plots for applicable plot types.
#'    Default: TRUE.
#'
#' @return results$output: List of evaluation measure plots.
#'
#' @examples
#' \dontrun{
#'   mt_plots_evaluate_ml(ml_name="mci_lasso", plot_measures = c("spec", "sens", "ppv", "f1"), per_fold_plots=F)
#' }
#'
#' @author KC
#'
#' @import ggplot2
#'
#' @export
mt_plots_ml_evaluate <- function(D,
                                 ml_name,
                                 plot_measures = c("spec", "sens", "ppv"),
                                 cutoff = 0.5,
                                 per_fold_plots = TRUE){

  PLOT_MEASURE_ABBR = c("spec", "sens", "ppv", "acc", "f1")

  # Validate arguments
  if("SummarizedExperiment" %in% class(D) == F) stop("D must be of class SummarizedExperiment!")
  if(missing(ml_name)) stop("A value must be provided for argument: ml_name.")

  # retrieve output from machine learning function
  #   extract predictions / test indices per fold
  res_list <- maplet:::mti_get_ml_res_by_name(D, ml_name)
  response_type <- res_list$output$response_type
  pred_list <- res_list$output$pred
  test_idx_list <- res_list$output$test_idx
  object_list <- res_list$output$object
  if(length(pred_list) > 20) maplet:::mti_logwarning("More than 20 folds detected. Plot generation may take a long time.")

  # get num_folds and labels from colData
  num_folds <- res_list$output$num_folds
  labels <- res_list$output$response_col %>% colData(D)[,.]

  if(response_type == "binary"){
    # ensure response column is of type factor - will cause errors if not
    pos_idx <- which(unique(labels)==res_list$output$pos_class)
    labels <- factor(labels, c(unique(labels)[pos_idx], unique(labels)[-pos_idx]))

    # if all indices unique, generate AUC and evaluation measure plots for all folds
    make_all_thresh_plots = FALSE
    all_test_indices <- test_idx_list %>% unname() %>% unlist()
    if(any(duplicated(all_test_indices))){
      maplet:::mti_logwarning("Indices are not unique. Varying threshold plots will not be generated.")
    }else{
      all_unique_samples = TRUE
    }

    # format data plotting, with columns: Class_Prob, Prediction, Label, Fold, Index
    pred_df <- purrr::map_df(pred_list, ~as.data.frame(.x), .id="Fold")
    test_idx_df <- purrr::map_df(test_idx_list, ~as.data.frame(.x), .id="Fold") %>% dplyr::select(-Fold)
    folds_df <- dplyr::bind_cols(pred_df, test_idx_df) %>%
      dplyr::rename(Class_Prob = 2, Index = 3) %>%
      dplyr::mutate(Prediction = as.numeric(Class_Prob > cutoff)) %>%
      dplyr::mutate(Label = as.numeric(labels[Index])-1) %>%
      dplyr::select(Class_Prob, Prediction, Label, Fold, Index)

    # Call plot generating functions
    # confusion matrix tile plot(s)
    confusion_matrix_plots <- make_confusion_matrix_plots(df = folds_df, num_folds = num_folds,
                                                          per_fold_plots = per_fold_plots, labels=labels)
    # dot plots include: AUC per fold, Evaluation Measure per fold, Folds per Evaluation Measure
    dot_plots <- make_dot_plots(df = folds_df, num_folds=num_folds, all_unique_samples = all_unique_samples,
                                plot_measures = plot_measures, cutoff = cutoff)
    # append lists
    plots <- c(confusion_matrix_plots, dot_plots)

    # only produced when indices are unique
    if(all_unique_samples){
      # threshold plots include:
      thresh_plots <- make_threshold_plots(df = folds_df, num_folds = num_folds, plot_measures = plot_measures,
                                           cutoff = cutoff)
      plots <- c(plots, thresh_plots)
    }
  }else if(response_type == "continuous"){
    #df_list <- lapply(1:num_folds, function(x){cbind(test_idx=test_idx_list[[x]], fold=x, label=labels[test_idx_list[[x]]], pred_list[[x]])})
    #fold_df <- do.call(rbind.data.frame, df_list)

    # Call plot generating functions
    # mean squared-error line plots
    mse_plots <- create_mse_plots(num_folds = num_folds, pred_list = pred_list, labels=labels, test_idx_list=test_idx_list, object_list = object_list)
    glmnet_mse_plots <- extract_glmnet_mse_plots(object_list = object_list)

    plots <- c(mse_plots, glmnet_mse_plots)

    # label-prediction correlation plots
  }else{
    stop("Value for 'response_type' must be one of: 'binary' or 'continuous'.")
  }


  # add status information & save plots
  logtxt <- glue::glue("Machine learning evaluation plots. Evalutaiton measures used ", glue::glue_collapse(plot_measures, ", "))
  funargs <- maplet:::mti_funargs()
  D %<>%
    maplet:::mti_generate_result(
      funargs = funargs,
      logtxt = logtxt,
      output = plots,
      output2 = NULL
    )

  # return
  D

}

make_confusion_matrix_plots <- function(df, num_folds, per_fold_plots, labels){

  cm_plots <- list()

  cm_table <- MLmetrics::ConfusionMatrix(y_pred = df$Prediction, y_true = df$Label) %>%
    as.data.frame %>% unlist() %>% unname()
  TClass <- as.factor(levels(labels)[cm_table[1:4]])
  PClass <- as.factor(levels(labels)[cm_table[5:8]])
  Y <- cm_table[9:12]
  cm_df <- data.frame(TClass, PClass, Y)

  p <- ggplot(data =  cm_df, mapping = aes(x = PClass, y = ordered(TClass, levels=rev(levels(TClass))))) +
    geom_tile(aes(fill = Y), colour = "white") +
    geom_text(aes(label = sprintf("%1.0f", Y)), size=6, vjust = 1) +
    xlab("Predicted") +
    ylab("True") +
    scale_fill_gradient(low = "#619CFF", high = "#F8766D") +
    theme_bw() + theme(legend.position = "none",
                       text = element_text(size=15)) +
    scale_x_discrete(position = "top")

  cm_plots$cm_all_folds <- p

  if(per_fold_plots){
    cm_per_fold <- lapply(1:num_folds, function(i, fdf = df){
      f <- glue::glue("Fold{i}")
      fold_df <- dplyr::filter(fdf, Fold == f)
      cm <- MLmetrics::ConfusionMatrix(fold_df$Prediction, fold_df$Label) %>%
        as.data.frame() %>%unlist() %>% unname()
      TClass <- as.factor(levels(labels)[cm_table[1:4]])
      PClass <- as.factor(levels(labels)[cm_table[5:8]])
      Y <- cm[9:12]
      cm_fold_df <- data.frame(TClass, PClass, Y, Fold=rep(f, 4))
    })

    cm_per_fold_df <- do.call(rbind, cm_per_fold)

    p <- ggplot(data =  cm_per_fold_df, mapping = aes(x = PClass, y = ordered(TClass, levels=rev(levels(TClass))))) +
      geom_tile(aes(fill = Y), colour = "white") +
      geom_text(aes(label = sprintf("%1.0f", Y)), vjust = 1) +
      xlab("Predicted") +
      ylab("True") +
      scale_fill_gradient(low = "#619CFF", high = "#F8766D") +
      theme_bw() + theme(legend.position = "none") +
      scale_x_discrete(position = "top") +
      facet_wrap(~Fold)

    cm_plots$cm_per_fold <- p
  }

  cm_plots

}


make_dot_plots <- function(df, num_folds,  all_unique_samples, plot_measures, cutoff){

  ### ------ Format Data ------ ###
  # for each fold, calculate the six evaluation measures
  #  combine list into data frame (folds as rows, evaluation measures as columns)
  em_per_fold <- lapply(1:num_folds, function(i, fdf = df){
    f <- glue::glue("Fold{i}")
    fold_df <- dplyr::filter(fdf, Fold == f)
    em_df <- list(Fold = f,
                  AUC = MLmetrics::AUC(fold_df$Prediction, fold_df$Label) %>% round(2),
                  sens = MLmetrics::Sensitivity(fold_df$Prediction, fold_df$Label) %>% round(2),
                  spec = MLmetrics::Specificity(fold_df$Prediction, fold_df$Label) %>% round(2),
                  f1 = MLmetrics::F1_Score(fold_df$Prediction, fold_df$Label) %>% round(2),
                  acc = MLmetrics::Accuracy(fold_df$Prediction, fold_df$Label) %>% round(2),
                  ppv = MLmetrics::Precision(fold_df$Prediction, fold_df$Label) %>% round(2))
  }) %>% do.call(rbind.data.frame, .)

  # Add an "Average of Folds" row with mean of measures
  avg_of_folds <- dplyr::summarise_all(em_per_fold, mean) %>%
    round(2) %>%
    dplyr::mutate(Fold = "Average") # replace 'NA' value generated for Fold column
  em_per_fold %<>% rbind(avg_of_folds)

  # if all samples are unique (e.g. in classic cross-validation), calculate the six evaulation measures
  #   for all samples together and include "All Samples" row in data frame
  if(all_unique_samples){
    all_samples_em <- list(Fold = "All",
                           AUC = MLmetrics::AUC(df$Prediction, df$Label) %>% round(2),
                           sens = MLmetrics::Sensitivity(df$Prediction, df$Label) %>% round(2),
                           spec = MLmetrics::Specificity(df$Prediction, df$Label) %>% round(2),
                           f1 = MLmetrics::F1_Score(df$Prediction, df$Label) %>% round(2),
                           acc = MLmetrics::Accuracy(df$Prediction, df$Label) %>% round(2),
                           ppv = MLmetrics::Precision(df$Prediction, df$Label) %>% round(2))
    em_per_fold %<>% rbind(all_samples_em)
  }

  ### ------ Make Plots ------ ###
  dot_plots <- list()

  # plot AUC per fold
  p <- ggplot(data = em_per_fold, aes(x = 1, y = AUC)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(shape=Fold, color=Fold), size=5, stroke=2, position = position_dodge2(w = 0.1)) +
    {if(all_unique_samples){
      scale_shape_manual(name = "Fold",
                         labels = c("All", "Average", paste0("Fold", 1:num_folds)),
                         values=c(17, 4, rep(16, num_folds)))
    }else{
      scale_shape_manual(name = "Fold",
                         labels = c("Average", paste0("Fold", 1:num_folds)),
                         values=c(4, rep(16, num_folds)))
    }} +
    ggtitle("AUC (Per Fold)") +
    theme(text = element_text(size=15),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  dot_plots$auc_per_fold <- p

  # plot fold per evaluation measure
  p <- reshape2::melt(em_per_fold, id=c("Fold")) %>%
    dplyr::filter(variable %in% plot_measures) %>%
    ggplot(., aes(x=factor(variable), y=value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(shape=Fold, color=Fold), size=5, stroke=2, position = position_dodge2(w = 0.2)) +
    {if(all_unique_samples){
      scale_shape_manual(name = "Fold",
                         labels = c("All", "Average", paste0("Fold", 1:num_folds)),
                         values=c(17, 4, rep(16, num_folds)))
    }else{
      scale_shape_manual(name = "Fold",
                         labels = c("Average", paste0("Fold", 1:num_folds)),
                         values=c(4, rep(16, num_folds)))
    }} +
    ylab("Value per Fold") +
    ggtitle(paste0("Evaluation Measures (Per Fold), cutoff: ", round(cutoff,2))) +
    theme(text = element_text(size=15))
  dot_plots$fold_per_em <- p

  # plot evaluation measures per fold
  p <- dplyr::select(em_per_fold, c("Fold", all_of(plot_measures))) %>%
    reshape2::melt(id=c("Fold")) %>%
    dplyr::filter(grepl("Fold", Fold)) %>%
    ggplot(., aes(x= factor(Fold), y=value, color=variable)) +
    geom_point(size=5) +
    xlab("Fold") +
    ylab("Value") +
    ggtitle(paste0("Evaluation Measures (Per Fold), cutoff: ", round(cutoff,2))) +
    theme(text = element_text(size=15),
          axis.title.x=element_blank()) +
    guides(color=guide_legend(title="Measure"))
  dot_plots$em_per_fold <- p

  # return the three dot plots
  dot_plots

}


make_threshold_plots <- function(df, num_folds, plot_measures, cutoff){

  thresh_plots <- list()

  measures_list = maplet:::mti_get_measures_list(trueclass=df$Label, predprob=df$Class_Prob)

  # (1) ROC curve
  roc_data=data.frame(FPR=rev(1-measures_list$specvals),
                      sensitivity=rev(measures_list$sensvals))
  roc_plot <- ggplot(data=roc_data, aes(x=FPR, y=sensitivity)) +
    geom_line(size=2) + geom_abline(intercept=0, slope=1, colour="blue") +
    ggtitle(sprintf("AUC: %.3f (All Folds)", measures_list$AUC)) + xlab("FPR") + ylab("TPR") +
    xlim(0,1) + ylim(0,1) +
    theme(text = element_text(size=15))
  thresh_plots$roc <- roc_plot

  # (2) Evaluation measures all samples
  measures_data <- data.frame(sens=measures_list$sensvals,
                              spec=measures_list$specvals,
                              ppv=measures_list$ppvvals,
                              acc=measures_list$accvals,
                              f1=measures_list$f1vals,
                              thresholds=measures_list$thresholds)
  measures_data <- measures_data[,c(plot_measures, "thresholds")]
  measures_plot_data <- reshape2::melt(measures_data, id=c("thresholds"))
  measures_plot <- ggplot(measures_plot_data, aes(x=thresholds, y=value, color=variable)) + geom_line(size=2) +
    xlab("Thresholds (values)") +
    ylab("Evaluation Measure Values") +
    labs(color='Measures') +
    geom_vline(xintercept = cutoff) +
    ggtitle(paste0("Evaluaiton Measures (All Folds), cutoff: ", round(cutoff,2))) +
    theme(text = element_text(size=15))
  thresh_plots$em_thresh_all <- measures_plot


  # (3) Evaluation measures per fold
  em_per_fold <- lapply(1:num_folds, function(i, fdf = df){
    f <- glue::glue("Fold{i}")
    fold_df <- dplyr::filter(fdf, Fold == f)
    fold_measures_list <- maplet:::mti_get_measures_list(trueclass = fold_df$Label, predprob = fold_df$Class_Prob)
    em_df <- data.frame(Fold = rep(f, length(fold_measures_list$sensvals)),
                  sens = fold_measures_list$sensvals,
                  spec = fold_measures_list$specvals,
                  f1 = fold_measures_list$f1vals,
                  acc = fold_measures_list$accvals,
                  ppv = fold_measures_list$ppvvals,
                  thresholds = fold_measures_list$thresholds) %>%
      dplyr::select(c(all_of(plot_measures), "thresholds", "Fold")) %>%
      reshape2::melt(id=c("thresholds", "Fold"))
  })

  # em_plots <- lapply(1:num_folds, function(i, em_list = em_per_fold){
  #
  #   em_df <- em_list[[i]]
  #   p <- ggplot(em_df, aes(x=thresholds, y=value, color=variable)) + geom_line(size=2) +
  #     xlab("Thresholds") +
  #     ylab("Evaluation Measure Values") +
  #     labs(color='Measures') +
  #     ggtitle(paste0("Evaluation Measures: Fold ", i)) +
  #     geom_vline(xintercept = cutoff) +
  #     theme(text = element_text(size=15))
  #
  # })
  # thresh_plots$measures_per_fold <- em_plots

  # (4) Evaluation Measures Per Fold - Faceted
  em_all <- do.call(rbind.data.frame, em_per_fold)
  p <- ggplot(em_all, aes(x=thresholds, y=value, color=variable)) + geom_line(size=2) +
    xlab("Thresholds") +
    ylab("Evaluation Measure Values") +
    labs(color='Measures') +
    ggtitle(paste0("Evaluation Measures Per Fold ")) +
    geom_vline(xintercept = cutoff) +
    theme(text = element_text(size=15)) +
    facet_wrap(~Fold)


  thresh_plots$facet_em_line <- p

  thresh_plots

}

create_mse_plots <- function(num_folds, pred_list, labels, test_idx_list, object_list){

  # not giving me same results as glmnet plot - why?
  mse_plots <- lapply(1:num_folds, function(x){

    obj <- object_list[[x]]
    til <- test_idx_list[[x]]
    pl <- pred_list[[x]]

    sqrd_err <- (labels[til] - pl)^2
    cvm <- colMeans(sqrd_err)
    cvsd <- apply(sqrd_err, 2, sd) / sqrt(nrow(sqrd_err))
    cvup <- cvm+cvsd
    cvlo <- cvm+cvsd

    plot_df <- data.frame(log_lambda=log(obj$lambda), mean=cvm, se=cvsd)

    g <- ggplot(plot_df, aes(x=log_lambda,y=mean)) +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
      geom_point(color="red") +
      geom_vline(xintercept = log(obj$lambda.min), linetype="dashed") +
      geom_vline(xintercept = log(obj$lambda.1se), linetype="dashed") +
      xlab("Log(lambda)") +
      ylab("Mean-Squared Error") +
      ggtitle(glue::glue("Fold {x}")) +
      theme_minimal()

  })

}

extract_glmnet_mse_plots <- function(object_list){

  glmnet_mse_plots <- list()

  tmp <- lapply(1:length(object_list), function(x){

    obj <- object_list[[x]]

    png(glue::glue("fold{x}.png"))
    dev.control('enable')
    plot(obj)
    title(glue::glue("Fold {x}", line=2))
    p = recordPlot()
    dev.off()
    #file.remove("tmp.png")

    glmnet_mse_plots[[x]] <- p

  })

  glmnet_mse_plots

}

create_cor_plots <- function(object){

  print("This function is empty!")

}
