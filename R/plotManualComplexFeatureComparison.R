
plotManualComplexFeatureComparison <- function(mapped_features,
                                        traces,
                                        manual_stats,
                                        calibration,
                                        complexID,
                                        plot_peak=TRUE,
                                        plot_monomer=TRUE,
                                        log=FALSE,
                                        plot_apex=TRUE,
                                        plot_in_complex_estimate=TRUE,
                                        TP=FALSE) {

  TP_result_ids <- as.numeric(unlist(strsplit(manual_stats$TP_result_ids,split=";")))
  FP_result_ids <- as.numeric(unlist(strsplit(manual_stats$FP_result_ids,split=";")))
  FN_manual_ids <- as.numeric(unlist(strsplit(manual_stats$FN_manual_ids,split=";")))

  features <- subset(mapped_features, complex_id == complexID | complex_id.manual == complexID)
  peptides <- unique(unlist(strsplit(features$subunits_annotated, split = ";")))
  #peptides <- traces$trace_annotation[complex_id == complexID]$id
  traces <- subset(traces,trace_ids = peptides)
  traces.long <- toLongFormat(traces$traces)
  setkey(traces.long, id)

  TP_idx=which(features$results_id %in% TP_result_ids)
  FP_idx=which(features$results_id %in% FP_result_ids)
  FN_idx=which(features$manual_id %in% FN_manual_ids)

  features$status = NA
  if (length(TP_idx) > 0){
    features$status[TP_idx]="TP"
  }
  if (length(FP_idx) > 0){
    features$status[FP_idx]="FP"
  }
  if (length(FN_idx) > 0){
    features$status[FN_idx]="FN"
  }

  proteinName = complexID

  features$status = as.factor(features$status)
  group.colors <- c("TP" = "green3", "FP" = "firebrick1", "FN"="orange")
  group.colors <- group.colors[which(names(group.colors) %in% features$status)]

  #grid.newpage()

  if (TP==FALSE){
    p <- ggplot(traces.long) +
      geom_line(aes_string(x='fraction', y='intensity', color='id'), na.rm = TRUE) +
      ggtitle(proteinName) +
      xlab('fraction') +
      ylab('intensity') +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      theme(plot.margin = unit(c(1,.5,.5,.5),"cm")) +
      theme(plot.title = element_text(vjust=19,size=12)) +
      guides(color=FALSE)
    p <- p + geom_rect(data=features,aes(xmin = left_pp, xmax = right_pp, ymin = 0, ymax = Inf),fill="blue",alpha = 0.15, na.rm = TRUE)
    p <- p + geom_rect(data=features,aes(xmin = left_pp.manual, xmax = right_pp.manual, ymin = 0, ymax = Inf),fill="green3",alpha=0.15, na.rm = TRUE)
    p <- p + geom_vline(data=features,aes(xintercept=apex), linetype="solid",colour="blue",alpha=0.5, na.rm = TRUE)
    p <- p + geom_vline(data=features,aes(xintercept=apex.manual), linetype="solid",colour="green3",alpha=0.5, na.rm = TRUE)
    #print(p)
  }

  if (TP){
    p <- ggplot(traces.long) +
      geom_line(aes_string(x='fraction', y='intensity', color='id'), na.rm = TRUE) +
      ggtitle(proteinName) +
      xlab('fraction') +
      ylab('intensity') +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      theme(plot.margin = unit(c(1,.5,.5,.5),"cm")) +
      theme(plot.title = element_text(vjust=19, size=12)) +
      guides(color=FALSE)
    p <- p + geom_rect(data=features,aes(xmin = left_pp, xmax = right_pp, ymin = 0, ymax = Inf,fill=status),alpha = 0.15, na.rm = TRUE)
    p <- p + geom_rect(data=features,aes(xmin = left_pp.manual, xmax = right_pp.manual, ymin = 0, ymax = Inf,fill=status),alpha=0.15, na.rm = TRUE)
    p <- p + geom_vline(data=features,aes(xintercept=apex), colour="grey2", linetype="dashed", alpha=0.5,size=1.5, na.rm = TRUE)
    p <- p + geom_vline(data=features,aes(xintercept=apex.manual), colour="green4", linetype="solid",alpha=0.5,size=1.5, na.rm = TRUE)
    p <- p + scale_fill_manual(values=group.colors)
    print(p)
  }

}
