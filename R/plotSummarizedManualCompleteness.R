#' plotSummarizedManualCompleteness
#' @description plotSummarizedManualCompleteness
#' @param hypotheses data.table with complex hypotheses
#' @param protTraces traces object of type proteins
#' @param PDF logical default = TRUE
#' @export
plotSummarizedManualCompleteness <- function(hypotheses,protTraces,PDF=TRUE){
  manual_annotation <- manual_corum_annotation
  manual_annotation <- manual_annotation[grep("DECOY",manual_annotation$complex_id,invert=TRUE)]
  setorder(manual_annotation,complex_id,-completeness)
  manual_annotation_best <- unique(manual_annotation,by="complex_id")
  manual_annotation_best_min50 <- subset(manual_annotation_best,completeness>=0.5)
  manual_annotation_best_lower50 <- subset(manual_annotation_best,completeness<0.5)

  proteins_in_hypotheses <- unique(hypotheses$protein_id)
  proteins_in_traces <- unique(protTraces$traces$id)
  hypotheses[,annotated:=1]
  hypotheses[,detected:=ifelse(protein_id %in% proteins_in_traces, 1, 0)]
  #hypotheses[,protein_collapsed := paste(protein_id,collapse=";"),by=complex_id]
  hypotheses[,annotated_collapsed := sum(annotated),by=complex_id]
  hypotheses[,detected_collapsed := sum(detected),by=complex_id]
  hypotheses[,ms_completeness := detected_collapsed/annotated_collapsed,by=complex_id]
  unique_hypotheses <- unique(hypotheses,by="complex_id")
  unique_hypotheses_50 <- subset(unique_hypotheses,ms_completeness >= 0.5)

  completenessSummary <- data.table(
    name=c("no co-elution","co-elution\n(< 50% complete)","co-elution\n(>= 50% complete)"),
    count=c(sum(!(unique_hypotheses_50$complex_id %in% manual_annotation_best$complex_id)),
    sum(unique_hypotheses_50$complex_id %in% manual_annotation_best_lower50$complex_id),
    sum(unique_hypotheses_50$complex_id %in% manual_annotation_best_min50$complex_id)
    )
  )

  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  if(PDF){pdf("manualComplexCompletenessPie.pdf")}
    print(pie(x=completenessSummary$count,labels=paste0(completenessSummary$name,"\n",completenessSummary$count),col=cbPalette[1:nrow(completenessSummary)]))
  if(PDF){dev.off()}

  if(PDF){pdf("manualComplexCompletenessScatter.pdf")}
    p <- ggplot(data=manual_annotation_best,aes(x=n_subunits_annotated,y=n_subunits,colour=completeness)) +
      geom_point() +
      geom_abline(intercept=0,slope=1) +
      geom_abline(intercept=0,slope=0.5, linetype=2) +
      annotate("text",x=70,y=80, label="100%", angle = 45) +
      annotate("text",x=70,y=9.5, label="50%", angle = 23) +
      scale_x_log10(name = "N subunits in hypothesis", breaks = c(1,10,100), limits = c(1,200)) +
      scale_y_log10(name = "N subunits observed co-eluting", breaks = c(1,10,100), limits = c(1,200)) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()) +
      labs(colour = "completeness\n") +
      theme(legend.justification=c(0,1), legend.position=c(0.05,0.75))
    print(p)
  if(PDF){dev.off()}
}
