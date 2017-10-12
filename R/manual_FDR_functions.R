#' estimateGridSearchComplexCentricManualFDR
#' @description estimateGridSearchComplexCentricManualFDR.
#' @param complex_features_list list containing complex feature finding results for different parameter sets.
#' @param manual_annotation data.table with manual annotations.
#' @return List with stats - complex centric, no feature mapping performed
#' @export
estimateGridSearchComplexCentricManualFDR<- function(complex_features_list,manual_annotation){
  x=lapply(complex_features_list,complexCentricManualFDR,manual=manual_annotation,grid_search_list=TRUE)
  x_names = names(x[[1]])
  y = as.data.table(t(setDT(x)))[,lapply(.SD,unlist)]
  names(y)=x_names
  y[]
}

#' complexCentricManualFDR
#' @description complexCentricManualFDR.
#' @param auto data.table containing complex feature finding results
#' @param manual data.table with manual annotations.
#' @param grid_search_list logical if grid search results are used
#' @return List with stats - complex centric, no feature mapping performed
#' @export
complexCentricManualFDR <- function(auto,manual,grid_search_list=FALSE){
  high.confidence.manual <- subset(manual, confidence=="High")
  if("complex_id" %in% names(auto)){
    detected.complexes <- unique(auto$complex_id)
    true.complexes <- unique(manual$complex_id)
    high.confidence.true.complexes <- unique(high.confidence.manual$complex_id)
  } else if ("protein_id" %in% names(auto)) {
    detected.complexes <- unique(auto$protein_id)
    true.complexes <- unique(manual$protein_id)
    high.confidence.true.complexes <- unique(high.confidence.manual$protein_id)
  } else {
    message("error in input data: neither complex_id nor protein_id column found")
  }

  TP <- sum(detected.complexes %in% true.complexes)
  FP <- sum(!(detected.complexes %in% true.complexes))
  FDR <- FP / (TP + FP)

  TP.high.confidence <- sum(detected.complexes %in% high.confidence.true.complexes)

  if(grid_search_list){
    list(FDR=FDR,TP=TP,TP_high=TP.high.confidence,corr=unique(auto$corr),window=unique(auto$window),rt_height=unique(auto$rt_height),smoothing_length=unique(auto$smoothing_length),peak_corr_cutoff=unique(auto$peak_corr_cutoff),completeness_cutoff=unique(auto$completeness_cutoff),n_subunits_cutoff=unique(auto$n_subunits_cutoff))
  } else {
    list(FDR=FDR,TP=TP,TP_high=TP.high.confidence)
  }
}

#' estimateGridSearchFeatureCentricManualFDR
#' @description estimateGridSearchFeatureCentricManualFDR.
#' @param complex_features_list list containing complex feature finding results for different parameter sets.
#' @param manual_annotation data.table with manual annotations.
#' @param parallelized logical
#' @param n_cores numeric number of cores to use (if parallelized == TRUE)
#' @return List with stats - feature centric, all features were mapped between automated and manual
#' @export
estimateGridSearchFeatureCentricManualFDR<- function(complex_features_list,manual_annotation,parallelized=FALSE,n_cores=1){
  if (parallelized) {
    cl <- snow::makeCluster(n_cores)
    # setting a seed is absolutely crutial to ensure reproducible results!!!!!!!!!!!!!!!!!!!
    clusterSetRNGStream(cl,123)
    doSNOW::registerDoSNOW(cl)
    clusterEvalQ(cl,library(SECprofiler))
    clusterEvalQ(cl,library(SECaddons))
    #clusterExport(cl, list("featureCentricManualFDR","manualFeatureMapping","estimate_errors","resolveDoubleAssignments"))
    x <- parLapply(cl,complex_features_list,fun=featureCentricManualFDR,manual_features=manual_annotation,grid_search_list=TRUE)
    stopCluster(cl)
  } else {
    x=lapply(complex_features_list,featureCentricManualFDR,manual_features=manual_annotation,grid_search_list=TRUE)
  }
  x <- lapply(x, "[[", 2)
  x_names = names(x[[1]])
  y = as.data.table(t(setDT(x)))[,lapply(.SD,unlist)]
  names(y)=x_names
  y[]
}

#' featureCentricManualFDR
#' @description featureCentricManualFDR.
#' @import data.table
#' @param detected_features data.table containing complex feature finding results
#' @param manual_features data.table with manual annotations.
#' @param grid_search_list logical if grid search results are used
#' @return List with stats - feature centric, all features were mapped between automated and manual
#' @export
featureCentricManualFDR <- function(detected_features, manual_features,grid_search_list=FALSE){
  if("complex_name" %in% names(detected_features)){
    detected_features[ , count := .N, by = complex_name]
    detected_features <- detected_features[order(-rank(count), complex_name)]
  } else if ("protein_id" %in% names(detected_features)) {
    detected_features[ , count := .N, by = protein_id]
    detected_features <- detected_features[order(-rank(count), protein_id)]
  } else {
    message("error in input data: neither complex_id nor protein_id column found")
  }
  results <- detected_features
  results[,results_id := seq_len(.N)]
  results$manual_id <- NA
  for (i in seq_len(nrow(results))) {
    results$manual_id[i] <- manualFeatureMapping(results[i],dist_cutoff=3,manual_annotation=manual_features)
  }
  results[,manual_id := as.integer(manual_id)]
  results[,manual_id_mapp := manual_id]
  manual_features[,manual_id := as.integer(manual_id)]
  results_merged <- merge(results,manual_features,by="manual_id",all=TRUE,suffixes=c("",".manual"))

  # don't allow same manual id to be assigned to two features
  results_merged[,apex_dist:= abs(apex-apex.manual)]
  results_merged[,left_pp_dist:= abs(left_pp-left_pp.manual)]
  results_merged[,right_pp_dist:= abs(right_pp-right_pp.manual)]
  results_merged[,peak_dist := (0.5*apex_dist)+(0.25*left_pp_dist)+(0.25*right_pp_dist)]
  # @TODO test if correctly done
  unique_manual_ids <- unique(results_merged$manual_id_mapp[which((!is.na(results_merged$manual_id_mapp)) & (results_merged$manual_id_mapp != 0))])
  result_ids_to_unassign <-  unlist(lapply(unique_manual_ids,resolveDoubleAssignments,mapped_features=results_merged))
  if (length(result_ids_to_unassign) > 0) {
    results_merged$manual_id_mapp[which(results_merged$results_id %in% result_ids_to_unassign)]=0
    results_merged$manual_id[which(results_merged$results_id %in% result_ids_to_unassign)]=0
  }
  stats <- estimate_errors(results_merged)
  if(grid_search_list){
    stats <- c(stats,corr=unique(detected_features$corr),window=unique(detected_features$window),rt_height=unique(detected_features$rt_height),smoothing_length=unique(detected_features$smoothing_length),peak_corr_cutoff=unique(detected_features$peak_corr_cutoff),completeness_cutoff=unique(detected_features$completeness_cutoff),n_subunits_cutoff=unique(detected_features$n_subunits_cutoff))
  }
  list(data=results_merged,stats=stats)
}

#' manualFeatureMapping
#' @description manualFeatureMapping.
#' @param feature feature
#' @param dist_cutoff dist_cutoff
#' @param manual_annotation manual_annotation
#' @export
manualFeatureMapping <- function(feature,dist_cutoff,manual_annotation){
  if("complex_id" %in% names(feature)){
    manual_features <- subset(manual_annotation,complex_id==feature$complex_id)
  } else if ("protein_id" %in% names(feature)) {
    manual_features <- subset(manual_annotation,protein_id==feature$protein_id)
  } else {
    message("error in input data: neither complex_id nor protein_id column found")
  }

  if (nrow(manual_features) == 0) {
    match = 0
  } else {
    apex_dist <- abs(feature$apex - manual_features$apex)
    if (min(apex_dist)[1] > dist_cutoff) {
      match = 0
    } else {
      sel_min_apex_dist <- which(apex_dist==min(apex_dist))
      if (length(sel_min_apex_dist) == 1) {
        match <- manual_features$manual_id[sel_min_apex_dist]
      } else {
        manual_features <- manual_features[sel_min_apex_dist]
        left_dist <- abs(feature$left_pp - manual_features$left_pp)
        right_dist <- abs(feature$right_pp - manual_features$right_pp)
        boundary_dist <- left_dist + right_dist
        sel_min_boundary_dist <- which(boundary_dist==min(boundary_dist))
        if (length(sel_min_boundary_dist) == 1) {
          match <- manual_features$manual_id[sel_min_boundary_dist]
        } else {
          message("two manual features found")
          match <- manual_features$manual_id[sel_min_boundary_dist][1]
        }
      }
    }
  }
  match
}

#' estimate_errors
#' @description estimate_errors.
#' @param table table
#' @export
estimate_errors <- function(table){
  FP_result_ids=table$results_id[which(table$manual_id_mapp == 0)]
  FN_manual_ids=table$manual_id[which(is.na(table$manual_id_mapp))]
  FN_high_manual_ids=table$manual_id[which(is.na(table$manual_id_mapp) & (table$confidence == "High"))]
  TP_result_ids=table$results_id[which((!is.na(table$manual_id_mapp)) & (table$manual_id_mapp != 0))]
  TP_high_result_ids = table$results_id[which((!is.na(table$manual_id_mapp)) & (table$manual_id_mapp != 0) & (table$confidence == "High"))]

  FP=length(FP_result_ids)
  FN=length(FN_manual_ids)
  FN_high=length(FN_high_manual_ids)
  TP=length(TP_result_ids)
  TP_high = length(TP_high_result_ids)

  FDR=FP/(TP+FP)

  # compare peak apex and boundary precision
  apex_r2 = summary(lm(apex~apex.manual, data = table))$adj.r.squared
  left_pp_r2 = summary(lm(left_pp~left_pp.manual, data = table))$adj.r.squared
  right_pp_r2 = summary(lm(right_pp~right_pp.manual, data = table))$adj.r.squared
  summary_r2 = (0.5*apex_r2)+(0.25*left_pp_r2)+(0.25*right_pp_r2)

  list(
    TP=TP,
    TP_high=TP_high,
    FP=FP,
    FN=FN,
    FN_high=FN_high,
    FDR=FDR,
    summary_r2=summary_r2,
    FP_result_ids=paste(FP_result_ids,collapse=";"),
    TP_result_ids=paste(TP_result_ids,collapse=";"),
    FN_manual_ids=paste(FN_manual_ids,collapse=";")
  )
}

#' resolveDoubleAssignments
#' @description resolveDoubleAssignments
#' @param id feature_id
#' @param mapped_features data.table with all mapped features
#' @return data.table with mached detected and manual features
#' @export
resolveDoubleAssignments <- function(id,mapped_features){
  features <- subset(mapped_features,manual_id_mapp==id)
  if (nrow(features) > 1) {
    keep <- features$results_id[which(features$peak_dist==(min(features$peak_dist)))[1]]
    unassign <- features$results_id[which(features$results_id != keep)]
    unassign
  }
}

#' plotManualVsDecoyFDR
#' @description plotManualVsDecoyFDR.
#' @param stats mergen manual and automatic stats
#' @param FDR_cutoff numeric default = 0.1
#' @param colour_parameter character string default = "completeness_cutoff"
#' @param PDF logical default=TRUE
#' @param name character string default = "ManualVsDecoyFDR_plot.pdf"
#' @return List with stats - feature centric, all features were mapped between automated and manual
#' @export
plotManualVsDecoyFDR <- function(stats,FDR_cutoff=0.1,colour_parameter="completeness_cutoff",PDF=TRUE,name="ManualVsDecoyFDR_plot.pdf"){
  if(PDF){pdf(name)}
  pl <- ggplot(data=stats,aes(x=FDR.auto,y=FDR.manual,colour=get(colour_parameter))) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    scale_x_continuous(breaks=seq(0,1,0.1),limits=c(0,1),minor_breaks=NULL) +
    scale_y_continuous(breaks=seq(0,1,0.1),limits=c(0,1),minor_breaks=NULL) +
    labs(colour = paste0(eval(colour_parameter),"\n")) +
    theme_bw()
  print(pl)
  if(PDF){dev.off()}
}
