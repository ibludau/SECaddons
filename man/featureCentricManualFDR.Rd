% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manual_FDR_functions.R
\name{featureCentricManualFDR}
\alias{featureCentricManualFDR}
\title{featureCentricManualFDR}
\usage{
featureCentricManualFDR(detected_features, manual_features,
  grid_search_list = FALSE)
}
\arguments{
\item{detected_features}{data.table containing complex feature finding results}

\item{manual_features}{data.table with manual annotations.}

\item{grid_search_list}{logical if grid search results are used}
}
\value{
List with stats - feature centric, all features were mapped between automated and manual
}
\description{
featureCentricManualFDR.
}
