% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRSWGS.R
\name{data_summary}
\alias{data_summary}
\title{Print data summary for numeric variables using gtsummary}
\usage{
data_summary(data, variables = NULL, by = NULL, add.pval = T)
}
\arguments{
\item{data}{data frame with numeric variables}

\item{variables}{character vector of numeric variable names to include in the summary}

\item{by}{character vector of variable names to group by}

\item{add.pval}{logical, if TRUE, add p-values to the summary table}
}
\value{
a gtsummary object
}
\description{
Print data summary for numeric variables using gtsummary
}
\author{
Charlotte Brault
}
