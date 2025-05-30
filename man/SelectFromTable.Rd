% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRSWGS.R
\name{SelectFromTable}
\alias{SelectFromTable}
\title{Interactively select rows from a table with filters and download}
\usage{
SelectFromTable(data, vars.inc = NULL, vars.dec = NULL, output.cols = NULL)
}
\arguments{
\item{data}{data frame with numeric variables}

\item{vars.inc}{character vector of variable names sought to be increased to include in the color gradient}

\item{vars.dec}{character vector of variable names sought to be decreased to include in the color gradient}

\item{output.cols}{character vector of variable names to include in the output file, if not provided, all columns will be included}
}
\value{
a shiny app with a table and a download button
}
\description{
Interactively select rows from a table with filters and download
}
\author{
Charlotte Brault
}
