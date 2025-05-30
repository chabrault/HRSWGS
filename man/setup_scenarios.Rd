% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HRSWGS.R
\name{setup_scenarios}
\alias{setup_scenarios}
\title{Setup prediction scenarios}
\usage{
setup_scenarios(
  myPheno,
  scenario,
  envs.train = NULL,
  envs.pred = NULL,
  ignore.genos = NULL,
  traits = NULL,
  genos = NULL,
  prop.CVS1 = 0.8
)
}
\arguments{
\item{myPheno}{data frame containing at least the following columns: gid (genotype identifier), env, value (trait phenotype)}

\item{scenario}{one of knEnv, knLoc.knYr, knLoc.nYr, nLoc.knYr, nLoc.nYr.
knEnv is known environment, including a sample of observations from the testing environment, as well as all other environments, knLoc.knYr is known location, known year, knLoc.nYr is known location / unknown year, nLoc.knYr is unknown location, known year, nLoc.nYr is unknown location / unknown year.}

\item{envs.train}{(optional) restriction of environments to be in training set, if not provided all environments are considered}

\item{envs.pred}{(optional) restriction of environments to be in testing set, if not provided all environments are considered}

\item{ignore.genos}{(optional) character vector of genotypes to exclude from testing set, while including them in training set, typically check genotypes}

\item{traits}{(optional) character vector for trait column(s)}

\item{genos}{(optional) character vector of restricted genotype list, typically genotyped lines.}

\item{prop.CVS1}{numeric, proportion of genotypes in testing set for scenario knEnv, default is 0.8.}
}
\value{
list with one lement per testing (environment) set.
Within each testing set, list of training environment, testing environment, genotypes in training and testing sets, and phenotypic data in training and testing sets.
}
\description{
Setup prediction scenarios
}
\author{
Charlotte Brault
}
