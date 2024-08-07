% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_func.R
\name{predict.ssnbayes}
\alias{predict.ssnbayes}
\title{Performs spatio-temporal prediction in R using an ssnbayes object from a fitted model.}
\usage{
\method{predict}{ssnbayes}(
  object = object,
  ...,
  path = path,
  obs_data = obs_data,
  pred_data = pred_data,
  net = net,
  nsamples = nsamples,
  addfunccol = addfunccol,
  locID_pred = locID_pred,
  chunk_size = chunk_size,
  seed = seed
)
}
\arguments{
\item{object}{A stanfit object returned from ssnbayes}

\item{...}{Other parameters}

\item{path}{Path with the name of the SpatialStreamNetwork object}

\item{obs_data}{The observed data frame}

\item{pred_data}{The predicted data frame}

\item{net}{(optional) Network from the SSN object}

\item{nsamples}{The number of samples to draw from the posterior distributions. (nsamples <= iter)}

\item{addfunccol}{The variable used for spatial weights}

\item{locID_pred}{(optional) the location id for the predictions. Used when the number of pred locations is large.}

\item{chunk_size}{(optional) the number of locID to make prediction from}

\item{seed}{(optional) A seed for reproducibility}
}
\value{
A data frame with the location (locID), time point (date), plus the MCMC draws from the posterior from 1 to the number of iterations.
The locID0 column is an internal consecutive location ID (locID) produced in the predictions, starting at max(locID(observed data)) + 1. It is used internally in the way predictions are made in chunks.
}
\description{
It will take an observed and a prediction data frame.
It requires the same number of observation/locations per day.
It requires location id (locID) and points id (pid).
The locID are unique for each site.
The pid is unique for each observation.
Missing values are allowed in the response but not in the covariates.
}
\details{
The returned data frame is melted to produce a long dataset. See examples.
Currently, the predict() function produces predictions for normal random variables. However, this can be easily transformed in to counts (Poisson distributed) and presence/absence (binomial distributed).
}
\examples{
\donttest{
#require('SSNdata')
#clear_preds <- readRDS(system.file("extdata/clear_preds.RDS", package = "SSNdata"))
#clear_preds$y <- NA
#pred <- predict(object = fit_ar,
#                 path = path,
#                 obs_data = clear,
#                 pred_data = clear_preds,
#                 net = 2,
#                 nsamples = 100, # numb of samples from the posterior
#                 addfunccol = 'afvArea', # var for spatial weights
#                 locID_pred = locID_pred,
#                 chunk_size = 60)
}
}
\author{
Edgar Santos-Fernandez
}
