% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_func.R
\name{ssnbayes}
\alias{ssnbayes}
\title{Fits a mixed linear regression model using Stan}
\usage{
ssnbayes(
  formula = formula,
  data = data,
  path = path,
  time_method = time_method,
  space_method = space_method,
  iter = 3000,
  warmup = 1500,
  chains = 3,
  refresh = max(iter/100, 1),
  net = 1,
  addfunccol = addfunccol,
  loglik = FALSE,
  ppd = FALSE,
  seed = seed
)
}
\arguments{
\item{formula}{A formula as in lm()}

\item{data}{A long data frame containing the locations, dates, covariates and the response variable. It has to have the locID and date. No missing values are allowed in the covariates.
The order in this data.fame MUST be: spatial locations (1 to S) at time t=1, then locations (1 to S) at t=2 and so on.}

\item{path}{Path with the name of the SpatialStreamNetwork object}

\item{time_method}{A list specifying the temporal structure (ar = Autorregressive; var = Vector autorregression) and coumn in the data with the time variable.}

\item{space_method}{A list defining if use or not of an SSN object and the spatial correlation structure. The second element is the spatial covariance structure. A 3rd element is a list with the lon and lat for Euclidean distance models.}

\item{iter}{Number of iterations}

\item{warmup}{Warm up samples}

\item{chains}{Number of chains}

\item{refresh}{Sampler refreshing rate}

\item{net}{The network id (optional). Used when the SSN object contains multiple networks.}

\item{addfunccol}{Variable to compute the additive function. Used to compute the spatial weights.}

\item{loglik}{Logic parameter denoting if the loglik will be computed by the model.}

\item{ppd}{Produce the posterior predictive distribution}

\item{seed}{(optional) A seed for reproducibility}
}
\value{
A list with the model fit

It returns a ssnbayes object (similar to stan returns). It includes the formula used to fit the model. The output can be transformed into the stanfit class using class(fits) <- c("stanfit").
}
\description{
It requires the same number of observation/locations per day.
It requires location id (locID) and points id (pid).
The locID are unique for each site.
The pid is unique for each observation.
Missing values are allowed in the response but not in the covariates.
}
\details{
Missing values are not allowed in the covariates and they must be imputed before using ssnbayes(). Many options can be found in https://cran.r-project.org/web/views/MissingData.html
The pid in the data has to be consecutive from 1 to the number of observations.
Users can use the SpatialStreamNetwork created with the SSN package. This will provide the spatial stream information used to compute covariance matrices. If that is the case, the data has
to have point ids (pid) matching the ones in SSN distance matrices, so that a mapping can occur.
}
\examples{
\dontrun{
#options(mc.cores = parallel::detectCores())
# Import SpatialStreamNetwork object
#path <- system.file("extdata/clearwater.ssn", package = "SSNbayes")
#n <- SSN2::ssn_import(path, predpts = "preds", overwrite  = TRUE)
## Imports a data.frame containing observations and covariates
#clear <- readRDS(system.file("extdata/clear_obs.RDS", package = "SSNbayes"))
#fit_ar <- ssnbayes(formula = y ~ SLOPE + elev + h2o_area + air_temp + sin + cos,
#                   data = clear,
#                   path = path,
#                   time_method = list("ar", "date"),
#                   space_method = list('use_ssn', c("Exponential.taildown")),
#                   iter = 2000,
#                   warmup = 1000,
#                   chains = 3,
#                   net = 2, # second network on the ssn object
#                   addfunccol='afvArea')
#space_method options examples
#use list('no_ssn', 'Exponential.Euclid', c('lon', 'lat')) if no ssn object is available
}
}
\author{
Edgar Santos-Fernandez
}
