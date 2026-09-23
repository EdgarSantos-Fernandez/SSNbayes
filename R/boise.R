#' Boise River spatio-temporal stream-network dataset
#'
#' @description
#' Example stream-network data from the Boise River Basin case study in
#' Santos-Fernandez et al. (2022). The example has two parts, both stored
#' under \code{inst/extdata}:
#' \itemize{
#'   \item a spatial stream network directory, \code{boise.ssn};
#'   \item an observation data frame, \code{boise_obs.RDS}.
#' }
#' \code{boise_obs.RDS} is the observation data used in
#' Santos-Fernandez et al. (2022), \emph{Computational Statistics & Data
#' Analysis}. \code{boise.ssn} is the spatial stream network supplied
#' with those observations. It contains network 3, the 42 sites used in
#' that paper, and other river networks as well.
#'
#' @details
#' \code{boise.ssn} is the directory passed to \code{SSN2::ssn_import()}.
#' It contains:
#' \itemize{
#'   \item \code{edges}: 32,852 stream segments (LINESTRING);
#'   \item \code{sites}: 226 observed sites (POINT);
#'   \item \code{preds}: 50,694 prediction locations (POINT);
#'   \item \code{netID1.dat} through \code{netID8.dat};
#'   \item \code{binaryID.db};
#'   \item \code{Variable_Definitions.xlsx}, a workbook of shapefile
#'     attribute definitions for \code{sites}, \code{preds}, and
#'     \code{edges}.
#' }
#' The prediction-site dataset is named \code{preds}. Pass
#' \code{predpts = "preds"} to import it. \code{binaryID.db} is already
#' present, so the default \code{overwrite = FALSE} is used and that file
#' is left unchanged.
#'
#' The stream features use the projected coordinate reference system stored
#' in the shapefiles (Albers, GNLCC, meters).
#'
#' \code{boise_obs.RDS} is the data frame used in the Computational
#' Statistics \& Data Analysis (CSDA) paper: daily
#' records for the 42 locations on network 3 (\code{netID} 3), on each of
#' 87 dates (3,654 rows). Every row in this file has \code{netID} 3.
#' Santos-Fernandez et al. (2022) describe these 87 dates as a systematic
#' 21-day subsample of mean daily temperatures measured from 2010-12-01 to
#' 2015-12-01. In this file the dates run from 2010-12-01 to 2015-11-11.
#'
#' \code{boise.ssn} contains other river networks in addition to network 3.
#' In \code{sites}, network 3 has the 42 sites used in the paper.
#' Networks 1, 2, and 5 contain the other observed sites (68, 105, and 11).
#' The 6,422 prediction locations in the paper are the \code{preds}
#' features with \code{netID} 3. That shapefile has 50,694 features in
#' total because it also includes prediction locations on other networks.
#'
#' The response in the case study is daily mean stream temperature
#' (\code{temp}, degrees Celsius). The published mean function uses stream
#' slope (\code{slope}), elevation (\code{elev}), cumulative drainage area
#' (\code{drainage}), daily air temperature (\code{air_temp}), and the
#' first Fourier pair (\code{sin}, \code{cos}).
#'
#' @format
#' \code{boise_obs.RDS} is a data frame with 3,654 rows and 49 columns.
#'
#' Columns used in the published analysis, or identical to a defined
#' site attribute:
#' \describe{
#'   \item{date}{Date. Observation date. Eighty-seven dates from
#'     2010-12-01 to 2015-11-11. Each of the 42 \code{locID} values
#'     occurs once on every date.}
#'   \item{locID}{Numeric location identifier. 42 locations.}
#'   \item{temp}{Numeric daily mean stream temperature in degrees Celsius.
#'     Present on all 2,620 rows with \code{dataset == "train"} and
#'     missing on all 1,034 rows with \code{dataset == "test"}.}
#'   \item{dataset}{Character. \code{"train"} (2,620 rows) or
#'     \code{"test"} (1,034 rows).}
#'   \item{air_temp}{Numeric daily air temperature in degrees Celsius,
#'     the air-temperature covariate in the case study.}
#'   \item{slope}{Numeric stream-slope covariate used in the case study.
#'     It differs from \code{SLOPE} by at most 0.0005.}
#'   \item{elev}{Integer elevation in meters. Identical to \code{Elev_M}.}
#'   \item{drainage}{Numeric cumulative drainage area in square kilometers.
#'     Identical to \code{CUMDRAINAG}.}
#'   \item{sin, cos}{Numeric first Fourier terms. Within each date,
#'     \code{sin} takes one value and \code{cos} takes one value.}
#'   \item{lat, long}{Numeric latitude and longitude in decimal degrees.}
#'   \item{aug_temp}{Numeric. Identical to \code{S1_93_11}.}
#'   \item{flow}{Numeric. Differs from \code{MA_HistCMS} by at most
#'     0.00005.}
#'   \item{pid}{Integer unique point identifier.}
#'   \item{netID}{Factor network identifier. Every row is network 3.
#'     The factor also has unused levels 1, 2, and 5.}
#'   \item{river}{Character. Every row is \code{"Boise_river".}}
#'   \item{stream_name}{Character stream name (30 streams).}
#'   \item{Stream}{Character stream name, as defined for the \code{sites}
#'     shapefile.}
#'   \item{Elev_M}{Integer elevation in meters.}
#'   \item{CBSP_ID}{Numeric site identifier.}
#'   \item{SLOPE}{Numeric stream-reach slope, in percent.}
#'   \item{CUMDRAINAG}{Numeric cumulative drainage area in square
#'     kilometers.}
#'   \item{Air_Aug}{Numeric mean August air temperature in degrees Celsius.
#'     One value is stored in this table (19.15). The workbook describes
#'     it as a basin-wide August mean.}
#'   \item{Flow_Aug}{Numeric mean August stream discharge. One value is
#'     stored in this table (4.15).}
#'   \item{S1_93_11}{Numeric historical composite: 19-year average August
#'     mean stream temperature for 1993-2011.}
#'   \item{MA_HistCMS}{Numeric mean annual flow from the historic NorWeST
#'     model, in cubic meters per second.}
#'   \item{Lat, Long_}{Numeric latitude and longitude.}
#'   \item{NEAR_FID}{Integer. FID of the line segment the site was snapped
#'     to.}
#'   \item{NEAR_DIST}{Numeric distance the site was moved when it was
#'     snapped to the network.}
#'   \item{NEAR_X, NEAR_Y}{Numeric coordinates of the snapped site.}
#'   \item{NEAR_ANGLE}{Numeric angle from the site to the nearest point on
#'     the line segment.}
#'   \item{rid}{Integer. Edge the site was snapped to.}
#'   \item{ratio}{Numeric position of the site along that edge.}
#'   \item{upDist}{Numeric upstream distance from the network outlet to
#'     the site.}
#'   \item{afvArea}{Numeric additive function value.}
#' }
#'
#' The variable-definitions workbook in \code{boise.ssn} is the source
#' for the \code{sites} attribute descriptions above.
#'
#' These columns are also present. Their names are listed without a
#' further definition: \code{UTM_Xcoord}, \code{UTM_Ycoord},
#' \code{locID_old}, \code{locID_con}, \code{LOCID}, \code{temp2},
#' \code{locID_backup}, \code{id}, \code{temp_backup},
#' \code{temp_pred_lm}, and \code{temp_pred_lm2}.
#' \code{locID_backup} is not the same vector as \code{locID}.
#' \code{temp_pred_lm} and \code{temp_pred_lm2} are missing on every
#' training row and present on every test row.
#'
#' @references
#' Santos-Fernandez, E., Ver Hoef, J. M., Peterson, E. E., McGree, J.,
#' Isaak, D. J., and Mengersen, K. (2022). Bayesian spatio-temporal
#' models for stream networks. \emph{Computational Statistics & Data
#' Analysis}, 170, 107446.
#' \doi{10.1016/j.csda.2022.107446}
#'
#' @examples
#' library(SSNbayes)
#'
#' # Path to the Boise stream-network directory
#' boise_path <- system.file(
#'   "extdata",
#'   "boise.ssn",
#'   package = "SSNbayes"
#' )
#'
#' # Import the stream network and the preds prediction sites
#' boise <- SSN2::ssn_import(
#'   boise_path,
#'   predpts = "preds"
#' )
#'
#' # Load the observation data
#' boise_obs <- readRDS(
#'   system.file(
#'     "extdata",
#'     "boise_obs.RDS",
#'     package = "SSNbayes"
#'   )
#' )
#'
#' head(boise_obs)
#' @name boise
NULL
