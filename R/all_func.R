#' Collapses a SpatialStreamNetwork object into a data frame
#'
#' @param ssn An S4 SpatialStreamNetwork object created with SSN2 package.
#' @param par A spatial parameter such as the computed_afv (additive function value).
#' @return A data frame with the lat and long of the line segments in the network. The column line_id refers to the ID of the line.
#' @importFrom dplyr arrange
#' @importFrom plyr .
#' @importFrom sf st_read st_geometry_type as_Spatial st_union
#' @export
#' @details The parameters (par) has to be present in the observed data frame via ssn_get_data(n, name = "obs"). More details of the argument par can be found in the additive.function() from SSN .
#' @examples
#' \donttest{
#' #require("SSN2")
#' #path <- system.file("extdata/clearwater.ssn", package = "SSNbayes")
#' #ssn <- SSN2::ssn_import(path, predpts = "preds", overwrite  = TRUE)
#' #t.df <- collapse(ssn, par = 'afvArea')}


collapse <- function(ssn, par = 'afvArea'){
  slot <- NULL
  df_all <- NULL
  line_id <- NULL

  df0 <- ssn[[1]] %>% st_union
  df0 <- df0[[1]]

  for (i in 1:length(df0)){
    df <- data.frame(df0[i])
    df$slot <- i #ssn@lines[[i]]@ID

    vec <- as.data.frame(ssn[[1]][i, par])

    df$computed_afv <- vec[,1]#ssn[[2]][i, par]

    df$line_id <- as.numeric(as.character(df$slot))
    df_all<- rbind(df, df_all)
    }
  df_all <-  dplyr::arrange(df_all, line_id)
  #df_all$slot <- NULL
  names(df_all)[names(df_all) == 'computed_afv'] <- par
  df_all
}





#' Creates a list containing the stream distances and weights
#'
#' @param path Path to the files
#' @param net (optional) A network from the SSN2 object
#' @param addfunccol (optional) A parameter to compute the spatial weights
#' @return A list of matrices
#' @importFrom dplyr mutate %>% distinct left_join case_when
#' @importFrom plyr .
#' @importFrom rstan stan
#' @importFrom stats dist
#' @importFrom SSN2 ssn_get_data
#' @importFrom sf st_coordinates
#' @importFrom methods is
#' @export
#' @examples
#' \donttest{
#' path <- system.file("extdata/clearwater.ssn", package = "SSNbayes")
#' mat_all <- dist_weight_mat(path, net = 2, addfunccol='afvArea')
#' }


dist_weight_mat <- function(path = path, net = 1, addfunccol='addfunccol'){
  pid <- NULL
  n  <- SSN2::ssn_import(path, overwrite  = TRUE)
  obs_data <- SSN2::ssn_get_data(n, name = "obs")

  if(is(n) == 'SSN'){
    xy <- st_coordinates(obs_data)
    colnames(xy) <- c('NEAR_X', 'NEAR_Y')
    obs_data <- cbind(obs_data, xy)
  }
  obs_data <- obs_data %>% as.data.frame()


  # creating distance matrices
  D <- readRDS(paste0(path, '/distance/obs/dist.net', net, '.RData')) # distance between observations

  # total distance
  H <- D + base::t(D)

  obs_data <- dplyr::filter(obs_data, pid %in% colnames(H)) %>% as.data.frame()
  # NB replace here by the variable used for spatial weights
  afv <- obs_data[c('locID', addfunccol)] %>% distinct()

  # codes from SSN glmssn
  nsofar <- 0
  dist.junc <- matrix(0, nrow = length(afv[, 1]), ncol = length(afv[,1]))
  distmat <- D
  ni <- length(distmat[1,])
  ordpi <- order(as.numeric(rownames(distmat)))
  dist.junc[(nsofar + 1):(nsofar + ni), (nsofar + 1):(nsofar + ni)] <-
    distmat[ordpi, ordpi, drop = FALSE]
  b.mat <- pmin(dist.junc, base::t(dist.junc))
  dist.hydro <- as.matrix(dist.junc + base::t(dist.junc))
  flow.con.mat <- 1 - (b.mat > 0) * 1

  n.all <- ni
  # weights matrix
  w.matrix <- sqrt(pmin(outer(afv[, addfunccol],rep(1, times = n.all)),
                        base::t(outer(afv[, addfunccol],rep(1, times = n.all) ))) /
                     pmax(outer(afv[, addfunccol],rep(1, times = n.all)),
                          base::t(outer(afv[, addfunccol], rep(1, times = n.all))))) *
    flow.con.mat

  # Euclidean distance

  #obs_data_coord <- data.frame(n@obspoints@SSNPoints[[1]]@point.coords)
  #obs_data_coord$locID <- factor(1:nrow(obs_data_coord))
  #obs_data <- obs_data %>% left_join(obs_data_coord, by = c('locID'))
  obs_data$point <- 'Obs'

  obs_data$coords.x1 <- obs_data$NEAR_X
  obs_data$coords.x2 <- obs_data$NEAR_Y

  #coor <- n@obspoints@SSNPoints[[1]]@point.coords
  coor <- obs_data[, c('NEAR_X', 'NEAR_Y')] #NB: check
  e <- coor %>%
    dist(., method = "euclidean", diag = FALSE, upper = FALSE) %>% as.matrix()

  list(e = e, D = D, H = H, w.matrix = w.matrix, flow.con.mat = flow.con.mat)
}





#' Creates a list of distances and weights between observed and prediction sites
#'
#' @param path Path with the name of the SpatialStreamNetwork object
#' @param net (optional) A network from the SpatialStreamNetwork object
#' @param addfunccol (optional) A parameter to compute the spatial weights
#' @return A list of matrices
#' @importFrom dplyr mutate %>% distinct left_join case_when bind_rows
#' @importFrom plyr .
#' @importFrom rstan stan
#' @importFrom stats dist
#' @importFrom SSN2 ssn_get_data
#' @importFrom methods is
#' @export
#' @description The output matrices are symmetric except the hydrologic distance matrix D.
#' @examples
#' \dontrun{
#' path <- system.file("extdata/clearwater.ssn", package = "SSNbayes")
#' mat_all_pred <- dist_weight_mat_preds(path, net = 2, addfunccol='afvArea')}


dist_weight_mat_preds <- function(path = path, net = 1, addfunccol = 'addfunccol'){
  netID <- NULL
  locID <- NULL
  pid <- NULL

  #Some weird error from SSN
  #n <- SSN2::ssn_import(path = path, predpts = "preds", overwrite  = TRUE)
  t <- try(n <- SSN2::ssn_import(path, predpts = "preds", overwrite  = TRUE),
           silent = TRUE)
  if("try-error" %in% class(t)){n <- SSN2::ssn_import(paste0(getwd(), path), predpts = "preds", overwrite  = TRUE)}


  obs_data <- SSN2::ssn_get_data(n, "obs")
  pred_data <- SSN2::ssn_get_data(n, name = "preds")


  obs_data <- dplyr::filter(obs_data, netID == net ) %>% arrange(locID)
  pred_data <- dplyr::filter(pred_data, netID == net ) %>% arrange(locID)

  obs_data$locID_backup <- obs_data$locID
  pred_data$locID_backup <- pred_data$locID

  if(is(n) == 'SSN'){
    xy <- st_coordinates(obs_data)
    colnames(xy) <- c('NEAR_X', 'NEAR_Y')
    obs_data <- cbind(obs_data, xy)

    xy2 <- st_coordinates(pred_data)
    colnames(xy2) <- c('NEAR_X', 'NEAR_Y')
    pred_data <- cbind(pred_data, xy2)
  }
  obs_data <- obs_data %>% as.data.frame()
  pred_data <- pred_data %>% as.data.frame()

  if(file.exists(paste0(path, '/distance/preds')) == FALSE) stop("no distance matrix available between predictions. Please, use SSN2::ssn_create_distmat(n, predpts = 'preds', overwrite =TRUE, among_predpts  = TRUE)")

  # creating distance matrices
  doo <- readRDS(paste0(path, '/distance/obs/dist.net', net, '.RData'))

  if(file.exists(paste0(path, '/distance/preds/dist.net', net, '.a.RData')) == FALSE) stop("no distance matrix available between observations and predictions. Please, use SSN2::ssn_create_distmat(n, predpts = 'preds', overwrite =TRUE, among_predpts  = TRUE)")
  dop <- readRDS(paste0(path, '/distance/preds/dist.net', net, '.a.RData')) # distance between observations
  dpo <- readRDS(paste0(path, '/distance/preds/dist.net', net, '.b.RData'))
  dpp <- readRDS(paste0(path, '/distance/preds/dist.net', net, '.RData'))

  D <- rbind(cbind(doo, dop), cbind(dpo, dpp))


  H <- D + base::t(D)


  pred_data <- dplyr::filter(pred_data, pid %in% colnames(H))
  # NB replace here by the variable used for spatial weights
  afv1 <-  obs_data[c('locID', addfunccol)] %>% distinct() %>% as.data.frame()
  afv2 <- pred_data[c('locID', addfunccol)] %>% distinct() %>% as.data.frame()

  afv <- dplyr::bind_rows(afv1, afv2)

  # codes from SSN::glmssn
  nsofar <- 0
  dist.junc <- matrix(0, nrow = length(afv[, 1]), ncol = length(afv[,1]))
  distmat <- D
  ni <- length(distmat[1,])
  ordpi <- order(as.numeric(rownames(distmat)))
  dist.junc[(nsofar + 1):(nsofar + ni), (nsofar + 1):(nsofar + ni)] <-
    distmat[ordpi, ordpi, drop = FALSE]
  b.mat <- pmin(dist.junc, base::t(dist.junc))
  colnames(b.mat) <- colnames(D)
  rownames(b.mat) <- rownames(D)
  dist.hydro <- as.matrix(dist.junc + base::t(dist.junc))
  flow.con.mat <- 1 - (b.mat > 0) * 1
  colnames(b.mat) <- colnames(D)
  rownames(b.mat) <- rownames(D)

  n.all <- ni
  # weights matrix
  w.matrix <- sqrt(pmin(outer(afv[, addfunccol],rep(1, times = n.all)),
                        base::t(outer(afv[, addfunccol],rep(1, times = n.all) ))) /
                     pmax(outer(afv[, addfunccol],rep(1, times = n.all)),
                          base::t(outer(afv[, addfunccol], rep(1, times = n.all))))) *
    flow.con.mat

  # Euclidean distance

  obs_data_coord <- obs_data[, c('NEAR_X', 'NEAR_Y', 'locID')] #NB: check

  #obs_data_coord <- data.frame(n@obspoints@SSNPoints[[1]]@point.coords)
  obs_data_coord$locID <- factor(1:nrow(obs_data_coord))
  obs_data_coord$locID <- as.numeric(as.character(obs_data_coord$locID))

  obs_data$locID <- as.numeric(factor(obs_data$locID))
  #obs_data <- obs_data %>% left_join(obs_data_coord, by = c('locID'))
  obs_data$point <- 'Obs'

  pred_data_coord <- pred_data[, c('NEAR_X', 'NEAR_Y', 'locID')] #data.frame(n@predpoints@SSNPoints[[1]]@point.coords)
  pred_data_coord$locID <- factor( (max(as.numeric(as.character(obs_data_coord$locID))) + 1)
                                   :(nrow(pred_data_coord)+ (max(as.numeric(as.character(obs_data_coord$locID))) )))
  pred_data_coord$locID <- as.numeric(as.character(pred_data_coord$locID))

  pred_data$locID <- as.numeric(factor(pred_data$locID)) + max(obs_data$locID )
  #pred_data <- pred_data %>% left_join(pred_data_coord, by = c('locID'))
  pred_data$point <- 'preds'

  #all_data <- rbind(obs_data[c('coords.x1', 'coords.x2')],
  #                  pred_data[c('coords.x1', 'coords.x2')])

  all_data <- dplyr::bind_rows(obs_data[c('NEAR_X', 'NEAR_Y')] %>% as.data.frame(),
                               pred_data[c('NEAR_X', 'NEAR_Y')] %>% as.data.frame())

  all_data$coords.x1 <- all_data$NEAR_X
  all_data$coords.x2 <- all_data$NEAR_Y

  e <- all_data %>%
    dplyr::select('coords.x1', 'coords.x2') %>%
    dist(., method = "euclidean", diag = FALSE, upper = FALSE) %>% as.matrix()
  colnames(e) <- colnames(D)
  rownames(e) <- rownames(D)

  list(e = e, D = D, H = H, w.matrix = w.matrix, flow.con.mat = flow.con.mat)
}







#' A simple modeling function using a formula and data
#'
#' @param formula A formula as in lm()
#' @param data A data.frame containing the elements specified in the formula
#' @return A list of matrices
#' @importFrom stats model.matrix model.response
#' @export
#' @author Jay ver Hoef
#' @examples
#' options(na.action='na.pass')
#' data("iris")
#' out_list = mylm(formula = Petal.Length ~ Sepal.Length + Sepal.Width, data = iris)


mylm <- function(formula, data) {
  # get response as a vector
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  y <- as.vector(model.response(mf, "numeric"))
  # create design matrix
  X <- model.matrix(formula, data)
  # return a list of response vector and design matrix
  return(list(y = y, X = X))
}





#' Fits a mixed linear regression model using Stan
#'
#' It requires the same number of observation/locations per day.
#' It requires location id (locID) and points id (pid).
#' The locID are unique for each site.
#' The pid is unique for each observation.
#' Missing values are allowed in the response but not in the covariates.
#'
#' @param path Path with the name of the SpatialStreamNetwork object
#' @param formula A formula as in lm()
#' @param data A long data frame containing the locations, dates, covariates and the response variable. It has to have the locID and date. No missing values are allowed in the covariates.
#' The order in this data.fame MUST be: spatial locations (1 to S) at time t=1, then locations (1 to S) at t=2 and so on.
#' @param space_method A list defining if use or not of an SSN object and the spatial correlation structure. The second element is the spatial covariance structure. A 3rd element is a list with the lon and lat for Euclidean distance models.
#' @param time_method A list specifying the temporal structure (ar = Autorregressive; var = Vector autorregression) and coumn in the data with the time variable.
#' @param iter Number of iterations
#' @param warmup Warm up samples
#' @param chains Number of chains
#' @param refresh Sampler refreshing rate
#' @param net The network id (optional). Used when the SSN object contains multiple networks.
#' @param addfunccol Variable to compute the additive function. Used to compute the spatial weights.
#' @param loglik Logic parameter denoting if the loglik will be computed by the model.
#' @param ppd Produce the posterior predictive distribution
#' @param seed (optional) A seed for reproducibility
#' @return A list with the model fit
#' @details Missing values are not allowed in the covariates and they must be imputed before using ssnbayes(). Many options can be found in https://cran.r-project.org/web/views/MissingData.html
#' The pid in the data has to be consecutive from 1 to the number of observations.
#' Users can use the SpatialStreamNetwork created with the SSN package. This will provide the spatial stream information used to compute covariance matrices. If that is the case, the data has
#' to have point ids (pid) matching the ones in SSN distance matrices, so that a mapping can occur.
#' @return It returns a ssnbayes object (similar to stan returns). It includes the formula used to fit the model. The output can be transformed into the stanfit class using class(fits) <- c("stanfit").
#' @export
#' @importFrom dplyr mutate %>% distinct left_join case_when
#' @importFrom plyr .
#' @importFrom rstan stan
#' @importFrom stats dist
#' @author Edgar Santos-Fernandez
#' @examples
#'\dontrun{
#'#options(mc.cores = parallel::detectCores())
#'# Import SpatialStreamNetwork object
#'#path <- system.file("extdata/clearwater.ssn", package = "SSNbayes")
#'#n <- SSN2::ssn_import(path, predpts = "preds", overwrite  = TRUE)
#'## Imports a data.frame containing observations and covariates
#'#clear <- readRDS(system.file("extdata/clear_obs.RDS", package = "SSNbayes"))
#'#fit_ar <- ssnbayes(formula = y ~ SLOPE + elev + h2o_area + air_temp + sin + cos,
#'#                   data = clear,
#'#                   path = path,
#'#                   time_method = list("ar", "date"),
#'#                   space_method = list('use_ssn', c("Exponential.taildown")),
#'#                   iter = 2000,
#'#                   warmup = 1000,
#'#                   chains = 3,
#'#                   net = 2, # second network on the ssn object
#'#                   addfunccol='afvArea')

#' #space_method options examples
#' #use list('no_ssn', 'Exponential.Euclid', c('lon', 'lat')) if no ssn object is available
#'}

ssnbayes <- function(formula = formula,
                     data = data,
                     path = path,
                     time_method = time_method, # list("ar", "date")
                     space_method = space_method, #list('use_ssn', 'Exponential.tailup'),
                     iter = 3000,
                     warmup = 1500,
                     chains = 3,
                     refresh = max(iter/100, 1),
                     net = 1,
                     addfunccol = addfunccol,
                     loglik = FALSE,
                     ppd = FALSE,
                     seed = seed
){

  # checks
  if(missing(time_method)){
    stop("Need to define the method (ar or var) and the column associated with time")
  }

  if(length(time_method) == 1){
    stop("Need to specify the column in the the data with the time variable")
  }

  time_points <- time_method[[2]]

  #if('date' %in% names(data) == FALSE) stop("There is no column date on the data. Please, set a column called date with the time")
  if('locID' %in% names(data) == FALSE) stop("There is no column locID on the data. Please, set a column called locID with the observation locations")

  if(missing(seed)) seed <- sample(1:1E6,1,replace=TRUE)

  if(!missing(space_method)){
    print('using SSN object')
    if(space_method[[1]] == 'use_ssn'){
      ssn_object <- TRUE


      if(length(space_method) > 1){
        if(space_method[[2]] %in% c("Exponential.tailup", "LinearSill.tailup" , "Spherical.tailup" ,
                                    "Exponential.taildown" ,"LinearSill.taildown" ,"Spherical.taildown",
                                    "Exponential.Euclid") == FALSE) {stop("Need to specify one or more of the following covariance matrices: Exponential.tailup, LinearSill.tailup , Spherical.tailup ,
		Exponential.taildown, LinearSill.taildown, Spherical.taildown or Exponential.Euclid")}
        CorModels <- space_method[[2]]
      }
      if(length(space_method) == 1){
        CorModels <- "Exponential.tailup"
        print('using an Exponential.tailup model')
      }

    }
    if(space_method[[1]] == 'no_ssn'){
      print('no SSN object defined')
      ssn_object <- FALSE
      if(space_method[[2]] %in% c("Exponential.Euclid") == FALSE) {stop("Need to specify Exponential.Euclid")}
      # when using Euclidean distance, need to specify the columns with lon and lat.
      if(length(space_method) < 3){ stop("Please, specify the columns in the data frame with the longitude and latitude (c('lon', 'lat'))") }

      data$lon <- data[,names(data) == space_method[[3]][1]]
      data$lat <- data[,names(data) == space_method[[3]][2]]
      CorModels <- space_method[[2]]
    }


  }

  if(missing(space_method)) {space_method <- 'no_ssn'; ssn_object <- FALSE; CorModels <- "Exponential.Euclid" }# if missing use Euclidean distance

 # if(loglik == T & ppd == T) {stop("Need to specify either loglik = T or ppd = T; or both = F ")}


  # Cov
  cor_tu <- case_when(CorModels == "Exponential.tailup" ~ 1,
                      CorModels == "LinearSill.tailup" ~ 2,
                      CorModels == "Spherical.tailup" ~ 3,
                      TRUE ~ 5)
  cor_tu <- sort(cor_tu)[1]

  cor_td <- case_when(CorModels == "Exponential.taildown" ~ 1,
                      CorModels == "LinearSill.taildown" ~ 2,
                      CorModels == "Spherical.taildown" ~ 3,
                      TRUE ~ 5)
  cor_td <- sort(cor_td)[1]


  cor_ed <- case_when(CorModels == "Exponential.Euclid" ~ 1,
                      TRUE ~ 5)
  #CorModels == "Spherical.Euclid" ~ 2, #NB: to be implemented
  #CorModels == "Gaussian.Euclid" ~ 3)
  cor_ed <- sort(cor_ed)[1]

  cor_re <- case_when(CorModels == "RE1" ~ 1,
                      TRUE ~ 5)
  cor_re <- sort(cor_re)[1]



  data_com <-  'data {
    int<lower=1> N;
    int<lower=1> K;
    int<lower=1> T;
    matrix[N,K] X[T] ; // real X[N,K,T]; //

    int<lower = 0> N_y_obs; // number observed values
    int<lower = 0> N_y_mis; // number missing values

    int<lower = 1> i_y_obs[N_y_obs] ;  //[N_y_obs,T]
    int<lower = 1> i_y_mis[N_y_mis] ;  // N_y_mis,T]

    vector[N_y_obs] y_obs;    //matrix[N_y_obs,1] y_obs[T];

    matrix[N, N] W ; // spatial weights
    matrix[N, N] h ; // total hydrological dist
    matrix[N, N] I ; // diag matrix

    matrix[N, N] D ; // downstream hydrological dist matrix
    matrix[N, N] flow_con_mat; // flow conected matrix

    matrix[N, N] e ; // Euclidean dist mat
    real<lower=1> alpha_max ;

  }'


  param_com <- '
  parameters {
    vector[K] beta;
    real<lower=0> sigma_nug;

    vector[N_y_mis] y_mis;//declaring the missing y
  '


  param_phi_ar <- '
    real <lower=-1, upper = 1> phi; // NB
  '
  param_phi_var <- '
    vector<lower=-1, upper = 1> [N] phi  ; // vector of autoregresion pars
  '



  param_tu <- '
    real<lower=0> sigma_tu;
    real<lower=0> alpha_tu;
'

  param_td <- '
    real<lower=0> sigma_td; // sd of tail-down
    real<lower=0> alpha_td; // range of the tail-down model
'

  param_ed <- '
    real<lower=0> sigma_ed;
    real<lower=0> alpha_ed; // range of the Euclidean dist model
'

  param_re <- '
    real<lower=0> sigma_RE1;
'

  tparam_com <- '
  transformed parameters {
    vector[N * T] y;
    vector[N] Y[T];

    vector[N] epsilon[T]; // error term
    vector[N] mu[T]; // mean

   real<lower=0> var_nug; // nugget

   matrix[N, N] C_tu; //tail-up cov
   matrix[N, N] C1; //tail-up cov
   matrix[N, N] Ind; //tail-up indicator

   matrix[N, N] C_td; //tail-down cov
   matrix[N, N] Ind2; //tail-down indicator
   matrix[2,1] iji;

   matrix[N, N] C_ed ;// Euclidean cov

   matrix[N, N] C_re ;// random effect cov
   matrix[N, N] RE1; // random effect 1

   '

  tparam_tu <- '
   // tail up exponential
   real<lower=0> var_tu; // parsil tail-down
  '

  tparam_td <- '
    real<lower=0> var_td; // parsil tail-down
 '

  tparam_ed <- '
    real<lower=0> var_ed; //  Euclidean dist var
'

  tparam_re <- '
    real<lower=0> var_RE1; // Random effect 1
'


  tparam_com2 <- '
    y[i_y_obs] = y_obs;
    y[i_y_mis] = y_mis;
    for (t in 1:T){
      Y[t] = y[((t - 1) * N + 1):(t * N)];
    }

    var_nug = sigma_nug ^ 2; // variance nugget
        mu[1] = X[1] * beta;
    epsilon[1] = Y[1] - mu[1];

'

  tparam_com_ar <- '
        for (t in 2:T){
        mu[t] = X[t] * beta;
        epsilon[t] = Y[t] - mu[t];
        mu[t] = mu[t] + phi * epsilon[t-1]; //
    }

  '

  tparam_com_var <- '
        for (t in 2:T){
        mu[t] = X[t] * beta;
        epsilon[t] = Y[t] - mu[t];
        mu[t] = mu[t] + phi .* epsilon[t-1]; // element wise mult two vectors
    }
  '


  tparam_tu2_exp <- '
    // tail up exponential
    var_tu = sigma_tu ^ 2; // variance tail-up
      C1 = var_tu * exp(- 3 * h / alpha_tu); // tail up exponential model
    C_tu = C1 .* W; // Hadamard (element-wise) product
'

  tparam_tu2_lin <- '
     //Tail-up linear-with-sill model
     var_tu = sigma_tu ^ 2; // variance tail-up
      for (i in 1:N) {
        for (j in 1:N) {
        Ind[i,j] = (h[i,j] / alpha_tu) <= 1 ? 1 : 0; // indicator
        }
      }
      C1 = var_tu * (1 - (h / alpha_tu)) .* Ind ; //Tail-up linear-with-sill model
    C_tu = C1 .* W; // Hadamard (element-wise) product
'

  tparam_tu2_sph <- '
      // Tail-up spherical model
      var_tu = sigma_tu ^ 2; // variance tail-up
      for (i in 1:N) {// Tail-up spherical model
        for (j in 1:N) {
          Ind[i,j] = (h[i,j] / alpha_tu) <= 1 ? 1 : 0; // indicator
        }
      }
    C1 = var_tu * (1 - (1.5 * h / alpha_tu) + (h .* h .* h / (2 * alpha_tu ^ 3))) .* Ind ; // Tail-up spherical model
    C_tu = C1 .* W; // Hadamard (element-wise) product
'
  #tail-up models end


  #tail-down models start

  # tail-down exponential
  tparam_td2_exp <- '
    var_td= sigma_td ^ 2; // variance tail-down
     	for (i in 1:N) {// Tail-down exponential model
          for (j in 1:N) {
            if(flow_con_mat[i,j] == 1){ // if points are flow connected
               C_td[i,j] = var_td * exp(- 3 * h[i,j] / alpha_td);
            }
            else{// if points are flow unconnected
              C_td[i,j] = var_td * exp(- 3 * (D[i,j] + D[j,i]) / alpha_td);
            }
          }
        }

'

  #Tail-down linear-with-sill model
  tparam_td2_lin <- '
    var_td= sigma_td ^ 2; // variance tail-down
        for (i in 1:N) {// Tail-down linear-with-sill model
          for (j in 1:N) {
            if(flow_con_mat[i,j] == 1){ // if points are flow connected
              Ind2[i,j] = (h[i,j] / alpha_td) <= 1 ? 1 : 0; // indicator
              C_td[i,j] = var_td * (1 - (h[i,j] / alpha_td)) .* Ind2[i,j] ; //Tail-up linear-with-sill model
            }
            else{// if points are flow unconnected
              iji[1,1] = D[i,j];
              iji [2,1] = D[j,i];
              Ind2[i,j] = (max(iji) / alpha_td) <= 1 ? 1 : 0; // indicator
              C_td[i,j] = var_td * (1 - (max(iji) / alpha_td)) * Ind2[i,j] ;
            }
          }
        }

  '

  #tail-down spherical model
  tparam_td2_sph <- '
    var_td= sigma_td ^ 2; // variance tail-down
        for (i in 1:N) {// tail-down spherical model
          for (j in 1:N) {
            if(flow_con_mat[i,j] == 1){ // if points are flow connected
              Ind2[i,j] = (h[i,j] / alpha_td) <= 1 ? 1 : 0; // indicator
              C_td[i,j] = var_td * (1 - (1.5 * h[i,j] / alpha_td) + ( (h[i,j] ^ 3) / (2 * alpha_td ^ 3))) * Ind2[i,j];
            }
            else{// if points are flow unconnected
              iji[1,1] = D[i,j];
              iji [2,1] = D[j,i];
              Ind2[i,j] = (max(iji) / alpha_td) <= 1 ? 1 : 0; // indicator
              C_td[i,j] = var_td * (1 - (1.5 * min(iji) / alpha_td) + ( max(iji)/(2 * alpha_td)   )) * (1 - (max(iji) / alpha_td) ) ^ 2 * Ind2[i,j];
            }
          }
        }

'

  #tail-down models end


  tparam_ed2 <- '
	  //Euclidean distance models start
    var_ed = sigma_ed ^ 2; // var Euclidean dist
      C_ed = var_ed * exp(- 3 * e / alpha_ed); // exponential model
    //Euclidean distance models end
'

  tparam_re2 <- '
    // random effect
    var_RE1 = sigma_RE1 ^ 2;
      C_re = var_RE1 * RE1;

'

  model_com <- '
  model {
    for (t in 1:T){
      target += multi_normal_cholesky_lpdf(Y[t] | mu[t], cholesky_decompose(C_tu + C_td + C_re + C_ed + var_nug * I + 1e-6) );
    }

    sigma_nug ~ uniform(0,50); // cauchy(0,1) prior nugget effect
    phi ~ uniform(-1, 1); // or can use phi ~ normal(0.5,0.3); //NB informative

'

  model_tu <- '
    sigma_tu ~ uniform(0,100);  // or cauchy(0,2) prior sd  tail-up model
    alpha_tu ~ uniform(0, alpha_max);
'

  model_td <- '
    sigma_td ~ uniform(0,100); // sd tail-down
    alpha_td ~ uniform(0, alpha_max);
'

  model_ed <- '
    sigma_ed ~ uniform(0,100); // sd Euclidean dist
    alpha_ed ~ uniform(0, alpha_max); // Euclidean dist range
'

  model_re <- '
    sigma_RE1 ~ uniform(0,5);
'

  gen_quant_ll <- '
  generated quantities {
   vector[T] log_lik;
     for (t in 1:T){
       log_lik[t] = multi_normal_cholesky_lpdf(Y[t]|mu[t],
        cholesky_decompose(C_tu + C_td + C_re + C_ed + var_nug * I + 1e-6) );
      }
  }
  '

  gen_quant_ppd <- '
  generated quantities {
   vector[N] y_pred[T];
     for (t in 1:T){
     y_pred[t] = multi_normal_cholesky_rng( mu[t], cholesky_decompose(C_tu + C_td + C_re + C_ed + var_nug * I + 1e-6) );
      }
  }
  '


  ssn_ar <- paste(
    data_com,

    param_com,

    if(cor_tu %in% 1:3) param_tu,
    if(cor_td %in% 1:3) param_td,
    if(cor_ed %in% 1:3) param_ed,
    if(cor_re %in% 1:3) param_re,

    if(time_method[[1]] == 'ar') param_phi_ar,
    if(time_method[[1]] == 'var') param_phi_var,

    '}',

    tparam_com,
    if(cor_tu %in% 1:3)tparam_tu,
    if(cor_td %in% 1:3)tparam_td,
    if(cor_ed %in% 1:3)tparam_ed,
    if(cor_re %in% 1:3) tparam_re,
    tparam_com2,


    if(time_method[[1]] == 'ar') tparam_com_ar,
    if(time_method[[1]] == 'var') tparam_com_var,



    case_when(cor_tu == 1 ~ tparam_tu2_exp,
              cor_tu == 2 ~ tparam_tu2_lin,
              cor_tu == 3 ~ tparam_tu2_sph,
              cor_tu >= 4 | cor_tu <= 0 ~ 'C_tu = rep_matrix(0, N, N);'),


    case_when(cor_td == 1 ~ tparam_td2_exp,
              cor_td == 2 ~ tparam_td2_lin,
              cor_td == 3 ~ tparam_td2_sph,
              cor_td >= 4 | cor_td <= 0 ~ 'C_td = rep_matrix(0, N, N);'),

    case_when(cor_ed == 1 ~ tparam_ed2,
              cor_ed >= 2 | cor_ed <= 0 ~ 'C_ed = rep_matrix(0, N, N);'),

    case_when(cor_re == 1 ~ tparam_re2,
              cor_re >= 2 | cor_re <= 0 ~ 'C_re = rep_matrix(0, N, N);'),

    '}',
    model_com,
    if(cor_tu %in% 1:3)model_tu,
    if(cor_td %in% 1:3)model_td,
    if(cor_ed %in% 1:3)model_ed,
    if(cor_re %in% 1:3) model_re,
    '}',

    if(loglik == TRUE) gen_quant_ll,
    if(ppd == TRUE) gen_quant_ppd
  )

  `%notin%` <- Negate(`%in%`)

  pars <- c(
    case_when(cor_tu %in% 1:3 ~ c('var_tu', 'alpha_tu'),
              cor_tu %notin% 1:3 ~ ""),

    case_when(cor_td %in% 1:3 ~ c('var_td', 'alpha_td'),
              cor_td %notin% 1:3 ~ ""),

    case_when(cor_ed %in% 1:3 ~ c('var_ed', 'alpha_ed'),
              cor_ed %notin% 1:3 ~ ""),

    case_when(cor_re %in% 1:3 ~ c('var_re', 'alpha_re'),
              cor_re %notin% 1:3 ~ ""),

    if(loglik == TRUE) 'log_lik',
    if(ppd == TRUE) 'y_pred',

    'var_nug',
    'beta',
    'phi',
    'y'
  )
  pars <- pars[pars != '']


  # data part
  old <- options()        # old options
  on.exit(options(old)) 	# reset once exit the function

  options(na.action='na.pass') # to preserve the NAs
  out_list <- mylm(formula = formula, data = data) # produces the design matrix

  response <- out_list$y # response variable
  design_matrix <- out_list$X # design matrix

  obs_data <- data

  ndays <- length(unique(obs_data[, names(obs_data) %in% time_points] ))


  N <- nrow(obs_data)/ndays #nobs


  nobs <- nrow(obs_data)/ndays #nobs

  obs_data$date_num <- as.numeric(factor(obs_data[, names(obs_data) %in% time_points]	))

  resp_var_name <- gsub("[^[:alnum:]]", " ", formula[2])
  obs_data$y <- obs_data[,names(obs_data) %in% resp_var_name]


  # array structure
  X <- design_matrix #cbind(1,obs_data[, c("X1", "X2", "X3")]) # design matrix

    # NB: this array order is Stan specific
  Xarray <- aperm(array( c(X), dim=c(N, ndays, ncol(X)) ),c(2, 1, 3))

  y_obs <- response[!is.na(response)]

  # index for observed values
  i_y_obs <- obs_data[!is.na(obs_data$y),]$pid

  # index for missing values
  i_y_mis <- obs_data[is.na(obs_data$y),]$pid

  if(ssn_object == TRUE){ # the ssn object exist?
    mat_all <- dist_weight_mat(path = path, net = net, addfunccol = addfunccol)
  }

  if(ssn_object == FALSE){ # the ssn object does not exist- purely spatial

    first_date <- unique(obs_data[, names(obs_data) %in% time_points])[1]

    di <- dist(obs_data[obs_data$date == first_date, c('lon', 'lat')], #data$date == 1
               method = "euclidean",
               diag = FALSE,
               upper = FALSE) %>% as.matrix()
    mat_all <-  list(e = di, D = di, H = di, w.matrix = di, flow.con.mat = di)
  }


  data_list <- list(N = N, # obs + preds  points
                    T = ndays, # time points
                    K = ncol(X),  # ncol of design matrix
                    y_obs = y_obs,# y values in the obs df

                    N_y_obs = length(i_y_obs),  #nrow(i_y_obs) numb obs points
                    N_y_mis = length(i_y_mis), #nrow(i_y_mis) numb preds points

                    i_y_obs = i_y_obs, # index of obs points
                    i_y_mis = i_y_mis, # index of preds points

                    X = Xarray, # design matrix
                    mat_all = mat_all,
                    alpha_max = 4 * max(mat_all$H) ) # a list with all the distance/weights matrices


  data_list$e = data_list$mat_all$e #Euclidean dist
  #for tail-up
  data_list$h = data_list$mat_all$H # total stream distance
  data_list$W = data_list$mat_all$w.matrix # spatial weights

  #for tail-down
  data_list$flow_con_mat = data_list$mat_all$flow.con.mat #flow connected matrix
  data_list$D = data_list$mat_all$D #downstream hydro distance matrix

  #RE1 = RE1mm # random effect matrix

  data_list$I = diag(1, nrow(data_list$W), nrow(data_list$W))  # diagonal matrix

  ini <- function(){list(var_nug =  .1)}

  fit <- rstan::stan(model_code = ssn_ar,
                     model_name = "ssn_ar",
                     data = data_list,
                     pars = pars,
                     iter = iter,
                     warmup = warmup,
                     init = ini,
                     chains = chains,
                     verbose = FALSE,
                     seed = seed,
                     refresh = refresh
  )
  attributes(fit)$formula <- formula

  class(fit) <- 'ssnbayes'

  fit
}





#' Performs spatio-temporal prediction in R using an ssnbayes object from a fitted model.
#'
#' It will take an observed and a prediction data frame.
#' It requires the same number of observation/locations per day.
#' It requires location id (locID) and points id (pid).
#' The locID are unique for each site.
#' The pid is unique for each observation.
#' Missing values are allowed in the response but not in the covariates.
#'
#' @param object A stanfit object returned from ssnbayes
#' @param ... Other parameters
#' @param path Path with the name of the SpatialStreamNetwork object
#' @param obs_data The observed data frame
#' @param pred_data The predicted data frame
#' @param net (optional) Network from the SSN object
#' @param nsamples The number of samples to draw from the posterior distributions. (nsamples <= iter)
#' @param addfunccol The variable used for spatial weights
#' @param chunk_size (optional) the number of locID to make prediction from
#' @param locID_pred (optional) the location id for the predictions. Used when the number of pred locations is large.
#' @param seed (optional) A seed for reproducibility
#' @return A data frame with the location (locID), time point (date), plus the MCMC draws from the posterior from 1 to the number of iterations.
#' The locID0 column is an internal consecutive location ID (locID) produced in the predictions, starting at max(locID(observed data)) + 1. It is used internally in the way predictions are made in chunks.
#' @details The returned data frame is melted to produce a long dataset. See examples.
#' Currently, the predict() function produces predictions for normal random variables. However, this can be easily transformed in to counts (Poisson distributed) and presence/absence (binomial distributed).
#' @export
#' @importFrom dplyr mutate %>% distinct left_join case_when
#' @importFrom plyr .
#' @importFrom rstan stan
#' @importFrom stats dist
#' @author Edgar Santos-Fernandez
#' @examples
#' \donttest{
#'#require('SSNdata')
#'#clear_preds <- readRDS(system.file("extdata/clear_preds.RDS", package = "SSNdata"))
#'#clear_preds$y <- NA
#'#pred <- predict(object = fit_ar,
#'#                 path = path,
#'#                 obs_data = clear,
#'#                 pred_data = clear_preds,
#'#                 net = 2,
#'#                 nsamples = 100, # numb of samples from the posterior
#'#                 addfunccol = 'afvArea', # var for spatial weights
#'#                 locID_pred = locID_pred,
#'#                 chunk_size = 60)
#'}

predict.ssnbayes <- function(object = object,
                             ...,
                             path = path,
                             obs_data = obs_data,
                             pred_data = pred_data,
                             net = net,
                             nsamples = nsamples, # number of samples to use from the posterior in the stanfit object
                             addfunccol = addfunccol, # variable used for spatial weights
                             locID_pred = locID_pred,
                             chunk_size = chunk_size,
                             seed = seed) {

  stanfit <- object
  formula <- as.formula(attributes(stanfit)$formula)
  obs_resp <- obs_data[,gsub("\\~.*", "", formula)[2]]
  if( any( is.na(obs_resp) )) {stop("Can't have missing values in the response in the observed data. You need to impute them before")}


  out <- pred_ssnbayes(object = object,
                       path = path,
                       obs_data = obs_data,
                       pred_data = pred_data,
                       net = net,
                       nsamples = nsamples, # number of samples to use from the posterior in the stanfit object
                       addfunccol = addfunccol, # variable used for spatial weights
                       locID_pred = locID_pred,
                       chunk_size = chunk_size,
                       seed = seed)
  out
}







#' Performs spatio-temporal prediction in R using an ssnbayes object from a fitted model.
#'
#' It will take an observed and a prediction data frame.
#' It requires the same number of observation/locations per day.
#' It requires location id (locID) and points id (pid).
#' The locID are unique for each site.
#' The pid is unique for each observation.
#' Missing values are allowed in the response but not in the covariates.
#'
#' @param object A stanfit object returned from ssnbayes
#' @param ... Other parameters
#' @param use_osm_ssn Use a supplied osm_ssn instead a SpatialStreamNetwork object.
#' @param osm_ssn The osm_ssn to be used.
# #' @param family A description of the response distribution and link function to be used in the model. Must be one of: 'gaussian', 'binomial', 'poisson".
#' @param path If not using an osm_ssn, path with the name of the SpatialStreamNetwork object.
#' @param obs_data The observed data frame
#' @param pred_data The predicted data frame
#' @param net (optional) Network from the SSN object
#' @param nsamples The number of samples to draw from the posterior distributions. (nsamples <= iter)
#' @param addfunccol If not using an osm_ssn, the variable used for spatial weights
#' @param chunk_size (optional) the number of locID to make prediction from
#' @param locID_pred (optional) the location id for the predictions. Used when the number of pred locations is large.
#' @param seed (optional) A seed for reproducibility
#' @return A data frame with the location (locID), time point (date), plus the MCMC draws from the posterior from 1 to the number of iterations.
#' The locID0 column is an internal consecutive location ID (locID) produced in the predictions, starting at max(locID(observed data)) + 1. It is used internally in the way predictions are made in chunks.
#' @details The returned data frame is melted to produce a long dataset. See examples.
#' @export
#' @importFrom dplyr mutate %>% distinct left_join case_when
#' @importFrom plyr .
#' @importFrom rstan stan
#' @importFrom stats dist
#' @author Edgar Santos-Fernandez
#' @examples
#' \donttest{
#'# require('SSNdata')
#'# require('sf')
#'#
#'# clear <- readRDS(system.file("extdata/clear_obs.RDS", package = "SSNdata"))
#'#
#'# clear_osm_ssn <- generate_osm_ssn(clear, "long", "lat", root_loc = 12, plot_network = TRUE)
#'#
#'# formula = y ~ SLOPE + elev + air_temp + sin + cos
#'#
#'# family = "gaussian"
#'#
#'# # Note - missing data must be imputed before using the predict() function.
#'# # This can be done using:
#'# data_impute <- mtsdi::mnimput(
#'#   formula = temp ~ SLOPE + elev + h2o_area + air_temp + sin + cos,
#'#   dataset = clear,
#'#   eps = 1e-3,
#'#   ts = FALSE,
#'#   method = "glm"
#'# )$filled.dataset
#'#
#'# clear <- left_join(clear, data_impute, by = c("elev", "air_temp", "sin", "cos", "SLOPE"))
#'#
#'# clear$y <- clear$temp.y
#'#
#'# fit_ar <- ssnbayes2(formula = formula,
#'#                     data = clear,
#'#                     osm_ssn = clear_osm_ssn,
#'#                     family = family,
#'#                     time_method = list("ar", "date"),
#'#                     space_method = list('use_osm_ssn', c("Exponential.taildown")),
#'#                     iter = 2000,
#'#                     warmup = 1000,
#'#                     chains = 3,
#'#                     cores = 3)
#'#
#'# # Get the coordinate reference system of the near_X and near_Y variables
#'# data_crs <- system.file("extdata/clearwater.ssn", package = "SSNbayes") %>%
#'#   SSN2::ssn_import(predpts = "preds", overwrite  = TRUE) %>%
#'#   SSN2::ssn_get_data(name = "preds") %>%
#'#   st_crs
#'#
#'# # Load in the predictions, and convert to 4326 CRS
#'# clear_preds <- readRDS(system.file("extdata/clear_preds.RDS", package = "SSNbayes")) %>%
#'#   st_as_sf(coords = c("NEAR_X","NEAR_Y"),
#'#            crs = data_crs) %>%
#'#   st_transform(crs = 4326)
#'#
#'# # Extract the long and lat columns
#'# xy <- st_coordinates(clear_preds)
#'#
#'# colnames(xy) <- c("long", "lat")
#'#
#'# Merge back in with the prediction dataset
#'# clear_preds <- cbind(clear_preds, xy) %>%
#'#   data.frame() %>%
#'#   select(-geometry)
#'#
#'# same_names <- intersect(names(clear), names(clear_preds))
#'#
#'# # Combine all observed and prediction sites into 1 dataframe
#'# all_sites <- rbind(clear %>% select(same_names),
#'#                    data.frame(clear_preds) %>% select(same_names)
#'#                    )
#'#
#'#
#'# clear_preds_osm_ssn <- generate_osm_ssn(sensor_data = all_sites,
#'#                                         lon_name = "long", lat_name = "lat",
#'#                                         root_loc = 12,
#'#                                         plot_network = TRUE,
#'#                                         gen_pred_sites = FALSE)
#'#
#'#
#'#
#'# locs <- clear_preds_osm_ssn$dist_mat_all$e %>% colnames
#'#
#'#
#'# clear_krig <- clear %>%
#'#   filter(locID %in% locs)
#'#
#'# clear_krig_preds <- clear_preds %>%
#'#   filter(locID %in% locs)
#'#
#'#
#'# preds <- predict(object = fit_ar,
#'#                  use_osm_ssn = TRUE,
#'#                  osm_ssn = clear_preds_osm_ssn,
#'#                  obs_data = clear_krig,
#'#                  pred_data = clear_krig_preds,
#'#                  seed = seed,
#'#                  nsamples = 25,
#'#                  chunk_size = length(unique(clear_krig_preds$locID))
#'#                  )
#'#
#'# # Condense data to posterior point estimates
#'# ys <- reshape2::melt(preds, id.vars = c('locID0', 'locID', 'date'), value.name ='y')
#'# ys$iter <- gsub("[^0-9.-]", "", ys$variable)
#'# ys$variable <- NULL
#'#
#'# ys <- data.frame(ys) %>% dplyr::group_by(date, locID, locID0) %>%
#'#   dplyr::summarise("sd" = sd(y, na.rm=T),
#'#                    "y_pred" = mean(y, na.rm=T))
#'#
#'# ys <- dplyr::arrange(ys, locID)
#'#'}

predict.ssnbayes2 <- function(object = object,
                              ...,
                              use_osm_ssn = TRUE,
                              osm_ssn = osm_ssn,
                              # family = family,
                              path = path,
                              obs_data = obs_data,
                              pred_data = pred_data,
                              net = net,
                              nsamples = nsamples, # number of samples to use from the posterior in the stanfit object
                              addfunccol = addfunccol, # variable used for spatial weights
                              locID_pred = locID_pred,
                              chunk_size = chunk_size,
                              seed = seed) {

  if(!use_osm_ssn & missing(path)){
    stop("Using SSN, yet no path is supplied. Please provide a path for SSN.")
  }

  if(use_osm_ssn == TRUE & missing(osm_ssn)){
    stop("Using osm_ssn, yet no osm_ssn is supplied. Please provide an osm_ssn.")
  }


  stanfit <- object
  formula <- as.formula(attributes(stanfit)$formula)
  obs_resp <- obs_data[,gsub("\\~.*", "", formula)[2]]
  if( any( is.na(obs_resp) )) {stop("Can't have missing values in the response in the observed data. You need to impute them before")}

  if(!use_osm_ssn){
    out <- pred_ssnbayes(object = object,
                         use_osm_ssn = FALSE,
                         path = path,
                         obs_data = obs_data,
                         pred_data = pred_data,
                         net = net,
                         nsamples = nsamples, # number of samples to use from the posterior in the stanfit object
                         addfunccol = addfunccol, # variable used for spatial weights
                         locID_pred = locID_pred,
                         chunk_size = chunk_size,
                         seed = seed)
  }
  else{
    out <- pred_ssnbayes(object = object,
                         use_osm_ssn = TRUE,
                         osm_ssn = osm_ssn,
                         # family = family,
                         obs_data = obs_data,
                         pred_data = pred_data,
                         nsamples = nsamples, # number of samples to use from the posterior in the stanfit object
                         locID_pred = locID_pred,
                         chunk_size = chunk_size,
                         seed = seed)
  }
  out
}






#' Internal function used to perform spatio-temporal prediction in R using a stanfit object from ssnbayes()
#'
#' Use predict.ssnbayes() instead.
#' It will take an observed and a prediction data frame.
#' It requires the same number of observation/locations per day.
#' It requires location id (locID) and points id (pid).
#' The locID are unique for each site.
#' The pid is unique for each observation.
#' Missing values are allowed in the response but not in the covariates.
#'
#' @param object A stanfit object returned from ssnbayes
#' @param mat_all_preds A list with the distance/weights matrices
#' @param nsamples The number of samples to draw from the posterior distributions. (nsamples <= iter)
#' @param start (optional) The starting location id
#' @param chunk_size (optional) the number of locID to make prediction from
#' @param obs_data The observed data frame
#' @param pred_data The predicted data frame
#' @param net (optional) Network from the SSN object
#' @param seed (optional) A seed for reproducibility
#' @return A data frame
#' @export
#' @importFrom dplyr mutate %>% distinct left_join case_when
#' @importFrom plyr .
#' @importFrom rstan stan
#' @importFrom stats dist
#' @importFrom stats as.formula
#' @author Edgar Santos-Fernandez

krig <- function(object = object,
                 mat_all_preds = mat_all_preds,
                 nsamples = 10,
                 start = 1,
                 chunk_size = 50,
                 obs_data = obs_data,
                 pred_data = pred_data,
                 net = net,
                 seed = seed){

  if(missing(seed)) seed <- sample(1:1E6,1,replace=TRUE)
  set.seed(seed)

  stanfit <- object
  formula <- as.formula(attributes(stanfit)$formula)

  phi <- rstan::extract(stanfit, pars = 'phi')$phi

  samples <- sample(1:nrow(phi), nsamples, replace = FALSE)
  phi <- phi[samples]

  betas <- rstan::extract(stanfit, pars = 'beta')$beta
  betas <- betas[samples,]

  var_nug <- rstan::extract(stanfit, pars = 'var_nug')$var_nug
  var_nug <- var_nug[samples]

  var_td <- rstan::extract(stanfit, pars = 'var_td')$var_td #NB
  var_td <- var_td[samples]

  alpha_td <- rstan::extract(stanfit, pars = 'alpha_td')$alpha_td #NB
  alpha_td <- alpha_td[samples]

  locID_obs <- sort(unique(obs_data$locID))

  pred_data$temp <- NA

  locID_pred <- sort(unique(pred_data$locID))

  pred_data$locID0 <- as.numeric(factor(pred_data$locID)) # conseq locID
  pred_data$locID0 <- pred_data$locID0 + length(locID_obs) # adding numb of locID in obs dataset

  locID_pred0 <- sort(unique(pred_data$locID0))

  # obs data frame. no missing in temp

  locID_pred_1 <- locID_pred0[start:(start + chunk_size - 1)] # NB

  pred_data_1 <- pred_data[pred_data$locID0 %in% locID_pred_1,]

  N <- nrow(pred_data_1)

  options(na.action='na.pass') # to preserve the NAs
  out_list <- mylm(formula = formula, data = obs_data) # produces the design matrix

  response_obs <- out_list$y # response variable
  design_matrix_obs <- out_list$X # design matrix

  out_list_pred <- mylm(formula = formula, data = pred_data_1)
  X_pred <- as.matrix(out_list_pred$X)


  locID_pred2 <- unique(pred_data_1$pid)
  locID_obs2 <- sort(unique(obs_data$pid))

  locIDs <- unique(obs_data$locID)
  NlocIDs <- length( unique(obs_data$locID))

  n_obs <- NlocIDs
  n_pred <- unique(pred_data_1$locID)

  X_obs2 <- design_matrix_obs
  X_pred2 <- X_pred

  Y2 <- response_obs


  mat_all_preds_1 <- lapply(mat_all_preds,
                            function(x){x[c(locID_obs, locID_pred_1),
                                          c(locID_obs, locID_pred_1)  ] })


  mat_all_preds_1 <- lapply(mat_all_preds_1, function(x){colnames(x) = (1:ncol(x));
  rownames(x) = (1:nrow(x)); x  })


  dim(mat_all_preds_1$H)

  h <- mat_all_preds_1$H # hydro distance matrix
  D <- mat_all_preds_1$D
  fc <- mat_all_preds_1$flow.con.mat

  total_numb_points <- length(c(locID_obs, locID_pred_1))

  niter <- nsamples

  # time

  t <- length(unique(obs_data$date))

  Ypred <- matrix(NA, (nrow(pred_data_1)/t)*t, niter)

  for(k in 1:niter){

    C_t <- (phi[k] ^ abs(outer(1:t,1:t,"-")))/(1-( phi[k]^2) )
    inv_t <- solve(C_t)

    C_td <- matrix(NA, total_numb_points, total_numb_points)
    for(i in 1:total_numb_points){
      for(j in 1:total_numb_points){
        C_td[i,j] <- ifelse(fc[i,j] == 1,
                            var_td[k] * exp(- 3 * h[i,j] / alpha_td[k]),
                            var_td[k] * exp(- 3 * (D[i,j] + D[j,i]) / alpha_td[k])
        )
      }
    }
    C_td <- C_td + var_nug[k] * diag(total_numb_points)

    Coo_td <- C_td[1:n_obs,1:n_obs]
    Cop_td <- C_td[1:n_obs,(n_obs+1):(n_obs+chunk_size)]


    # Separable covariance space-time matrix

    Coo_all <- kronecker(C_t, Coo_td, FUN = "*")
    inv_all  <- chol2inv(chol(Coo_all))
    Cop_all <- kronecker(C_t, Cop_td, FUN = "*")

    mu = X_obs2 %*% betas[k,]
    mu_pred = X_pred2 %*% betas[k,]
    Ypred[,k] = mu_pred + t(Cop_all) %*%
      inv_all %*% (Y2 - mu)

  }

  pred_data_1$ypred_all <- apply(Ypred, 1,mean)

  cbind(pred_data_1[,c('locID0','locID','date')], data.frame(Ypred))

}














#' Generates an Open Street Maps Spatial Stream Network.
#'
#' It will take a dataframe containing sensor locations, with a coulmn `locID` denoting location IDs (can be string or numeric).
#' Requires a root sensor location, in which all sensor locations are connected to, and flow towards.
#'
#' @param sensor_data A dataframe containing sensor locations and a column called `locID` denoting the location IDs.
#' @param lon_name The column name containing the longitude of the sensor.
#' @param lat_name The column name containing the latitude of the sensor.
#' @param root_loc A sensor location which all water flows towards.
#' @param plot_network A bool indicating whether to plot the network. Useful for determining potential errors. Will be saved in object.
#' @param gen_pred_sites A bool indicating whether to generate new prediction sites.
#' @param num_pred_sites Required if `gen_pred_sites=TRUE`. An integer indicating the number of prediction sites to generate.
#' @return An `osm_ssn` object, containing: the graph network as an `sfnetwork`; the distance adjacency matrices; if `plot_network=TRUE`, an interactive leaflet plot, if `gen_pred_sites=TRUE`, an `sf` object containing the locations of the predictions sites.
#' @details Assumes longitude and latitude data are in the 4326 coordinate reference system. If they are not, data will need to be converted beforehand using `st_transform()`. Predictions on sites using gen_pred_sites is under development specially for tail-down models.
#' @export
#' @importFrom dplyr %>% filter left_join
#' @importFrom sfnetworks activate
#' @importFrom sf st_as_sf st_coordinates st_as_sfc
#' @importFrom leaflet leaflet addPolylines addTiles addMarkers addCircleMarkers addAwesomeMarkers makeAwesomeIcon
#' @author Sean Francis
#' @examples
#' \donttest{
#'# require('SSNdata')
#'# clear <- readRDS(system.file("extdata/clear_obs.RDS", package = "SSNdata"))
#'#
#'# # Generate osm ssn
#'# clear_osm_ssn <- generate_osm_ssn(clear,
#'#                                   "long", "lat",
#'#                                   root_loc = 12,
#'#                                   plot_network = TRUE,
#'#                                   gen_pred_sites = TRUE,
#'#                                   num_pred_sites = 20)
#'# # Show plot
#'# clear_osm_ssn$plot
#'}
generate_osm_ssn <- function(sensor_data,
                             lon_name, lat_name,
                             root_loc,
                             plot_network = TRUE,
                             gen_pred_sites = FALSE,
                             num_pred_sites = NA
){

  type <- NULL
  locID <- NULL
  weight <- NULL


  if('locID' %in% names(sensor_data) == FALSE){
    stop("There is no column locID on the data. Please, set a column called locID with the sensor locations")
  }

  if(gen_pred_sites & is.na(num_pred_sites)){
    stop("Numbr of prediction sites not supplied. Please supply a number.")
  }

  sensors_formatted <- format_sensor_data(sensor_data,
                                          lon_name,
                                          lat_name)

  river_line_network <- generate_river_network(sensors_formatted,
                                               root_loc)


  river_graph_network <- convert_network_to_graph(river_line_network,
                                                  sensors_formatted)

  net <- river_graph_network %>%
    activate("nodes")

  network_df <- st_as_sf(net, "nodes") %>%
    cbind(st_coordinates(st_as_sf(net, "nodes"))) %>%
    filter(type == "sensor")

  missing_locations <- setdiff(sensors_formatted$locID, network_df$locID)

  if(length(missing_locations) > 0){
    warning(
      paste0("The following location ID could not be blended in to the OSM SSN: ", missing_locations,
             " and will be excluded. Please plot network for clarification.\n")
    )
  }

  # Filter out missing sensors
  sensors_formatted_filtered <- sensors_formatted %>%
    filter(!locID %in% missing_locations)

  if(nrow(sensors_formatted_filtered) < nrow(sensors_formatted)){
    # Recalculated graph network with missing sensors
    river_graph_network <- convert_network_to_graph(river_line_network,
                                                    sensors_formatted_filtered)
  }

  dist_mat_all <- gen_distance_matrices(river_graph_network,
                                        sensors_formatted_filtered,
                                        root_loc)


  osm_ssn <- list(
    graph_network = river_graph_network,
    dist_mat_all = dist_mat_all
  )


  if(gen_pred_sites){

    predictions <- generate_prediction_locations(river_line_network,
                                                 num_pred_sites)

    all_sensors <- rbind(sensors_formatted_filtered,
                         predictions)

    river_graph_network_preds <- convert_network_to_graph(river_line_network,
                                                          all_sensors)

    dist_mat_sensors_preds <- gen_distance_matrices(river_graph_network_preds,
                                                    all_sensors,
                                                    root_loc = root_loc)

    net <- river_graph_network_preds %>%
      activate("nodes")

    network_df <- st_as_sf(net, "nodes") %>%
      cbind(st_coordinates(st_as_sf(net, "nodes"))) %>%
      filter(type == "sensor")

    pred_sites <- cbind(predictions, st_coordinates(predictions))

    names(pred_sites)[names(pred_sites) == "X"] <- lon_name
    names(pred_sites)[names(pred_sites) == "Y"] <- lat_name

    osm_ssn$dist_mat_all_preds <- dist_mat_sensors_preds
    osm_ssn$pred_sites <- pred_sites
  }

  if(plot_network){

    orig_loc <- cbind(sensors_formatted, st_coordinates(sensors_formatted))

    plot <- st_as_sf(net, "edges") %>%
      mutate(length = as.numeric(weight)) %>%
      st_as_sfc() %>%
      leaflet() %>%
      addPolylines() %>%
      addTiles() %>%
      addCircleMarkers(data = orig_loc,
                       color = "red",
                       label = ~paste0("Original location: ", locID))

    if(gen_pred_sites){

      network_df_pred <- network_df %>%
        filter(!locID %in% sensors_formatted_filtered$locID)

      network_df_sensors <- network_df %>%
        filter(locID %in% sensors_formatted_filtered$locID)

      plot <- plot %>%
        addAwesomeMarkers(
          data = network_df_pred,
          label = ~locID,
          icon = makeAwesomeIcon(icon = "circle",
                                 # iconColor = "white",
                                 markerColor = "green",
                                 library = "fa")
        )  %>%
        addMarkers(
          data = network_df_sensors,
          label = ~ifelse(is.na(locID),
                          yes = "Missing from network",
                          no = locID)
        )

    }else{

      plot <- plot %>%
        addMarkers(
          data = network_df,
          label = ~ifelse(is.na(locID),
                          yes = "Missing from network",
                          no = locID)
        )
    }

    osm_ssn$plot <- plot
  }


  class(osm_ssn) <- "osm_ssn"

  return(osm_ssn)
}




#' Fits a mixed linear regression model using Stan. This is an updated version of ssnbayes()
#'
#' It requires the same number of observation/locations per day.
#' It requires location id (locID) and points id (pid).
#' The locID are unique for each site.
#' The pid is unique for each observation.
#' Missing values are allowed in the response and in the covariates.
#' Missing values are imputed using the `mtsdi` pacakge, with a `glm` family.
#'
#' @param formula A formula as in lm()
#' @param family A description of the response distribution and link function to be used in the model. Must be one of: 'gaussian', 'binomial', 'poisson".
#' @param data A long data frame containing the locations, dates, covariates and the response variable. It has to have the locID and date. No missing values are allowed in the covariates.
#' The order in this data.fame MUST be: spatial locations (1 to S) at time t=1, then locations (1 to S) at t=2 and so on.
#' @param osm_ssn An `osm_ssn` generated using `generate_osm_ssn()`.
#' @param path File path with the name of the SpatialStreamNetwork object
#' @param space_method A list, first element must be one of 'use_ssn', 'use_osm_ssn', 'no_ssn. Whether to use SSN or osm ssn, and the spatial correlation structure. The second element is the spatial covariance structure. A 3rd element is a list with the lon and lat for Euclidean distance models.
#' @param time_method A list specifying the temporal structure (ar = Autorregressive; var = Vector autorregression) and column in the data with the time variable.
#' @param iter Number of iterations
#' @param warmup Warm up samples
#' @param chains Number of chains
#' @param cores Number of cores
#' @param refresh Sampler refreshing rate
#' @param net The network id (optional). Used when the SSN object cotains multiple networks.
#' @param addfunccol Variable to compute the additive function. Used to compute the spatial weights.
#' @param loglik Logic parameter denoting if the loglik will be computed by the model.
#' @param seed (optional) A seed for reproducibility
#' @return A list with the model fit
#' @details Missing values on the covariates and response can be imputed by default using mtsdi::mnimput(). We strongly recommend the documentation for this function be read before use.
#' The pid in the data has to be consecutive from 1 to the number of observations.
#' Users can use the SpatialStreamNetwork created with the SSN package. This will provide the spatial stream information used to compute covariance matrices. If that is the case, the data has
#' to have point ids (pid) matching the ones in SSN distance matrices, so that a mapping can occur.
#' @return It returns a ssnbayes object (similar to stan returns). It includes the formula used to fit the model. The output can be transformed into the stanfit class using class(fits) <- c("stanfit").
#' @export
#' @importFrom dplyr mutate %>% distinct left_join case_when filter
#' @importFrom plyr .
#' @importFrom rstan stan_model sampling
#' @importFrom stats dist
#' @importFrom mtsdi mnimput
#' @author Edgar Santos-Fernandez
#' @examples
#'\dontrun{
#'#options(mc.cores = parallel::detectCores())
#'#require(SSNdata)
#'#formula = temp ~ SLOPE + elev + h2o_area + air_temp + sin + cos
#'## Imports a data.frame containing observations and covariates
#'#clear <- readRDS(system.file("extdata/clear_obs.RDS", package = "SSNdata"))
#'#family <-  "gaussian"
#'
#'# # If using osm_ssn:
#'# # Generate osm_ssn
#'# clear_osm_ssn <- generate_osm_ssn(clear, "long", "lat", root_loc = 12, plot_network = TRUE)
#'# fit_ar <- ssnbayes2(formula,
#'#                     data = clear,
#'#                     osm_ssn = clear_osm_ssn,
#'#                     family = family,
#'#                     time_method = list("ar", "date"),
#'#                     space_method = list('use_osm_ssn', c("Exponential.taildown")),
#'#                     iter = 2000,
#'#                     warmup = 1000,
#'#                     chains = 3,
#'#                     cores = 3)
#'#
#'#
#'#
#'# If not using osm_ssn, and instead using SSN:
#'# Import SpatialStreamNetwork object
#'# path <- system.file("extdata/clearwater.ssn", package = "SSNbayes")
#'# fit_ar <- ssnbayes2(formula = formula,
#'#                     data = clear,
#'#                     path = path,
#'#                     family = family,
#'#                     time_method = list("ar", "date"),
#'#                     space_method = list('use_ssn', c("Exponential.taildown")),
#'#                     iter = 2000,
#'#                     warmup = 1000,
#'#                     chains = 3,
#'#                     cores = 3,
#'#                     net = 2, # second network on the ssn object
#'#                     addfunccol='afvArea')
#'#
#'#
#'#
#'# #space_method options examples
#'# #use list('no_ssn', 'Exponential.Euclid', c('lon', 'lat')) if no ssn object is available
#'}

ssnbayes2 <- function(formula = formula,
                      family = family,
                      data = data,
                      osm_ssn = osm_ssn,
                      path = NA,
                      time_method = time_method, # list("ar", "date")
                      space_method = space_method, #list('use_osm_ssn', 'Exponential.tailup'),
                      iter = 3000,
                      warmup = 1500,
                      chains = 3,
                      cores = 3,
                      refresh = max(iter/100, 1),
                      net = NA,
                      addfunccol = NA,
                      loglik = FALSE,
                      seed = seed){

  locID <- NULL

  ssn_type <- space_method[[1]]

  # checks
  if(! ssn_type %in% c("use_ssn", "use_osm_ssn", "no_ssn")){
    stop("First argument of space_method must be one of, 'use_ssn', 'use_osm_ssn', 'no_ssn'")
  }


  if(ssn_type == "use_ssn" & missing(path)){
    stop("Using SSN, yet no path is supplied. Please provide a path for SSN.")
  }

  if(ssn_type == "use_osm_ssn" & missing(osm_ssn)){
    stop("Using osm_ssn, yet no osm_ssn is supplied. Please provide an osm_ssn.")
  }



  if(missing(time_method)){
    stop("Need to define the method (ar or var) and the column associated with time")
  }

  if(length(time_method) == 1){
    stop("Need to specify the column in the the data with the time variable")
  }

  if(missing(family)){
    stop("Need to specify the response distribution with the family variable")
  }

  if(!family %in% c("gaussian", "binomial", "poisson")){
    stop("family must be one of: 'gaussian', 'binomial', 'poisson")
  }

  if("data.table" %in% class(data)){
    message("Data is of class data.table; converting to data.frame...")
    data <- as.data.frame(data)
  }



  time_points <- time_method[[2]]

  #if('date' %in% names(data) == FALSE) stop("There is no column date on the data. Please, set a column called date with the time")
  if('locID' %in% names(data) == FALSE) stop("There is no column locID on the data. Please, set a column called locID with the observation locations")

  if(missing(seed)) seed <- sample(1:1E6,1,replace=TRUE)

  if(!missing(space_method)){

    ssn_type <- space_method[[1]]

    if(ssn_type != 'no_ssn'){

      if(ssn_type == 'use_ssn'){ message('Using SSN object...') }
      else{ message('Using osm ssn object...') }

      ssn_object <- TRUE


      if(length(space_method) > 1){
        if(space_method[[2]] %in% c("Exponential.tailup", "LinearSill.tailup" , "Spherical.tailup" ,
                                    "Exponential.taildown" ,"LinearSill.taildown" ,"Spherical.taildown",
                                    "Exponential.Euclid") == FALSE) {stop("Need to specify one or more of the following covariance matrices: Exponential.tailup, LinearSill.tailup , Spherical.tailup ,
		Exponential.taildown, LinearSill.taildown, Spherical.taildown or Exponential.Euclid")}
        CorModels <- space_method[[2]]
      }
      if(length(space_method) == 1){
        CorModels <- "Exponential.tailup"
        message('Using an Exponential.tailup model...')
      }

    }
    else if(ssn_type == 'no_ssn'){
      message('No SSN object defined')
      ssn_object <- FALSE
      if(space_method[[2]] %in% c("Exponential.Euclid") == FALSE) {stop("Need to specify Exponential.Euclid")}
      # when using Euclidean distance, need to specify the columns with lon and lat.
      if(length(space_method) < 3){ stop("Please, specify the columns in the data frame with the longitude and latitude (c('lon', 'lat'))") }

      data$lon <- data[,names(data) == space_method[[3]][1]]
      data$lat <- data[,names(data) == space_method[[3]][2]]
      CorModels <- space_method[[2]]
    }


  }

  if(missing(space_method)) {space_method <- 'no_ssn'; ssn_object <- FALSE; CorModels <- "Exponential.Euclid" }# if missing use Euclidean distance




  # Cov
  cor_tu <- case_when(CorModels == "Exponential.tailup" ~ 1,
                      CorModels == "LinearSill.tailup" ~ 2,
                      CorModels == "Spherical.tailup" ~ 3,
                      TRUE ~ 5)
  cor_tu <- sort(cor_tu)[1]

  cor_td <- case_when(CorModels == "Exponential.taildown" ~ 1,
                      CorModels == "LinearSill.taildown" ~ 2,
                      CorModels == "Spherical.taildown" ~ 3,
                      TRUE ~ 5)
  cor_td <- sort(cor_td)[1]


  cor_ed <- case_when(CorModels == "Exponential.Euclid" ~ 1,
                      TRUE ~ 5)
  #CorModels == "Spherical.Euclid" ~ 2, #NB: to be implemented
  #CorModels == "Gaussian.Euclid" ~ 3)
  cor_ed <- sort(cor_ed)[1]

  cor_re <- case_when(CorModels == "RE1" ~ 1,
                      TRUE ~ 5)
  cor_re <- sort(cor_re)[1]


  data_initial <- 'data {
      int<lower=1> N;
      int<lower=1> K;
      int<lower=1> T;
      matrix[N,K] X[T] ; // real X[N,K,T]; //

      //int<lower = 0> N_y_obs; // number observed values
      //int<lower = 0> N_y_mis; // number missing values

      //int<lower = 1> i_y_obs[N_y_obs] ;  //[N_y_obs,T]
      //int<lower = 1> i_y_mis[N_y_mis] ;  // N_y_mis,T]

      matrix[N, N] W ; // spatial weights
      matrix[N, N] h ; // total hydrological dist
      matrix[N, N] I ; // diag matrix

      matrix[N, N] D ; // downstream hydrological dist matrix
      matrix[N, N] flow_con_mat; // flow conected matrix

      matrix[N, N] e ; // Euclidean dist mat
      real<lower=1> alpha_max ;
      '


  if (family == "gaussian"){
    data_com <- paste0(data_initial,
                       '
                       vector[N*T] y_obs;
                       }'
    )
  }else if (family == "poisson") {
    data_com <- paste0(data_initial,
                       '
                       int<lower=0> y_obs[N*T];
                       }'
    )
  }else{
    data_com <- paste0(data_initial,
                       '
                       int<lower=0, upper=1> y_obs[N*T];
                       }
                       '
    )
  }


  tdata_com <- paste0(
    '
    transformed data {
      vector[N] y_vec[T];

      for(t in 1:T){
        y_vec[t] = to_vector(y_obs[((t - 1) * N + 1):(t * N)]);
      }
    }

    '
  )


  param_com <- '
    parameters {
      vector[K] beta;
      real<lower=0> sigma_nug;

      vector[N] eta;

      //vector[N] y_mis;//declaring the missing y
    '


  param_phi_ar <- '
    real <lower=-1, upper = 1> phi; // NB
  '
  param_phi_var <- '
    vector<lower=-1, upper = 1> [T] phi  ; // vector of autoregresion pars
  '



  param_tu <- '
    real<lower=0> sigma_tu;
    real<lower=0, upper=alpha_max> alpha_tu;
    '

  param_td <- '
    real<lower=0> sigma_td; // sd of tail-down
    real<lower=0, upper=alpha_max> alpha_td; // range of the tail-down model
    '

  param_ed <- '
    real<lower=0> sigma_ed;
    real<lower=0, upper=alpha_max> alpha_ed; // range of the Euclidean dist model
    '

  param_re <- '
    real<lower=0> sigma_RE1;
    '


  tparam_com <- '
    transformed parameters {

      vector[N] epsilon[T]; // error term
      vector[N] mu[T]; // mean
      real<lower=0> var_nug; // nugget

       matrix[N, N] C_tu; //tail-up cov
       matrix[N, N] C1; //tail-up cov
       matrix[N, N] Ind; //tail-up indicator

       matrix[N, N] C_td; //tail-down cov
       matrix[N, N] Ind2; //tail-down indicator
       matrix[2,1] iji;

       matrix[N, N] C_ed ;// Euclidean cov

       matrix[N, N] C_re ;// random effect cov
       matrix[N, N] RE1; // random effect 1

     '

  tparam_tu <- '
    // tail up exponential
    real<lower=0> var_tu; // parsil tail-down
    '

  tparam_td <- '
    real<lower=0> var_td; // parsil tail-down
    '

  tparam_ed <- '
    real<lower=0> var_ed; //  Euclidean dist var
    '

  tparam_re <- '
    real<lower=0> var_RE1; // Random effect 1
    '


  tparam_com2 <- '

    var_nug = sigma_nug ^ 2; // variance nugget
    mu[1] = X[1] * beta;
    epsilon[1] = y_vec[1] - mu[1];
    '

  tparam_com_ar <- '
    for (t in 2:T){
      mu[t] = X[t] * beta;
      epsilon[t] = y_vec[t] - mu[t];
      mu[t] = mu[t] + phi * epsilon[t-1];
    }
    '

  tparam_com_var <- '
    for (t in 2:T){
      mu[t] = X[t] * beta;
      epsilon[t] = y_vec[t] - mu[t];
      mu[t] = mu[t] + phi .* epsilon[t-1]; // element wise mult two vectors
    }
    '


  tparam_tu2_exp <- '
    // tail up exponential
    var_tu = sigma_tu ^ 2; // variance tail-up
    C1 = var_tu * exp(- 3 * h / alpha_tu); // tail up exponential model
    C_tu = C1 .* W; // Hadamard (element-wise) product
    '

  tparam_tu2_lin <- '
    //Tail-up linear-with-sill model
    var_tu = sigma_tu ^ 2; // variance tail-up
    for (i in 1:N) {
      for (j in 1:N) {
        Ind[i,j] = (h[i,j] / alpha_tu) <= 1 ? 1 : 0; // indicator
      }
    }
    C1 = var_tu * (1 - (h / alpha_tu)) .* Ind ; //Tail-up linear-with-sill model
    C_tu = C1 .* W; // Hadamard (element-wise) product
    '

  tparam_tu2_sph <- '
    // Tail-up spherical model
    var_tu = sigma_tu ^ 2; // variance tail-up
    for (i in 1:N) {// Tail-up spherical model
      for (j in 1:N) {
        Ind[i,j] = (h[i,j] / alpha_tu) <= 1 ? 1 : 0; // indicator
      }
    }
    C1 = var_tu * (1 - (1.5 * h / alpha_tu) + (h .* h .* h / (2 * alpha_tu ^ 3))) .* Ind ; // Tail-up spherical model
    C_tu = C1 .* W; // Hadamard (element-wise) product
  '
  #tail-up models end


  #tail-down models start

  # tail-down exponential
  tparam_td2_exp <- '
    var_td= sigma_td ^ 2; // variance tail-down
     	for (i in 1:N) {// Tail-down exponential model
          for (j in 1:N) {
            if(flow_con_mat[i,j] == 1){ // if points are flow connected
               C_td[i,j] = var_td * exp(- 3 * h[i,j] / alpha_td);
            }
            else{// if points are flow unconnected
              C_td[i,j] = var_td * exp(- 3 * (D[i,j] + D[j,i]) / alpha_td);
            }
          }
        }

'

  #Tail-down linear-with-sill model
  tparam_td2_lin <- '
    var_td= sigma_td ^ 2; // variance tail-down
        for (i in 1:N) {// Tail-down linear-with-sill model
          for (j in 1:N) {
            if(flow_con_mat[i,j] == 1){ // if points are flow connected
              Ind2[i,j] = (h[i,j] / alpha_td) <= 1 ? 1 : 0; // indicator
              C_td[i,j] = var_td * (1 - (h[i,j] / alpha_td)) .* Ind2[i,j] ; //Tail-up linear-with-sill model
            }
            else{// if points are flow unconnected
              iji[1,1] = D[i,j];
              iji [2,1] = D[j,i];
              Ind2[i,j] = (max(iji) / alpha_td) <= 1 ? 1 : 0; // indicator
              C_td[i,j] = var_td * (1 - (max(iji) / alpha_td)) * Ind2[i,j] ;
            }
          }
        }

  '

  #tail-down spherical model
  tparam_td2_sph <- '
    var_td= sigma_td ^ 2; // variance tail-down
        for (i in 1:N) {// tail-down spherical model
          for (j in 1:N) {
            if(flow_con_mat[i,j] == 1){ // if points are flow connected
              Ind2[i,j] = (h[i,j] / alpha_td) <= 1 ? 1 : 0; // indicator
              C_td[i,j] = var_td * (1 - (1.5 * h[i,j] / alpha_td) + ( (h[i,j] ^ 3) / (2 * alpha_td ^ 3))) * Ind2[i,j];
            }
            else{// if points are flow unconnected
              iji[1,1] = D[i,j];
              iji [2,1] = D[j,i];
              Ind2[i,j] = (max(iji) / alpha_td) <= 1 ? 1 : 0; // indicator
              C_td[i,j] = var_td * (1 - (1.5 * min(iji) / alpha_td) + ( max(iji)/(2 * alpha_td)   )) * (1 - (max(iji) / alpha_td) ) ^ 2 * Ind2[i,j];
            }
          }
        }

'

  #tail-down models end


  tparam_ed2 <- '
	  //Euclidean distance models start
    var_ed = sigma_ed ^ 2; // var Euclidean dist
      C_ed = var_ed * exp(- 3 * e / alpha_ed); // exponential model
    //Euclidean distance models end
   '

  tparam_re2 <- '
    // random effect
    var_RE1 = sigma_RE1 ^ 2;
      C_re = var_RE1 * RE1;

    '

  tparam_chol <- paste0(
    if(family == "gaussian"){
      '
      matrix[N, N] L_chol = cholesky_decompose(C_tu + C_td + C_re + C_ed);
      '
    }else{
      '
      matrix[N, N] L_chol = cholesky_decompose(C_tu + C_td + C_re + C_ed + var_nug * I);
      '
    },
    # if(family == "gaussian"){
    #   '
    #   matrix[N, N] L_chol = cholesky_decompose(C_tu + C_td + C_re + C_ed);
    #   '
    # }else{
    #   '
    #   matrix[N, N] L_chol = cholesky_decompose(C_td + var_nug * I);
    #   '
    # },

    '
    vector[N] gp_function = L_chol * eta;

    vector[N*T] mu2; // a vector of length N * T

    for (i in 1:T) {
      mu2[((i - 1) * N + 1):(i * N)] = mu[i] + gp_function;
    }

  '
  )

  if(family == "gaussian"){
    target <- 'normal(mu2[t], var_nug);'
  }else if (family == "poisson"){
    target <- 'poisson_log(mu2[t]);'
  }else{
    target <- 'bernoulli_logit(mu2[t]);'
  }

  model_com <- paste0(
    '
  model {
    for (t in 1:N*T){
      y_obs[t] ~ ', target, '
    }
    // Multiplier for non-centred GP parameterisation
    eta ~ std_normal();

    sigma_nug ~ uniform(0,50); // cauchy(0,1) prior nugget effect
    phi ~ uniform(-1, 1); // or can use phi ~ normal(0.5,0.3); //NB informative
  '
  )

  model_tu <- '
    sigma_tu ~ uniform(0,100);  // or cauchy(0,2) prior sd  tail-up model
    alpha_tu ~ uniform(0, alpha_max);
    '

  model_td <- '
    sigma_td ~ uniform(0,100); // sd tail-down
    alpha_td ~ uniform(0, alpha_max);
    '

  model_ed <- '
    sigma_ed ~ uniform(0,100); // sd Euclidean dist
    alpha_ed ~ uniform(0, alpha_max); // Euclidean dist range
    '

  model_re <- '
    sigma_RE1 ~ uniform(0,5);
    '

  if(loglik == TRUE){
    loglik_code <- paste0(
      '
    vector[N*T] log_lik;

    for (t in 1:N*T){
      log_lik[t] = ', target, '
    }'
    )
  }else{
    loglik_code = ''
  }

  if(family == "gaussian"){
    target_fitted_vals <- 'normal_rng(mu2[t], var_nug);'
  }else if (family == "poisson"){
    target_fitted_vals <- 'poisson_log_rng(mu2[t]);'
  }else{
    target_fitted_vals <- 'bernoulli_logit_rng(mu2[t]);'
  }

  gen_quant <- paste0(
    '
    generated quantities {',
    loglik_code,'
      vector[N*T] fitted_values;

      for (t in 1:N*T) {
        fitted_values[t] = ', target_fitted_vals,'
      }
     }
    '
  )

  ssn_ar <- paste(
    data_com,

    tdata_com,

    param_com,

    if(cor_tu %in% 1:3) param_tu,
    if(cor_td %in% 1:3) param_td,
    if(cor_ed %in% 1:3) param_ed,
    if(cor_re %in% 1:3) param_re,

    if(time_method[[1]] == 'ar') param_phi_ar,
    if(time_method[[1]] == 'var') param_phi_var,

    '}',

    tparam_com,
    if(cor_tu %in% 1:3) tparam_tu,
    if(cor_td %in% 1:3) tparam_td,
    if(cor_ed %in% 1:3) tparam_ed,
    if(cor_re %in% 1:3) tparam_re,
    tparam_com2,


    if(time_method[[1]] == 'ar') tparam_com_ar,
    if(time_method[[1]] == 'var') tparam_com_var,



    case_when(cor_tu == 1 ~ tparam_tu2_exp,
              cor_tu == 2 ~ tparam_tu2_lin,
              cor_tu == 3 ~ tparam_tu2_sph,
              cor_tu >= 4 | cor_tu <= 0 ~ 'C_tu = rep_matrix(0, N, N);'),


    case_when(cor_td == 1 ~ tparam_td2_exp,
              cor_td == 2 ~ tparam_td2_lin,
              cor_td == 3 ~ tparam_td2_sph,
              cor_td >= 4 | cor_td <= 0 ~ 'C_td = rep_matrix(0, N, N);'),

    case_when(cor_ed == 1 ~ tparam_ed2,
              cor_ed >= 2 | cor_ed <= 0 ~ 'C_ed = rep_matrix(0, N, N);'),

    case_when(cor_re == 1 ~ tparam_re2,
              cor_re >= 2 | cor_re <= 0 ~ 'C_re = rep_matrix(0, N, N);'),

    tparam_chol,

    '}',
    model_com,
    if(cor_tu %in% 1:3)model_tu,
    if(cor_td %in% 1:3)model_td,
    if(cor_ed %in% 1:3)model_ed,
    if(cor_re %in% 1:3)model_re,
    '}',

    gen_quant
  )

  `%notin%` <- Negate(`%in%`)

  pars <- c(
    case_when(cor_tu %in% 1:3 ~ c('var_tu', 'alpha_tu'),
              cor_tu %notin% 1:3 ~ ""),

    case_when(cor_td %in% 1:3 ~ c('var_td', 'alpha_td'),
              cor_td %notin% 1:3 ~ ""),

    case_when(cor_ed %in% 1:3 ~ c('var_ed', 'alpha_ed'),
              cor_ed %notin% 1:3 ~ ""),

    case_when(cor_re %in% 1:3 ~ c('var_re', 'alpha_re'),
              cor_re %notin% 1:3 ~ ""),

    if(loglik == TRUE) 'log_lik',

    'var_nug',
    'beta',
    'phi',
    'fitted_values',
    'eta'
  )

  pars <- pars[pars != '']


  # data part
  old <- options()        # old options
  on.exit(options(old)) 	# reset once exit the function

  options(na.action='na.pass') # to preserve the NAs

  # If using an osm_ssn, filter the data to the distance matrices column names
  #  - in the event some sensor locations were not able to be snapped to the osm ssn
  if(space_method[[1]] == "use_osm_ssn"){
    sensor_names <- osm_ssn$dist_mat_all$e %>% colnames()

    data <- data %>%
      filter(locID %in% sensor_names)
  }

  out_list <- mylm(formula = formula, data = data) # produces the design matrix

  if (anyNA(out_list$y)){

    if(interactive()){
      toImpute = readline(prompt = "There is data missing in your response variable. Do you want to impute? [Y/n]") %>%
        tolower()
      warning("Please read function documentation for imputation assumptions\n")

      while(toImpute %notin% c('y', 'n')){
        toImpute = readline(prompt = "Please enter either 'Y' or 'n'") %>%
          tolower()
        if(toImpute == ''){message("Function exited");invokeRestart("abort")}
      }

      if(toImpute == 'n'){
        stop("Please supply data with non-missing response values")
      }
    }

    # # Impute missing data
    # require(mtsdi)

    message("Imputing missing data")

    data_impute <- mtsdi::mnimput(
      formula = formula,
      dataset = data,
      eps = 1e-3,
      ts = FALSE,
      method = "glm"
    )$filled.dataset

    # Extract the response variable name from the formula
    response_var <- formula[[2]]

    # Ensure for pois, data is still discrete
    if(family == "poisson"){
      data_impute[[response_var]] <- round(data_impute[[response_var]], 0)
    }

    if(family == "binomial"){
      data_impute[[response_var]] <- ifelse(data_impute[[response_var]] > 0.5,
                                            yes = 1,
                                            no = 0)
    }

    # Redefine the out_list object
    out_list <- mylm(formula = formula, data = data_impute) # produces the design matrix
  }


  response <- out_list$y # response variable



  design_matrix <- out_list$X # design matrix

  obs_data <- data

  ndays <- length(unique(obs_data[, names(obs_data) %in% time_points] ))


  N <- nrow(obs_data)/ndays #nobs


  nobs <- nrow(obs_data)/ndays #nobs

  obs_data$date_num <- as.numeric(factor(obs_data[, names(obs_data) %in% time_points]	))

  resp_var_name <- gsub("[^[:alnum:]]", " ", formula[2])
  obs_data$y <- obs_data[,names(obs_data) %in% resp_var_name]


  # array structure
  X <- design_matrix #cbind(1,obs_data[, c("X1", "X2", "X3")]) # design matrix

  # NB: this array order is Stan specific
  Xarray <- aperm(array( c(X), dim=c(N, ndays, ncol(X)) ),c(2, 1, 3))

  y_obs <- response #[!is.na(response)]
  #
  #
  # Y <- matrix(nrow = N, ncol = T)
  #
  # # Loop through each t from 1 to T
  # for (t in 1:T) {
  #   # Extract the segment of y corresponding to Y[t] and assign it to the t-th column of Y
  #   Y[, t] <- y_obs[((t - 1) * N + 1):(t * N)]
  # }

  # # index for observed values
  # i_y_obs <- obs_data[!is.na(obs_data$y),]$pid
  #
  # # index for missing values
  # i_y_mis <- obs_data[is.na(obs_data$y),]$pid

  if(space_method[[1]] == "use_osm_ssn"){
    mat_all <- osm_ssn$dist_mat_all
  }
  else{
    message("Using supplied SSN")
    if(ssn_object == TRUE){ # the ssn object exist?
      mat_all <- dist_weight_mat(path = path, net = net, addfunccol = addfunccol)
    }

    if(ssn_object == FALSE){ # the ssn object does not exist- purely spatial

      first_date <- unique(obs_data[, names(obs_data) %in% time_points])[1]

      di <- dist(obs_data[obs_data$date == first_date, c('lon', 'lat')], #data$date == 1
                 method = "euclidean",
                 diag = FALSE,
                 upper = FALSE) %>% as.matrix()
      mat_all <-  list(e = di, D = di, H = di, w.matrix = di, flow.con.mat = di)
    }
  }


  data_list <- list(N = N, # obs + preds  points
                    T = ndays, # time points
                    K = ncol(X),  # ncol of design matrix
                    y_obs = y_obs,# y values in the obs df

                    # N_y_obs = length(i_y_obs),  #nrow(i_y_obs) numb obs points
                    # N_y_mis = length(i_y_mis), #nrow(i_y_mis) numb preds points
                    #
                    # i_y_obs = i_y_obs, # index of obs points
                    # i_y_mis = i_y_mis, # index of preds points

                    X = Xarray, # design matrix
                    mat_all = mat_all,
                    # Below is original code
                    # alpha_max = 4 * max(mat_all$H)
                    alpha_max = 2 * max(mat_all$H) # THIS IS CHANGED FROM ORIGINAL VERSION
  ) # a list with all the distance/weights matrices



  data_list$e = data_list$mat_all$e #Euclidean dist
  #for tail-up
  data_list$h = data_list$mat_all$H # total stream distance
  data_list$W = data_list$mat_all$w.matrix # spatial weights


  #for tail-down
  data_list$flow_con_mat = data_list$mat_all$flow.con.mat #flow connected matrix
  data_list$D = data_list$mat_all$D #downstream hydro distance matrix

  #RE1 = RE1mm # random effect matrix

  data_list$I = diag(1, nrow(data_list$e), nrow(data_list$e))  # diagonal matrix

  if(family %in% c("poisson", "binomial")){
    ini = 0
  }else{
    ini <- function(){list(var_nug =  .1)}
  }



  message("Compiling model ...")

  compiled_model <- rstan::stan_model(
    model_code = ssn_ar,
    model_name = "ssn_ar"
  )

  message("Done!")

  message("Sampling model ...")

  fit <- rstan::sampling(object = compiled_model,
                         data = data_list,
                         pars = pars,
                         iter = iter,
                         warmup = warmup,
                         init = ini,
                         chains = chains,
                         cores = cores,
                         verbose = FALSE,
                         seed = seed,
                         refresh = refresh)

  message("Done!")

  attributes(fit)$formula <- formula

  class(fit) <- 'ssnbayes2'

  return(fit)
}





#' Internal function used to perform spatio-temporal prediction in R using a stanfit object from ssnbayes()
#'
#' Use predict.ssnbayes() instead.
#' It will take an observed and a prediction data frame.
#' It requires the same number of observation/locations per day.
#' It requires location id (locID) and points id (pid).
#' The locID are unique for each site.
#' The pid is unique for each observation.
#' Missing values are allowed in the response but not in the covariates.
#'
#' @param object A stanfit object returned from ssnbayes
#' @param mat_all_preds A list with the distance/weights matrices
#' @param nsamples The number of samples to draw from the posterior distributions. (nsamples <= iter)
#' @param start (optional) The starting location id
#' @param chunk_size (optional) the number of locID to make prediction from
#' @param obs_data The observed data frame
#' @param pred_data The predicted data frame
#' @param net (optional) Network from the SSN object
#' @param seed (optional) A seed for reproducibility
# #' @param family A description of the response distribution and link function to be used in the model. Must be one of: 'gaussian', 'binomial', 'poisson".
#' @return A data frame
#' @export
#' @importFrom dplyr mutate %>% distinct left_join case_when
#' @importFrom plyr .
#' @importFrom rstan stan
#' @importFrom stats dist
#' @importFrom stats as.formula
#' @author Edgar Santos-Fernandez

krig2 <- function(object = object,
                  mat_all_preds = mat_all_preds,
                  nsamples = 10,
                  start = 1,
                  chunk_size = 50,
                  obs_data = obs_data,
                  pred_data = pred_data,
                  net = net,
                  seed = seed
                  # ,
                  # family = family
){

  if(missing(seed)) seed <- sample(1:1E6,1,replace=TRUE)
  set.seed(seed)

  stanfit <- object
  formula <- as.formula(attributes(stanfit)$formula)

  phi <- rstan::extract(stanfit, pars = 'phi')$phi

  samples <- sample(1:nrow(phi), nsamples, replace = FALSE)
  # phi <- phi[samples]
  phi <- rbeta(samples, 9,1)

  betas <- rstan::extract(stanfit, pars = 'beta')$beta
  betas <- betas[samples,]

  var_nug <- rstan::extract(stanfit, pars = 'var_nug')$var_nug
  var_nug <- var_nug[samples]

  var_td <- rstan::extract(stanfit, pars = 'var_td')$var_td #NB
  var_td <- var_td[samples]

  alpha_td <- rstan::extract(stanfit, pars = 'alpha_td')$alpha_td #NB
  alpha_td <- alpha_td[samples]


  locID_obs <- sort(unique(obs_data$locID))

  if(is.numeric(locID_obs)){
    if(length(setdiff(locID_obs, 1:length(locID_obs))) > 1){
      warning("Numeric locID is not consecutive.")
      message("Converting locID to consecutive.")

      locID_obs <- as.numeric(factor(locID_obs)) # conseq locID
    }
  }


  pred_data$temp <- NA

  locID_pred <- sort(unique(pred_data$locID))

  pred_data$locID0 <- as.numeric(factor(pred_data$locID)) # conseq locID
  pred_data$locID0 <- pred_data$locID0 + length(locID_obs) # adding numb of locID in obs dataset

  locID_pred0 <- sort(unique(pred_data$locID0))

  # obs data frame. no missing in temp

  locID_pred_1 <- locID_pred0[start:(start + chunk_size - 1)] # NB

  pred_data_1 <- pred_data[pred_data$locID0 %in% locID_pred_1,]

  N <- nrow(pred_data_1)

  options(na.action='na.pass') # to preserve the NAs
  out_list <- mylm(formula = formula, data = obs_data) # produces the design matrix

  response_obs <- out_list$y # response variable
  design_matrix_obs <- out_list$X # design matrix

  out_list_pred <- mylm(formula = formula, data = pred_data_1)
  X_pred <- as.matrix(out_list_pred$X)


  locID_pred2 <- unique(pred_data_1$pid)
  locID_obs2 <- sort(unique(obs_data$pid))

  locIDs <- unique(obs_data$locID)
  NlocIDs <- length( unique(obs_data$locID))

  n_obs <- NlocIDs
  n_pred <- unique(pred_data_1$locID)

  X_obs2 <- design_matrix_obs
  X_pred2 <- X_pred

  Y2 <- response_obs


  mat_all_preds_1 <- lapply(mat_all_preds,
                            function(x){x[c(locID_obs, locID_pred_1),
                                          c(locID_obs, locID_pred_1)  ] })


  mat_all_preds_1 <- lapply(mat_all_preds_1, function(x){colnames(x) = (1:ncol(x));
  rownames(x) = (1:nrow(x)); x  })


  dim(mat_all_preds_1$H)

  h <- mat_all_preds_1$H # hydro distance matrix
  D <- mat_all_preds_1$D
  fc <- mat_all_preds_1$flow.con.mat

  total_numb_points <- length(c(locID_obs, locID_pred_1))

  niter <- nsamples

  # time

  t <- length(unique(obs_data$date))

  Ypred <- matrix(NA, (nrow(pred_data_1)/t)*t, niter)

  for(k in 1:niter){

    C_t <- (phi[k] ^ abs(outer(1:t,1:t,"-")))/(1-( phi[k]^2) )
    inv_t <- solve(C_t)

    C_td <- matrix(NA, total_numb_points, total_numb_points)
    for(i in 1:total_numb_points){
      for(j in 1:total_numb_points){
        C_td[i,j] <- ifelse(fc[i,j] == 1,
                            var_td[k] * exp(- 3 * h[i,j] / alpha_td[k]),
                            var_td[k] * exp(- 3 * (D[i,j] + D[j,i]) / alpha_td[k])
        )
      }
    }
    C_td <- C_td + var_nug[k] * diag(total_numb_points)

    Coo_td <- C_td[1:n_obs,1:n_obs]
    Cop_td <- C_td[1:n_obs,(n_obs+1):(n_obs+chunk_size)]


    # Separable covariance space-time matrix

    Coo_all <- kronecker(C_t, Coo_td, FUN = "*")
    inv_all  <- chol2inv(chol(Coo_all))
    Cop_all <- kronecker(C_t, Cop_td, FUN = "*")

    mu = X_obs2 %*% betas[k,]
    mu_pred = X_pred2 %*% betas[k,]
    Ypred[,k] = mu_pred + t(Cop_all) %*%
      inv_all %*% (Y2 - mu)

  }

  pred_data_1$ypred_all <- apply(Ypred, 1,mean)

  cbind(pred_data_1[,c('locID0','locID','date')], data.frame(Ypred))

  # # if family .... rounding / convert to 0, 1
  # if(family == "poisson"){
  #
  # }
}








##### Internal Functions #####


#' Internal function used to perform spatio-temporal prediction in R using a stanfit object from ssnbayes()
#'
#' Use predict.ssnbayes() instead.
#' It will take an observed and a prediction data frame.
#' It requires the same number of observation/locations per day.
#' It requires location id (locID) and points id (pid).
#' The locID are unique for each site.
#' The pid is unique for each observation.
#' Missing values are allowed in the response but not in the covariates.
#'
#' @param object A stanfit object returned from ssnbayes
#' @param use_osm_ssn Use a supplied osm_ssn instead a SpatialStreamNetwork object.
#' @param osm_ssn The osm_ssn to be used.
# #' @param family A description of the response distribution and link function to be used in the model. Must be one of: 'gaussian', 'binomial', 'poisson".
#' @param path If not using an osm_ssn, path with the name of the SpatialStreamNetwork object.
#' @param obs_data The observed data frame
#' @param pred_data The predicted data frame
#' @param net (optional) Network from the SSN object
#' @param nsamples The number of samples to draw from the posterior distributions. (nsamples <= iter)
#' @param addfunccol If not using an osm_ssn, the variable used for spatial weights.
#' @param chunk_size (optional) the number of locID to make prediction from
#' @param locID_pred (optional) the location id for the predictions. Used when the number of pred locations is large.
#' @param seed (optional) A seed for reproducibility
#' @return A data frame
#' @export
#' @importFrom dplyr mutate %>% distinct left_join case_when
#' @importFrom plyr .
#' @importFrom rstan stan
#' @importFrom stats dist
#' @author Edgar Santos-Fernandez
#' @examples
#' #pred <- pred_ssnbayes(path = path,
#' #obs_data = clear,
#' #stanfit = fit_ar,
#' #pred_data = preds,
#' #net = 2,
#' #nsamples = 100, # number of samples to use from the posterior in the stanfit object
#' #addfunccol = 'afvArea') # variable used for spatial weights



pred_ssnbayes <- function(
    object = object,
    use_osm_ssn = TRUE,
    osm_ssn = osm_ssn,
    # family = family,
    path = path,
    obs_data = obs_data,
    pred_data = pred_data,
    net = 1,
    nsamples = 100, # number of samples to use from the posterior in the stanfit object
    addfunccol = 'afvArea', # variable used for spatial weights
    locID_pred = locID_pred, # location ID of the points to predict
    chunk_size = chunk_size,
    seed = seed
){

  stanfit <- object
  class(object) <- 'stanfit'

  if(missing(seed)) seed <- sample(1:1E6,1,replace=TRUE)
  set.seed(seed)

  if(missing(use_osm_ssn)) use_osm_ssn <- FALSE

  if(!use_osm_ssn){
    mat_all_preds <- dist_weight_mat_preds(path = path,
                                           net = net,
                                           addfunccol = addfunccol)
  }
  else{
    mat_all_preds <- osm_ssn$dist_mat_all
  }

  # the row and col names is not conseq
  rownames(mat_all_preds$e) <- 1:(nrow(mat_all_preds$e))
  colnames(mat_all_preds$e) <- 1:(nrow(mat_all_preds$e))

  rownames(mat_all_preds$D) <- 1:(nrow(mat_all_preds$D))
  colnames(mat_all_preds$D) <- 1:(nrow(mat_all_preds$D))

  rownames(mat_all_preds$H) <- 1:(nrow(mat_all_preds$H))
  colnames(mat_all_preds$H) <- 1:(nrow(mat_all_preds$H))

  rownames(mat_all_preds$w.matrix) <- 1:(nrow(mat_all_preds$w.matrix))
  colnames(mat_all_preds$w.matrix) <- 1:(nrow(mat_all_preds$w.matrix))

  rownames(mat_all_preds$flow.con.mat) <- 1:(nrow(mat_all_preds$flow.con.mat))
  colnames(mat_all_preds$flow.con.mat) <- 1:(nrow(mat_all_preds$flow.con.mat))

  obs_points <- length(unique(obs_data$locID))
  pred_points <- length(unique(pred_data$locID))

  if(missing(chunk_size)) chunk_size <- pred_points

  out_all <- NULL

  is <- ceiling(pred_points/chunk_size)

  for(j in 1:is){
    start <- ((j - 1) * chunk_size + 1)

    chunk_size <- ifelse(j != is, chunk_size, pred_points - (j - 1) * chunk_size)

    if(!use_osm_ssn){
      out <- krig(object = object,
                  mat_all_preds = mat_all_preds,
                  nsamples = nsamples,
                  start = start,
                  chunk_size = chunk_size,
                  obs_data = obs_data,
                  pred_data = pred_data,
                  net = net)
    }else{
      out <- krig2(object = object,
                   mat_all_preds = mat_all_preds,
                   nsamples = nsamples,
                   start = start,
                   chunk_size = chunk_size,
                   obs_data = obs_data,
                   pred_data = pred_data
                   # ,
                   # family = family
      )
    }

    out_all <- rbind(out_all, out)
  }

  data.frame(out_all)

}







#' Converts a dataframe containing sensor data in to an `sf` object
#'
#' @param sensor_data A dataframe containing sensor locations and a column called `locID` denoting the location IDs.
#' @param lon_name The column name containing the longitude of the sensor.
#' @param lat_name The column name containing the latitude of the sensor.
#' @return An `sf` dataframe containing sensor locations and `locID`.
#' @importFrom dplyr all_of %>% distinct select
#' @importFrom sf st_as_sf
#' @author Sean Francis
format_sensor_data <- function(sensor_data, lon_name, lat_name){
  sensor_data %>%
    select(all_of(c("locID", lat_name, lon_name))) %>%
    distinct %>%
    st_as_sf(coords = c(lon_name,lat_name),
             crs = 4326) %>%
    return()
}


#' Takes an `sf` object as input, and queries Open Street Maps for rivers and streams within the bounds of the sensor locations.
#'
#' @param sensor_locations An `sf` dataframe containing sensor locations and `locID`.
#' @param root_loc A sensor location which all water flows towards.
#' @return An `sf` object of `MULTILINESTRING`(s) containing the streams and rivers, cropped to the sensor location, with all lines not containing the outlet/root_loc removed.
#' @details This function queries Open Street Maps using `osmdata` package in R. Crops the data to the bounding box of the sensor locations. Removes all lines not connected to the outlet/root_loc.
#' @importFrom sf st_bbox st_touches st_crop st_nearest_feature
#' @importFrom osmdata opq add_osm_feature osmdata_sf
#' @importFrom dplyr %>% group_by summarise
#' @importFrom igraph graph_from_adj_list components
#' @importFrom purrr pluck
#' @author Sean Francis
generate_river_network <- function(sensor_locations, root_loc){

  locID <- NULL

  # Get the bounding box of the sensor locations
  sensors_bbox <- st_bbox(sensor_locations)

  buffer <- 0.05

  # Increase the size of the bounding box
  sensors_bbox[1] <- sensors_bbox[1] - buffer
  sensors_bbox[2] <- sensors_bbox[2] - buffer
  sensors_bbox[3] <- sensors_bbox[3] + buffer
  sensors_bbox[4] <- sensors_bbox[4] + buffer

  message("Searching OSM for river data...")
  # Query OSM
  lines <- ((opq(bbox = sensors_bbox,
                 timeout = 2000,
                 memsize = 10e8) %>%
               add_osm_feature(key = 'waterway', value = c('river', 'stream')) %>%
               osmdata_sf())$osm_lines)
  message("Done!")

  # Crop river network to the bounding box
  lines <- st_crop(lines, sensors_bbox)

  # Concatenate `LINESTRING` into connected `MULTILINESTRING`(s).
  concat <- (lines %>%
               st_touches() %>%
               graph_from_adj_list() %>%
               components())$membership

  river_network <- group_by(
    lines,
    section = as.character({{concat}})
  ) %>%
    summarise()

  root_loc_point <- sensor_locations %>%
    filter(locID == root_loc)


  touching_index <- st_nearest_feature(root_loc_point, river_network)

  # figure out which multilinestring contains root_loc, remove all others
  # assumes all sensors are connected by flow
  river_network_main <- river_network[touching_index, ]

  return(river_network_main)
}


#' Combines two `sf` objects, the sensor locations and river network, in to an `sfnetwork`.
#'
#' @param river_network An `sf` object of `MULTILINESTRING`(s) containing the streams and rivers, cropped to the sensor locations.
#' @param sensor_locations An `sf` dataframe containing sensor locations and `locID`.
#' @return An `sfnetwork` graph containing the river network and sensors.
#' @details This function "snaps" sensor locations to the nearest point on the network.
#' Additionally, the sensor locations are stored in the graph, `g`, and can be accessed via `V(g)$type`, as "sensor".
#' @importFrom sf st_zm st_geometry st_cast st_sfc st_crs
#' @importFrom dplyr %>% mutate
#' @importFrom sfnetworks as_sfnetwork st_network_blend
#' @importFrom igraph V
#' @author Sean Francis
convert_network_to_graph <- function(river_network,
                                     sensor_locations
                                     # # If elevation variable is not given
                                     #, elev_var = NA
){

  # if(is.na(elev_var)){
  #   sensor_locations <- elevatr::get_elev_point(sensor_locations,
  #                                               prj = st_crs(river_network),
  #                                               src = "epqs") %>%
  #     arrange(desc(elevation))
  # }

  geom <- NULL

  network <-
    river_network %>%
    st_zm(geom) %>%
    st_geometry() %>%
    st_cast("LINESTRING") %>%
    # back to sfc
    st_sfc(crs = st_crs(river_network)) %>%
    # directed graph
    as_sfnetwork(
      length_as_weight = TRUE, # Check what this argument does
      directed = FALSE) %>%
    # Add sensors as nodes - this does the snapping to lines automatically
    st_network_blend(sensor_locations)



  network <- network %>%
    mutate(
      "node_ID" = 1:length(network)
    )

  num_sensors <- nrow(sensor_locations)

  len_net <- length(igraph::V(network))

  num_nodes_not_needed <- len_net-num_sensors


  igraph::V(network)$type <- c(rep("node", times = num_nodes_not_needed),
                               rep("sensor", times = num_sensors))

  return(network)

}

# gen_river_network_df <- function(network, ){
#
#   network_df <- as.data.frame(network,
#                               what = "vertices") %>%
#     filter(type == "sensor") %>%
#     st_as_sf
#
#   network_df <- network_df %>%
#     cbind(st_coordinates(network_df)) %>%
#     arrange()
#
#   return(network_df)
# }


#' Creates a list containing the stream distances and weights
#'
#' @param network An `sfnetwork` containing the river network and sensor locations.
#' @param sensor_locations An `sf` dataframe containing sensor locations and `locID`.
#' @param root_loc A sensor location which all water flows towards.
#' @return A list of matrices.
#' @importFrom dplyr %>%
#' @importFrom igraph distances shortest_paths
#' @importFrom sf st_coordinates
#' @importFrom geosphere distHaversine
#' @importFrom stats rbeta
# #' @importMethodsFrom sfnetworks as.data.frame
#' @author Sean Francis
gen_distance_matrices <- function(network,
                                  sensor_locations,
                                  root_loc){

  network_df_all <- as.data.frame(network,
                                  what = "vertices")

  root <- which(network_df_all$locID == root_loc)

  num_sensors <- nrow(sensor_locations)

  sensor_indexes <- which(!is.na(network_df_all$locID))

  num_nodes <- nrow(network_df_all)

  names <- network_df_all$locID %>%
    unique %>%
    stats::na.omit()

  ordered_names <- sensor_locations$locID %>%
    as.character


  directed_dist <- matrix(0,
                          nrow = num_nodes,
                          ncol = num_nodes)

  w <- matrix(0,
              nrow = num_sensors,
              ncol = num_sensors)


  dimnames(w) <- list(names, names)

  w <- w[ordered_names,
         ordered_names]

  message("Tail-up model cannot be computed. Please ensure you use tail-down for model fitting.")

  hydro_distance <- distances(network,
                              v = sensor_indexes,
                              to = sensor_indexes,
                              mode = "out")

  dimnames(hydro_distance) <- list(names, names)

  hydro_distance <- hydro_distance[ordered_names,
                                   ordered_names]




  for (sensor in sensor_indexes){

    if(sensor == root){
      directed_dist[root, root] <- 0
    }else{

      # Get shortest path
      path <- (shortest_paths(network,
                              from = sensor,
                              to = root,
                              mode = 'out')$vpath
      )[[1]] %>% as.numeric


      path <- path[path %in% sensor_indexes]


      tmp_directed_distances <- distances(network,
                                          v = sensor,
                                          to = path,
                                          mode = "out")

      dimnames(tmp_directed_distances) <- list(sensor, path)

      # Get the above and put in matrix
      directed_dist[sensor, path] <- tmp_directed_distances

    }
  }

  directed_distances <- directed_dist[sensor_indexes, sensor_indexes]

  dimnames(directed_distances) <- list(names, names)

  # reorder
  directed_distances <- directed_distances[ordered_names,
                                           ordered_names]


  D <- matrix(0,num_sensors, num_sensors)

  dimnames(D) <- list(ordered_names,
                      ordered_names)

  e <- matrix(0,num_sensors, num_sensors)

  dimnames(e) <- list(ordered_names,
                      ordered_names)


  # Now we need to iterate through the matrix
  for (i in row.names(D)){
    for (j in colnames(D)){

      # Calculate e while we're iterating through the matrices
      tmp_e <- cbind(sensor_locations,
                     st_coordinates(sensor_locations))

      point_i <- c(tmp_e$X[tmp_e$locID == i],
                   tmp_e$Y[tmp_e$locID == i])

      point_j <- c(tmp_e$X[tmp_e$locID == j],
                   tmp_e$Y[tmp_e$locID == j])

      e[i, j] <- distHaversine(point_i, point_j)

      # Calculate D
      if(i != j){

        node_i <- which(network_df_all$locID == i)
        node_j <- which(network_df_all$locID == j)

        # Get the path from
        path_i <- (shortest_paths(network,
                                  from = node_i,
                                  to = root,
                                  mode = 'out')$vpath
        )[[1]] %>% as.numeric


        # Get the path to
        path_j <- (shortest_paths(network,
                                  from = node_j,
                                  to = root,
                                  mode = 'out')$vpath
        )[[1]] %>% as.numeric

        # Find the common node
        intersecting_nodes <- intersect(path_i, path_j)


        if (length(intersecting_nodes) == 0){
          D[i, j] <- 0
        }else{
          # Find the closest node to the "from" node
          node_dist <- distances(network,
                                 v = node_i,
                                 to = intersecting_nodes)

          D[i, j] <- node_dist[1, which(node_dist == min(node_dist))]
        }

      }

      # Set all instances of "Inf" in hydro_distance to 0
      if (hydro_distance[i, j] == Inf) {hydro_distance[i, j] = 0}
    }
  }

  # Set flow connectivity to directed distance
  flow.con.mat <- directed_distances

  for (i in 1:(num_sensors)){
    for (j in 1:(num_sensors)){

      flow.con.mat[i, j] = ifelse(flow.con.mat[i, j] == 0,
                                  yes = 0,
                                  no = 1)

      if(i == j){flow.con.mat[i, j] <- 1}

    }
  }

  # Convert to symmetric according to top triangular
  flow.con.mat = flow.con.mat + t(flow.con.mat) - 2*diag(diag(flow.con.mat)) + diag(1,nrow=dim(flow.con.mat)[1])


  # total distance
  # H <- D + base::t(D)

  H <- hydro_distance

  # e_df <- network_df_all %>%
  #   filter(type == "sensor") %>%
  #   st_as_sf()

  # e <- sensor_locations %>%
  #   cbind(st_coordinates(sensor_locations)) %>%
  #   as.data.frame() %>%
  #   select(X, Y) %>%
  #   dist(., method = "euclidean", diag = FALSE, upper = FALSE) %>%
  #   as.matrix()
  #
  # dimnames(e) <- dimnames(hydro_distance)

  return(
    list(e = e,
         D = D,
         H = H,
         w.matrix = w, # to complete later
         flow.con.mat = flow.con.mat)
  )
}




#' Creates an `sf` dataframe containing equally spaced points on a river network
#'
#' @param river_network An `sf` object of `MULTILINESTRING`(s) containing the streams and rivers, cropped to the sensor locations.
#' @param num_pred_sites An integer indicating the number of prediction sites to generate.
#' @return An `sf` dataframe of `POINT`(s) containing prediction sites.
#' @importFrom dplyr %>% rename filter mutate row_number
#' @importFrom sf st_sample st_is_empty st_as_sf st_geometry_type st_cast
#' @author Sean Francis
generate_prediction_locations <- function(river_network,
                                          num_pred_sites){

  sampled_points <- st_sample(river_network,
                              type = "regular",
                              size = num_pred_sites)

  x <- NULL
  geomtype <- NULL

  # Remove empty polygons
  sampled_points <- sampled_points[!st_is_empty(sampled_points)]

  d <- st_as_sf(sampled_points) %>%
    rename("geometry" = x)


  d$geomtype <- st_geometry_type(d)

  d1 <- d %>%
    filter(geomtype == "POINT")

  suppressWarnings({
    d2 <- d %>%
      filter(geomtype == "MULTIPOINT") %>%
      st_cast("POINT")
  })

  df <- rbind(d1, d2) %>%
    # Keep column location
    rename("locID" = geomtype) %>%
    mutate("locID" = paste0("pred_site_", row_number()))

  rownames(df) <- NULL


  return(df)

}










#' Spatio-Temporal Attention Regression Network (STARN)
#'
#' Fits a deep learning spatio-temporal model using LSTM layers with attention mechanisms and spatial embeddings.
#' This model is designed for detecting anomalies in high-frequency sensor data by capturing both spatial
#' dependencies (via distance-based embeddings) and temporal dynamics (via recurrent layers with attention).
#' Uncertainty is quantified via Monte Carlo dropout during inference.
#'
#' @param formula A formula specifying the response and predictors, e.g., \code{y ~ x1 + x2}.
#' @param data A data frame containing the spatio-temporal sensor data.
#' @param spat_matrix A spatial distance matrix between sensor locations (square and symmetric).
#' @param time_col Name of the column indicating temporal ordering (unquoted or quoted).
#' @param loc_col Name of the column indicating spatial location IDs (unquoted or quoted).
#' @param timesteps Length of the temporal window used to construct sequences for the LSTM model (default: 11).
#' @param spat_embedding_dim Dimension of the spatial embedding extracted via classical MDS (default: 10).
#' @param time_method Method used for temporal modeling (currently supports "lstm").
#' @param space_method Method used for spatial integration (currently supports "attention").
#' @param batch_size Training batch size (default: 64).
#' @param epochs Number of training epochs (default: 100).
#' @param validation_split Proportion of training data used for validation (default: 0.2).
#' @param n_iter Number of forward passes with dropout for uncertainty estimation (default: 100).

#' @importFrom keras keras_model compile fit
#' @importFrom tibble tibble
#' @importFrom dplyr pull ungroup
#' @importFrom stats predict sd quantile setNames

#' @return A list with components:
#' \describe{
#'   \item{\code{model}}{The trained Keras model object.}
#'   \item{\code{history}}{Training history returned by Keras.}
#'   \item{\code{data}}{The input data frame augmented with predictions, prediction intervals, and anomaly indicators.}
#' }
#'
#' @details
#' The function uses bidirectional LSTMs followed by an attention mechanism to extract temporal features,
#' and incorporates spatial information via embeddings obtained through classical multidimensional scaling (MDS).
#' Uncertainty is quantified using Monte Carlo dropout, and anomalies are flagged when observed values fall
#' outside the estimated 90% prediction interval.
#'
#' @references
#' Santos-Fernandez, E., et al. (2025). *New Bayesian and deep learning spatio-temporal models can reveal anomalies in sensor data more effectively*. (in review).
#'
#' @examples
#' \dontrun{
#' # Load example data from package
#' obs_path <- system.file("extdata", "obs_data.RDS", package = "SSNbayes")
#' dm_path <- system.file("extdata", "Cov_mat.rds", package = "SSNbayes")
#'
#' obs_data <- readRDS(obs_path)
#' dm <- readRDS(dm_path)
#' spat_matrix <- dm$H
#'
#' # Fit the STARN model
#' fit_starn <- starn(
#'   formula = y_impute ~ X1 + X2 + X3,
#'   data = obs_data,
#'   time_col = 'date',
#'   loc_col = 'locID',
#'   spat_matrix = spat_matrix,
#'   timesteps = 11,
#'   time_method = "lstm",
#'   space_method = "attention",
#'   batch_size = 64,
#'   epochs = 5,
#'   validation_split = 0.2
#' )
#' }
#'
#' @export


starn <- function(formula,
                  data,
                  spat_matrix,
                  time_col,
                  loc_col,
                  timesteps = 11,
                  spat_embedding_dim = 10,
                  time_method = "lstm",
                  space_method = "attention",
                  batch_size = 64,
                  epochs = 100,
                  validation_split = 0.2,
                  n_iter = 100) {

  # Parse response and predictors from formula
  response <- all.vars(formula)[1]
  predictors <- all.vars(formula)[-1]
  variables <- c(predictors, response)

  # Rename columns for flexibility based on time_col and loc_col arguments
  data <- data %>%
    rename(date = {{time_col}}, locID = {{loc_col}})

  # Prepare windowed data
  create_window_array <- function(data, window_size) {
    t(sapply(1:(length(data) - window_size + 1), function(x)
      data[x:(x + window_size - 1)]))
  }

  windowed_data <- data %>%
    group_by(locID) %>%
    arrange(date) %>%
    do({
      windows <- lapply(variables, function(var) {
        create_window_array(pull(., var), timesteps)
      })
      tibble(
        locID = unique(.$locID),
        windows = list(windows)
      )
    }) %>%
    ungroup()

  # Prepare input and output data
  n_samples <- sum(sapply(windowed_data$windows, function(x) nrow(x[[1]])))
  x_windows <- array(0, dim = c(n_samples, timesteps, length(variables) - 1))
  y_adjusted <- numeric(n_samples)
  locIDs <- numeric(n_samples)

  sample_index <- 1
  for (i in seq_along(windowed_data$windows)) {
    loc_windows <- windowed_data$windows[[i]]
    n_loc_samples <- nrow(loc_windows[[1]])
    for (j in 1:(length(variables) - 1)) {
      x_windows[sample_index:(sample_index + n_loc_samples - 1), , j] <- loc_windows[[j]]
    }
    y_adjusted[sample_index:(sample_index + n_loc_samples - 1)] <- loc_windows[[length(variables)]][, timesteps]
    locIDs[sample_index:(sample_index + n_loc_samples - 1)] <- windowed_data$locID[i]
    sample_index <- sample_index + n_loc_samples
  }

  unique_locIDs <- sort(unique(locIDs))
  locID_to_int <- setNames(seq_along(unique_locIDs), unique_locIDs)
  locIDs_int <- locID_to_int[as.character(locIDs)]
  num_locIDs <- length(unique_locIDs)

  # Split into train and test sets
  x_train <- x_windows
  y_train <- y_adjusted
  locIDs_train <- locIDs_int

  # Spatial embeddings using MDS
  embedding_dim <- spat_embedding_dim
  locID_embeddings_matrix <- cmdscale(spat_matrix, k = embedding_dim)
  rownames(locID_embeddings_matrix) <- as.character(unique_locIDs)
  locID_embeddings_train <- locID_embeddings_matrix[as.character(locIDs_train), ]

  # Define custom dropout layer
  DropoutAlways <- R6::R6Class(
    "DropoutAlways",
    inherit = KerasLayer,
    public = list(
      rate = NULL,
      seed = NULL,
      initialize = function(rate, seed = NULL) {
        self$rate <- rate
        self$seed <- seed
      },
      call = function(inputs, mask = NULL) {
        K <- backend()
        K$dropout(inputs, self$rate, seed = self$seed)
      }
    )
  )

  layer_dropout_always <- function(object, rate, seed = NULL, name = NULL, trainable = TRUE) {
    create_layer(DropoutAlways, object, list(rate = rate, seed = seed, name = name, trainable = trainable))
  }

  # Model architecture
  features <- length(predictors)
  input_series <- layer_input(shape = c(timesteps, features), name = 'series_input')
  input_embedding <- layer_input(shape = c(embedding_dim), name = 'embedding_input')

  lstm1 <- input_series %>%
    bidirectional(layer_lstm(units = 64, return_sequences = TRUE)) %>%
    layer_dropout_always(rate = 0.3)
  lstm2 <- lstm1 %>%
    bidirectional(layer_lstm(units = 64, return_sequences = TRUE)) %>%
    layer_dropout_always(rate = 0.3)

  attention_layer <- function(x) {
    a <- layer_dense(units = 1, activation = "tanh")(x)
    a_probs <- layer_activation(a, activation = "softmax")
    output_attention_mul <- layer_multiply(list(x, a_probs))
    output <- output_attention_mul %>%
      layer_lambda(f = function(z) k_sum(z, axis = 2))
    output
  }

  attention <- lstm2 %>%
    attention_layer()

  concat <- layer_concatenate(list(attention, input_embedding))
  dense1 <- concat %>%
    layer_dense(units = 32, activation = "relu") %>%
    layer_dropout_always(rate = 0.2)
  output <- dense1 %>%
    layer_dense(units = 1)

  model <- keras_model(inputs = list(series_input = input_series, embedding_input = input_embedding), outputs = output)

  # Compile and fit the model
  model %>% compile(optimizer = optimizer_adam(learning_rate = 0.001), loss = "mse")
  early_stopping <- callback_early_stopping(patience = 10, restore_best_weights = TRUE)
  reduce_lr <- callback_reduce_lr_on_plateau(factor = 0.2, patience = 5, min_lr = 1e-6)

  history <- model %>% fit(
    x = list(series_input = x_train, embedding_input = as.matrix(locID_embeddings_train)),
    y = y_train,
    epochs = epochs,
    batch_size = batch_size,
    validation_split = validation_split,
    callbacks = list(early_stopping, reduce_lr)
  )

  # Prediction with dropout
  predict_with_dropout <- function(model, x_input, n_iter) {
    preds <- replicate(n_iter, {
      predict(model, x_input)
    }, simplify = "array")

    preds <- array(preds, dim = c(dim(preds)[1], n_iter))
    preds_mean <- rowMeans(preds)
    preds_sd <- apply(preds, 1, sd)
    pred_low <- apply(preds, 1, function(x) quantile(x, 0.05))
    pred_up <- apply(preds, 1, function(x) quantile(x, 0.95))

    list(mean = preds_mean, sd = preds_sd, pred_low = pred_low, pred_up = pred_up)
  }

  K <- backend()
  K$set_learning_phase(1)  # Set learning phase to training

  x_train_input <- list(series_input = x_train, embedding_input = as.matrix(locID_embeddings_train))
  results_train <- predict_with_dropout(model, x_train_input, n_iter = n_iter)

  K$set_learning_phase(0)  # Reset learning phase to inference

  data$dataset <- 'test'
  data[1:nrow(x_train), ]$dataset <- 'train'

  data <- data %>% arrange(locID)
  data$pred <- data$pred_low <- data$pred_up <- NA

  data[data$date %in% timesteps:max(unique(data$date)), ]$pred = results_train$mean
  data[data$date %in% timesteps:max(unique(data$date)), ]$pred_low = results_train$pred_low
  data[data$date %in% timesteps:max(unique(data$date)), ]$pred_up = results_train$pred_up

  # Indicate anomalies
  data$anom_pred_ind_k <- ifelse(data$yobs < data$pred_low | data$yobs > data$pred_up, 1, 0)
  data$anom_pred_ind_k <- factor(data$anom_pred_ind_k, levels = c(1, 0))

  list(model = model, history = history, data = data)
}






#' Bayesian Autoregressive Spatio-Temporal model (BARST)
#'
#' Fits a Bayesian spatio-temporal model with autoregressive temporal structure and stream-network-based spatial modeling using Stan. The method automatically selects spatial knots based on iterative RMSE performance and returns predictions with uncertainty and anomaly indicators.
#'
#' @param formula A formula specifying the response and covariates, e.g., \code{yobs ~ X1 + X2 + X3}.
#' @param data A data frame with the full space-time dataset, including covariates and response.
#' @param path Path to the folder containing SSN or spatial network objects and distance matrices.
#' @param time_method A list specifying temporal structure, e.g., \code{list("ar", "date")}.
#' @param space_method A list specifying spatial covariance function and structure, e.g., \code{list("use_ssn", c("Exponential.taildown"))}.
#' @param iter Number of iterations per Stan run (default 3000).
#' @param warmup Number of warmup iterations (default 1500).
#' @param chains Number of Stan chains (default 3).
#' @param refresh Frequency of Stan progress printout (default 100).
#' @param loglik Logical, whether to return pointwise log-likelihood (default \code{FALSE}).
#' @param ppd Logical, whether to return posterior predictive distributions (default \code{TRUE}).
#' @param seed Random seed (default: 123).
#' @param total_iter Number of RMSE-based refinement iterations (default: 10).
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{fit}}{The final Stan model object.}
#'   \item{\code{pred}}{Data frame with posterior predictions, intervals, and anomaly indicators.}
#'   \item{\code{data_list}}{List of Stan inputs used for the final model.}
#'   \item{\code{replic}}{Matrix of knot selections across iterations.}
#' }
#'
#' @export
barst <- function(formula = yobs ~ X1 + X2 + X3,
                  data,
                  path,
                  time_method = list("ar", "date"),
                  space_method = list("use_ssn", c("Exponential.taildown")),
                  iter = 3000,
                  warmup = 1500,
                  chains = 3,
                  refresh = 100,
                  loglik = FALSE,
                  ppd = TRUE,
                  seed = 123,
                  total_iter = 10) {

  stopifnot(requireNamespace("rstan"))
  stopifnot(requireNamespace("dplyr"))
  stopifnot(requireNamespace("Matrix"))
  stopifnot(requireNamespace("stats"))


  obs_data <- readRDS(file.path(path, "obs_data.rds"))
  dm <- readRDS(file.path(path, "dist_mat.rds"))
  h <- dm$H
  fc <- dm$flow.con.mat

  # Add some missing values (if desired)
  obs_data[c(1000, 1500), "yobs"] <- NA

  # Create sampling probabilities
  fc_prob <- rowSums(fc) / sum(fc)

  n <- length(unique(obs_data$locID))
  time_points <- length(unique(obs_data$date))
  m <- 30
  set.seed(seed)

  locs <- sort(sample(unique(obs_data$locID), m, replace = FALSE, prob = fc_prob))

  D_star <- as.matrix(h[locs, locs])
  D_site_star <- as.matrix(h[, locs])

  options(na.action = 'na.pass')
  out_list <- mylm(formula = formula, data = obs_data)
  response <- out_list$y
  design_matrix <- out_list$X

  i_y_obs <- obs_data[!is.na(obs_data$yobs), ]$pid
  i_y_mis <- obs_data[is.na(obs_data$yobs), ]$pid
  y_obs <- response[!is.na(response)]

  # Priors
  K <- ncol(design_matrix)
  beta_mean <- rep(0, K)
  beta_sd <- rep(1, K)

  hyperpars <- list(
    eta = c(1, 1),
    alpha = c(5, 1),
    phi = c(0.5, 0.41),
    sigma_nug = c(1, 1),
    sigma = c(1, 1)
  )

  w_z_mean <- rep(1, m)
  w_z_sd <- rep(1, m)
  e_z_mean <- rep(1, n)
  e_z_sd <- rep(1, n)

  Xarray <- aperm(array(design_matrix, dim = c(n, time_points, K)), c(2, 1, 3))

  # Initialize replicate matrix
  replic <- matrix(0, m, total_iter + 2)
  replic[, 1] <- locs

  # Construct data list
  data_list <- list(
    obs_data = obs_data, T = time_points, replic = replic, niter = total_iter, iter_num = 1,
    n = n, m = m, y_obs = y_obs, N_y_obs = length(i_y_obs), N_y_mis = length(i_y_mis),
    i_y_obs = i_y_obs, i_y_mis = i_y_mis, D = as.matrix(h), K = K, X = Xarray,
    beta_mean = beta_mean, beta_sd = beta_sd,
    eta_mean = hyperpars$eta[1], eta_sd = hyperpars$eta[2],
    alpha_mean = hyperpars$alpha[1], alpha_sd = hyperpars$alpha[2],
    phi_mean = hyperpars$phi[1], phi_sd = hyperpars$phi[2],
    sigma_nug_mean = hyperpars$sigma_nug[1], sigma_nug_sd = hyperpars$sigma_nug[2],
    sigma_mean = hyperpars$sigma[1], sigma_sd = hyperpars$sigma[2],
    w_z_mean = w_z_mean, w_z_sd = w_z_sd,
    e_z_mean = e_z_mean, e_z_sd = e_z_sd
  )

  # Stan model
  stan_code <- readChar(system.file("stan", "barst.stan", package = "SSNbayes"), nchars = 1e6)

  for (jj in 10:total_iter) {
    message("Running iteration: ", jj)

    data_list$iter_num <- jj

    current_iter <- if (jj == total_iter) iter else 100
    current_warmup <- if (jj == total_iter) warmup else 50

    fit <- rstan::stan(
      model_code = stan_code,
      data = data_list,
      iter = current_iter,
      warmup = current_warmup,
      chains = chains,
      refresh = refresh
    )

  }

  list(fit = fit, pred = ypred, data_list = data_list, replic = replic)
}


