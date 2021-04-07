#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Collapses a SSN object into a data frame
#'
#' @param t Path to a SSN object
#' @return A data frame
#' @importFrom dplyr arrange
#' @importFrom plyr .
#' @export
#' @examples
#' #t.df <- collapse(t)
#'

collapse <- function(t){
  slot <- NULL
  df_all <- NULL
  for (i in 1:length(t@lines)){
    df <- data.frame(t@lines[[i]]@Lines[[1]]@coords)
    df$slot <- t@lines[[i]]@ID
    df$computed_afv <-  t@data$afvArea[i] #afvArea
    df_all<- rbind(df, df_all)

    df_all$slot <- as.numeric(as.character(df_all$slot))
  }
  df_all <-  dplyr::arrange(df_all, slot)
  df_all
}





#' Creates a list of distances and weights
#'
#' @param path Path to the files
#' @param net (optional) A network from the SSN object
#' @param addfunccol (optional) A parameter to compute the spatial weights
#' @return A list of matrices
#' @importFrom dplyr mutate %>% distinct left_join case_when
#' @importFrom plyr .
#' @importFrom SSN importSSN getSSNdata.frame
#' @importFrom rstan stan
#' @importFrom stats dist
#' @export
#' @examples
#' #mat_all <- dist_wei_mat(path)
#'


dist_wei_mat <- function(path = path, net = 1, addfunccol='addfunccol'){
  pid <- NULL
  n  <- importSSN(path, o.write = TRUE)
  obs_data <- getSSNdata.frame(n, "Obs")

  # creating distance matrices
  D <- readRDS(paste0(path, '/distance/obs/dist.net', net, '.RData')) # distance between observations

  # total distance
  H <- D + base::t(D)

  obs_data <- dplyr::filter(obs_data, pid %in% colnames(H))
  # NB replace here by the variable used for spatial weights
  afv <- obs_data[c('locID', addfunccol)] %>% distinct()

  # codes from SSN::glmssn
  nsofar <- 0
  dist.junc <- matrix(0, nrow = length(afv[, 1]), ncol = length(afv[,1]))
  distmat <- D
  ni <- length(distmat[1,])
  ordpi <- order(as.numeric(rownames(distmat)))
  dist.junc[(nsofar + 1):(nsofar + ni), (nsofar + 1):(nsofar + ni)] <-
    distmat[ordpi, ordpi, drop = F]
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

  coor <- n@obspoints@SSNPoints[[1]]@point.coords
  e <- coor %>%
    dist(., method = "euclidean", diag = FALSE, upper = FALSE) %>% as.matrix()

  #e <- obs_data %>%
  #  dplyr::select('coords.x1', 'coords.x2') %>%
  #  dist(., method = "euclidean", diag = FALSE, upper = FALSE) %>% as.matrix()

  list(e = e, D = D, H = H, w.matrix = w.matrix, flow.con.mat = flow.con.mat)
}




#' Creates a list of distances and weights between observed and prediction sites
#'
#' @param path Path with the name of the SSN object
#' @param net (optional) A network from the SSN object
#' @param addfunccol (optional) A parameter to compute the spatial weights
#' @return A list of matrices
#' @importFrom dplyr mutate %>% distinct left_join case_when
#' @importFrom plyr .
#' @importFrom SSN importSSN getSSNdata.frame
#' @importFrom rstan stan
#' @importFrom stats dist
#' @export
#' @description The output matrices are symmetric except the hydrologic distance matrix D.
#' @examples
#' #mat_all <- dist_wei_mat_preds(path, net = 1, addfunccol = 'addfunccol')
#'
#'

dist_wei_mat_preds <- function(path = path, net = 1, addfunccol = 'addfunccol'){
  netID <- NULL
  locID <- NULL
  pid <- NULL

  n  <- importSSN(path, predpts = 'preds', o.write = TRUE)

  obs_data <- getSSNdata.frame(n, "Obs")
  pred_data <- getSSNdata.frame(n, "preds")

  obs_data <- dplyr::filter(obs_data, netID == net ) %>% arrange(locID)
  pred_data <- dplyr::filter(pred_data, netID == net ) %>% arrange(locID)

  obs_data$locID_backup <- obs_data$locID
  pred_data$locID_backup <- pred_data$locID

  if(file.exists(paste0(path, '/distance/preds')) == F) stop("no distance matrix available between predictions. Please, use createDistMat(ssn_object, predpts = 'preds', o.write=TRUE, amongpreds = T)")

  # creating distance matrices
  doo <- readRDS(paste0(path, '/distance/obs/dist.net', net, '.RData'))

  if(file.exists(paste0(path, '/distance/preds/dist.net', net, '.a.RData')) == F) stop("no distance matrix available between observations and predictions. Please, use createDistMat(ssn_object, predpts = 'preds', o.write=TRUE, amongpreds = T)")
  dop <- readRDS(paste0(path, '/distance/preds/dist.net', net, '.a.RData')) # distance between observations
  dpo <- readRDS(paste0(path, '/distance/preds/dist.net', net, '.b.RData'))
  dpp <- readRDS(paste0(path, '/distance/preds/dist.net', net, '.RData'))

  D <- rbind(cbind(doo, dop), cbind(dpo, dpp))

  colnames(D)

  H <- D + base::t(D)

  print(dim(D))

  pred_data <- dplyr::filter(pred_data, pid %in% colnames(H))
  # NB replace here by the variable used for spatial weights
  afv1 <-  obs_data[c('locID', addfunccol)] %>% distinct()
  afv2 <- pred_data[c('locID', addfunccol)] %>% distinct()

  afv <- rbind(afv1, afv2)

  # codes from SSN::glmssn
  nsofar <- 0
  dist.junc <- matrix(0, nrow = length(afv[, 1]), ncol = length(afv[,1]))
  distmat <- D
  ni <- length(distmat[1,])
  ordpi <- order(as.numeric(rownames(distmat)))
  dist.junc[(nsofar + 1):(nsofar + ni), (nsofar + 1):(nsofar + ni)] <-
    distmat[ordpi, ordpi, drop = F]
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

  obs_data_coord <- data.frame(n@obspoints@SSNPoints[[1]]@point.coords)
  obs_data_coord$locID <- factor(1:nrow(obs_data_coord))
  obs_data_coord$locID <- as.numeric(as.character(obs_data_coord$locID))

  obs_data$locID <- as.numeric(factor(obs_data$locID))
  obs_data <- obs_data %>% left_join(obs_data_coord, by = c('locID'))
  obs_data$point <- 'Obs'

  pred_data_coord <- data.frame(n@predpoints@SSNPoints[[1]]@point.coords)
  pred_data_coord$locID <- factor( (max(as.numeric(as.character(obs_data_coord$locID))) + 1)
                                   :(nrow(pred_data_coord)+ (max(as.numeric(as.character(obs_data_coord$locID))) )))
  pred_data_coord$locID <- as.numeric(as.character(pred_data_coord$locID))

  pred_data$locID <- as.numeric(factor(pred_data$locID)) + max(obs_data$locID )
  pred_data <- pred_data %>% left_join(pred_data_coord, by = c('locID'))
  pred_data$point <- 'preds'

  all_data <- rbind(obs_data[c('coords.x1', 'coords.x2')],
                    pred_data[c('coords.x1', 'coords.x2')])

  e <- all_data %>%
    dplyr::select('coords.x1', 'coords.x2') %>%
    dist(., method = "euclidean", diag = FALSE, upper = FALSE) %>% as.matrix()
  colnames(e) <- colnames(D)
  rownames(e) <- rownames(D)

  list(e = e, D = D, H = H, w.matrix = w.matrix, flow.con.mat = flow.con.mat)
}

#mat_all_preds <- dist_wei_mat_preds(path = path,  net = 2, addfunccol='afvArea')





#' A simple modeling function using a formula and data
#'
#' @param formula A formula as in lm()
#' @param data A data.frame containing the elements specified in the formula
#' @return A list of matrices
#' @importFrom stats model.matrix model.response
#' @export
#' @author Jay ver Hoef
#' @examples
#' #options(na.action='na.pass')
#' #out_list = mylm(formula = y ~ X1 + X2 + X3, data = data)
#'

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




#' Fits a mixed regression model using Stan
#' It requires the same number of observation/locations per day.
#' It requires location id (locID) and points id (pid).
#' The locID are unique for each site.
#' The pid is unique for each observation.
#' Missing values are allowed in the response but not in the covariates.
#'
#' @param path Path with the name of the SSN object
#' @param formula A formula as in lm()
#' @param data A long data frame containing the locations, dates, covariates and the response variable. It has to have the locID and date. No missing values are allowed in the covariates.
#' @param space_method A list defining if use or not of an SSN object and the spatial correlation structure. The second element is the spatial covariance structure. A 3rd element is a a list with the lon and lat for Euclidean distance models.
#' @param time_method A list specifying the temporal structure (ar = Autorregressive; var = Vector autorregression) and coumn in the data with the time variable.
#' @param iter Number of iterations
#' @param warmup Warm up samples
#' @param chains Number of chains
#' @param refresh Sampler refreshing rate
#' @param net The network id (optional). Used when the SSN object cotains multiple networks.
#' @param addfunccol Variable to compute the additive function. Used to compute the spatial weights.
#' @param loglik Logic parameter denoting if the loglik will be computed by the model.
#' @param seed (optional) A seed for reproducibility
#' @return A list with the model fit
#' @export
#' @importFrom dplyr mutate %>% distinct left_join case_when
#' @importFrom plyr .
#' @importFrom SSN importSSN getSSNdata.frame
#' @importFrom rstan stan
#' @importFrom stats dist
#' @examples
#' #fit_td <- ssnbayes(formula = 'y ~ X1 + X2 + X3',
#' #                    data = data, ssn = ssn,
#' #                    space_method = list("use_ssn", "Exponential.taildown"),
#' #                    time_method = list("ar", "date"),
#' #                    iter = 3000,
#' #                    warmup = 1500,
#' #                    chains = 3)



# list('use_ssn', 'Exponential.tailup')
# list('no_ssn', 'Exponential.Euc', c('lon', 'lat'))

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
                     loglik = F,
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

  #if('date' %in% names(data) == F) stop("There is no column date on the data. Please, set a column called date with the time")
  if('locID' %in% names(data) == F) stop("There is no column locID on the data. Please, set a column called locID with the observation locations")

  if(missing(seed)) seed <- sample(1:1E6,1,replace=T)

  if(!missing(space_method)){
	print('using SSN object')
   if(space_method[[1]] == 'use_ssn'){
	  ssn_object <- T


	  if(length(space_method) > 1){
		 if(space_method[[2]] %in% c("Exponential.tailup", "LinearSill.tailup" , "Spherical.tailup" ,
		"Exponential.taildown" ,"LinearSill.taildown" ,"Spherical.taildown",
		"Exponential.Euclid") == F) {stop("Need to specify one or more of the following covariance matrices: Exponential.tailup, LinearSill.tailup , Spherical.tailup ,
		Exponential.taildown, LinearSill.taildown, Spherical.taildown or Exponential.Euc")}
		CorModels <- space_method[[2]]
		}
	  if(length(space_method) == 1){
		CorModels <- "Exponential.tailup"
		print('using an Exponential.tailup model')
	  }

  }
  if(space_method[[1]] == 'no_ssn'){
  print('no SSN object defined')
  ssn_object <- F
    if(space_method[[2]] %in% c("Exponential.Euc") == F) {stop("Need to specify Exponential.Euc")}
  # when using Euclidean distance, need to specify the columns with long and lat.
  if(length(space_method) < 3){ stop("Please, specify the columns in the data frame with the latitude and longitude (c('lon', 'lat'))") }

  data$lon <- data[,names(data) == space_method[[3]][1]]
  data$lat <- data[,names(data) == space_method[[3]][2]]
  CorModels <- space_method[[2]]
  }


  }

  if(missing(space_method)) {space_method <- 'no_ssn'; ssn_object <- F; CorModels <- "Exponential.Euc" }# if missing use Euclidean distance




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
  #CorModels == "Spherical.Euclid" ~ 2,
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
      //y[ i_y_obs[,t] ] = y_obs[((t - 1) * N_y_obs + 1):(t * N_y_obs)];
      //y[ i_y_mis[,t] ] = y_mis[((t - 1) * N_y_mis + 1):(t * N_y_mis)];
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
    //phi ~ uniform(-1, 1); //
    phi ~ normal(0.5,0.3); //NB informative
'

  model_tu <- '
    sigma_tu ~ uniform(0,100);  // cauchy(0,2) prior sd  tail-up model
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


  gen_quant <- '
  generated quantities {
   vector[T] log_lik;
   // vector[N] log_lik[T];
     for (t in 1:T){
       log_lik[t] = multi_normal_cholesky_lpdf(Y[t]|mu[t],
        cholesky_decompose(C_tu + C_td + C_re + C_ed + var_nug * I + 1e-6) );
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


    #ifelse(cor_tu %in% 1:3, tparam_tu2, 'C_tu = rep_matrix(0, N, N);'),
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

    if(loglik == T) gen_quant
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

    if(loglik == T) 'log_lik',

    'var_nug',
    'beta',
    'phi',
    'y'
  )
  pars <- pars[pars != '']


  # data part
  options(na.action='na.pass') # to preserve the NAs
  out_list <- mylm(formula = formula, data = data) # produces the design matrix

  response <- out_list$y # response variable
  design_matrix <- out_list$X # design matrix

  obs_data <- data
  #ndays <- length(unique(obs_data$date))

  ndays <- length(unique(obs_data[, names(obs_data) %in% time_points] ))


  N <- nrow(obs_data)/ndays #nobs


  nobs <- nrow(obs_data)/ndays #nobs

  #obs_data$date_num <- as.numeric(factor(obs_data$date))

obs_data$date_num <- as.numeric(factor(obs_data[, names(obs_data) %in% time_points]	))


  resp_var_name <- gsub("[^[:alnum:]]", " ", formula[2])
  obs_data$y <- obs_data[,names(obs_data) %in% resp_var_name]

  #train <- nrow(obs_data[obs_data$date_num==1 & !is.na(obs_data$y),])
  #test <- nrow(obs_data[obs_data$date_num==1 & is.na(obs_data$y),])



  # array structure
  X <- design_matrix #cbind(1,obs_data[, c("X1", "X2", "X3")]) # design matrix
  #Xarray <- abind(matrix(X[,1],nrow = ndays, ncol = nobs, byrow = T),
  #                matrix(X[,2],nrow = ndays, ncol = nobs, byrow = T),
  #                matrix(X[,3],nrow = ndays, ncol = nobs, byrow = T),
  #                matrix(X[,4],nrow = ndays, ncol = nobs, byrow = T), along = 3)

  # NB: this array order is Stan specific
  Xarray <- aperm(array( c(X), dim=c(N, ndays, ncol(X)) ),c(2, 1, 3))


  y_obs <- response[!is.na(response)]#obs_data[!is.na(obs_data$y),]$y
  #Yarray <- array(matrix(ys, nrow = ndays, ncol = train, byrow = T),
  #                dim = c(ndays, train))

  # index for observed values
  i_y_obs <- obs_data[!is.na(obs_data$y),]$pid
  # i_y_obs <- matrix(i_y_obs, nrow = train , ncol = ndays, byrow = F)


  # index for missing values
  i_y_mis <- obs_data[is.na(obs_data$y),]$pid
  #  i_y_mis <- matrix(i_y_mis, nrow = test, ncol = ndays, byrow = F)



  if(ssn_object == T){ # the ssn object exist?
    mat_all <- dist_wei_mat(path = path, net = net, addfunccol = addfunccol)
  }

  if(ssn_object == F){ # the ssn object does not exist- purely spatial

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

  # phi_ini <- ifelse(time_method == "ar", 0.5, rep(0.5,N))


  ini <- function(){list(var_nug =  .1#, phi = phi_ini #phi= rep(0.5,N)
                         #y = rep( mean(obs_data$temp, na.rm = T),T*N)
  )}


  fit <- rstan::stan(model_code = ssn_ar,
                     model_name = "ssn_ar",
                     data = data_list,
                     pars = pars,
                     iter = iter,
                     warmup = warmup,
                     init = ini,
                     chains = chains,
                     verbose = F,
                     seed = seed,
                     refresh = refresh
  )
  attributes(fit)$formula <- formula

  fit
}



#' Performs spatial prediction in R using a stanfit object from ssnbayes()
#' It will take an observed and a prediction data frame.
#' It requires the same number of observation/locations per day.
#' It requires location id (locID) and points id (pid).
#' The locID are unique for each site.
#' The pid is unique for each observation.
#' Missing values are allowed in the response but not in the covariates.
#'
#' @param stanfit A stanfit object returned from ssnbayes
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
#' @importFrom SSN importSSN getSSNdata.frame
#' @importFrom rstan stan
#' @importFrom stats dist
#' @importFrom stats as.formula

krig <- function(stanfit = stanfit,
                 mat_all_preds = mat_all_preds,
                 nsamples = 10,
                 start = 1,
                 chunk_size = 50,
                 obs_data = obs_data,
                 pred_data = pred_data,
                 net = net,
				 seed = seed){

  if(missing(seed)) seed <- sample(1:1E6,1,replace=T)
  set.seed(seed)

  formula <- as.formula(attributes(stanfit)$formula)

  phi <- rstan::extract(stanfit, pars = 'phi')$phi

  samples <- sample(1:nrow(phi), nsamples, replace = F)
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

  locID_pred <- sort(unique(pred_data$locID)) #6422

  pred_data$locID0 <- as.numeric(factor(pred_data$locID)) # conseq locID
  pred_data$locID0 <- pred_data$locID0 + length(locID_obs) # adding numb of locID in obs dataset

  locID_pred0 <- sort(unique(pred_data$locID0))

  # obs data frame. no missing in temp

  locID_pred_1 <- locID_pred0[start:(start + chunk_size - 1)] # NB

  pred_data_1 <- pred_data[pred_data$locID0 %in% locID_pred_1,]

  #dim(pred_data_1)

  N <- nrow(pred_data_1)

  #   formula = temp ~ SLOPE + elev + CUMDRAINAG + air_temp + sin + cos

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
  #print('kriging')

  for(k in 1:niter){

    #print(k)

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



#' Performs spatial prediction in R using a stanfit object from ssnbayes()
#' It will take an observed and a prediction data frame.
#' It requires the same number of observation/locations per day.
#' It requires location id (locID) and points id (pid).
#' The locID are unique for each site.
#' The pid is unique for each observation.
#' Missing values are allowed in the response but not in the covariates.
#'
#' @param path Path with the name of the SSN object
#' @param obs_data The observed data frame
#' @param stanfit A stanfit object returned from ssnbayes
#' @param pred_data The predicted data frame
#' @param net (optional) Network from the SSN object
#' @param nsamples The number of samples to draw from the posterior distributions. (nsamples <= iter)
#' @param addfunccol The variable used for spatial weights
#' @param chunk_size (optional) the number of locID to make prediction from
#' @param locID_pred (optional) the location id for the predictions. Used when the number of pred locations is large.
#' @param seed (optional) A seed for reproducibility

#' @return A data frame
#' @export
#' @importFrom dplyr mutate %>% distinct left_join case_when
#' @importFrom plyr .
#' @importFrom SSN importSSN getSSNdata.frame
#' @importFrom rstan stan
#' @importFrom stats dist
#' @examples
#' #pred <- pred_ssnbayes(path = path,
#' #obs_data = clear,
#' #stanfit = fit_ar,
#' #pred_data = preds,
#' #net = 2,
#' #nsamples = 100, # number of samples to use from the posterior in the stanfit object
#' #addfunccol = 'afvArea') # variable used for spatial weights



pred_ssnbayes <- function(
  path = path,
  obs_data = obs_data,
  stanfit = stanfit,
  pred_data = pred_data,
  net = 1,
  nsamples = 100, # number of samples to use from the posterior in the stanfit object
  addfunccol = 'afvArea', # variable used for spatial weights
  locID_pred = locID_pred, # location ID of the points to predict
  #ssn_object = T # if it will use an SSN object
  #CorModels = "Exponential.tailup",
  chunk_size = chunk_size,
  seed = seed
){

  if(missing(seed)) seed <- sample(1:1E6,1,replace=T)
  set.seed(seed)


  mat_all_preds <- dist_wei_mat_preds(path = path,
                                      net = net,
                                      addfunccol = addfunccol)

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

  #chunk_size <- length(unique(obs_data$locID))

  if(missing(chunk_size)) chunk_size <- pred_points

  out_all <- NULL

  is <- ceiling(pred_points/chunk_size)

  for(j in 1:is){ #NB
    print('krig_funct')
    print(j)
    start <- ((j - 1) * chunk_size + 1)

    chunk_size <- ifelse(j != is, chunk_size, pred_points - (j - 1) * chunk_size)

    out <- krig(stanfit = stanfit,
                mat_all_preds = mat_all_preds,
                nsamples = nsamples,
                start = start,
                chunk_size = chunk_size,
                obs_data = obs_data,
                pred_data = pred_data,
                net = 2)

    out_all <- rbind(out_all, out)
  }

  data.frame(out_all)

}

