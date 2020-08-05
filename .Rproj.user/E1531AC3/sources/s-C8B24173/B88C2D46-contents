# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

#' Creates a list of distances and weights
#'
#' @param path Path to the files
#' @return A list of matrices
#' @importFrom dplyr mutate %>% distinct left_join case_when
#' @importFrom plyr .
#' @importFrom SSN importSSN getSSNdata.frame
#' @importFrom rstan stan
#' @importFrom stats dist
#' @export
#' @examples
dist_wei_mat <- function(path = path){

  n  <- importSSN(path, o.write = TRUE)
  obs_data <- getSSNdata.frame(n, "Obs")

  # creating distance matrices
  D <- readRDS(paste0(path, '/distance/obs/dist.net1.RData')) # distance between observations

  # total distance
  H <- D + base::t(D)

  # NB replace here by the variable used for spatial weights
  afv <- obs_data[c('locID', 'addfunccol')] %>% distinct()

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
  w.matrix <- sqrt(pmin(outer(afv[, 'addfunccol'],rep(1, times = n.all)),
                        base::t(outer(afv[, 'addfunccol'],rep(1, times = n.all) ))) /
                     pmax(outer(afv[, 'addfunccol'],rep(1, times = n.all)),
                          base::t(outer(afv[, 'addfunccol'], rep(1, times = n.all))))) *
    flow.con.mat

  # Euclidean distance

  obs_data_coord <- data.frame(n@obspoints@SSNPoints[[1]]@point.coords)
  obs_data_coord$locID <- factor(1:nrow(obs_data_coord))

  obs_data <- obs_data %>% left_join(obs_data_coord, by = c('locID'))
  obs_data$point <- 'Obs'


  e <- obs_data %>%
    dplyr::select('coords.x1', 'coords.x2') %>%
    dist(., method = "euclidean", diag = FALSE, upper = FALSE) %>% as.matrix()

  list(e = e, D = D, H = H, w.matrix = w.matrix, flow.con.mat = flow.con.mat)
}

#mat_all <- dist_wei_mat(path)


#' Fits the model using Stan
#'
#' @param formula A formula as in lm()
#' @param data A list
#' @param CorModels Correlation structure
#' @param iter Number of iterations
#' @param warmup Warm up samples
#' @param chains Number of chains
#' @param refresh Sampler refreshing rate
#' @return A list with the fit
#' @export
#' @importFrom dplyr mutate %>% distinct left_join case_when
#' @importFrom plyr .
#' @importFrom SSN importSSN getSSNdata.frame
#' @importFrom rstan stan
#' @importFrom stats dist
#' @examples
#' #fit_td <- ssnbayes(formula = 'y ~ X1 + X2 + X3',
#' #                     data = data,
#' #                    CorModels = "Exponential.taildown",
#' #                     iter = 3000,
#' #                     warmup = 1500,
#' #                    chains = 3)

ssnbayes <- function(formula = formula,
                     data = data,
                     CorModels = "Exponential.tailup",
                     iter = 3000,
                     warmup = 1500,
                     chains = 3,
                     refresh = max(iter/100, 1)){

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

    int<lower = 1> i_y_obs[N_y_obs,T] ;
    int<lower = 1> i_y_mis[N_y_mis,T] ;

    vector[N_y_obs * T] y_obs;    //matrix[N_y_obs,1] y_obs[T];

    matrix[N, N] W ; // spatial weights
    matrix[N, N] h ; // total hydrological dist
    matrix[N, N] I ; // diag matrix

    matrix[N, N] D ; // downstream hydrological dist matrix
    matrix[N, N] flow_con_mat; // flow conected matrix

    matrix[N, N] e ; // Euclidean dist mat

  }'


  param_com <- '
  parameters {
    vector[K] beta;
    real<lower=0> sigma_nug;

    vector[N_y_mis * T] y_mis;//declaring the missing y
    real phi;
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

   real<lower=0> alpha_max;

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
    alpha_max = 4 * max(h);
    for (t in 1:T){
      y[ i_y_mis[,t] ] = y_mis[((t - 1) * N_y_mis + 1):(t * N_y_mis)];
      y[ i_y_obs[,t] ] = y_obs[((t - 1) * N_y_obs + 1):(t * N_y_obs)];
      Y[t] = y[((t - 1) * N + 1):(t * N)];
    }

    var_nug = sigma_nug ^ 2; // variance nugget
        mu[1] = X[1] * beta;
    epsilon[1] = Y[1] - mu[1];

    for (t in 2:T){
        mu[t] = X[t] * beta;
        epsilon[t] = Y[t] - mu[t];
        mu[t] = mu[t] + phi * epsilon[t-1];
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
      target += multi_normal_cholesky_lpdf(Y[t] | mu[t], cholesky_decompose(C_tu + C_td + C_re + C_ed + var_nug * I) );
    }

    sigma_nug ~ cauchy(0,1); // prior nugget effect
    phi ~ uniform(-1, 1); //
'

  model_tu <- '
    sigma_tu ~ cauchy(0,2);  // prior sd  tail-up model
    alpha_tu ~ uniform(0, alpha_max);
'

  model_td <- '
    sigma_td ~ cauchy(0,2); // sd tail-down
    alpha_td ~ uniform(0, alpha_max);
'

  model_ed <- '
    sigma_ed ~ cauchy(0,2); // sd Euclidean dist
    alpha_ed ~ uniform(0, alpha_max); // Euclidean dist range
'

  model_re <- '
    sigma_RE1 ~ uniform(0,5);
'


  ssn_ar1 <- paste(
    data_com,

    param_com,

    if(cor_tu %in% 1:3) param_tu,
    if(cor_td %in% 1:3) param_td,
    if(cor_ed %in% 1:3) param_ed,
    if(cor_re %in% 1:3) param_re,
    '}',

    tparam_com,
    if(cor_tu %in% 1:3)tparam_tu,
    if(cor_td %in% 1:3)tparam_td,
    if(cor_ed %in% 1:3)tparam_ed,
    if(cor_re %in% 1:3) tparam_re,
    tparam_com2,

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
    '}'
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
    'var_nug',
    'beta',
    'phi',
    'y'
  )
  pars <- pars[pars != '']


  data$e = data$mat_all$e #Euclidean dist
  #for tail-up
  data$h = data$mat_all$H # total stream distance
  data$W = data$mat_all$w.matrix # spatial weights

  #for tail-down
  data$flow_con_mat = data$mat_all$flow.con.mat #flow connected matrix
  data$D = data$mat_all$D #downstream hydro distance matrix

  #RE1 = RE1mm # random effect matrix

  data$I = diag(1, nrow(data$W), nrow(data$W))  # diagonal matrix

  fit <- stan(model_code = ssn_ar1,
              model_name = "ssn_ar1",
              data = data,
              pars = pars,
              iter = iter,
              warmup = warmup,
              chains = chains,
              verbose = F,
              #seed = 22,
              refresh = refresh
  )

  fit
}

