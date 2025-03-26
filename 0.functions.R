
#basic Rao's Q function

raoq <- function(props, dist){
  c(crossprod(props, dist %*% props))
  }


#computerized resampling function

resample <- function(df, steps, length){
  if (!is.numeric(df[[1]]) & !is.character(df[[2]])){
    stop("Check column types")
  }
  if (any(steps>length)){
    stop("maximum step must be less than length")
  }
  if (ncol(df)!=2){
    stop("data must have exactly two columns")
  }
  copy <- df
  copy[1] <- df[1]+length
  df <- rbind(df, copy) #concatenates two transects to allow "wrapping around"
  output <- vector(mode = "list", length = length(steps))
  names(output) <- steps
  for (j in seq_along(steps)){
    resampled <- data.frame()
    for (i in seq_len(length)){
    range <- i:(i+steps[j]-1)
    records <- df[df[[1]] %in% range, 2, drop = T]
    extend <- data.frame(rep(i, length(records)), records)
    resampled <- rbind(resampled, extend)
    }
    names(resampled) <- names(df)
    output[[j]] <- resampled
  }
  return(output)
}


#Q decomposition and equivalent numbers

qdecomp <- function(data, distances, species){
  cont <- table(factor(data$species, levels = species), data$X) |> #contingency table of abundances
    proportions(margin = 2) |> 
    as.data.frame.matrix() |> 
    select(where(~ !any(is.nan(.)))) #drops empty samples
  alpha <- sapply(cont, function(x) raoq(x, distances)) |> #alpha Q
    mean()
  P_i <- apply(cont, 1, mean) #regional relative abundances as mean of local ab.
  gamma <- raoq(P_i, distances) #gamma Q
  beta_add <- gamma - alpha #additive beta
  beta_rel <- beta_add / gamma #additive beta normalized by gamma
  E_alpha <- 1/(1-alpha) #equivalent alpha
  E_gamma <- 1/(1-gamma) #equivalent gamma
  E_beta_mult <- E_gamma / E_alpha #multiplicative eq. beta
  E_beta_add <- E_gamma - E_alpha #additive eq. beta
  E_beta_add_prop_norm <- (E_beta_add / E_gamma) / (1 - 1 / ncol(cont)) #proportional additive beta, normalized on sample n
  taxon <- matrix(1, nrow = dim(distances)[1], ncol = dim(distances)[2])-
    diag(1, dim(distances)[1], dim(distances)[2]) #matrix of taxonomic distances (0/1)
  alpha_sp <- sapply(cont, function(x) raoq(x, taxon)) |> #taxonomic alpha
    mean()
  gamma_sp <- raoq(P_i, taxon) #taxonomic gamma
  beta_sp_add <- gamma_sp - alpha_sp #additive tax. beta
  beta_sp_rel <- beta_sp_add / gamma_sp #normalized tax. beta
  E_alpha_sp <- 1/(1-alpha_sp) #equivalent tax. alpha
  E_gamma_sp <- 1/(1-gamma_sp) #equivalent tax. gamma
  E_beta_sp <- E_gamma_sp / E_alpha_sp #equivalent tax. beta
  clustering <- beta_rel / beta_sp_rel #beta Q clustering
  uniqueness <- gamma / gamma_sp #gamma Q uniqueness
  U_gamma_star <- E_gamma / E_gamma_sp #gamma E uniqueness
  redundancy_a <- E_alpha_sp / E_alpha #alpha E redundancy, inverse of Ualpha*
  output <- data.frame(alpha, gamma, beta_add, beta_rel,
                       E_alpha, E_gamma, E_beta_mult, E_beta_add,
                       E_beta_add_prop_norm, clustering, uniqueness,
                       alpha_sp, gamma_sp, U_gamma_star, redundancy_a,
                       beta_sp_add, E_alpha_sp, E_gamma_sp, E_beta_sp)
  return(output)
}

#distribution plot function

histdensity <- function(data, var, bins = 10, adjust = 1)  {
  data |> 
    ggplot(aes(x = {{var}}))+
    geom_histogram(aes(y=after_stat(density)), color = "white", bins=bins)+
    geom_density(aes(y=after_stat(density)), color = "red", adjust = 1)
}

#Box-Cox transformation

boxcox <- function(x) {
  if (any(x == 0)) x <- x + 1

  boxcox_result <- MASS::boxcox(x ~ 1)
  lambda <- boxcox_result$x[which.max(boxcox_result$y)]
  
  message("Box-Cox transform with lambda = ", lambda)
  
  if (lambda == 0) {
    x <- log(x)
  } else {
    x <- (x ^ lambda - 1) / lambda
  }
  return(x)
}
