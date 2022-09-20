calculate_amp_phi <- function (a_cos, b_sin) {
  amp <- unname(sqrt(a_cos^2 + b_sin^2))
  phi <- unname(atan2(b_sin, a_cos) %% (2*pi))
  cbind(amp = amp, phi = phi)
}

calculate_ci_amp_phi <- function (amp, a_cos, b_sin, fit_res_ssr, R, deg_f) {
  # if (det(ssx) == 0 | kappa(ssx) > 0.5/.Machine$double.eps)
  #   return(cbind(amp = NA, phi = NA))
  ## manually compute variance-covariance matrix (faster here than vcov())
  ssxinv <- chol2inv(R)
  ## need only rows 2 and 3 corresponding to cos and sin, according to
  ## the delta rule
  ssxinvab <- ssxinv[2:3, 2:3]
  ## Jacobians stored row-wise in a matrix
  Ja <- cbind(ifelse(amp > 0, a_cos/amp, 0), 
              ifelse(amp > 0, b_sin/amp, 0))
  Jp <- cbind(ifelse(amp > 0, -b_sin/amp^2, 0), 
              ifelse(amp > 0, a_cos/amp^2, 0))
  ## dimension of non-matrix vector x will automatically be adjusted by R
  ## to conform with matrix multiplication, see ?matmult
  var_a <- apply(Ja, 1, function (x) x %*% ssxinvab %*% x) * fit_res_ssr/deg_f
  var_p <- apply(Jp, 1, function (x) x %*% ssxinvab %*% x) * fit_res_ssr/deg_f
  ##ciquant <- qnorm(0.025, lower.tail = FALSE)
  ciquant <- 2
  ci_amp <- sqrt(as.numeric(var_a))*ciquant
  ci_phi <- sqrt(as.numeric(var_p))*ciquant
  
  cbind(amp = ci_amp, phi = ifelse(ci_phi < pi, ci_phi, pi))
  
}

calculate_ci_amp_phi_r <- function (amp, a_cos, b_sin, ssx, thefit) {
  # if (det(ssx) == 0 | kappa(ssx) > 0.5/.Machine$double.eps)
  #   return(cbind(amp = NA, phi = NA))
  vcovmat <- vcov(thefit)[2:3, 2:3]
  Ja <- cbind(ifelse(amp > 0, a_cos/amp, 0), 
              ifelse(amp > 0, b_sin/amp, 0))
  Jp <- cbind(ifelse(amp > 0, -b_sin/amp^2, 0), 
              ifelse(amp > 0, a_cos/amp^2, 0))
  var_a <- Ja %*% vcovmat %*% t(Ja)
  var_p <- Jp %*% vcovmat %*% t(Jp)
  ##ciquant <- qnorm(0.025, lower.tail = FALSE)
  ciquant <- 2
  ci_amp <- sqrt(as.numeric(var_a))*ciquant
  ci_phi <- sqrt(as.numeric(var_p))*ciquant
  
  cbind(amp = ci_amp, phi = ifelse(ci_phi < pi, ci_phi, pi))
  
}

noncentral_f_test <- function (Fval, deg_f1, deg_f2, X, a_over_sigma) {
  # where X is the N x 3 design matrix.
  X1 <- X[, -(2:3), drop = FALSE]
  ## model matrix harmonics
  X2 <- X[, 2:3]
  ## see e.g., A.C. Davison, Statistical Models, p. 367
  H1 <- X1 %*% solve(crossprod(X1)) %*% t(X1)
  Z2 <- (diag(dim(X2)[1]) - H1) %*% X2
  
  ## smallest eigenvalue of Z2 (equal to N/2 if balanced)
  mineig <- min(eigen(crossprod(Z2), 
                      only.values = TRUE, symmetric = TRUE)$values)
  
  ## p value according to noncentral F distribution
  delsq <- a_over_sigma^2*mineig
  ## note that this function could be useful also in the robust context; the
  ## 3*pi correction (see below) is not necessarily needed here, since this
  ## function would in the context of semiparametric reduction of dispersion be
  ## used only with upper tail (test against weak or zero rhythmicity compound
  ## null), for which higher ncp makes the test slightly more conservative
  stats::pf(Fval, df1 = deg_f1, df2 = deg_f2, ncp = delsq)
  
}

noncentral_chisq_test <- function (x2, deg_f, X, a_over_sigma) {
  # where X is the N x 3 design matrix.
  X1 <- X[, -(2:3), drop = FALSE]
  ## model matrix harmonics
  X2 <- X[, 2:3]
  ## see e.g., A.C. Davison, Statistical Models, p. 367
  H1 <- X1 %*% solve(crossprod(X1)) %*% t(X1)
  Z2 <- (diag(dim(X2)[1]) - H1) %*% X2
  
  ## smallest eigenvalue of Z2 (equal to N/2 if balanced)
  mineig <- min(eigen(crossprod(Z2), only.values = TRUE)$values)
  
  ## p value according to noncentral F distribution
  ## see Hettmansperger & McKean (2011) page 203 for factor 3/pi
  delsq <- a_over_sigma^2*mineig*3/pi
  pchisq(x2, df = deg_f, ncp = delsq)
  
}


# normalize time series matrix (no NAs) -----------------------------------
normalize_ts_matrix <- function (inputts, inputtime, 
                                 norm.pol, norm.pol.degree) {
  if (norm.pol) {
    if (length(inputtime) < (1 + norm.pol.degree)) {
      stop(paste("Normalization polynomial degree too high for the number",
                 "of time points"))
    }
    trendfit <- lm(inputts ~ poly(inputtime, norm.pol.degree, raw = TRUE))
    trend.ts <- fitted(trendfit)
    trend.coef <- coef(trendfit)
    if(any(zapsmall(trend.ts) == 0)) {
      stop(paste("Normalization using polynomial failed (zero-crossing);", 
                 "try other normalization settings"))
    }
    return(list(norm.ts = inputts/trend.ts, norm.w = trend.coef, 
                norm.vals = trend.ts))
  }
  else {
    tsmeans <- colMeans(inputts)
    if(any(zapsmall(tsmeans) == 0)) {
      stop(paste("Some time series have zero means, normalization failed.",
                 "Try other normalization settings, or leave these data out."))
    }
    return(list(norm.ts = scale(inputts, FALSE, tsmeans), norm.w = tsmeans))
  }
}



# harmonic regression matrix (no NAs) -------------------------------------
harmonic_regression_matrix <- function (inputts, inputtime, Tau,
                                        a_over_sigma) {
  
  ## check time series length
  if (length(inputtime) < 3) {
    stop(paste("These time series are too short for a meaningful analysis.",
               "At least 3 time points are needed."))
  }
  
  ## matrix fit of the unrestricted model (harmonic regression)
  inputts.fit <- lm(inputts ~ 1 + cos(2*pi/Tau*inputtime) + 
                      sin(2*pi/Tau*inputtime), x = TRUE)
  
  ## covariance matrix of independent variables
  ssx <- zapsmall(crossprod(inputts.fit$x))
  # if (det(ssx) == 0 | kappa(ssx) > 0.5/.Machine$double.eps)
  #   return(cbind(amp = NA, phi = NA))
  ## refrain from parameter estimation if the design matrix is bad
  if (det(ssx) == 0 || 
      (log10(kappa(ssx)) > (-log10(.Machine$double.eps) - 4))) {
    stop(paste("The time points are so unfortunately spaced that a phase and",
               "amplitude determination is impossible."))
  }
  
  ## fitted values, possibly coerce to matrix
  fit.vals <- as.matrix(fitted(inputts.fit))
  
  ## coefficients, amplitudes, phases
  coeffs <- t(coef(inputts.fit))
  pars <- as.data.frame(calculate_amp_phi(coeffs[, 2], coeffs[, 3]))
  if (!any(duplicated(colnames(inputts)))) {
    rownames(pars) <- colnames(inputts)
  }
  
  ## if more than 3 time points are available, confidence intervals and p-values
  ## can be computed.
  if (length(inputtime) == 3) {
    
    pvals <- NA
    ci <- NA
    fit.res.ssr <- NA
    pvals_son <- NA
    # pvals_son_null_weak <- NA
    
  } else {
    
    ## sum squared residual of the restricted and unrestricted models
    # inputts.ssr <- apply(inputts, 2, var)*(length(inputtime) - 1)
    
    # if (is.matrix(residuals(inputts.fit)))
    #   fit.res.ssr <- colSums(residuals(inputts.fit)^2)
    # ## handle also the case with one single input ts
    # else
    #   fit.res.ssr <- sum(residuals(inputts.fit)^2)
    
    fit.res.ssr <- deviance(inputts.fit)
    
    inputts.fit.summaries <- summary(inputts.fit)
    names(inputts.fit.summaries) <- 
      gsub("Response ", "", names(inputts.fit.summaries))
    
    ## f-statistic and pvalues
    
    fstats <- as.data.frame(t(sapply(inputts.fit.summaries, 
                                     function (x) x$fstatistic)))
    pvals <- with(fstats, stats::pf(value, numdf, dendf, lower.tail = FALSE))
    names(pvals) <- names(inputts.fit.summaries)
    
    pvals_son <- 
      with(fstats, 
           noncentral_f_test(value, numdf, dendf, inputts.fit$x, a_over_sigma))
    names(pvals_son) <- names(inputts.fit.summaries)
    
    # pvals_son_null_weak <- 1 - pvals_son
    
    ## distance to upper confidence interval limit
    ci <- as.data.frame(calculate_ci_amp_phi(pars$amp, 
                                             coeffs[, 2], coeffs[, 3],
                                             fit.res.ssr,
                                             qr.R(inputts.fit$qr),
                                             length(inputtime) - 3))
    if (!any(duplicated(colnames(inputts))))
      rownames(ci) <- colnames(inputts)
    
    
  }
  
  deg_f <- length(inputtime) - 3
  
  ## return values
  list(fit.vals = fit.vals,
       pars = pars, pvals = pvals, 
       ci = ci, coeffs = coeffs[, 2:3], 
       ssr = fit.res.ssr, df = deg_f, sigma_hat = sqrt(fit.res.ssr/deg_f),
       ssx = ssx,
       # pvals_son_null_weak = pvals_son_null_weak,
       pvals_son = pvals_son)
  
}


harmonic_regression_matrix_nuisance <- function (inputts, inputtime, Tau,
                                                 nuisance_f,
                                                 a_over_sigma) {
  
  ## "move" the formula to a new formula with local environment
  nuisance_f_local <- stats::formula(deparse(nuisance_f))
  parent.env(environment(nuisance_f_local)) <- environment(nuisance_f)
  
  nuisance_dim <- ncol(stats::model.matrix(nuisance_f))
  ## check for enough degrees of freedom
  deg_f <- length(inputtime) - (nuisance_dim + 2)
  if (deg_f < 0) {
    stop(paste("Too few time points.  Unable to continue.  Try fewer nuisance",
               "variables."))
  }
  
  ## matrix fit of the restricted model (polynomial)
  rest.fit <- stats::lm(stats::update(nuisance_f_local, inputts ~ .))
  
  ## matrix fit of the unrestricted model (harmonic regression)
  unrest.fit <- stats::update(rest.fit, . ~ 
                                cos(2*pi/Tau*inputtime) + 
                                sin(2*pi/Tau*inputtime) + .,
                              x = TRUE)
  ## possibly coerce fitted values to matrix
  fit.vals <- as.matrix(stats::fitted(unrest.fit))
  
  ## refrain from parameter estimation if the design matrix is bad
  ssx <- zapsmall(crossprod(unrest.fit$x))
  if (det(ssx) == 0 || 
      (log10(kappa(ssx)) > (-log10(.Machine$double.eps) - 4))) {
    stop(paste("The time points are so unfortunately spaced that a phase and",
               "amplitude determination is impossible."))
  }
  
  ## coefficients, amplitudes, phases
  coeffs <- t(stats::coef(unrest.fit))
  pars <- as.data.frame(calculate_amp_phi(coeffs[, 2], coeffs[, 3]))
  if (!any(duplicated(colnames(inputts)))) {
    rownames(pars) <- colnames(inputts)
  }
  
  if (deg_f == 0) {
    
    pvals <- NA
    ci <- NA
    unrest.ssr <- NA
    pvals_son <- NA
    # pvals_son_null_weak <- NA
    
  } else {
    
    unrest.ssr <- stats::deviance(unrest.fit)
    rest.ssr <- stats::deviance(rest.fit)
    
    ## f-statistic and pvalues (anova() is not useful here)
    fstats <- ((rest.ssr - unrest.ssr)/2) / (unrest.ssr/deg_f)
    pvals <- stats::pf(fstats, 2, deg_f, lower.tail = FALSE)
    
    pvals_son <- noncentral_f_test(fstats, 2, deg_f, unrest.fit$x, 
                                   a_over_sigma)
    names(pvals_son) <- names(unrest.ssr)
    
    # pvals_son_null_weak <- 1 - pvals_son
    
    
    ## distance to upper confidence interval limit 
    ci <- as.data.frame(calculate_ci_amp_phi(pars$amp, 
                                             coeffs[, 2], coeffs[, 3],
                                             unrest.ssr, qr.R(unrest.fit$qr),
                                             deg_f))
    if (!any(duplicated(colnames(inputts)))) {
      rownames(ci) <- colnames(inputts)
    }
    
  }
  
  
  ## return values
  list(fit.vals = fit.vals,
       pars = pars, pvals = pvals, 
       ci = ci, coeffs = coeffs, 
       ssr = unrest.ssr, df = deg_f, sigma_hat = sqrt(unrest.ssr/deg_f), 
       ssx = ssx,
       # pvals_son_null_weak = pvals_son_null_weak,
       pvals_son = pvals_son)
  
}


# normalize time series vector (with NAs) ---------------------------------
normalize_one_ts <- function (inputts, inputtime, 
                              norm.pol, norm.pol.degree) {
  n.non.na <- length(which(!is.na(inputts)))
  if (norm.pol) {
    if (n.non.na < (1 + norm.pol.degree)) {
      return(list(norm.ts = rep(NA, length(inputtime)), 
                  norm.w  = rep(NA, 1 + norm.pol.degree),
                  norm.vals = rep(NA, length(inputtime))))
    }
    trendfit <- lm(inputts ~ poly(inputtime, norm.pol.degree, raw = TRUE),
                   na.action = na.exclude)
    trend.ts <- fitted(trendfit)
    trend.coef <- coef(trendfit)
    if(any(zapsmall(trend.ts) == 0, na.rm = TRUE)) {
      warning(paste("Zero-crossing of normalization polynomial.", 
                    "One of the time series is left out the fitting procedure"))
      return(list(norm.ts = rep(NA, length(inputts)), norm.w = trend.coef,
                  norm.vals = trend.ts))
    }
    return(list(norm.ts = inputts/trend.ts, norm.w = trend.coef,
                norm.vals = trend.ts))
  }
  else {
    if (n.non.na == 0) {
      return(list(norm.ts = inputts, norm.w = NA))
    }
    tsmean <- mean(inputts, na.rm = TRUE)
    if(abs(tsmean) < 10^-getOption("digits")) {
      warning(paste("Mean value zero of a time series.  It is left out of the",
                    "fitting procedure"))
      return(list(norm.ts = rep(NA, length(inputts)), norm.w = trend.ts))
    }
    return(list(norm.ts = inputts/tsmean, norm.w = tsmean))
  }
}

## this is only for the case that multiplicative normalization is done for 
## polynomials.  "mean" normalization is done as part of the main robust fit.
normalize_one_ts_r_pol <- function (inputts, inputtime, norm.pol.degree,
                                    robust_scores) {
  n.non.na <- length(which(!is.na(inputts)))
  if (n.non.na < (1 + norm.pol.degree)) {
    return(list(norm.ts = rep(NA, length(inputtime)), 
                norm.w  = rep(NA, 1 + norm.pol.degree),
                norm.vals = rep(NA, length(inputtime))))
  }
  trendfit <- try(Rfit::rfit(inputts ~ poly(inputtime, norm.pol.degree, 
                                            raw = TRUE),
                             na.action = na.exclude, scores = robust_scores),
                  silent = TRUE)
  if (inherits(trendfit, "try-error")) {
    warning(paste("Failure of robust normalization trend fitting.", 
                  "One of the time series is left out the fitting procedure"))
    return(list(norm.ts = rep(NA, length(inputts)), 
                norm.w = rep(NA, 1 + norm.pol.degree),
                norm.vals = rep(NA, length(inputtime))))
  }
  if (n.non.na < length(inputts)) {
    trend.ts <- rep(NA, length(inputts))
    trend.ts[!is.na(inputts)] <- fitted(trendfit)
  } else {
    trend.ts <- fitted(trendfit)
  }
  trend.coef <- coef(trendfit)
  if(any(zapsmall(trend.ts) == 0, na.rm = TRUE)) {
    warning(paste("Zero-crossing of normalization polynomial.", 
                  "One of the time series is left out the fitting procedure"))
    return(list(norm.ts = rep(NA, length(inputts)), norm.w = trend.coef,
                norm.vals = trend.ts))
  }
  
  list(norm.ts = inputts/trend.ts, norm.w = trend.coef,
       norm.vals = trend.ts)
  
}


# harmonic regression one time series (with NAs) --------------------------
fit_one_harmonic <- function (inputts, inputtime, Tau, a_over_sigma) {
  
  n.non.na <- length(which(!is.na(inputts)))
  if (n.non.na < 3) {
    return(list(pars = c(amp = NA, phi = NA),
                coeffs = rep(NA, 3), ci = c(amp = NA, phi = NA),
                fit.vals = rep(NA, length(inputtime)),
                ssr = NA, deg_f = NA, sigma_hat = NA, 
                ssx = NA,
                pval = NA,
                # pval_son_null_weak = NA,
                pval_son = NA))
  }
  ## harmonic regression
  inputts.fit <- lm(inputts ~ (1 + cos(2*pi/Tau*inputtime) + 
                                 sin(2*pi/Tau*inputtime)),
                    na.action = na.exclude, x = TRUE) 
  
  ## refrain from parameter estimation if the design matrix is bad
  ssx <- zapsmall(crossprod(inputts.fit$x))
  if (det(ssx) == 0 || 
      (log10(kappa(ssx)) > (-log10(.Machine$double.eps) - 4))) {
    return(list(pars = c(amp = NA, phi = NA),
                coeffs = rep(NA, 3), ci = c(amp = NA, phi = NA),
                fit.vals = rep(NA, length(inputtime)),
                ssr = NA, deg_f = NA, sigma_hat = NA,
                ssx = ssx,
                pval = NA,
                # pval_son_null_weak = NA,
                pval_son = NA))
  }
  
  fit.vals <- fitted(inputts.fit)
  
  coeffs <- coef(inputts.fit)
  pars <- calculate_amp_phi(coeffs[2], coeffs[3])
  
  if (n.non.na == 3) {
    
    pval <- NA
    ci <- c(amp = NA, phi = NA)
    fit.res.ssr <- NA
    pval_son <- NA
    # pval_son_null_weak <- NA
    
  } else {
    
    ## sum squared residual of the unrestricted model; fitted values
    fit.res.ssr <- deviance(inputts.fit)
    
    ## f-statistic and pvalues
    inputts_fstat <- summary(inputts.fit)$fstatistic
    pval <- unname(
      stats::pf(inputts_fstat["value"], 
                inputts_fstat["numdf"], inputts_fstat["dendf"],
                lower.tail = FALSE)
    )
    
    pval_son <- 
      unname(noncentral_f_test(inputts_fstat["value"], 
                               inputts_fstat["numdf"], 
                               inputts_fstat["dendf"], 
                               inputts.fit$x, a_over_sigma))
    # pval_son_null_weak <- 1 - pval_son
    
    ## distance to upper confidence interval limit 
    ci <- calculate_ci_amp_phi(pars[, "amp"], coeffs[2], coeffs[3],
                               fit.res.ssr, qr.R(inputts.fit$qr),
                               n.non.na - 3)
    
    
  }
  
  deg_f <- n.non.na - 3
  
  list(pars = pars,
       coeffs = coeffs, ci = ci,
       fit.vals = fit.vals,
       ssr = fit.res.ssr, deg_f = deg_f, sigma_hat = sqrt(fit.res.ssr/deg_f),
       ssx = ssx,
       pval = pval,
       # pval_son_null_weak <- 1 - pval_son
       pval_son = pval_son)
  
}


fit_one_harmonic_nuisance <- function (inputts, inputtime, Tau,
                                       nuisance_f, a_over_sigma) {
  
  nuisance_f_local <- stats::formula(deparse(nuisance_f))
  parent.env(environment(nuisance_f_local)) <- environment(nuisance_f)
  
  n.non.na <- length(which(!is.na(inputts)))
  nuisance_dim <- ncol(stats::model.matrix(nuisance_f))
  ## check for enough degrees of freedom
  deg_f <- n.non.na - (nuisance_dim + 2)
  if (deg_f < 0) {
    return(list(pars = c(amp = NA, phi = NA),
                coeffs = rep(NA, nuisance_dim + 2), ci = c(amp = NA, phi = NA),
                fit.vals = rep(NA, length(inputtime)),
                ssr = NA, deg_f = NA, sigma_hat = NA,
                ssx = NA,
                pval = NA,
                # pval_son_null_weak = NA,
                pval_son = NA))
  }
  
  ## fit of the restricted model (nuisance_f)
  rest.fit <- stats::lm(stats::update(nuisance_f_local, inputts ~ .), 
                        na.action = na.exclude)
  
  ## matrix fit of the unrestricted model (harmonic regression)
  unrest.fit <- stats::update(rest.fit, . ~ 
                                cos(2*pi/Tau*inputtime) + 
                                sin(2*pi/Tau*inputtime) + .,
                              x = TRUE, na.action = na.exclude)
  
  ## refrain from parameter estimation if the design matrix is bad
  ssx <- zapsmall(crossprod(unrest.fit$x))
  if (det(ssx) == 0 || 
      (log10(kappa(ssx)) > (-log10(.Machine$double.eps) - 4))) {
    return(list(pars = c(amp = NA, phi = NA),
                coeffs = rep(NA, nuisance_dim + 2), ci = c(amp = NA, phi = NA),
                fit.vals = rep(NA, length(inputtime)),
                ssr = NA, deg_f = NA, sigma_hat = NA,
                ssx = ssx, 
                pval = NA,
                # pval_son_null_weak = NA,
                pval_son = NA))
  }
  
  fit.vals <- fitted(unrest.fit)
  
  ## coefficients, amplitudes, phases
  coeffs <- coef(unrest.fit)
  pars <- calculate_amp_phi(coeffs[2], coeffs[3])
  
  if (deg_f == 0) {
    
    pval <- NA
    ci <- c(amp = NA, phi = NA)
    unrest.ssr <- NA
    pval_son <- NA
    # pval_son_null_weak <- NA
    
  } else {
    
    ## sum squared residual of the restricted and unrestricted models
    rest.ssr <- stats::deviance(rest.fit)
    unrest.ssr <- stats::deviance(unrest.fit)
    
    ## f-statistic and pvalues
    test <- stats::anova(rest.fit, unrest.fit)
    pval <- test$`Pr(>F)`[2]
    
    pval_son <- 
      unname(noncentral_f_test(test$F[2], 
                               test$Df[2], 
                               test$Res.Df[2], unrest.fit$x, a_over_sigma))
    # pval_son_null_weak <- 1 - pval_son
    
    
    ## distance to upper confidence interval limit (Halberg 1967)
    ci <- calculate_ci_amp_phi(pars[, "amp"], coeffs[2], coeffs[3],
                               unrest.ssr, qr.R(unrest.fit$qr), deg_f)
    
  }    
  
  list(pars = pars,
       coeffs = stats::coef(unrest.fit), ci = ci,
       fit.vals = fit.vals,
       ssr = unrest.ssr, deg_f = deg_f, sigma_hat = sqrt(unrest.ssr/deg_f),
       ssx = ssx,
       # pval_son_null_weak = pval_son_null_weak,
       pval = pval, pval_son = pval_son)
}


# robust harmonic regression one time series (with NAs) -------------------
fit_one_harmonic_r <- function (inputts, inputtime, Tau, normalize = FALSE,
                                a_over_sigma, robust_scores) {
  
  n.non.na <- length(which(!is.na(inputts)))
  if (n.non.na < 3) {
    return(list(pars = c(amp = NA, phi = NA),
                coeffs = rep(NA, 3), ci = c(amp = NA, phi = NA),
                fit.vals = rep(NA, length(inputtime)),
                ssr = NA, deg_f = NA, sigma_hat = NA,
                ssx = NA,
                pval = NA,
                # pval_son_null_weak = NA,
                pval_son = NA,
                mean = NA))
  }
  
  # inputtime <- times
  # inputts <- rhythm.data[5, ]
  # Tau <- 24
  
  inputts.fit <- try(Rfit::rfit(inputts ~ (1 + cos(2*pi/Tau*inputtime) + 
                                             sin(2*pi/Tau*inputtime)),
                                na.action = na.exclude,
                                scores = robust_scores), silent = TRUE)
  if (inherits(inputts.fit, "try-error")) {
    warning(paste("Robust fitting procedure failed in one case. NAs", 
                  "are reported for this case."))
    return(list(pars = c(amp = NA, phi = NA),
                coeffs = rep(NA, 3), ci = c(amp = NA, phi = NA),
                fit.vals = rep(NA, length(inputtime)),
                ssr = NA, deg_f = NA, sigma_hat = NA,
                ssx = NA,
                pval = NA,
                # pval_son_null_weak = NA,
                pval_son = NA,
                mean = NA))
  }
  
  
  ## refrain from parameter estimation if the design matrix is bad
  ssx <- zapsmall(crossprod(inputts.fit$x))
  if (det(ssx) == 0 || 
      (log10(kappa(ssx)) > (-log10(.Machine$double.eps) - 4))) {
    return(list(pars = c(amp = NA, phi = NA),
                coeffs = rep(NA, 3), ci = c(amp = NA, phi = NA),
                fit.vals = rep(NA, length(inputtime)),
                ssr = NA, deg_f = NA, sigma_hat = NA,
                ssx = NA,
                pval = NA,
                # pval_son_null_weak = NA,
                pval_son = NA,
                mean = NA))
  }
  
  ## workaround for a bug in rfit(), where NAs are not propagated by fitted() 
  ## although na.exclude() is used
  if (n.non.na < length(inputts)) {
    fit.vals <- rep(NA, length(inputts))
    fit.vals[!is.na(inputts)] <- fitted(inputts.fit)
  } else {
    fit.vals <- fitted(inputts.fit)
  }
  
  
  coeffs <- coef(inputts.fit)
  pars <- calculate_amp_phi(coeffs[2], coeffs[3])
  
  ## compute robust mean by averaging over the design matrix 
  ## (see for instance Davison pp. 383--384)
  mean_r <- as.numeric(crossprod(coeffs, colMeans(inputts.fit$x)))
  
  if (n.non.na == 3) {
    
    pval <- NA
    ci <- c(amp = NA, phi = NA)
    fit.res.ssr <- NA
    sigma_hat <- NA
    pval_son <- NA
    # pval_son_null_weak <- NA,
    if (normalize) {
      pars[, "amp"] <- pars[, "amp"]/mean_r
      fit.vals <- fit.vals/mean_r
    }
    
    
  } else {
    
    ## p value
    rfit_summary <- try(
      Rfit::summary.rfit(inputts.fit, overall.test = "drop"),
      silent = TRUE)
    
    
    ## there may be conditions for which summary.rfit() fails
    if (inherits(rfit_summary, "try-error")) {
      warning(paste("The robust testing procedure against the null hypothesis", 
                    "did not converge for one",
                    "sample.  NA is reported for this case"))
      pval <- NA
      pval_son <- NA
      # pval_son_null_weak <- NA
      
      #   return(list(pars = c(amp = NA, phi = NA),
      #               coeffs = rep(NA, 3), ci = c(amp = NA, phi = NA),
      #               fit.vals = rep(NA, length(inputtime)),
      #               ssr = NA, deg_f = NA, pval = NA,
      #               pval_son = NA))
    } else {
      pval <- as.numeric(rfit_summary$droppval)
      # NOTE: summary.rfit returns F statistic, i.e. with reduction in 
      # dispersion in numerator being divided by additional 2 degrees of 
      # freedom. Therefore multiply by 2 to get the chi-square statistic
      df_ssq_red <- length(inputts.fit$coefficients) - 1
      pval_son <- 
        noncentral_chisq_test(as.numeric(rfit_summary$dropstat)*df_ssq_red, 
                              df_ssq_red,
                              inputts.fit$x,
                              a_over_sigma)
      # pval_son_null_weak <- 1 - 
      #   noncentral_f_test(as.numeric(rfit_summary$dropstat), 
      #                     df_ssq_red,
      #                     length(inputts.fit$y) - inputts.fit$qrx1$rank,
      #                     inputts.fit$x,
      #                     a_over_sigma)
      
    }
    
    ## sum squared residual of the unrestricted model; fitted values
    fit.res.ssr <- sum(residuals(inputts.fit)^2, na.rm = TRUE)
    
    ## sigma_hat, see Hettmansperger & McKean (2011) page 203
    sigma_hat <- inputts.fit$tauhat*sqrt(3/pi)
    
    ## distance to upper confidence interval limit 
    ci <- calculate_ci_amp_phi_r(pars[, "amp"], coeffs[2], coeffs[3],
                                 ssx, inputts.fit)
    
    if (normalize) {
      ## if we normalize, we normalize here, to avoid an extra earlier 
      ## run of Rfit just to compute the robust mean.
      fit.res.ssr <- fit.res.ssr/(mean_r^2)
      sigma_hat <- sigma_hat/mean_r
      ci[, "amp"] <- ci[, "amp"]/mean_r
      pars[, "amp"] <- pars[, "amp"]/mean_r
      fit.vals <- fit.vals/mean_r
      coeffs[2:3] <- coeffs[2:3]/mean_r
    }
    
  }
  
  list(pars = pars,
       coeffs = coeffs, ci = ci,
       mean = mean_r,
       fit.vals = fit.vals,
       ssr = unname(fit.res.ssr), deg_f = (n.non.na - 3), 
       sigma_hat = unname(sigma_hat),
       ssx = ssx,
       # pval_son_null_weak = pval_son_null_weak,
       pval = pval, pval_son = pval_son)
  
}


fit_one_harmonic_nuisance_r <- function (inputts, inputtime, Tau,
                                         nuisance_f, normalize = FALSE, 
                                         a_over_sigma, robust_scores) {
  
  nuisance_f_local <- stats::formula(deparse(nuisance_f))
  parent.env(environment(nuisance_f_local)) <- environment(nuisance_f)
  
  n.non.na <- length(which(!is.na(inputts)))
  nuisance_dim <- ncol(stats::model.matrix(nuisance_f))
  ## check for enough degrees of freedom
  deg_f <- n.non.na - (nuisance_dim + 2)
  if (deg_f < 0) {
    return(list(pars = c(amp = NA, phi = NA),
                coeffs = rep(NA, nuisance_dim + 2), ci = c(amp = NA, phi = NA),
                fit.vals = rep(NA, length(inputtime)),
                ssr = NA, deg_f = NA, sigma_hat = NA,
                ssx = NA,
                pval = NA,
                # pval_son_null_weak = NA, 
                pval_son = NA,
                mean = NA))
  }
  
  ## fit of the restricted model (nuisance_f)
  rest_f <- stats::update(nuisance_f_local, inputts ~ .)
  rest.fit <- try(
    Rfit::rfit(rest_f, na.action = na.exclude, scores = robust_scores), 
    silent = TRUE
  )
  
  ## matrix fit of the unrestricted model (harmonic regression)
  unrest_f <- stats::update(rest_f, ~ cos(2*pi/Tau*inputtime) + 
                              sin(2*pi/Tau*inputtime) + .)
  unrest.fit <- try(
    Rfit::rfit(unrest_f, na.action = na.exclude, scores = robust_scores), 
    silent = TRUE
  )
  if (inherits(rest.fit, "try-error") || inherits(unrest.fit, "try-error")) {
    warning(paste("Robust fitting procedure failed in one case. NAs", 
                  "are reported for this case."))
    return(list(pars = c(amp = NA, phi = NA),
                coeffs = rep(NA, nuisance_dim + 2), ci = c(amp = NA, phi = NA),
                fit.vals = rep(NA, length(inputtime)),
                ssr = NA, deg_f = NA, sigma_hat = NA,
                ssx = NA,
                pval = NA,
                # pval_son_null_weak = NA, 
                pval_son = NA,
                mean = NA))
  }
  
  ## refrain from parameter estimation if the design matrix is bad
  ssx <- zapsmall(crossprod(unrest.fit$x))
  if (det(ssx) == 0 || 
      (log10(kappa(ssx)) > (-log10(.Machine$double.eps) - 4))) {
    return(list(pars = c(amp = NA, phi = NA),
                coeffs = rep(NA, nuisance_dim + 2), ci = c(amp = NA, phi = NA),
                fit.vals = rep(NA, length(inputtime)),
                ssr = NA, deg_f = NA, sigma_hat = NA,
                ssx = ssx,
                pval = NA,
                # pval_son_null_weak = NA, 
                pval_son = NA,
                mean = NA))
  }
  
  ## workaround for a bug in rfit(), where NAs are not propagated by fitted() 
  ## although na.exclude() is used
  if (n.non.na < length(inputts)) {
    fit.vals <- rep(NA, length(inputts))
    fit.vals[!is.na(inputts)] <- stats::fitted(unrest.fit)
  } else {
    fit.vals <- stats::fitted(unrest.fit)
  }
  
  
  ## coefficients, amplitudes, phases
  coeffs <- stats::coef(unrest.fit)
  pars <- calculate_amp_phi(coeffs[2], coeffs[3])
  
  ## compute robust mean with by averaging over the design matrix 
  ## (see for instance Davison pp. 383--384)
  mean_r <- as.numeric(crossprod(coeffs, colMeans(unrest.fit$x)))
  
  if (deg_f == 0) {
    
    pval <- NA
    ci <- c(amp = NA, phi = NA)
    unrest.ssr <- NA
    sigma_hat <- NA
    pval_son <- NA
    # pval_son_null_weak <- NA 
    if (normalize) {
      pars[, "amp"] <- pars[, "amp"]/mean_r
      fit.vals <- fit.vals/mean_r
    }
    
  } else {
    
    unrest.ssr <- sum(stats::residuals(unrest.fit)^2, na.rm = TRUE)
    
    ## sigma_hat, see Hettmansperger & McKean (2011) page 203
    sigma_hat <- unrest.fit$tauhat*sqrt(3/pi)
    
    
    ## p value
    rfit_testresult <- try(Rfit::drop.test(unrest.fit, rest.fit),
                           silent = TRUE)
    
    ## there may be conditions for which rfit() fails
    if (inherits(rfit_testresult, "try-error")) {
      warning(paste("The robust testing procedure against the null hypothesis", 
                    "did not converge for one",
                    "sample.  NA is reported for this case"))
      pval <- NA
      pval_son <- NA
      # pval_son_null_weak <- NA 
      
      #   return(list(pars = c(amp = NA, phi = NA),
      #               coeffs = rep(NA, 3), ci = c(amp = NA, phi = NA),
      #               fit.vals = rep(NA, length(inputtime)),
      #               ssr = NA, deg_f = NA, pval = NA,
      #               pval_son = NA))
    } else {
      pval <- as.numeric(rfit_testresult$p.value)
      pval_son <- 
        noncentral_chisq_test(as.numeric(rfit_testresult$F)*rfit_testresult$df1, 
                              rfit_testresult$df1,
                              unrest.fit$x,
                              a_over_sigma)
      # pval_son_null_weak <- 1 - 
      #   noncentral_f_test(as.numeric(rfit_testresult$F), 
      #                     rfit_testresult$df1,
      #                     rfit_testresult$df2,
      #                     unrest.fit$x,
      #                     a_over_sigma)
      
    }
    
    
    ## distance to upper confidence interval limit (Halberg 1967)
    ci <- calculate_ci_amp_phi_r(pars[, "amp"], coeffs[2], coeffs[3],
                                 ssx, unrest.fit)
    
    if (normalize) {
      ## if we normalize, we normalize here, to avoid an extra earlier
      ## run of Rfit just to compute the robust mean.
      unrest.ssr <- unrest.ssr/(mean_r^2)
      sigma_hat <- sigma_hat/mean_r
      ci[, "amp"] <- ci[, "amp"]/mean_r
      pars[, "amp"] <- pars[, "amp"]/mean_r
      fit.vals <- fit.vals/mean_r
      coeffs <- coeffs/mean_r
    }
    
    
  }    
  
  list(pars = pars,
       coeffs = coeffs, ci = ci,
       mean = mean_r,
       fit.vals = fit.vals,
       ssr = unname(unrest.ssr), deg_f = deg_f, sigma_hat = unname(sigma_hat),
       ssx = ssx,
       # pval_son_null_weak = pval_son_null_weak, 
       pval = pval, pval_son = pval_son)
  
}



# harmonic regression one time series at a time, with NAs -----------------
harmonic_regression_nas <- function (inputts, inputtime, Tau, robust = FALSE,
                                     normalize = FALSE, 
                                     a_over_sigma,
                                     n_cores, robust_scores = robust_scores) {
  
  if (length(inputtime) < 3) {
    stop(paste("These time series are too short for a meaningful analysis.",
               "At least 3 time points are needed"))
  }
  
  # For now, we don't allow future.apply because of unresolved bugs
  # if (n_cores == 1L || !requireNamespace("future.apply", quietly = TRUE)) {
  if (!robust) {
    result.list <- apply(inputts, 2, fit_one_harmonic, inputtime, Tau,
                         a_over_sigma)
  } else {
    result.list <- apply(inputts, 2, fit_one_harmonic_r, inputtime, Tau,
                         normalize = normalize,
                         a_over_sigma, robust_scores)
  }
  # } else {
  #   future::plan(future::multisession,
  #                workers = min(n_cores, future::availableCores()))
  #   if (!robust) {
  #     result.list <- 
  #       future.apply::future_apply(inputts, 2, fit_one_harmonic, inputtime, 
  #       Tau, a_over_sigma)
  #   } else {
  #     result.list <- 
  #       future.apply::future_apply(inputts, 2, fit_one_harmonic_r, 
  #                                  inputtime, Tau, normalize = normalize,
  #                                  a_over_sigma, robust_scores)
  #   }
  # }
  
  coeffs <- t(sapply(result.list, "[[", "coeffs"))
  fit.vals <- sapply(result.list, "[[", "fit.vals")
  ssr <- sapply(result.list, "[[", "ssr")
  deg_f <- sapply(result.list, "[[", "deg_f")
  sigma_hat <- sapply(result.list, "[[", "sigma_hat")
  ssx <- lapply(result.list, "[[", "ssx")
  pvals <- sapply(result.list, "[[", "pval")
  pvals_son <- sapply(result.list, "[[", "pval_son")
  # pvals_son_null_weak <- 
  #   sapply(result.list, "[[", "pval_son_null_weak")
  pars <- t(sapply(result.list, "[[", "pars"))
  colnames(pars) <- c("amp", "phi")
  if (!any(duplicated(colnames(inputts)))) {
    rownames(pars) <- colnames(inputts)
  }
  pars <- as.data.frame(pars)
  ci <- t(sapply(result.list, "[[", "ci"))
  colnames(ci) <- c("amp", "phi")
  if (!any(duplicated(colnames(inputts)))) {
    rownames(ci) <- colnames(inputts)
  }
  ci <- as.data.frame(ci)
  
  return_val <- 
    list(fit.vals = fit.vals,
         pars = pars, pvals = pvals,
         ci = ci, coeffs = coeffs, ssr = ssr, df = deg_f, sigma_hat = sigma_hat,
         ssx = ssx, 
         # pvals_son_null_weak = pvals_son_null_weak,
         pvals_son = pvals_son)
  
  if (robust) {
    c(list(means = sapply(result.list, "[[", "mean")), return_val)
  } else {
    return_val
  }
  
}


harmonic_regression_nas_nuisance <- function (inputts, inputtime, Tau,
                                              nuisance_f, robust = FALSE,
                                              normalize = FALSE, a_over_sigma,
                                              n_cores, robust_scores =
                                                robust_scores) {
  
  nuisance_dim <- ncol(stats::model.matrix(nuisance_f))
  ## check for enough degrees of freedom
  deg_f <- length(inputtime) - (nuisance_dim + 2)
  if (deg_f < 0) {
    stop(paste("Too few time points. Unable to continue. Try fewer nuisance",
               "variables."))
  }
  
  # if (n_cores == 1L || !requireNamespace("future.apply", quietly = TRUE)) {
  if (!robust) {
    result.list <- apply(inputts, 2, fit_one_harmonic_nuisance, inputtime, 
                         Tau, nuisance_f, a_over_sigma)
  } else {
    result.list <- apply(inputts, 2, fit_one_harmonic_nuisance_r, inputtime, 
                         Tau, nuisance_f, normalize = normalize, a_over_sigma,
                         robust_scores)
  }
  # } else {
  #   future::plan(future::multisession, 
  #                workers = min(n_cores, future::availableCores()))
  #   if (!robust) {
  #     result.list <- 
  #       future.apply::future_apply(inputts, 2, fit_one_harmonic_nuisance, 
  #                                  inputtime, Tau, nuisance_f, a_over_sigma)
  #   } else {
  #     result.list <- 
  #       future.apply::future_apply(inputts, 2, fit_one_harmonic_nuisance_r, 
  #                                  inputtime, Tau, nuisance_f, 
  #                                  normalize = normalize, a_over_sigma,
  #                                  robust_scores)
  #   }
  # }
  
  coeffs <- t(sapply(result.list, "[[", "coeffs"))
  fit.vals <- sapply(result.list, "[[", "fit.vals")
  ssr <- sapply(result.list, "[[", "ssr")
  deg_f <- sapply(result.list, "[[", "deg_f")
  sigma_hat <- sapply(result.list, "[[", "sigma_hat")
  ssx <- lapply(result.list, "[[", "ssx")
  pvals <- sapply(result.list, "[[", "pval")
  pvals_son <- sapply(result.list, "[[", "pval_son")
  # pvals_son_null_weak <- 
  #   sapply(result.list, "[[", "pval_son_null_weak")
  pars <- t(sapply(result.list, "[[", "pars"))
  colnames(pars) <- c("amp", "phi")
  if (!any(duplicated(colnames(inputts)))) {
    rownames(pars) <- colnames(inputts)
  }
  pars <- as.data.frame(pars)
  ci <- t(sapply(result.list, "[[", "ci"))
  colnames(ci) <- c("amp", "phi")
  if (!any(duplicated(colnames(inputts)))) {
    rownames(ci) <- colnames(inputts)
  }
  ci <- as.data.frame(ci)
  
  return_val <- 
    list(fit.vals = fit.vals,
         pars = pars, pvals = pvals,
         ci = ci, coeffs = coeffs, ssr = ssr, df = deg_f, sigma_hat = sigma_hat,
         ssx = ssx, 
         # pvals_son_null_weak = pvals_son_null_weak,
         pvals_son = pvals_son)
  
  if (robust) {
    c(list(means = sapply(result.list, "[[", "mean")), return_val)
  } else {
    return_val
  }
  
}

# Documentation noncentral test for (stronger) rhythmicity (not used at the
# moment):
# \code{pvals_son_null_weak} \tab Vector of p-values according to a 
# noncentral F-test for rhythmicity against a compound null hypothesis of
# arrhythmicity or weak rhythmicity, defined by the cutoff 


# @param n_cores Numeric or integer, determines whether parallel processing
#   will be attempted. This requires the \code{future.apply} package installed,
#   and only applies if there are NAs in the data, and/or \code{robust = TRUE}.
#   If \code{n_cores} is greater than the number of available cores, the number
#   of available cores is used instead. If \code{n_cores = 1L} (default),
#   and/or if the \code{future.apply} is not available, no parallel processing
#   is performed.


# harmonic regression documentation ---------------------------------------
#' Harmonic Regression
#'
#' Estimates amplitudes and phases along with confidence intervals and p-values 
#' from a set of time series that may oscillate with a specified period. A 
#' model, per default \deqn{y = m + a cos(\omega t) + b sin(\omega t),} is 
#' fitted to the time series.  This model is equivalent to the model \deqn{m + c
#' cos(\omega t - \phi),} with amplitude \eqn{c = \sqrt(a^2 + b^2)} and phase 
#' \eqn{\phi = atan2(b, a)}. P-values for \eqn{c > 0} (more precisely: either 
#' \eqn{a} or \eqn{b > 0} ) are computed by an F-test.  Confidence intervals for
#' the amplitudes and phases are computed by a linear error propagation 
#' approximation.
#' 

#' The default setting is that the time series are normalized with their mean
#' values.  Optionally a polynomial of degree 1 or more is first fitted to each
#' time series, whereupon the original time series are normalized by dividing
#' with the fitted values at each point, thus trends in a fold-change sense are
#' assumed.  Another option is to include nuisance variables, in which case the
#' same model plus any nuisance variables: \eqn{y = m + a cos(\omega t) + b
#' sin(\omega t) + (...) } is fitted to the (possibly normalized) data.  In this
#' case, returned p-values still only concern the alternative \eqn{c > 0} as
#' defined above.

#'
#' There is a robust option, in which a robust method rather than linear least
#' squares version is used.  This is highly recommended if there is a suspicion
#' that the data may contain outliers. For this, the package 'Rfit' is used. The
#' method allows hypothesis testing and confidence interval calculation similar
#' to what is poissible with least squares.  See Kloke and McKean for more
#' information.  When using normalization here, the robust intercept rather than
#' mean values are used and reported (see values below). Note that under
#' relatively rare circumstances, the optimization procedure used by Rfit may
#' fail.  In that case, NAs are reported for p values and estimates.  A warning
#' for each such occurrence is given.

#'
#' Values returned include normalized time series (if normalization is
#' performed), normalization weights (means or polynomial coefficients if
#' polynomial normalization is used), fitted normalized curves, fitted
#' non-normalized curves, a data frame of amplitudes and phases (in radians),
#' p-values according to an F-test (Halberg 1967), a data frame of approximately
#' 1.96 standard deviations for the amplitude and phase estimates, a matrix of
#' coefficients a and b and possibly c,... , the sum square resuduals after the
#' fit for each time series, and the covariance matrix for all independent
#' variables, starting with (\eqn{1}, \eqn{cos(\omega t)}, and \eqn{sin(\omega
#' t)}), and followed by possible nuisance variables. The latter can be used in
#' post-processing e.g. to obtain individual p-values for coefficients by
#' t-tests.

#' @param inputts Matrix of time series.  Rows correspond to time points,
#'   columns to samples.  If a vector is provided, it is coerced to a matrix.
#' @param inputtime Vector of the time points corresponding to the row in the
#'   time series matrix.
#' @param Tau Scalar giving the oscillation period to estimate and test for.
#' @param normalize Boolean, set to \code{TRUE} if normalization is to be
#'   performed (default).  Unless \code{norm.pol = TRUE}, normalization is
#'   performed by dividing with the mean.
#' @param norm.pol Boolean, set to \code{TRUE} if a polynomial should be fitted
#'   to each time series, and used for normalization.  In this case, each point
#'   in a time series will be divided by the value of the fitted polynomial.
#'   Defaults to \code{FALSE}.
#' @param norm.pol.degree Scalar indicating the polynomial degree for the
#'   normalization (ignored if \code{norm.pol = FALSE}).
#' @param nuisance_f A formula object, representing any nuisance variables to be
#'   used in addition to the intercept (the estimated overall mean). The user is
#'   responsible for ensuring that included variables are defined. If
#'   \code{normalize = TRUE}, normalization by the intercept is performed. Thus
#'   users should ensure that nuisance variables are centered, or if factors are
#'   used, that their contrasts are set so that the intercept can be interpreted
#'   as an overall mean.
#' @param robust Boolean should a robust semi-parametric method from the package
#'   "Rfit" be used?
#' @param a_over_sigma Numeric, for non-central F test for arrhythmicity, this
#'   is the definition of the threshold amplitude divided by noise standard
#'   deviation below which arrhythmicity is defined.  Note that "arrhythmicity"
#'   here always refers to the absence of a rhythmic component with period Tau:
#'   Rhythms with other periods may well be present.
#' @param n_cores Unused at this stage
#' @param robust_scores Object of class \code{scores}, see the documentation of
#'   package \code{Rfit}, and the textbook (Kloke and McKean, 2014). Defaults to
#'   Wilcoxon scores, which are recommended for unknown error distributions. If
#'   a strongly right (left)-skewed error distribution is suspected,
#'   \code{Rfit::bentscores1} (\code{Rfit::bentscores3}) may be more
#'   appropriate, respectively. For heavy-tailed distributions in both
#'   directions, \code{Rfit::bentscores4} are recommended.
#'   
#' @return A list containing: \tabular{ll}{
#'  \code{means} \tab Vector (if \code{norm.pol = FALSE}) or matrix (otherwise) 
#'  of the means or robust means, or coefficients of the fitted polynomial used
#'  for normalization.  If using trend elimination AND (robust) mean
#'  normalization, then the mean (median) is used and reported.  If using
#'  polynomial normalization, the coefficients are reported, in which case the
#'  user is left to interprete the meaning of the normalization and computed
#'  amplitudes. If using trend elimination without normalization the median
#'  \code{(robust = TRUE)} or mean (otherwise) is reported. \cr
#'  \code{normts} \tab Matrix of mean-scaled or normalized-by-polynomial time 
#'  series, same dimensionality as \code{inputts}. \cr
#'  \code{fit.vals} \tab Matrix of model fitted values to \code{inputts}. \cr
#'  \code{norm.fit.vals} \tab Matrix of model fitted values to the normalized
#'  (trend eliminated or mean scaled) time series. \cr
#'  \code{pars} \tab Data frame of estimated amplitudes and phases (in radians,
#'  between 0 and \eqn{2\pi}). \cr
#'  \code{pvals} \tab Vector of p-values according to an F-test of the model fit
#'  against a restricted model (mean-centering only). \cr
#'  \code{ci} \tab Data frame of one-sided approximative 95\% (\eqn{2\sigma})
#'  confidence intervals for the estimated amplitudes and phases. \cr
#'  \code{coeffs} \tab Matrix of estimated model parameters \eqn{a} and \eqn{b}.
#'  \cr
#'  \code{ssr} \tab Vector of sum square residuals for the model fits. These are
#'  the normalized residuals if \code{normalize = TRUE}, equivalent to the 
#'  squared coefficient of variation. \cr
#'  \code{df} \tab Scalar if \code{inputts} does not contain \code{NA}s and
#'  Vector otherwise, representing the degrees of freedom of the residual from
#'  the fit. \cr
#'  \code{sigma_hat} \tab Scalar estimate of the residual standard deviation.
#'  For information on the robust estimate (when \code{robust = TRUE}), see
#'  Hettmansperger & McKean (2011) page 203. \cr
#'  \code{ssx} \tab Matrix (3 times 3, if \code{inputts} does not contain 
#'  \code{NA}s, a list of such matrices, one for each time series, otherwise) of
#'  covariances for the dependent variables corresponding to (\eqn{m}, \eqn{a
#'  cos(\omega t)}, and \eqn{b sin(\omega t)}, respecively). \cr
#'  \code{pvals_son} \tab Vector of p-values according to a noncentral F-
#'  test for arrhythmicity against a compound null hypothesis of rhythmicity, 
#'  defined by the cutoff \code{a_over_sigma}. If \code{robust = TRUE}, this is 
#'  a noncentral chi-square test, which is slightly more conservative. \cr
#' }
#' 
#' @examples
#' # load example data
#' data(rna_polya)
#' rna_mat <- t(as.matrix(rna_polya[1:200, -1]))
#' head(rownames(rna_mat))
#' 
#' # extract time points
#' timepoints <- as.numeric(
#'   regmatches(names(rna_polya)[-1],
#'   regexpr("\\d+", names(rna_polya)[-1])))
#'   
#' # perform harmonic regression
#' hreg <- HarmonicRegression::harmonic_regression(
#'   inputts = rna_mat, inputtime = timepoints)
#'
#' # inspect results
#' hreg
#' summary(hreg)
#' 
#' # inspect p values
#' summary(hreg$pvals)
#' summary(hreg$pvals_son)

#' @importFrom stats coef fitted lm median na.exclude pf residuals var
#'   vcov anova update deviance pnorm pt pchisq
#' @export

#' @references Halberg F, Tong YL, Johnson EA: Circadian System Phase -- An
#'   Aspect of Temporal Morphology; Procedures and Illustrative Examples. in:
#'   The Cellular Aspects of Biorhythms, Springer 1967. \cr Hettmansperger TP
#'   and McKean, JW: Robust nonparametric statistical methods. CRC Press 2011.
#'   \cr Kloke J and McKean JW: Nonparametric Statistical Methods Using R. CRC
#'   Press 2014. \cr Davison, AC: Statistical Models. Cambridge University Press
#'   2003.


# harmonic regression main function ---------------------------------------
harmonic_regression <- function (inputts, inputtime, Tau = 24,
                                 normalize = TRUE, norm.pol = FALSE, 
                                 norm.pol.degree = 1,
                                 nuisance_f = NULL, robust = FALSE,
                                 a_over_sigma = 1,
                                 n_cores = 1L,
                                 robust_scores = Rfit::wscores) {
  
  ## check that input data come as numerics
  if (!is.numeric(inputts) || !is.numeric(inputtime)) {
    stop("Input data (inputts, inputtime) must be numeric")
  }
  
  ## try to coerce input to matrix, if needed
  if (is.vector(inputts)) {
    inputts <- as.matrix(inputts)
  }
  
  ## check series lengths
  if (nrow(inputts) != length(inputtime)) {
    stop(paste("Length of time series (inputts):", nrow(inputts), " and time",
               "points (inputtime):", length(inputtime), "do not match."))
  }
  
  if (norm.pol.degree < 1) {
    stop(paste("Polynomial for normalization must ",
               "be of degree 1 or more. Aborting."))
  }
  
  if (!is.null(nuisance_f) &&
      attr(stats::terms(nuisance_f), "response") != 0L) {
    stop("nuisance_f must be given and formulated without response variables")
  }
  
  ## can't accept ts objects with NAs
  if (any(is.na(inputts)) && ("ts" %in% class(inputts))) {
    stop(paste("Time series object with NAs was supplied.  This is not yet",
               "supported; please supply data containing NAs as a plain", 
               "matrix.  Aborting."))
  }
  
  if (robust) {
    
    # inputtime <- times
    # inputts <- t(rhythm.data)
    # Tau <- 24
    
    if (normalize) {
      
      if (norm.pol) {
        ## in this special case, user is responsible for sensible interpretation
        ## of the normalization
        # if (n_cores == 1L || !requireNamespace("future.apply", 
        #                                        quietly = TRUE)) {
        norm.ts.list <- apply(inputts, 2, normalize_one_ts_r_pol,
                              inputtime, norm.pol.degree,
                              robust_scores)
        # } else {
        #   future::plan(future::multisession,
        #                workers = min(n_cores, future::availableCores()))
        #   norm.ts.list <- 
        #     future.apply::future_apply(inputts, 2, normalize_one_ts_r_pol,
        #                                inputtime, norm.pol.degree,
        #                                robust_scores)
        # }
        norm.ts <- sapply(norm.ts.list, "[[", "norm.ts")
        norm.w <- sapply(norm.ts.list, "[[", "norm.w")
        norm.vals <- sapply(norm.ts.list, "[[", "norm.vals")
        ## normalization just performed; use normalize = FALSE below
        if (!is.null(nuisance_f)) {
          results <- harmonic_regression_nas_nuisance(norm.ts, inputtime,
                                                      Tau, nuisance_f,
                                                      robust = TRUE,
                                                      normalize = FALSE,
                                                      a_over_sigma = 
                                                        a_over_sigma,
                                                      n_cores = n_cores,
                                                      robust_scores =
                                                        robust_scores)
        } else { 
          results <- harmonic_regression_nas(norm.ts, inputtime, Tau,
                                             robust = TRUE, normalize = FALSE,
                                             a_over_sigma = a_over_sigma,
                                             n_cores = n_cores,
                                             robust_scores = robust_scores)
        }
        results <- append(results, list(norm.fit.vals = results$fit.vals),
                          after = 1)
        ## delete the overall mean returned by the regression, we use the
        ## normalization polynomial coefficients here instead.
        results$means <- NULL
        results <- c(list(means = norm.w), results)
        results <- append(results, list(normts = norm.ts),
                          after = 1)
        results$fit.vals <- results$norm.fit.vals*norm.vals
        
      } else {
        ## in this case, normalization will be performed from within
        ## fit_one_harmonic.*.r for efficiency reasons
        norm.ts <- inputts
        if (!is.null(nuisance_f)) {
          results <- harmonic_regression_nas_nuisance(norm.ts, inputtime,
                                                      Tau, nuisance_f,
                                                      robust = TRUE,
                                                      normalize = TRUE,
                                                      a_over_sigma = 
                                                        a_over_sigma,
                                                      n_cores = n_cores,
                                                      robust_scores =
                                                        robust_scores)
          
        } else {
          results <- harmonic_regression_nas(norm.ts, inputtime, Tau,
                                             robust = TRUE, normalize = TRUE,
                                             a_over_sigma = a_over_sigma,
                                             n_cores = n_cores,
                                             robust_scores = robust_scores)
          ## reset the intercept to 1 for consistency with classic harmonic
          ## regression
          results$coeffs[, 1] <- 1
        }
        
        results <- append(results, list(norm.fit.vals = results$fit.vals),
                          after = 1)
        ## we have to do the real ts normalization here:
        results <- append(results, 
                          list(normts = sweep(norm.ts, 2, results$means, "/")),
                          after = 1)
        results$fit.vals <- sweep(results$norm.fit.vals, 2, results$means, "*")
        
      }
      
      ## end if (normalize)
      
    } else {
      if (!is.null(nuisance_f)) {
        results <- harmonic_regression_nas_nuisance(inputts, inputtime, 
                                                    Tau, nuisance_f,
                                                    robust = TRUE,
                                                    normalize = FALSE,
                                                    a_over_sigma = a_over_sigma,
                                                    n_cores = n_cores,
                                                    robust_scores =
                                                      robust_scores)
      } else {
        results <- harmonic_regression_nas(inputts, inputtime, Tau,
                                           robust = TRUE, normalize = FALSE,
                                           a_over_sigma = a_over_sigma,
                                           n_cores = n_cores,
                                           robust_scores = robust_scores)
      }
    }
    
    ## end if (robust) 
    
  } else {
    
    ## check if there are NAs
    if (any(is.na(inputts))) {
      
      ## if NAs; call the slow versions handling NAs separately for each time
      ## series
      if (normalize) {
        # if (n_cores == 1L || !requireNamespace("future.apply", 
        #                                        quietly = TRUE)) {
          norm.ts.list <- apply(inputts, 2, normalize_one_ts,
                                inputtime, norm.pol, norm.pol.degree)
        # } else {
        #   future::plan(future::multisession,
        #                workers = min(n_cores, future::availableCores()))
        #   norm.ts.list <- 
        #     future.apply::future_apply(inputts, 2, normalize_one_ts,
        #                                inputtime, norm.pol, norm.pol.degree)
        # }
        norm.ts <- sapply(norm.ts.list, "[[", "norm.ts")
        norm.w <- sapply(norm.ts.list, "[[", "norm.w")
        if (!is.null(nuisance_f)) {
          results <- harmonic_regression_nas_nuisance(norm.ts, inputtime,
                                                      Tau, nuisance_f,
                                                      robust = FALSE,
                                                      normalize = FALSE,
                                                      a_over_sigma = 
                                                        a_over_sigma,
                                                      n_cores = n_cores)
        } else {
          results <- harmonic_regression_nas(norm.ts, inputtime, Tau,
                                             robust = FALSE,
                                             normalize = FALSE,
                                             a_over_sigma = a_over_sigma,
                                             n_cores = n_cores)
        }
        ## results$fit.vals are really fits to normalized time series
        ## fit.vals will be recalculated below
        results <- append(results, list(norm.fit.vals = results$fit.vals),
                          after = 1)
        results <- c(list(means = norm.w), results)
        results <- append(results, list(normts = norm.ts),
                          after = 1)
        if (norm.pol) {
          norm.vals <- sapply(norm.ts.list, "[[", "norm.vals")
          results$fit.vals <- results$norm.fit.vals*norm.vals
        } else {
          results$fit.vals <- sweep(results$norm.fit.vals, 2, norm.w, "*")
        }
      } else {
        if (!is.null(nuisance_f)) {
          results <- harmonic_regression_nas_nuisance(inputts, inputtime, 
                                                      Tau, nuisance_f,
                                                      robust = FALSE,
                                                      normalize = FALSE,
                                                      a_over_sigma = 
                                                        a_over_sigma,
                                                      n_cores = n_cores)
        } else {
          results <- harmonic_regression_nas(inputts, inputtime, Tau,
                                             robust = FALSE,
                                             normalize = FALSE,
                                             a_over_sigma = a_over_sigma,
                                             n_cores = n_cores)
        }
        results <- c(list(means = colMeans(inputts, na.rm = TRUE)), results)
      }
      
    } else {
      
      ## use fast vectorized versions
      
      if (normalize) {
        ## normalization
        norm.ts.list <- normalize_ts_matrix(inputts, inputtime, 
                                            norm.pol, norm.pol.degree)
        if (!is.null(nuisance_f)) {
          results <- harmonic_regression_matrix_nuisance(norm.ts.list$norm.ts,
                                                         inputtime, Tau, 
                                                         nuisance_f,
                                                         a_over_sigma = 
                                                           a_over_sigma)
        } else {
          results <- harmonic_regression_matrix(norm.ts.list$norm.ts, 
                                                inputtime, Tau,
                                                a_over_sigma = a_over_sigma)
        }
        results <- append(results, list(norm.fit.vals = results$fit.vals),
                          after = 1)
        results <- c(list(means = norm.ts.list$norm.w), results)
        results <- append(results, list(normts = norm.ts.list$norm.ts),
                          after = 1)
        ## compute non-normalized fitted values
        if (norm.pol) {
          results$fit.vals <- results$norm.fit.vals*norm.ts.list$norm.vals
        } else {
          results$fit.vals <- sweep(results$norm.fit.vals, 2, 
                                    norm.ts.list$norm.w, "*")
        }
        
      } else {
        if (!is.null(nuisance_f)) {
          results <- harmonic_regression_matrix_nuisance(inputts,
                                                         inputtime, Tau, 
                                                         nuisance_f,
                                                         a_over_sigma = 
                                                           a_over_sigma)
        } else {
          results <- harmonic_regression_matrix(inputts, inputtime, Tau,
                                                a_over_sigma = a_over_sigma)
        }
        results <- c(list(means = colMeans(inputts)), results)
      }
      
    }
    
  }
  
  structure(results, class = "hregm", normalized = normalize)
  
}



#' @rdname harmonic_regression
#' @examples
#' 
#' # Old usage (package HarmonicRegression)
#' hreg <- HarmonicRegression::harmonic.regression(
#'   inputts = rna_mat, inputtime = timepoints)
#' @export
harmonic.regression <- harmonic_regression


#' @export
print.hregm <- function (x, fdr = 0.1, amp = 0.15, ...) {
  
  n_ts <- length(x$pvals)
  cat(paste("\n\tHarmonic regression results for", n_ts, "time series.\n"))
  
  if (requireNamespace("qvalue", quietly = TRUE) && 
      length(x$pvals) >= 11L) {
    qobj <- qvalue::qvalue(x$pvals, pi0.method = "bootstrap")
    p_adj <- qobj$qvalues
    cat(paste("\tEstimated proportion rhythmic time series:", 
              paste0(round(1 - qobj$pi0, digits = 2)*100, "%\n")))
  } else {
    p_adj <- stats::p.adjust(x$pvals, method = "BH")
  }
  
  n_fdr <- length(which(p_adj <= fdr))
  n_fdr_amp <- length(which(p_adj <= fdr & x$pars$amp >= amp))
  
  cat(paste0("\n\tAt FDR ", fdr, ":\n\t", round(n_fdr/n_ts, digits = 2)*100,
             "% (", n_fdr, "/", n_ts, ") rhythmic time series\n"))
  
  if (attr(x, "normalized")) {
    cat(paste0("\t", round(n_fdr_amp/n_ts, digits = 2)*100,
               "% (", n_fdr_amp, "/", n_ts, ") rhythmic time series with ",
               "relative amplitude >= ", amp, "\n\n"))
  } else {
    cat("\n\n")
  }
  
  invisible(x)
  
}

#' @export
summary.hregm <- function (object, ...) {
  results <- data.frame(
    means = object$means,
    amplitudes = object$pars$amp,
    amplitudes_ci = object$ci$amp,
    phi = object$pars$phi,
    phi_ci = object$ci$phi,
    sigma_hat = object$sigma_hat,
    pvalues = object$pvals,
    pvalues_son = object$pvals_son
  )
  rownames(results) <- rownames(object$pars)
  results
}

#' @export
coef.hregm <- function (object, ...) {
  object$coeffs
}

#' @export
fitted.hregm <- function (object, normalized = FALSE, ...) {
  if (normalized && length(object) == 13L) {
    object$norm.fit.vals
  } else if (normalized) {
    stop("No time series normalization was performed")
  } else {
    object$fit.vals
  }
}




# deviance.hregm <- function (object, normalized = FALSE) {
#   if (normalized && length(object) == 13L) {
#     object$ssr
#   } else if (normalized) {
#     stop("No time series normalization was performed")
#   } else if (length(object) == 13L) {
#     object$ssr*object$means^2
#   } else {
#     object$ssr
#   }
# }




redistribute_pvals <- function (pvals, breakpoint = 0.9, next_breakpoint = 0.8,
                                ref_low = 0.6, ref_high = 0.7) {
  candidate_ind <- which(pvals > breakpoint)
  ## compute predicted number of p values in the interval 
  ## 1 >= p > breakpoint, based on the interval
  ## ref_high > p > ref_low
  n_unif_reference <-
    round(length(which(pvals > ref_low & pvals < ref_high)) *
            (1 - breakpoint)/(ref_high - ref_low))
  ## stop if number of p values > breakpoint is already lower than the 
  ## predicted number
  if (length(candidate_ind) <= n_unif_reference) {
    stop("redistribute_pvals(): No right tail peak")
  }
  ## otherwhise: redistribute a number of p values > breakpoint; the number
  ## exceeding the predicted number
  to_replace <- sample(candidate_ind,
                       length(candidate_ind) - n_unif_reference,
                       replace = FALSE)
  ## they get redistributed into an interval 1 >= p > next_breakpoint.
  replace(pvals, to_replace, runif(length(to_replace),
                                   min = next_breakpoint, max = 1))
}




#' Regularize p Values for Size Estimation of Non-null (H1) Population
#'
#' This algorithm assumes a bimodal p value distribution with one peak towards p
#' = 0 and one peak towards p = 1. The latter peak violates assumptions of
#' algorithms generating q values, e.g., \code{qvalue::pi0est()} off the qvalue
#' package. The p values can be redistributed so that the p value population
#' behind the peak towards p = 1 is artificially redistributed to a uniform
#' distribution, removing the peak. This algorithm works with a vector of
#' decreasing breakpoints. It redistributes the p values successively by
#' considering breakpoints one at a time. P values larger than the current
#' breakpoint are replaced as uniformly distributed p values on an interval
#' between the next breakpoint in the breakpoints vector and 1. It uses a
#' reference interval to determine the expected number of p values between the
#' breakpoint and 1. Only a number of p values exceeding this reference number
#' are redistributed.
#'
#'
#' @param pvals Vector of p values
#' @param breakpoints Break points for successive redistribution of p values.
#'   Must be decreasing. A final breakpoint of 0 is implicitly assumed for the
#'   final step
#' @param ref_low Lower limit of reference interval
#' @param ref_high Upper limit of reference interbal
#' @param suppress_message Should the note of caution be suppressed? (Default
#'   \code{FALSE})
#'
#' @return Vector of regularized p values. NOTE: These p values must only be
#'   used for estimating the size of the non-null (H1) population, e.g., with
#'   \code{qvalue::pi0est()} in the qvalue package.
#'
#' @importFrom stats runif
#' @importFrom utils head tail
#'
#' @export
#'
#' @examples
#' ## generate example p values from beta distributions, exhibiting a
#' ## bimodal distribution
#' set.seed(123)
#' pvals <- c(stats::rbeta(1000, 0.5, 5), stats::rbeta(300, 5, 0.5))
#' hist(pvals)
#' hist(regularize_pvals(pvals))
regularize_pvals <- function (pvals, breakpoints = c(0.99, 0.98, 0.9, 0.8),
                              ref_low = 0.6, ref_high = 0.7,
                              suppress_message = FALSE) {
  if (!suppress_message) {
    message(paste("Regularized (pseudo) p values must only be used for",
                  "estimating the size of the non-null (H1) population, e.g.,",
                  "with qvalue::pi0est() in the qvalue package.",
                  "They must not be used as p values for any other purposes."))
  }
  breakpoint <- utils::head(breakpoints, 1)
  breakpoints <- utils::tail(breakpoints, -1)
  if (length(breakpoints) == 0L) {
    next_breakpoint = 0
  } else {
    next_breakpoint = utils::head(breakpoints, 1)
  }
  pvals_next <- try(redistribute_pvals(pvals, breakpoint = breakpoint,
                                       next_breakpoint = next_breakpoint,
                                       ref_low = ref_low, 
                                       ref_high = ref_high),
                    silent = TRUE)
  if (inherits(pvals_next, "try-error")) {
    return(pvals)
  } else if (length(breakpoints) == 0L) {
    return(pvals_next)
  } else {
    return(
      regularize_pvals(pvals_next, breakpoints = breakpoints,
                       ref_low = ref_low, ref_high = ref_high,
                       suppress_message = suppress_message)
    )
  }
}


#' Convert absolute amplitude in log2 scale to relative amplitude in linear
#' scale (vectorized)
#'
#' @param log_amp Logarithmic amplitude(s)
#'
#' @return Relative amplitude(s), which are amplitudes divided by means 
#' @export
#'
#' @examples
#' log2_amp_to_relamp(c(1, 2))
log2_amp_to_relamp <- function (log_amp) {
  (2^(2*log_amp) - 1)/(2^(2*log_amp) + 1)
}

# Data sets documentation -------------------------------------------------

#' Menet et al. RNA-Seq Data
#' 
#' Data originating from a molecular biology research group.  Objective was to 
#' study circadian (rhythmic with 24 hr period) transcriptional activity in 
#' mouse liver.  To detect rhythmicity in these data, harmonic regression may be
#' used.
#' 
#' @format a data frame with nascent RNA-seq data
#' @details The nascent-seq data are thought to reflect transcriptional
#'   activities. Raw data was collected from the supplementary material of the
#'   original publication cited below. In that study, samples were collected at
#'   6 different times of day in two biological replicates (fields ZT* in the
#'   data frame). The field 'mgi_symbol' refers to MGI gene names.
#' @source Jerome S Menet, Joseph Rodriguez, Katharine C Abruzzi, Michael 
#'   Rosbash. Nascent-Seq reveals novel features of mouse circadian 
#'   transcriptional regulation. eLife 1:e00011 (2012).
#' @examples 
#' data(rna_nasc)
#' nasc_t <- seq(0, 44, 4)
#' plot(nasc_t, rna_nasc["Arntl", -1], type = "b")
#' @name rna_nasc
NULL
# "rna_nasc"


#' Menet et al. RNA-Seq Data
#' 
#' Data originating from a molecular biology research group.  Objective was to 
#' study circadian (rhythmic with 24 hr period) transcriptional activity in 
#' mouse liver.  To detect rhythmicity in these data, harmonic regression may be
#' used.
#' 
#' @format a data frame with poly(A)+ RNA-seq data

#' @details The poly(A)+-seq data represent mature mRNA abundances.  mRNA
#'   abundance reflects how strongly a gene is expressed. Raw data was collected
#'   from the supplementary material of the original publication cited below. In
#'   that study, samples were collected at 6 different times of day in two
#'   biological replicates (fields ZT* in the data frame). The field
#'   'mgi_symbol' refers to MGI gene names.

#' @source Jerome S Menet, Joseph Rodriguez, Katharine C Abruzzi, Michael 
#'   Rosbash. Nascent-Seq reveals novel features of mouse circadian 
#'   transcriptional regulation. eLife 1:e00011 (2012).
#' @examples 
#' data(rna_polya)
#' polya_t <- seq(2, 46, 4)
#' plot(polya_t, rna_polya["Arntl", -1], type = "b")
#' @name rna_polya
NULL
# "rna_polya"









