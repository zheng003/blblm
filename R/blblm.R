#' @import purrr
#' @import furrr
#' @import future
#' @importFrom vroom vroom_lines
#' @importFrom readr read_csv
#' @importFrom utils capture.output
#' @import stats
#' @importFrom magrittr %>%
#' @aliases NULL
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' Linear Regression with Little Bag of Bootstraps
#'
#'
#' @param formula formula
#' @param file dataframe or a list of file of datasets
#' @param m number of parts data split into
#' @param B number of bootstraps
#' @param cl number of clusters to use
#'
#' @return list
#' @export
blblm <- function(formula, file, m = 10, B = 5000, cl = 1) {
  if (class(file)=="data.frame") {
    data_list <- split_data(file, m)
    n = nrow(file)
  }
  else {
    data_list <- map(file, read_csv,col_types= "dd")
    n = length(vroom_lines(file))-length(file)
  }
  if (cl !=1) {suppressWarnings(plan(multiprocess, workers = cl))}
  estimates <- map(
    data_list,
    ~ lm_each_subsample(formula = formula, data = ., n, B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}

#' Generalized linear models with Little Bag of Bootstraps
#'
#'
#' @param formula formula
#' @param family a description of the error distribution and link function to be used in the model.
#' @param file dataframe or a list of file of datasets
#' @param m number of parts data split into
#' @param B number of bootstraps
#' @param cl number of clusters to use
#'
#' @return list
#' @export
blbglm <- function(formula, family = gaussian, file, m = 10, B = 5000, cl = 1) {
  if (class(file)=="data.frame") {
    data_list <- split_data(file, m)
    n = nrow(file)
  }
  else {
    data_list <- map(file, read_csv,col_types= "dd")
    n = length(vroom_lines(file))-length(file)
  }
  if (cl !=1) {suppressWarnings(plan(multiprocess, workers = cl))}
  estimates <- map(
    data_list,
    ~ glm_each_subsample(formula = formula, family = family, data = ., n, B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blbglm"
  invisible(res)
}



#' split data into m parts of approximated equal sizes
#' @param data datafrmae
#' @param m number of parts data split into
#'
#' @return list
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}



#' compute the estimates
#' @param formula fomula: V1~V2
#' @param data dataframe
#' @param n total number of rows in data
#' @param B number of bootstraps
#'
#' @return list
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}

#' compute the estimates
#' @param formula fomula: V1~V2
#' @param data dataframe
#' @param family a description of the error distribution and link function to be used in the model.
#' @param n total number of rows in data
#' @param B number of bootstraps
#'
#' @return list
glm_each_subsample <- function(formula, family, data, n, B) {
  replicate(B, glm_each_boot(formula, family, data, n), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset
#' @param formula fomula: V1~V2
#' @param data dataframe
#' @param n total number of rows in data
#'
#' @return model
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}

#' compute the regression estimates for a blb dataset
#' @param formula fomula: V1~V2
#' @param family a description of the error distribution and link function to be used in the model.
#' @param data dataframe
#' @param n total number of rows in data
#'
#' @return model
glm_each_boot <- function(formula, family, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  glm1(formula, family, data, freqs)
}

#' estimate the regression estimates based on given the number of repetitions
#' @param formula fomula: V1~V2
#' @param data dataframe
#' @param freqs multinomially distributed random number vectors
#'
#' @return list
lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}

#' estimate the regression estimates based on given the number of repetitions
#' @param formula fomula: V1~V2
#' @param family a description of the error distribution and link function to be used in the model.
#' @param data dataframe
#' @param freqs multinomially distributed random number vectors
#'
#' @return list
glm1 <- function(formula, family, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- glm(formula, family, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' compute the coefficients from fit
#' @param fit lm model
#'
#' @return list
blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit
#' @param fit model
#'
#' @return double
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}



#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}

#' @export
#' @method print blbglm
print.blbglm <- function(x, ...) {
  cat("blbglm model:", capture.output(x$formula))
  cat("\n")
}


#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' @export
#' @method sigma blbglm
sigma.blbglm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#' @export
#' @method coef blbglm
coef.blbglm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}


#' @export
#' @method confint blbglm
confint.blbglm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms.formula(object$formula, data = new_data), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


#' @export
#' @method predict blbglm
predict.blbglm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms.formula(object$formula, data = new_data), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
               apply(1, mean_lwr_upr, level = level) %>%
               t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
