#' @title Pure and Proper Class Cover Catch Digraph Resampled Classifier
#'
#' @description \code{juggle} fits a resampled Pure and Proper Class Cover Catch
#' Digraph (PCCCD) classification model.
#'
#' @param x feature matrix or dataframe.
#' @param y class factor variable.
#' @param n_model an integer. Number of bagged iterations.
#' @param replace a bool. Should replacement be used in data sampling, if true, uses bootstrap sampling.
#' @param min_proportion a numeric. Minimum proportion of samples to use in each model, used only if replace is FALSE.
#' @param tau a numeric between \eqn{(0,1]}. Tau parameter for PCCCD.
#' @param verbose a bool. Should progress be printed to console.
#'
#' @details
#' Resampling framework for PCCCD.
#'
#' @return an object of "cccd_classifier" which includes:
#'  \item{i_dominant_list}{dominant sample indexes.}
#'  \item{x_dominant_list}{dominant samples from feature matrix, x}
#'  \item{radii_dominant_list}{Radiuses of the circle for dominant samples}
#'  \item{class_names}{class names}
#'  \item{k_class}{number of classes}
#'  \item{proportions}{proportions each class covered}
#'  \item{tau}{tau parameter for P-CCCD}
#'
#' @author Jordan Eckert
#'
#' @importFrom  stats runif
#'
#' @examples
#'
#' set.seed(123)
#' n <- 1000
#' x1 <- runif(n, 1, 10)
#' x2 <- runif(n, 1, 10)
#' x <- cbind(x1, x2)
#' y <- as.factor(ifelse(3 < x1 & x1 < 7 & 3 < x2 & x2 < 7, "A", "B"))
#'
#' # testing the performance
#' i_train <- sample(1:n, round(n*0.8))
#'
#' x_train <- x[i_train,]
#' y_train <- y[i_train]
#'
#' x_test <- x[-i_train,]
#' y_test <- y[-i_train]
#'
#' m_pcccd <- pcccd(x = x_train, y = y_train, tau = 1)
#' pred <- classify_pcccd(pcccd = m_pcccd, newdata = x_test)
#'
#' juggle_pcccd <- juggle(x = x_train, y = y_train, n_model = 30)
#' pred2 <- classify_juggle(object = juggle_pcccd, newdata = x_test)
#'
#' # confusion matrix
#' table(y_test, pred)
#' table(y_test, pred2)
#'
#' # test accuracy
#' sum(y_test == pred)/nrow(x_test)
#' sum(y_test == pred2)/nrow(x_test)
#'
#' @rdname juggle
#' @export

juggle <-
  function(x,
           y,
           n_model = 30,
           replace = TRUE,
           min_proportion = 0.5,
           tau = 1,
           verbose = FALSE) {

    if (!is.matrix(x) & !is.data.frame(x)) {
      stop("x must be a matrix or data.frame")
    }

    if (is.data.frame(x)) {
      x <- as.matrix(x)
    }

    prop_sample = ifelse(replace, 1, min_proportion)

    class_names <- levels(y)
    k_class <- length(class_names)
    n <- nrow(x)
    p <- ncol(x)
    n_var = ncol(x)

    m_list <- list()
    var_list <- list()
    for (i in 1:n_model) {
      classes <- unique(y)
      guaranteed_points_index <- unlist(lapply(classes, function(cls) {
        sample(which(y == cls), size = 2)
      }))

      i_sample <- sample(1:n, size = (n * prop_sample) - req_points, replace = replace)
      i_sample <- c(i_sample, guaranteed_points_index)

      x_boot <- x[i_sample, , drop = FALSE]
      y_boot <- y[i_sample]

      m_list[[i]] <-
        pcccd(
          x = x_boot,
          y = y_boot,
          tau = tau
        )
      var_list[[i]] <- i_var
      if (verbose) {
        cat("\r", round(i/n_model*100, 2), "%")
      }
    }
    results <- list(m_list = m_list,
                    var_list = var_list,
                    k_class = k_class,
                    n_model = n_model,
                    class_names = class_names)
    class(results) <- "resampled_pcccd"
    return(results)
  }

#' @title  Resampled Pure and Proper Class Cover Catch Digraph Prediction
#'
#' @description \code{classify_juggles} makes prediction using
#' \code{resampled_pcccd} object.
#'
#' @param object a \code{resampled_pcccd} object
#' @param newdata newdata as matrix or dataframe.
#' @param type "pred" or "prob". Default is "pred". "pred" is class estimations,
#'  "prob" is \eqn{n\times k} matrix of class probabilities.
#'
#' @return a vector of class predictions (if type is "pred") or a \eqn{n\times p}
#' matrix of class probabilities (if type is "prob").
#'
#' @author Jordan Eckert jpe0018@auburn.edu
#'
#' @importFrom  stats predict
#'
#' @examples
#' n <- 1000
#' x1 <- runif(n, 1, 10)
#' x2 <- runif(n, 1, 10)
#' x <- cbind(x1, x2)
#' y <- as.factor(ifelse(3 < x1 & x1 < 7 & 3 < x2 & x2 < 7, "A", "B"))
#'
#' # testing the performance
#' i_train <- sample(1:n, round(n*0.8))
#'
#' x_train <- x[i_train,]
#' y_train <- y[i_train]
#'
#' x_test <- x[-i_train,]
#' y_test <- y[-i_train]
#'
#' juggle_pcccd <- juggle(x = x_train, y = y_train, n_model = 30)
#' pred2 <- classify_juggle(object = juggle_pcccd, newdata = x_test)
#'
#' # confusion matrix
#' table(y_test, pred2)
#'
#' # test accuracy
#' sum(y_test == pred2)/nrow(x_test)
#'
#' @rdname classify_juggle
#' @export

classify_juggle <- function(object, newdata, type = "pred") {
  m_list <- object$m_list
  var_list <- object$var_list
  k_class <- object$k_class
  n_model <- object$n_model
  class_names <- object$class_names

  x <- newdata
  n <- nrow(x)

  M_votes <- matrix(data = 0, nrow = n, ncol = k_class)

  for (i in 1:n_model) {
    x_selected <- x[, var_list[[i]], drop = FALSE]
    prob <-
      classify_pcccd(pcccd = m_list[[i]],
                     newdata = x_selected,
                     type = "prob")

    for (j in 1:n) {
      M_votes[j, which.max(prob[j,])] <-
        M_votes[j, which.max(prob[j,])] + 1
    }
  }

  prob <- t(apply(M_votes, 1, function(m)
    m / sum(m)))

  if (type == "prob") {
    colnames(prob) <- class_names
    return(prob)
  }
  if (type == "pred") {
    pred <-
      factor(class_names[max.col(prob)], levels = class_names, labels = class_names)
    return(pred)
  }
}
