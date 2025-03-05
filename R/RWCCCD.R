#' @title  Random Walk Class Cover Catch Digraph Classifier
#'
#' @description \code{rwcccd} fits a
#' Random Walk Class Cover Catch Digraph (RWCCCD) classification model.
#'
#' @param x feature matrix or dataframe.
#' @param y class factor variable.
#' @param method "default" or "balanced".
#' @param m penalization parameter. Takes value in \eqn{[0,\infty)}.
#' @param proportion proportion of covered samples. A real number between \eqn{(0,1]}.
#' @param partial_ordering \code{TRUE} or \code{FALSE} Default is \code{FALSE} \code{TRUE} uses partial
#' ordering in determining dominant points. It orders incompletely but faster.
#' Only for \code{rwcccd_classifier_2}.
#'
#' @details
#' Random Walk Class Cover Catch Digraphs (RWCCD) are determined by calculating
#' \eqn{T_{\text{target}}} score for each class as target class as
#'
#' \deqn{
#' T_{\text{target}}=R_{\text{target}}(r_{\text{target}})-\frac{r_{\text{target}}n_u}{2d_m(x)}.
#' }
#'
#' Here, \eqn{r_{\text{target}}} is radius and determined by maximum
#' \eqn{R_{\text{target}}(r) - P_{\text{target}}(r)} calculated for each target sample.
#' \eqn{R_{\text{target}}(r)} is
#' \deqn{
#'   R_{\text{target}}(r):=
#'   w_{target}|{z\in X^{\text{target}}_{n_{\text{target}}}:d(x^{\text{target}},z)\leq r}| -
#'   w_{non-target}|{z\in X^{\text{non-target}}_{n_{\text{non-target}}}:d(x^{\text{target}},z)\leq r}|
#' }
#'
#' and \eqn{P_{\text{target}}(r)} is
#' \deqn{
#'   P_{\text{target}}(r) = m\times d(x^{\text{target}},z)^p.
#' }
#' \eqn{m=0} removes penalty. \eqn{w_{target}=1} for default and
#' \eqn{w_{target}=n_{\text{target}/n_{\text{non-target}}}} for balanced method.
#' \eqn{n_u} is the number of uncovered samples in the current iteration and
#' \eqn{d_m(x)} is \eqn{\max{d(x^{\text{target}},x^{\text{uncovered}})}}.
#'
#' This method is more robust to noise compared to PCCCD However, balls covers
#' classes improperly and \eqn{r = 0} can be selected.
#'
#' For detail, please refer to Priebe et al. (2001), Priebe et al. (2003),
#' and Manukyan and Ceyhan (2016).
#'
#' @return a rwcccd_classifier object
#'  \item{i_dominant_list}{dominant sample indexes.}
#'  \item{x_dominant_list}{dominant samples from feature matrix, x}
#'  \item{radii_dominant_list}{Radiuses of the circle for dominant samples}
#'  \item{class_names}{class names}
#'  \item{k_class}{number of classes}
#'  \item{proportions}{proportions each class covered}
#'
#' @author Jordan Eckert
#'
#' @importFrom  RANN nn2
#' @importFrom  Rfast Dist
#' @importFrom  Rfast dista
#'
#' @references
#' Priebe, C. E., DeVinney, J., & Marchette, D. J. (2001). On the distribution
#' of the domination number for random class cover catch digraphs. Statistics &
#' Probability Letters, 55(3), 239–246. https://doi.org/10.1016/s0167-7152(01)00129-8
#'
#' Priebe, C. E., Marchette, D. J., DeVinney, J., & Socolinsky, D. A. (2003).
#' Classification Using Class Cover Catch Digraphs. Journal of Classification,
#' 20(1), 3–23. https://doi.org/10.1007/s00357-003-0003-7
#'
#' Manukyan, A., & Ceyhan, E. (2016). Classification of imbalanced data with a
#' geometric digraph family. Journal of Machine Learning Research, 17(1),
#' 6504–6543. https://jmlr.org/papers/volume17/15-604/15-604.pdf
#'
#' @examples
#'
#' n <- 500
#' x1 <- runif(n, 1, 10)
#' x2 <- runif(n, 1, 10)
#' x <- cbind(x1, x2)
#' y <- as.factor(ifelse(3 < x1 & x1 < 7 & 3 < x2 & x2 < 7, "A", "B"))
#'
#' # dataset
#' m_rwcccd_1 <- rwcccd(x = x, y = y, method = "default", m = 1)
#'
#' plot(x, col = y, asp = 1, main = "default")
#' # dominant samples of second class
#' x_center <- m_rwcccd_1$x_dominant_list[[2]]
#' # radii of balls for second class
#' radii <- m_rwcccd_1$radii_dominant_list[[2]]
#'
#' # balls
#' for (i in 1:nrow(x_center)) {
#'   xx <- x_center[i, 1]
#'   yy <- x_center[i, 2]
#'   r <- radii[i]
#'   theta <- seq(0, 2*pi, length.out = 100)
#'   xx <- xx + r*cos(theta)
#'   yy <- yy + r*sin(theta)
#'   lines(xx, yy, type = "l", col = "green")
#' }
#'
#' # dataset
#' m_rwcccd_2 <- rwcccd(x = x, y = y, method = "balanced", m = 1, proportion = 0.5)
#'
#' plot(x, col = y, asp = 1, main = "balanced, proportion = 0.5")
#' # dominant samples of second class
#' x_center <- m_rwcccd_2$x_dominant_list[[2]]
#' # radii of balls for second class
#' radii <- m_rwcccd_2$radii_dominant_list[[2]]
#'
#' # balls
#' for (i in 1:nrow(x_center)) {
#'   xx <- x_center[i, 1]
#'   yy <- x_center[i, 2]
#'   r <- radii[i]
#'   theta <- seq(0, 2*pi, length.out = 100)
#'   xx <- xx + r*cos(theta)
#'   yy <- yy + r*sin(theta)
#'   lines(xx, yy, type = "l", col = "green")
#' }
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
#' m_rwcccd <- rwcccd(x = x_train, y = y_train, method = "balanced")
#' pred <- classify_rwcccd(rwcccd = m_rwcccd, newdata = x_test, e = 1)
#'
#' # confusion matrix
#' table(y_test, pred)
#'
#' # accuracy
#' sum(y_test == pred)/nrow(x_test)
#'
#' @rdname rwcccd
#' @export

rwcccd <- function(
    x,
    y,
    method = "default",
    m = 1,
    proportion = 0.99,
    partial_ordering = FALSE) {

  if (proportion < 0 | proportion > 1) {
    stop("proportion must be in range [0,1]")
  }

  if (!is.matrix(x) & !is.data.frame(x)) {
    stop("x must be a matrix or data.frame")
  }

  if (!is.factor(y)) {
    stop("y must be a factor")
  }

  if (!(method %in% c("default", "balanced"))) {
    stop("method must be 'default' or 'balanced'")
  }

  if (m < 0) {
    stop("m must be bigger then zero")
  }

  if (!isTRUE(partial_ordering) & !isFALSE(partial_ordering)) {
    stop("partial_ordering must be TRUE or FALSE")
  }

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  p <- ncol(x)
  n <- nrow(x)
  class_names <- levels(y)
  k_class <- length(class_names)

  i_dominant_list <- vector(mode = "list", length = k_class)
  x_dominant_list <- vector(mode = "list", length = k_class)
  radii_dominant_list <- vector(mode = "list", length = k_class)
  proportions <- c()
  T_score_list <- vector(mode = "list", length = k_class)

  if (partial_ordering) {
    f_order <- function(x) {
      Rfast::Order(x, partial = TRUE)
    }
  } else {
    f_order <- order
  }

  for (i in 1:k_class) {
    i_main  <- which(y == class_names[i])
    x_main  <- x[i_main,]
    x_other <- x[-i_main,]
    n_main  <- nrow(x_main)
    n_other <- nrow(x_other)

    n_main <- nrow(x_main)

    dist_main2other <- Rfast::dista(xnew = x_main, x = x_other)
    dist_main2main  <- Rfast::Dist(x = x_main)
    dist_main2all   <- cbind(dist_main2main, dist_main2other)

    w <- switch(method,
                default = c(rep(1, n_main), rep(-1, n_other)),
                balanced = c(rep(n_other/n_main, n_main), rep(-1, n_other))
    )

    cover_main    <- rep(0, n_main)
    cover_all     <- rep(0, n)
    code_dominant <- rep(0, n_main)
    radii         <- rep(NA, n_main)
    T_score       <- rep(NA, n_main)
    ww <- w
    cover_proportion <- 0

    while (Inf) {
      T_scores <- rep(NA, n_main)
      rs <- rep(NA, n_main)
      i_potentialdom <- which(code_dominant == 0 & cover_main == 0)

      if (length(i_potentialdom) == 0) {
        break
      }

      n_uncovered <- sum(cover_main == 0)

      for (j in i_potentialdom) {
        dd <- dist_main2all[j,]
        ii_sort <- f_order(dd)
        dd_sorted <- dd[ii_sort]
        w_sorted <- ww[ii_sort]
        R_all <- cumsum(w_sorted)
        P_all <- m*dd_sorted^p
        R_penalized_all <- R_all - P_all
        i_R_penalized_max <- which.max(R_penalized_all)
        rs[j] <- dd_sorted[i_R_penalized_max]
        T_scores[j] <- R_all[i_R_penalized_max] - rs[j]/max(dist_main2main[j,]) * (n_uncovered/2)
      }

      i_selected <- which.max(T_scores)
      T_score_selected <- T_scores[i_selected]
      r_selected <- rs[i_selected]
      i_cover_main <- which(dist_main2main[i_selected,] < r_selected)

      if (length(i_cover_main) > 0) {
        cover_main[i_cover_main] <- 1
      }

      i_cover_all <- which(dist_main2all[i_selected,] < r_selected)

      if (length(i_cover_all) > 0) {
        cover_all[i_cover_all] <- 1
        ww[i_cover_all] <- 0
      }

      code_dominant[i_selected] <- 1
      radii[i_selected] <- r_selected
      T_score[i_selected] <- T_score_selected
      cover_proportion <- sum(cover_main)/n_main

      if(is.na(cover_proportion)) {
        cover_proportion <- 0
        break
      }

      if (cover_proportion >= proportion) {
        break
      }

      if (sum(code_dominant) == n_main) {
        break
      }

      if (sum(cover_main) == n_main) {
        break
      }
    }

    i_dominant_list[[i]] <- which(code_dominant == 1)
    radii_dominant_list[[i]] <- radii[i_dominant_list[[i]]]
    x_dominant_list[[i]] <- x_main[i_dominant_list[[i]],,drop = FALSE]
    T_score_list[[i]] <- T_score[i_dominant_list[[i]]]
    proportions[i] <- cover_proportion
  }

  results <- list(
    i_dominant_list = i_dominant_list,
    x_dominant_list = x_dominant_list,
    radii_dominant_list = radii_dominant_list,
    T_score_list = T_score_list,
    class_names = class_names,
    k_class = k_class,
    proportions = proportions
  )
  class(results) <- "rwcccd_classifier"
  return(results)
}

#' @title  Random Walk Class Cover Catch Digraph Prediction
#'
#' @description \code{classify_rwcccd} makes prediction using
#' \code{rwcccd_classifier} object.
#'
#' @param rwcccd a \code{rwcccd_classifier} object
#' @param newdata newdata as matrix or dataframe.
#' @param type "pred" or "prob". Default is "pred". "pred" is class estimations,
#'  "prob" is \eqn{n\times k} matrix of class probabilities.
#' @param e 0 or 1. Default is 0. Penalty based on \eqn{T} scores in
#' \code{rwcccd_classifier} object.
#'
#' @details
#' Estimations are based on nearest dominant neighbor.
#' \code{e} argument is used to penalize estimations based on \eqn{T} scores in
#' \code{rwcccd_classifier} object.
#'
#' For detail, please refer to Priebe et al. (2001), Priebe et al. (2003),
#' and Manukyan and Ceyhan (2016).
#'
#' @return a vector of class predictions (if type is "pred") or a \eqn{n\times p}
#' matrix of class probabilities (if type is "prob").
#'
#' @author Jordan Eckert
#'
#' @references
#' Priebe, C. E., DeVinney, J., & Marchette, D. J. (2001). On the distribution
#' of the domination number for random class cover catch digraphs. Statistics &
#' Probability Letters, 55(3), 239–246. https://doi.org/10.1016/s0167-7152(01)00129-8
#'
#' Priebe, C. E., Marchette, D. J., DeVinney, J., & Socolinsky, D. A. (2003).
#' Classification Using Class Cover Catch Digraphs. Journal of Classification,
#' 20(1), 3–23. https://doi.org/10.1007/s00357-003-0003-7
#'
#' Manukyan, A., & Ceyhan, E. (2016). Classification of imbalanced data with a
#' geometric digraph family. Journal of Machine Learning Research, 17(1),
#' 6504–6543. https://jmlr.org/papers/volume17/15-604/15-604.pdf
#'
#' @importFrom proxy dist
#' @importFrom Rfast colMins
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
#' m_rwcccd <- rwcccd(x = x_train, y = y_train)
#' pred <- classify_rwcccd(rwcccd = m_rwcccd, newdata = x_test, e = 1)
#'
#' # confusion matrix
#' table(y_test, pred)
#'
#' # test accuracy
#' sum(y_test == pred)/nrow(x_test)
#'
#' @rdname classify_rwcccd
#' @export

classify_rwcccd <- function(rwcccd, newdata, type = "pred", e = 0) {

  x_dominant_list <- rwcccd$x_dominant_list
  radii_dominant_list <- rwcccd$radii_dominant_list
  T_score_list <- rwcccd$T_score_list
  class_names <- rwcccd$class_names
  k_class <- rwcccd$k_class

  x <- newdata

  if (e < 0 | e > 1) {
    stop("e must be [0,1]")
  }

  if (!(type %in% c("pred", "prob"))) {
    stop("method must be 'pred' or 'prob'")
  }

  if (!is.matrix(x) & !is.data.frame(x)) {
    stop("newdata must be a matrix or data.frame")
  }

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  n <- nrow(x)
  dist_prop <- matrix(data = NA, nrow = n, ncol = k_class)

  for (i in 1:k_class) {
    dist_x2dom <- as.matrix(proxy::dist(x_dominant_list[[i]], x))
    prop_x2dom <- (dist_x2dom/radii_dominant_list[[i]])^(T_score_list[[i]]^e)
    dist_prop[,i] <- Rfast::colMins(prop_x2dom, value = TRUE)
  }

  prob <- 1 - t(apply(dist_prop, 1, function(m) m/sum(m)))

  if (type == "prob") {
    colnames(prob) <- class_names
    return(prob)
  }

  if (type == "pred") {
    pred <- factor(class_names[max.col(prob)], levels = class_names, labels = class_names)
    return(pred)
  }
}
