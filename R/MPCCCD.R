#' @title  Pure and Proper Class Cover Catch Digraph Prediction
#'
#' @description \code{classify_mpcccd} makes prediction using \code{mpcccd} object. The decision rule uses a modified scoring method based on the cardinality of the dominating ball. The hyperparameter e acts as a weight for the cardinality of the dominating ball in the scoring method.
#'
#' @param object a \code{pcccd_classifier} object
#' @param newdata newdata as matrix or dataframe.
#' @param type "pred" or "prob". Default is "pred". "pred" is class estimations,
#'  "prob" is \eqn{n\times k} matrix of class probabilities.
#' @param e value between 0 and 1 that determines the effect of the cardinality of the neighbors on classification.
#' When \code{e = 0}, the prediction is same as original PCCCD, when \code{e = 1}, cardinality is fully incorporated.
#' Default is 01
#'
#' @details
#' Estimations are based on nearest dominant neighbor in radius unit.
#'
#' For detail, please refer to Priebe et al. (2001), Priebe et al. (2003),
#' and Manukyan and Ceyhan (2016).
#'
#' @return a vector of class predictions (if type is "pred") or a \eqn{n\times p}
#' matrix of class probabilities (if type is "prob").
#'
#' @author Fatih Saglam, saglamf89@gmail.com
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
#' m_pcccd <- mpcccd(x = x_train, y = y_train)
#' pred <- classify_mpcccd(object = m_pcccd, newdata = x_test, e = .7)
#'
#' # confusion matrix
#' table(y_test, pred)
#'
#' # test accuracy
#' sum(y_test == pred)/nrow(x_test)
#'
#' @rdname classify_mpcccd
#' @export

classify_mpcccd <- function(object, newdata, type = "pred", e = 1) {
  x_dominant_list <- object$x_dominant_list
  radii_dominant_list <- object$radii_dominant_list
  class_names <- object$class_names
  k_class <- object$k_class
  cardinality <- object$cardinality

  x <- newdata
  n <- nrow(x)

  dist_prop <- matrix(data = NA, nrow = n, ncol = k_class)

  for (i in 1:k_class) {
    dist_x2dom <- as.matrix(proxy::dist(x_dominant_list[[i]], x))
    prop_x2dom <- (dist_x2dom/radii_dominant_list[[i]])^(cardinality[[i]]^e)
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

#' @title Pure-CCCD Classifier
#'
#' @description Fits a Pure Class Cover Catch Digraph (P-CCCD) classification model. Based on \code{rcccd} package function.
#'
#' @param x Either a matrix or data frame containing the training set.
#' @param y A factor containing the class labels for the training set. Note:
#' must be the same length as the number of rows in \code{x}.
#' @param tau A real number between \eqn{(0,1]} controlling radii of the ball. Choice does not affect the digraph,
#' but can affect classification performance by constructing bigger balls. Default is 1.
#'
#' @details
#' Framework for PCCCD. PCCCD determines target class dominant points
#' set \eqn{S} and their circular cover area by determining balls
#' \eqn{B(x^{\text{target}}, r_{i})} with radii using minimum amount of
#' dominant point which satisfies \eqn{X^{\text{non-target}}\cap \bigcup_{i}
#' B_{i} = \varnothing} (pure) and \eqn{X^{\text{target}}\subset \bigcup_{i}
#' B_{i}} (proper). The calculation of the radii r is based off
#'
#' \eqn{r(x) := (1- \tau)d(x, l(x)) + \tau d(x, u(x))} where
#' \eqn{u(x)} is the argmin distance between x and the non-target points and
#' \eqn{l(x)} is the argmax distance between x and the other target points such that their distance is less than \eqn{u(x)}.
#'
#' This guarantees that balls of target class never covers any non-target
#' samples (pure) and balls cover all target samples (proper).
#'
#' Multiclass classification is handled using a one-versus-rest approach where the
#' non-target class is taken to be the union of all classes that are not the target class.
#'
#' For further detail, please refer to Priebe et al. (2001), Priebe et al. (2003),
#' and Manukyan and Ceyhan (2016).
#'
#' Much of the framework is based off the \code{rcccd} package with minor changes.
#' Used for Jordan Eckert's dissertation research.
#'
#' @return an object of "mpcccd" which includes:
#'  \item{i_dominant_list}{dominant sample indexes.}
#'  \item{x_dominant_list}{dominant samples from feature matrix, x}
#'  \item{radii_dominant_list}{Radiuses of the circle for dominant samples}
#'  \item{cardinality}{cardinality of each neighborhood for dominant samples}
#'  \item{class_names}{class names}
#'  \item{k_class}{number of classes}
#'
#' @author Jordan Eckert
#'
#' @import igraph
#' @importFrom cccd dominate
#'
#' @rdname mpcccd
#' @export

mpcccd <- function(x, y, tau = 1) {
  if (!is.matrix(x) & !is.data.frame(x)) {
    stop("x must be a matrix or data.frame")
  }

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if(tau <= 0 || tau > 1) {
    stop("Tau must be in range (0,1]")
  }

  y <- as.factor(y)
  class_names <- levels(y)
  k_class <- length(class_names)

  i_dominant_list <- vector(mode = "list", length = k_class)
  x_dominant_list <- vector(mode = "list", length = k_class)
  radii_dominant_list <- vector(mode = "list", length = k_class)
  cardinality <- vector(mode = "list", length = k_class)

  for (i in 1:k_class) {
    i_main <- which(y == class_names[i])
    x_main <- x[i_main,]
    x_other <- x[-i_main,]

    n_main <- nrow(x_main)
    dataf <- rbind(x_main, x_other)
    ind <- 1:n_main

    ddata <- distance(dataf, 2)

    R <- apply(ddata[-ind, ind], 2, min)
    M <- matrix(as.integer(ddata[ind, ind] <= R), length(R))
    in.dxx <- ddata[ind,ind]*M
    Rnew <- apply(in.dxx , 1, max)
    R <- Rnew*(1-tau)+R*tau
    M <- matrix(as.integer(ddata[ind,ind] <=R), length(R))
    diag(M) <- 0

    g <- graph_from_adjacency_matrix(M, mode = "directed", diag = FALSE)
    m_dominant <- dominate(g)

    temp <- degree(g, mode = "out")

    i_dominant_list[[i]] <- m_dominant
    x_dominant_list[[i]] <- x_main[m_dominant,,drop = FALSE]
    cardinality[[i]] <- ifelse(temp[m_dominant] > 0, temp[m_dominant] + 1, temp[m_dominant])
    radii_dominant_list[[i]] <- R[m_dominant]
  }

  results <- list(
    i_dominant_list = i_dominant_list,
    x_dominant_list = x_dominant_list,
    radii_dominant_list = radii_dominant_list,
    cardinality = cardinality,
    class_names = class_names,
    k_class = k_class)

  class(results) <- "mpcccd"
  return(results)
}


#' @title Distance function used
#'
#' @import intervals
#' @import geometry
#'
#' @rdname mpcccd
#' @export

distance <- function(data,method=1)
{
  m=method
  ds <- function(x,y,m)
  {
    if(m!=Inf)
    {
      result <- apply(y,1,function(z){
        temp <- sum(abs(z-x)^m)^(1/m)
        return(temp)
      })
    }
    else
    {
      result <- apply(y,1,function(z){
        temp <- max(abs(z-x))
        return(temp)
      })
    }
  }

  result <- apply(data,1,function(x){
    return(ds(x,data,m))
  })

  return(result)
}



