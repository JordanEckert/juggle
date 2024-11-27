#' @title Greedy Domination Algorithm
#'
#' @description \code{dominate_greedy_matrix} finds the MDS of a cover. Uses greedy
#' domination algorithm from West (2000).
#'
#' @param A adjacency matrix of the cover.
#'
#' @return a list which includes:
#'  \item{i_dominant_list}{dominant sample indexes.}
#'
#' @author Jordan Eckert
#'
#' @references
#' Manukyan, A., & Ceyhan, E. (2016). Classification of imbalanced data with a
#' geometric digraph family. Journal of Machine Learning Research, 17(1),
#' 6504–6543. https://jmlr.org/papers/volume17/15-604/15-604.pdf
#'
#' @rdname dominate_greedy_matrix
#' @export

dominate_greedy_matrix <- function(A)
{

  if (!is.matrix(A)) {
    A <- matrix(A, nrow = 1)
  }

  S <- NULL
  n <- nrow(A)
  covered <- rep(FALSE, n)
  card_index <- c()

  while (!all(covered)) {
    domination_scores <- apply(A, 1, sum)
    i <- which.max(domination_scores)
    covered[A[i, ]==TRUE] <- TRUE
    S <- c(S, i)
    card_index <- c(card_index, domination_scores[i])
    A[, covered==TRUE] <- FALSE
  }

  return(list(i_dominant_list = S, cardinality = card_index))
}

#' @title Pure and Proper Class Cover Catch Digraph Classifier
#'
#' @description \code{pcccd} fits a Pure and Proper Class Cover Catch
#' Digraph (PCCCD) classification model.
#'
#' @param x feature matrix or dataframe.
#' @param y class factor variable.
#' @param proportion proportion of covered samples. A real number between \eqn{(0,1]}.
#' 1 by default. Smaller numbers results in less dominant samples.
#' @param tau radius hyperparameter. Default is \code{.Machine$double.eps}.
#'
#' @details
#' Multiclass framework for PCCCD. PCCCD determines target class dominant points
#' set \eqn{S} and their circular cover area by determining balls
#' \eqn{B(x^{\text{target}}, r_{i})} with radii r using minimum amount of
#' dominant point which satisfies \eqn{X^{\text{non-target}}\cap \bigcup_{i}
#' B_{i} = \varnothing} (pure) and \eqn{X^{\text{target}}\subset \bigcup_{i}
#' B_{i}} (proper).
#'
#' This guarantees that balls of target class never covers any non-target
#' samples (pure) and balls cover all target samples (proper).
#'
#' For detail, please refer to Priebe et al. (2001), Priebe et al. (2003),
#' and Manukyan and Ceyhan (2016).
#'
#' @return an object of "cccd_classifier" which includes:
#'  \item{i_dominant_list}{dominant sample indexes.}
#'  \item{x_dominant_list}{dominant samples from feature matrix, x}
#'  \item{radii_dominant_list}{Radiuses of the circle for dominant samples}
#'  \item{class_names}{class names}
#'  \item{k_class}{number of classes}
#'  \item{proportions}{proportions each class covered}
#'  \item{tau}{radius hyperparameter}
#'
#' @author Jordan Eckert
#'
#' @importFrom  RANN nn2
#' @importFrom  Rfast Dist
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
#' n <- 1000
#' x1 <- runif(n, 1, 10)
#' x2 <- runif(n, 1, 10)
#' x <- cbind(x1, x2)
#' y <- as.factor(ifelse(3 < x1 & x1 < 7 & 3 < x2 & x2 < 7, "A", "B"))
#'
#' m_pcccd <- pcccd(x = x, y = y, tau = 1)
#'
#' # dataset
#' plot(x, col = y, asp = 1)
#'
#' # dominant samples of first class
#' x_center <- m_pcccd$x_dominant_list[[1]]
#'
#' # radii of balls for first class
#' radii <- m_pcccd$radii_dominant_list[[1]]
#'
#' # balls
#' for (i in 1:nrow(x_center)) {
#' xx <- x_center[i, 1]
#' yy <- x_center[i, 2]
#' r <- radii[i]
#' theta <- seq(0, 2*pi, length.out = 100)
#' xx <- xx + r*cos(theta)
#' yy <- yy + r*sin(theta)
#' lines(xx, yy, type = "l", col = "green")
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
#' m_pcccd <- pcccd(x = x_train, y = y_train, tau = 1)
#' pred <- classify_pcccd(pcccd = m_pcccd, newdata = x_test)
#'
#' # confusion matrix
#' table(y_test, pred)
#'
#' # test accuracy
#' sum(y_test == pred)/nrow(x_test)
#'
#' # Example 2
#' data("iris")
#'
#' index <- sample(1:nrow(iris), nrow(iris) * 0.7)
#'
#' x <- iris[index, 1:4]
#' y <- iris[index, 5]
#' test_data <- iris[-index, 1:4]
#' test_label <- iris[-index, 5]
#'
#' model <- pcccd(x, y, tau = .5)
#' pred2 <- classify_pcccd(pcccd = model, newdata = test_data)
#'
#' sum(test_label == pred2)/nrow(test_data)
#'
#' @rdname pcccd
#' @export

pcccd <- function(x, y, proportion = 1, tau = 1) {

  if (proportion < 0 | proportion > 1) {
    stop("proportion must be in range [0,1]")
  }

  if (!is.matrix(x) & !is.data.frame(x)) {
    stop("x must be a matrix or data.frame")
  }

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (!is.factor(y)) {
    stop("y must be a factor")
  }

  if (tau <= 0 | tau > 1) {
    stop("tau must be in range (0,1]")
  }

  class_names <- levels(y)
  k_class <- length(class_names)

  i_dominant_list <- vector(mode = "list", length = k_class)
  x_dominant_list <- vector(mode = "list", length = k_class)
  radii_dominant_list <- vector(mode = "list", length = k_class)
  cardinality <- vector(mode = "list", length = k_class)
  proportions <- c()

  for (i in 1:k_class) {
    i_main <- which(y == class_names[i])
    x_main <- x[i_main,]
    x_other <- x[-i_main,]

    n_main <- nrow(x_main)

    dist_main2other <- RANN::nn2(data = x_other, query = x_main, k = 1)$nn.dist
    dist_main2main <- Rfast::Dist(x = x_main)

    M <- dist_main2main < c(dist_main2other)
    M <- matrix(as.numeric(M), n_main)

    cover <- rep(0, n_main)

    m_dominant <- dominate_greedy_matrix(M)

    i_dominant_list[[i]] <- m_dominant$i_dominant
    x_dominant_list[[i]] <- x_main[m_dominant$i_dominant,,drop = FALSE]
    radii_dominant_list[[i]] <- (1-tau)*which.max(dist_main2main[m_dominant$i_dominant_list, m_dominant$i_dominant_list]) + tau*dist_main2other[m_dominant$i_dominant,]
    proportions[i] <- m_dominant$cover_proportion
    cardinality[[i]] <- m_dominant$cardinality
  }

  results <- list(
    i_dominant_list = i_dominant_list,
    x_dominant_list = x_dominant_list,
    radii_dominant_list = radii_dominant_list,
    class_names = class_names,
    k_class = k_class,
    proportions = proportions,
    cardinality = cardinality,
    tau = tau
  )
  class(results) <- "pcccd_classifier"
  return(results)
}

#' @title  Pure and Proper Class Cover Catch Digraph Prediction
#'
#' @description \code{classify_pcccd} makes prediction using \code{pcccd_classifier} object.
#'
#' @param pcccd a \code{pcccd_classifier} object
#' @param newdata newdata as matrix or dataframe.
#' @param type "pred" or "prob". Default is "pred". "pred" is class estimations,
#'  "prob" is \eqn{n\times k} matrix of class probabilities.
#' @details
#' Estimations are based on nearest dominant neighbor.
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
#' m_pcccd <- pcccd(x = x_train, y = y_train, tau = 1)
#' pred <- classify_pcccd(pcccd = m_pcccd, newdata = x_test)
#'
#' # confusion matrix
#' table(y_test, pred)
#'
#' # test accuracy
#' sum(y_test == pred)/nrow(x_test)
#'
#' @rdname classify_pcccd
#' @export


classify_pcccd <- function(pcccd, newdata, type = "pred") {
  if (!(type %in% c("pred", "prob"))) {
    stop("type must be 'pred' or 'prob'")
  }

  if (!is.matrix(newdata) & !is.data.frame(newdata)) {
    stop("newdata must be a matrix or data.frame")
  }

  x_dominant_list <- pcccd$x_dominant_list
  radii_dominant_list <- pcccd$radii_dominant_list
  class_names <- pcccd$class_names
  k_class <- pcccd$k_class

  x <- newdata
  n <- nrow(x)

  dist_prop <- matrix(data = NA, nrow = n, ncol = k_class)

  for (i in 1:k_class) {
    dist_x2dom <- as.matrix(proxy::dist(x_dominant_list[[i]], x))
    prop_x2dom <- dist_x2dom/radii_dominant_list[[i]]
    dist_prop[,i] <- Rfast::colMins(prop_x2dom, value = TRUE)
  }
  prob <- t(apply(1/(dist_prop + 1), 1, function(m) m/sum(m)))

  if (type == "prob") {
    colnames(prob) <- class_names
    return(prob)
  }
  if (type == "pred") {
    pred <- factor(class_names[max.col(prob)], levels = class_names, labels = class_names)
    return(pred)
  }
}

#' @title  Modified Pure and Proper Class Cover Catch Digraph Prediction Rule
#'
#' @description \code{modified_classify_pcccd} makes prediction using \code{pcccd_classifier} object.
#' We use \eqn{\rho = (\frac{d(x,z)}{r(x)})^{|N(x)|^e} } as a dissimlarity metric where \eqn{e} acts
#' as a hyperparameter for the effect of cardinality of the neighbors on classification.
#'
#' @param pcccd a \code{pcccd_classifier} object
#' @param newdata newdata as matrix or dataframe.
#' @param type "pred" or "prob". Default is "pred". "pred" is class estimations,
#'  "prob" is \eqn{n\times k} matrix of class probabilities.
#' @param e value between 0 and 1 that determines the effect of the cardinality of the neighbors on classification.
#' When \code{e = 0}, the prediction is same as original PCCCD, when \code{e = 1}, cardinality is fully incorporated.
#' Default is 0.
#'
#' @return a vector of class predictions (if type is "pred") or a \eqn{n\times p}
#' matrix of class probabilities (if type is "prob").
#'
#' @author Jordan Eckert
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
#' m_pcccd <- pcccd(x = x_train, y = y_train)
#' pred <- modified_classify_pcccd(pcccd = m_pcccd, newdata = x_test, e = 1)
#'
#' # confusion matrix
#' table(y_test, pred)
#'
#' # test accuracy
#' sum(y_test == pred)/nrow(x_test)
#'
#' @rdname modified_classify_pcccd
#' @export

modified_classify_pcccd <- function(pcccd, newdata, type = "pred", e = 0) {
  x_dominant_list <- pcccd$x_dominant_list
  radii_dominant_list <- pcccd$radii_dominant_list
  class_names <- pcccd$class_names
  k_class <- pcccd$k_class
  cardinality <- pcccd$cardinality

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
