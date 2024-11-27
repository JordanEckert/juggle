#' @title AdaCover Digraph Classifier
#'
#' @description Ensemble boosting method using FlexiBall classifiers as weak learners.
#' FlexiBalls are built using a similar framework as P-CCCD classifiers, but with
#' stochastic selection of covering balls and a parameter that controls the number
#' of balls used to cover the data. Boosted weights are used in the stochastic selection.
#' Multiclass classification uses SAMME framework.
#'
#' @param x Training data frame or matrix without labels
#' @param y Training labels
#' @param tau P-CCCD radii parameter between 0 and 1. Default is 1.
#' @param test_data Test data frame or matrix without labels
#' @param num_iter Number of iterations for the AdaBoost algorithm. Default is 10.
#' @param num_balls Cardinality of the covering ball set. Default is 1.
#'
#' @return A \code{list} containing the following attributes:
#' \describe{
#'      \item{final_predictions}{The final predictions on the test data}
#'      \item{alphas}{The alpha value at each stage of boosting}
#' }
#' @author Jordan Eckert
#' @examples
#'
#' data("iris")
#'
#' index <- sample(1:nrow(iris), nrow(iris) * 0.7)
#'
#' x <- iris[index, 1:4]
#' y <- iris[index, 5]
#' test_data <- iris[-index, 1:4]
#'
#' model <- acdc(x, y, tau = .5, test_data, num_iter = 5, num_balls = 1)
#'
#' sum(model$final_predictions == iris[-index, 5])/nrow(iris[-index,])
#'
#' @rdname acdc
#' @export

acdc <- function(x,
                 y,
                 tau = 1,
                 test_data,
                 num_iter = 10,
                 num_balls = 1){

  if (!is.matrix(x) & !is.data.frame(x)) {
    stop("x must be a matrix or data.frame")
  }

  if (!is.numeric(num_iter) || !is.numeric(num_balls)) {
    stop("Number of iterations and size of MDS must be numeric")
  }

  if (num_balls > nrow(x)) {
    stop("Cardinality of dominating set must not be bigger than training data size")
  }

  if (tau <= 0 || tau > 1) {
    stop("Tau must be in (0, 1]")
  }

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  # Store classifiers and weights
  classifiers <- list()
  alphas <- list()
  predictions <- data.frame(matrix(ncol = num_iter, nrow = nrow(x)))

  # Initialize weights
  weights <- rep(1/nrow(x), nrow(x))

  # Factor for class labels
  y <- as.factor(y)

  # Number of classes
  class_names <- levels(y)
  k_class <- length(class_names)

  for(i in 1:num_iter){

    classifier <- weighted_pcccd(x, y, weights = weights, num_balls = num_balls, tau = tau)
    classifiers[[i]] <- classifier

    # Predict
    pred <- classify_pcccd(classifier, x, type = "pred")
    predictions[,i] <- pred

    # Calculate error (adds epsilon to avoid division by zero in alpha calc)
    err <- sum(weights * as.numeric(pred != y))/sum(weights) + .Machine$double.eps

    # Calculate alpha
    alpha <- (1/(k_class))*log((1-err)/err) + log(k_class - 1)
    alphas[[i]] <- alpha

    # Update weights
    weights <- weights * exp(alpha * (pred != y))
  }

  # Predict on test data
  test_pred <- data.frame(matrix(ncol = num_iter, nrow = nrow(test_data)))
  weighted_pred <- list()

  for(i in 1:num_iter){
    test_pred <- classify_pcccd(classifiers[[i]], test_data, type = "prob")
    weighted_pred[[i]] <- test_pred * alphas[[i]]
  }

  final_pred <- Reduce(`+`, weighted_pred)
  final_predictions <- apply(final_pred, 1, function(x) colnames(final_pred)[which.max(x)])

  return(list(final_predictions = final_predictions, alphas = alphas))
}

#' @title Stochastic selection for MDS
#'
#' @description Finds the minimum dominating set of a digraph greedily using AdaBoosted weights
#'
#' @param M Adjacency matrix of the digraph
#' @param num_balls Final cardinality of the MDS
#' @param weights Weights of the vertices
#'
#' @return A \code{list} containing the following attributes:
#' \describe{
#'      \item{i_dominant_list}{index of the dominate set vertices}
#' }
#' @author Jordan Eckert
#' @export

dominate_stochastic_matrix <- function(M = M, num_balls,
                                       weights = weights_main){
  if (!is.matrix(M)) {
    M <- matrix(M, nrow = 1)
  }

  S <- NULL
  n <- nrow(M)
  covered <- rep(FALSE, n)

  i <- sample(1:n, num_balls, replace = TRUE, prob = weights)
  S <- c(S, i)

  return(list(i_dominant_list = S))
}

#' @title Weighted Pure and Proper Class Cover Catch Digraph Classifier
#'
#' @description \code{weighted_pcccd} fits a Pure and Proper Class Cover Catch
#' Digraph (PCCCD) classification model using weights of points to greedily determine MDS
#'
#' @param x feature matrix or dataframe.
#' @param y class factor variable.
#' @param weights Weighting vector. Default is equal weighting.
#' Smaller numbers results in less likely to be used in dominating algorithm.
#' @param num_balls Number of minimum dominating balls to be used in the cover
#' @param tau PCCCD radii parameter between 0 and 1. Default is 1.
#'
#' @return an object of "cccd_classifier" which includes:
#'  \item{i_dominant_list}{dominant sample indexes.}
#'  \item{x_dominant_list}{dominant samples from feature matrix, x}
#'  \item{radii_dominant_list}{Radiuses of the circle for dominant samples}
#'  \item{weights}{boosted weights used to be returned}
#'  \item{class_names}{class names}
#'  \item{k_class}{number of classes}
#'  \item{proportions}{defunct parameter}
#'
#' @author Jordan Eckert jordan@eckert.network
#'
#' @importFrom  RANN nn2
#' @importFrom  Rfast Dist
#'
#' @rdname weighted_pcccd
#' @export

weighted_pcccd <- function(x = x, y = y, weights = weights,
                           num_balls, tau){

  # Number of classes
  class_names <- levels(y)
  k_class <- length(class_names)

  # Store model information for each class
  i_dominant_list <- vector(mode = "list", length = k_class)
  x_dominant_list <- vector(mode = "list", length = k_class)
  radii_dominant_list <- vector(mode = "list", length = k_class)

  for (i in 1:k_class) {
    # Separate into target/nontarget
    i_main <- which(y == class_names[i])
    x_main <- x[i_main,]
    x_other <- x[-i_main,]

    n_main <- nrow(x_main)

    # Create adjacency matrix
    dist_main2other <- RANN::nn2(data = x_other, query = x_main, k = 1)$nn.dist
    dist_main2main <- Rfast::Dist(x = x_main)

    # Incorporate tau
    M <- dist_main2main <= c(tau*dist_main2other)
    M <- matrix(as.numeric(M), n_main)

    # Normalize weights_main
    weights_main <- weights[i_main]/sum(weights[i_main])

    # Stochastic Cover Selection
    m_dominant <- dominate_stochastic_matrix(M, num_balls, weights_main)

    # Storing class cover information
    i_dominant_list[[i]] <- m_dominant$i_dominant_list
    x_dominant_list[[i]] <- x_main[m_dominant$i_dominant_list,,drop = FALSE]
    radii_dominant_list[[i]] <- tau * dist_main2other[m_dominant$i_dominant,]
    weights[i_main] <- weights_main

  }

  # Normalize after weights have been used in model
  weights <- weights/sum(weights)

  results <- list(
    i_dominant_list = i_dominant_list,
    x_dominant_list = x_dominant_list,
    radii_dominant_list = radii_dominant_list,
    weights = weights,
    class_names = class_names,
    k_class = k_class,
    proportions = proportions,
    tau = tau)

  class(results) <- "pcccd_classifier"
  return(results)
}
