#' @title Creates region of proposed overlap
#'
#' @description \code{find_hyperplane_boundary_points} uses SVM boundary and decision values to find the points between two parallel hyperplanes.
#'
#' @param svm_model SVM model with decision values
#' @param x_pair training data of the proposed overlap region
#' @param y_pair labels of the proposed overlap region
#' @param delta_1 minimum distance for the first hyperplane
#' @param delta_2 maximum distance for the second hyperplane
#'
#' @details
#' Fits a support vector machine using a Gaussian kernel to find the discrimination hyperplane.
#' Two parallel hyperplanes \eqn{H_1} and {H_2} are fit such that \eqn{d(H_1, H_2) = \delta_1 + \delta_2}.
#' Points inside region are considered in overlap region, if such region exists.
#'
#' @return a list which includes:
#' \item{indices}{index of points in overlap region}
#' \item{points}{coordinates of points in overlap region}
#' \item{labels}{labels of points in overlap region}
#'
#' @author Jordan Eckert
#'
#' @importFrom  stats runif
#' @import e1071
#'
#' @rdname find_hyperplane_boundary_points
#' @export

# Function to find points between parallel hyperplanes
find_hyperplane_boundary_points <- function(svm_model, x_pair, y_pair, delta_1, delta_2) {
  # Get decision values for all points
  decision_values <- predict(svm_model, as.matrix(x_pair), decision.values = TRUE)
  decision_values_abs <- abs(attr(decision_values, "decision.values"))

  # Find points between delta_1 and delta_2
  boundary_points_indices <- which(
    decision_values_abs >= min(delta_1, delta_2) &
      decision_values_abs <= max(delta_1, delta_2)
  )

  # Return a list with indices and the corresponding points
  return(list(
    indices = boundary_points_indices,
    points = x_pair[boundary_points_indices, ],
    labels = y_pair[boundary_points_indices]
  ))
}

#' @title Aggregates final decision rule
#'
#' @description \code{aggregate_predictions} finds the most frequent class in each row
#'
#' @param matrix_data matrix of predictions
#'
#' @details
#' Maximum is determined by which.max in event of tie
#'
#' @return an list which includes:
#'  \item{predictions}{aggregated final predictions}
#'
#' @author Jordan Eckert
#'
#' @importFrom  stats runif
#'
#' @rdname aggregate_predictions
#' @export

# Aggregate predictions to give final prediction
aggregate_predictions <- function(matrix_data) {
  # Get the number of rows
  n_rows <- nrow(matrix_data)

  # Initialize a vector to store the final predictions for each row
  predictions <- character(n_rows)

  # Loop through each row
  for (i in 1:n_rows) {
    # Get the current row
    row <- matrix_data[i, ]

    # Calculate the frequency of each factor in the row
    factor_counts <- table(row)

    # Find the most frequent factor
    # `which.max` gives the index of the maximum count
    most_frequent_factor <- names(factor_counts)[which.max(factor_counts)]

    # Store the prediction for this row
    predictions[i] <- most_frequent_factor
  }

  return(predictions)
}

#' @title Support Vector Hybrid Class Cover Catch Digraph Classifier
#'
#' @description \code{svcccd} fits a support vector hybrid class cover catch digraph
#' classifier.
#'
#' @param x feature matrix or dataframe.
#' @param y class factor variable.
#' @param test_data Test data frame or matrix without labels
#' @param gamma Parameter to control Gaussian SVM kernel
#' @param e Parameter to RW-CCCD
#' @param cost boolean. If TRUE, uses cost weighted SVM to determine boundary.
#'
#' @details
#' Fits a support vector machine using a Gaussian kernel to find a discrimination hyperplane.
#' Two parallel hyperplanes \eqn{H_1} and \eqn{H_2} are fit such that \eqn{d(H_1, H_2) = \delta}.
#' Points inside region are classified using RW-CCCD, if such region exists. If not, classification
#' proceeds using SVM.
#'
#' @return an object of "svcccd" which includes:
#'  \item{final_predictions}{final predictions}
#'
#' @author Jordan Eckert
#'
#' @importFrom  stats runif
#' @import e1071
#'
#' @examples
#' library(e1071)
#' data("iris")
#' set.seed(1234)
#' index <- sample(1:nrow(iris), nrow(iris) * 0.7)
#'
#' x <- iris[index, 1:4]
#' y <- iris[index, 5]
#'
#' test_data <- iris[-index, 1:4]
#' test_label <- iris[-index, 5]
#'
#' gamma <- .01
#' e <- 1
#'
#' model <- svcccd(x, y, test_data, gamma, e)
#' sum(model$final_predictions == test_label) / length(test_label)
#' @rdname svcccd
#' @export

svcccd <- function(x, y, test_data, gamma = .01, e = 1, cost = TRUE){

  if (!is.matrix(x) & !is.data.frame(x)) {
    stop("x must be a matrix or data.frame")
  }

  if (!is.numeric(gamma)) {
    stop("gamma must be a numeric")
  }

  if (!is.factor(y)) {
    stop("y must be a factor")
  }

  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  # Initialize variables
  class_names <- levels(y)
  k_class <- length(class_names)
  pairwise_combinations <- combn(class_names, 2, simplify = FALSE)
  predictions <- as.data.frame(matrix(NA, nrow = nrow(test_data),
                                      ncol = length(pairwise_combinations)))

  for(i in seq_along(pairwise_combinations)){

    # Subset data based on current class pair
    class_pair <- pairwise_combinations[[i]]
    x_pair <- x[y %in% class_pair,]
    y_pair <- droplevels(y[y %in% class_pair])

    if(cost == TRUE){
      wt.t <- table(y_pair)

      inv.wt.t <- 1/wt.t
      inv.wt.t <- inv.wt.t/min(inv.wt.t)

      svm_model <- svm(x = x_pair, y = y_pair,
                       kernel = "radial",
                       gamma = gamma,
                       cost = inv.wt.t,
                       class.weights = wt.t,
                       decision.values = TRUE)
    } else {
      svm_model <- svm(x = x_pair, y = y_pair,
                       kernel = "radial",
                       gamma = gamma,
                       decision.values = TRUE)
    }

    # See where SVM misclassifies on training data
    svm_labels <- predict(svm_model, as.matrix(x_pair))
    mislabels <- which(y_pair != svm_labels)

    # Edge case of no mislabels, just use SVM predictions
    if(length(mislabels) == 0){
      predictions[,i] <- predict(svm_model, as.matrix(test_data))

      next
    }

    # Constructing delta
    distances <- svm_model$decision.values[mislabels]
    delta_1 <- min(distances)
    delta_2 <- max(distances)

    # Add buffer if degenerate
    if (abs(delta_2 - delta_1) < 1e-6) {
      delta_1 <- delta_1 - 1e-3
      delta_2 <- delta_2 + 1e-3
    }

    # Find points between parallel hyperplanes
    boundary_points <- find_hyperplane_boundary_points(svm_model,
                                                       x_pair, y_pair,
                                                       delta_1, delta_2)

    # If boundary points are empty, just use SVM
    if(length(boundary_points$indices) == 0 || length(levels(droplevels(boundary_points$labels))) == 1){
      predictions[,i] <- predict(svm_model, as.matrix(test_data))

      next
    }

    # If only one unique observation per class exists, use SVM
    if(any(table(boundary_points$indices)) == 1){
      predictions[,i] <- predict(svm_model, as.matrix(test_data))

      next
    }

    # Fit RW-CCCD on boundary points
    rw_model <- rwcccd(x = as.matrix(boundary_points$points),
                       y = boundary_points$labels)

    # Location of test_data in the feature space
    test_data_decision_values <- predict(svm_model,
                                         as.matrix(test_data),
                                         decision.values = TRUE)
    predictions[,i] <- test_data_decision_values

    # Find points between parallel hyperplanes
    boundary_points_test <- find_hyperplane_boundary_points(svm_model, test_data, test_data_decision_values, delta_1, delta_2)

    # If no points are in boundary_points_test, use SVM
    if(length(boundary_points_test$indices) == 0){
      predictions[,i] <- predict(svm_model, as.matrix(test_data))

      next
    }

    # Predict RW-CCCD on boundary_points_test
    rw_labels <- classify_rwcccd(rw_model, boundary_points_test$points, e = e)
    predictions[boundary_points_test$indices, i] <- rw_labels
  }


  final_predictions <- as.factor(aggregate_predictions(as.matrix(predictions)))

  return(list(final_predictions = final_predictions))
}

