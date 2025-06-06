#' @title Generate datasets from Unit Square
#'
#' @description \code{generate_unit_square} generates a data frame with two classes
#' from a unit square with varying overlap, that is \eqn{X_0 \sim U(0,1)^d} and
#' \eqn{X_1 \sim U(\Delta, 1 + \Delta)^d}.
#'
#' @param n a numeric value. The size of class \eqn{X_0}.
#' @param m a numeric value. The size of class \eqn{X_1}.
#' @param d a numeric value. The dimension of the data.
#' @param Delta a numeric value. The overlap parameter between the two classes.
#'
#' @return a data frame with two classes.
#'
#' @author Jordan Eckert
#'
#' @examples
#' n <- 50
#' m <- 25
#' d <- 2
#' Delta <- .5
#'
#' x <- generate_unit_square(n, m, d, Delta)
#'
#' plot(x[,c(1:2)], col = x$Class, asp = 1)
#'
#' @rdname generate_unit_square
#' @export

generate_unit_square <- function(n, m, d, Delta) {
  # Simulate data from first class
  class1_data <- matrix(runif(n * d), nrow = n)

  # Simulate data from second class
  class2_data <- matrix(runif(m * d, min = 0 + Delta, max = 1 + Delta), nrow = m)

  # Combine data from both classes
  data <- rbind(class1_data, class2_data)

  # Create numeric labels
  labels_numeric <- c(rep(1, n), rep(2, m))

  # Convert numeric labels to factor
  labels_factor <- factor(labels_numeric, levels = c(1, 2), labels = c("Class 1", "Class 2"))

  # Combine data and labels into a data frame
  df <- data.frame(data)
  df$Class <- labels_factor

  # Return data frame
  return(df)
}

#' @title Generate datasets from disjointed Unit Square
#'
#' @description \code{generate_disjoint_square} generates a data frame with two classes
#' from a unit square with varying overlap, but disjointed along the first feature,
#' that is \eqn{X_0 \sim U(0,1)^{d}} and \eqn{X_1 \sim U(1,2) \times U(\Delta, 1 + \Delta)^{d-1}}.
#'
#' @param n a numeric value. The size of class \eqn{X_0}.
#' @param m a numeric value. The size of class \eqn{X_1}.
#' @param d a numeric value. The dimension of the data.
#' @param Delta a numeric value. The overlap parameter between the two classes.
#'
#' @return a data frame with two classes.
#'
#' @author Jordan Eckert
#'
#' @examples
#' library(plotly)
#' n <- 100
#' m <- 25
#' d <- 3
#' Delta <- 0
#'
#' data <- generate_disjoint_square(n, m, d, Delta)
#' plot_ly(x = data$X1, y=data$X2, z=data$X3, type="scatter3d", mode="markers", color=data$Class)
#'
#' @rdname generate_disjoint_square
#' @export

generate_disjoint_square <- function(n, m, d, Delta){
  # Simulate data from first class
  class1_data <- matrix(runif(n * d), nrow = n)

  # Simulate data from second class
  class2_data1 <- matrix(runif(m, min = 1, max = 2))
  class2_data2 <- matrix(runif(m * (d-1), min = 0 + Delta, max = 1 + Delta), nrow = m)
  class2_data <- cbind(class2_data1, class2_data2)

  # Combine data from both classes
  data <- rbind(class1_data, class2_data)

  # Create numeric labels
  labels_numeric <- c(rep(1, n), rep(2, m))

  # Convert numeric labels to factor
  labels_factor <- factor(labels_numeric, levels = c(1, 2), labels = c("Class 1", "Class 2"))

  # Combine data and labels into a data frame
  df <- data.frame(data)
  df$Class <- labels_factor

  # Return data frame
  return(df)
}

#' @title Generate datasets from Normal distributions
#'
#' @description \code{generate_normal_data} generates a data frame with two classes
#' from Normal distributions for each feature with varying overlap,
#' that is \eqn{X_0 \sim N(0,1)^{d}} and \eqn{X_1 \sim N(\Delta, 1)^{d}} where
#' \eqn{\Delta} is a random direction some specified Mahalanobis distance away.
#'
#' @param n a numeric value. The size of class \eqn{X_0}.
#' @param m a numeric value. The size of class \eqn{X_1}.
#' @param d a numeric value. The dimension of the data.
#' @param dist a numeric value. The Mahalanobis distance between the two classes.
#'
#' @return a data frame with two classes.
#'
#' @author Jordan Eckert
#'
#' @examples
#' library(plotly)
#' n <- 100
#' m <- 25
#' d <- 3
#' dist <- 1
#'
#' data <- generate_normal_data(n, m, d, dist)
#' plot_ly(x = data$X1, y=data$X2, z=data$X3, type="scatter3d", mode="markers", color=data$Class)
#'
#' @rdname generate_normal_data
#' @export

generate_normal_data <- function(n, m, d, dist) {
  # Simulate data for the first class (centered at origin)
  class1_data <- matrix(rnorm(n * d, mean = 0, sd = 1), nrow = n)

  # Generate a random direction for the second class
  direction <- rnorm(d)
  direction <- direction / sqrt(sum(direction^2)) # Normalize to unit vector

  # Calculate the mean vector for the second class
  class2_mean <- direction * dist

  # Simulate data for the second class
  class2_data <- matrix(rnorm(m * d, mean = 0, sd = 1), nrow = m)
  class2_data <- sweep(class2_data, 2, class2_mean, `+`) # Add the mean vector

  # Combine data from both classes
  data <- rbind(class1_data, class2_data)

  # Create numeric labels
  labels_numeric <- c(rep(1, n), rep(2, m))

  # Convert numeric labels to factor
  labels_factor <- factor(labels_numeric, levels = c(1, 2), labels = c("Class 1", "Class 2"))

  # Combine data and labels into a data frame
  df <- data.frame(data)
  df$Class <- labels_factor

  # Return data frame
  return(df)
}
