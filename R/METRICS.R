#' @title 5 x 2 CV F-Test
#'
#' @description \code{cv52test} calculates the 5 x 2 CV F-Test as well as the
#' 10 t-tests for the difference between the two models.
#'
#' @param m1 A matrix of the first model's predictions.
#' @param m2 A matrix of the second model's predictions.
#'
#' @return a vector where the first value is the p-value for the F-test
#' and the next 10 values are the p-values for the t-tests.
#'
#' @author Jordan Eckert
#'
#' @rdname cv52test
#' @export

cv52test <- function(m1,m2)
{
  md <- m1-m2
  sumpi <- sum(md^2)
  mpi <- apply(md,1,mean)
  sdpi <- sum((md-mpi)^2)
  Fres <- sumpi/(2*sdpi)
  Ftest <- 1 - pf(Fres,10,5)
  tres <- md/sqrt(sdpi/5)
  ttest <- 1 - pt(tres,5)
  return(c(Ftest,ttest))
}

#' @title AUC of the ROC curve for two class setting.
#'
#' @description \code{auc} calculates AUC of the ROC curve.
#' It evaluates the model's ability to distinguish between positive and negative
#' classes across all classification thresholds.
#'
#' Best use case: Evaluating and comparing the overall discriminative power of models.
#'
#' @param predict A vector of predicted class labels. Must be labeled as 1 and 2.
#' @param actual A vector of actual class labels. Must be labeled as 1 and 2.
#'
#' @return AUC of the ROC curve. Values range from 0 to 1
#'
#' @author Jordan Eckert
#'
#' @examples
#' n <- 1000
#' x1 <- runif(n, 1, 10)
#' x2 <- runif(n, 1, 10)
#' x <- cbind(x1, x2)
#' y <- as.factor(ifelse(3 < x1 & x1 < 7 & 3 < x2 & x2 < 7, "A", "B"))
#'
#' #' # testing the performance
#' i_train <- sample(1:n, round(n*0.8))
#'
#' x_train <- x[i_train,]
#' y_train <- y[i_train]
#'
#' x_test <- x[-i_train,]
#' y_test <- y[-i_train]
#'
#' m_pcccd <- pcccd(x = x_train, y = y_train, tau = 1)
#' p_pred <- classify_pcccd(pcccd = m_pcccd, newdata = x_test)
#'
#' auc(as.numeric(as.factor(p_pred)), as.numeric(as.factor(y_test)))
#'
#' @rdname auc
#' @export

auc <- function(predict,actual){

  if(!is.vector(predict)) predict <- as.vector(predict)
  if(!is.vector(actual)) actual <- as.vector(actual)

  if(length(unique(actual)) > 2)
    stop("you need at most 2 classes")
  if(length(unique(predict)) == 1)
    return(0.5)

  ind <- which(predict==2)
  TP  <- sum(actual[ind]==2)
  FP  <- sum(actual[ind]==1)
  FN  <- sum(actual[-ind]==2)
  TN  <- sum(actual[-ind]==1)
  rate1 <- TP/(TP+FN)
  rate2 <- FP/(FP+TN)
  return((1+rate1-rate2)/2)
}

#' @title F1 - Score
#'
#' @description \code{f1_score} calculates F1 score for a model's prediction.
#' Focuses on the balance between false positives (Precision) and false negatives (Recall).
#' Harmonic mean of precision and recall. Places more emphasis on correctly identifying the minority class.
#' Does not directly account for true negatives. No symmetry between classes.
#'
#' Best use case: Use F1 Score when dealing with imbalanced datasets and focusing on the positive class.
#'
#' @param predict A vector of predicted class labels. Must be labeled as 1 and 2.
#' @param actual A vector of actual class labels. Must be labeled as 1 and 2.
#'
#' @return F1 score. Values range from 0 to 1.
#'
#' @examples
#' n <- 1000
#' x1 <- runif(n, 1, 10)
#' x2 <- runif(n, 1, 10)
#' x <- cbind(x1, x2)
#' y <- as.factor(ifelse(3 < x1 & x1 < 7 & 3 < x2 & x2 < 7, "A", "B"))
#'
#' #' # testing the performance
#' i_train <- sample(1:n, round(n*0.8))
#'
#' x_train <- x[i_train,]
#' y_train <- y[i_train]
#'
#' x_test <- x[-i_train,]
#' y_test <- y[-i_train]
#'
#' m_pcccd <- pcccd(x = x_train, y = y_train, tau = 1)
#' p_pred <- classify_pcccd(pcccd = m_pcccd, newdata = x_test)
#'
#' f1_score(as.numeric(as.factor(p_pred)), as.numeric(as.factor(y_test)))
#'
#'
#' @author Jordan Eckert
#'
#' @rdname f1_score
#' @export

f1_score <- function(predict, actual) {

  if(!is.vector(predict)) predict <- as.vector(predict)
  if(!is.vector(actual)) actual <- as.vector(actual)

  if(length(unique(actual)) > 2)
    stop("you need at most 2 classes")

  TP <- sum(predict == 2 & actual == 2, na.rm = TRUE)
  FP <- sum(predict == 2 & actual == 1, na.rm = TRUE)
  FN <- sum(predict == 1 & actual == 2, na.rm = TRUE)

  precision <- ifelse((TP + FP) == 0, 0, TP / (TP + FP))
  recall <- ifelse((TP + FN) == 0, 0, TP / (TP + FN))

  f1 <- ifelse((precision + recall) == 0, 0, 2 * (precision * recall) / (precision + recall))

  return(f1)
}

#' @title Geometric Mean
#'
#' @description \code{g_mean} calculates geometric mean of sensitivity (recall) and specificity.
#' Indicates the balance between correctly identifying positive and negative classes.
#' Accounts for both sensitivity and specificity, ensuring no bias toward any class.
#' Less commonly used compared to AUC and F1 Score, making it harder to compare across studies.
#'
#' Best use case: When maintaining performance for both classes is equally important, especially in imbalanced datasets.
#'
#' @param predict A vector of predicted class labels. Must be labeled as 1 and 2.
#' @param actual A vector of actual class labels. Must be labeled as 1 and 2.
#'
#' @return Geometric Mean. Values range from 0 to 1
#'
#' @author Jordan Eckert
#'
#' @examples
#' n <- 1000
#' x1 <- runif(n, 1, 10)
#' x2 <- runif(n, 1, 10)
#' x <- cbind(x1, x2)
#' y <- as.factor(ifelse(3 < x1 & x1 < 7 & 3 < x2 & x2 < 7, "A", "B"))
#'
#' #' # testing the performance
#' i_train <- sample(1:n, round(n*0.8))
#'
#' x_train <- x[i_train,]
#' y_train <- y[i_train]
#'
#' x_test <- x[-i_train,]
#' y_test <- y[-i_train]
#'
#' m_pcccd <- pcccd(x = x_train, y = y_train, tau = 1)
#' p_pred <- classify_pcccd(pcccd = m_pcccd, newdata = x_test)
#'
#' g_mean(as.numeric(as.factor(p_pred)), as.numeric(as.factor(y_test)))
#'
#' @rdname g_mean
#' @export

g_mean <- function(predict, actual) {

  if(!is.vector(predict)) predict <- as.vector(predict)
  if(!is.vector(actual)) actual <- as.vector(actual)

  if(length(unique(actual)) > 2)
    stop("you need at most 2 classes")

  TP <- sum(predict == 2 & actual == 2, na.rm = TRUE)
  TN <- sum(predict == 1 & actual == 1, na.rm = TRUE)
  FP <- sum(predict == 2 & actual == 1, na.rm = TRUE)
  FN <- sum(predict == 1 & actual == 2, na.rm = TRUE)

  sensitivity <- ifelse((TP + FN) == 0, 0, TP / (TP + FN))
  specificity <- ifelse((TN + FP) == 0, 0, TN / (TN + FP))

  g_mean <- sqrt(sensitivity * specificity)

  return(g_mean)
}

#' @title Matthew's Correlation Coefficient
#'
#' @description \code{mcc} calculates Matthew's Correlation Coefficient.
#' MCC is essentially the correlation coefficient between the predicted and actual classifications.
#' Unlike F1 Score, which focuses on the positive class,
#' MCC evaluates the correlation between true and predicted labels,
#' ensuring both classes are treated symmetrically.
#' Can be undefined if the confusion matrix contains all zeros for certain
#' categories (e.g., no true positives or no false positives).
#'
#' Best use case: Both positive and negative outcomes are equally important.
#'
#' @param predict A vector of predicted class labels. Must be labeled as 1 and 2.
#' @param actual A vector of actual class labels. Must be labeled as 1 and 2.
#'
#' @return Matthew's Correlation Coefficient. Values range from -1 to 1.
#' 1 indicates perfect prediction, 0 indicates random prediction, and -1 indicates total disagreement.
#'
#' @author Jordan Eckert
#'
#' @examples
#' n <- 1000
#' x1 <- runif(n, 1, 10)
#' x2 <- runif(n, 1, 10)
#' x <- cbind(x1, x2)
#' y <- as.factor(ifelse(3 < x1 & x1 < 7 & 3 < x2 & x2 < 7, "A", "B"))
#'
#' #' # testing the performance
#' i_train <- sample(1:n, round(n*0.8))
#'
#' x_train <- x[i_train,]
#' y_train <- y[i_train]
#'
#' x_test <- x[-i_train,]
#' y_test <- y[-i_train]
#'
#' m_pcccd <- pcccd(x = x_train, y = y_train, tau = 1)
#' p_pred <- classify_pcccd(pcccd = m_pcccd, newdata = x_test)
#'
#' mcc(as.numeric(as.factor(p_pred)), as.numeric(as.factor(y_test)))
#'
#' @rdname mcc
#' @export

mcc <- function(predict, actual) {

  # Ensure inputs are vectors
  if (!is.vector(predict)) predict <- as.vector(predict)
  if (!is.vector(actual)) actual <- as.vector(actual)

  # Check for binary classification
  classes <- sort(unique(actual))
  if (length(classes) != 2) stop("You need exactly 2 classes")

  # Recode labels to 0 (neg) and 1 (pos) using actual classes
  predict_bin <- as.numeric(predict == classes[2])
  actual_bin  <- as.numeric(actual == classes[2])

  # Handle NA predictions by treating them as wrong
  na_idx <- is.na(predict_bin)
  predict_bin[na_idx] <- 1 - actual_bin[na_idx]  # force to be incorrect

  # Compute confusion matrix components
  TP <- sum(predict_bin == 1 & actual_bin == 1)
  TN <- sum(predict_bin == 0 & actual_bin == 0)
  FP <- sum(predict_bin == 1 & actual_bin == 0)
  FN <- sum(predict_bin == 0 & actual_bin == 1)

  # Avoid integer overflow
  TP <- as.double(TP)
  TN <- as.double(TN)
  FP <- as.double(FP)
  FN <- as.double(FN)

  # Calculate MCC
  numerator <- (TP * TN) - (FP * FN)
  denominator <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

  mcc <- ifelse(denominator == 0, 0, numerator / denominator)
  return(mcc)
}
