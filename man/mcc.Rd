% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/METRICS.R
\name{mcc}
\alias{mcc}
\title{Matthew's Correlation Coefficient}
\usage{
mcc(predict, actual)
}
\arguments{
\item{predict}{A vector of predicted class labels. Must be labeled as 1 and 2.}

\item{actual}{A vector of actual class labels. Must be labeled as 1 and 2.}
}
\value{
Matthew's Correlation Coefficient. Values range from -1 to 1.
1 indicates perfect prediction, 0 indicates random prediction, and -1 indicates total disagreement.
}
\description{
\code{mcc} calculates Matthew's Correlation Coefficient.
MCC is essentially the correlation coefficient between the predicted and actual classifications.
Unlike F1 Score, which focuses on the positive class,
MCC evaluates the correlation between true and predicted labels,
ensuring both classes are treated symmetrically.
Can be undefined if the confusion matrix contains all zeros for certain
categories (e.g., no true positives or no false positives).

Best use case: Both positive and negative outcomes are equally important.
}
\examples{
n <- 1000
x1 <- runif(n, 1, 10)
x2 <- runif(n, 1, 10)
x <- cbind(x1, x2)
y <- as.factor(ifelse(3 < x1 & x1 < 7 & 3 < x2 & x2 < 7, "A", "B"))

#' # testing the performance
i_train <- sample(1:n, round(n*0.8))

x_train <- x[i_train,]
y_train <- y[i_train]

x_test <- x[-i_train,]
y_test <- y[-i_train]

m_pcccd <- pcccd(x = x_train, y = y_train, tau = 1)
p_pred <- classify_pcccd(pcccd = m_pcccd, newdata = x_test)

mcc(as.numeric(as.factor(p_pred)), as.numeric(as.factor(y_test)))

}
\author{
Jordan Eckert
}
