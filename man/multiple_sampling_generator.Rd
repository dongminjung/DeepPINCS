\name{multiple_sampling_generator}
\alias{multiple_sampling_generator}

\title{
Generator function for multiple inputs
}
\description{
This is a generator function that yields batches of data with multiple inputs.
}


\usage{
multiple_sampling_generator(X_data, Y_data = NULL, batch_size,
    shuffle = TRUE)
}


\arguments{
    \item{X_data}{
list of multiple inputs
}
    \item{Y_data}{
targets (default: NULL)
}
    \item{batch_size}{
batch size
}
    \item{shuffle}{
whether to shuffle the data or not (default: TRUE)
}
}


\value{
generator for "keras::fit_generator" or "keras::predict_genertor"
}


\author{
Dongmin Jung
}


\examples{
X_data <- c(list(matrix(rnorm(200), ncol = 2)),
    list(matrix(rnorm(200), ncol = 2)))
Y_data <- matrix(rnorm(100), ncol = 1)
multiple_sampling_generator(X_data, Y_data, 32)
}


