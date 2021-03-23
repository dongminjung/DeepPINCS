\name{metric_f1_score}
\alias{metric_f1_score}

\title{
F1-score
}

\description{
The F1-score is a metric combining precision and recall. It is typically used instead of accuracy in the case of severe class imbalance in the dataset. The higher the values of F1-score, the better the validation of the model.
}

\author{
Dongmin Jung
}

\seealso{
keras::k_equal, keras::k_sum, tensorflow::tf
}

\references{
Kubben, P., Dumontier, M., & Dekker, A. (2019). Fundamentals of clinical data science. Springer.

Mishra, A., Suseendran, G., & Phung, T. N. (Eds.). (2020). Soft Computing Applications and Techniques in Healthcare. CRC Press.
}

\examples{
layer_multi_linear <- Layer(
    classname <- "MultiLinear", 
    initialize <- function(units, ...) {
        super()$`__init__`(...)
        self$units <- units
    },
    build <- function(input_shape) {
        self$kernel <- self$add_weight(shape = shape(input_shape[3], self$units))
        self$bias <- self$add_weight(shape = shape(self$units))
    },
    call <- function(inputs, ...) {
        keras::k_dot(inputs, self$kernel) + self$bias
    },
    get_config <- function() {
        list(
            units <- self$units,
            name <- self$name
        )
    }
)

max_atoms <- 10
feature_dim <- 24
inputA <- keras::layer_input(shape = c(max_atoms, max_atoms))
inputX <- keras::layer_input(shape = c(max_atoms, feature_dim))
x <- keras::layer_dot(c(inputA, inputX), axes = 0)
x  <- x \%>\% 
    layer_multi_linear(units = 10) \%>\%
    keras::layer_activation(activation = "relu")
x <- keras::layer_global_average_pooling_1d()(x) \%>\% 
    keras::layer_dense(units = 1, activation = "sigmoid")
model <- keras::keras_model(c(inputA, inputX), x) \%>\%
    keras::compile(loss = "mean_squared_error",
        optimizer = keras::optimizer_adam(),
        metrics = custom_metric("F1_score", metric_f1_score))

X <- example_bioassay[1:100, 1]
graph_structure_node_feature <- get_graph_structure_node_feature(X, 10)
X <- list(graph_structure_node_feature$A_pad,
    graph_structure_node_feature$X_pad)
y <- example_bioassay[1:100, 4]
model \%>\% fit(X, y, epochs = 2, batch_size = 16)
}
