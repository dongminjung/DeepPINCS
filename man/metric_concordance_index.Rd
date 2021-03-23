\name{metric_concordance_index}
\alias{metric_concordance_index}

\title{
Concordance index
}

\description{
The concordance index or c-index can be seen as one of the model performance metrics. It represents a good fit of the model.
}

\author{
Dongmin Jung
}

\seealso{
keras::k_cast, keras::k_equal, keras::k_sum, tensorflow::tf
}

\references{
Kose, U., & Alzubi, J. (2020). Deep learning for cancer diagnosis. Springer.
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
    keras::layer_dense(units = 1, activation = "linear")
model <- keras::keras_model(c(inputA, inputX), x) \%>\%
    keras::compile(loss = "mean_squared_error",
        optimizer = keras::optimizer_adam(),
        metrics = custom_metric("concordance_index",
            metric_concordance_index))

X <- example_bioassay[1:100, 1]
graph_structure_node_feature <- get_graph_structure_node_feature(X, 10)
X <- list(graph_structure_node_feature$A_pad,
    graph_structure_node_feature$X_pad)
y <- example_bioassay[1:100, 4]
model \%>\% fit(X, y, epochs = 2, batch_size = 16)
}
