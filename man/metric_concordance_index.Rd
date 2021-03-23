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
compound_length_seq <- 50
compound_embedding_dim <- 16
protein_embedding_dim <- 16
protein_length_seq <- 100

mlp_cnn_cpi <- fit_cpi(
    smiles = example_cpi[1:100, 1],
    AAseq = example_cpi[1:100, 2], 
    outcome = example_cpi[1:100, 3],
    compound_type = "sequence",
    compound_length_seq = compound_length_seq,
    compound_embedding_dim = compound_embedding_dim,
    protein_length_seq = protein_length_seq,
    protein_embedding_dim = protein_embedding_dim,
    net_args = list(
    compound = "mlp_in_out",
    compound_args = list(
            fc_units = c(10),
            fc_activation = c("relu")),
        protein = "cnn_in_out",
        protein_args = list(
            cnn_filters = c(32),
            cnn_kernel_size = c(3),
            cnn_activation = c("relu"),
            fc_units = c(10),
            fc_activation = c("relu")),
        fc_units = c(1),
        fc_activation = c("sigmoid"),
        loss = "binary_crossentropy",
        optimizer = keras::optimizer_adam(),
        metrics = custom_metric("concordance_index",
            metric_concordance_index)),
    epochs = 2,
    batch_size = 16)
}
