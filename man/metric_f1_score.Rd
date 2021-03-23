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
        metrics = custom_metric("F1_score",
            metric_f1_score)),
    epochs = 2,
    batch_size = 16)
}
