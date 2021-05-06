\name{cpi_model}
\alias{fit_cpi}
\alias{predict_cpi}


\title{
Deep learning model fitting and prediction for compound-protein interactions
}

\description{
The model for compound-protein interactions (CPI) takes the pair of SMILES strings of compounds and amino acid sequences (one letter amino acid code) of proteins as input. They are fed into the compound and protein encoders, respectively, and then these encoders are concatenated. Due to the combination of compound and protein encoders, there are many kinds of CPI models. However, the graph neural network such as the graph concolutional network (GCN) is only available for compounds. We need to select one of types of compounds. For graph and fingerprint, the SMILES sequences are not used for encoders, because the information of graph or fingerprint is extracted from the SMILES sequenes and then it is fed into encoders. For sequence, the unigram is used as default, but the n-gram is available only for proteins. Since the CPI model needs some arguments of encoders, we may have to match the names of such arguments.
}

\usage{
fit_cpi(smiles = NULL, AAseq = NULL, outcome,
        convert_canonical_smiles = TRUE,
        compound_type = NULL, compound_max_atoms,
        compound_length_seq, protein_length_seq,
        compound_embedding_dim, protein_embedding_dim,
        protein_ngram_max = 1, protein_ngram_min = 1,
        smiles_val = NULL, AAseq_val = NULL, outcome_val = NULL,
        net_args = list(
            compound,
            compound_args,
            protein,
            protein_args,
            fc_units = c(1),
            fc_activation = c("linear"), ...),
        net_names = list(
            name_compound_max_atoms = NULL,
            name_compound_feature_dim = NULL,
            name_compound_fingerprint_size = NULL,
            name_compound_embedding_layer = NULL,
            name_compound_length_seq = NULL,
            name_compound_num_tokens = NULL,
            name_compound_embedding_dim = NULL,
            name_protein_length_seq = NULL,
            name_protein_num_tokens = NULL,
            name_protein_embedding_dim = NULL),
        preprocessor_only = FALSE,
        preprocessing = list(
            outcome = NULL,
            outcome_val = NULL,
            convert_canonical_smiles = NULL,
            canonical_smiles = NULL,
            compound_type = NULL,
            compound_max_atoms = NULL,
            compound_A_pad = NULL,
            compound_X_pad = NULL,
            compound_A_pad_val = NULL,
            compound_X_pad_val = NULL,
            compound_fingerprint = NULL,
            compound_fingerprint_val = NULL,
            smiles_encode_pad = NULL,
            smiles_val_encode_pad = NULL,
            compound_lenc = NULL,
            compound_length_seq = NULL,
            compound_num_tokens = NULL,
            compound_embedding_dim = NULL,
            AAseq_encode_pad = NULL,
            AAseq_val_encode_pad = NULL,
            protein_lenc = NULL,
            protein_length_seq = NULL,
            protein_num_tokens = NULL,
            protein_embedding_dim = NULL,
            protein_ngram_max = NULL,
            protein_ngram_min = NULL),
        batch_size, use_generator = FALSE,
        validation_split = 0, ...)

predict_cpi(modelRes, smiles = NULL, AAseq = NULL,
            preprocessing = list(
                canonical_smiles = NULL,
                compound_A_pad = NULL,
                compound_X_pad = NULL,
                compound_fingerprint = NULL,
                smiles_encode_pad = NULL,
                AAseq_encode_pad = NULL),
            use_generator = FALSE,
            batch_size = NULL)
}

\arguments{
    \item{smiles}{
SMILES strings, each column for the element of a pair (default: NULL)
}
    \item{AAseq}{
amino acid sequences, each column for the element of a pair (default: NULL)
}
    \item{outcome}{
a variable that indicates how strong two molecules interact with each other or whether there is an interaction between them
}
    \item{convert_canonical_smiles}{
SMILES strings are converted to canonical SMILES strings if TRUE (default: TRUE)
}
    \item{compound_type}{
"graph", "fingerprint" or "sequence"
}
    \item{compound_max_atoms}{
maximum number of atoms for compounds
}
    \item{compound_length_seq}{
length of compound sequence
}
    \item{protein_length_seq}{
length of protein sequence
}
    \item{compound_embedding_dim}{
dimension of the dense embedding for compounds
}
    \item{protein_embedding_dim}{
dimension of the dense embedding for proteins
}
    \item{protein_ngram_max}{
maximum size of an n-gram for protein sequences (default: 1)
}
    \item{protein_ngram_min}{
minimum size of an n-gram for protein sequences (default: 1)
}
    \item{smiles_val}{
SMILES strings for validation (default: NULL)
}
    \item{AAseq_val}{
amino acid sequences for validation (default: NULL)
}
    \item{outcome_val}{
outcome for validation (default: NULL)
}
    \item{net_args}{
list of arguments for compound and protein encoder networks and for fully connected layer
\itemize{
    \item compound : encoder network for compounds
    \item compound_args : arguments of compound encoder
    \item protein : encoder network for proteins
    \item protein_args : arguments of protein encoder
    \item fc_units : dimensionality of the output space in the fully connected layer (default: 1)
    \item fc_activation : activation of the fully connected layer (default: "linear")
    \item ... : arguments of "keras::compile" but for object
}
}
    \item{net_names}{
list of names of arguments used in both the CPI model and encoder networks, names are set to NULL as default
\itemize{
    \item name_compound_max_atoms : corresponding name for the maximum number of atoms in the compound encoder, "max_atoms" if NULL
    \item name_compound_feature_dim : corresponding name for the dimension of node features in the compound encoder, "feature_dim" if NULL
    \item name_compound_fingerprint_size : corresponding name for the length of a fingerprint in the compound encoder, "fingerprint_size" if NULL
    \item name_compound_embedding_layer : corresponding name for the use of the embedding layer in the compound encoder, "embedding_layer" if NULL
    \item name_compound_length_seq : corresponding name for the length of sequences in the compound encoder, "length_seq" if NULL
    \item name_compound_num_tokens : corresponding name for the total number of distinct strings in the compound encoder, "num_tokens" if NULL
    \item name_compound_embedding_dim : corresponding name for dimension of the dense embedding in the compound encoder, "embedding_dim" if NULL
    \item name_protein_length_seq : corresponding name for the length of sequences in the protein encoder, "length_seq" if NULL
    \item name_protein_num_tokens : corresponding name for the total number of distinct strings in the protein encoder, "num_tokens" if NULL
    \item name_protein_embedding_dim : corresponding name for dimension of the dense embedding in the protein encoder, "embedding_dim" if NULL
}
}
    \item{preprocessor_only}{
model is not fitted after preprocessing if TRUE (default: FALSE)
}
    \item{preprocessing}{
list of preprocessed results for "fit_cpi" or "predict_cpi", they are set to NULL as default
\itemize{
    \item outcome : outcome variable
    \item outcome_val : outcome variable for validation
    \item convert_canonical_smiles : canonical representation used for preprocessing if TRUE
    \item canonical_smiles : canonical representation of SMILES
    \item compound_type : "graph", "fingerprint" or "sequence"
    \item compound_max_atoms : maximum number of atoms for compounds
    \item compound_A_pad : padded or turncated adjacency matrix of compounds
    \item compound_X_pad : padded or turncated node features of compounds
    \item compound_A_pad_val : padded or turncated adjacency matrix for validation
    \item compound_X_pad_val : padded or turncated node features for validation
    \item compound_fingerprint : fingerprint of compounds
    \item compound_fingerprint_val : fingerprint for validation
    \item smiles_encode_pad : encoded SMILES sequence which is padded or truncated
    \item smiles_val_encode_pad : encoded SMILES sequence for validation
    \item compound_lenc : encoded labels for characters of SMILES strings
    \item compound_length_seq : length of compound sequence
    \item compound_num_tokens : total number of characters of compounds
    \item compound_embedding_dim : dimension of the dense embedding for compounds
    \item AAseq_encode_pad : encoded amino acid sequence which is padded or truncated
    \item AAseq_val_encode_pad : encoded amino acid sequence for validation
    \item protein_lenc : encoded labels for characters of amino acid sequenes
    \item protein_length_seq : length of protein sequence
    \item protein_num_tokens : total number of characters of proteins
    \item protein_embedding_dim : dimension of the dense embedding for proteins
    \item protein_ngram_max : maximum size of an n-gram for protein sequences
    \item protein_ngram_min : minimum size of an n-gram for protein sequences
    \item removed_smiles : index for removed smiles while checking
    \item removed_AAseq : index for removed AAseq while checking
    \item removed_smiles_val : index for removed smiles of validation
    \item removed_AAseq_val : index for removed AAseq of validation
}
}
    \item{batch_size}{
batch size
}
    \item{use_generator}{
use data generator if TRUE (default: FALSE)
}
    \item{validation_split}{
proportion of validation data, it is ignored when there is a validation set (default: 0)
}
    \item{modelRes}{
result of the "fit_cpi"
}
    \item{...}{
additional parameters for the "keras::fit" or "keras::fit_generator"
}
}


\value{
model
}


\author{
Dongmin Jung
}


\seealso{
keras::compile, keras::fit, keras::fit_generator, keras::layer_dense, keras::keras_model, purrr::pluck, webchem::is.smiles
}


\examples{
compound_max_atoms <- 50
protein_embedding_dim <- 16
protein_length_seq <- 100
gcn_cnn_cpi <- fit_cpi(
    smiles = example_cpi[1:100, 1],
    AAseq = example_cpi[1:100, 2],
    outcome = example_cpi[1:100, 3],
    compound_type = "graph",
    compound_max_atoms = compound_max_atoms,
    protein_length_seq = protein_length_seq,
    protein_embedding_dim = protein_embedding_dim,
    net_args = list(
        compound = "gcn_in_out",
        compound_args = list(
            gcn_units = c(128, 64),
            gcn_activation = c("relu", "relu"),
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
        metrics = "accuracy"),
    epochs = 2, batch_size = 16)
pred <- predict_cpi(gcn_cnn_cpi, example_cpi[101:110, 1], example_cpi[101:110, 2])

gcn_cnn_cpi2 <- fit_cpi(
    preprocessing = gcn_cnn_cpi$preprocessing,
    net_args = list(
        compound = "gcn_in_out",
        compound_args = list(
            gcn_units = c(128, 64),
            gcn_activation = c("relu", "relu"),
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
        metrics = "accuracy"),
    epochs = 2, batch_size = 16)
pred <- predict_cpi(gcn_cnn_cpi2, preprocessing = pred$preprocessing)
}
