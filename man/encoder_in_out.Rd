\name{encoder_in_out}
\alias{gcn_in_out}
\alias{rnn_in_out}
\alias{cnn_in_out}
\alias{mlp_in_out}


\title{
Input and output tensors of encoders
}

\description{
The graph convolutional network (GCN), recurrent neural network (RNN), convolutional neural network (CNN), and multilayer perceptron (MLP) are used as encoders. The last layer of the encoders is the fully connected layer. The units and activation can be vectors and the length of the vectors represents the number of layers.
}

\usage{
gcn_in_out(max_atoms, feature_dim, gcn_units, gcn_activation,
    fc_units, fc_activation)

rnn_in_out(length_seq, fingerprint_size, embedding_layer = TRUE,
    num_tokens, embedding_dim, rnn_type, rnn_bidirectional,
    rnn_units, rnn_activation, fc_units, fc_activation)

cnn_in_out(length_seq, fingerprint_size, embedding_layer = TRUE,
    num_tokens, embedding_dim, cnn_filters, cnn_kernel_size, cnn_activation,
    fc_units, fc_activation)

mlp_in_out(length_seq, fingerprint_size, embedding_layer = TRUE,
    num_tokens, embedding_dim, fc_units, fc_activation)
}

\arguments{
    \item{max_atoms}{
maximum number of atoms for gcn
}
    \item{feature_dim}{
dimension of atom features for gcn
}
    \item{gcn_units}{
dimensionality of the output space in the gcn layer
}
    \item{gcn_activation}{
activation of the gcn layer
}
    \item{fingerprint_size}{
the length of a fingerprint
}
    \item{embedding_layer}{
use the embedding layer if TRUE (default: TRUE)
}
    \item{embedding_dim}{
a non-negative integer for dimension of the dense embedding
}
    \item{length_seq}{
length of input sequences
}
    \item{num_tokens}{
total number of distinct strings
}
    \item{cnn_filters}{
dimensionality of the output space in the cnn layer
}
    \item{cnn_kernel_size}{
length of the 1D convolution window in the cnn layer
}
    \item{cnn_activation}{
activation of the cnn layer
}
    \item{rnn_type}{
"lstm" or "gru"
}
    \item{rnn_bidirectional}{
use the bidirectional wrapper for rnn if TRUE
}
    \item{rnn_units}{
dimensionality of the output space in the rnn layer
}
    \item{rnn_activation}{
activation of the rnn layer
}
    \item{fc_units}{
dimensionality of the output space in the fully connected layer
}
    \item{fc_activation}{
activation of the fully connected layer
}
}

\value{
input and output tensors of encoders
}

\author{
Dongmin Jung
}

\seealso{
keras::layer_activation, keras::bidirectional, keras::layer_conv_1d, keras::layer_dense, keras::layer_dot, keras::layer_embedding, keras::layer_global_average_pooling_1d, keras::layer_input, keras::layer_lstm, keras::layer_gru, keras::layer_flatten
}

\examples{
gcn_in_out(max_atoms = 50,
    feature_dim = 50,
    gcn_units = c(128, 64),
    gcn_activation = c("relu", "relu"),
    fc_units = c(10),
    fc_activation = c("relu"))
}
