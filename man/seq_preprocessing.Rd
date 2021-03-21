\name{seq_preprocessing}
\alias{seq_preprocessing}


\title{
Preprocessing for SMILES strings and amino acid sequences
}

\description{
Preprocessing helps make the data suitable for the model depending on the type of data the preprocessing works upon. Preprocessing is more time consuming for text data. The adjacency matrix and node feature, fingerprint, or string data are preprocessed from sequences.
}

\usage{
seq_preprocessing(smiles = NULL,
    AAseq = NULL,
    type,
    convert_canonical_smiles,
    max_atoms,
    length_seq,
    lenc = NULL,
    ngram_max = 1,
    ngram_min = 1)

}

\arguments{
    \item{smiles}{
SMILES strings (default: NULL)
}
    \item{AAseq}{
amino acid sequences (default: NULL)
}
    \item{type}{
"graph", "fingerprint" or "sequence"
}
    \item{convert_canonical_smiles}{
SMILES strings are converted to canonical SMILES strings if TRUE
}
    \item{max_atoms}{
maximum number of atoms for compounds
}
    \item{length_seq}{
length of compound or protein sequence
}
    \item{lenc}{
encoded labels for characters of SMILES strings or amino acid sequenes (default: NULL)
}
    \item{ngram_max}{
maximum size of an n-gram for protein sequences (default: 1)
}
    \item{ngram_min}{
minimum size of an n-gram for protein sequences (default: 1)
}
}

\value{
    \item{canonical_smiles}{
canonical representation of SMILES
}
    \item{convert_canonical_smiles}{
canonical representation is used or not
}
    \item{A_pad}{
padded or turncated adjacency matrix of compounds if type is "graph"
}
    \item{X_pad}{
padded or turncated node features of compounds if type is "graph"
}
    \item{fp}{
fingerprint of compounds if type is "fingerprint"
}
    \item{sequences_encode_pad}{
encoded sequences which are padded or truncated
}
    \item{lenc}{
encoded labels for characters of SMILES strings or amino acid sequenes
}
    \item{length_seq}{
length of compound or protein sequence
}
    \item{num_tokens}{
total number of characters of compounds or proteins
}
}

\author{
Dongmin Jung
}

\references{
Dey, N., Wagh, S., Mahalle, P. N., & Pathan, M. S. (Eds.). (2019). Applied machine learning for smart data analysis. CRC Press.
}

\examples{
seq_preprocessing(smiles = cbind(example_cpi[1, 1]),
    type = "fingerprint",
    convert_canonical_smiles = TRUE)
}
