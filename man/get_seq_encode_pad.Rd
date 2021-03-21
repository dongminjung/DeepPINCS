\name{get_seq_encode_pad}
\alias{get_seq_encode_pad}


\title{
Vectorization of characters of strings
}

\description{
A vectorization of characters of strings is necessary. Vectorized characters are padded or truncated.
}

\usage{
get_seq_encode_pad(sequences, length_seq, ngram_max = 1, ngram_min = 1,
    lenc = NULL)
}

\arguments{
    \item{sequences}{
SMILE strings or amino acid sequences
}
    \item{length_seq}{
length of input sequences
}
    \item{ngram_max}{
maximum size of an n-gram (default: 1)
}
    \item{ngram_min}{
minimum size of an n-gram (default: 1)
}
    \item{lenc}{
encoded labels for characters, LableEncoder object fitted by "CatEncoders::LabelEncoder.fit" (default: NULL)
}
}

\value{
    \item{sequences_encode_pad}{
for each SMILES string, an encoded sequence which is padded or truncated
}
    \item{lenc}{
encoded labels for characters
}
    \item{num_token}{
total number of characters
}
}

\author{
Dongmin Jung
}

\seealso{
CatEncoders::LabelEncoder.fit, CatEncoders::transform, keras::pad_sequences,
stringdist::qgrams, tokenizers::tokenize_ngrams
}

\examples{
get_seq_encode_pad(example_cpi[1, 2], 10)
}
