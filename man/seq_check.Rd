\name{seq_check}
\alias{seq_check}


\title{
Check SMILES strings and amino acid sequences
}

\description{
In real-world cases, most of the data are not complete and contains incorrect values, missing values, and so on. Thus, there may be invalid sequences in the data. This function can find such sequences and remove them from the data. For SMILES strings, the function "webchem::is.smiles" is used. A valid amino acid sequence means a string that only contains capital letters of an alphabet.
}

\usage{
seq_check(smiles = NULL, AAseq = NULL, outcome = NULL)
}

\arguments{
    \item{smiles}{
SMILES strings (default: NULL)
}
    \item{AAseq}{
amino acid sequences (default: NULL)
}
    \item{outcome}{
a variable that indicates how strong two molecules interact with each other or whether there is an interaction between them (default: NULL)
}
}

\value{
valid sequences
}

\author{
Dongmin Jung
}

\seealso{
webchem::is.smiles
}

\references{
Dey, N., Wagh, S., Mahalle, P. N., & Pathan, M. S. (Eds.). (2019). Applied machine learning for smart data analysis. CRC Press.
}

\examples{
seq_check(smiles = example_cpi[1, 1], outcome = example_cpi[1, 3])
}
