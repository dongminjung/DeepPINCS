\name{get_fingerprint}
\alias{get_fingerprint}


\title{
Molecular fingerprint of compounds from SMILES strings
}

\description{
A molecular fingerprint is a way of encoding the structural features of a molecule. The most common type of fingerprint is a sequence of ones and zeros. Fingerprints are special kinds of descriptors that characterize a molecule and its properties as a binary bit vector that represents the presence or absence of particular substructure in the molecule. For such a fingerprint, the Chemistry Development Kit (CDK) is used as a cheminformatics tool.
}

\usage{
get_fingerprint(smiles, ...)
}

\arguments{
    \item{smiles}{
SMILES strings
}
    \item{...}{
arguments for "rcdk::get.fingerprint" but for molecule
}
}

\value{
a fingerprint of a compound
}

\author{
Dongmin Jung
}

\seealso{
rcdk::get.fingerprint, rcdk::parse.smiles
}

\references{
Balakin, K. V. (2009). Pharmaceutical data mining: approaches and applications for drug discovery. Wiley.
}

\examples{
get_fingerprint(example_cpi[1, 1])
}
