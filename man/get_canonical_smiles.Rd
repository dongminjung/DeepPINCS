\name{get_canonical_smiles}
\alias{get_canonical_smiles}


\title{
Convert SMILES strings to canonical SMILES strings
}

\description{
There may be many different ways to construct the SMILES string for a given molecule. A canonical representation is a unique ordering of the atoms for a given molecular graph.
}

\usage{
get_canonical_smiles(smiles)
}

\arguments{
    \item{smiles}{
SMILES strings
}
}

\value{
canonical representation of SMILES
}

\author{
Dongmin Jung
}

\seealso{
rcdk::parse.smile, rcdk::get.smiles, rcdk::smiles.flavors
}

\references{
Leach, A. R., & Gillet, V. J. (2007). An introduction to chemoinformatics. Springer.
}

\examples{
get_canonical_smiles(example_cpi[1, 1])
}
