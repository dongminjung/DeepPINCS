\name{get_graph_structure_node_feature}
\alias{get_graph_structure_node_feature}


\title{
Graph structure and node features from SMILES strings
}

\description{
In molecular graph representations, nodes represent atoms and edges represent bonds. For molecular features, the Chemistry Development Kit (CDK) is used as a cheminformatics tool. The degree of an atom in the graph representation and the atomic symbol and implicit hydrogen count for an atom are used as molecular features.
}

\usage{
get_graph_structure_node_feature(smiles, max_atoms,
    element_list = c(
        "C", "N", "O", "S", "F", "Si", "P", "Cl",
        "Br", "Mg", "Na", "Ca", "Fe",  "Al", "I",
        "B", "K", "Se", "Zn", "H", "Cu", "Mn"))
}

\arguments{
    \item{smiles}{
SMILES strings
}
    \item{max_atoms}{
maximum number of atoms
}
    \item{element_list}{
list of atom symbols
}
}

\value{
    \item{A_pad}{
a padded or turncated adjacency matrix for each SMILES string
}
    \item{X_pad}{
a padded or turncated node features for each SMILES string
}
    \item{feature_dim}{
dimension of node features
}
    \item{element_list}{
list of atom symbols
}
}

\author{
Dongmin Jung
}

\seealso{
matlab::padarray, purrr::chuck, rcdk::get.adjacency.matrix, rcdk::get.atoms,
rcdk::get.hydrogen.count, rcdk::get.symbol rcdk::parse.smiles
}

\references{
Balakin, K. V. (2009). Pharmaceutical data mining: approaches and applications for drug discovery. Wiley.
}

\examples{
get_graph_structure_node_feature(example_cpi[1, 1], 10)
}
