% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tiplength.R
\name{tiplength}
\alias{tiplength}
\title{Extracts the length of a tip from a tree}
\usage{
tiplength(tree, tipname)
}
\arguments{
\item{tree}{a phylogenetic tree (as a \code{phylo} object)}

\item{tipname}{the tip name, as a character string, or a numeric index}
}
\value{
The tip length (as a \code{double}).
}
\description{
\code{tiplength} returns the length of a tip of a tree, given either (a) the tip name or
(b) the index of the tip (in the order of the tip labels in the tree).
}
\author{
Simon Frost (\email{sdwfrost@gmail.com})
}

