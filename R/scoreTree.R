#' Scores an observed tree against a known simulation. Scores are the quantile the observed tree metric is in
#'
#' \code{scoreTree} returns the scoring of an observed tree
#' Each row is a sequential pruning, each column a metric from allMetrics.
#'
#' @param obs_tree a phylogenetic tree to score (as a \code{phylo} object)
#' @param pruned a set of trees from prunePermute to score against (as a \code{list} of \code{phylo} objects)
#' @return An \code{array} of scores.
#' @author John Lees (\email{lees.john6@@gmail.com})
#' @export
#' @examples
#' sim_tree <- simulate()
#' sim_tree2 <- simulate()
#' pruned <- prunePermute(sim_tree$tree)
#' obs_tree <- sim_tree2$tree
#' scoreTree(obs_tree, pruned)

scoreTree <- function(obs_tree, pruned_trees)
{
  # metrics calculated for simulated and observed passed trees
  # TODO will the sim_metrics be pre-calculated and passed? Hopefully, but
  # assume not for now
  sim_metrics <- lapply(pruned$seqPrunes, lapply, allMetrics)
  real_metrics <- lapply(timeprune(obs_tree)$trees, allMetrics)

  # Loop over each metric within each pruning
  scores <- array(dim=c(length(real_metrics),length(real_metrics[[1]])))
  for (prunedTree in 1:length(real_metrics))
  {
    for (metric in 1:length(real_metrics[[1]]))
    {
      # This gives the quantile the observed score is in
      scoreRange <- rep(NA, length=length(sim_metrics[[1]]))
      for (permutation in 1:length(sim_metrics[[1]]))
      {
        scoreRange[permutation] <- sim_metrics[[prunedTree]][[permutation]][[metric]]
      }
      scoreRange <- sort(scoreRange)
      scores[prunedTree,metric] <-
        sum(scoreRange < real_metrics[[prunedTree]][[metric]])/length(scoreRange)
    }
  }

  return(scores)
}

# Wrapper for branchLengthMetrics and imbalanceMetrics
allMetrics <- function(tree)
{
  all_metrics <- c(imbalanceMetrics(tree), branchLengthMetrics(tree))
  return(all_metrics)
}
