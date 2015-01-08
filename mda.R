mda <- function(x, level=0.99) {
  # Multiples Discriminants Analysis for n classes.
  # Assumes classes are defined in the first column of x and
  # named as consecutive positive integers (base+1, base+2, ..., base+nclass).
  group <- x[, 1]
  base <- min(group) - 1
  nclass <- length(unique(group))
  sets <- list()
  for (c in 1:nclass) sets[[c]] <- x[group == c + base, -1]
  ms <- list()
  for (c in 1:nclass) ms[[c]] <- colMeans(sets[[c]])
  m <- colMeans(x[, -1])
  nfeat <- ncol(x[, -1])
  Sinter <- matrix(0, ncol=nfeat, nrow=nfeat)
  for (c in 1:nclass) Sinter <- Sinter + (ms[[c]] - m) %*% t(ms[[c]] - m)
  Sintra <- matrix(0, ncol=nfeat, nrow=nfeat)
  for (c in 1:nclass) {
    for (i in 1:ncol(sets[[c]])) {
      Sintra <- Sintra + (sets[[c]][i, ] - ms[[c]]) %*% t(sets[[c]][i, ] - ms[[c]])
    }
  }
  W <- solve(Sintra) %*% Sinter
  decomp <- eigen(W)
  nvec <- sum(cumsum(Re(decomp$values)) / sum(Re(decomp$values)) < level)
  if (nvec < 1) nvec <- 1
  Wtr <- Re(decomp$vectors[, 1:nvec])
  return(cbind(x[, 1], x[, -1] %*% Wtr))
}
