mda <- function(y, x) {
  # Multiples Discriminants Analysis for n classes.
  # Classes must be defined in y and
  # named as consecutive positive integers (base+1, base+2, ..., base+nclass).
  base <- min(y) - 1
  nclass <- length(unique(y))
  sets <- list()
  for (c in 1:nclass) sets[[c]] <- x[y == c + base, -1]
  ms <- list()
  for (c in 1:nclass) ms[[c]] <- colMeans(sets[[c]])
  m <- colMeans(x)
  nfeat <- ncol(x)
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
  values <- Re(decomp$values)
  vectors <- Re(decomp$vectors)
  return(list("values"=values, "vectors"=vectors, "x"=x %*% vectors))
}
