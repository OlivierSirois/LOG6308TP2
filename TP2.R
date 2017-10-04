m = as.matrix(read.table("citeseer.rtable"))
library("expm")
d <- 0.85
pr <- rep(1, dim(m)[1])

sum.powers.matrix <- function(m, n) {
  powers <- c(1:n)
  res <- Reduce('+', lapply(powers, function(x) m %^% x))
  res[res > 1] <- 1
  
  return(res)
}

remove.autoreferences <- function(m) {
  res <- m
  diag(res) <- 0
}

page.rank.until.stab <- function(m, d, pr){
  pr.next <- page.rank(m, d, pr)
  while(abs(mean(pr.next-pr)) > 0.0001){
    pr <- pr.next
    pr.next <- page.rank(m,d,pr)
  }
  return (pr)
}

page.rank <- function(m, d, pr){
  denum <- (pr/colSums(m))
  denum[denum == Inf] <- 0
  (pr <- (1-d)/3 + (d * (m %*%denum)))
  #print(pr[450])
  return(pr)
}
pr <- page.rank.until.stab(m, d, pr)
print(pr)




pagerank.iteration <- function(refs, n, d, pr) {
  # Number of articles
  n.articles <- dim(refs)[1]
  
  # n-level references
  m <- sum.powers.matrix(refs, n)
  
  # Compute PageRank
  pr.res <- (1-d)/n.articles + (d * (m %*% (pr/colSums(m))))
  
  return(pr.res)
}


m <- matrix(c(0,1,1,0,0,1,1,0,0),3)
pr <- rep(1,3)

m2 <- sum.powers.matrix(m,2)