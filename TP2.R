# Load libraries
library("expm")

sum.powers.matrix <- function(m, n) {
    powers <- c(1:n)
    res <- Reduce('+', lapply(powers, function(x) m %^% x))
    res[res > 1] <- 1

    return(res)
}

remove.autoreferences <- function(m) {
    res <- m
}

pagerank.iteration <- function(refs, n, d, pr) {
    # Number of articles
    n.articles <- dim(refs)[1]

    # n-level references
    m <- sum.powers.matrix(refs, n)

    # Compute PageRank
    pr.res <- (1-d)/n.articles + (d * (m %*% (pr/colSums(m))))

    return(pr.res)
}

# Read source file
data <- read.table("citeseer.rtable")

# Cast data to matrix
references <- as.matrix(data)

m <- matrix(c(0,1,1,0,0,1,1,0,0),3)
pr <- rep(1,3)

pr <- pagerank.iteration(m, 2, 0.85, pr)
print(pr)

pr <- pagerank.iteration(m, 2, 0.85, pr)
print(pr)

pr <- pagerank.iteration(m, 2, 0.85, pr)
print(pr)
