m = read.table("citeseer.rtable")
d <- 0.85
pr <- rep(1, dim(m)[1])
(pr <- (1-d)/3 + (d * (m %*% (pr/colSums(m)))))