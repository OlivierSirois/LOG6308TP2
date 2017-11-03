m = as.matrix(read.table("citeseer.rtable", check.names=F))
library("expm")
library("dplyr")
d <- 0.85
pr <- rep(1, dim(m)[1])

## Cosinus entre un vecteur v et chaque colonne dela matrice m
# cosinus.vm <- function(v,m) {
#   # On on met tous nos valeurs de NA a 0, sinon on va avoir des problemes de calculs avec des matrices sparse
#   m[is.na(m)] <- 0
#   v[is.na(v)] <- 0 ;
#   # On calcule le cosinus entre le vecteur V et les colonnes de la matrice m en utilisant la formule vu en classe
#   (v %*% m)/(sqrt(colSums(m^2)) * sqrt(sum(v^2)))
# }
#
# #Correlation entre la rangee v (v = index) et chaque colonne de la matrice m
# corr.vm <- function(v,m) {
#   # on centre les valeurs de la matrice m en fonction de la moyenne, on la renomme m.centre
#   v.i <- rowMeans(m[,], na.rm=T)
#   # on enleve les NA
#   m[is.na(m)] <- 0
#   m.centre <- m - v.i
#   # on centre le vecteur v en fonction de sa moyenne
#   v.index <- v
#   v.index[is.na(v)] <- 0
#   v.index <- v.index - mean(v, na.rm=T)
#   # on retourne ensuite un vecteur correspondant entre le vecteur v et sa correlation avec chaque rangers de m
#   return( (v.index%*%t(m.centre))/(sqrt(sum(v.index^2) * rowSums(m.centre^2))))
# }

#cette fonction fais la somme de toutes les puissance de la matrice m allant de 0 jusqu'a n
sum.powers.matrix <- function(m, n) {
  powers <- c(1:n)
  res <- Reduce('+', lapply(powers, function(x) m %^% x))
  res[res > 1] <- 1

  return(res)
}

#on enleve toutes les autoreferences dans la matrice referentielle. Ce probleme survient lorsqu'on fait une addition des différentes puissances de matrices.
#ce comportement est a évité car sa donne de l'importance artificielle à notre cible lorsqu'elle est référencé elle même par sa référence.
remove.autoreferences <- function(m) {
  res <- m
  diag(res) <- 0
}

# Cette fonction exécute l'algorithme PageRank jusqu'à ce que c'est valeurs soit stabilisés. On définie une stabilit. lorsque l'erreur moyenne absolue est moins de .0001
page.rank.until.stab <- function(m, d, pr){
  pr.next <- page.rank(m, d, pr)
  while(abs(mean(pr.next-pr)) > 0.0001){
    pr <- pr.next
    pr.next <- page.rank(m,d,pr)
  }
  return (pr)
}

#Cette fonction fait une itérations de l'algorithme page rank. Nous avons ajouter une petite modifications dans le cas ou la somme d'un colonne est de 0 (c'est à dire, une
# article non référencé) on remplace ensuite la valeur qui sera égale a Inf, par un 0
page.rank <- function(m, d, pr){
  denum <- (pr/colSums(m))
  denum[denum == Inf] <- 0
  (pr <- (1-d)/3 + (d * (m %*%denum)))
  #print(pr[450])
  return(pr)
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
# # On calcule notre domaine S comme étant tout les article sont référencés par notre origine. Pour faire cela, on regarde ceux qui on une valeur positive dans notre matrice référentielle
# S <- which(m["422908",]==1)
# # Pour le domaine S prime, nous voulons aussi rajouter les références des références. Pour faire cela, nous faisont que prendre la somme des deux première puissance de la
# # matrice référentielle. C'est a dire, la matrice référentielle elle-même et la deuxième puissance.
# S.prime <- which(sum.powers.matrix(m,2)["422908",]==1)
#
#
# #On calcule le pagerank de tout nos articles (pas très long)
# pr <- page.rank.until.stab(m, d, pr)
# #on place les rankings PageRank de notre domaine dans un vecteur
# S.rankings <- pr[S]
# #pour le domaine S prime
# S.prime.rankings <- pr[S.prime]
# #On crée un dataframe avec nos données
# S.dat <- data.frame(S.rankings, S)
# S.dat$article = rownames(S.dat)
# #même chose pour le domaine S prime
# S.prime.dat <- data.frame(S.prime.rankings, S.prime)
# S.prime.dat$article = rownames(S.prime.dat)
# #on effectue un trie sur nos valeurs, on sort celles qui sont les plus hautes en premier
# S.best <- S.dat %>% select(S.rankings, S, article) %>% arrange(desc(S.rankings, arr.ind=T))
# #même chose pour le domaine S prime
# S.prime.best <- S.prime.dat %>% select(S.prime.rankings, S.prime, article) %>% arrange(desc(S.prime.rankings, arr.ind=T))
#
# ## on calcule nos coefficients de correlation et du cosinus
#
# corr.ratings <- corr.vm(m["422908", ], m[S,])
# cos.ratings <- cosinus.vm(m["422908",], t(m[S,]))
#
# # on remplace les NaN par 0
# corr.ratings[is.nan(corr.ratings)] <- 0
# cos.ratings[is.nan(cos.ratings)] <- 0
#
# # on prends nos etiquettes
# labels.corr <- colnames(corr.ratings)
# labels.cos <- colnames(cos.ratings)
#
# # on cree nos dataframes
# df.corr <- data.frame(corr = as.vector(corr.ratings), article = labels.corr)
# df.cos <- data.frame(cos = as.vector(cos.ratings), article = labels.cos)
#
# #on trie
# df.best.cos <- df.cos %>% select(cos, article) %>% arrange(desc(cos, arr.ind=T))
# df.best.corr <- df.corr %>% select(corr, article) %>% arrange(desc(corr, arr.ind=T))

#print(cl)

# ------------------- Cross-validation: item-item method -----------------------

#Correlation entre la rangee v (v = index) et chaque colonne de la matrice m
corr.vm <- function(v,m) {
  # on centre les valeurs de la matrice m en fonction de la moyenne, on la renomme m.centre
  v.i <- rowMeans(m[,], na.rm=T)
  # on enleve les NA
  m[is.na(m)] <- 0
  m.centre <- m - v.i
  # on centre le vecteur v en fonction de sa moyenne
  v.index <- v
  v.index[is.na(v)] <- 0
  v.index <- v.index - mean(v, na.rm=T)
  # on retourne ensuite un vecteur correspondant entre le vecteur v et sa correlation avec chaque rangers de m
  return( (v.index%*%t(m.centre))/(sqrt(sum(v.index^2) * rowSums(m.centre^2))))
}

cosinus.vm <- function(v,m) {
  # On on met tous nos valeurs de NA a 0, sinon on va avoir des problemes de calculs avec des matrices sparse
  m[is.na(m)] <- 0
  v[is.na(v)] <- 0 ;
  # On calcule le cosinus entre le vecteur V et les colonnes de la matrice m en utilisant la formule vu en classe
  denom <- sqrt(colSums(m^2)) * sqrt(sum(v^2))
  denom[denom == 0] <- 1
  res <- (v %*% m) / denom
  # if (denom == 0)
  # {
  #     if (colSums(m^2) == sum(v^2))
  #     {
  #         res <- 1
  #     }
  #     else
  #     {
  #         res <- 0
  #     }
  # }
  # else
  # {
  #
  # }



  return(res)
}

min.nindex <- function(m, n=11) {
  i <- order(m)
  return(i[1:n])
}

max.nindex <- function(m, n=5) {
    i <- order(m, decreasing=TRUE)
    return(i[1:n])
}

# predict <- function(refs, article) {
#
#   dist.article <- sqrt(colSums(refs[,article] - refs)^2)
#
#   i.distance.article <- min.nindex(dist.article)[-1]
#
#   # print(i.distance.article)
#   print(sapply(i.distance.article, function(x) cosinus.vm(refs[, article], matrix(refs[,x]))))

  # corr.articles <- sapply(i.distance.article, function(x) corr.vm(refs[,article], refs[,x]))
  #
  # print(corr.articles)


  # # Vecteur contenant les distances entre Star Trek et les autres films
  # distance.450 <- sqrt(colSums(ratings[,450] - ratings)^2)
  #
  #
  # # Calcul des 20 voisins les plus proches
  # n.voisins <- 20 + 1
  # votes.communs <- colSums((ratings[,450] * ratings) > 0) # nombre de votes communs
  # #print(which(votes.communs==0))
  # i.distance.450 <- min.nindex(distance.450, n.voisins)
  # # votes.communs[i.distance.450]
  #
  # i.distance.450 <- i.distance.450[i.distance.450 != 450]
  # # Moyenne des votes par film (sans les NA)
  # i.mean.item <- matrix(colMeans(m[], na.rm=TRUE))
  # i.450.mean <- mean(m[,450], na.rm=T)
  #
  # # on sort nos predictions en utilisant la formule apprise dans le cours
  # res <- sapply(unique(1:943), function(x) predict.vote(t(m[,i.distance.450]), i.450.mean, x, rowMeans(t(m[,i.distance.450]), na.rm = T), cosinus.vm(m[,450], m[,i.distance.450])))
  #
  # return(res)
# }



predict.value <- function(refs, user, item)
{
    items.avg <- colMeans(refs, na.rm=T)
    # print(items.avg)
    items.similarity <- cosinus.vm(refs[,item], refs)
    distance.item <- sqrt(colSums(refs[,item] - refs, na.rm=T)^2)
    n.voisins <- 20 + 1
    most.similar <- max.nindex(items.similarity, n.voisins)[-item]

    correction.factor <- 1 / sum(items.similarity[most.similar])
    if (correction.factor == Inf)
    {
        correction.factor <- 1
    }

    value <- items.avg[item] + correction.factor *
             sum(sapply(most.similar, function(x)
                items.similarity[x] * (refs[user, x] - items.avg[x])), na.rm=T)

    return(value)
}

rmse <- function(results, target)
{
    # print(results^2 - target^2)
    res <- sqrt(sum(abs(results^2 - target^2)) / dim(target)[1]^2)

    return(res)
}

# Proportion de références de test (ici : 10%)
cross.validation.factor <- 0.05

# Nombre d'articles dans la base
nb.articles <- dim(m)[1]

# Sélection aléatoire des indices des références de la base de test
test.refs.row.indices <- sample(1:nb.articles, round(cross.validation.factor * nb.articles))
test.refs.col.indices <- sample(1:nb.articles, round(cross.validation.factor * nb.articles))

# Séparation de la base d'entraînement et de la base de test
m.test <- m[test.refs.row.indices, test.refs.col.indices]
m.training <- m
m.training[test.refs.row.indices, test.refs.col.indices] <- NA

# print(predict(as.matrix(m.training), test.refs.col.indices[1]))
# predict.value(as.matrix(m.training), test.refs.row.indices[1], test.refs.col.indices[1])
results <-sapply(test.refs.row.indices, function(x)
            sapply(test.refs.col.indices, function(y)
                predict.value(as.matrix(m.training), x, y)))

print(rmse(as.matrix(results), as.matrix(m.test)))

# print(m.test)

# print("Test: ")
# m.test
#
# print("Training:")
# m.training



# cos.ratings <- sapply(test.articles.indices, function(x) cosinus.vm(m[x,], m[, -x]))
# cos.ratings
