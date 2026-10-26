# 32_cart_surrogate.R - Comparison: does a surrogate CART read the forest as
# well as the explainable ensemble tree?
#
# The permutation test of 04_e2tree.R says that the structural fidelity of the
# explanatory tree is significant, but the null distribution has a tiny spread
# (sd about 4e-5), so significance alone says little. The real question is
# whether an ordinary regression tree, grown on the response and unaware of
# the forest, would give the same reading. If so, the apparatus is not needed.
#
# Same training set, same dissimilarity matrix of the forest D, same nLoI
# metric, same tree size. Only the induction criterion changes:
#   - e2tree : partition grown to reproduce D;
#   - CART   : partition grown to predict y.
#
# eValidation() weights co-membership with the internal 'prob' of the e2tree
# nodes, which an rpart tree does not have. To compare like with like, the
# binary convention is used for both: Ps[i,j] = 1 if i and j fall in the same
# leaf, 0 otherwise. The partitions are thus compared, which is exactly what
# distinguishes the two methods. The native value of e2tree (0.083, with the
# weights) is in T6 and is not comparable with these two.
#
# Output: T47_cart_surrogate.csv
source("R/00_setup.R")
suppressPackageStartupMessages({library(e2tree); library(rpart)})

su <- readRDS(file.path(PATHS$models, "e2tree_setup.rds"))
e2 <- readRDS(file.path(PATHS$models, "e2tree_fit.rds"))
S  <- su$S; fml <- su$formula
D  <- e2$D; df <- e2$tree$tree
n  <- nrow(S)

# ---- nLoI under the binary convention, replicating eValidation --------------
# eValidation: ord <- hclust(as.dist(D))$order ; O <- sqrt(1 - D[ord,ord]) ;
#              O_hat <- sqrt(Ps[ord,ord]) ; loi(O, O_hat)
ord <- hclust(as.dist(D))$order
O   <- sqrt(pmax(1 - D[ord, ord], 0))

nloi_from_membership <- function(memb) {
  Ps <- matrix(0, n, n)
  for (g in unique(memb)) { idx <- which(memb == g); Ps[idx, idx] <- 1 }
  diag(Ps) <- 1
  loi(O, sqrt(pmin(pmax(Ps[ord, ord], 0), 1)))$nloi
}

# leaf membership of the e2tree tree
memb_e2 <- integer(n)
tn <- df$node[df$terminal]
for (k in tn) memb_e2[unlist(df$obs[df$node == k])] <- k
n_leaves_e2 <- length(tn)
stopifnot(all(memb_e2 > 0))

# ---- CART pruned to the same size -------------------------------------------
set.seed(SEED)
big <- rpart(fml, data = S, method = "anova",
             control = rpart.control(cp = 0, minbucket = 20, xval = 0))
tab <- as.data.frame(big$cptable); tab$leaves <- tab$nsplit + 1
ok  <- tab[tab$leaves <= n_leaves_e2, ]
cart <- prune(big, cp = ok$CP[which.max(ok$leaves)])
memb_cart <- as.integer(factor(cart$where))
n_leaves_cart <- length(unique(memb_cart))

message(sprintf("e2tree: %d leaves | pruned CART: %d leaves",
                n_leaves_e2, n_leaves_cart))

# ---- Main metric: correlation between partition and proximity ----------------
# Under the binary convention nLoI depends on the distribution of leaf sizes:
# a random partition into 23 equal groups gets 0.058, better than the e2tree
# tree (0.083), only because the tree has very different leaves (from 46 to
# 729 units). The comparison therefore uses the point-biserial correlation
# between the forest proximity, 1 - D, and the indicator of belonging to the
# same leaf: it is symmetric, does not depend on size and measures whether the
# partition reproduces the proximity structure of the forest.
prox_corr <- function(memb) {
  P <- 1 - D                       # proximity of the ensemble
  sx <- sum(P) - sum(diag(P)); sxx <- sum(P^2) - sum(diag(P)^2)
  m  <- n * (n - 1)                # ordered off-diagonal pairs
  sy <- 0; sxy <- 0
  for (g in unique(memb)) {
    idx <- which(memb == g); k <- length(idx)
    if (k < 2) next
    blk <- P[idx, idx]
    sy  <- sy  + k * (k - 1)
    sxy <- sxy + sum(blk) - sum(diag(blk))
  }
  syy <- sy                        # binary indicator: y^2 = y
  num <- sxy / m - (sx / m) * (sy / m)
  den <- sqrt((sxx / m - (sx / m)^2) * (syy / m - (sy / m)^2))
  num / den
}

nloi_e2   <- nloi_from_membership(memb_e2)
nloi_cart <- nloi_from_membership(memb_cart)
cor_e2    <- prox_corr(memb_e2)
cor_cart  <- prox_corr(memb_cart)

# ---- Reference: random partition of the same size ------------------------------
# It tells how much of the fidelity comes from having 23 groups and how much
# from having chosen them well. The sizes are those of the e2tree tree, so the
# control gains nothing from being more balanced.
set.seed(SEED)
sizes <- as.integer(table(memb_e2))
rand_rep <- replicate(10, {
  mm <- integer(n); pool <- sample(n)
  s <- 1; for (j in seq_along(sizes)) {
    mm[pool[s:(s + sizes[j] - 1)]] <- j; s <- s + sizes[j] }
  c(nloi = nloi_from_membership(mm), cor = prox_corr(mm)) })
nloi_rand <- mean(rand_rep["nloi", ]); cor_rand <- mean(rand_rep["cor", ])

r2 <- function(o, p) 1 - sum((o - p)^2) / sum((o - mean(o))^2)
out <- data.frame(
  summariser = c("e2tree (induced on the ensemble structure)",
                 "CART surrogate (grown on the response)",
                 "random partition, same leaf sizes"),
  leaves    = c(n_leaves_e2, n_leaves_cart, n_leaves_e2),
  prox_corr = round(c(cor_e2, cor_cart, cor_rand), 4),
  nLoI_bin  = round(c(nloi_e2, nloi_cart, nloi_rand), 4))
out$in_sample_R2 <- c(NA, round(r2(S[[TARGET]], predict(cart, S)), 3), NA)

write.csv(out, file.path(PATHS$tables, "T47_cart_surrogate.csv"),
          row.names = FALSE)
cat("\n#### e2tree vs surrogate CART (binary convention, same D) ####\n")
print(out, row.names = FALSE)
message("32_cart_surrogate OK -> T47_cart_surrogate.csv")
