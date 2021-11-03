#' @useDynLib mtest
#' @importFrom Rcpp sourceCpp
NULL

#' p-val calculation to compare mulinomial trials. This function is R only,
#' use `m.test` for a faster alternative.
#'
#'
#'
#' @param experiments List of multinomial trials. Each
#'   trial is codified as a vector with the number of results in each category.
#'   The number and order of categories must be the same
#'   in every experiment.
#'   In a
#'   \emph{binomial} experiment, this would be \code{c(number_of_successes,
#'   number_of_failures)}. The function also accepts a matrix (one column per
#'   experiment).
#'
#' @return p-value according to the null hypothesis that the probability of each
#'   outcome is the same in every experiment (two-sided)
#' @export
#'
#' @examples
#' bernPval(list(c(8, 2), c(4, 7)))
bernPval <- function(experiments=list()) {
  mex <- .toMatrix(experiments)
  r <- bernDist(experiments)
  result <- .getCVal(r, mex)
  return(result)
}

#' Get the probabilities of each possible result in a set of multinomial trials
#'
#' This function is called by \link{bernPval}, and therefore is not usually
#' called independently. It is exported so that interested users can check the
#' distribution. To also compute the cumulative probabilities, the result must
#' be sent to \link{calCumul}. A given result can be retrieved with
#' \link{probeDist}
#'
#' @inheritParams bernPval
#'
#' @return List containing each possible result with the same number of
#'   experiments, tries and outcomes as the input. Each element in the list
#'   contains a matricial description of the result (\code{desc}), the total
#'   number of times each outcome was measured (\code{sdesc}) and the absolute
#'   probability of the result (\code{val}).
#'
#' @export
#'
#' @examples
#' r <- bernDist(list(c(8, 2), c(4, 7)))
bernDist <- function(experiments=list()){
  mex <- .toMatrix(experiments)
  df <- ncol(mex) - 1
  startDesc <- .zero(mex)
  n <- (sum(unlist(experiments)))
  denom <- prod(unlist(Map({function(x) n + x}, c(1:df))))
  v <- factorial(df) / denom
  sumDesc <- colSums(startDesc)
  l <- .bdist(list(val=v, desc=startDesc, sdesc=sumDesc))
  r <- l[order(sapply(l, function(x) x$val))]
  return(r)
}


.getCVal <- function(l, m){
  cummVal <- 1
  subs <- 0
  lastVal <- l[[length(l)]]$val
  foundIt <- F
  i <- length(l)
  while (!foundIt && i > 0){
    cVal <- l[[i]]$val
    if (!isTRUE(all.equal(cVal, lastVal))){
      lastVal = cVal
      cummVal <- cummVal - subs
      subs <- 0
    }
    l[[i]]$cval <- cummVal
    subs <- subs + cVal
    if (.mcomp(m, l[[i]]$desc)){
      foundIt <- T
    }
    i <- i - 1
  }
  return(cummVal)
}



#' Check a given result in a distribution
#'
#' Once a distribution of results has been calculated with \link{bernDist} and
#' possibly \link{calCumul}, this function returns information on a particular
#' result
#'
#' @inheritParams calCumul
#' @inheritParams bernDist
#'
#' @return Data calculated for \code{experiment} in \code{l}. This includes the
#'   matricial description of the result (\code{desc}), the total number of
#'   times each outcome was measured (\code{sdesc}) and the absolute probability
#'   of the result (\code{val}). If \link{calCumul} has been run, it will also
#'   show the cumulative probability (\code{cval}). If the \code{experiment} is
#'   not found in \code{l}, the function returns \code{NULL}.
#' @export
#'
#' @examples
#' r <- bernDist(list(c(8, 2), c(4, 7)))
#' s <- calCumul(r)
#' t <- probeDist(s, list(c(8, 2), c(4, 7)))
probeDist <- function(l, experiments){
  m <- .toMatrix(experiments)
  foundIt <- F
  i <- length(l)
  result <- NULL
  while (!foundIt && i > 0){
    if (.mcomp(m, l[[i]]$desc)){
      foundIt <- T
      result <- l[[i]]
    }
    i <- i - 1
  }
  return(result)
}

.mcomp <- function(m1, m2){
  result <- F
  if (ncol(m1) != ncol(m2) || nrow(m1) != nrow(m2)){
    return(result)
  }
  m <- ncol(m1)
  n <- nrow(m1)
  for (i in 1:n){
    for (j in 1:m){
      if (m1[i,j] != m2[i,j]){
        return(result)
      }
    }
  }
  result <- T
  return(result)
}

#' Calculate the cumulative probability of each result
#'
#' The cumulative probability is the sum of the probability of a given result
#' plus the probabilities of every other equally or less likely result. For a
#' given result, this is the \code{p-val} returned by \link{bernPval}.
#'
#'
#' @param l List containing each possible result with the same number of
#'   experiments, tries and outcomes as the original input. Usually, this will
#'   be the output of \link{bernDist}.
#'
#' @return List with the same elements and fields as the input, plus a
#'   \code{cval} field with the cumulative probability.
#'
#' @export
#'
#' @examples
#' r <- bernDist(list(c(8, 2), c(4, 7)))
#' s <- calCumul(r)
calCumul <- function(l){
  cummVal <- 1
  subs <- 0
  lastVal <- l[[length(l)]]$val
  for (i in length(l):1){
    cVal <- l[[i]]$val
    if (!isTRUE(all.equal(cVal, lastVal))){
      lastVal = cVal
      cummVal <- cummVal - subs
      subs <- 0
    }
    l[[i]]$cval <- cummVal
    subs <- subs + cVal
  }
  return(l)
}

.bdist <- function(l){
  result <- list()
  tt <- .slice(l)
  h <- tt[[1]]
  t <- tt[[2]]
  sd <- l$sdesc
  if (is.null(t)){
    return(lorder(l))
  }
  lt <- .bdist(t)
  result <- .unslice(h, lt)
  return(result)
}

.nOutcomes <- function(experiments){
  nd <- max(unlist((Map({function(l) length(l)}, experiments))))
  return(nd)
}

.toMatrix <- function(experiments){
  r <- experiments
  if (is.matrix(r)){
    r <- t(r)
  }
  else{
    nd <- .nOutcomes(experiments)
    r <- matrix(unlist(experiments, recursive = T), ncol = nd, byrow = T)
  }
  return(r)
}


.zero <- function(experiments){
  result <- matrix(unlist(Map({function(x) rep(0, length(x))}, experiments)), ncol = ncol(experiments), byrow = T)
  i <- 1
  for (xp in 1:nrow(experiments)){
    ns <- sum(unlist(experiments[xp,]))
    result[i,1] <- ns
    i = i + 1
  }
  return(result);
}


.lappend <- function(l, v){
  result <- c(l, list(v))
  return(result)
}

.lhead <- function(l, n=1){
  result <- NULL
  if (length(l) > (n-1)){
    result <- l[1:n]
  }
  return(result)
}

.ltail <- function(l, nt=1){
  result <- NULL
  if (length(l) > nt){
    result <- l[(nt+1):length(l)]
  }
  return(result)
}

.mhead <- function(m, n=1){
  result <- NULL
  if (!is.null(nrow(m)) && nrow(m) > (n - 1)){
    result <- m[1:n,]
  }
  return(result)
}

.mtail <- function(m, n=1){
  result <- NULL
  if (!is.null(nrow(m)) && nrow(m) > n){
    result <- m[(n+1):nrow(m),]
  }
  return(result)
}

.head <- function(l, n=1){
  result <- NULL
  if (!is.null(l) && !is.null(nrow(l))){
    result <- .mhead(l, n)
  }
  else{
    result <- .lhead(l, n)
  }
  return(result)
}

.tail <- function(l, n=1){
  result <- NULL
  if (!is.null(l) && !is.null(nrow(l))){
    result <- .mtail(l, n)
  }
  else{
    result <- .ltail(l, n)
  }
  return(result)
}

.dice <- function(l, n=1){
  tl <- NULL
  if (!is.null(l) && !is.null(l$desc)){
    d <- .ltail(as.vector(l$desc, mode="integer"), n)
    if (!is.null(d)){
      tl <- list(val=l$val, desc=d, sdesc=.ltail(l$sdesc, n))
    }
  }
  result <- list(list(val=1, desc=.lhead(as.vector(l$desc, mode="integer"), n), sdesc=.lhead(l$sdesc, n)),
                 tl)
  return(result)
}

.slice <- function(l, n=1){
  tl <- NULL
  if (!is.null(l) && !is.null(l$desc)){
    d <- .mtail(l$desc, n)
    if (!is.null(d)){
      tl <- list(val=l$val, desc=d, sdesc=l$sdesc)
    }
  }
  result <- list(list(val=1, desc=.mhead(l$desc, n), sdesc=l$sdesc),
                 tl)
  return(result)
}

.undice <- function(l1, l2){
  result <- list(val=l2$val, desc=append(l1$desc, l2$desc), sdesc=append(l1$sdesc, l2$sdesc))
  return(result)
}

.unslice <- function(h, l2){
  result <- list()
  for (t in l2){
    th <- h
    th$val <- t$val
    th$sdesc <- t$sdesc
    l1 <- lorder(th)
    for (l in l1){
      d <- rbind(l$desc, t$desc)
      #d <- .redim(d)
      r <- list(val=l$val, desc=d, sdesc=l$sdesc)
      result <- .lappend(result, r)
    }
  }
  return(result)
}

.redim <- function(m){
  n <- ncol(m)
  o <- nrow(m)
  dimnames(m) <- list(1:o, 1:n)
  return(m)
}

ddesc <- function(l){
  for (n in l){
    print(n$desc)
  }
}

lorder <- function(v){
  lv <- v
  tt <- .dice(lv)
  h <- tt[[1]]
  t <- tt[[2]]
  result <- list(lv)
  if (is.null(t)){
    return(result)
  }
  vh <- h$desc[1]
  while (vh > 0){
    t$val <- t$val * h$desc*(1+t$sdesc[1]) / ((1 + t$desc[1])*h$sdesc)
    t$desc[1] <- t$desc[1] + 1
    t$sdesc[1] <- t$sdesc[1] + 1
    h$sdesc <- h$sdesc - 1
    n <- lorder(t)
    for (l in n){
      h$desc[1] <- vh - 1
      r <- .undice(h, l)
      result <- .lappend(result, r)
    }
    vh <- vh - 1
  }
  return(result)
}


#' @export
bpval <- function(experiments){
  mex <- .toMatrix(experiments)
  df <- ncol(mex) - 1
  st <- .zero(mex)
  n <- sum(mex)
  denom <- prod(unlist(Map({function(x) n + x}, c(1:df))))
  v <- factorial(df) / denom
  su <- colSums(st)
  for (i in 1:nrow(mex)){
    for (j in 2:ncol(mex)){
      while (st[i, j] < mex[i, j]){
        v <- v * st[i, 1] * (1 + su[j]) / ((1 + st[i, j]) * su[1])
        st[i, 1] <- st[i, 1] - 1
        st[i, j] <- st[i, j] + 1
        su[1] <- su[1] - 1
        su[j] <- su[j] + 1
      }
    }
  }
  return(v)
}

#' @export
lbpval <- function(experiments){
  mex <- .toMatrix(experiments)
  df <- ncol(mex) - 1
  st <- .zero(mex)
  n <- sum(mex)
  denom <- sum(log((unlist(Map({function(x) n + x}, c(1:df))))))
  v <- sum(log(1:df)) - denom
  su <- colSums(st)
  for (i in 1:nrow(mex)){
    for (j in 2:ncol(mex)){
      while (st[i, j] < mex[i, j]){
        v <- v + log(st[i, 1]) + log(1 + su[j]) - log((1 + st[i, j])) -log(su[1])
        st[i, 1] <- st[i, 1] - 1
        st[i, j] <- st[i, j] + 1
        su[1] <- su[1] - 1
        su[j] <- su[j] + 1
      }
    }
  }
  return(v)
}


binomial.coef <- function(n, k){
  result <- 1
  for (i in 1:k){
    result <- result * (n - k + i) / i
  }
  return(result)
}

.tzero <- function(experiments){
  result <- matrix(unlist(Map({function(x) rep(0, length(x))}, experiments)), ncol = ncol(experiments), byrow = T)
  n1 <- sum(unlist(experiments[1,]))
  n2 <- sum(unlist(experiments[2,]))
  result[1, 2] <- n1
  result[2, 1] <- n2
  return(result);
}

#' @export
tbpval <- function(experiments){
  mex <- .toMatrix(experiments)
  st <- .tzero(mex)
  n <- sum(mex)
  denom <- (n + 1) * (n + 2) * binomial.coef(n, st[1, 2])
  v <- 1 / denom
  su <- colSums(st)
  ns <- rowSums(st)
  vo <- 1 / ((n + 2) * binomial.coef(n+1, st[1, 2]+1))
  vo <- vo * (st[1, 2]+1) * (1 + su[1]) / ((1 + st[1, 1]) * (su[2]+1))
  while (st[1, 1] < mex[1, 1]){
    v <- v + vo / (ns[1] + 1)
    st[1, 1] <- st[1, 1] + 1
    st[1, 2] <- st[1, 2] - 1
    su[1] <- su[1] + 1
    su[2] <- su[2] - 1
    vo <- vo * (st[1, 2]+1) * (1 + su[1]) / ((1 + st[1, 1]) * (su[2]+1))
  }
  vo <- vo * ((ns[2] + 1) * (st[1, 1] + 1) * (su[2] + 1)) / ((ns[1] + 1) * (st[2, 2] + 1) * (su[1] + 1))
  while (st[2, 2] < mex[2, 2]){
    v <- v + vo / (ns[2] + 1)
    st[2, 1] <- st[2, 1] - 1
    st[2, 2] <- st[2, 2] + 1
    su[1] <- su[1] - 1
    su[2] <- su[2] + 1
    vo <- vo * (st[2, 1]+1) * (1 + su[2]) / ((1 + st[2, 2]) * (su[1]+1))
  }
  return(v)
}

#' One-sided p-val calculation to compare binomial trials
#'
#' This function only takes two experiments with two outcomes, usually called
#' \emph{sucess} and \emph{failure}. The null hypothesis states that the
#' probability of success in the first experiment is higher than the
#' probability of success in the second experiment. Therefore, the order
#' of the experiments is important.
#'
#' @inheritParams bernPval
#'
#' @return p-value for the one-sided null hypothesis.
#' @export
#'
#' @examples
#' # H0: p(success in exp1) > p(success in exp2)
#' os.m.test(list(c(10, 8), c(10, 3)))
os.m.test <- function(experiments){
  mex <- .toMatrix(experiments)
  cutoff <- tbpval(experiments)
  if (cutoff == 0){
    message(paste('The probability of the result is too low. Probable overflow, ',
            'the result is unreliable.', sep = ''))
  }
  st <- .tzero(mex)
  n <- sum(mex)
  denom <- (n + 1) * (n + 2) * binomial.coef(n, st[1, 2])
  v <- 1 / denom
  su <- colSums(st)
  ns <- rowSums(st)
  vo <- 1 / ((n + 2) * binomial.coef(n+1, st[1, 2]+1))
  vo <- vo * (st[1, 2]+1) * (1 + su[1]) / ((1 + st[1, 1]) * (su[2]+1))
  result <- 2 * tmtest(st, v, cutoff, ns, su, vo)
  return(result)
}

#' @export
lm.test <- function(experiments){
  mex <- .toMatrix(experiments)
  cutoff <- lbpval(experiments)
  df <- ncol(mex) - 1
  nr <- nrow(mex)
  st <- .zero(mex)
  su <- colSums(st)
  n <- sum(mex)
  denom <- sum(log(unlist(Map({function(x) n + x}, c(1:df)))))
  v <- sum(log(1:df)) - denom
  r <- 0
  if (v <= cutoff){
    r <- v
  }
  result <- lmtest(list(val=v, desc=st, sdesc=su, r=0), cutoff, 0, 0)
  r <- sumlog(r, result$r)
  if (nr > 1){
    for (rw in 2:nr){
      t <- lmtest(list(val=v, desc=st, sdesc=su, r=0), cutoff,
                 rw - 1, 0)
      r <- sumlog(r, t$r)
    }
  }
  return(r)
}

#' p-val calculation to compare multinomial trials
#'
#' @inheritParams bernPval
#'
#' @return p-value for the null hypothesis (all underlying outcome probabilities
#' are the same in every experiment)
#' @export
#'
#' @examples
#' m.test(list(c(8, 2), c(4, 7)))
m.test <- function(experiments){
  mex <- .toMatrix(experiments)
  cutoff <- bpval(experiments)
  if (cutoff == 0){
    #message(paste('The probability of the result is too low. Probable overflow, ',
    #              'the result is unreliable.', sep = ''))
    return(lm.test(experiments))
  }
  df <- ncol(mex) - 1
  nr <- nrow(mex)
  st <- .zero(mex)
  su <- colSums(st)
  n <- sum(mex)
  denom <- prod(unlist(Map({function(x) n + x}, c(1:df))))
  v <- factorial(df) / denom
  r <- 0
  if (v <= cutoff){
    r <- v
  }
  result <- mtest(list(val=v, desc=st, sdesc=su, r=0), cutoff, 0, 0)
  r <- r + result$r
  if (nr > 1){
    for (rw in 2:nr){
      t <- mtest(list(val=v, desc=st, sdesc=su, r=0), cutoff,
                 rw - 1, 0)
      r <- r + t$r
    }
  }
  return(r)
}

modrunif <- function(nc, low, high, dround = 0){
  result <- runif(nc, low, high);
  if (dround > 0){
    result <- round(result, dround)
  }
  return(result);
}

#' @export
.mc <- function(n, nc, nr, algo=1, NRUNIF=0){
  result <- matrix(nrow = n*nr, ncol=nc)
  if (algo == 1){
    for (r in 1:n){
      w <- modrunif(nc, 0, 1, NRUNIF)
      w <- w / sum(w)
      for (i in 1:nr){
        result[(r-1) * nr + i,] <- w
      }
    }
  }
  if (algo == 2){
    df <- nc - 1
    for (i in 1:n){
      rord <- sample(1:nc, nc, replace = F)
      fst <- rord[1:df]
      lst <- rord[length(rord)]
      t <- 1
      for (d in fst){
        tt <- as.numeric(modrunif(1, 0, t, NRUNIF))
        for (j in 1:nr){
          result[nr * (i-1) + j, d] <- tt
        }
        t <- t - tt
      }
      for (j in 1:nr){
        result[nr * (i-1) + j, lst] <- t
      }
    }
  }
  if (algo == 3){
    df <- nc - 1
    for (i in 1:n){
      rord <- sample(1:nc, nc, replace = F)
      fst <- rord[1:df]
      lst <- rord[length(rord)]
      t <- 1
      for (d in fst){
        tt <- as.numeric(modrunif(1, 0, t, NRUNIF))
        for (j in 1:nr){
          result[nr * (i-1) + j, d] <- min(tt, t-tt)
        }
        t <- max(tt, t - tt)
      }
      for (j in 1:nr){
        result[nr * (i-1) + j, lst] <- t
      }
    }
  }
  if (algo == 4){
    for (i in 1:n){
      rord <- c(sort(rep(modrunif(nc-1, 0, 1, NRUNIF)), decreasing = F), 1)
      r2 <- rord
      for (j in 2:length(r2)){
        r2[j] <- rord[j] - rord[j-1]
      }
      for (j in 1:nr){
        result[nr * (i-1) + j, ] <- r2
      }
    }
  }
  if (algo == 5){
    for (i in 1:n){
      rord <- c(sort(rep(modrunif(nr, 0, 1, NRUNIF)), decreasing = T), 1)
      r2 <- rord
      for (j in 1:nr){
        result[nr * (i-1) + j, ] <- c(r2[j], 1-r2[j])
      }
    }
  }
  if (algo == 6){
    for (i in 1:n){
      rord <- c(rep(modrunif(nr, 0, 1, NRUNIF)), 1)
      r2 <- rord
      for (j in 1:nr){
        result[nr * (i-1) + j, ] <- c(r2[j], 1-r2[j])
      }
    }
  }
  return(result)
}


#' Monte Carlo simulation of multinomial process
#'
#' @inheritParams bernPval
#' @param n Number of simulations to run
#' @param algo Algorithm for generating the pseudorandom matrix. The two
#' most relevant algorithms are `4` (two-sided) and `5` (one-sided).
#' @param NRUNIF number of decimal places in the pseudo-random numbers.
#' Only for debugging purposes.
#'
#' @return List of results, as described in \link{bernDist}
#' @section Details:
#' Algorithms `1-4` assume that the probability of each outcome is the same
#' in every experiment, and therefore the alternative hypothesis
#' is that at least one of the probabilities is different (either higer or
#' lower) in at least one experiment. Algorithm `5` only works on `2x2`
#' tables and assumes that the probability of success in the first experiment
#' is higher than in the second experiment. Each algorithm divides the interval
#' `[0, 1]` in as many sub-intervals as outcomes (`nc`), with the first
#' sub-interval `[0, p_1]` and the last `(p_nc-1, 1]`:
#' 1. Generate `nc - 1` pseudo-random numbers and divide by their sum. The
#' resulting values are biased towards `0.5`.
#' 2. Generate a pseudo-random number `p1` between 0 and 1. Then recursively
#' generate numbers between `1-sum(p1, p2, ...)` until `nc - 1` numbers have
#' been provided. Assign every sub-interval to a pseudo-randomly chosen
#' outcome. The resulting values are biased towards `0` if `nc > 2`.
#' 3. Like `2`, but the sub-interval that is sub-divided is always the largest
#' remaining. Less biased than `2` if `nc > 2`.
#' 4. Generate `nc - 1` pseudo-random numbers in `[0, 1]`, sort them and
#' use them to divide the interval. This is the default algorithm.
#' 5. Same as `4`, but one-sided and only for `2x2` tables. The sub-intervals
#' are generated twice, and the result with a lower probability of success is
#' assigned to the second experiment.
#' 6. Same as `5`, but results are not sorted.
#' @export
#'
#' @examples
#' mcm(list(c(4, 5), c(6, 7)), n=1000, algo=4)
mcm <- function(experiments=list(), n=10000, algo=4, NRUNIF=0){
  mex <- .toMatrix(experiments)
  startDesc <- .zero(mex)
  r1 <- list()
  rp <- .mc(n, ncol(mex), nrow(mex), algo = algo, NRUNIF = NRUNIF)
  for (i in 1:n){
    mr <- matrix(c(0), ncol = ncol(mex), nrow = nrow(mex))
    nr <- nrow(startDesc)
    for (r in 1:nr){
      w <- nr * (i-1) + r
      t <- startDesc[r, 1]
      for (x in 1:t){
        p <- as.numeric(modrunif(1, 0, 1, NRUNIF))
        d <- 1
        s <- rp[w, 1]
        while (d < ncol(mex) && s < p){
          d <- d + 1
          v <- rp[w, d]
          s <- s + v
        }
        mr[r, d] <- mr[r, d] + 1
      }
    }
    repr <- paste(as.vector(mr), collapse = '_')
    if (is.null(r1[[repr]])){
      r1[[repr]]$n <- 0
      r1[[repr]]$val <- 0
      r1[[repr]]$desc <- mr
    }
    r1[[repr]]$n <- r1[[repr]]$n + 1
    r1[[repr]]$val <- r1[[repr]]$val + 1 / n
  }
  r2 <- r1[order(sapply(r1, function(x) x$val))]
  return(r2)
}



