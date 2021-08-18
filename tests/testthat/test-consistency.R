#Matricial input

l4e <- list(c(1, 0), c(0, 2), c(3, 0), c(1, 3))
tk <- matrix(c(1, 0, 0, 2, 3, 0, 1, 3), ncol = 4, byrow = F)
rtk <- bernPval(tk)

test_that("Example p-val is correct with matricial input", {
  expect_equal(rtk, 0.01890332, tolerance=1e-6)
})


#Floating-point
r <- bernDist(list(c(10,0)))
s <- calCumul(r)


test_that("One trial: Every result has the same probability", {
  for (i in 2:length(s)){
    expect_equal(s[[i - 1]]$val, s[[i]]$val)
  }
})

test_that("One trial: Every result has a cval of 1", {
  for (i in 1:length(s)){
    expect_equal(s[[i]]$cval, 1)
  }
})



# General
r_4_trials <- bernDist(l4e)
s_4_trials <- calCumul(r_4_trials)
pval <- bernPval(l4e)

r_4_outcomes <- bernDist(list(c(1, 0, 2, 0), c(3, 0, 4, 0)))


test_that("Number of elements is correct with four trials", {
  expect_equal(length(r_4_trials), 120)
})

test_that("Values sum 1 with four trials", {
  expect_equal(sum(unlist(Map({function(x) x$val}, r_4_trials))), 1, tolerance=1e-6)
})

test_that("Number of elements is correct with four outcomes", {
  expect_equal(length(r_4_outcomes), 2400)
})

test_that("Values sum 1 with four outcomes", {
  expect_equal(sum(unlist(Map({function(x) x$val}, r_4_outcomes))), 1, tolerance=1e-6)
})

test_that("Populate sets last element cval to 1", {
  expect_equal(s_4_trials[[length(s_4_trials)]]$cval, 1)
})

test_that("Populate respects simmetry in first two elements", {
  expect_equal(s_4_trials[[1]]$cval, s_4_trials[[2]]$cval)
  expect_equal(s_4_trials[[1]]$cval, 2 * s_4_trials[[1]]$val, tolerance=1e-6)
})

test_that("Example p-val is correct", {
  expect_equal(pval, 0.01890332, tolerance=1e-6)
})

test_that("Example one-sided p-val is correct", {
  expect_equal(os.m.test(list(c(10, 30), c(30, 50))), 0.006068082, tolerance=1e-7)
})


test_that("bernPval gives the same result as bernDist plus calCumul", {
  el <- probeDist(s_4_trials, l4e)
  expect_equal(pval, el$cval, tolerance=1e-6)
})

#Incompatibility between testthat and Rcpp
test_that("bernDist gives same value as m.test", {
  el <- probeDist(s_4_trials, l4e)
  expect_equal(m.test(l4e), el$cval, tolerance=1e-6)
})

test_that("Sum of val is lower than or equal to cval", {
  c <- 0
  a <- T
  for (i in 1:length(s_4_trials)){
    c <- c + s_4_trials[[i]]$val
    if (c > (s_4_trials[[i]]$cval + 1e7)){
      a <- F
    }
  }
  expect_true(a)
})


