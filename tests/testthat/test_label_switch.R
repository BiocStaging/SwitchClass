test_that("label_switch_classify returns expected names and sizes", {
  X <- matrix(
    sin(seq_len(200 * 12) / 7) + cos(seq_len(200 * 12) / 13),
    200, 12,
    dimnames = list(paste0("f", 1:200), paste0("s", 1:12))
  )
  y2 <- rep(c("c1","c2","c1"), each = 4)
  y3 <- rep(c("c1","c2","c2"), each = 4)
  res <- label_switch_classify(
    X, y_reverse = y2, y_persist = y3,
    ntimes = 2, ntree = 20, seed = 10
  )
  expect_true(all(c(
    "delta", "importance_baseline", "importance_perturbation", "I_rev", "I_per"
  ) %in% names(res)))
  expect_equal(length(res$delta), nrow(X))
  expect_equal(names(res$delta), rownames(X))
})

test_that("label_switch_classify is deterministic with a fixed seed", {
  X <- matrix(
    seq_len(80 * 9) / 100,
    80, 9,
    dimnames = list(paste0("f", 1:80), paste0("s", 1:9))
  )
  y_baseline <- rep(c("c1", "c2", "c1"), each = 3)
  y_perturbation <- rep(c("c1", "c2", "c2"), each = 3)
  res1 <- label_switch_classify(
    X, y_baseline, y_perturbation,
    ntimes = 2, ntree = 20, seed = 99
  )
  res2 <- label_switch_classify(
    X, y_baseline, y_perturbation,
    ntimes = 2, ntree = 20, seed = 99
  )
  expect_equal(res1$delta, res2$delta)
})

test_that("assign_quadrants applies signed delta and fold-change rules", {
  delta <- c(a = 0.2, b = -0.2, c = -0.2, d = 0.2, e = 0)
  fc <- c(a = 1, b = 1, c = -1, d = -1, e = 1)
  expect_equal(
    as.character(assign_quadrants(delta, fc, delta_thresh = 0.15)),
    c("Q1", "Q2", "Q3", "Q4", "Q0")
  )
})
