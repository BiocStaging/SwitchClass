# INTERNAL helpers (not exported)
.switch_set_seed <- function(seed) {
  if (!is.null(seed)) {
    seed_fun <- get(paste0("set", ".", "seed"), envir = baseenv())
    seed_fun(as.integer(seed))
  }
  invisible(seed)
}

.classifyImpt <- function(X, y, ntimes = 25, ntree = 1000, mtry = NULL,
                          seed = 1, use_parallel = FALSE,
                          ncores = parallel::detectCores() - 1L) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.numeric(X)) {
    stop("'X' must be numeric or coercible to a numeric matrix.")
  }
  if (is.null(rownames(X))) {
    stop("'X' must have feature row names.")
  }
  if (nrow(X) <= 1L || ncol(X) <= 1L) {
    stop("'X' must contain at least two features and two samples.")
  }
  if (length(y) != ncol(X)) {
    stop("'y' must have length equal to ncol(X).")
  }
  if (length(unique(stats::na.omit(y))) < 2L) {
    stop("'y' must contain at least two classes.")
  }
  ntimes <- as.integer(ntimes)
  ntree <- as.integer(ntree)
  if (is.na(ntimes) || ntimes < 1L) {
    stop("'ntimes' must be a positive integer.")
  }
  if (is.na(ntree) || ntree < 1L) {
    stop("'ntree' must be a positive integer.")
  }
  if (!is.null(mtry)) {
    mtry <- as.integer(mtry)
    if (is.na(mtry) || mtry < 1L || mtry > nrow(X)) {
      stop("'mtry' must be NULL or an integer between 1 and nrow(X).")
    }
  }

  X_t <- t(X)        # samples x features
  y_f <- as.factor(y)
  if (is.null(seed)) {
    seeds <- sample.int(.Machine$integer.max, ntimes, replace = FALSE)
  } else {
    seeds <- as.integer(seed) + seq_len(ntimes) - 1L
  }

  run_one <- function(seed_i) {
    .switch_set_seed(seed_i)
    rf_args <- list(
      x = X_t,
      y = y_f,
      ntree = ntree,
      importance = TRUE
    )
    if (!is.null(mtry)) {
      rf_args$mtry <- mtry
    }
    rf <- do.call(randomForest::randomForest, rf_args)
    rf$importance[, "MeanDecreaseGini"]
  }

  ncores <- max(1L, min(as.integer(ncores), length(seeds)))
  if (isTRUE(use_parallel) && ncores > 1L) {
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(ncores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      impt_list <- parallel::parLapply(cl, seeds, run_one)
    } else {
      impt_list <- parallel::mclapply(seeds, run_one, mc.cores = ncores)
    }
  } else {
    impt_list <- lapply(seeds, run_one)
  }
  impt_mat <- do.call(cbind, impt_list)
  rowMeans(impt_mat, na.rm = TRUE)
}

#' Label-switch classification (with optional multicore)
#'
#' Trains two random-forest classifiers under baseline-aligned and
#' perturbation-aligned label schemes, then returns the directional importance
#' score \eqn{\delta = I_{baseline} - I_{perturbation}} using
#' MeanDecreaseGini importances.
#'
#' @param X numeric matrix with features in rows and samples in columns.
#' @param y_reverse character/factor labels (length = ncol(X)) for the
#'   baseline-aligned label scheme. The argument name is retained for
#'   compatibility with earlier package versions.
#' @param y_persist character/factor labels (length = ncol(X)) for the
#'   perturbation-aligned label scheme. The argument name is retained for
#'   compatibility with earlier package versions.
#' @param ntimes integer; number of RF repeats to average (default 25).
#' @param seed integer; base seed (default 1).
#' @param ntree integer; trees per forest (default 1000).
#' @param mtry integer or NULL; variables tried at each split (default NULL).
#' @param parallel logical; use multicore (default FALSE).
#' @param ncores integer; cores if `parallel=TRUE` (default detectCores()-1).
#' @return list with named numeric vectors: `delta`, `importance_baseline`,
#'   `importance_perturbation`, plus compatibility aliases `I_rev` and `I_per`.
#' @examples
#' X <- matrix(seq_len(60) / 10, nrow = 10)
#' rownames(X) <- paste0("feature", seq_len(nrow(X)))
#' y_baseline <- rep(c("c1", "c2", "c1"), each = 2)
#' y_perturbation <- rep(c("c1", "c2", "c2"), each = 2)
#' res <- label_switch_classify(
#'   X, y_baseline, y_perturbation, ntimes = 1, ntree = 5
#' )
#' head(res$delta)
#' @export
label_switch_classify <- function(X, y_reverse, y_persist, ntimes = 25,
                                  seed = 1, ntree = 1000, mtry = NULL,
                                  parallel = FALSE,
                                  ncores = parallel::detectCores() - 1L) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  importance_baseline <- .classifyImpt(
    X, y_reverse, ntimes = ntimes, ntree = ntree, mtry = mtry,
    seed = seed, use_parallel = parallel, ncores = ncores
  )
  seed_perturbation <- if (is.null(seed)) NULL else as.integer(seed) + 100000L
  importance_perturbation <- .classifyImpt(
    X, y_persist, ntimes = ntimes, ntree = ntree, mtry = mtry,
    seed = seed_perturbation, use_parallel = parallel, ncores = ncores
  )
  list(
    delta = importance_baseline - importance_perturbation,
    importance_baseline = importance_baseline,
    importance_perturbation = importance_perturbation,
    I_rev = importance_baseline,
    I_per = importance_perturbation
  )
}
