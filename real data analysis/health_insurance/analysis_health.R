# Healthcare PPI: repeat the labeling step 1000 times on the same full dataset.
# Y = insurance status; f = supplied prediction; R = income.
# TPR/FPR use R > alpha. AUC gives half credit to ties.

root <- "results/real_data"

cfg <- list(
  B = 1000,
  seed = 20250903,
  lambda = 0.004,        # Section 6.1 ratio n/N
  shift_n = 1200,        # approximate number labeled under covariate shift
  gamma = 0.25,          # income effect in the labeling model
  alpha = 50000,
  eps = 1e-5,            # diagnostic lower bound; propensities are not clipped
  level = 0.95,
  data = file.path(root, "data"),
  out = file.path(root, "healthcare_ppi_lambda_0.004")
)

# Download and read the data 

is_zip <- function(path) {
  if (!file.exists(path)) return(FALSE)
  con <- file(path, "rb")
  on.exit(close(con))
  identical(readBin(con, "raw", 2), charToRaw("PK"))
}

download <- function(path) {
  if (is_zip(path)) return(invisible(path))
  id <- "1RjWsnq-gMngRFRj22DvezcdCVl2MxAIX"
  url <- paste0("https://drive.usercontent.google.com/download?id=", id,
                "&export=download&confirm=t")
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".npz")
  on.exit(unlink(tmp))
  tryCatch(download.file(url, tmp, mode = "wb"), error = function(e) NULL)
  
  
  if (!is_zip(tmp)) {
    if (!requireNamespace("googledrive", quietly = TRUE)) {
      install.packages("googledrive", repos = "https://cloud.r-project.org")
    }
    googledrive::drive_deauth()
    googledrive::drive_download(googledrive::as_id(id), tmp, overwrite = TRUE)
  }
  if (!is_zip(tmp)) stop("The download is not an NPZ archive.")
  stopifnot(file.copy(tmp, path, overwrite = TRUE))
  invisible(path)
}

read_npy <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  readBin(con, "raw", 6)                     
  version <- as.integer(readBin(con, "raw", 2))
  size <- if (version[1] == 1) 2 else 4
  len <- sum(as.numeric(readBin(con, "raw", size)) * 256^(0:(size - 1)))
  header <- rawToChar(readBin(con, "raw", len))
  field <- function(pattern) {
    regmatches(header, regexec(pattern, header, perl = TRUE))[[1]][2]
  }
  type <- field("['\"]descr['\"]\\s*:\\s*['\"]([^'\"]+)['\"]")
  shape <- field("['\"]shape['\"]\\s*:\\s*\\(([^)]*)\\)")
  dims <- trimws(strsplit(shape, ",", fixed = TRUE)[[1]])
  dims <- as.numeric(dims[nzchar(dims)])
  fortran <- field("['\"]fortran_order['\"]\\s*:\\s*(True|False)") == "True"
  n <- prod(dims)
  if (!type %in% c("<i8", "<f4", "<f8")) stop("Unexpected data type: ", type)
  
  if (type == "<i8") {
    bytes <- readBin(con, "raw", 8 * n)
    stopifnot(length(bytes) == 8 * n)
    bytes <- matrix(as.numeric(bytes), nrow = 8)
    low <- colSums(bytes[1:4, , drop = FALSE] * 256^(0:3))
    high <- colSums(bytes[5:8, , drop = FALSE] * 256^(0:3))
    high <- high - (high >= 2^31) * 2^32
    value <- low + high * 2^32
  } else {
    value <- readBin(con, "double", n, size = as.integer(substr(type, 3, 3)),
                     endian = "little")
    stopifnot(length(value) == n)
  }
  if (length(dims) == 2) {
    value <- matrix(value, dims[1], dims[2], byrow = !fortran)
  }
  value
}

load_data <- function(path) {
  tmp <- tempfile("npy_", tmpdir = dirname(path))
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))
  files <- unzip(path, list = TRUE)$Name
  files <- files[match(c("Y.npy", "Yhat.npy", "X.npy"), basename(files))]
  unzip(path, files = files, exdir = tmp, junkpaths = TRUE)
  a <- lapply(file.path(tmp, basename(files)), read_npy)
  d <- data.frame(Y = as.numeric(a[[1]]), f = as.numeric(a[[2]]), R = a[[3]][, 1])
  stopifnot(all(d$Y %in% 0:1), all(is.finite(d$R)),
            all(is.finite(d$f) & d$f >= 0 & d$f <= 1))
  d
}

# TPR, FPR, AUC and their influence functions 

mid_cdf <- function(R, w) {
  keep <- w > 0
  r <- R[keep]
  w <- w[keep]
  ord <- order(r)
  r <- r[ord]
  w <- w[ord]
  groups <- rle(r)
  ends <- cumsum(groups$lengths)
  cw <- cumsum(w)
  mass <- cw[ends] - c(0, cw[ends[-length(ends)]])
  below <- c(0, cumsum(mass))[findInterval(R, groups$values) + 1]
  equal <- mass[match(R, groups$values)]
  equal[is.na(equal)] <- 0
  (below - 0.5 * equal) / sum(mass)
}


metrics <- function(R, Y, alpha, w = rep(1, length(Y))) {
  H <- as.numeric(R > alpha)
  w1 <- w * Y
  w0 <- w * (1 - Y)
  stopifnot(sum(w1) > 0, sum(w0) > 0)
  m1 <- sum(w1) / sum(w)
  m0 <- sum(w0) / sum(w)
  r1 <- sum(w1 * H) / sum(w)
  r0 <- sum(w0 * H) / sum(w)
  
  tpr <- sum(w1 * H) / sum(w1)
  fpr <- sum(w0 * H) / sum(w0)
  q0 <- mid_cdf(R, w0)
  q1 <- 1 - mid_cdf(R, w1)
  auc <- sum(w1 * q0) / sum(w1)
  
  phi <- cbind(
    TPR = H * Y / m1 - r1 * Y / m1^2,
    FPR = H * (1 - Y) / m0 - r0 * (1 - Y) / m0^2,
    AUC = Y * (q0 - auc) / m1 + (1 - Y) * (q1 - auc) / m0
  )
  list(theta = c(TPR = tpr, FPR = fpr, AUC = auc), phi = phi)
}


col_var <- function(x) {
  v <- apply(x, 2, var)
  replace(v, !is.finite(v) | v <= 0, NA_real_)
}

# No shift: original Section 6.1 PPI 

no_shift <- function(d, C, full, alpha) {
  n <- sum(C)
  N <- nrow(d) - n
  lab <- d[C == 1, ]
  y <- metrics(lab$R, lab$Y, alpha)
  f <- metrics(lab$R, lab$f, alpha)
  vy <- col_var(y$phi)
  vf <- col_var(full$phi)
  cyf <- diag(cov(y$phi, f$phi))
  omega <- cyf / vf
  delta <- full$theta - f$theta
  
  # Section 6.1 uses lambda = n/N.
  v <- (vy - cyf^2 / ((1 + n / N) * vf)) / n
  v[v < 0] <- NA_real_
  list(base = y$theta, ppi = y$theta + omega * delta,
       sd_base = sqrt(vy / n), sd_ppi = sqrt(v),
       omega = omega, rectifier = delta,
       pi_min = NA_real_, pi_max = NA_real_, condition = NA_real_,
       components = list(theta_y_lab = y$theta, theta_f_lab = f$theta,
                         theta_f_full = full$theta, omega = omega,
                         lambda_n_over_N = n / N))
}

#  Shift: estimate pi and adjust the IF using Section 4 

design <- function(R) {
  data.frame(income = (R - mean(R)) / sd(R))
}


adjust <- function(phi, C, p, Z, info) {
  b <- crossprod(Z, phi * (C * (1 - p) / p)) / length(C)
  (Z * (C - p)) %*% solve(info, b)
}

shift <- function(d, C, full, alpha, eps) {
  M <- nrow(d)
  fit <- glm(C ~ income, data = design(d$R), family = binomial())
  p <- as.numeric(fitted(fit))
  if (!fit$converged || any(p < eps | p > 1 - eps)) {
    stop("Check the fitted labeling model and propensity overlap.")
  }
  Z <- model.matrix(fit)
  info <- crossprod(Z, Z * (p * (1 - p))) / M
  w <- C / p
  Y <- numeric(M)                     
  Y[C == 1] <- d$Y[C == 1]           
  y <- metrics(d$R, Y, alpha, w)
  f <- metrics(d$R, d$f, alpha, w)
  
  A <- w * y$phi - adjust(y$phi, C, p, Z, info)
  Lf <- w * f$phi - adjust(f$phi, C, p, Z, info)
  B <- full$phi - Lf
  omega <- -diag(cov(A, B)) / col_var(B)
  delta <- full$theta - f$theta
  ppi_if <- A + sweep(B, 2, omega, "*")
  
  list(base = y$theta, ppi = y$theta + omega * delta,
       sd_base = sqrt(col_var(A) / M), sd_ppi = sqrt(col_var(ppi_if) / M),
       omega = omega, rectifier = delta, pi_hat = p,
       pi_min = min(p), pi_max = max(p), condition = kappa(info),
       beta = coef(summary(fit)),
       components = list(theta_y_lab = y$theta, theta_f_lab = f$theta,
                         theta_f_full = full$theta, omega = omega))
}



result_rows <- function(fit, scenario, ref, n, N, level) {
  methods <- if (scenario == "no_shift") {
    c("Labeled only", "PPI Section 6.1")
  } else c("IPW labeled, estimated pi", "PPI Section 4, estimated pi")
  est <- c(fit$base, fit$ppi)
  se <- c(fit$sd_base, fit$sd_ppi)
  truth <- rep(ref, 2)
  z <- qnorm((1 + level) / 2)
  data.frame(
    scenario = scenario, method = rep(methods, each = 3),
    metric = rep(names(ref), 2), estimate = unname(est),
    all_labeled_estimate = unname(truth), bias = unname(est - truth),
    sd = unname(se), lower = unname(est - z * se), upper = unname(est + z * se),
    omega = c(rep(NA_real_, 3), unname(fit$omega)),
    rectifier = c(rep(NA_real_, 3), unname(fit$rectifier)),
    n_labeled = n, n_unlabeled = N, unlabeled_over_labeled = N / n,
    lambda_n_over_N_in_theory = n / N,
    pi_hat_min = fit$pi_min, pi_hat_max = fit$pi_max,
    propensity_condition_number = fit$condition
  )
}

summarize <- function(rows) {
  groups <- split(rows, interaction(rows$scenario, rows$method, rows$metric,
                                    drop = TRUE))
  avg <- function(x) {
    x <- x[is.finite(x)]
    if (length(x)) mean(x) else NA_real_
  }
  ans <- lapply(groups, function(x) {
    est <- x$estimate[is.finite(x$estimate)]
    covered <- is.finite(x$lower) & is.finite(x$upper) &
      is.finite(x$all_labeled_estimate)
    emp_var <- if (length(est) > 1) var(est) else NA_real_
    
    data.frame(
      scenario = x$scenario[1], method = x$method[1], metric = x$metric[1],
      all_labeled_estimate = x$all_labeled_estimate[1], repetitions = nrow(x),
      mean_estimate = avg(x$estimate), mean_bias = avg(x$bias),
      mean_absolute_error = avg(abs(x$bias)), rmse = sqrt(avg(x$bias^2)),
      mean_analytic_sd = avg(x$sd),
      
      mean_analytic_variance = avg(x$sd^2),
      empirical_variance = emp_var,
      empirical_sd = sqrt(emp_var),
      coverage = if (any(covered)) {
        mean(x$lower[covered] <= x$all_labeled_estimate[covered] &
               x$all_labeled_estimate[covered] <= x$upper[covered])
      } else NA_real_,
      mean_ci_width = avg(x$upper - x$lower),
      mean_omega = avg(x$omega), mean_n_labeled = mean(x$n_labeled),
      mean_n_unlabeled = mean(x$n_unlabeled),
      mean_unlabeled_over_labeled = mean(x$unlabeled_over_labeled),
      valid_estimates = sum(is.finite(x$estimate)),
      valid_analytic_sds = sum(is.finite(x$sd)),
      valid_intervals = sum(covered),
      invalid_estimates = sum(!is.finite(x$estimate)),
      invalid_analytic_sds = sum(!is.finite(x$sd)),
      invalid_intervals = sum(!covered)
    )
  })
  ans <- do.call(rbind, ans)
  ans <- ans[order(ans$scenario, ans$method, ans$metric), ]
  rownames(ans) <- NULL
  ans
}


run <- function(d, cfg, md5 = NA_character_) {
  dir.create(cfg$out, recursive = TRUE, showWarnings = FALSE)
  csv <- function(x, name) write.csv(x, file.path(cfg$out, name), row.names = FALSE)
  M <- nrow(d)
  n <- as.integer(round(cfg$lambda * M / (1 + cfg$lambda)))
  if (n < 2 || n >= M) stop("lambda gives an invalid labeled sample size.")
  N <- M - n
  lambda_actual <- n / N
  
  
  ref_fit <- metrics(d$R, d$Y, cfg$alpha)
  ref <- ref_fit$theta
  ref_table <- data.frame(metric = names(ref), all_labeled_estimate = unname(ref),
                          analytic_sd_all_labeled = sqrt(col_var(ref_fit$phi) / M),
                          n_total = M)
  csv(ref_table, "all_labeled_reference.csv")
  full <- metrics(d$R, d$f, cfg$alpha)   
  
  
  X <- design(d$R)
  p_shift <- cfg$shift_n / M
  g0 <- qlogis(p_shift)
  pi_true <- plogis(g0 + cfg$gamma * X$income)
  if (any(pi_true <= cfg$eps | pi_true >= 1 - cfg$eps)) {
    stop("The labeling design does not satisfy the requested overlap.")
  }
  seeds <- data.frame(replication = seq_len(cfg$B),
                      seed_no_shift = cfg$seed + 2 * (seq_len(cfg$B) - 1),
                      seed_shift = cfg$seed + 2 * (seq_len(cfg$B) - 1) + 1)
  csv(seeds, "replication_seeds.csv")
  settings <- list(config = cfg, dataset_md5 = md5, n_total = M, seeds = seeds,
                   design = list(n = n, N = N,
                                 target_shift_labels = cfg$shift_n,
                                 shift_base_rate = p_shift, intercept = g0,
                                 lambda_requested = cfg$lambda,
                                 lambda_n_over_N = lambda_actual,
                                 unlabeled_over_labeled = N / n,
                                 expected_shift_labels = sum(pi_true),
                                 pi_range = range(pi_true), pi_mean = mean(pi_true)),
                   rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
                   reference = ref_table, session_info = sessionInfo())
  saveRDS(settings, file.path(cfg$out, "analysis_settings.rds"))
  
  message("No shift: n = ", n, ", N = ", N,
          ", requested lambda = ", cfg$lambda,
          ", actual n/N = ", signif(lambda_actual, 7), ".")
  message("Shift: E(number labeled) = ", signif(sum(pi_true), 7),
          " under the simple income-only pi model.")
  
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  results <- vector("list", cfg$B)
  for (b in seq_len(cfg$B)) {
    
    set.seed(seeds$seed_no_shift[b])
    C0 <- integer(M)
    C0[sample.int(M, n, replace = FALSE)] <- 1L
    
    
    set.seed(seeds$seed_shift[b])
    C1 <- rbinom(M, 1, pi_true)
    observed <- d
    observed$Y[C0 == 0] <- NA_real_
    fit0 <- no_shift(observed, C0, full, cfg$alpha)
    observed <- d
    observed$Y[C1 == 0] <- NA_real_
    fit1 <- shift(observed, C1, full, cfg$alpha, cfg$eps)
    rows <- rbind(
      result_rows(fit0, "no_shift", ref, n, N, cfg$level),
      result_rows(fit1, "estimated_pi_shift", ref, sum(C1), M - sum(C1), cfg$level)
    )
    rows <- data.frame(replication = b,
                       labeling_seed = rep(c(seeds$seed_no_shift[b], seeds$seed_shift[b]),
                                           each = 6), rows)
    results[[b]] <- rows
   
    write.table(rows, file.path(cfg$out, "replication_level_results.csv"),
                sep = ",", row.names = FALSE, col.names = b == 1, append = b > 1,
                na = "NA")
    
    if (b == 1) {
      csv(rows, "point_estimates.csv")
      csv(rows[grepl("^PPI", rows$method), ], "ppi_point_estimates.csv")
      csv(data.frame(row_id = seq_len(M), C_no_shift = C0,
                     C_estimated_pi_shift = C1, pi_hat_shift = fit1$pi_hat),
          "replication_1_assignments.csv")
      first <- list(results = rows, no_shift = fit0$components,
                    estimated_pi = fit1$components, propensity_coefficients = fit1$beta)
    }
    if (b == 1 || b %% 10 == 0 || b == cfg$B) {
      saveRDS(list(completed = b, results = results[seq_len(b)], metadata = settings),
              file.path(cfg$out, "checkpoint.rds"))
      message("Completed ", b, " / ", cfg$B)
    }
  }
  
  results <- do.call(rbind, results)
  rownames(results) <- NULL
  summary <- summarize(results)
  csv(summary, "operating_characteristics.csv")
  ans <- list(settings = settings, point_estimates = first$results,
              first_analysis = first, reference = ref_table,
              results = results, summary = summary)
  saveRDS(ans, file.path(cfg$out, "healthcare_ppi_1000_repetitions.rds"))
  print(summary[, c("scenario", "method", "metric", "mean_estimate",
                    "mean_bias", "mean_analytic_sd", "empirical_sd",
                    "coverage")], row.names = FALSE)
  invisible(ans)
}



dir.create(cfg$data, recursive = TRUE, showWarnings = FALSE)
file <- file.path(cfg$data, "census_healthcare.npz")
download(file)
dat <- load_data(file)
analysis <- run(dat, cfg, md5 = unname(tools::md5sum(file)))
