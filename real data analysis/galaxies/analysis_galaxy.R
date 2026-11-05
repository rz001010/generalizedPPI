# Galaxies data: no-covariate-shift PPI for TPR, FPR, and AUC.
#
# Y    = human label indicating whether a galaxy has spiral arms
# f    = computer-vision prediction score
# R    = f, used as the continuous score for TPR, FPR, and AUC
#
# TPR and FPR use R > 0.6. AUC gives half credit to tied scores.
# The full cohort is fixed; only the labeled subset changes across repetitions.

library(tidyverse)
library(grid)



root <- "results/real_data"

cfg <- list(
  B = 1000,
  seed = 20250903,
  lambda = 0.08,                 # lambda = n/N
  cutoff = 0.6,
  level = 0.95,
  data_dir = file.path(root, "data"),
  out_dir = file.path(root, "galaxies_ppi_lambda_0.08")
)

dir.create(cfg$data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

# Download and read the NPZ file 

is_zip <- function(path) {
  if (!file.exists(path)) return(FALSE)
  con <- file(path, "rb")
  on.exit(close(con))
  identical(readBin(con, "raw", 2), charToRaw("PK"))
}

download_galaxies <- function(path) {
  if (is_zip(path)) return(invisible(path))
  
  id <- "1pDLQesPhbH5fSZW1m4aWC-wnJWnp1rGV"
  url <- paste0(
    "https://drive.usercontent.google.com/download?id=", id,
    "&export=download&confirm=t"
  )
  tmp <- tempfile(tmpdir = dirname(path), fileext = ".npz")
  on.exit(unlink(tmp))
  tryCatch(download.file(url, tmp, mode = "wb"), error = function(e) NULL)
  
  if (!is_zip(tmp)) {
    if (!requireNamespace("googledrive", quietly = TRUE)) {
      install.packages("googledrive", repos = "https://cloud.r-project.org")
    }
    googledrive::drive_deauth()
    googledrive::drive_download(
      googledrive::as_id(id), tmp, overwrite = TRUE
    )
  }
  
  if (!is_zip(tmp)) stop("The downloaded file is not a valid NPZ archive.")
  file.copy(tmp, path, overwrite = TRUE)
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
  
  if (type == "<i8") {
    bytes <- matrix(as.numeric(readBin(con, "raw", 8 * n)), nrow = 8)
    low <- colSums(bytes[1:4, , drop = FALSE] * 256^(0:3))
    high <- colSums(bytes[5:8, , drop = FALSE] * 256^(0:3))
    high <- high - (high >= 2^31) * 2^32
    value <- low + high * 2^32
  } else {
    value <- readBin(
      con, "double", n,
      size = as.integer(substr(type, 3, 3)), endian = "little"
    )
  }
  
  if (length(dims) == 2) {
    value <- matrix(value, dims[1], dims[2], byrow = !fortran)
  }
  value
}

read_galaxies <- function(path) {
  tmp <- tempfile("galaxies_", tmpdir = dirname(path))
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))
  
  files <- unzip(path, list = TRUE)$Name
  files <- files[match(c("Y.npy", "Yhat.npy"), basename(files))]
  unzip(path, files = files, exdir = tmp, junkpaths = TRUE)
  
  data.frame(
    Y = as.numeric(read_npy(file.path(tmp, "Y.npy"))),
    f = as.numeric(read_npy(file.path(tmp, "Yhat.npy")))
  )
}

# TPR, FPR, AUC, and influence functions 


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

metrics <- function(R, Y, cutoff, w = rep(1, length(Y))) {
  H <- as.numeric(R > cutoff)
  w1 <- w * Y
  w0 <- w * (1 - Y)
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

# Section 6.1 PPI estimator 

ppi_fit <- function(d, C, full_prediction, cutoff) {
  n <- sum(C)
  N <- nrow(d) - n
  lab <- d[C == 1, ]
  
  y <- metrics(lab$R, lab$Y, cutoff)
  f <- metrics(lab$R, lab$f, cutoff)
  vy <- col_var(y$phi)
  vf <- col_var(full_prediction$phi)
  cyf <- diag(cov(y$phi, f$phi))
  
  omega <- cyf / vf
  correction <- full_prediction$theta - f$theta
  ppi <- y$theta + omega * correction
  
  # Analytic variance from Section 6.1, with lambda = n/N.
  ppi_var <- (vy - cyf^2 / ((1 + n / N) * vf)) / n
  ppi_var[ppi_var <= 0] <- NA_real_
  
  list(
    labeled = y$theta,
    ppi = ppi,
    sd_labeled = sqrt(vy / n),
    sd_ppi = sqrt(ppi_var),
    omega = omega
  )
}



run_analysis <- function(d, cfg) {
  M <- nrow(d)
  n <- round(cfg$lambda * M / (1 + cfg$lambda))
  N <- M - n
  z <- qnorm((1 + cfg$level) / 2)
  
  # Yhat is both the prediction f and the continuous score R.
  d$R <- d$f
  reference_fit <- metrics(d$R, d$Y, cfg$cutoff)
  reference <- reference_fit$theta
  full_prediction <- metrics(d$R, d$f, cfg$cutoff)
  
  reference_table <- data.frame(
    metric = names(reference),
    all_labeled_estimate = unname(reference),
    analytic_sd_all_labeled = sqrt(col_var(reference_fit$phi) / M),
    n_total = M
  )
  write.csv(
    reference_table,
    file.path(cfg$out_dir, "all_labeled_reference.csv"),
    row.names = FALSE
  )
  
  seeds <- data.frame(
    replication = seq_len(cfg$B),
    seed = cfg$seed + 2 * (seq_len(cfg$B) - 1)
  )
  write.csv(
    seeds, file.path(cfg$out_dir, "replication_seeds.csv"),
    row.names = FALSE
  )
  
  results <- vector("list", cfg$B)
  for (b in seq_len(cfg$B)) {
    set.seed(seeds$seed[b])
    C <- integer(M)
    C[sample.int(M, n)] <- 1L
    fit <- ppi_fit(d, C, full_prediction, cfg$cutoff)
    
    estimate <- c(fit$labeled, fit$ppi)
    sd <- c(fit$sd_labeled, fit$sd_ppi)
    truth <- rep(reference, 2)
    
    results[[b]] <- data.frame(
      replication = b,
      seed = seeds$seed[b],
      method = rep(c("Labeled only", "PPI Section 6.1"), each = 3),
      metric = rep(names(reference), 2),
      estimate = unname(estimate),
      all_labeled_estimate = unname(truth),
      bias = unname(estimate - truth),
      sd = unname(sd),
      lower = unname(estimate - z * sd),
      upper = unname(estimate + z * sd),
      omega = c(rep(NA_real_, 3), unname(fit$omega)),
      n_labeled = n,
      n_unlabeled = N,
      lambda_actual = n / N
    )
    
    if (b == 1) {
      write.csv(
        results[[b]], file.path(cfg$out_dir, "point_estimates.csv"),
        row.names = FALSE
      )
      write.csv(
        data.frame(row_id = seq_len(M), labeled = C),
        file.path(cfg$out_dir, "replication_1_assignments.csv"),
        row.names = FALSE
      )
    }
  }
  
  results <- bind_rows(results)
  write.csv(
    results, file.path(cfg$out_dir, "replication_level_results.csv"),
    row.names = FALSE
  )
  
  summary <- results %>%
    group_by(method, metric) %>%
    summarise(
      all_labeled_estimate = first(all_labeled_estimate),
      repetitions = n(),
      mean_estimate = mean(estimate, na.rm = TRUE),
      mean_bias = mean(bias, na.rm = TRUE),
      mean_analytic_sd = mean(sd, na.rm = TRUE),
      empirical_sd = sd(estimate, na.rm = TRUE),
      coverage = mean(lower <= all_labeled_estimate &
                        all_labeled_estimate <= upper, na.rm = TRUE),
      mean_omega = mean(omega, na.rm = TRUE),
      valid_estimates = sum(is.finite(estimate)),
      valid_analytic_sds = sum(is.finite(sd)),
      n_labeled = first(n_labeled),
      n_unlabeled = first(n_unlabeled),
      lambda_actual = first(lambda_actual),
      .groups = "drop"
    ) %>%
    mutate(
      metric = factor(metric, levels = c("TPR", "FPR", "AUC")),
      method = factor(method, levels = c("Labeled only", "PPI Section 6.1"))
    ) %>%
    arrange(method, metric)
  
  write.csv(
    summary, file.path(cfg$out_dir, "operating_characteristics.csv"),
    row.names = FALSE
  )
  
  ratio <- summary %>%
    select(metric, method, mean_analytic_sd) %>%
    pivot_wider(names_from = method, values_from = mean_analytic_sd) %>%
    mutate(sd_ratio = `PPI Section 6.1` / `Labeled only`)
  
  write.csv(
    ratio, file.path(cfg$out_dir, "sd_ratio.csv"),
    row.names = FALSE
  )
  
  analysis <- list(
    settings = cfg,
    reference = reference_table,
    results = results,
    summary = summary,
    sd_ratio = ratio
  )
  saveRDS(
    analysis,
    file.path(
      cfg$out_dir,
      paste0("galaxies_ppi_", cfg$B, "_repetitions.rds")
    )
  )
  
  message(
    "M = ", M, ", n = ", n, ", N = ", N,
    ", actual lambda = ", signif(n / N, 7)
  )
  print(summary, n = Inf)
  print(ratio, n = Inf)
  analysis
}



data_file <- file.path(cfg$data_dir, "galaxies.npz")
download_galaxies(data_file)
galaxies <- read_galaxies(data_file)
analysis <- run_analysis(galaxies, cfg)


message("Results saved in: ", cfg$out_dir)
