# Mean estimated-SE ratios for four PPI estimators.
# Ratio = mean of the replication-level SE ratios.

library(ggplot2)


# Directories and settings 

package_results_dir <- "results/noncovshift_mean"
own_ppi_results_dir <- "results/noncovshift_ppi"

output_dir <- file.path(package_results_dir, "four_ppi_ratio_comparison")
cache_dir <- file.path(output_dir, "scenario_cache")
final_csv <- file.path(output_dir, "mean_se_ratio_four_ppi_estimators.csv")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

n <- 1000
lambda_vals <- c(0.8, 0.5, 0.25, 0.1)
models <- c("all", "ideal", "noise", "bad", "pure_noise")
progress_every <- 250

reuse_completed_csv <- TRUE
# Set TRUE only if underlying simulation results have changed.
# This overrides reuse_completed_csv and rebuilds all scenario summaries.
rebuild_cache <- FALSE

method_levels <- c(
  "Original PPI", "PPI++", "Calibrated PPI", "Proposed PPI estimator"
)
package_method_fields <- c(
  "Original PPI" = "mean_sd_ppi_i",
  "PPI++" = "mean_sd_ppipp_i",
  "Calibrated PPI" = "mean_sd_calppi_i"
)
own_ppi_method_field <- c(
  "Proposed PPI estimator" = "mean_sd_ppi_i"
)


status_message <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
          paste0(..., collapse = ""))
  flush.console()
}

lambda_text <- function(x) {
  format(x, trim = TRUE, scientific = FALSE)
}

# Accept older CSVs/caches without recalculating their numerical results.
standardize_results <- function(result) {
  result$method <- as.character(result$method)
  old_name <- tolower(trimws(result$method)) == "our construction ppi"
  result$method[which(old_name)] <- "Proposed PPI estimator"
  if (!"metric" %in% names(result)) result$metric <- "mean"
  result[, c("model", "lambda", "metric", "method",
             "mean_ratio", "valid_replications"), drop = FALSE]
}

load_replication_result <- function(sim_file, required_fields) {
  e <- new.env(parent = emptyenv())
  object_names <- load(sim_file, envir = e)
  candidates <- unique(c("res.i", object_names))
  candidates <- candidates[candidates %in% object_names]
  
  for (object_name in candidates) {
    object <- e[[object_name]]
    if (is.list(object) && all(required_fields %in% names(object))) {
      return(object)
    }
  }
  
  stop("Required mean-SE fields not found in: ", sim_file,
       "\nRequired: ", paste(required_fields, collapse = ", "))
}




calculate_scenario_ratios <- function(
  root_dir, model, lambda_value, method_fields, source_label) {
  
  scenario <- paste0("f_", model, "_n_", n,
                     "_lambda_", lambda_text(lambda_value))
  folder_path <- file.path(root_dir, scenario)
  cache_file <- file.path(cache_dir, paste0(source_label, "_", scenario, ".rds"))
  description <- paste(names(method_fields), collapse = ", ")
  
  if (file.exists(cache_file) && !rebuild_cache) {
    status_message("Using cache: ", description, " / ", scenario)
    return(standardize_results(readRDS(cache_file)))
  }
  
  status_message("Finding files: ", description, " / ", folder_path)
  sim_files <- list.files(
    folder_path,
    pattern = "^res\\.[0-9]+\\.(RData|rda)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(sim_files) == 0L) {
    stop("No numbered replication files found in: ", folder_path)
  }
  
  status_message("Reading ", length(sim_files), " files: ", description)
  ratio_values <- matrix(
    NA_real_, nrow = length(sim_files), ncol = length(method_fields),
    dimnames = list(NULL, names(method_fields))
  )
  required_fields <- unique(c("mean_sd_n_i", unname(method_fields)))
  
  for (i in seq_along(sim_files)) {
    result <- load_replication_result(sim_files[i], required_fields)
    labeled_sd <- as.numeric(result[["mean_sd_n_i"]])
    
    for (method_name in names(method_fields)) {
      method_sd <- as.numeric(result[[method_fields[[method_name]]]])
      if (length(labeled_sd) == 1L && length(method_sd) == 1L) {
        ratio_values[i, method_name] <- method_sd / labeled_sd
      }
    }
    
    if (i %% progress_every == 0L || i == length(sim_files)) {
      status_message("  Read ", i, "/", length(sim_files),
                     " files: ", description, " / ", scenario)
    }
  }
  
  scenario_result <- data.frame(
    model = model,
    lambda = lambda_value,
    metric = "mean",
    method = names(method_fields),
    mean_ratio = colMeans(ratio_values, na.rm = TRUE),
    valid_replications = colSums(!is.na(ratio_values)),
    stringsAsFactors = FALSE
  )
  saveRDS(scenario_result, cache_file)
  status_message("Saved scenario cache: ", description, " / ", scenario)
  scenario_result
}


# Obtain summaries

get_mean_results <- function() {
  if (reuse_completed_csv && !rebuild_cache && file.exists(final_csv)) {
    status_message("Using completed mean-ratio CSV: ", final_csv)
    result <- standardize_results(read.csv(
      final_csv, stringsAsFactors = FALSE, check.names = FALSE
    ))
    return(result[
      result$metric == "mean" & result$lambda %in% lambda_vals &
        result$model %in% models,
      , drop = FALSE
    ])
  }
  
  rows <- list()
  scenario_number <- 0L
  total_scenarios <- length(models) * length(lambda_vals)
  status_message("Reading summaries for ", total_scenarios, " scenarios")
  
  for (lambda_value in lambda_vals) {
    for (model in models) {
      scenario_number <- scenario_number + 1L
      scenario_start <- Sys.time()
      status_message("Scenario ", scenario_number, "/", total_scenarios,
                     ": f=", model, ", lambda=", lambda_text(lambda_value))
      
      rows[[length(rows) + 1L]] <- calculate_scenario_ratios(
        package_results_dir, model, lambda_value,
        package_method_fields, "package_methods"
      )
      rows[[length(rows) + 1L]] <- calculate_scenario_ratios(
        own_ppi_results_dir, model, lambda_value,
        own_ppi_method_field, "own_ppi"
      )
      
      write.csv(do.call(rbind, rows),
                file.path(output_dir, "mean_ratio_results_in_progress.csv"),
                row.names = FALSE)
      status_message("Finished scenario ", scenario_number, "/", total_scenarios,
                     " in ", round(as.numeric(difftime(
                       Sys.time(), scenario_start, units = "mins"
                     )), 2), " minutes")
    }
  }
  
  result <- do.call(rbind, rows)
  write.csv(result, final_csv, row.names = FALSE)
  status_message("Saved mean-ratio summary: ", final_csv)
  result
}

overall_start <- Sys.time()
ratio_results <- get_mean_results()


# Four-panel mean-ratio plot 

model_labels <- c(
  all = expression(RF * "(" * Y * plain("~") * X * ")"),
  ideal = expression(Pr(Y == 1 ~ "|" ~ X)),
  noise = expression(RF * "(" * Y * plain("~") * X[4] + X[5] * ")"),
  bad = expression(RF * "(" * Y * plain("~") * R + X[1] + X[4] + X[5] * ")"),
  pure_noise = "Unif[0.01,0.99]"
)
method_colors <- c(
  "Original PPI" = "salmon",
  "PPI++" = "turquoise3",
  "Calibrated PPI" = "orchid",
  "Proposed PPI estimator" = "chartreuse4"
)
method_shapes <- c(
  "Original PPI" = 16,
  "PPI++" = 17,
  "Calibrated PPI" = 15,
  "Proposed PPI estimator" = 18
)

plot_data <- ratio_results
plot_data$model <- factor(plot_data$model, levels = models)
plot_data$method <- factor(plot_data$method, levels = method_levels)
plot_data$lambda <- factor(
  plot_data$lambda, levels = lambda_vals,
  labels = paste0("lambda == ", lambda_vals)
)

status_message("Creating the mean-ratio plot")

ratio_plot <- ggplot(
  plot_data,
  aes(x = model, y = mean_ratio, color = method, shape = method,
      group = method)
) +
  
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray45",
             linewidth = 0.6) +
  
  geom_point(size = 4, position = position_dodge(width = 0.55)) +
  facet_wrap(~lambda, nrow = 1, labeller = label_parsed) +
  scale_x_discrete(labels = model_labels, drop = FALSE) +
  scale_color_manual(values = method_colors, breaks = method_levels) +
  scale_shape_manual(values = method_shapes, breaks = method_levels) +
  labs(
    x = "Prediction model",
    y = expression(paste("mean(",
                         hat(sigma)[plain(PPI)] / hat(sigma)[plain(lab)], ")")),
    color = "PPI estimator",
    shape = "PPI estimator"
  ) +
  theme_bw(base_size = 20) +
  theme(
    text = element_text(face = "plain"),
    axis.text = element_text(size = 16),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.title = element_text(size = 20),
    axis.title.y = element_text(margin = margin(r = 16)),
    strip.text = element_text(size = 18, face = "plain"),
    panel.spacing = grid::unit(1, "lines"),
    legend.position = "bottom",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 18, face = "plain"),
    legend.key.width = grid::unit(0.8, "cm"),
    
    plot.margin = margin(t = 18, r = 20, b = 18, l = 30, unit = "pt")
  ) +
  guides(
    color = guide_legend(nrow = 1, override.aes = list(size = 4.5)),
    shape = guide_legend(nrow = 1)
  )

plot_file <- file.path(output_dir, "mean_se_ratio_four_ppi_estimators.png")
ggsave(
  filename = plot_file, plot = ratio_plot,
  width = 18, height = 6.2, units = "in", dpi = 300, bg = "white"
)

status_message("Saved: ", plot_file)
status_message("Complete. Total time: ", round(as.numeric(difftime(
  Sys.time(), overall_start, units = "mins"
)), 2), " minutes")
