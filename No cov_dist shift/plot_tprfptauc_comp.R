# Boxplots of replication-level estimated SEs for TPR, FPR, and AUC.
# Two methods: efficient projection and proposed PPI.


library(ggplot2)


efficient_results_dir <- "results/noncovshift_tprfprauc_eff"

own_ppi_results_dir <- "results/noncovshift_ppi"

output_dir <- file.path(
  efficient_results_dir,
  "sd_comparison_tpr_fpr_auc"
)

cache_dir <- file.path(
  output_dir,
  "scenario_cache_se_boxplots"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

n <- 1000
lambda_vals <- c(0.8, 0.5, 0.25, 0.1)
models <- c("all", "ideal", "noise", "bad", "pure_noise")
metrics <- c("tpr", "fpr", "auc")
progress_every <- 250

# Keep FALSE when changing only the appearance of the figure.
# Set TRUE only after adding or replacing underlying simulation results.
rebuild_cache <- FALSE

combined_cache_file <- file.path(
  output_dir,
  paste0("se_boxplot_inputs_n_", n, ".rds")
)

cache_settings <- list(
  efficient_results_dir = efficient_results_dir,
  own_ppi_results_dir = own_ppi_results_dir,
  n = n,
  lambda_vals = lambda_vals,
  models = models,
  metrics = metrics
)


status_message <- function(...) {
  message(
    "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
    paste0(..., collapse = "")
  )
  flush.console()
}

scenario_name <- function(model, lambda_value) {
  paste0(
    "f_", model,
    "_n_", n,
    "_lambda_",
    format(lambda_value, trim = TRUE, scientific = FALSE)
  )
}

read_replication <- function(sim_file, fields) {
  e <- new.env(parent = emptyenv())
  object_names <- load(sim_file, envir = e)
  
  candidates <- unique(c("res.i", object_names))
  candidates <- candidates[candidates %in% object_names]
  
  for (object_name in candidates) {
    object <- e[[object_name]]
    
    if (is.list(object) && all(fields %in% names(object))) {
      return(
        vapply(
          fields,
          function(field) {
            value <- as.numeric(object[[field]])
            
            if (length(value) != 1L) {
              stop(
                "Expected one value for ", field,
                " in: ", sim_file
              )
            }
            
            value
          },
          numeric(1)
        )
      )
    }
  }
  
  stop(
    "Required fields not found in: ", sim_file,
    "\nRequired: ", paste(fields, collapse = ", ")
  )
}



read_scenario <- function(root_dir, source, model, lambda_value) {
  scenario <- scenario_name(model, lambda_value)
  folder_path <- file.path(root_dir, scenario)
  
  cache_file <- file.path(
    cache_dir,
    paste0(source, "_", scenario, ".rds")
  )
  
  fields <- paste0(metrics, "_sd_", source, "_i")
  
  result_columns <- c(
    "model", "lambda", "metric", "series", "sd_value"
  )
  
  description <- if (source == "eff") {
    "Efficient projection"
  } else {
    "Proposed PPI estimator"
  }
  
  if (file.exists(cache_file) && !rebuild_cache) {
    cached <- readRDS(cache_file)
    
    # Reject old caches 
    valid_cache <- identical(cached$source_folder, folder_path) &&
      identical(cached$fields, fields) &&
      is.data.frame(cached$result) &&
      identical(names(cached$result), result_columns)
    
    if (valid_cache) {
      status_message(
        "Using replication-level cache: ",
        description, " / ", scenario
      )
      return(cached$result)
    }
    
    status_message(
      "Old summary/incompatible cache cannot form boxplots: ",
      description, " / ", scenario
    )
  }
  
  sim_files <- list.files(
    folder_path,
    pattern = "^res\\.[0-9]+\\.(RData|rda)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(sim_files) == 0L) {
    stop("No numbered replication files found in: ", folder_path)
  }
  
  status_message(
    "Reading ", length(sim_files), " replication files: ",
    description, " / ", scenario
  )
  
  values <- matrix(
    NA_real_,
    nrow = length(sim_files),
    ncol = length(fields),
    dimnames = list(NULL, fields)
  )
  
  for (i in seq_along(sim_files)) {
    values[i, ] <- read_replication(sim_files[i], fields)
    
    if (i %% progress_every == 0L || i == length(sim_files)) {
      status_message(
        "  ", description, ": ", i, "/", length(sim_files),
        " files; ", scenario
      )
    }
  }
  
  
  result <- do.call(
    rbind,
    lapply(metrics, function(metric) {
      se_values <- values[, paste0(metric, "_sd_", source, "_i")]
      valid <- is.finite(se_values)
      
      if (!any(valid)) {
        stop(
          "No finite SEs for ",
          description, " / ", metric, " / ", scenario
        )
      }
      
      if (any(!valid)) {
        status_message(
          "Using ", sum(valid), "/", length(valid),
          " finite SEs: ",
          description, " / ", metric, " / ", scenario
        )
      }
      
      data.frame(
        model = model,
        lambda = lambda_value,
        metric = metric,
        series = paste0(source, "_est"),
        sd_value = se_values[valid],
        stringsAsFactors = FALSE
      )
    })
  )
  
  saveRDS(
    list(
      source_folder = folder_path,
      fields = fields,
      result = result
    ),
    cache_file
  )
  
  status_message(
    "Saved replication-level cache: ",
    description, " / ", scenario
  )
  
  result
}


# Efficient-only table from saved scenario summaries 
#
# bias: mean point estimate minus saved target.
# sd: empirical SD of efficient-estimator point estimates.
# coverage_95: saved 95% CI coverage, on the 0-to-1 scale.

read_efficient_summary <- function(model, lambda_value) {
  scenario <- scenario_name(model, lambda_value)
  folder_path <- file.path(efficient_results_dir, scenario)
  
  summary_file <- list.files(
    folder_path,
    pattern = "^res\\.summary\\.(RData|rda)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(summary_file) != 1L) {
    stop("Expected one res.summary file in: ", folder_path)
  }
  
  status_message(
    "Reading efficient-estimator summary: ",
    scenario
  )
  
  fields <- unlist(lapply(metrics, function(metric) {
    paste0(
      metric,
      c(
        "_target_0",
        "_est_eff",
        "_empirical_sd_eff",
        "_cov_eff"
      )
    )
  }))
  
  saved <- read_replication(summary_file, fields)
  
  data.frame(
    n = n,
    lambda = lambda_value,
    model = model,
    estimand = toupper(metrics),
    
    bias = unname(
      saved[paste0(metrics, "_est_eff")] -
        saved[paste0(metrics, "_target_0")]
    ),
    
    sd = unname(
      saved[paste0(metrics, "_empirical_sd_eff")]
    ),
    
    coverage_95 = unname(
      saved[paste0(metrics, "_cov_eff")]
    ),
    
    stringsAsFactors = FALSE
  )
}


# Load combined cache, or build it once 

overall_start <- Sys.time()
plot_inputs <- NULL

if (file.exists(combined_cache_file) && !rebuild_cache) {
  cached <- readRDS(combined_cache_file)
  
  if (
    identical(cached$format, "raw_se_boxplots_v1") &&
    identical(cached$settings, cache_settings)
  ) {
    plot_inputs <- cached
    
    status_message(
      "Using combined cache; skipping replication and summary files: ",
      combined_cache_file
    )
  }
}

if (is.null(plot_inputs)) {
  rows <- list()
  table_rows <- list()
  scenario_number <- 0L
  total_scenarios <- length(models) * length(lambda_vals)
  
  for (lambda_value in lambda_vals) {
    for (model in models) {
      scenario_number <- scenario_number + 1L
      
      status_message(
        "Scenario ", scenario_number, "/", total_scenarios,
        ": f=", model, ", lambda=", lambda_value
      )
      
      rows[[length(rows) + 1L]] <- read_scenario(
        efficient_results_dir,
        "eff",
        model,
        lambda_value
      )
      
      rows[[length(rows) + 1L]] <- read_scenario(
        own_ppi_results_dir,
        "ppi",
        model,
        lambda_value
      )
      
      table_rows[[length(table_rows) + 1L]] <- read_efficient_summary(
        model,
        lambda_value
      )
    }
  }
  
  plot_inputs <- list(
    format = "raw_se_boxplots_v1",
    settings = cache_settings,
    sd_results = do.call(rbind, rows),
    efficient_table = do.call(rbind, table_rows)
  )
  
  saveRDS(plot_inputs, combined_cache_file)
  
  status_message(
    "Saved combined cache for future plotting: ",
    combined_cache_file
  )
  

  write.csv(
    plot_inputs$sd_results,
    file.path(output_dir, "tpr_fpr_auc_estimated_se_replications.csv"),
    row.names = FALSE
  )
}

sd_results <- plot_inputs$sd_results
efficient_table <- plot_inputs$efficient_table

table_file <- file.path(
  output_dir,
  "efficient_bias_sd_coverage95.csv"
)

write.csv(
  efficient_table,
  table_file,
  row.names = FALSE
)

print(
  efficient_table,
  row.names = FALSE,
  digits = 4
)

status_message(
  "Saved efficient-estimator table: ",
  table_file
)

if (any(!is.finite(as.matrix(
  efficient_table[, c("bias", "sd", "coverage_95")]
)))) {
  status_message(
    "Some saved summaries are nonfinite; retained as NA/NaN. ",
    "Inspect those scenarios before interpreting coverage."
  )
}


# PLOTTING SECTION

model_labels <- c(
  all = expression(RF * "(" * Y * plain("~") * X * ")"),
  ideal = expression(Pr(Y == 1 ~ "|" ~ X)),
  noise = expression(
    RF * "(" * Y * plain("~") * X[4] + X[5] * ")"
  ),
  bad = expression(
    RF * "(" * Y * plain("~") * R + X[1] + X[4] + X[5] * ")"
  ),
  pure_noise = "Unif[0.01,0.99]"
)

series_levels <- c("eff_est", "ppi_est")

series_labels <- c(
  eff_est = "Efficient projection",
  ppi_est = "Proposed PPI estimator"
)

series_colors <- c(
  eff_est = "#0072B2",
  ppi_est = "#D55E00"
)

series_offsets <- c(
  eff_est = -0.18,
  ppi_est =  0.18
)

plot_data <- sd_results

plot_data$model <- factor(
  plot_data$model,
  levels = models
)

plot_data$metric <- factor(
  plot_data$metric,
  levels = metrics,
  labels = c("TPR", "FPR", "AUC")
)

plot_data$lambda <- factor(
  plot_data$lambda,
  levels = lambda_vals,
  labels = paste0("lambda == ", lambda_vals)
)

plot_data$series <- factor(
  plot_data$series,
  levels = series_levels
)

plot_data$x_position <- as.numeric(plot_data$model) +
  unname(series_offsets[as.character(plot_data$series)])

status_message("Creating the stacked estimated-SE boxplots")

sd_plot <- ggplot(
  plot_data,
  aes(
    x = x_position,
    y = sd_value,
    fill = series,
    group = interaction(model, series, drop = TRUE)
  )
) +
  geom_boxplot(
    color = "gray20",
    width = 0.28,
    position = "identity",
    coef = 1.5,
    linewidth = 0.5,
    outlier.size = 0.8,
    outlier.alpha = 0.35
  ) +
  facet_grid(
    metric ~ lambda,
    scales = "free_y",
    switch = "y",
    labeller = labeller(lambda = label_parsed)
  ) +
  scale_x_continuous(
    breaks = seq_along(models),
    labels = model_labels,
    limits = c(0.5, length(models) + 0.5),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.08, 0.12))
  ) +
  scale_fill_manual(
    name = NULL,
    values = series_colors,
    breaks = series_levels,
    labels = series_labels,
    drop = FALSE
  ) +
  labs(
    x = "Prediction model",
    y = "Estimated standard error",
    caption = paste0(
      "Boxes: interquartile range; line: median; ",
      "whiskers: 1.5-IQR rule; points: outliers."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    text = element_text(size = 20, face = "plain"),
    axis.text = element_text(
      size = 16,
      face = "plain",
      color = "gray20"
    ),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      vjust = 1,
      size = 16,
      face = "plain"
    ),
    axis.title = element_text(
      size = 20,
      face = "plain"
    ),
    strip.text = element_text(
      size = 18,
      face = "plain"
    ),
    strip.text.y.left = element_text(
      angle = 0,
      face = "plain"
    ),
    strip.background = element_rect(
      fill = "gray96",
      color = "gray75"
    ),
    strip.placement = "outside",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(
      color = "gray92",
      linewidth = 0.3
    ),
    panel.spacing = grid::unit(0.8, "lines"),
    legend.position = "bottom",
    legend.text = element_text(
      size = 18,
      face = "plain"
    ),
    legend.key = element_blank(),
    legend.key.width = grid::unit(0.7, "cm"),
    plot.caption = element_text(
      size = 16,
      face = "plain",
      hjust = 0
    ),
    plot.margin = margin(10, 15, 10, 10)
  ) +
  guides(
    fill = guide_legend(nrow = 1)
  )


# Save figure 

png_file <- file.path(
  output_dir,
  "tpr_fpr_auc_estimated_se_boxplots.png"
)

ggsave(
  filename = png_file,
  plot = sd_plot,
  width = 17,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)

status_message("Saved figure: ", png_file)

status_message(
  "Complete. Total time: ",
  round(
    as.numeric(
      difftime(Sys.time(), overall_start, units = "mins")
    ),
    2
  ),
  " minutes"
)