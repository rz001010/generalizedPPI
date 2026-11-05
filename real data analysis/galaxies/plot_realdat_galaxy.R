# Stacked distribution plot for the Galaxies PPI analysis.


library(tidyverse)
library(grid)


root <- "results/real_data"

result_dir <- file.path(
  root,
  "galaxies_ppi_lambda_0.08"
)

res <- read_csv(
  file.path(result_dir, "replication_level_results.csv"),
  show_col_types = FALSE
)


metric_names <- c("TPR", "FPR", "AUC")

r_colors <- palette.colors(8, "Okabe-Ito")

method_colors <- c(
  "PPI" = unname(r_colors[3]),
  "Labeled data" = unname(r_colors[4])
)

truth_color <- "red3"

res <- res %>%
  mutate(
    metric = factor(
      metric,
      levels = metric_names
    ),
    estimator = if_else(
      str_detect(method, "^PPI"),
      "PPI",
      "Labeled data"
    ),
    estimator = factor(
      estimator,
      levels = c("PPI", "Labeled data")
    )
  )


ridge_data <- function(x, value) {
  
  pieces <- x %>%
    group_split(metric, estimator) %>%
    map(function(d) {
      
      values <- d[[value]]
      
      den <- density(
        values,
        n = 500,
        adjust = 1.05,
        na.rm = TRUE
      )
      
      base <- if_else(
        as.character(d$estimator[1]) == "PPI",
        2,
        1
      )
      
      height <- 0.45 * den$y / max(den$y)
      mean_x <- mean(values, na.rm = TRUE)
      
      list(
        curve = tibble(
          metric = d$metric[1],
          estimator = d$estimator[1],
          x = den$x,
          base = base,
          top = base + height
        ),
        mark = tibble(
          metric = d$metric[1],
          estimator = d$estimator[1],
          mean = mean_x,
          base = base,
          top = base + approx(
            den$x,
            height,
            xout = mean_x,
            rule = 2
          )$y,
          xmin = min(den$x),
          xmax = max(den$x)
        )
      )
    })
  
  list(
    curve = map_dfr(pieces, "curve"),
    mark = map_dfr(pieces, "mark")
  )
}


ridge_plot <- function(
    density_data,
    reference = NULL,
    row_label,
    show_legend = FALSE,
    show_metric = TRUE,
    show_mean = TRUE,
    show_baseline = TRUE) {
  
  p <- ggplot(
    density_data$curve,
    aes(
      x = x,
      fill = estimator,
      color = estimator
    )
  ) +
    geom_ribbon(
      aes(
        ymin = base,
        ymax = top
      ),
      alpha = 0.72,
      color = NA
    ) +
    geom_line(
      aes(y = top),
      linewidth = 0.65
    ) +
    facet_wrap(
      ~metric,
      nrow = 1,
      scales = "free_x"
    ) +
    scale_fill_manual(
      values = method_colors,
      breaks = c("PPI", "Labeled data"),
      labels = c(
        "PPI",
        "Labeled data estimator"
      )
    ) +
    scale_color_manual(
      values = method_colors,
      breaks = c("PPI", "Labeled data"),
      labels = c(
        "PPI",
        "Labeled data estimator"
      )
    ) +
    scale_x_continuous(
      expand = expansion(
        mult = c(0.06, 0.22)
      )
    ) +
    scale_y_continuous(
      limits = c(0.90, 2.84),
      expand = c(0, 0)
    ) +
    labs(
      x = NULL,
      y = row_label,
      fill = NULL,
      color = NULL
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(alpha = 0.75)
      ),
      color = "none"
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = if (show_legend) "bottom" else "none",
      legend.text = element_text(size = 14),
      legend.key.width = unit(1.5, "lines"),
      legend.margin = margin(0, 0, 1, 0),
      strip.background = element_blank(),
      strip.text = if (show_metric) {
        element_text(
          face = "plain",
          size = 16
        )
      } else {
        element_blank()
      },
      axis.title.y = element_text(
        face = "bold",
        size = 15,
        angle = 0,
        vjust = 0.5,
        lineheight = 0.9,
        margin = margin(r = 8)
      ),
      axis.text.x = element_text(size = 13),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.spacing.x = unit(1.5, "lines"),
      plot.margin = margin(7, 24, 2, 10)
    ) +
    coord_cartesian(clip = "off")
  
  # Horizontal baselines
  if (show_baseline) {
    p <- p +
      geom_segment(
        data = density_data$mark,
        aes(
          x = xmin,
          xend = xmax,
          y = base,
          yend = base
        ),
        color = "black",
        linewidth = 0.35,
        inherit.aes = FALSE
      )
  }
  
  # Black vertical marks for the mean across replications
  if (show_mean) {
    p <- p +
      geom_segment(
        data = density_data$mark,
        aes(
          x = mean,
          xend = mean,
          y = base,
          yend = top
        ),
        color = "black",
        linewidth = 0.8,
        inherit.aes = FALSE
      )
  }
  
  # All-labeled reference values
  if (!is.null(reference)) {
    p <- p +
      geom_vline(
        data = reference,
        aes(
          xintercept = all_labeled_estimate
        ),
        color = truth_color,
        linewidth = 1.05,
        inherit.aes = FALSE
      ) +
      geom_text(
        data = reference,
        aes(
          x = all_labeled_estimate,
          y = 2.67,
          label = sprintf(
            "%.3f",
            all_labeled_estimate
          )
        ),
        color = truth_color,
        hjust = -0.03,
        size = 4.7,
        inherit.aes = FALSE
      )
  }
  
  p
}


reference <- res %>%
  distinct(
    metric,
    all_labeled_estimate
  )

point_estimates <- ridge_plot(
  ridge_data(res, "estimate"),
  reference = reference,
  row_label = "Point\nestimates",
  show_legend = FALSE,
  show_metric = TRUE,
  show_mean = TRUE,
  show_baseline = TRUE
)

standard_errors <- ridge_plot(
  ridge_data(res, "sd"),
  row_label = "Standard error\nestimates",
  show_legend = TRUE,
  show_metric = FALSE,
  show_mean = FALSE,
  show_baseline = FALSE
)


top <- ggplotGrob(point_estimates)
bottom <- ggplotGrob(standard_errors)

common_widths <- unit.pmax(
  top$widths,
  bottom$widths
)

top$widths <- common_widths
bottom$widths <- common_widths


output_file <- file.path(
  result_dir,
  "stacked_estimates_and_standard_errors.png"
)

png(
  output_file,
  width = 3000,
  height = 1150,
  res = 240,
  bg = "white"
)

grid.newpage()

layout <- grid.layout(
  2,
  1,
  heights = unit(
    c(1.05, 1.18),
    "null"
  )
)

pushViewport(
  viewport(layout = layout)
)

pushViewport(
  viewport(layout.pos.row = 1)
)

grid.draw(top)

popViewport()

pushViewport(
  viewport(layout.pos.row = 2)
)

grid.draw(bottom)

popViewport(2)

dev.off()

message(
  "Saved stacked plot in: ",
  output_file
)