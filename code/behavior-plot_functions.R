
## LINEAR PLOTS

single_linear_plot <- function(linear_array, type_id, metric, type_color, ylims, y_label){

  plot <- ggplot(linear_array[linear_array$type == type_id,],
                 aes(x = ap_coords, y = .data[[metric]], group = type)) +
    geom_point(aes(fill = type), pch = 21, colour = "black") +
    scale_colour_manual(values = type_color) +
    scale_fill_manual(values = type_color) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous(limits = c(-2.6, -1.2), expand = c(.01,.01)) +
    scale_y_continuous(limits = ylims) +
    ylab(y_label) + xlab("AP location (mm from bregma)") +
    geom_smooth(aes(fill = type, colour = type), alpha = 0.2, method = "lm", fullrange = TRUE, se = TRUE) +
    theme_classic() +
    theme(axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
          axis.line.x = element_blank(),
          legend.position = "none", legend.title = element_blank())

  return(plot)
}

dual_linear_plot <- function(linear_array, metric, ylims, y_label){

  plot <- ggplot(linear_array, aes(x = ap_coords, y = .data[[metric]], group = type)) +
    geom_point(aes(fill = type), pch = 21, colour = "black") +
    scale_colour_manual(values = c(chr2_blue, eyfp_green)) +
    scale_fill_manual(values = c(chr2_blue, eyfp_green)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous(limits = c(-2.6, -1.2), expand = c(.01,.01)) +
    scale_y_continuous(limits = ylims) +
    ylab(y_label) + xlab("AP location (mm from bregma)") +
    geom_smooth(aes(fill = type, colour = type), alpha = 0.2, method = "lm", fullrange = TRUE, se = TRUE) +
    theme_classic() +
    theme(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"),
          legend.title = element_blank())

  return(plot)
}

## HEATMAPS

get_heatmap <- function(array){

  baseline_array <- lapply(array, function(x) lapply(x, function(y) y[["baseline"]][,c("XCoordinate", "YCoordinate")]))
  treatment_array <- lapply(array, function(x) lapply(x, function(y) y[["treatment"]][,c("XCoordinate", "YCoordinate")]))

  baseline_array <- lapply(baseline_array, bind_rows)
  treatment_array <- lapply(treatment_array, bind_rows)

  baseline_image <- lapply(baseline_array, function(x) matrix(data = 0, nrow = 591, ncol = 591,
                                                              dimnames = list(as.character(round(seq(from = -29.5, to = 29.5, by = 0.1),
                                                                                                 digits = 1)),
                                                                              as.character(round(seq(from = -29.5, to = 29.5, by = 0.1),
                                                                                                 digits = 1)))))
  treatment_image <- lapply(treatment_array, function(x) matrix(data = 0, nrow = 591, ncol = 591,
                                                                dimnames = list(as.character(round(seq(from = -29.5, to = 29.5, by = 0.1),
                                                                                                   digits = 1)),
                                                                                as.character(round(seq(from = -29.5, to = 29.5, by = 0.1),
                                                                                                   digits = 1)))))
  # equalize proportions
  for (i in 1:length(baseline_image)){
    for (j in 1:length(rownames(baseline_array[[i]]))){
      baseline_image[[i]][as.character(round(baseline_array[[i]]$XCoordinate[j], digits = 1)),
                     as.character(round(baseline_array[[i]]$YCoordinate[j], digits = 1))] <- baseline_image[[i]][as.character(round(baseline_array[[i]]$XCoordinate[j], digits = 1)),
                       as.character(round(baseline_array[[i]]$YCoordinate[j], digits = 1))] + (1 / length(rownames(baseline_array[[i]])))
    }
  }


  for (i in 1:length(treatment_image)){
    for (j in 1:length(rownames(treatment_array[[i]]))){
      treatment_image[[i]][as.character(round(treatment_array[[i]]$XCoordinate[j], digits = 1)),
                          as.character(round(treatment_array[[i]]$YCoordinate[j], digits = 1))] <- treatment_image[[i]][as.character(round(treatment_array[[i]]$XCoordinate[j], digits = 1)),
                            as.character(round(treatment_array[[i]]$YCoordinate[j], digits = 1))] + (1 / length(rownames(treatment_array[[i]])))
    }
  }

  combined_image <- mapply(function(x, y) x - y, treatment_image, baseline_image, SIMPLIFY = FALSE)
  combined_image <- lapply(combined_image, function(x) im(x))
  combined_image <- lapply(combined_image, function(x) blur(x, sigma = 7))
  combined_image <- lapply(combined_image, function(x) x$v)
  combined_image <- lapply(combined_image, function(x) reshape2::melt(x))
  combined_image <- lapply(combined_image, function(x) {colnames(x) <- c("x", "y", "z"); x})

  return(combined_image)

}

get_heatmap_epm <- function(array, slot){

  array <- lapply(array, function(x) lapply(x, function(y) y[[slot]][,c("XCoordinate", "YCoordinate")]))

  array <- lapply(array, bind_rows)

  image <- lapply(array, function(x) matrix(data = 0, nrow = 741, ncol = 741,
                                            dimnames = list(as.character(round(seq(from = -37.0, to = 37.0, by = 0.1), digits = 1)),
                                                            as.character(round(seq(from = -37.0, to = 37.0, by = 0.1), digits = 1)))))

  # equalize proportions
  for (i in 1:length(image)){
    for (j in 1:length(rownames(array[[i]]))){
      image[[i]][as.character(round(array[[i]]$XCoordinate[j], digits = 1)),
                 as.character(round(array[[i]]$YCoordinate[j], digits = 1))] <- image[[i]][as.character(round(array[[i]]$XCoordinate[j], digits = 1)),
                                                                                           as.character(round(array[[i]]$YCoordinate[j], digits = 1))] + 1
    }
  }

  image <- lapply(image, function(x) im(x))
  image <- lapply(image, function(x) blur(x, sigma = 12))
  image <- lapply(image, function(x) x$v)
  image <- lapply(image, function(x) reshape2::melt(x))
  image <- lapply(image, function(x) {colnames(x) <- c("x", "y", "z"); x})

  return(image)

}

get_heatmap_epm_diff <- function(array){

  baseline_array <- lapply(array, function(x) lapply(x, function(y) y[["baseline_all"]][,c("XCoordinate", "YCoordinate")]))
  treatment_array <- lapply(array, function(x) lapply(x, function(y) y[["treatment"]][,c("XCoordinate", "YCoordinate")]))

  baseline_array <- lapply(baseline_array, bind_rows)
  treatment_array <- lapply(treatment_array, bind_rows)

  baseline_image <- lapply(baseline_array, function(x) matrix(data = 0, nrow = 741, ncol = 741,
                                                     dimnames = list(as.character(round(seq(from = -37.0, to = 37.0, by = 0.1), digits = 1)),
                                                                     as.character(round(seq(from = -37.0, to = 37.0, by = 0.1), digits = 1)))))

  treatment_image <- lapply(treatment_array, function(x) matrix(data = 0, nrow = 741, ncol = 741,
                                                      dimnames = list(as.character(round(seq(from = -37.0, to = 37.0, by = 0.1), digits = 1)),
                                                                      as.character(round(seq(from = -37.0, to = 37.0, by = 0.1), digits = 1)))))

  # equalize proportions
  for (i in 1:length(baseline_image)){
    for (j in 1:length(rownames(baseline_array[[i]]))){
      baseline_image[[i]][as.character(round(baseline_array[[i]]$XCoordinate[j], digits = 1)),
                          as.character(round(baseline_array[[i]]$YCoordinate[j], digits = 1))] <- baseline_image[[i]][as.character(round(baseline_array[[i]]$XCoordinate[j], digits = 1)),
                                                                                                                      as.character(round(baseline_array[[i]]$YCoordinate[j], digits = 1))] + (1 / length(rownames(baseline_array[[i]])))
    }
  }


  for (i in 1:length(treatment_image)){
    for (j in 1:length(rownames(treatment_array[[i]]))){
      treatment_image[[i]][as.character(round(treatment_array[[i]]$XCoordinate[j], digits = 1)),
                           as.character(round(treatment_array[[i]]$YCoordinate[j], digits = 1))] <- treatment_image[[i]][as.character(round(treatment_array[[i]]$XCoordinate[j], digits = 1)),
                                                                                                                         as.character(round(treatment_array[[i]]$YCoordinate[j], digits = 1))] + (1 / length(rownames(treatment_array[[i]])))
    }
  }

  combined_image <- mapply(function(x, y) x - y, treatment_image, baseline_image, SIMPLIFY = FALSE)
  combined_image <- lapply(combined_image, function(x) im(x))
  combined_image <- lapply(combined_image, function(x) blur(x, sigma = 12))
  combined_image <- lapply(combined_image, function(x) x$v)
  combined_image <- lapply(combined_image, function(x) reshape2::melt(x))
  combined_image <- lapply(combined_image, function(x) {colnames(x) <- c("x", "y", "z"); x})

  return(combined_image)

}

get_heatmap_openfield <- function(array, slot){

  array <- lapply(array, function(x) lapply(x, function(y) y[[slot]][,c("XCoordinate", "YCoordinate")]))

  array <- lapply(array, bind_rows)

  image <- lapply(array, function(x) matrix(data = 0, nrow = 275, ncol = 275,
                  dimnames = list(as.character(round(seq(from = -13.7, to = 13.7, by = 0.1), digits = 1)),
                                  as.character(round(seq(from = -13.7, to = 13.7, by = 0.1), digits = 1)))))

  # equalize proportions
  for (i in 1:length(image)){
    for (j in 1:length(rownames(array[[i]]))){
      image[[i]][as.character(round(array[[i]]$XCoordinate[j], digits = 1)),
                          as.character(round(array[[i]]$YCoordinate[j], digits = 1))] <- image[[i]][as.character(round(array[[i]]$XCoordinate[j], digits = 1)),
                                                                                                    as.character(round(array[[i]]$YCoordinate[j], digits = 1))] + 1
    }
  }

  image <- lapply(image, function(x) im(x))
  image <- lapply(image, function(x) blur(x, sigma = 7))
  image <- lapply(image, function(x) x$v)
  image <- lapply(image, function(x) reshape2::melt(x))
  image <- lapply(image, function(x) {colnames(x) <- c("x", "y", "z"); x})

  return(image)

}

get_heatmap_openfield_diff <- function(array){

  baseline_array <- lapply(array, function(x) lapply(x, function(y) y[["baseline_all"]][,c("XCoordinate", "YCoordinate")]))
  treatment_array <- lapply(array, function(x) lapply(x, function(y) y[["treatment_all"]][,c("XCoordinate", "YCoordinate")]))

  baseline_array <- lapply(baseline_array, bind_rows)
  treatment_array <- lapply(treatment_array, bind_rows)

  baseline_image <- lapply(baseline_array, function(x) matrix(data = 0, nrow = 275, ncol = 275,
                                                     dimnames = list(as.character(round(seq(from = -13.7, to = 13.7, by = 0.1), digits = 1)),
                                                                     as.character(round(seq(from = -13.7, to = 13.7, by = 0.1), digits = 1)))))

  treatment_image <- lapply(treatment_array, function(x) matrix(data = 0, nrow = 275, ncol = 275,
                                                      dimnames = list(as.character(round(seq(from = -13.7, to = 13.7, by = 0.1), digits = 1)),
                                                                      as.character(round(seq(from = -13.7, to = 13.7, by = 0.1), digits = 1)))))

  # equalize proportions
  for (i in 1:length(baseline_image)){
    for (j in 1:length(rownames(baseline_array[[i]]))){
      baseline_image[[i]][as.character(round(baseline_array[[i]]$XCoordinate[j], digits = 1)),
                          as.character(round(baseline_array[[i]]$YCoordinate[j], digits = 1))] <- baseline_image[[i]][as.character(round(baseline_array[[i]]$XCoordinate[j], digits = 1)),
                                                                                                                      as.character(round(baseline_array[[i]]$YCoordinate[j], digits = 1))] + (1 / length(rownames(baseline_array[[i]])))
    }
  }


  for (i in 1:length(treatment_image)){
    for (j in 1:length(rownames(treatment_array[[i]]))){
      treatment_image[[i]][as.character(round(treatment_array[[i]]$XCoordinate[j], digits = 1)),
                           as.character(round(treatment_array[[i]]$YCoordinate[j], digits = 1))] <- treatment_image[[i]][as.character(round(treatment_array[[i]]$XCoordinate[j], digits = 1)),
                                                                                                                         as.character(round(treatment_array[[i]]$YCoordinate[j], digits = 1))] + (1 / length(rownames(treatment_array[[i]])))
    }
  }

  combined_image <- mapply(function(x, y) x - y, treatment_image, baseline_image, SIMPLIFY = FALSE)
  combined_image <- lapply(combined_image, function(x) im(x))
  combined_image <- lapply(combined_image, function(x) blur(x, sigma = 7))
  combined_image <- lapply(combined_image, function(x) x$v)
  combined_image <- lapply(combined_image, function(x) reshape2::melt(x))
  combined_image <- lapply(combined_image, function(x) {colnames(x) <- c("x", "y", "z"); x})

  return(combined_image)

}

plot_heatmap <- function(heatmap_array){

  limits <- max(abs(heatmap_array$z)) * c(-1, 1)

  plot <- ggplot(heatmap_array, aes(x, y, fill = z)) +
    geom_raster() +
    #geom_rect(aes(xmin = 0, xmax = max(.data$x), ymin = 0, ymax = min(.data$y)),
    #    fill = "grey", alpha = 0.03) +
    #geom_segment(aes(x = 0, xend = -Inf, y = 295.5, yend = 295.5), linetype = 2, colour = chr2_blue) +
    #geom_segment(aes(y = 0, yend = -Inf, x = 295.5, xend = 295.5), linetype = 2, colour = chr2_blue) +
    geom_vline(xintercept = 295.5, linetype = 2, colour = chr2_blue) +
    geom_hline(yintercept = 295.5, linetype = 2, colour = chr2_blue) +
    #scale_fill_gradientn(colours = rev(brewer.pal(n = 100, name = "RdBu")), limit = limits) +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", p1 = 0.4, p2 = 0.4) +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    theme(panel.background = element_blank(), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

  return(plot)
}

plot_heatmap_epm <- function(heatmap_array){

  limits <- max(abs(heatmap_array$z)) * c(-1, 1)

  plot <- ggplot(heatmap_array, aes(x, y, fill = z)) +
    geom_raster() +
    annotate("rect", xmin = 0, xmax = 345, ymin = 0, ymax = 345, fill = "white", colour = NA, size = 0.5) +
    annotate("rect", xmin = 395, xmax = 741, ymin = 0, ymax = 345, fill = "white", colour = NA, size = 0.5) +
    annotate("rect", xmin = 0, xmax = 345, ymin = 395, ymax = 741, fill = "white", colour = NA, size = 0.5) +
    annotate("rect", xmin = 395, xmax = 741, ymin = 395, ymax = 741, fill = "white", colour = NA, size = 0.5) +
    annotate("rect", xmin = 345, xmax = 395, ymin = 0, ymax = 345, colour = "blue", fill = NA, size = 0.5) +
    annotate("rect", xmin = 345, xmax = 395, ymin = 395, ymax = 741, colour = "blue", fill = NA, size = 0.5) +
    annotate("rect", ymin = 345, ymax = 395, xmin = 0, xmax = 345, colour = "red", fill = NA, size = 0.5) +
    annotate("rect", ymin = 345, ymax = 395, xmin = 395, xmax = 741, colour = "red", fill = NA, size = 0.5) +
    scale_fill_viridis() +
    theme(panel.background = element_blank(), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

  return(plot)
}

plot_heatmap_epm_diff <- function(heatmap_array){

  limits <- max(abs(heatmap_array$z)) * c(-1, 1)

  plot <- ggplot(heatmap_array, aes(x, y, fill = z)) +
    geom_raster() +
    annotate("rect", xmin = 0, xmax = 345, ymin = 0, ymax = 345, fill = "white", colour = NA, size = 0.5) +
    annotate("rect", xmin = 395, xmax = 741, ymin = 0, ymax = 345, fill = "white", colour = NA, size = 0.5) +
    annotate("rect", xmin = 0, xmax = 345, ymin = 395, ymax = 741, fill = "white", colour = NA, size = 0.5) +
    annotate("rect", xmin = 395, xmax = 741, ymin = 395, ymax = 741, fill = "white", colour = NA, size = 0.5) +
    annotate("rect", xmin = 345, xmax = 395, ymin = 0, ymax = 345, colour = "purple", fill = NA, size = 0.5) +
    annotate("rect", xmin = 345, xmax = 395, ymin = 395, ymax = 741, colour = "purple", fill = NA, size = 0.5) +
    annotate("rect", ymin = 345, ymax = 395, xmin = 0, xmax = 345, colour = "darkgreen", fill = NA, size = 0.5) +
    annotate("rect", ymin = 345, ymax = 395, xmin = 395, xmax = 741, colour = "darkgreen", fill = NA, size = 0.5) +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", p1 = 0.4, p2 = 0.4) +
    theme(panel.background = element_blank(), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

  return(plot)
}

plot_heatmap_openfield <- function(heatmap_array){

  limits <- max(abs(heatmap_array$z)) * c(-1, 1)

  plot <- ggplot(heatmap_array, aes(x, y, fill = z)) +
    geom_raster() +
    geom_vline(xintercept = 68, linetype = 2) + geom_vline(xintercept = 206, linetype = 2) +
    geom_hline(yintercept = 68, linetype = 2) + geom_hline(yintercept = 206, linetype = 2) +
    scale_fill_viridis() +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    theme(panel.background = element_blank(), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

  return(plot)
}

plot_heatmap_openfield_diff <- function(heatmap_array){

  limits <- max(abs(heatmap_array$z)) * c(-1, 1)

  plot <- ggplot(heatmap_array, aes(x, y, fill = z)) +
    geom_raster() +
    geom_vline(xintercept = 68, linetype = 2) + geom_vline(xintercept = 206, linetype = 2) +
    geom_hline(yintercept = 68, linetype = 2) + geom_hline(yintercept = 206, linetype = 2) +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", p1 = 0.4, p2 = 0.4) +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    theme(panel.background = element_blank(), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

  return(plot)
}

## BAR GRAPHS

bg_twogroup <- function(summary_array, metric, group_color, axis_limits, treatment_color){

  colnames(summary_array)[4] <- "metric"
  summary <- summary_array %>%
    group_by(group, time) %>%
    summarize(sem = sd(metric)/sqrt(n()), metric = mean(metric))

  ggplot(summary_array, aes(x = time, y = metric, fill = as.character(group), colour = as.character(group))) +
    geom_rect(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
              fill = c("white", treatment_color, "white", treatment_color), data = summary, alpha = 0.4, colour = NA) +
    geom_bar(stat = "summary", position = "dodge", fill = NA) +
    geom_errorbar(aes(ymin = metric - sem, ymax = metric + sem), data = summary, width = .5, position = position_dodge(width = 0.9)) +
    geom_line(aes(group = id)) +
    geom_point(aes(group = id), shape = 21, color = 'black') +
    scale_y_continuous(expand = c(0,0), limits = axis_limits) +
    scale_fill_manual(values = rev(c(control_color, group_color))) +
    scale_colour_manual(values = rev(c(control_color, group_color))) +
    geom_hline(yintercept = 0) +
    xlab("") +
    ylab(metric) +
    facet_wrap(~ group) +
    theme_classic() +
    theme(legend.position="none",
          axis.ticks.x=element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          strip.background = element_blank())

}

bg_twogroup_nobar <- function(summary_array, metric, group_color, axis_limits, treatment_color, chance_level){

  colnames(summary_array)[4] <- "metric"
  summary <- summary_array %>%
    group_by(group, time) %>%
    summarize(sem = sd(metric)/sqrt(n()), metric = mean(metric))

  ggplot(summary_array, aes(x = time, y = metric, fill = as.character(group), colour = as.character(group))) +
    geom_rect(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
              fill = c("white", treatment_color, "white", treatment_color), data = summary, alpha = 0.4, colour = NA) +
    geom_line(aes(group = id)) +
    geom_point(aes(group = id), shape = 21, color = 'black') +
    geom_errorbar(aes(ymin = metric - sem, ymax = metric + sem), data = summary, position = position_dodge(width = 0.9), colour = "black") +
    geom_errorbar(aes(ymin = metric, ymax = metric), data = summary, position = position_dodge(width = 0.9)) +
    #geom_point(stat = "summary", position = "dodge", colour = "black", shape = 21, size = 6) +
    scale_y_continuous(expand = c(0,0), limits = axis_limits) +
    scale_fill_manual(values = rev(c(control_color, group_color))) +
    scale_colour_manual(values = rev(c(control_color, group_color))) +
    geom_hline(yintercept = chance_level, linetype = 2) +
    xlab("") +
    ylab(metric) +
    facet_wrap(~ group) +
    theme_classic() +
    theme(legend.position="none",
          axis.ticks.x=element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          strip.background = element_blank())

}

bg_diff <- function(summary_array, metric, group_color, axis_limits, treatment_color){

  colnames(summary_array)[3] <- "metric"

  summary <- summary_array %>%
    group_by(group) %>%
    summarize(sem = sd(metric)/sqrt(n()), metric = mean(metric))

  ggplot(summary_array, aes(x = group, y = metric, fill = as.character(group), colour = as.character(group))) +
    geom_rect(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
              fill = c("white", treatment_color), data = summary, alpha = 0.4, colour = NA) +
    geom_rect(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
              fill = "white", data = summary, alpha = 0.4, colour = NA) +
    geom_bar(stat = "summary", position = "dodge", fill = NA) +
    geom_errorbar(aes(ymin = metric - sem, ymax = metric + sem), data = summary, width = .5, position = position_dodge(width = 0.9)) +
    geom_point(aes(group = id), shape = 21, color = 'black') +
    scale_y_continuous(expand = c(0,0), limits = axis_limits) +
    scale_fill_manual(values = c(group_color, control_color)) +
    scale_colour_manual(values = c(group_color, control_color)) +
    geom_hline(yintercept = 0) +
    xlab("") +
    ylab(metric) +
    theme_classic() +
    theme(legend.position="none",
          axis.ticks.x=element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          strip.background = element_blank())

}

bg_epm_prepost <- function(summary_array, metric, group_color, axis_limits, treatment_color){

  colnames(summary_array)[4] <- "metric"
  summary <- summary_array %>%
    group_by(group, time) %>%
    summarize(sem = sd(metric)/sqrt(n()), metric = mean(metric))

  ggplot(summary_array, aes(x = time, y = metric, fill = as.character(group))) +
    geom_rect(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
              fill = c("white", treatment_color, "white", treatment_color, "white", treatment_color), data = summary, alpha = 0.4) +
    geom_line(aes(group = group, colour = group), stat = "summary") +
    geom_errorbar(aes(ymin = metric - sem, ymax = metric + sem, colour = group), data = summary, width = .1) +
    geom_point(aes(group = group), shape = 21, size = 3, colour = 'black', stat = "summary") +
    scale_y_continuous(expand = c(0,0), limits = axis_limits) +
    scale_x_discrete(labels = c("0-5", "5-10", "10-15")) +
    scale_fill_manual(values = rev(c(control_color, group_color))) +
    scale_colour_manual(values = c(control_color, group_color)) +
    xlab("Time (min)") +
    ylab(metric) +
    theme_classic() +
    theme(legend.position="none",
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          strip.background = element_blank())

}

bg_epm_mins <- function(subset_array, metric, group_color, axis_limits, treatment_color){

  summary_array <- bind_rows(subset_array, .id = "group")

  colnames(summary_array)[4] <- "metric"
  summary <- summary_array %>%
    group_by(group, mins) %>%
    summarize(sem = sd(metric)/sqrt(n()), metric = mean(metric))

  ggplot(summary_array, aes(x = mins, y = metric, fill = as.character(group))) +
    geom_rect(xmin = 5, xmax = 10, ymin = -Inf, ymax = Inf, fill = rep(c(treatment_color, "white"), times = 15), data = summary, alpha = 0.1) +
    geom_line(aes(group = group, colour = group), stat = "summary") +
    geom_errorbar(aes(ymin = metric - sem, ymax = metric + sem, colour = group), data = summary, width = .5) +
    geom_point(aes(group = group), shape = 21, size = 2.7, colour = 'black', stat = "summary") +
    scale_y_continuous(expand = c(0,0), limits = axis_limits) +
    scale_x_continuous(expand = c(0,0), breaks = seq(0, 15, 5), limits = c(0, 16)) +
    scale_fill_manual(values = c(control_color, group_color)) +
    scale_colour_manual(values = c(control_color, group_color)) +
    xlab("Time (min)") +
    ylab(metric) +
    theme_classic() +
    theme(legend.position="none",
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          strip.background = element_blank())

}

bg_openfield_5x <- function(summary_array, metric, group_color, axis_limits, treatment_color){

  colnames(summary_array)[4] <- "metric"
  summary <- summary_array %>%
    group_by(group, time) %>%
    summarize(sem = sd(metric)/sqrt(n()), metric = mean(metric))

  ggplot(summary_array, aes(x = time, y = metric, fill = as.character(group))) +
    geom_rect(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
              fill = c("white", treatment_color, "white", treatment_color, "white", treatment_color,
                       "white", treatment_color, "white", treatment_color), data = summary, alpha = 0.4) +
    geom_rect(xmin = 3.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
              fill = c("white", treatment_color, "white", treatment_color, "white", treatment_color,
                       "white", treatment_color, "white", treatment_color), data = summary, alpha = 0.4) +
    geom_line(aes(group = group, colour = group), stat = "summary") +
    geom_errorbar(aes(ymin = metric - sem, ymax = metric + sem, colour = group), data = summary, width = .1) +
    geom_point(aes(group = group), shape = 21, size = 3, colour = 'black', stat = "summary") +
    scale_y_continuous(expand = c(0,0), limits = axis_limits) +
    scale_x_discrete(labels = c("0-5", "5-10", "10-15", "15-20", "20-25")) +
    scale_fill_manual(values = rev(c(control_color, group_color))) +
    scale_colour_manual(values = c(control_color, group_color)) +
    xlab("Time (min)") +
    ylab(metric) +
    theme_classic() +
    theme(legend.position="none",
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          strip.background = element_blank())

}

bg_oft_mins <- function(subset_array, metric, group_color, axis_limits, treatment_color){

  summary_array <- bind_rows(subset_array, .id = "group")

  colnames(summary_array)[4] <- "metric"
  summary <- summary_array %>%
    group_by(group, mins) %>%
    summarize(sem = sd(metric)/sqrt(n()), metric = mean(metric))

  ggplot(summary_array, aes(x = mins, y = metric, fill = as.character(group))) +
    geom_rect(xmin = 5, xmax = 10, ymin = -Inf, ymax = Inf,
              fill = rep(c(treatment_color, "white"), times = 25), data = summary, alpha = 0.4) +
    geom_rect(xmin = 15, xmax = 20, ymin = -Inf, ymax = Inf,
              fill = rep(c(treatment_color, "white"), times = 25), data = summary, alpha = 0.4) +
    geom_line(aes(group = group, colour = group), stat = "summary") +
    geom_errorbar(aes(ymin = metric - sem, ymax = metric + sem, colour = group), data = summary, width = .5) +
    geom_point(aes(group = group), shape = 21, size = 1.7, colour = 'black', stat = "summary") +
    scale_y_continuous(expand = c(0,0), limits = axis_limits) +
    scale_x_continuous(expand = c(0,0), breaks = seq(0, 25, 5), limits = c(0, 26)) +
    scale_fill_manual(values = c(control_color, group_color)) +
    scale_colour_manual(values = c(control_color, group_color)) +
    xlab("Time (min)") +
    ylab(metric) +
    theme_classic() +
    theme(legend.position="none",
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          strip.background = element_blank())

}

## PATHS

plot_path <- function(x, path_color) {

  ggplot(x$treatment, aes(x = XCoordinate, y = YCoordinate)) +
    #geom_rect(aes(xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0),
    #          fill = "#E3FCFF") +
    geom_path(colour = path_color) +
    geom_hline(yintercept = 0, linetype = 2, colour = chr2_blue) +
    geom_vline(xintercept = 0, linetype = 2, colour = chr2_blue) +
    theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
}

## STRIPES

discrete_stripe <- function(x, colname, colorlist){

  plot <- ggplot(x, aes(Time, 60)) + geom_tile(aes(fill = .data[[colname]])) +
                 scale_fill_manual(values = colorlist) + ylab(" ") +
                 geom_vline(xintercept = 600) +
                 theme(panel.background = element_blank(),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                       axis.title.y=element_text(size=30,face="bold"),
                       axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                       legend.position = "none") + scale_x_continuous(expand = c(0,0))
  return(plot)
}

continuous_stripe <- function(x, colname){

  plot <- ggplot(x, aes(Time, 60)) + geom_tile(aes(fill = .data[[colname]])) +
    scale_fill_gradientn(colours = matlab.like(10)) + ylab(" ") +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_text(size=30,face="bold"),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          legend.position = "none") + scale_x_continuous(expand = c(0,0))
  return(plot)
}

velocity_stripe <- function(x){

  x$Time <- (x$Time / 60)
  plot <- ggplot(x, aes(Time, Velocity)) +
                 geom_segment(x = x$Time, xend = x$Time, y = 0, yend = max(x$Velocity), aes(colour = x$Freezing)) +
                 geom_path() +
                 scale_colour_manual(values = c('lightcoral', 'white')) +
                 scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0), breaks = pretty_breaks()) +
                 theme(legend.position="none", panel.background = element_blank(),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black")) + xlab("Time (min)") + ylab("Velocity (cm/s)")
  return(plot)
}

## QUAD PLOTS (TO BE USED?)

quad_series <- function(x){

  quad <- cbind(x$Time,
                integer(length = length(rownames(x))),
                integer(length = length(rownames(x))),
                integer(length = length(rownames(x))),
                integer(length = length(rownames(x))))
  quad[,1] <- quad[,1] - quad[1,1] + .25 # remove any potential issues from start time if non-zero
  quad[,1] <- as.character(quad[,1])
  colnames(quad) <- c("Time", "LL", "LR", "UL", "UR")
  quad_avg <- quad
  quad_list <- c("LL", "LR", "UL", "UR")

  for (i in 2:5) {
    for (j in 1:length(rownames(x))){
      if (x[j, 4] == quad_list[i-1]){
        quad[j,i] <- 1
      }
    }
    quad_avg[,i] <- cummean(quad[,i])
    quad[,i] <- quad[,i]
    quad[,i] <- cumsum(quad[,i])
  }

  quad <- as.data.frame(quad) %>% pivot_longer(!Time)
  quad_avg <- as.data.frame(quad_avg) %>% pivot_longer(!Time)

  quad <- cbind(quad, quad_avg$value)
  colnames(quad) <- c("Time", "Quadrant", "Cumulative", "Average")
  quad$Time <- as.numeric(quad$Time)
  quad[,3:4] <- cbind(as.numeric(quad$Cumulative)*.25, as.numeric(quad$Average)) # format for graphs
  quad$Average <- quad$Average*100

  return(quad)
}

quad_plot <- function(x, colname){

  if (colname == "Cumulative"){
    y_label <- "Quadrant Occupancy (sec)"
  } else {
    y_label <- "Quadrant Occupancy (%)"
  }

  plot <- ggplot(x, aes(Time / 60, .data[[colname]], colour = Quadrant)) + geom_line() +
                 scale_colour_manual(values = c("cornflowerblue", "red", "blue", "blue4")) +
                 scale_y_continuous(expand = c(0,0)) +
                 scale_x_continuous(expand = c(0,0), breaks = pretty_breaks()) +
                 theme(legend.position="bottom", panel.background = element_blank(),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black"), legend.key = element_rect(fill = NA, color = NA)) +
                 xlab("Time (min)") + ylab(y_label)
  return(plot)
}

## SERIES PLOTS (TO BE USED?)

cont_series_plot <- function(x, colname, y_label){

  series_array <- cbind.data.frame(x$Time, cummean(x[[colname]]))
  colnames(series_array) <- c("Time", "Metric")

  plot <- ggplot(series_array, aes(Time / 60, Metric)) + geom_line() +
                 scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0), breaks = pretty_breaks()) +
                 theme(legend.position="none", panel.background = element_blank(),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black")) +
                 xlab("Time (min)") + ylab(y_label)
  return(plot)
}

pi_series_plot <- function(x){

  series_array <- cbind.data.frame(x$Time, x$Quadrant)
  colnames(series_array) <- c("Time", "PI")
  series_array$Time <- series_array$Time - series_array$Time[1] + .25 # remove any potential issues from start time if non-zero

  series_array$PI <- cummean(ifelse (series_array$PI == "LR", 1, 0)) # get cumulative means
  series_array$PI <- (series_array$PI - 0.25) / 0.25

  plot <- ggplot(series_array, aes(Time / 60, PI)) + geom_line() +
                 scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0), breaks = pretty_breaks()) +
                 theme(legend.position="none", panel.background = element_blank(),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       axis.line = element_line(colour = "black")) +
                 xlab("Time (min)") + ylab("Performance Index")
  return(plot)
}

## LEGENDS/PARAMETERS

points_legend <- function(summary_array, metric, group_color, axis_limits, treatment_color, direction){

  colnames(summary_array)[4] <- "metric"
  summary <- summary_array %>%
    group_by(group, time) %>%
    summarize(sem = sd(metric)/sqrt(n()), metric = mean(metric))

  plot <- ggplot(summary_array, aes(x = time, y = metric, fill = as.character(group), colour = as.character(group))) +
    geom_rect(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
              fill = c("white", treatment_color, "white", treatment_color), data = summary, alpha = 0.4, colour = NA) +
    geom_bar(stat = "summary", position = "dodge", fill = NA) +
    geom_errorbar(aes(ymin = metric - sem, ymax = metric + sem), data = summary, width = .5, position = position_dodge(width = 0.9)) +
    geom_line(aes(group = id)) +
    geom_point(aes(group = id), shape = 21, color = 'black') +
    scale_y_continuous(expand = c(0,0), limits = axis_limits) +
    scale_fill_manual(values = rev(c(control_color, group_color))) +
    scale_colour_manual(values = rev(c(control_color, group_color))) +
    geom_hline(yintercept = 0) +
    xlab("") +
    ylab(metric) +
    facet_wrap(~ group) +
    theme_classic() +
    theme(legend.direction = direction, legend.position = "top",
          axis.ticks.x=element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          strip.background = element_blank())

  plot <- cowplot::get_legend(plot)

  return(plot)
}

plot_heatmap_legend <- function(heatmap_array, direction){

  limits <- max(abs(heatmap_array$z)) * c(-1, 1)

  plot <- ggplot(heatmap_array, aes(x, y, fill = z)) +
    geom_raster() +
    #geom_rect(aes(xmin = 0, xmax = max(.data$x), ymin = 0, ymax = min(.data$y)),
    #    fill = "grey", alpha = 0.03) +
    #geom_segment(aes(x = 0, xend = -Inf, y = 295.5, yend = 295.5), linetype = 2, colour = chr2_blue) +
    #geom_segment(aes(y = 0, yend = -Inf, x = 295.5, xend = 295.5), linetype = 2, colour = chr2_blue) +
    geom_vline(xintercept = 295.5, linetype = 2, colour = chr2_blue) +
    geom_hline(yintercept = 295.5, linetype = 2, colour = chr2_blue) +
    #scale_fill_gradientn(colours = rev(brewer.pal(n = 100, name = "RdBu")), limit = limits) +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", p1 = 0.4, p2 = 0.4) +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    theme(panel.background = element_blank(), legend.direction = direction,
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

  plot <- cowplot::get_legend(plot)

  return(plot)
}

plot_heatmap_parameters <- function(heatmap_array, par1, par2, col_scale){

  limits <- max(abs(heatmap_array$z)) * c(-1, 1)

  plot <- ggplot(heatmap_array, aes(x, y, fill = z)) +
    geom_raster() +
    #geom_rect(aes(xmin = 0, xmax = max(.data$x), ymin = 0, ymax = min(.data$y)),
    #    fill = "grey", alpha = 0.03) +
    #geom_segment(aes(x = 0, xend = -Inf, y = 295.5, yend = 295.5), linetype = 2, colour = chr2_blue) +
    #geom_segment(aes(y = 0, yend = -Inf, x = 295.5, xend = 295.5), linetype = 2, colour = chr2_blue) +
    geom_vline(xintercept = 295.5, linetype = 2, colour = chr2_blue) +
    geom_hline(yintercept = 295.5, linetype = 2, colour = chr2_blue) +
    #scale_fill_gradientn(colours = rev(brewer.pal(n = 100, name = "RdBu")), limit = limits) +
    scale_fill_continuous_diverging(palette = col_scale, p1 = par1, p2 = par2) +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    theme(panel.background = element_blank(), legend.position = "none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

  return(plot)
}

plot_heatmap_openfield_legend <- function(heatmap_array, direction){

  limits <- max(abs(heatmap_array$z)) * c(-1, 1)

  plot <- ggplot(heatmap_array, aes(x, y, fill = z)) +
    geom_raster() +
    #geom_rect(aes(xmin = 0, xmax = max(.data$x), ymin = 0, ymax = min(.data$y)),
    #    fill = "grey", alpha = 0.03) +
    geom_vline(xintercept = 68, linetype = 2) + geom_vline(xintercept = 206, linetype = 2) +
    geom_hline(yintercept = 68, linetype = 2) + geom_hline(yintercept = 206, linetype = 2) +
    scale_fill_viridis() +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.direction = direction,
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

  plot <- cowplot::get_legend(plot)

  return(plot)
}

plot_heatmap_openfield_legend_diff <- function(heatmap_array, direction){

  limits <- max(abs(heatmap_array$z)) * c(-1, 1)

  plot <- ggplot(heatmap_array, aes(x, y, fill = z)) +
    geom_raster() +
    #geom_rect(aes(xmin = 0, xmax = max(.data$x), ymin = 0, ymax = min(.data$y)),
    #    fill = "grey", alpha = 0.03) +
    geom_vline(xintercept = 68, linetype = 2) + geom_vline(xintercept = 206, linetype = 2) +
    geom_hline(yintercept = 68, linetype = 2) + geom_hline(yintercept = 206, linetype = 2) +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", p1 = 0.4, p2 = 0.4) +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.direction = direction,
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

  plot <- cowplot::get_legend(plot)

  return(plot)
}

plot_heatmap_epm_legend <- function(heatmap_array, direction){

  limits <- max(abs(heatmap_array$z)) * c(-1, 1)

  plot <- ggplot(heatmap_array, aes(x, y, fill = z)) +
    geom_raster() +
    annotate("rect", xmin = 0, xmax = 345, ymin = 0, ymax = 345, fill = "white", colour = NA, size = 0.5) +
    annotate("rect", xmin = 395, xmax = 741, ymin = 0, ymax = 345, fill = "white", colour = NA, size = 0.5) +
    annotate("rect", xmin = 0, xmax = 345, ymin = 395, ymax = 741, fill = "white", colour = NA, size = 0.5) +
    annotate("rect", xmin = 395, xmax = 741, ymin = 395, ymax = 741, fill = "white", colour = NA, size = 0.5) +
    annotate("rect", xmin = 345, xmax = 395, ymin = 0, ymax = 345, colour = "blue", fill = NA, size = 0.5) +
    annotate("rect", xmin = 345, xmax = 395, ymin = 395, ymax = 741, colour = "blue", fill = NA, size = 0.5) +
    annotate("rect", ymin = 345, ymax = 395, xmin = 0, xmax = 345, colour = "red", fill = NA, size = 0.5) +
    annotate("rect", ymin = 345, ymax = 395, xmin = 395, xmax = 741, colour = "red", fill = NA, size = 0.5) +
    scale_fill_viridis() +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.direction = direction,
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

  plot <- cowplot::get_legend(plot)

  return(plot)
}

plot_heatmap_epm_legend_diff <- function(heatmap_array, direction){

  limits <- max(abs(heatmap_array$z)) * c(-1, 1)

  plot <- ggplot(heatmap_array, aes(x, y, fill = z)) +
    geom_raster() +
    annotate("rect", xmin = 0, xmax = 345, ymin = 0, ymax = 345, fill = "white", colour = NA, size = 0.5) +
    annotate("rect", xmin = 395, xmax = 741, ymin = 0, ymax = 345, fill = "white", colour = NA, size = 0.5) +
    annotate("rect", xmin = 0, xmax = 345, ymin = 395, ymax = 741, fill = "white", colour = NA, size = 0.5) +
    annotate("rect", xmin = 395, xmax = 741, ymin = 395, ymax = 741, fill = "white", colour = NA, size = 0.5) +
    annotate("rect", xmin = 345, xmax = 395, ymin = 0, ymax = 345, colour = "purple", fill = NA, size = 0.5) +
    annotate("rect", xmin = 345, xmax = 395, ymin = 395, ymax = 741, colour = "purple", fill = NA, size = 0.5) +
    annotate("rect", ymin = 345, ymax = 395, xmin = 0, xmax = 345, colour = "darkgreen", fill = NA, size = 0.5) +
    annotate("rect", ymin = 345, ymax = 395, xmin = 395, xmax = 741, colour = "darkgreen", fill = NA, size = 0.5) +
    scale_fill_continuous_diverging(palette = "Blue-Red 3", p1 = 0.4, p2 = 0.4) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.direction = direction,
          axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

  plot <- cowplot::get_legend(plot)

  return(plot)
}
