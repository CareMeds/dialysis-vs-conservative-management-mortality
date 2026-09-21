create_ps_distribution_plot <- function(data,
                                        PS_varname,
                                        trt_varname,
                                        weights = NULL,
                                        PS_title_hist,
                                        PS_title_scaled_hist,
                                        titleSize = 22,
                                        TextSize = 26,
                                        xlab = TRUE,
                                        ylab = TRUE,
                                        palette = c("blue", "darkgreen"),
                                        x_axis_text = "Propensity score") {
  # extract ps and trt
  data <- copy(data)
  data$ps <- data[, PS_varname, with = FALSE][[1]]
  data$trt <- data[, trt_varname, with = FALSE][[1]]
  
  # calculate unadjusted AUC
  AUC <- pROC::auc(pROC::roc(
    predictor = data$ps,
    response = data$trt,
    quiet = TRUE
  ))
  
  # create histogram for propensity scores
  # !! added to make sure to use weights from data
  hist <- ggplot2::ggplot(data, ggplot2::aes(
    x = ps,
    fill = as.factor(trt),
    weight = !!weights
  )) +
    ggplot2::geom_histogram(
      alpha = 0.6,
      position = "dodge",
      binwidth = 0.01,
      boundary = 0
    ) +
    ggplot2::scale_x_continuous(x_axis_text, limits = c(0, 1)) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::theme(
      text = ggplot2::element_text(size = TextSize),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = titleSize
      )
    ) +
    ggplot2::ggtitle(PS_title_hist)
  
  # create scaled density for propensity scores
  scaled_hist <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = ps,
      fill = as.factor(trt),
      weight = !!weights,
      y = ggplot2::after_stat(scaled)
    )
  ) +
    ggplot2::geom_density(alpha = 0.6, bw = 0.05) +
    ggplot2::scale_x_continuous(x_axis_text, limits = c(0, 1)) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::theme(
      text = ggplot2::element_text(size = TextSize),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = titleSize
      )
    ) +
    ggplot2::ggtitle(PS_title_scaled_hist)
  
  if (!xlab) {
    hist <- hist +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
    scaled_hist <- scaled_hist +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
  }
  
  if (!ylab) {
    hist <- hist +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
    scaled_hist <- scaled_hist +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank()
      )
  }
  
  return(list(
    hist = hist,
    scaled_hist = scaled_hist,
    AUC = AUC
  ))
}

# create love plot to compare SMD before and after weighting
love_plot <- function(SMDs_dt,
                      SMD_names,
                      plotColors,
                      xlab_title = "Standardized mean difference",
                      xmax,
                      titleSize = 22,
                      TextSize = 26) {
  # make one love plot
  love.plot <- ggplot2::ggplot(data = SMDs_dt, ggplot2::aes(y = factor(rownames(SMDs_dt), levels = rev(
    rownames(SMDs_dt)
  )))) +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid") +
    ggplot2::geom_vline(xintercept = 0.1, linetype = "dashed") +
    ggplot2::geom_point(x = SMDs_dt[, SMD_names[1]],
                        colour = plotColors[1],
                        size = 4) +
    ggplot2::geom_point(x = SMDs_dt[, SMD_names[2]],
                        colour = plotColors[2],
                        size = 4) +
    ggplot2::xlab(xlab_title) +
    ggplot2::ylab("") +
    ggplot2::scale_x_continuous(limits = c(0, max(max(SMDs_dt), xmax)), breaks =
                                  seq(0, max(max(SMDs_dt), xmax), 0.1)) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = titleSize
      ),
      text = ggplot2::element_text(size = TextSize),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank()
    )
  
  if (length(SMD_names) == 4) {
    love.plot <- love.plot +
      ggplot2::geom_point(x = SMDs_dt[, SMD_names[3]],
                          colour = plotColors[1],
                          size = 4) +
      ggplot2::geom_point(x = SMDs_dt[, SMD_names[4]],
                          colour = plotColors[2],
                          size = 4)
  }
  return(love.plot)
}

# create Kaplan-Meier plot with effect measures
create_KM_plot <- function(data,
                           trt_var,
                           event_var,
                           time2event_var,
                           w_meth,
                           out_est, 
                           horizon,
                           unit = "months",
                           manual_colors,
                           trt_labels = c("trt=0", "trt=1")) {
  # extract data from out_est
  plot_data <- data.table(
    time = out_est$est_full$time,
    strata = factor(c(rep(0, length(out_est$est_full$time)), 
                      rep(1, length(out_est$est_full$time))), 
                    levels = c(0, 1)),
    surv = c(1 - out_est$est_full$R0, 1 - out_est$est_full$R1),
    lower = c(
      1 - out_est$est_CI_full |> 
        dplyr::filter(name == "conf.low") |> 
        dplyr::pull(R0),
      1 - out_est$est_CI_full |>
        dplyr::filter(name == "conf.low") |> 
        dplyr::pull(R1)
    ),
    upper = c(
      1 - out_est$est_CI_full |> 
        dplyr::filter(name == "conf.high") |> 
        dplyr::pull(R0),
      1 - out_est$est_CI_full |> 
        dplyr::filter(name == "conf.high") |> 
        dplyr::pull(R1)
    )
  )
  
  # get location for censoring, and risk table
  model_formula <- as.formula(paste0(
    "survival::Surv(",
    time2event_var,
    ", ",
    event_var,
    ") ~ ",
    trt_var
  ))
  
  # define weights
  if (w_meth=="unweighted"){
    weights <- rep(1, nrow(data))
  } else{
    weights <- data[[paste0("sw_", w_meth)]]
  }
  
  # Create a list of arguments for the survfit call
  keep <- weights > 0
  fit_args <- list(
    formula = model_formula,
    data = data[keep],
    robust = TRUE,
    weights = weights[keep]
  )
  
  # Use do.call to ensure the function call is constructed correctly
  KM_fit <- do.call(survival::survfit, fit_args)
  KM_curve <- summary(KM_fit, times = unique(plot_data$time))
  
  # censoring data
  censor_dt <- data.table(
    time = KM_curve$time,
    strata = factor(as.numeric(KM_curve$strata)-1, 
                    levels = c(0, 1)),
    surv = KM_curve$surv,
    n.censor = KM_curve$n.censor
  )
  censor_data <- censor_dt[n.censor>0]
  
  # create plot
  KM_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = time,
      y = surv,
      color = strata,
      fill = strata,
      linetype = strata
    )
  ) +
    ggplot2::geom_line(linewidth = 1) +
    pammtools::geom_stepribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.2, # Transparency for the shading
      color = NA   # Remove the outline from the ribbon itself
    ) +
    ggplot2::geom_point(data = censor_data,  # add censoring
                        ggplot2::aes(x = time, y = surv),
                        shape = 3, 
                        size = 2,
                        stroke = 0.8,
                        show.legend = FALSE) +
    ggplot2::labs(
      x = paste0("Time (", unit, ")"),
      y = "Survival probability (%)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black"),
      axis.ticks = ggplot2::element_line(color = "black"),
      text = ggplot2::element_text(size = 14),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.1),
      labels = seq(0, 100, by = 10),
      expand = c(0, 0)
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, (horizon+1), by = 365 / 2), # break every 6 months
      labels = function(x)
        round(x / 30.5) # days -> months
    ) +
    ggplot2::scale_color_manual(labels = trt_labels,
                                values = manual_colors) +
    ggplot2::scale_fill_manual(guide = "none",
                               values = manual_colors) +
    ggplot2::scale_linetype_manual(guide = "none", 
                                   values = c("solid", "solid"))
  
  # find numbers at risk at each time point of interest
  KM_table <- summary(KM_fit, times = seq(0, horizon, by = 365 / 2))
  risk_table <- data.table(
    time = KM_table$time / 365,
    strata = factor(as.numeric(KM_table$strata)-1, 
                    levels = c(0, 1)),
    n.risk = round(KM_table$n.risk)  # after weighting might not be integer
  )
  
  # create colors labels using HTML
  colored_labels <- paste0("<span style='color:", manual_colors[1:2], "'>", trt_labels, "</span>")
  
  # create table
  KM_table <- ggplot2::ggplot(risk_table, 
                              ggplot2::aes(x = time, 
                                           y = strata, 
                                           label = n.risk)) +
    ggplot2::geom_text() +
    ggplot2::scale_x_continuous(limits = c(0, horizon / 365),
                                breaks = seq(0, horizon / 365, by = 0.5)) +
    ggplot2::scale_y_discrete(labels = colored_labels) + 
    ggplot2::theme_void() +
    ggplot2::theme(
      axis.text.y = ggtext::element_markdown(hjust = 1)
    )
  
  return(list(KM_plot = KM_plot,
              KM_table = KM_table))
}

# combine histograms on propensity score distributions
combine_PS_plots <- function(plot_list) {
  lapply(seq_along(plot_list), function(i) {
    plot_list[[i]]$hist + plot_list[[i]]$scaled_hist
  })
}

create_forest_plot_all_measures <- function(dt,
                                            print_metrics = c("N", "Nonoverlap", "Imbalance")) {
  # Loop over each measure to generate forest plots
  plot_list <- list()
  for (measure in c("RD", "dRMST", "HR")) {
    # create label
    dt[, paste0(measure, "_label") :=
         fmt_ci(get(measure),
                get(paste0(measure, "_lower")),
                get(paste0(measure, "_upper")),
                digits = ifelse(measure == "RD", 1, 2))]
    
    # fix the order
    dt$analysis_name <- factor(dt$analysis_name, levels = rev(dt$analysis_name))
    
    #------------------------------------------------------------
    # Base header
    #------------------------------------------------------------
    header_table <- ggplot2::ggplot(data.frame(y = 0), ggplot2::aes(y = y)) +
      ggplot2::geom_text(
        x = 0,
        label = "Subgroup",
        hjust = 0,
        vjust = 0,
        fontface = "bold"
      ) +
      ggplot2::geom_text(
        x = 1,
        label = ifelse(
          measure == "RD",
          "2-year RD\n% (95% CI)",
          ifelse(
            measure == "dRMST",
            "2-year \u0394RMST\nmonths (95% CI)",
            "2-year HR\n(95% CI)"
          )
        ),
        hjust = 1,
        vjust = 0,
        fontface = "bold"
      ) +
      ggplot2::theme_void() +
      ggplot2::scale_y_continuous(limits = c(0, 1))
    
    
    #------------------------------------------------------------
    # Base table
    #------------------------------------------------------------
    table <- ggplot2::ggplot(dt, ggplot2::aes(y = analysis_name)) +
      ggplot2::geom_text(ggplot2::aes(x = 0, label = analysis_name), hjust = 0) +
      ggplot2::geom_text(ggplot2::aes(x = 1, label = .data[[paste0(measure, "_label")]]), hjust = 1) +
      ggplot2::theme_void() +
      ggplot2::geom_segment(
        ggplot2::aes(
          x = 0,
          xend = 1,
          y = as.numeric(analysis_name) + 0.45,
          yend = as.numeric(analysis_name) + 0.45
        ),
        color = "gray90"
      )
    
    
    #------------------------------------------------------------
    # Dynamic metric columns
    #------------------------------------------------------------
    metric_positions <- get_metric_positions(length(print_metrics))
    names(metric_positions) <- print_metrics
    metric_labels    <- c(N = "Patients,\nNo.",
                          Nonoverlap = "Non overlap,\nNo.",
                          Imbalance = "SMD > 0.1,\nNo.")
    
    for (m in print_metrics) {
      # Header label
      header_table <- header_table +
        ggplot2::geom_text(
          x = metric_positions[m],
          label = metric_labels[m],
          hjust = 0,
          vjust = 0,
          fontface = "bold"
        )
      
      # Table column
      table <- table +
        ggplot2::geom_text(
          data = dt,
          mapping = ggplot2::aes(y = analysis_name, label = .data[[m]]),
          x = metric_positions[m],
          hjust = 0,
          inherit.aes = FALSE
        )
    }
    
    # header plot
    header_forest <- ggplot2::ggplot(data.frame(y = 0), ggplot2::aes(y = y)) +
      ggplot2::geom_text(
        x = 0.4,
        label = ifelse(measure == "dRMST", "Favor conservative", "Favor dialysis"),
        hjust = 1,
        vjust = 0,
        fontface = "bold"
      ) +
      ggplot2::geom_text(
        x = 0.6,
        label = ifelse(measure == "dRMST", "Favor dialysis", "Favor conservative"),
        hjust = 0,
        vjust = 0,
        fontface = "bold"
      ) +
      ggplot2::theme_void() +
      ggplot2::scale_y_continuous(limits = c(0, 1))
    
    # forest plot
    if (measure == "HR") {
      center <- 1
    } else {
      center <- 0  # default fallback
    }
    
    x_dev <- max(abs(dt[[paste0(measure, "_lower")]] - center), abs(dt[[paste0(measure, "_upper")]] - center), na.rm = TRUE)
    x_min <- ifelse(measure == "RD", -75, ifelse(measure == "dRMST", -14, 0))
    x_max <- ifelse(measure == "RD", 75, ifelse(measure == "dRMST", 14, 2))
    x_break <- ifelse(measure == "RD", 10, ifelse(measure == "dRMST", 2, 0.2))
    forest <- ggplot2::ggplot(dt, ggplot2::aes(y = analysis_name)) +
      ggplot2::geom_point(ggplot2::aes(x = .data[[measure]]),
                          shape = 15,
                          size = 2) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          xmin = .data[[paste0(measure, "_lower")]], 
          xmax = .data[[paste0(measure, "_upper")]],
          y = analysis_name
        ),
        width = 0,
        orientation = "y"
      ) + 
      ggplot2::geom_vline(
        xintercept = center,
        linetype = "dashed",
        color = "gray50"
      ) +
      # x-axis line in padding area
      ggplot2::geom_hline(
        yintercept = 0.5,
        # below the first row (padding area)
        color = "black"
      ) +
      ggplot2::coord_cartesian(clip = "off") +  # allow drawing in the padding area
      ggplot2::scale_x_continuous(
        limits = c(x_min, x_max),
        breaks = seq(x_min, x_max, x_break),
        labels = seq(x_min, x_max, x_break)
      ) +
      ggplot2::theme(
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(0, 0, -10, 0)
      )
    
    plot_list[[paste0("header_table_", measure)]] <- header_table
    plot_list[[paste0("table_", measure)]] <- table
    plot_list[[paste0("header_forest_", measure)]] <- header_forest
    plot_list[[paste0("forest_", measure)]] <- forest
  }
  
  # Save plot to results folder
  blank_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  combined_plot <- cowplot::plot_grid(
    plot_list$header_table_RD,
    plot_list$header_forest_RD,
    plot_list$table_RD,
    plot_list$forest_RD,
    plot_list$header_table_dRMST,
    plot_list$header_forest_dRMST,
    plot_list$table_dRMST,
    plot_list$forest_dRMST,
    plot_list$header_table_HR,
    plot_list$header_forest_HR,
    plot_list$table_HR,
    plot_list$forest_HR,
    blank_plot,
    blank_plot,
    ncol = 2,
    nrow = 7,
    rel_heights = c(0.4, 1, 0.4, 1, 0.4, 1, 0.1),
    rel_widths = c(1, 0.5)
  ) +
    ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA))
  
  return(list(combined_plot = combined_plot, dt = dt))
}

# metric positions of forest plots
get_metric_positions <- function(n) {
  if (n == 1) {
    return(c(0.50))
  } else if (n == 2) {
    return(c(0.40, 0.60))
  } else if (n == 3) {
    return(c(0.37, 0.49, 0.65))
  } else {
    stop("Only 1, 2, or 3 metrics supported.")
  }
}

# create effect plot of benefit versus risk
effect_plot <- function(estimates_df,
                        effect_modifier,
                        y_middle,
                        measure,
                        y_min_RD = -70,
                        y_max_RD,
                        y_min_RR = -1,
                        y_max_RR,
                        y_min_dRMST,
                        y_max_dRMST = 10,
                        y_min_HR = 0,
                        y_max_HR,
                        show_favor_annotation = TRUE) {
  # add padding on y-axis
  padding_y <- ifelse(measure == "HR" | measure == "RR", 0.1, ifelse(measure == "RD", 5, 1))
  
  # Define y-axis breaks based on measure
  if (measure == "HR") {
    y_min <- y_min_HR
    y_max <- y_max_HR
    y_breaks <- sort(c(seq(y_min, y_max, by = 0.2), y_middle))
  } else if (measure == "RD") {
    y_min <- y_min_RD
    y_max <- y_max_RD
    y_breaks <- sort(c(seq(y_min, y_max, by = 10), y_middle))
  } else if (measure == "RR") {
    y_min <- y_min_RR
    y_max <- y_max_RR
    y_breaks <- sort(c(seq(y_min, y_max, by = 0.2), y_middle))
  } else {
    y_min <- y_min_dRMST
    y_max <- y_max_dRMST
    y_breaks <- sort(c(seq(y_min, y_max, by = 1), y_middle))
  }
  
  # x-axis minimum
  if (min(estimates_df$effect_modifier_range) < 1) {
    x_min <- 0
    x_min_text <- 0.05
  } else{
    x_min <- 60
    x_min_text <- 61.5
  }
  
  # label describing the effect measure (now used as plot title, not y-axis label)
  # dRMST uses a plotmath expression so the delta (Δ) symbol always renders,
  # regardless of font/device (unicode escapes can silently drop on some devices)
  measure_label <- if (measure == "RD") {
    "Risk difference in %"
  } else if (measure == "dRMST") {
    bquote(bold(Delta * "RMST in months"))
  } else if (measure == "RR") {
    "Risk ratio"
  } else {
    "Hazard ratio"
  }
  
  # effect plot
  plot <- ggplot2::ggplot(estimates_df,
                          ggplot2::aes(x = effect_modifier_range, y = get(measure))) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = get(paste0(measure, "_lower")),
                                      ymax = get(paste0(measure, "_upper"))),
                         alpha = 0.3) +
    ggplot2::geom_hline(
      yintercept = y_middle,
      linetype = "dashed",
      color = "black"
    ) +
    ggplot2::scale_y_continuous(limits = c(y_min, y_max), breaks = y_breaks) +
    ggplot2::labs(x = ifelse(effect_modifier == "age", "Age in years", "Predicted 2-year mortality risk (%)"),
                  y = NULL) +
    ggplot2::ggtitle(measure_label) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      text = ggplot2::element_text(size = 18),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.ticks.x = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.5)
    )
  
  # add "favor dialysis" arrow + text indicator (optional)
  if (show_favor_annotation) {
    plot <- plot +
      ggplot2::annotate(
        "text",
        x = x_min_text,
        y = ifelse(measure == "dRMST", y_middle + padding_y, y_middle - padding_y),
        label = "Favor dialysis",
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 6
      ) +
      ggplot2::annotate(
        "segment",
        x = x_min,
        xend = x_min,
        y = ifelse(measure == "dRMST", y_middle + padding_y, y_middle - padding_y),
        yend = ifelse(measure == "dRMST", y_max, y_min),
        arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
        color = "black"
      )
  }
  
  # reverse y-axis for dRMST
  if (measure == "dRMST") {
    plot <- plot +
      ggplot2::scale_y_reverse(limits = c(y_max, y_min),
                               breaks = rev(y_breaks))
  }
  
  return(plot)
}

# histogram of effect modifier
create_histogram_stratified <- function(dt, var_name, trt_name, manual_colors) {
  # histogram plot
  hist_stratified <- ggplot2::ggplot(dt, ggplot2::aes(x = get(var_name), fill = as.factor(get(trt_name)))) +
    ggplot2::geom_histogram(
      binwidth = ifelse(var_name == "pred_risk", 0.01, 1),
      color = NA,
      position = "dodge"
    ) +
    ggplot2::labs(
      x = ifelse(
        var_name == "pred_risk",
        "Predicted 2-year mortality risk (%)",
        "Age in years"
      ),
      y = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = c("0" = manual_colors[1], "1" = manual_colors[2]),
      labels = c("Conservative", "Dialysis")
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid = ggplot2::element_blank(),
      text = ggplot2::element_text(size = 18),
      axis.ticks.x = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.ticks.y = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.5),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.direction = "horizontal"
    )
  
  # histogram x-axis limits
  if (var_name == "age") {
    hist_stratified <- hist_stratified +
      ggplot2::scale_x_continuous(breaks = seq(60, 100, 5),
                                  labels = seq(60, 100, 5)) +
      ggplot2::coord_cartesian(xlim = c(60, 100))
  } else if (var_name == "pred_risk") {
    hist_stratified <- hist_stratified +
      ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.1),
                                  labels = seq(0, 100, 10)) + # scales::percent_format(accuracy = 1)) +
      ggplot2::coord_cartesian(xlim = c(0, 1))
  }
  
  return(hist_stratified)
}

calibration_plot <- function(out_measures){
  # compute calibration plot data
  calibration_data <- data.frame(
    risk = out_measures$pseudos$risk,
    observed = out_measures$smooth_pseudos$fit,
    lower = out_measures$smooth_pseudos$fit - 1.96 * out_measures$smooth_pseudos$se.fit,
    upper = out_measures$smooth_pseudos$fit + 1.96 * out_measures$smooth_pseudos$se.fit
  )
  calibration_data$lower <- pmax(0, calibration_data$lower)
  calibration_data$upper <- pmin(1, calibration_data$upper)
  
  # create plot
  cal_plot <- ggplot2::ggplot(calibration_data, 
                              ggplot2::aes(x = risk, y = observed)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = "steelblue",
      alpha = 0.2
    ) +
    ggplot2::geom_line(linewidth = 0.75, alpha = 0.8) +
    ggplot2::annotate(
      "segment",
      x = 0,
      y = 0,
      xend = 1,
      yend = 1,
      linetype = "dashed",
      color = "gray40"
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(y = "Observed Risks") +
    ggthemes::theme_clean() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.background = ggplot2::element_rect(colour = NA),
      legend.position = "bottom",
      plot.subtitle = ggplot2::element_text(size = 10),
      panel.border = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank()
    )
  
  # histogram of predicted risk
  cal_hist <- ggplot2::ggplot(calibration_data, ggplot2::aes(x = risk)) +
    ggplot2::geom_histogram(
      binwidth = 0.01,
      fill = "steelblue",
      color = "white"
    ) +
    ggthemes::theme_clean() +
    ggplot2::coord_cartesian(c(0, 1)) +
    ggplot2::labs(title = NULL, y = "Count", x = "Predicted two-year mortality risk") +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      plot.background = ggplot2::element_blank()
    )
  
  return(list(cal_plot = cal_plot,
              cal_hist = cal_hist))
}

# --- Helper: extract named metrics from compute_measures() output --------
extract_metrics <- function(m) {
  c(Intercept = m$Intercept, Slope = m$Slope, AUC = m$AUC)
}

# --- Helpers to subset metrics by cohort prefix -------------------------
cohort_apparent <- function(prefix) {
  c(Intercept = apparent[[paste0(prefix, "_Intercept")]],
    Slope     = apparent[[paste0(prefix, "_Slope")]],
    AUC       = apparent[[paste0(prefix, "_AUC")]])
}

cohort_boots <- function(prefix) {
  data.frame(
    Intercept = orig_boots[[paste0(prefix, "_Intercept")]],
    Slope     = orig_boots[[paste0(prefix, "_Slope")]],
    AUC       = orig_boots[[paste0(prefix, "_AUC")]]
  )
}

cohort_optimism <- function(prefix) {
  c(Intercept = optimism[[paste0(prefix, "_Intercept")]],
    Slope     = optimism[[paste0(prefix, "_Slope")]],
    AUC       = optimism[[paste0(prefix, "_AUC")]])
}

# --- Annotation builder -------------------------------------------------
build_annotation <- function(apparent_vals, boot_vals, optimism_vals) {
  corr  <- function(metric) apparent_vals[metric] - optimism_vals[metric]
  ci_lo <- function(metric) quantile(boot_vals[[metric]], 0.025) - optimism_vals[metric]
  ci_hi <- function(metric) quantile(boot_vals[[metric]], 0.975) - optimism_vals[metric]
  
  paste0(
    "Calibration intercept ", fmt_ci(corr("Intercept"), ci_lo("Intercept"), ci_hi("Intercept"), digits = 3),
    "\nCalibration slope ",   fmt_ci(corr("Slope"),     ci_lo("Slope"),     ci_hi("Slope"), digits = 3),
    "\nAUC ",                 fmt_ci(corr("AUC"),       ci_lo("AUC"),       ci_hi("AUC"), digits = 3)
  )
}

annotate_cal_plot <- function(cal_plot_obj, apparent_vals, boot_vals, optimism_vals) {
  cal_plot_obj$cal_plot +
    ggplot2::annotate(
      "text",
      x     = 0.01,
      y     = 0.975,
      label = build_annotation(apparent_vals, boot_vals, optimism_vals),
      hjust = 0,
      vjust = 1,
      size  = 3
    )
}

save_cal_plot <- function(cal_plot_obj, annotated_plot, filename) {
  ggplot2::ggsave(
    plot     = annotated_plot / cal_plot_obj$cal_hist +
      patchwork::plot_layout(heights = c(3, 1)),
    filename = file.path(results_path, "Supplemental", filename),
    width    = 5,
    height   = 5,
    dpi      = 300
  )
}

plot_metric <- function(df,
                        y_var,
                        y_lab,
                        x_lab = "",
                        x_lim = c(0, 1),
                        y_lim = c(0, 1)) {
  if (y_var == "dRMST"){
    y_lab <- "\u0394RMST (months)" 
  }
  
  if (y_var == "RD"){
    y_lab <- "RD (%)" 
  }
  
  plot <- ggplot2::ggplot(df, ggplot2::aes(x = mortality_risk, y = .data[[y_var]])) +
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
    ggplot2::scale_x_continuous(limits = c(x_lim[1]-0.05, x_lim[2]+0.05), 
                                breaks = seq(x_lim[1], x_lim[2], 0.1)) +
    ggplot2::labs(y = y_lab, x = x_lab) +
    ggplot2::theme_classic()
  
  if (y_var == "RD"){
    plot <- plot + ggplot2::scale_y_continuous(
      limits = c(y_lim[1], y_lim[2]),
      breaks = seq(y_lim[1], y_lim[2], 0.1),
      labels = function(x) paste0(x * 100)
    )
  } else{
    plot <- plot + ggplot2::ylim(y_lim)
  }
  
  return(plot)
}