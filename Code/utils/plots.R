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
        round(x / 30.44) # days -> months
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

# Horizontal stacked bars, one per treatment arm, of the mean time (months) spent in
# each component over the follow-up (e.g. total time at home, total time in-center,
# death = the remaining time), so every bar adds up to total_months. Numbers (estimate
# and 95% CI) are written inside the segments; segments narrower than small_threshold
# months get their label outside the bar, with a short leader line. Dashed lines connect
# the segment boundaries between neighbouring bars.
# bar_dt: one row per arm x component with columns arm (label), arm_order (1 = top bar),
# component, value, lower, upper (all in months). component_colors: named vector whose
# names are the components, in the order they are stacked (left to right).
# Optional column connector_group: components that share a group (e.g. two parts of the
# same segment in one bar) are treated as one segment for the dashed connectors, so bars
# that split a segment can still be connected to bars that do not. Defaults to component.
# A component may be absent from some arms.
make_time_bar_figure <- function(bar_dt,
                                 component_colors,
                                 total_months = 24,
                                 small_threshold = 1.8,
                                 text_size = 3.3,
                                 bar_height = 0.4,
                                 legend_title = "Mean time over 2 years",
                                 x_title = "Months") {
  dt <- data.table::copy(data.table::as.data.table(bar_dt))
  dt[, component := factor(component, levels = names(component_colors))]
  data.table::setorder(dt, arm_order, component)
  
  if (!"connector_group" %in% names(dt)) dt[, connector_group := as.character(component)]
  
  # segment limits within each bar
  dt[, xmax := cumsum(value), by = arm_order]
  dt[, xmin := xmax - value]
  dt[, xmid := (xmin + xmax) / 2]
  
  # bar positions: arm_order 1 on top
  n_arm <- data.table::uniqueN(dt$arm_order)
  dt[, y := n_arm - arm_order + 1]
  h <- bar_height / 2
  
  # label text; dark text on light fills, white on dark fills
  dt[, label := sprintf("%.1f\n(%.1f, %.1f)", value, lower, upper)]
  lum <- function(col) {
    rgb <- grDevices::col2rgb(col) / 255
    0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
  }
  dt[, text_col := ifelse(lum(component_colors[as.character(component)]) < 0.55, "white", "black")]
  
  # labels inside the segment, or outside with a leader line for narrow segments:
  # above the top bar, below the others
  dt[, outside := value < small_threshold]
  dt[, direction := ifelse(arm_order == min(arm_order), 1, -1)]
  inside_dt  <- dt[outside == FALSE]
  outside_dt <- dt[outside == TRUE]
  outside_dt[, `:=`(y_edge  = y + direction * h,
                    y_label = y + direction * (h + 0.30))]
  
  # dashed connectors between the boundaries of neighbouring bars
  # (per connector_group: the right edge of the group's last segment in each bar; the last
  # group of a bar ends at total_months and gets no connector)
  group_order <- unique(dt[order(as.integer(component)), connector_group])
  boundaries <- dt[, .(boundary = max(xmax), y = y[1],
                       last_comp = max(as.integer(component))), by = .(arm_order, connector_group)]
  boundaries <- boundaries[connector_group != group_order[length(group_order)]]
  connect_dt <- merge(boundaries[, .(arm_order, y, boundary, connector_group)],
                      data.table::copy(boundaries[, .(arm_order, y, boundary, connector_group)])[, arm_order := arm_order - 1],
                      by = c("arm_order", "connector_group"),
                      suffixes = c("_top", "_bottom"))
  
  arms <- unique(dt[, .(y, arm)])
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = connect_dt,
      ggplot2::aes(x = boundary_top, xend = boundary_bottom,
                   y = y_top - h, yend = y_bottom + h),
      linetype = "dashed", colour = "grey60", linewidth = 0.3
    ) +
    ggplot2::geom_rect(
      data = dt,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = y - h, ymax = y + h, fill = component),
      colour = "white", linewidth = 0.4
    ) +
    ggplot2::geom_text(
      data = inside_dt,
      ggplot2::aes(x = xmid, y = y, label = label, colour = text_col),
      size = text_size, lineheight = 0.95
    ) +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_fill_manual(values = component_colors, name = legend_title) +
    # coord_cartesian (not scale limits) for the x range: a bar ending at total_months
    # plus a floating-point hair would be dropped as out of bounds, leaving it blank
    ggplot2::scale_x_continuous(breaks = seq(0, total_months, by = 3)) +
    ggplot2::scale_y_continuous(breaks = arms$y, labels = arms$arm) +
    ggplot2::coord_cartesian(xlim = c(-0.1, total_months + 0.5),
                             ylim = c(0.25, n_arm + 0.75),
                             expand = FALSE) +
    ggplot2::labs(x = x_title, y = NULL) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, title.position = "top")) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.text.y = ggplot2::element_text(face = "bold", colour = "black"),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank()
    )
  
  if (nrow(outside_dt) > 0) {
    p <- p +
      ggplot2::geom_segment(
        data = outside_dt,
        ggplot2::aes(x = xmid, xend = xmid, y = y_edge, yend = y_label - direction * 0.13),
        colour = "black", linewidth = 0.3
      ) +
      ggplot2::geom_text(
        data = outside_dt,
        ggplot2::aes(x = xmid, y = y_label, label = label),
        size = text_size, lineheight = 0.95, colour = "black"
      )
  }
  
  p
}

# Curly brace (accolade) as a path along x0..x1 next to a bar edge at y0. side = +1 draws it
# above the edge (tip pointing up), -1 below. rx (x units) and ry (y units) set the size of
# the rounded corners; narrow segments get a smaller rx.
brace_path <- function(x0, x1, y0, side, rx_max = 0.35, ry = 0.05, n = 12) {
  w  <- x1 - x0
  rx <- min(rx_max, w / 4)
  xm <- (x0 + x1) / 2
  th <- seq(0, pi / 2, length.out = n)
  lx <- c(x0 + rx - rx * cos(th), xm - rx + rx * sin(th))
  ly <- c(y0 + side * ry * sin(th), y0 + side * (2 * ry - ry * cos(th)))
  data.frame(x = c(lx, rev(2 * xm - lx)), y = c(ly, rev(ly)))
}

# Flexible version of make_time_bar_figure() with (1) numbers either inside the boxes or
# on curly braces next to the bars, (2) free placement of the bars, so extra (thinner)
# blocks can be added, and (3) per-bar choice of the side the labels go on.
# bar_dt: one row per bar x component, columns
#   y (numeric vertical position of the bar, larger = higher), y_label (axis label),
#   component, value, lower, upper (months), side (+1 labels above the bar, -1 below),
#   optional bar_height (default bar_height), optional connector_group (see
#   make_time_bar_figure; used by connectors = TRUE).
# legend_note: optional italic text centred just below the legend (at the top of the panel),
# e.g. to explain a lighter shade.
# legend_breaks / legend_labels: optional components (and their labels) to show in the legend;
# all components are still coloured. Use it to keep the legend small when several components
# share a colour (e.g. before / after parts of one component).
# component_colors: named vector, names = components in stacking order (left to right).
# label_style "boxes": numbers inside the segments (outside with a leader line when narrower
#   than small_threshold months); "braces": every segment gets a brace with the numbers.
# Neighbouring outside labels that would overlap are staggered over several levels.
# connectors: dashed lines between the segment boundaries of vertically neighbouring bars.
# extra_lines: optional data.frame (x, xend, y, yend) of additional dashed lines.
make_time_bar_figure2 <- function(bar_dt,
                                  component_colors,
                                  total_months = 24,
                                  label_style = c("boxes", "braces"),
                                  connectors = TRUE,
                                  extra_lines = NULL,
                                  small_threshold = 1.8,
                                  text_size = 3.3,
                                  bar_height = 0.4,
                                  min_label_dx = 3.0,
                                  level_step = 0.36,
                                  legend_title = "Mean time over 2 years",
                                  legend_nrow = 1,
                                  legend_breaks = NULL,
                                  legend_labels = legend_breaks,
                                  legend_note = NULL,
                                  x_title = "Months") {
  label_style <- match.arg(label_style)
  dt <- data.table::copy(data.table::as.data.table(bar_dt))
  dt[, component := factor(component, levels = names(component_colors))]
  if (!"bar_height" %in% names(dt)) dt[, bar_height := bar_height]
  if (!"connector_group" %in% names(dt)) dt[, connector_group := as.character(component)]
  data.table::setorder(dt, -y, component)
  
  # segment limits within each bar
  dt[, xmax := cumsum(value), by = y]
  dt[, xmin := xmax - value]
  dt[, xmid := (xmin + xmax) / 2]
  dt[, h := bar_height / 2]
  dt[, label := sprintf("%.1f\n(%.1f, %.1f)", value, lower, upper)]
  
  # dark text on light fills, white on dark fills
  lum <- function(col) {
    rgb <- grDevices::col2rgb(col) / 255
    0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
  }
  dt[, text_col := ifelse(lum(component_colors[as.character(component)]) < 0.55, "white", "black")]
  
  # which labels sit outside the bar
  dt[, outside := if (label_style == "braces") TRUE else value < small_threshold]
  
  # stagger overlapping outside labels (same bar, same side) over levels 0, 1, 2, ...
  dt[, level := 0L]
  for (key in unique(dt[outside == TRUE, paste(y, side)])) {
    idx <- which(dt$outside & paste(dt$y, dt$side) == key)
    idx <- idx[order(dt$xmid[idx])]
    last_x <- numeric(0)                       # x of the last label on each level
    for (i in idx) {
      lv <- 0L
      while (lv < length(last_x) && dt$xmid[i] - last_x[lv + 1] < min_label_dx) lv <- lv + 1L
      last_x[lv + 1] <- dt$xmid[i]
      data.table::set(dt, i, "level", lv)
    }
  }
  
  # geometry of the outside labels
  ry <- 0.05
  dt[, y_edge := y + side * (h + 0.04)]
  if (label_style == "braces") {
    dt[, y_tip := y_edge + side * 2 * ry]
  } else {
    dt[, y_tip := y + side * h]
  }
  dt[, lab_y := y_tip + side * ((if (label_style == "braces") 0.04 else 0.14) + level * level_step)]
  
  # braces
  brace_dt <- NULL
  if (label_style == "braces") {
    brace_dt <- data.table::rbindlist(lapply(seq_len(nrow(dt)), function(i) {
      b <- brace_path(dt$xmin[i], dt$xmax[i], dt$y_edge[i], dt$side[i], ry = ry)
      b$id <- i
      b
    }))
  }
  
  # leader lines from the bar / brace tip up to the label (when not directly adjacent)
  leader_dt <- dt[outside == TRUE & (label_style == "boxes" | level > 0)]
  
  # dashed connectors between boundaries of neighbouring bars
  ys <- sort(unique(dt$y), decreasing = TRUE)
  connect_dt <- NULL
  if (connectors && length(ys) > 1) {
    grp_order <- unique(dt[order(as.integer(component)), connector_group])
    bnd <- dt[, .(boundary = max(xmax), h = h[1]), by = .(y, connector_group)]
    bnd <- bnd[connector_group != grp_order[length(grp_order)]]
    bnd[, row := match(y, ys)]
    connect_dt <- merge(bnd, data.table::copy(bnd)[, row := row - 1],
                        by = c("row", "connector_group"), suffixes = c("_top", "_bottom"))
  }
  
  arms <- unique(dt[, .(y, y_label)])
  data.table::setorder(arms, y)
  
  # vertical room for the labels
  pad <- function(sd) {
    m <- dt[side == sd & outside == TRUE]
    if (nrow(m) == 0) return(0.35)
    0.35 + 0.25 + (if (label_style == "braces") 0.1 else 0) + max(m$level) * level_step
  }
  ylim <- c(min(dt$y - dt$h) - pad(-1) + 0.15,
            max(dt$y + dt$h) + pad(1) - 0.05 + (if (is.null(legend_note)) 0 else 0.4))
  
  p <- ggplot2::ggplot()
  if (!is.null(connect_dt) && nrow(connect_dt) > 0) {
    p <- p + ggplot2::geom_segment(
      data = connect_dt,
      ggplot2::aes(x = boundary_top, xend = boundary_bottom,
                   y = y_top - h_top, yend = y_bottom + h_bottom),
      linetype = "dashed", colour = "grey60", linewidth = 0.3)
  }
  if (!is.null(extra_lines)) {
    p <- p + ggplot2::geom_segment(
      data = extra_lines,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      linetype = "dashed", colour = "grey60", linewidth = 0.3)
  }
  p <- p +
    ggplot2::geom_rect(
      data = dt,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = y - h, ymax = y + h, fill = component),
      colour = "white", linewidth = 0.4)
  
  # numbers inside the segments (boxes style)
  inside_dt <- dt[outside == FALSE]
  if (nrow(inside_dt) > 0) {
    p <- p + ggplot2::geom_text(
      data = inside_dt,
      ggplot2::aes(x = xmid, y = y, label = label, colour = text_col),
      size = text_size, lineheight = 0.95) +
      ggplot2::scale_colour_identity()
  }
  
  # braces, leader lines and outside labels
  if (!is.null(brace_dt)) {
    p <- p + ggplot2::geom_path(data = brace_dt, ggplot2::aes(x = x, y = y, group = id),
                                colour = "black", linewidth = 0.35)
  }
  if (nrow(leader_dt) > 0) {
    p <- p + ggplot2::geom_segment(
      data = leader_dt,
      ggplot2::aes(x = xmid, xend = xmid, y = y_tip, yend = lab_y),
      colour = "black", linewidth = 0.3)
  }
  outside_dt <- dt[outside == TRUE]
  if (nrow(outside_dt) > 0) {
    p <- p + ggplot2::geom_text(
      data = outside_dt,
      ggplot2::aes(x = xmid, y = lab_y, label = label, vjust = ifelse(side > 0, 0, 1)),
      size = text_size, lineheight = 0.95, colour = "black")
  }
  
  if (!is.null(legend_note)) {
    p <- p + ggplot2::annotate("text", x = total_months / 2, y = ylim[2] - 0.05,
                               label = legend_note, fontface = "italic",
                               hjust = 0.5, vjust = 1, size = text_size)
  }
  
  p +
    ggplot2::scale_fill_manual(
      values = component_colors, name = legend_title,
      limits = names(component_colors)[names(component_colors) %in% as.character(dt$component)],
      breaks = if (is.null(legend_breaks)) ggplot2::waiver() else legend_breaks,
      labels = if (is.null(legend_breaks)) ggplot2::waiver() else legend_labels) +
    ggplot2::scale_x_continuous(breaks = seq(0, total_months, by = 3)) +
    ggplot2::scale_y_continuous(breaks = arms$y, labels = arms$y_label) +
    ggplot2::coord_cartesian(xlim = c(-0.1, total_months + 0.5), ylim = ylim, expand = FALSE) +
    ggplot2::labs(x = x_title, y = NULL) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = legend_nrow, title.position = "top")) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.text.y = ggplot2::element_text(face = "bold", colour = "black"),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank()
    )
}

# Builds the "months at home / in-center / death over 2 years" bar figure (Figure 2) from the
# bootstrap results and returns it as a ggplot object; saving is left to the caller. No
# objects from the global environment are used.
# Three bars (conservative management, choose dialysis, and an extra bar splitting the time at
# home of the dialysis arm into before / after the start of dialysis), linked by dashed lines;
# compact legend with a note on the lighter shade.
# Arguments
#   estimate_ci_dt   long table (trt, state, adjustment, estimate, ...) of the point estimates in
#                    DAYS, states as produced by compute_all_estimates() (incl. "Time at home",
#                    "In-center", "Days after dialysis", "In-center after dialysis")
#   bootstrap_dt     long table (b, trt, state, adjustment, value) of the bootstrap replicates, days
#   days_per_patient per-patient table with a trt column (used for the N per arm)
#   state_cols       names of the states; everything but "Death" counts as alive
#   state_colors     named colours "At home", "Hospitalization" (used for in-center) and "Death"
#   adjustment_level which adjustment level of the estimates to plot
#   days_per_month   days per month for the conversion to months (30.44 = 365.25 / 12)
#   total_months     length of the follow-up in months; death = total_months - months alive
create_time_bar_figure <- function(estimate_ci_dt,
                                   bootstrap_dt,
                                   days_per_patient,
                                   state_cols,
                                   state_colors,
                                   adjustment_level = "Adjusted for confounding and censoring",
                                   days_per_month = 30.44,
                                   total_months = 24) {
  # component names
  c_home_total  <- "Total time at home"
  c_home_before <- "Time at home before dialysis"
  c_home_after  <- "Time at home after dialysis"
  c_center      <- "Total time in-center"
  c_death       <- "Death"
  
  # colours: "after dialysis" is a lighter (opaque) version of the same colour (55% colour,
  # 45% white)
  lighten <- function(col) {
    grDevices::rgb(t(0.55 * grDevices::col2rgb(col) + 0.45 * 255), maxColorValue = 255)
  }
  green <- state_colors[["At home"]]
  blue  <- state_colors[["Hospitalization"]]
  # names = components in stacking order (left to right within a bar)
  colors_all <- c(
    setNames(green,          c_home_total),
    setNames(green,          c_home_before),
    setNames(lighten(green), c_home_after),
    setNames(blue,           c_center),
    setNames(state_colors[["Death"]], c_death)
  )
  bar_colors <- colors_all[c(c_home_before, c_home_after, c_home_total, c_center, c_death)]
  
  # estimates: one column per state (days) at the point estimate and per bootstrap replicate;
  # every component is a function of those columns, so its CI is taken within each replicate
  alive_states <- setdiff(state_cols, "Death")
  arms <- c("1", "0")
  pt_wide  <- data.table::dcast(
    estimate_ci_dt[adjustment == adjustment_level & as.character(trt) %in% arms],
    trt ~ state, value.var = "estimate")
  rep_wide <- data.table::dcast(
    bootstrap_dt[adjustment == adjustment_level & as.character(trt) %in% arms],
    b + trt ~ state, value.var = "value")
  pt_wide[, trt := as.character(trt)]
  rep_wide[, trt := as.character(trt)]
  
  # months for one component: f() takes a table with the state columns (days), returns days
  component_row <- function(trt_value, component, f) {
    p <- pt_wide[trt == trt_value]
    r <- rep_wide[trt == trt_value]
    data.table::data.table(
      trt = trt_value, component = component,
      value = f(p) / days_per_month,
      lower = quantile(f(r) / days_per_month, 0.025, na.rm = TRUE, names = FALSE),
      upper = quantile(f(r) / days_per_month, 0.975, na.rm = TRUE, names = FALSE))
  }
  f_home_total  <- function(d) d[["Time at home"]]
  f_center      <- function(d) d[["In-center"]]
  f_alive       <- function(d) rowSums(as.matrix(d[, ..alive_states]))
  f_death       <- function(d) total_months * days_per_month - f_alive(d)
  # after the start of dialysis: home = days after start - in-center days after start (includes
  # days at home while on HD/PD outside the in-center sessions/visits); before = total - after,
  # so before + after add up exactly to the total
  f_home_after  <- function(d) d[["Days after dialysis"]] - d[["In-center after dialysis"]]
  f_home_before <- function(d) f_home_total(d) - f_home_after(d)
  
  # the three components shared by all versions, for both arms
  base_dt <- data.table::rbindlist(lapply(arms, function(tv) rbind(
    component_row(tv, c_home_total, f_home_total),
    component_row(tv, c_center,     f_center),
    component_row(tv, c_death,      f_death)
  )))
  # dialysis arm: home time split before / after the start of dialysis
  home_split_dt <- rbind(component_row("1", c_home_before, f_home_before),
                         component_row("1", c_home_after,  f_home_after))
  
  # axis labels with the sample size
  n_dt <- data.table::as.data.table(days_per_patient)[, .(n = .N), by = .(trt = as.character(trt))]
  arm_label <- function(trt_value) {
    paste0(ifelse(trt_value == "1", "Choose dialysis", "Choose conservative\nmanagement"),
           "\n(N=", n_dt[trt == trt_value, n], ")")
  }
  
  # ---- layout: conservative management / choose dialysis / time at home split ---------
  # vertical positions (larger = higher): the two dialysis bars are close together, with more
  # room between them and the conservative management bar
  y_cm    <- 3.2
  y_dial  <- 2
  y_split <- 1.4
  bar_dt <- rbind(
    base_dt[, .(trt, component, value, lower, upper,
                y = ifelse(trt == "1", y_dial, y_cm),
                y_label = sapply(trt, arm_label))],
    home_split_dt[, .(trt, component, value, lower, upper,
                      y = y_split,
                      y_label = "Time at home:\nbefore vs. after\nstart of dialysis")]
  )
  bar_dt[, side := ifelse(y == y_cm, 1, -1)]
  # connector groups: the parts of the time at home count as one segment
  bar_dt[, connector_group := data.table::fcase(
    grepl("home",      component), "home",
    grepl("in-center", component), "center",
    default = "death")]
  make_time_bar_figure2(
    bar_dt           = bar_dt,
    component_colors = bar_colors,
    total_months     = total_months,
    label_style      = "boxes",
    connectors       = TRUE,
    # compact legend: only the three components; the lighter shade is explained in an italic
    # note below the legend
    legend_breaks    = c(c_home_total, c_center, c_death),
    legend_labels    = c("Time at home", "Time in-center", "Death"),
    legend_note      = "Lighter shade = time at home after start of dialysis",
    legend_nrow      = 1
  )
}