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
                           trt_labels = c("trt=0", "trt=1"),
                           show_censoring = TRUE,  # show censoring ticks?
                           censor_height = 0.03,   # tick height (in survival probability units)
                           censor_linewidth = 0.4, # tick line width in mm (increase if too faint)
                           censor_alpha = 0.9,     # tick transparency (1 = opaque)
                           label_width = 13,       # max. characters per line in risk-table labels
                           font_family = "sans") {
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
  
  # ---- shared style: same font, size and color for plot and table ----
  base_size <- 14
  txt_size  <- base_size * 0.8      # = size of the axis labels (theme_minimal: rel(0.8))
  txt_col   <- "black"              # one text color for axes, legend, title and table
  x_breaks  <- seq(0, horizon, by = 365 / 2)  # every 6 months
  x_expand  <- ggplot2::expansion(mult = 0.03)
  
  # legend inside the plot (lower left is empty); use 'inside' if ggplot2 >= 3.5
  legend_theme <- if (utils::packageVersion("ggplot2") >= "3.5.0") {
    ggplot2::theme(legend.position = "inside",
                   legend.position.inside = c(0.02, 0.03))
  } else {
    ggplot2::theme(legend.position = c(0.02, 0.03))
  }
  
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
    )
  
  # censoring ticks: thin and light (turn off with show_censoring = FALSE)
  if (show_censoring) {
    KM_plot <- KM_plot +
      # vertical segments instead of point symbols: the width is set directly
      # by linewidth, so it does not depend on the symbol/font rendering
      ggplot2::geom_segment(data = censor_data,
                            ggplot2::aes(x = time, xend = time,
                                         # clamp to 0-1 so no tick falls outside the y-scale limits
                                         y = pmax(surv - censor_height / 2, 0),
                                         yend = pmin(surv + censor_height / 2, 1)),
                            linewidth = censor_linewidth,
                            alpha = censor_alpha,
                            show.legend = FALSE)
  }
  
  KM_plot <- KM_plot +
    # y title at the top left (horizontal), same left edge as the risk table
    ggplot2::labs(
      x = paste0("Time (", unit, ")"),
      y = NULL,
      title = "Survival probability (%)"
    ) +
    ggplot2::theme_minimal(base_size = base_size, base_family = font_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = txt_size, color = txt_col,
                                         hjust = 0,
                                         margin = ggplot2::margin(b = 8)),
      plot.title.position = "plot",
      axis.text = ggplot2::element_text(color = txt_col),
      axis.title = ggplot2::element_text(color = txt_col),
      legend.text = ggplot2::element_text(color = txt_col, size = txt_size),
      legend.title = ggplot2::element_blank(),
      legend.justification = c(0, 0),
      legend.background = ggplot2::element_blank(),
      legend.key.height = grid::unit(0.9, "lines"),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.25),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.25),
      axis.ticks.length = grid::unit(2, "pt"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    legend_theme +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.1),
      labels = seq(0, 100, by = 10),
      expand = c(0, 0)
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = function(x) round(x / 30.44), # days -> months
      expand = x_expand
    ) +
    ggplot2::coord_cartesian(xlim = c(0, horizon)) +
    # reverse = TRUE: same order as the risk table (Dialysis on top, CM below)
    ggplot2::scale_color_manual(labels = trt_labels,
                                values = manual_colors,
                                guide = ggplot2::guide_legend(reverse = TRUE)) +
    ggplot2::scale_fill_manual(guide = "none",
                               values = manual_colors) +
    ggplot2::scale_linetype_manual(guide = "none", 
                                   values = c("solid", "solid"))
  
  # find numbers at risk at each time point of interest
  KM_table <- summary(KM_fit, times = x_breaks)
  risk_table <- data.table(
    time = KM_table$time,            # in days: same x-axis as the KM plot
    strata = factor(as.numeric(KM_table$strata)-1, 
                    levels = c(0, 1)),
    n.risk = round(KM_table$n.risk)  # after weighting might not be integer
  )
  
  # row labels: drop underscores, wrap long labels over two lines, color via HTML
  wrap_label <- function(x) {
    x <- gsub("_", " ", x)
    paste(strwrap(x, width = label_width), collapse = "<br>")
  }
  colored_labels <- paste0("<span style='color:", manual_colors[1:2], "'>",
                           vapply(trt_labels, wrap_label, character(1)),
                           "</span>")
  
  # create table
  KM_table <- ggplot2::ggplot(risk_table, 
                              ggplot2::aes(x = time, 
                                           y = strata, 
                                           label = n.risk)) +
    ggplot2::geom_text(size = txt_size / ggplot2::.pt,
                       color = txt_col,
                       family = font_family) +
    ggplot2::scale_x_continuous(breaks = x_breaks, expand = x_expand) +
    # extra space above/below so the rows are further apart
    ggplot2::scale_y_discrete(labels = colored_labels,
                              expand = ggplot2::expansion(add = 0.8)) + 
    ggplot2::coord_cartesian(xlim = c(0, horizon), clip = "off") +
    ggplot2::labs(title = "Number at risk") +
    ggplot2::theme_void(base_size = base_size, base_family = font_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = txt_size, color = txt_col,
                                         face = "bold", hjust = 0),
      plot.title.position = "plot",
      # left-aligned labels, with padding between label and numbers
      axis.text.y = ggtext::element_markdown(size = txt_size, color = txt_col,
                                             hjust = 0, halign = 0, lineheight = 1,
                                             margin = ggplot2::margin(r = 10)),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
  
  return(list(KM_plot = KM_plot,
              KM_table = KM_table))
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

# Horizontal stacked bars of the mean time (months) in each component over the follow-up
# (e.g. time at home, time in-center, death = the remaining time), so every bar adds up to
# total_months. Free placement of the bars (so extra bars can be added), numbers inside the
# boxes (bold estimate, CI below) or, for narrow segments, in one line outside the bar with a
# short leader line, optional curly braces, optional text under segments, and per-bar choice
# of the side the outside labels go on. Needs ggtext.
# bar_dt: one row per bar x component, columns
#   y (numeric vertical position of the bar, larger = higher), y_label (axis label, markdown),
#   component, value, lower, upper (months), side (+1 outside labels above the bar, -1 below),
#   optional bar_height, optional connector_group, optional connect (FALSE = this bar takes
#   no part in the dashed connectors).
# component_colors: named vector, names = components in stacking order (left to right).
# legend_title: NULL (default) = no legend title. legend_size: size of the legend relative to
#   the base font.
# braces_extra: optional data.frame (x0, x1, y, side): one brace spanning x0..x1 along the bar
#   edge at y (side = -1: below the edge, tip pointing down; +1: above).
# segment_labels: optional data.frame (x, y, label): plain text, top-aligned at y, centred at x.
# caption: optional small note below the figure.
# label_gap: free space (y units) between a leader line and its label.
# label_height: height (y units) reserved for a two-line label (one-line labels get 55%).
make_time_bar_figure <- function(bar_dt,
                                 component_colors,
                                 total_months = 24,
                                 braces_extra = NULL,
                                 segment_labels = NULL,
                                 small_threshold = 1.8,
                                 text_size = 3.3,
                                 base_size = 11,
                                 bar_height = 0.4,
                                 min_label_dx = 4,
                                 level_step = NULL,
                                 label_gap = 0.08,
                                 label_height = 0.09 * text_size,
                                 x_break_by = 6,
                                 axis_colour = "grey30",
                                 legend_title = NULL,
                                 legend_size = 1.15,
                                 legend_nrow = 1,
                                 legend_breaks = NULL,
                                 legend_labels = legend_breaks,
                                 caption = NULL,
                                 x_title = "Months") {
  if (is.null(level_step)) level_step <- label_height + 0.03
  dt <- data.table::copy(data.table::as.data.table(bar_dt))
  dt[, component := factor(component, levels = names(component_colors))]
  if (!"bar_height" %in% names(dt)) dt[, bar_height := bar_height]
  if (!"connector_group" %in% names(dt)) dt[, connector_group := as.character(component)]
  if (!"connect" %in% names(dt)) dt[, connect := TRUE]
  data.table::setorder(dt, -y, component)
  
  # segment limits within each bar
  dt[, xmax := cumsum(value), by = y]
  dt[, xmin := xmax - value]
  dt[, xmid := (xmin + xmax) / 2]
  dt[, h := bar_height / 2]
  
  # dark text on light fills, white on dark fills
  lum <- function(col) {
    rgb <- grDevices::col2rgb(col) / 255
    0.2126 * rgb[1, ] + 0.7152 * rgb[2, ] + 0.0722 * rgb[3, ]
  }
  dt[, text_col := ifelse(lum(component_colors[as.character(component)]) < 0.55, "white", "black")]
  
  # narrow segments get their label outside the bar
  dt[, outside := value < small_threshold]
  
  # labels: bold estimate + CI on the next line; narrow segments get one line outside the bar
  dt[, label := ifelse(outside,
                       sprintf("%.1f (%.1f, %.1f)", value, lower, upper),
                       sprintf("<b>%.1f</b><br>(%.1f, %.1f)", value, lower, upper))]
  dt[, lab_h := ifelse(outside, 0.55 * label_height, label_height)]
  
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
  dt[, y_tip := y + side * h]
  dt[, lab_y := y_tip + side * (0.2 + level * level_step)]
  # leader lines stop label_gap short of the label so they never touch the text
  dt[, leader_end := lab_y - side * label_gap]
  
  # extra braces (e.g. one brace under a whole segment)
  extra_brace_dt <- NULL
  if (!is.null(braces_extra) && nrow(braces_extra) > 0) {
    extra_brace_dt <- data.table::rbindlist(lapply(seq_len(nrow(braces_extra)), function(i) {
      b <- brace_path(braces_extra$x0[i], braces_extra$x1[i], braces_extra$y[i],
                      braces_extra$side[i], rx_max = 0.6, ry = 0.09)
      b$id <- i
      b
    }))
  }
  
  # leader lines from the bar / brace tip to the label
  leader_dt <- dt[outside == TRUE]
  
  # dashed connectors between boundaries of neighbouring (connecting) bars
  connect_dt <- NULL
  cdt <- dt[connect == TRUE]
  ys <- sort(unique(cdt$y), decreasing = TRUE)
  if (length(ys) > 1) {
    grp_order <- unique(dt[order(as.integer(component)), connector_group])
    bnd <- cdt[, .(boundary = max(xmax), h = h[1]), by = .(y, connector_group)]
    bnd <- bnd[connector_group != grp_order[length(grp_order)]]
    bnd[, row := match(y, ys)]
    connect_dt <- merge(bnd, data.table::copy(bnd)[, row := row - 1],
                        by = c("row", "connector_group"), suffixes = c("_top", "_bottom"))
  }
  
  arms <- unique(dt[, .(y, y_label)])
  data.table::setorder(arms, y)
  
  # panel limits: exactly as much room as the bars and their outside labels need
  out <- dt[outside == TRUE]
  lo <- min(c(dt$y - dt$h, out[side == -1, lab_y - lab_h]))
  hi <- max(c(dt$y + dt$h, out[side == 1,  lab_y + lab_h]))
  if (!is.null(segment_labels)) lo <- min(lo, segment_labels$y - 0.55 * label_height)
  ylim <- c(lo - 0.08, hi + 0.06)
  
  p <- ggplot2::ggplot()
  if (!is.null(connect_dt) && nrow(connect_dt) > 0) {
    p <- p + ggplot2::geom_segment(
      data = connect_dt,
      ggplot2::aes(x = boundary_top, xend = boundary_bottom,
                   y = y_top - h_top, yend = y_bottom + h_bottom),
      linetype = "dashed", colour = "grey70", linewidth = 0.3)
  }
  p <- p +
    ggplot2::geom_rect(
      data = dt,
      ggplot2::aes(xmin = xmin, xmax = xmax, ymin = y - h, ymax = y + h, fill = component),
      colour = "white", linewidth = 0.4)
  
  # numbers inside the segments
  no_pad <- grid::unit(rep(0, 4), "pt")
  inside_dt <- dt[outside == FALSE]
  if (nrow(inside_dt) > 0) {
    p <- p + ggtext::geom_richtext(
      data = inside_dt,
      ggplot2::aes(x = xmid, y = y, label = label, colour = text_col),
      size = text_size, lineheight = 0.95, fill = NA, label.colour = NA,
      label.padding = no_pad, label.margin = no_pad) +
      ggplot2::scale_colour_identity()
  }
  
  # braces, leader lines and outside labels
  if (!is.null(extra_brace_dt)) {
    p <- p + ggplot2::geom_path(data = extra_brace_dt, ggplot2::aes(x = x, y = y, group = id),
                                colour = "grey40", linewidth = 0.4)
  }
  if (nrow(leader_dt) > 0) {
    p <- p + ggplot2::geom_segment(
      data = leader_dt,
      ggplot2::aes(x = xmid, xend = xmid, y = y_tip, yend = leader_end),
      colour = "grey30", linewidth = 0.3)
  }
  if (nrow(out) > 0) {
    p <- p + ggtext::geom_richtext(
      data = out,
      ggplot2::aes(x = xmid, y = lab_y, label = label, vjust = ifelse(side > 0, 0, 1)),
      size = text_size, lineheight = 0.95, colour = "black", fill = NA, label.colour = NA,
      label.padding = no_pad, label.margin = no_pad)
  }
  if (!is.null(segment_labels)) {
    p <- p + ggplot2::geom_text(
      data = segment_labels, ggplot2::aes(x = x, y = y, label = label),
      vjust = 1, size = text_size, colour = "grey20")
  }
  
  p +
    ggplot2::scale_fill_manual(
      values = component_colors, name = legend_title,
      limits = names(component_colors)[names(component_colors) %in% as.character(dt$component)],
      breaks = if (is.null(legend_breaks)) ggplot2::waiver() else legend_breaks,
      labels = if (is.null(legend_breaks)) ggplot2::waiver() else legend_labels) +
    ggplot2::scale_x_continuous(breaks = seq(0, total_months, by = x_break_by)) +
    ggplot2::scale_y_continuous(breaks = arms$y, labels = arms$y_label) +
    ggplot2::coord_cartesian(xlim = c(-0.1, total_months + 0.1), ylim = ylim, expand = FALSE) +
    ggplot2::labs(x = x_title, y = NULL, caption = caption) +
    ggplot2::guides(fill = if (is.null(legend_title)) {
      ggplot2::guide_legend(nrow = legend_nrow)
    } else {
      ggplot2::guide_legend(nrow = legend_nrow, title.position = "top")
    }) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = if (is.null(legend_title)) ggplot2::element_blank() else
        ggplot2::element_text(face = "bold", hjust = 0.5),
      legend.text = ggplot2::element_text(size = ggplot2::rel(legend_size), colour = "grey20"),
      legend.key.size = grid::unit(1.2 * legend_size, "lines"),
      legend.spacing.x = grid::unit(8, "pt"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(4, "pt"),
      axis.text.y = ggtext::element_markdown(colour = "black", hjust = 1, lineheight = 1.1,
                                             margin = ggplot2::margin(r = 6)),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(colour = axis_colour),
      axis.title.x = ggplot2::element_text(colour = axis_colour),
      axis.line.x = ggplot2::element_line(colour = "grey50"),
      axis.ticks.x = ggplot2::element_line(colour = "grey50"),
      plot.caption = ggplot2::element_text(hjust = 0, colour = "grey30",
                                           size = ggplot2::rel(0.75)),
      plot.caption.position = "plot",
      plot.margin = ggplot2::margin(4, 10, 4, 4)
    )
}

# Builds the "months alive at home / alive in-center / dead over 2 years" bar figure (Figure 2)
# from the bootstrap results and returns it as a ggplot object; saving is left to the caller.
# No objects from the global environment are used. Needs ggtext.
# Layout: conservative management on top, dialysis below it (dashed lines connect the segment
# boundaries), and below that a bar splitting the time at home of the dialysis arm into before /
# after the start of dialysis, linked to the home segment of the dialysis bar by one curly brace.
# Arguments
#   estimate_ci_dt   long table (trt, state, adjustment, estimate, ...) of the point estimates in
#                    DAYS, states as produced by compute_all_estimates() (incl. "Time at home",
#                    "In-center", "Days after dialysis", "In-center after dialysis")
#   bootstrap_dt     long table (b, trt, state, adjustment, value) of the bootstrap replicates, days
#   days_per_patient per-patient table with a trt column (used for the n per arm)
#   state_cols       names of the states; everything but "Death" counts as alive
#   state_colors     named colours "At home", "Hospitalization" (used for in-center) and "Death"
#   adjustment_level which adjustment level of the estimates to plot
#   days_per_month   days per month for the conversion to months (30.44 = 365.25 / 12)
#   total_months     length of the follow-up in months; dead = total_months - months alive
#   arm_names        axis labels of the arms (markdown; <br> = line break), named "1" and "0"
#   legend_labels    legend labels: home, in-center, dead
#   text_size        size of the numbers in the figure (ggplot size, mm)
#   bar_height       height of the bars (y units)
#   bar_gap_main     space between the conservative management and dialysis bars
#   bar_gap_brace    space between the dialysis bar and the split bar (holds the brace)
#   legend_size      size of the legend relative to the base font
#   ci_caption       add a caption in the figure stating how the CIs were obtained (default
#                    FALSE: put this in the figure legend of the manuscript instead)
#   check_ci         print a check of where each point estimate lies within its bootstrap CI
#                    (also stored as attr(<plot>, "ci_check"))
create_time_bar_figure <- function(estimate_ci_dt,
                                   bootstrap_dt,
                                   days_per_patient,
                                   state_cols,
                                   state_colors,
                                   adjustment_level = "Adjusted for confounding and censoring",
                                   days_per_month = 30.44,
                                   total_months = 24,
                                   arm_names = c(`1` = "Choose dialysis",
                                                 `0` = "Choose conservative<br>management"),
                                   legend_labels = c("Alive at home", "Alive in-center", "Dead"),
                                   text_size = 4.2,
                                   bar_height = 0.62,
                                   bar_gap_main = 0.5,
                                   bar_gap_brace = 0.4,
                                   legend_size = 1.2,
                                   ci_caption = FALSE,
                                   check_ci = TRUE) {
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
  
  # months for one component: f() takes a table with the state columns (days), returns days.
  # Also returns the mean / median of the bootstrap replicates for the CI check
  component_row <- function(trt_value, component, f) {
    p <- pt_wide[trt == trt_value]
    r <- rep_wide[trt == trt_value]
    boot <- f(r) / days_per_month
    data.table::data.table(
      trt = trt_value, component = component,
      value = f(p) / days_per_month,
      lower = quantile(boot, 0.025, na.rm = TRUE, names = FALSE),
      upper = quantile(boot, 0.975, na.rm = TRUE, names = FALSE),
      boot_mean = mean(boot, na.rm = TRUE),
      boot_median = median(boot, na.rm = TRUE),
      n_boot = sum(!is.na(boot)))
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
  
  # ---- check: where does each point estimate sit within its percentile bootstrap CI? -----
  # position 0 = at the lower limit, 1 = at the upper limit, ~0.5 = centred. A point estimate
  # near one end means the bootstrap distribution is shifted relative to the estimate
  # (bootstrap bias), e.g. because the replicates use fewer imputations / resampled data
  ci_check <- rbind(base_dt, home_split_dt)[, .(
    arm = ifelse(trt == "1", "dialysis", "conservative"), component,
    estimate = value, lower, upper, boot_mean, boot_median,
    bias = boot_mean - value,
    position = (value - lower) / (upper - lower),
    n_boot)]
  ci_check[, flag := position < 0.15 | position > 0.85]
  if (check_ci) {
    cat("\nCheck of point estimates versus percentile bootstrap CIs (months):\n")
    print(ci_check[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)])
    if (any(ci_check$flag)) {
      message("Point estimate lies in the outer 15% of its bootstrap CI for: ",
              paste(ci_check[flag == TRUE, paste0(component, " (", arm, ")")], collapse = "; "),
              ". Check that estimate and replicates come from the same procedure (same number of ",
              "imputations, same estimator) before reporting.")
    }
  }
  
  # axis labels: bold name, sample size in grey below it
  n_dt <- data.table::as.data.table(days_per_patient)[, .(n = .N), by = .(trt = as.character(trt))]
  grey_txt <- function(x) paste0("<span style='color:#6B6B6B'>", x, "</span>")
  arm_label <- function(trt_value) {
    paste0("<b>", arm_names[[trt_value]], "</b><br>",
           grey_txt(paste0("N = ", format(n_dt[trt == trt_value, n], big.mark = ",", trim = TRUE))))
  }
  
  # ---- layout: conservative management / dialysis / time at home split ----------------
  y_split <- 1
  y_dial  <- y_split + bar_height + bar_gap_brace
  y_cm    <- y_dial  + bar_height + bar_gap_main
  h <- bar_height / 2
  bar_dt <- rbind(
    base_dt[, .(trt, component, value, lower, upper,
                y = ifelse(trt == "1", y_dial, y_cm),
                y_label = sapply(trt, arm_label),
                connect = TRUE)],
    home_split_dt[, .(trt, component, value, lower, upper,
                      y = y_split,
                      y_label = paste0("<b>Time at home</b><br>", grey_txt("choose dialysis group")),
                      connect = FALSE)]
  )
  bar_dt[, side := ifelse(y == y_cm, 1, -1)]
  # connector groups: the parts of the time at home count as one segment
  bar_dt[, connector_group := data.table::fcase(
    grepl("home",      component), "home",
    grepl("in-center", component), "center",
    default = "death")]
  
  # one brace under the home segment of the dialysis bar (pointing at the split bar), and the
  # names of the two parts under the split bar
  home_total_dial <- base_dt[trt == "1" & component == c_home_total, value]
  v_before <- home_split_dt[component == c_home_before, value]
  v_after  <- home_split_dt[component == c_home_after,  value]
  braces_extra <- data.frame(x0 = 0, x1 = home_total_dial, y = y_dial - h - 0.05, side = -1)
  segment_labels <- data.frame(
    x = c(v_before / 2, v_before + v_after / 2),
    y = y_split - h - 0.1,
    label = c("Before dialysis start", "After dialysis start"))
  
  n_boot <- max(ci_check$n_boot)
  caption <- if (ci_caption) {
    sprintf(paste0("Point estimates from the full cohort; 95%% confidence intervals are ",
                   "percentile intervals from %s bootstrap replicates."),
            format(n_boot, big.mark = ","))
  } else NULL
  
  p <- make_time_bar_figure(
    bar_dt           = bar_dt,
    component_colors = bar_colors,
    total_months     = total_months,
    braces_extra     = braces_extra,
    segment_labels   = segment_labels,
    bar_height       = bar_height,
    text_size        = text_size,
    base_size        = 13,
    x_break_by       = 6,
    legend_title     = NULL,
    legend_size      = legend_size,
    legend_breaks    = c(c_home_total, c_center, c_death),
    legend_labels    = legend_labels,
    legend_nrow      = 1,
    caption          = caption
  )
  attr(p, "ci_check") <- ci_check
  p
}