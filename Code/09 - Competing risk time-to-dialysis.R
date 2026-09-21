################################################################################
### Decision for dialysis versus conservative management
### PART - Competing risk analysis for time-to-dialysis
################################################################################

# set-up
rm(list = ls(all.names = TRUE))
knitr::opts_knit$set(root.dir = "P:/SCREAM2/SCREAM2_Research/Carolien Maas/")
set.seed(1)
setwd(
  "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/"
)
results_path <- "P:/SCREAM2/SCREAM2_Research/Carolien Maas/Project Dialysis versus Conservative Management/Results/"

# load libraries
library(mstate)
library(patchwork)
library(foreach)   # parallel computation
library(doRNG)     # handle parallel seeds
set.seed(1)        # set seed for parallel backend

# load functions
source("Code/utils/data_manipulation.R")
source("Code/utils/weighting.R")
source("Code/utils/competing_risk.R")
source("Code/utils/plots.R")
source("Code/utils/compute_absolute_relative_risks.R")

# load data
load("Data/merged_ckd.Rdata")
load("Data/cohort_with_prob.Rdata")

# Two patients have their KRT exactly at the time horizon, we can ignore warnings
baseline[LOPNR == 287508019, .(
  LOPNR,
  event_KRT_2y,
  time2event_KRT_2y,
  event_KRT_inf,
  time2event_KRT_inf,
  event_death_2y,
  time2event_death_2y,
  event_death_inf,
  time2event_death_inf
)]
baseline[LOPNR == 366058177, .(
  LOPNR,
  event_KRT_2y,
  time2event_KRT_2y,
  event_KRT_inf,
  time2event_KRT_inf,
  event_death_2y,
  time2event_death_2y,
  event_death_inf,
  time2event_death_inf
)]

# prepare data with patients who chose dialysis in long format
vars <- c(
  id_name,
  "event_KRT_inf",
  "time2event_KRT_inf",
  "event_death_inf",
  "time2event_death_inf",
  "sw_IPTW"
)
cmp_dt <- baseline[trt == 1, ..vars, with = FALSE]

# estimate state probabilities
cmp_result <- state_probabilities(dt = cmp_dt, horizon = horizon)

# create data frame for CI plot
state_prob_df <- data.frame(
  time = rep(cmp_result$times, 3) / 365,
  state_prob = c(
    cmp_result$state_prob$P11,
    cmp_result$state_prob$P12,
    cmp_result$state_prob$P13
  ),
  Event = rep(
    c("Event_free", "Dialysis", "Death"),
    each = nrow(cmp_result$state_prob)
  )
)

### Use bootstrapping for 95% CI
# 1. Make a bootstrap sample of the full baseline dt (dia and CCC)
# 2. Fit a IPTW model on each bootstrap sample
# 3. Fit Cox model using bootstrapped model and bootstrapped IPTW
state_prob_boot <- data.frame(time = state_prob_df$time)
RMST_boot <- data.frame(RMST = c("Event_free", "Dialysis", "Death"))
for (B in 1:n_bootstraps) {
  # create bootstrap sample
  bootstrap <- baseline[sample(1:nrow(baseline), replace = TRUE), ]
  
  # compute IPTW on bootstrap
  out_weights <- create_weights(
    data = bootstrap,
    id_name = id_name,
    model_PS = model_PS,
    w_meth = "IPTW",
    verbose = FALSE
  )
  bootstrap[["sw_IPTW"]] <- as.numeric(out_weights$data$w)
  
  # prepare data with patients who chose dialysis in long format
  cmp_boot_dt <- bootstrap[trt == 1, ..vars, with = FALSE]
  
  # estimate state probabilities
  cmp_boot_result <- state_probabilities(cmp_boot_dt, horizon)
  
  # save state probabilities for each bootstrap
  state_prob_boot[, paste0("boot_", B)] <-
    c(
      cmp_boot_result$state_prob$P11,
      cmp_boot_result$state_prob$P12,
      cmp_boot_result$state_prob$P13
    )
  RMST_boot[, paste0("boot_", B)] <-
    unlist(cmp_boot_result$RMST)
}
state_prob_df$state_prob_lower <- apply(state_prob_boot[, -1], 1, quantile, probs = 0.025)
state_prob_df$state_prob_upper <- apply(state_prob_boot[, -1], 1, quantile, probs = 0.975)
RMST_lower <- apply(RMST_boot[, -1], 1, quantile, probs = 0.025)
RMST_upper <- apply(RMST_boot[, -1], 1, quantile, probs = 0.975)

################################################################################
### SANITY CHECK
################################################################################
# RMST_death has to be the same as the KM RMST
fit <- survival::survfit(survival::Surv(time2event_death_2y, event_death_2y) ~ 1,
                         data = baseline[trt == 1, ])
RMST <- summary(fit, rmean = horizon)$table["rmean"]
# P13 RMST "decision -> death" should be 19.7
RMST / 30.5
# 19.71609
cmp_result$RMST$RMST_13
# 19.27663

################################################################################
### Create plot
################################################################################
# --- Left: State probability curves ---
state_colors <- c(manual_colors[3], manual_colors[2], manual_colors[1])
p1 <- ggplot2::ggplot(state_prob_df,
                      ggplot2::aes(
                        x = time,
                        y = state_prob,
                        fill = Event,
                        color = Event
                      )) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = state_prob_lower,
      ymax = state_prob_upper,
      group = Event,
      fill = Event
    ),
    alpha = 0.2,
    colour = NA
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Event_free" = state_colors[1],
      "Dialysis" = state_colors[2],
      "Death" = state_colors[3]
    )
  ) +
  ggplot2::scale_color_manual(
    labels = c("Alive dialysis-free", "Dialysis", "Death"),
    values = c(
      "Event_free" = state_colors[1],
      "Dialysis" = state_colors[2],
      "Death" = state_colors[3]
    ),
    breaks = c("Event_free", "Dialysis", "Death")
  ) +
  ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(fill = NA)),
                  fill = "none") +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                              breaks = seq(0, 1, by = 0.1)) +
  ggplot2::theme_minimal() +
  ggplot2::labs(x = "Time (years)", y = "Probability (%)", color = "Transition") +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = ggplot2::element_blank()
  ) +
  ggplot2::annotate(
    "text",
    x = 0.5,
    y = 1,
    label = paste0("At 2 years\nAlive dialysis-free\nDialysis\nDeath"),
    hjust = 0,
    vjust = 1,
    size = 3
  ) +
  ggplot2::annotate(
    "text",
    x = 1.1,
    y = 1,
    label = paste0(
      "%, 95% CI\n",
      fmt(
        as.numeric(
          state_prob_df |>
            dplyr::filter(Event == "Event_free" &
                            time == horizon / 365) |>
            dplyr::select(state_prob)
        ) * 100,
        1
      ),
      " (",
      fmt(
        as.numeric(
          state_prob_df |>
            dplyr::filter(Event == "Event_free" &
                            time == horizon / 365) |>
            dplyr::select(state_prob_lower)
        ) * 100,
        1
      ),
      ", ",
      fmt(
        as.numeric(
          state_prob_df |>
            dplyr::filter(Event == "Event_free" &
                            time == horizon / 365) |>
            dplyr::select(state_prob_upper)
        ) * 100,
        1
      ),
      ")\n",
      fmt(
        as.numeric(
          state_prob_df |>
            dplyr::filter(Event == "Dialysis" &
                            time == horizon / 365) |>
            dplyr::select(state_prob)
        ) * 100,
        1
      ),
      " (",
      fmt(
        as.numeric(
          state_prob_df |>
            dplyr::filter(Event == "Dialysis" &
                            time == horizon / 365) |>
            dplyr::select(state_prob_lower)
        ) * 100,
        1
      ),
      ", ",
      fmt(
        as.numeric(
          state_prob_df |>
            dplyr::filter(Event == "Dialysis" &
                            time == horizon / 365) |>
            dplyr::select(state_prob_upper)
        ) * 100,
        1
      ),
      ")\n",
      fmt(
        as.numeric(
          state_prob_df |>
            dplyr::filter(Event == "Death" &
                            time == horizon / 365) |>
            dplyr::select(state_prob)
        ) * 100,
        1
      ),
      " (",
      fmt(
        as.numeric(
          state_prob_df |>
            dplyr::filter(Event == "Death" &
                            time == horizon / 365) |>
            dplyr::select(state_prob_lower)
        ) * 100,
        1
      ),
      ", ",
      fmt(
        as.numeric(
          state_prob_df |>
            dplyr::filter(Event == "Death" &
                            time == horizon / 365) |>
            dplyr::select(state_prob_upper)
        ) * 100,
        1
      ),
      ")"
    ),
    hjust = 0,
    vjust = 1,
    size = 3
  ) +
  ggplot2::annotate(
    "text",
    x = 1.6,
    y = 1,
    label = paste0(
      "months, 95% CI\n",
      fmt(cmp_result$RMST$RMST_11, 1),
      " (",
      fmt(RMST_lower[1], 1),
      ", ",
      fmt(RMST_upper[1], 1),
      ")\n",
      fmt(cmp_result$RMST$RMST_12, 1),
      " (",
      fmt(RMST_lower[2], 1),
      ", ",
      fmt(RMST_upper[2], 1),
      ")\n",
      fmt(cmp_result$RMST$RMST_13, 1),
      " (",
      fmt(RMST_lower[3], 1),
      ", ",
      fmt(RMST_upper[3], 1),
      ")"
    ),
    hjust = 0,
    vjust = 1,
    size = 3
  )

# --- Right: Stacked state probabilities ---
# Long format for stacked plot
df_long <- state_prob_df[, 1:3] |>
  tidyr::pivot_longer(cols = state_prob,
                      names_to = "State",
                      values_to = "Probability") |>
  dplyr::mutate(Event = factor(Event, levels = c("Dialysis", "Event_free", "Death")))
p2 <- ggplot2::ggplot(df_long, ggplot2::aes(x = time, y = Probability, fill = Event)) +
  ggplot2::geom_area(color = "black",
                     linewidth = 0.2,
                     alpha = 0.8) +
  ggplot2::scale_fill_manual(
    values = c(
      "Event_free" = state_colors[1],
      "Dialysis" = state_colors[2],
      "Death" = state_colors[3]
    ),
    labels = c("Dialysis", "Alive dialysis-free", "Death")
  ) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                              breaks = seq(0, 1, by = 0.1)) +
  ggplot2::labs(x = "Time (years)") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    legend.position = "bottom",
    # move legend below plot
    legend.direction = "horizontal",
    # horizontal layout
    legend.title = ggplot2::element_blank()
  )

# --- Arrange plots side by side ---
ggplot2::ggsave(
  plot = gridExtra::grid.arrange(p1, p2, ncol = 2),
  filename = paste0(results_path, "Main/Figure_2.pdf"),
  width = 10,
  height = 5,
  dpi = 300
)
