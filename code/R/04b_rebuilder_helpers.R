# =============================================================
# SBAPC smoking-to-mortality replication code
# Rebuilder helpers for the sequential pipeline
# =============================================================

#' Recover Demographic Keys (Age, Period, Cohort)
#' @description Ensures a data frame has consistent age, period, and cohort columns.
#' Follows the identity: cohort = period - age.
recover_demographic_keys <- function(df) {
  if (is.null(df)) return(NULL)
  df <- tibble::as_tibble(df)
  
  # Ensure period is present
  if (!"period" %in% names(df)) stop("Missing 'period' column in demographic recovery.")
  
  # Try to find cohort
  ch_col <- intersect(c("cohort", "cohort_true", "cohort_ref"), names(df))[1]
  
  if (!is.na(ch_col)) {
    # If we have cohort, derive age
    df$cohort <- as.integer(as.character(df[[ch_col]]))
    if (!"age" %in% names(df)) {
      df$age <- as.integer(df$period - df$cohort)
    }
  } else if ("age" %in% names(df)) {
    # If we have age, derive cohort
    df$age <- as.integer(as.character(df$age))
    df$cohort <- as.integer(df$period - df$age)
  } else {
    stop("Cannot recover demographic keys: missing both age and cohort.")
  }
  
  return(df)
}

#' Get Mortality Kernel Weights
#' @description Annualizes 3 post-diagnosis death probabilities (0-1, 1-3, 3-5 years)
#' into a 6-lag kernel (w0...w5) using the midyear_uniform rule (Manuscript SM-A.3).
get_mortality_kernel_weights <- function(p01, p13, p35) {
  # Weights according to SM-A.3
  w <- numeric(6)
  w[1] <- 0.5 * p01                      # w0
  w[2] <- 0.5 * p01 + 0.25 * p13         # w1
  w[3] <- 0.5 * p13                      # w2
  w[4] <- 0.25 * p13 + 0.25 * p35        # w3
  w[5] <- 0.5 * p35                      # w4
  w[6] <- 0.25 * p35                     # w5
  return(w)
}

#' Calculate Incidence Sensitivity Coefficient (bz_hat)
#' @description Returns the sensitivity coefficient based on BETA_MODE.
get_incidence_sensitivity_coef <- function(beta_mode = "fixed_rr_offset", beta_P_eff = NULL) {
  return(1.0)
}

#' Apply Scenario Incidence Shift
#' @description Rebuilds incidence by applying the change in epidemiologic offset
#' scaled by the sensitivity coefficient.
apply_incidence_scenario_shift <- function(base_rate, 
                                           off_scen, 
                                           off_base, 
                                           bz_hat = 1.0) {
  # Delta Offset logic: exp(bz_hat * (off_scen - off_base))
  delta_off <- .safe_num(off_scen) - .safe_num(off_base)
  log_ratio <- exp(bz_hat * delta_off)
  
  new_rate <- pmax(base_rate * log_ratio, 1e-12)
  return(new_rate)
}

#' Rebuild Expected Deaths Surface (Convolution)
#' @description Implementation of the discrete convolution D_{a,t} = sum w_l * I_{a-l, t-l}.
#' Strictly cohort-consistent (SM-A.3).
rebuild_expected_deaths_surface <- function(inc_df, weights) {
  # inc_df must have sex, age, period, cases_hat
  # weights is a vector of length 6 (w0...w5)
  
  if (is.null(inc_df) || nrow(inc_df) == 0) return(NULL)
  
  # Ensure sorted for faster lookups if needed, though we'll use a join/shift approach
  inc_df <- recover_demographic_keys(inc_df)
  
  # Create shifted copies for convolution
  conv_list <- list()
  for (l in 0:5) {
    w_l <- weights[l+1]
    if (w_l <= 0) next
    
    # Lag incident cases along cohort lines: diagnosis at (a-l, t-l) leads to death at (a,t)
    shifted <- inc_df %>%
      dplyr::mutate(
        age = age + l,
        period = period + l
      ) %>%
      dplyr::mutate(cases_w = cases_hat * w_l) %>%
      dplyr::select(sex, age, period, cases_w)
    
    conv_list[[l+1]] <- shifted
  }
  
  # Aggregate shifted contributions
  D_at <- dplyr::bind_rows(conv_list) %>%
    dplyr::group_by(sex, age, period) %>%
    dplyr::summarise(mort_ext_deaths = sum(cases_w, na.rm = TRUE), .groups = "drop")
  
  return(D_at)
}
