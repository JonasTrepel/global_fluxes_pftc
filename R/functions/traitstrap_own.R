kurtosis <- function (x, na.rm = FALSE) 
{
  if (is.matrix(x)) 
    apply(x, 2, kurtosis, na.rm = na.rm)
  else if (is.vector(x)) {
    if (na.rm) 
      x <- x[!is.na(x)]
    n <- length(x)
    n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
  }
  else if (is.data.frame(x)) 
    sapply(x, kurtosis, na.rm = na.rm)
  else kurtosis(as.vector(x), na.rm = na.rm)
}

skewness <- function (x, na.rm = FALSE) 
{
  if (is.matrix(x)) 
    apply(x, 2, skewness, na.rm = na.rm)
  else if (is.vector(x)) {
    if (na.rm) 
      x <- x[!is.na(x)]
    n <- length(x)
    (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
  }
  else if (is.data.frame(x)) 
    sapply(x, skewness, na.rm = na.rm)
  else skewness(as.vector(x), na.rm = na.rm)
}

get_ci <- function(data, sd_mult = 1, ci = 0.95, which, parametric = TRUE) {
  if (isTRUE(parametric)) {
    if (which == "high") {
      return(mean(data) + sd(data) * sd_mult)
    }
    if (which == "low") {
      return(mean(data) - sd(data) * sd_mult)
    }
  } else {
    if (which == "high") {
      return(quantile(data, probs = (1 + ci) / 2, type = 1))
    }
    if (which == "low") {
      return(quantile(data, probs = (1 - ci) / 2, type = 1))
    }
  }
}

trait_np_bootstrap_own <- function (filled_traits, nrep = 100, sample_size = 200, raw = FALSE){
  if (isTRUE(raw)) {
    nrep <- 1
  }
  attrib <- attr(filled_traits, "attrib")
  value_col <- attrib$value_col
  bootstrap_moments <- list_rbind(map(seq_len(nrep), ~{
    raw_dist <- slice_sample(filled_traits, n = sample_size, 
                             replace = TRUE, weight_by = weight)
    if (raw) {
      return(raw_dist)
    } else {
      summarise(raw_dist, mean = mean(.data[[value_col]]), sd = sd(.data[[value_col]]), 
                variance = var(.data[[value_col]]), skewness = skewness(.data[[value_col]]), 
                kurtosis = kurtosis(.data[[value_col]]), .groups = "drop")
    }
  }, .id = "n"))
  attr(bootstrap_moments, "attrib") <- attrib
  class(bootstrap_moments) <- class(bootstrap_moments)[!class(bootstrap_moments) == 
                                                         "filled_trait"]
  return(bootstrap_moments)
}


trait_summarise_boot_moments_own <- function (bootstrap_moments, parametric = TRUE, sd_mult = 1, 
          ci = 0.95) 
{
  attrib <- attr(bootstrap_moments, "attrib")
  groups <- c(as.character(attrib$scale_hierarchy), attrib$trait_col, 
              attrib$other_col)
  if (!is.null(attrib$treatment_col)) {
    groups <- c(groups, paste0(attrib$treatment_col, "_comm"))
  }
  summ_bootstrap_moments <- summarise(rename(group_by(ungroup(bootstrap_moments), 
                                                      across(any_of(groups))), MEAN = "mean"), n = n(), mean = mean(.data$MEAN), 
                                      ci_low_mean = get_ci(data = .data$MEAN, sd_mult = sd_mult, 
                                                           ci = ci, which = "low", parametric = parametric), 
                                      ci_high_mean = get_ci(data = .data$MEAN, sd_mult = sd_mult, 
                                                            ci = ci, which = "high", parametric = parametric), 
                                      sdd = mean(.data$sd), ci_low_sd = get_ci(data = .data$sd, 
                                                                                      sd_mult = sd_mult, ci = ci, which = "low", parametric = parametric), 
                                      ci_high_sd = get_ci(data = .data$sd, sd_mult = sd_mult, 
                                                           ci = ci, which = "high", parametric = parametric), 
                                      var = mean(.data$variance), ci_low_var = get_ci(data = .data$variance, 
                                                                                      sd_mult = sd_mult, ci = ci, which = "low", parametric = parametric), 
                                      ci_high_var = get_ci(data = .data$variance, sd_mult = sd_mult, 
                                                           ci = ci, which = "high", parametric = parametric), 
                                      skew = mean(.data$skewness), ci_low_skew = get_ci(data = .data$skewness, 
                                                                                        sd_mult = sd_mult, ci = ci, which = "low", parametric = parametric), 
                                      ci_high_skew = get_ci(data = .data$skewness, sd_mult = sd_mult, 
                                                            ci = ci, which = "high", parametric = parametric), 
                                      kurt = mean(.data$kurtosis), ci_low_kurt = get_ci(data = .data$kurtosis, 
                                                                                        sd_mult = sd_mult, ci = ci, which = "low", parametric = parametric), 
                                      ci_high_kurt = get_ci(data = .data$kurtosis, sd_mult = sd_mult, 
                                                            ci = ci, which = "high", parametric = parametric)) %>% 
    rename(sd = sdd)
  return(summ_bootstrap_moments)
}
