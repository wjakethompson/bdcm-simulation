### Capture command line arguments ---------------------------------------------
keep_cond <- c(5, 54)
frst_cond <- as.integer(keep_cond[1])
last_cond <- as.integer(keep_cond[2])


### Setup ----------------------------------------------------------------------
# install.packages(c("remotes", "pak"))
library(remotes)
library(pak)

# Package dependencies are managed with the renv package. To install the package
# versions used for this simulation, run renv::restore(), which will restore the
# packages from renv.lock.

# renv::restore()

library(tidyverse)
library(here)
library(fs)
library(glue)
library(modelr)
library(portableParallelSeeds)
library(measr)


### Utility functions ----------------------------------------------------------
logit <- function(x) {
  log(x / (1 - x))
}
inv_logit <- function(x) {
  exp(x) / (1 + exp(x))
}
qmatrix_entries <- function(att, item, all_entries) {
  entry <- if (item %in% 1:3) {
    all_entries %>%
      filter(!!sym(att) == 1L, total == 1L) %>%
      select(-c(total, prob))
  } else {
    all_entries %>%
      filter(!!sym(att) == 1L) %>%
      slice_sample(n = 1, weight_by = prob) %>%
      select(-c(total, prob))
  }

  return(entry)
}
generate_items <- function(q_matrix, generate) {
  all_param <- get_parameters(q_matrix, item_id = "item_id", type = "lcdm")
  intercepts <- all_param %>%
    filter(class == "intercept") %>%
    mutate(value = runif(n(), min = -3.00, max = 0.60))
  maineffects <- all_param %>%
    filter(class == "maineffect") %>%
    mutate(value = runif(n(), min = 1.00, max = 5.00))

  interactions <- if (generate == "lcdm") {
    all_param %>%
      filter(class == "interaction") %>%
      left_join(maineffects %>%
                  slice_min(value, n = 1, by = item_id) %>%
                  mutate(min_value = -1 * value) %>%
                  select(item_id, min_value),
                by = "item_id") %>%
      mutate(value = map_dbl(min_value, ~runif(1, min = .x, max = 2.00))) %>%
      select(-min_value)
  } else if (generate == "dina") {
    all_param %>%
      filter(class == "interaction") %>%
      mutate(value = runif(n(), min = 0.00, max = 5.00))
  }
  item_params <- full_join(all_param,
                           bind_rows(intercepts, maineffects, interactions),
                           by = c("item_id", "class", "attributes", "coef"))

  if (generate == "dina") {
    item_params <- item_params %>%
      mutate(num_att = case_when(is.na(attributes) ~ 0L,
                                 TRUE ~ str_count(attributes, "__") + 1)) %>%
      mutate(max_att = max(num_att), .by = item_id) %>%
      filter((num_att == 0) | (num_att == max_att)) %>%
      select(-c(num_att, max_att))
  }

  return(item_params)
}
pimat <- function(attributes, item_params) {
  possible_profiles <- create_profiles(attributes = attributes) %>%
    rowid_to_column(var = "class_id")
  class_params <- possible_profiles %>%
    select(-class_id) %>%
    modelr::model_matrix(as.formula(paste0("~ .^", attributes))) %>%
    rowid_to_column(var = "class_id") %>%
    pivot_longer(-class_id, names_to = "param", values_to = "relevant") %>%
    mutate(class = case_when(param == "(Intercept)" ~ "intercept",
                             str_detect(param, ":") ~ "interaction",
                             TRUE ~ "maineffect"),
           attributes = case_when(class == "intercept" ~ NA_character_,
                                  TRUE ~ str_replace(param, ":", "__"))) %>%
    select(class_id, class, attributes, relevant)

  pi_mat <- crossing(item_id = unique(item_params$item_id),
                     class_id = unique(class_params$class_id)) %>%
    left_join(item_params, by = "item_id", relationship = "many-to-many") %>%
    left_join(class_params, by = c("class_id", "class", "attributes"),
              relationship = "many-to-one") %>%
    filter(relevant == 1) %>%
    summarize(log_odds = sum(value), .by = c(item_id, class_id))

  return(pi_mat)
}
generate_data <- function(attributes, items, sample_size, generate) {
  # generate q-matrix -----
  all_entries <- create_profiles(attributes = attributes) %>%
    rowwise() %>%
    mutate(total = sum(c_across(where(is.integer)))) %>%
    ungroup() %>%
    filter(between(total, 1, 2)) %>%
    mutate(prob = case_when(total == 1 ~ 0.1,
                            TRUE ~ (0.9 / (attributes - 1))))

  q_matrix <- crossing(att = paste0("att", seq_len(attributes)),
           item = seq_len(items)) %>%
    pmap_dfr(qmatrix_entries, all_entries = all_entries) %>%
    mutate(item_id = paste0("item_", sprintf("%02d", 1:n())),
           .before = 1)

  # generate true profiles -----
  profiles <- create_profiles(attributes = attributes) %>%
    rowid_to_column(var = "class_id") %>%
    slice_sample(n = sample_size, replace = TRUE) %>%
    rowid_to_column(var = "resp_id")

  # generate true item parameters -----
  item_params <- generate_items(q_matrix = q_matrix, generate = generate)
  pi_matrix <- pimat(attributes = attributes, item_params = item_params) %>%
    mutate(prob = map_dbl(log_odds, inv_logit))

  # generate data -----
  data <- profiles %>%
    select(resp_id, class_id) %>%
    left_join(pi_matrix, by = "class_id", relationship = "many-to-many") %>%
    mutate(rand = runif(n(), min = 0, max = 1),
           score = case_when(rand < prob ~ 1L,
                             rand >= prob ~ 0L)) %>%
    select(resp_id, item_id, score) %>%
    mutate(item_id = paste0("item_", sprintf("%02d", item_id))) %>%
    pivot_wider(names_from = item_id, values_from = score)

  # return data -----
  ret_list <- list(
    data = data,
    q_matrix = q_matrix,
    true_person = profiles,
    true_pi = pi_matrix
  )
  return(ret_list)
}
calc_kappa <- function(score1, score2, min_score, max_score) {
    if (missing(min_score)) {
    min_score <- min(min(score1), min(score2))
  }
  if (missing(max_score)) {
    max_score <- max(max(score1), max(score2))
  }

  score1 <- factor(score1, levels = min_score:max_score)
  score2 <- factor(score2, levels = min_score:max_score)

  # pairwise frequencies
  confusion.mat <- table(data.frame(score1, score2))
  confusion.mat <- confusion.mat / sum(confusion.mat)

  # get expected pairwise frequencies under independence
  histogram.a <- table(score1) / length(table(score1))
  histogram.b <- table(score2) / length(table(score2))
  expected.mat <- histogram.a %*% t(histogram.b)
  expected.mat <- expected.mat / sum(expected.mat)

  # get weights
  labels <- as.numeric(as.vector(names(table(score1))))
  weights <- outer(labels, labels, FUN = function(x, y) (x - y) ^ 2)

  # calculate kappa
  kappa <- 1 - sum(weights * confusion.mat) / sum(weights * expected.mat)
  kappa
}
agreement <- function(model, threshold, data) {
  measr_extract(model, "attribute_prob") %>%
    mutate(resp_id = as.integer(as.character(resp_id)),
           across(starts_with("att"), ~case_when(.x >= threshold ~ 1L,
                                                 TRUE ~ 0L))) %>%
    pivot_longer(cols = starts_with("att"), names_to = "att",
                 values_to = "est") %>%
    left_join(data$true_person %>%
                select(-class_id) %>%
                pivot_longer(cols = starts_with("att"), names_to = "att",
                             values_to = "true"),
              by = c("resp_id", "att")) %>%
    summarize(pct_cor = mean(est == true),
              kappa = calc_kappa(est, true, min_score = 0, max_score = 1),
              type1 = mean(est == 1 & true == 0),
              type2 = mean(est == 0 & true == 1),
              precision = sum(est == 1 & true == 1) / sum(est == 1))
}
avg_prop <- function(values) {
  max999 <- function(x) sign(x) * min(0.999, abs(x))
  min001 <- function(x) sign(x) * max(0.001, abs(x))
  values <- purrr::map_dbl(.x = values, .f = max999)
  values <- purrr::map_dbl(.x = values, .f = min001)

  r2z <- function(x) 0.5 * log((1 + x) / (1 - x))
  z2r <- function(x) (exp(2 * x) - 1) / (exp(2 * x) + 1)
  values <- purrr::map_dbl(.x = values, .f = r2z)
  values <- mean(values)
  values <- z2r(values)
  values
}


### Simulation functions -------------------------------------------------------
# cond <- 1
# rep <- 4
# attributes <- 2
# items <- 5
# sample_size <- 500
# generate <- "dina"
# stan_cores <- 2
# seeds <- read_rds(here("data", "random-seeds.rds"))
single_sim <- function(cond, rep, attributes, items, sample_size, generate,
                       stan_cores, seeds) {
  # print progress message ------
  message <- glue("Now running: Condition {cond}, Replication {rep}")
  print(message)

  # create output directory -----
  outdir <- here("output", "replications", glue("cond{sprintf('%02d', cond)}"),
                 glue("rep_{sprintf('%03d', rep)}"))
  if (!dir_exists(outdir)) dir_create(outdir)

  if (file_exists(glue("{outdir}/results.rds"))) {
    data <- read_rds(glue("{outdir}/data.rds"))
    results <- read_rds(glue("{outdir}/results.rds"))

    ret_list <- list(sim_data = data,
                     results = results)

    return(ret_list)
  }

  # set seed -----
  setSeeds(seeds, run = rep)
  useStream(cond)

  # generate data -----
  if (!file_exists(glue("{outdir}/data.rds"))) {
    data <- generate_data(attributes = attributes, items = items,
                          sample_size = sample_size, generate = generate)
    write_rds(data, glue("{outdir}/data.rds"))
  } else {
    data <- read_rds(glue("{outdir}/data.rds"))
  }

  # estimate model -----
  lcdm <- suppressWarnings(
    measr_dcm(data = data$data, qmatrix = data$q_matrix,
              resp_id = "resp_id", item_id = "item_id",
              type = "lcdm", method = "mcmc", backend = "rstan",
              iter = 2000, warmup = 1500, chains = 2,
              cores = stan_cores, refresh = 0,
              control = list(adapt_delta = 0.99),
              prior = c(prior(normal(-1.5, 1), class = "intercept"),
                        prior(lognormal(0, 1), class = "maineffect"),
                        prior(normal(0, 2), class = "interaction")))
  )

  dina <- suppressWarnings(
    measr_dcm(data = data$data, qmatrix = data$q_matrix,
              resp_id = "resp_id", item_id = "item_id",
              type = "dina", method = "mcmc", backend = "rstan",
              iter = 2000, warmup = 1500, chains = 2,
              cores = stan_cores, refresh = 0)
  )

  # calculate fit measures -----
  lcdm <- add_fit(lcdm, method = c("m2", "ppmc"), item_fit = NULL)
  lcdm <- suppressWarnings(add_criterion(lcdm, cirterion = c("loo", "waic")))
  lcdm <- add_reliability(lcdm)

  dina <- add_fit(dina, method = c("m2", "ppmc"), item_fit = NULL)
  dina <- suppressWarnings(add_criterion(dina, cirterion = c("loo", "waic")))
  dina <- add_reliability(dina)

  # model comparisons -----
  loo <- loo_compare(lcdm, dina, criterion = "loo")
  loo_sig <- abs(loo[2, "elpd_diff"]) > (loo[2, "se_diff"] * 2.5)
  loo_prefer <- loo %>%
    as_tibble(rownames = "model") %>%
    slice(1) %>%
    pull(model)

  waic <- loo_compare(lcdm, dina, criterion = "waic")
  waic_sig <- abs(waic[2, "elpd_diff"]) > (waic[2, "se_diff"] * 2.5)
  waic_prefer <- waic %>%
    as_tibble(rownames = "model") %>%
    slice(1) %>%
    pull(model)

  # classification accuracy -----
  lcdm <- add_respondent_estimates(lcdm)
  dina <- add_respondent_estimates(dina)

  lcdm50 <- agreement(lcdm, threshold = 0.5, data = data)
  lcdm80 <- agreement(lcdm, threshold = 0.8, data = data)
  dina50 <- agreement(dina, threshold = 0.5, data = data)
  dina80 <- agreement(dina, threshold = 0.8, data = data)

  # summarize results -----
  results <- tibble(
    lcdm_m2 = lcdm$fit$m2$m2,
    lcdm_m2_pval = lcdm$fit$m2$pval,
    lcdm_m2_type1 = lcdm$fit$m2$pval < .05,
    lcdm_m2_type2 = NA,
    dina_m2 = dina$fit$m2$m2,
    dina_m2_pval = dina$fit$m2$pval,
    dina_m2_type1 = ifelse(generate == "dina", dina$fit$m2$pval < .05, NA),
    dina_m2_type2 = ifelse(generate == "lcdm", dina$fit$m2$pval > .05, NA),
    lcdm_ppmc_ppp = lcdm$fit$ppmc$model_fit$raw_score$ppp,
    lcdm_ppmc_type1 = lcdm$fit$ppmc$model_fit$raw_score$ppp < .05,
    lcdm_ppmc_type2 = NA,
    dina_ppmc_ppp = dina$fit$ppmc$model_fit$raw_score$ppp,
    dina_ppmc_type1 = ifelse(generate == "dina",
                             dina$fit$ppmc$model_fit$raw_score$ppp < .05,
                             NA),
    dina_ppmc_type2 = ifelse(generate == "lcdm",
                             dina$fit$ppmc$model_fit$raw_score$ppp > .05,
                             NA),
    loo_prefer = loo_prefer,
    loo_sig = loo_sig,
    waic_prefer = waic_prefer,
    waic_sig = waic_sig,
    true_prefer = generate,
    lcdm_attr_pctcor_50    = lcdm50$pct_cor,
    lcdm_attr_kappa_50     = lcdm50$kappa,
    lcdm_attr_type1_50     = lcdm50$type1,
    lcdm_attr_type2_50     = lcdm50$type2,
    lcdm_attr_precision_50 = lcdm50$precision,
    lcdm_attr_pctcor_80    = lcdm80$pct_cor,
    lcdm_attr_kappa_80     = lcdm80$kappa,
    lcdm_attr_type1_80     = lcdm80$type1,
    lcdm_attr_type2_80     = lcdm80$type2,
    lcdm_attr_precision_80 = lcdm80$precision,
    dina_attr_pctcor_50    = dina50$pct_cor,
    dina_attr_kappa_50     = dina50$kappa,
    dina_attr_type1_50     = dina50$type1,
    dina_attr_type2_50     = dina50$type2,
    dina_attr_precision_50 = dina50$precision,
    dina_attr_pctcor_80    = dina80$pct_cor,
    dina_attr_kappa_80     = dina80$kappa,
    dina_attr_type1_80     = dina80$type1,
    dina_attr_type2_80     = dina80$type2,
    dina_attr_precision_80 = dina80$precision,
    lcdm_reli_acc = mean(lcdm$reliability$map_reliability$accuracy$acc),
    lcdm_reli_con = mean(lcdm$reliability$map_reliability$consistency$consist),
    dina_reli_acc = mean(dina$reliability$map_reliability$accuracy$acc),
    dina_reli_con = mean(dina$reliability$map_reliability$consistency$consist)
  ) %>%
    add_column(rep = rep, .before = 1) %>%
    add_column(cond = cond, .before = 1) %>%
    write_rds(glue("{outdir}/results.rds"), compress = "gz")

  # return results -----
  ret_list <- list(sim_data = data,
                   results = results)

  rm(list = setdiff(ls(), c("ret_list"))); gc()
  return(ret_list)
}

# cond <- 1
# attributes <- 2
# items <- 5
# sample_size <- 500
# generate <- "lcdm"
# reps <- 3
# stan_cores <- 2
# seeds <- read_rds(here("data", "random-seeds.rds"))
single_cond <- function(cond, attributes, items, sample_size, generate, reps,
                        stan_cores, seeds) {
  # create output directories -----
  outdir <- here("output", "replications",
                 glue("cond{sprintf('%02d', cond)}"))
  if (!dir_exists(outdir)) dir_create(outdir)

  if (length(reps) == 1) {
    data_list <- tibble(cond = cond,
                        rep = seq_len(reps),
                        attributes = attributes,
                        items = items,
                        sample_size = sample_size,
                        generate = generate)

    outlist <- glue("{outdir}/rep_{sprintf('%03d', seq_len(reps))}")
  } else {
    data_list <- tibble(cond = cond,
                        rep = as.integer(reps),
                        attributes = attributes,
                        items = items,
                        sample_size = sample_size,
                        generate = generate)

    outlist <- glue("{outdir}/rep_{sprintf('%03d', as.integer(reps))}")
  }
  walk(outlist, ~ if(!dir_exists(.)) dir_create(.))

  # run a condition -----
  all_cond <- pmap(.l = as.list(data_list), .f = single_sim,
                   stan_cores = stan_cores, seeds = seeds)

  # summarize condition results -----
  all_cond_reps <- all_cond %>%
    map_dfr(pluck, "results") %>%
    rename_with(~str_replace(.x, "ppp_type", "type"))

  m2_res <- all_cond_reps %>%
    select(contains("m2_pval")) %>%
    mutate(generate = generate) %>%
    pivot_longer(cols = -generate, names_to = "estimate", values_to = "pval",
                 names_pattern = "(.*)_m2_pval") %>%
    mutate(true_flag = estimate == "dina" & generate != "dina",
           obs_flag = pval < .05) %>%
    summarize(m2_type1 = sum(obs_flag & !true_flag) / sum(!true_flag),
              m2_type2 = sum(!obs_flag & true_flag) / sum(true_flag),
              m2_precision = sum(obs_flag & true_flag) / sum(obs_flag))

  ppmc_res <- all_cond_reps %>%
    select(contains("ppmc_ppp")) %>%
    mutate(generate = generate) %>%
    pivot_longer(cols = -generate, names_to = "estimate", values_to = "ppp",
                 names_pattern = "(.*)_ppmc_ppp") %>%
    mutate(true_flag = estimate == "dina" & generate != "dina",
           obs_flag = ppp < .05) %>%
    summarize(ppmc_type1 = sum(obs_flag & !true_flag) / sum(!true_flag),
              ppmc_type2 = sum(!obs_flag & true_flag) / sum(true_flag),
              ppmc_precision = sum(obs_flag & true_flag) / sum(obs_flag))

  comp_res <- all_cond_reps %>%
    select(matches("loo|waic")) %>%
    summarize(loo_correct = mean(loo_prefer == generate),
              loo_correct_same = mean(loo_prefer == generate | !loo_sig),
              waic_correct = mean(waic_prefer == generate),
              waic_correct_same = mean(waic_prefer == generate | !waic_sig))

  attr_agree <- all_cond_reps %>%
    select(matches("_attr_")) %>%
    summarize(across(everything(), mean))

  attr_reli <- all_cond_reps %>%
    select(matches("_reli_")) %>%
    summarize(across(everything(), mean))

  # return results -----
  all_cond_reps <- all_cond_reps %>%
    mutate(attributes = attributes, items = items,
           sample_size = sample_size, generate = generate,
           .after = cond)
  cond_summary <- bind_cols(m2_res, ppmc_res, comp_res, attr_agree,
                            attr_reli) %>%
    mutate(cond = cond, attributes = attributes, items = items,
           sample_size = sample_size, generate = generate,
           .before = 1)

  ret_list <- list(
    all_reps = all_cond_reps,
    summary = cond_summary
  )
  write_rds(ret_list, glue("{outdir}/results.rds"), compress = "gz")

  return(ret_list)
}

# conditions <- tibble(cond = 1:2, attributes = 2, items = 5, sample_size = 500,
#                      generate = c("lcdm", "dina"))
# reps <- 3
# stan_cores <- 2
# seeds <- read_rds(here("data", "random-seeds.rds"))
full_sim <- function(conditions, reps, stan_cores, seeds) {
  write_path <- here("output", "replications")
  if (!dir_exists(write_path)) dir_create(write_path)

  sim_results <- pmap(conditions, .f = single_cond, reps = reps,
                      stan_cores = stan_cores, seeds = seeds)

  write_rds(sim_results,
            file = here("output", "replications", "full_results.rds"),
            compress = "gz")

  ret_list <- list(
    all_reps = map_df(sim_results, "all_reps"),
    rep_summary = map_df(sim_results, "summary")
  )

  return(ret_list)
}


### Run simulation -------------------------------------------------------------
conditions <- crossing(attributes = 2:4,
                       items = c(5, 7, 10),
                       sample_size = c(500, 1000, 5000),
                       generate = c("lcdm", "dina")) %>%
  rowid_to_column(var = "cond")
reps <- 100

# Create seeds
if (file_exists(here("data", "random-seeds.rds"))) {
  projSeeds <- read_rds(here("data", "random-seeds.rds"))
} else {
  projSeeds <- seedCreator(nReps = reps,
                           streamsPerRep = nrow(conditions),
                           seed = 1213, file = here("data", "random-seeds.rds"))
}

# Filter conditions
if (length(keep_cond) == 0) {
  frst_cond <- min(conditions$cond)
  last_cond <- max(conditions$cond)
}
conditions <- conditions %>%
  filter(between(cond, frst_cond, last_cond))

# Run simulation
sim_results <- full_sim(conditions = conditions, reps = reps, stan_cores = 2,
                        seeds = projSeeds)

# Save output
write_rds(sim_results,
          file = here("output", "sim_results.rds"),
          compress = "gz")
write_rds(sim_results$all_reps,
          file = here("output", "all_reps.rds"),
          compress = "gz")
write_rds(sim_results$rep_summary,
          file = here("output", "rep_summary.rds"),
          compress = "gz")
