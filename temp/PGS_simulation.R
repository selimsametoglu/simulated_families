## ---------------------------
## PACKAGES
## ---------------------------
if (!requireNamespace("lme4", quietly = TRUE)) install.packages("lme4")
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")

library(lme4)
library(MASS)

set.seed(123)

## ---------------------------
## SIMULATE TRIO DATA
## ---------------------------
simulate_trio_data <- function(n_families,
                               beta_child   = 0.25,
                               beta_mother  = 0.10,
                               beta_father  = 0.10,
                               r_parent_pgs = 0.30,
                               resid_trait_sd = 1) {
  
  n_kids_per_family <- 2
  n_children <- n_families * n_kids_per_family
  
  # 1) Parental PGS (correlated)
  Sigma_parents <- matrix(c(1, r_parent_pgs,
                            r_parent_pgs, 1), nrow = 2)
  parent_pgs <- MASS::mvrnorm(
    n   = n_families,
    mu  = c(0, 0),
    Sigma = Sigma_parents
  )
  mother_pgs_family <- parent_pgs[, 1]
  father_pgs_family <- parent_pgs[, 2]
  
  # Expand to child level
  family_id  <- rep(1:n_families, each = n_kids_per_family)
  mother_pgs <- rep(mother_pgs_family, each = n_kids_per_family)
  father_pgs <- rep(father_pgs_family, each = n_kids_per_family)
  
  # 2) Child PGS
  mid_parent <- 0.5 * mother_pgs + 0.5 * father_pgs
  child_pgs_raw <- mid_parent + rnorm(n_children, 0, 1)
  
  # Standardize PGS
  z_child_pgs   <- as.numeric(scale(child_pgs_raw))
  z_mother_pgs  <- as.numeric(scale(mother_pgs))
  z_father_pgs  <- as.numeric(scale(father_pgs))
  
  # 3) Trait
  trait <- beta_child  * z_child_pgs +
    beta_mother * z_mother_pgs +
    beta_father * z_father_pgs +
    rnorm(n_children, 0, resid_trait_sd)
  
  data.frame(
    family_id    = family_id,
    z_child_pgs  = z_child_pgs,
    z_mother_pgs = z_mother_pgs,
    z_father_pgs = z_father_pgs,
    trait        = trait
  )
}

## ---------------------------
## POWER FOR MODEL 1 & MODEL 2
## ---------------------------
power_lmm_models <- function(n_families,
                             n_sims      = 1000,
                             alpha       = 0.05,
                             beta_child  = 0.25,
                             beta_mother = 0.10,
                             beta_father = 0.10) {
  
  p_child_m1   <- numeric(n_sims)
  p_child_m2   <- numeric(n_sims)
  p_mother_m2  <- numeric(n_sims)
  p_father_m2  <- numeric(n_sims)
  
  for (s in 1:n_sims) {
    dat <- simulate_trio_data(
      n_families  = n_families,
      beta_child  = beta_child,
      beta_mother = beta_mother,
      beta_father = beta_father
    )
    
    ## --- Model 1: trait ~ child PGS + (1 | family_id) ---
    m1 <- lmer(trait ~ z_child_pgs + (1 | family_id), data = dat)
    coefs1 <- summary(m1)$coefficients
    t_child_m1 <- coefs1["z_child_pgs", "t value"]
    p_child_m1[s] <- 2 * pnorm(abs(t_child_m1), lower.tail = FALSE)
    
    ## --- Model 2: trait ~ child + mother + father PGS + (1 | family_id) ---
    m2 <- lmer(trait ~ z_child_pgs + z_mother_pgs + z_father_pgs +
                 (1 | family_id), data = dat)
    coefs2 <- summary(m2)$coefficients
    t_child_m2  <- coefs2["z_child_pgs",  "t value"]
    t_mother_m2 <- coefs2["z_mother_pgs", "t value"]
    t_father_m2 <- coefs2["z_father_pgs", "t value"]
    
    p_child_m2[s]  <- 2 * pnorm(abs(t_child_m2),  lower.tail = FALSE)
    p_mother_m2[s] <- 2 * pnorm(abs(t_mother_m2), lower.tail = FALSE)
    p_father_m2[s] <- 2 * pnorm(abs(t_father_m2), lower.tail = FALSE)
  }
  
  data.frame(
    n_families      = n_families,
    power_m1_child  = mean(p_child_m1  < alpha, na.rm = TRUE),
    power_m2_child  = mean(p_child_m2  < alpha, na.rm = TRUE),
    power_m2_mother = mean(p_mother_m2 < alpha, na.rm = TRUE),
    power_m2_father = mean(p_father_m2 < alpha, na.rm = TRUE)
  )
}


## ---------------------------
## HELPER: R2 -> beta_child (Option 3)
## ---------------------------
beta_from_R2 <- function(R2) {
  sqrt(R2 / (1 - R2))
}



## ---------------------------
## RUN POWER FOR MULTIPLE R2s AND Ns
## ---------------------------

# Target R2 values for CHILD PGS
target_R2_child_vec <- c(0.15, 0.10, 0.05, 0.025)

# Grid of sample sizes (number of families)
n_families_grid <- c(200, 500, 1000, 2000)   # tweak as you like

# TRUE indirect effects for parents (set as you wish)
beta_mother_true <- 0.10
beta_father_true <- 0.10

# Number of simulation replicates per condition
n_sims <- 500   # increase to 1000+ for smoother estimates

results_all <- do.call(rbind, lapply(target_R2_child_vec, function(R2_child) {
  
  beta_child_true <- beta_from_R2(R2_child)
  
  do.call(rbind, lapply(n_families_grid, function(N) {
    out <- power_lmm_models(
      n_families  = N,
      n_sims      = n_sims,
      alpha       = 0.05,
      beta_child  = beta_child_true,
      beta_mother = beta_mother_true,
      beta_father = beta_father_true
    )
    
    # Add columns for scenario info
    out$target_R2_child <- R2_child
    out$beta_child_true <- beta_child_true
    out
  }))
}))

results_all


print(results_all)

# If you want a quick look:
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(dplyr)
library(tidyr)
library(ggplot2)

results_long <- results_all %>%
  pivot_longer(
    cols = starts_with("power_"),
    names_to = "effect",
    values_to = "power"
  )

ggplot(results_long,
       aes(x = n_families, y = power,
           color = effect, group = effect)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ target_R2_child, labeller = label_bquote(R^2 == .(target_R2_child))) +
  scale_x_continuous(breaks = n_families_grid) +
  ylim(0, 1) +
  theme_minimal() +
  labs(
    x = "Number of families (2 children per family)",
    y = "Estimated power",
    color = "Effect",
    title = "Power for Model 1 and Model 2 across child-PGS R² scenarios"
  )

