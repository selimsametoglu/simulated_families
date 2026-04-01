## ============================
## PACKAGES
## ============================
pkgs <- c("lme4", "MASS", "dplyr", "tidyr", "ggplot2")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(lme4)
library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(123)

## ============================
## OPTION 3: R2 -> beta
## beta = sqrt(R2/(1-R2)) when residual SD = 1 in the DGM
## ============================
beta_from_R2 <- function(R2) sqrt(R2 / (1 - R2))

## ============================
## SIMULATE TRIO + TWIN DATA (MZ/DZ)
## - One twin pair per family (2 children)
## - MZ/DZ affects within-pair correlation of child PGS
## - Trait generated from child/mother/father PGS
## ============================
simulate_trio_twin_data <- function(n_families,
                                    beta_child   = 0.25,
                                    beta_mother  = 0.10,
                                    beta_father  = 0.10,
                                    r_parent_pgs = 0.30,
                                    resid_trait_sd = 1,
                                    prop_MZ = 0.5,
                                    mz_twin_cor = 1.0,
                                    dz_twin_cor = 0.5) {
  
  n_kids_per_family <- 2
  n_children <- n_families * n_kids_per_family
  
  # Twin type per family
  twin_type_family <- ifelse(runif(n_families) < prop_MZ, "MZ", "DZ")
  
  # Parental PGS (correlated)
  Sigma_parents <- matrix(c(1, r_parent_pgs,
                            r_parent_pgs, 1), nrow = 2)
  parent_pgs <- MASS::mvrnorm(n_families, mu = c(0, 0), Sigma = Sigma_parents)
  
  mother_pgs_family <- parent_pgs[, 1]
  father_pgs_family <- parent_pgs[, 2]
  
  # Expand to child level
  family_id  <- rep(1:n_families, each = 2)
  child_id   <- rep(1:2, times = n_families)
  mother_pgs <- rep(mother_pgs_family, each = 2)
  father_pgs <- rep(father_pgs_family, each = 2)
  twin_type  <- rep(twin_type_family, each = 2)
  
  # Child PGS: mid-parent + correlated residuals (cor depends on MZ/DZ)
  child_pgs_raw <- numeric(n_children)
  
  for (f in 1:n_families) {
    idx <- which(family_id == f)
    mp  <- 0.5 * mother_pgs_family[f] + 0.5 * father_pgs_family[f]
    
    cor_tw <- if (twin_type_family[f] == "MZ") mz_twin_cor else dz_twin_cor
    Sigma_tw <- matrix(c(1, cor_tw,
                         cor_tw, 1), nrow = 2)
    
    e_pair <- MASS::mvrnorm(n = 1, mu = c(0, 0), Sigma = Sigma_tw)
    child_pgs_raw[idx] <- mp + e_pair
  }
  
  # Standardize PGS (good practice for interpretability)
  z_child_pgs   <- as.numeric(scale(child_pgs_raw))
  z_mother_pgs  <- as.numeric(scale(mother_pgs))
  z_father_pgs  <- as.numeric(scale(father_pgs))
  
  # Outcome
  trait <- beta_child  * z_child_pgs +
    beta_mother * z_mother_pgs +
    beta_father * z_father_pgs +
    rnorm(n_children, 0, resid_trait_sd)
  
  data.frame(
    family_id    = family_id,
    child_id     = child_id,
    twin_type    = factor(twin_type, levels = c("DZ", "MZ")),
    z_child_pgs  = z_child_pgs,
    z_mother_pgs = z_mother_pgs,
    z_father_pgs = z_father_pgs,
    trait        = trait
  )
}

## ============================
## POWER FUNCTION FOR MODEL 1 & MODEL 2
## Model 1: trait ~ child + (1|family)
## Model 2: trait ~ child + mother + father + (1|family)
## ============================
power_lmm_models_twins <- function(n_families,
                                   n_sims = 500,
                                   alpha = 0.05,
                                   beta_child  = 0.25,
                                   beta_mother = 0.10,
                                   beta_father = 0.10,
                                   prop_MZ = 0.5,
                                   mz_twin_cor = 1.0,
                                   dz_twin_cor = 0.5,
                                   r_parent_pgs = 0.30,
                                   resid_trait_sd = 1) {
  
  p_child_m1   <- numeric(n_sims)
  p_child_m2   <- numeric(n_sims)
  p_mother_m2  <- numeric(n_sims)
  p_father_m2  <- numeric(n_sims)
  
  for (s in 1:n_sims) {
    dat <- simulate_trio_twin_data(
      n_families     = n_families,
      beta_child     = beta_child,
      beta_mother    = beta_mother,
      beta_father    = beta_father,
      prop_MZ        = prop_MZ,
      mz_twin_cor    = mz_twin_cor,
      dz_twin_cor    = dz_twin_cor,
      r_parent_pgs   = r_parent_pgs,
      resid_trait_sd = resid_trait_sd
    )
    
    # Model 1
    m1 <- lmer(trait ~ z_child_pgs + (1 | family_id), data = dat, REML = FALSE)
    coefs1 <- summary(m1)$coefficients
    t_child_m1 <- coefs1["z_child_pgs", "t value"]
    p_child_m1[s] <- 2 * pnorm(abs(t_child_m1), lower.tail = FALSE)
    
    # Model 2
    m2 <- lmer(trait ~ z_child_pgs + z_mother_pgs + z_father_pgs + (1 | family_id),
               data = dat, REML = FALSE)
    coefs2 <- summary(m2)$coefficients
    t_child_m2  <- coefs2["z_child_pgs",  "t value"]
    t_mother_m2 <- coefs2["z_mother_pgs", "t value"]
    t_father_m2 <- coefs2["z_father_pgs", "t value"]
    
    p_child_m2[s]  <- 2 * pnorm(abs(t_child_m2),  lower.tail = FALSE)
    p_mother_m2[s] <- 2 * pnorm(abs(t_mother_m2), lower.tail = FALSE)
    p_father_m2[s] <- 2 * pnorm(abs(t_father_m2), lower.tail = FALSE)
  }
  
  tibble(
    n_families      = n_families,
    power_m1_child  = mean(p_child_m1  < alpha, na.rm = TRUE),
    power_m2_child  = mean(p_child_m2  < alpha, na.rm = TRUE),
    power_m2_mother = mean(p_mother_m2 < alpha, na.rm = TRUE),
    power_m2_father = mean(p_father_m2 < alpha, na.rm = TRUE)
  )
}

## ============================
## RUN POWER GRID (R2 scenarios x N scenarios)
## ============================

# Your requested R2 values for child PGS (Option 3)
target_R2_child_vec <- c(0.15, 0.10, 0.05, 0.025)

# Choose family sample sizes (each family contributes exactly 2 twins)
n_families_grid <- c(200, 500, 1000, 2000)

# Simulation settings
n_sims <- 500          # raise to 1000–3000 for smoother curves
alpha  <- 0.05

# Twin settings (tweak these)
prop_MZ     <- 0.50
mz_twin_cor <- 1.00
dz_twin_cor <- 0.50

# Parental indirect effects (choose what you expect)
beta_mother_true <- 0.10
beta_father_true <- 0.10

results_all <- bind_rows(lapply(target_R2_child_vec, function(R2_child) {
  
  beta_child_true <- beta_from_R2(R2_child)
  
  bind_rows(lapply(n_families_grid, function(N) {
    out <- power_lmm_models_twins(
      n_families  = N,
      n_sims      = n_sims,
      alpha       = alpha,
      beta_child  = beta_child_true,
      beta_mother = beta_mother_true,
      beta_father = beta_father_true,
      prop_MZ     = prop_MZ,
      mz_twin_cor = mz_twin_cor,
      dz_twin_cor = dz_twin_cor
    )
    
    out %>%
      mutate(
        target_R2_child = R2_child,
        beta_child_true = beta_child_true,
        prop_MZ = prop_MZ,
        mz_twin_cor = mz_twin_cor,
        dz_twin_cor = dz_twin_cor
      )
  }))
}))

print(results_all)

## ============================
## GRAPH 1: All effects, faceted by target R2
## ============================

results_long <- results_all %>%
  pivot_longer(
    cols = starts_with("power_"),
    names_to = "effect",
    values_to = "power"
  ) %>%
  mutate(
    effect = recode(effect,
                    power_m1_child  = "Model 1: Direct (child PGS)",
                    power_m2_child  = "Model 2: Direct (child PGS)",
                    power_m2_mother = "Model 2: Indirect (mother PGS)",
                    power_m2_father = "Model 2: Indirect (father PGS)"
    ),
    target_R2_child = factor(target_R2_child,
                             levels = c(0.15, 0.10, 0.05, 0.025),
                             labels = c("R²=0.15", "R²=0.10", "R²=0.05", "R²=0.025"))
  )

p_all <- ggplot(results_long, aes(x = n_families, y = power, group = effect, color = effect)) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ target_R2_child, nrow = 2) +
  scale_x_continuous(breaks = n_families_grid) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  labs(
    title = "Power curves for Model 1 and Model 2 (MZ/DZ twin families)",
    subtitle = paste0("Dashed line = 80% power | prop_MZ=", prop_MZ,
                      ", MZ cor=", mz_twin_cor, ", DZ cor=", dz_twin_cor,
                      " | sims per point=", n_sims),
    x = "Number of families (each family = 1 twin pair)",
    y = "Estimated power",
    color = "Effect"
  )

print(p_all)

## ============================
## GRAPH 2: Direct effect only (Model 1 vs Model 2)
## ============================

direct_only <- results_long %>%
  filter(effect %in% c("Model 1: Direct (child PGS)", "Model 2: Direct (child PGS)"))

p_direct <- ggplot(direct_only, aes(x = n_families, y = power, group = effect, color = effect)) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ target_R2_child, nrow = 2) +
  scale_x_continuous(breaks = n_families_grid) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  labs(
    title = "Power for the direct effect: Model 1 vs Model 2",
    x = "Number of families",
    y = "Estimated power",
    color = "Model"
  )

print(p_direct)

## ============================
## GRAPH 3: Indirect effects only (mother vs father)
## ============================

indirect_only <- results_long %>%
  filter(effect %in% c("Model 2: Indirect (mother PGS)", "Model 2: Indirect (father PGS)"))

p_indirect <- ggplot(indirect_only, aes(x = n_families, y = power, group = effect, color = effect)) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ target_R2_child, nrow = 2) +
  scale_x_continuous(breaks = n_families_grid) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  labs(
    title = "Power for indirect effects in Model 2 (mother vs father)",
    x = "Number of families",
    y = "Estimated power",
    color = "Effect"
  )

print(p_indirect)
