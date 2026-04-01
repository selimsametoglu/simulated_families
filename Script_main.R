# Power analyses # 
# a.k.a. PGS_simulation_v4 (with MZDZ)
# Selim Sametoglu #


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
## SIMULATE TWIN-TRIO DATA (MZ/DZ)
## - One twin pair per family (2 children)
## - Parental PGS correlated (assortative mating)
## - Child PGS correlated within twin pair (MZ=1, DZ=0.5 by default)
## - Trait depends on child PGS + mother PGS + father PGS (truth)
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
  
  n_children <- n_families * 2
  
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
  
  # Child PGS: mid-parent + correlated residuals within twin pair
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
  
  # Standardize PGS (so betas are in SD units)
  z_child_pgs   <- as.numeric(scale(child_pgs_raw))
  z_mother_pgs  <- as.numeric(scale(mother_pgs))
  z_father_pgs  <- as.numeric(scale(father_pgs))
  
  # Outcome generated from truth
  trait <- beta_child  * z_child_pgs +
    beta_mother * z_mother_pgs +
    beta_father * z_father_pgs +
    rnorm(n_children, 0, resid_trait_sd)
  
  dat <- data.frame(
    family_id    = family_id,
    child_id     = child_id,
    twin_type    = factor(twin_type, levels = c("DZ", "MZ")),
    z_child_pgs  = z_child_pgs,
    z_mother_pgs = z_mother_pgs,
    z_father_pgs = z_father_pgs,
    trait        = trait
  )
  
  # Shared parental predictor (strategy)
  dat$z_parent_mean <- as.numeric(scale((dat$z_mother_pgs + dat$z_father_pgs) / 2))
  
  dat
}

## ============================
## POWER FUNCTION (shared parental strategy)
## Model 1: trait ~ child + (1|family)
## Model 2: trait ~ child + parent_mean + (1|family)
##
## Power outputs:
## - direct effect in Model 1
## - direct effect in Model 2
## - indirect effect (parent_mean) in Model 2
## - LRT power for Model2 vs Model1 (chi-square test; df=1)
## ============================
power_shared_parent_strategy <- function(n_families,
                                         n_sims = 500,
                                         alpha = 0.008, # tweak alpha from here, too
                                         beta_child  = 0.25,
                                         beta_indirect_each_parent = 0.10,
                                         prop_MZ = 0.5,
                                         mz_twin_cor = 1.0,
                                         dz_twin_cor = 0.5,
                                         r_parent_pgs = 0.30,
                                         resid_trait_sd = 1) {
  
  p_child_m1 <- numeric(n_sims)
  p_child_m2 <- numeric(n_sims)
  p_ind_m2   <- numeric(n_sims)
  p_lrt      <- numeric(n_sims)
  
  for (s in 1:n_sims) {
    
    dat <- simulate_trio_twin_data(
      n_families     = n_families,
      beta_child     = beta_child,
      beta_mother    = beta_indirect_each_parent,
      beta_father    = beta_indirect_each_parent,
      prop_MZ        = prop_MZ,
      mz_twin_cor    = mz_twin_cor,
      dz_twin_cor    = dz_twin_cor,
      r_parent_pgs   = r_parent_pgs,
      resid_trait_sd = resid_trait_sd
    )
    
    # ML fits for LRT validity
    m1 <- lmer(trait ~ z_child_pgs + (1 | family_id), data = dat, REML = FALSE)
    m2 <- lmer(trait ~ z_child_pgs + z_parent_mean + (1 | family_id), data = dat, REML = FALSE)
    
    # Wald p-values (normal approx) for fixed effects
    c1 <- summary(m1)$coefficients
    c2 <- summary(m2)$coefficients
    
    t_child_m1 <- c1["z_child_pgs", "t value"]
    p_child_m1[s] <- 2 * pnorm(abs(t_child_m1), lower.tail = FALSE)
    
    t_child_m2 <- c2["z_child_pgs", "t value"]
    p_child_m2[s] <- 2 * pnorm(abs(t_child_m2), lower.tail = FALSE)
    
    t_ind_m2 <- c2["z_parent_mean", "t value"]
    p_ind_m2[s] <- 2 * pnorm(abs(t_ind_m2), lower.tail = FALSE)
    
    # Likelihood ratio test (chi-square, df=1)
    lrt <- anova(m1, m2)
    p_lrt[s] <- lrt$`Pr(>Chisq)`[2]
  }
  
  tibble(
    n_families      = n_families,
    power_m1_direct = mean(p_child_m1 < alpha, na.rm = TRUE),
    power_m2_direct = mean(p_child_m2 < alpha, na.rm = TRUE),
    power_m2_ind    = mean(p_ind_m2   < alpha, na.rm = TRUE),
    power_lrt_m2_vs_m1 = mean(p_lrt   < alpha, na.rm = TRUE)
  )
}

## ============================
## RUN POWER GRID
## - Direct R2 scenarios for child PGS
## - Indirect R2 scenarios for EACH parent (assumed equal)
## ============================

# Direct-effect scenarios (child PGS R² in isolation)
target_R2_child_vec <- c(0.15, 0.10, 0.05, 0.025, 0.015, 0.01, 0.005)

# Indirect-effect scenarios for each parent (mother = father)
target_R2_indirect_vec <- c(0.15, 0.10, 0.05, 0.025, 0.015, 0.01, 0.005)

# Sample sizes: number of families (each family contributes 2 twins)
n_families_grid <- c(200, 500, 1000, 2000)

# Simulation settings
n_sims <- 500
alpha  <- 0.008 # tweak the alpha from here

# Twin settings
prop_MZ     <- 0.50
mz_twin_cor <- 1.00
dz_twin_cor <- 0.50

results_all <- bind_rows(lapply(target_R2_child_vec, function(R2_child) {
  
  beta_child_true <- beta_from_R2(R2_child)
  
  bind_rows(lapply(target_R2_indirect_vec, function(R2_ind) {
    
    beta_ind_each_parent <- beta_from_R2(R2_ind)
    
    bind_rows(lapply(n_families_grid, function(N) {
      
      out <- power_shared_parent_strategy(
        n_families  = N,
        n_sims      = n_sims,
        alpha       = alpha,
        beta_child  = beta_child_true,
        beta_indirect_each_parent = beta_ind_each_parent,
        prop_MZ     = prop_MZ,
        mz_twin_cor = mz_twin_cor,
        dz_twin_cor = dz_twin_cor
      )
      
      out %>%
        mutate(
          target_R2_child    = R2_child,
          target_R2_indirect = R2_ind,
          beta_child_true    = beta_child_true,
          beta_ind_each_parent = beta_ind_each_parent
        )
    }))
  }))
}))

print(results_all)

## ============================
## PLOT GRAPHS
## ============================

results_long <- results_all %>%
  pivot_longer(
    cols = starts_with("power_"),
    names_to = "metric",
    values_to = "power"
  ) %>%
  mutate(
    metric = recode(metric,
                    power_m1_direct      = "Model 1: Total genetic effect Wald test",
                    power_m2_direct      = "Model 2: Direct genetic effect Wald test",
                    power_m2_ind         = "Model 2: Indirect genetic effect (parent mean) Wald test",
                    power_lrt_m2_vs_m1   = "Model 2 vs Model 1: LRT (χ², df=1)"
    ),
    child_R2_lab = factor(
      target_R2_child,
      levels = c(0.005, 0.01, 0.015, 0.025, 0.05, 0.10, 0.15),
      labels = c("Total R²=0.005",
                 "Total R²=0.01",
                 "Total R²=0.015",
                 "Total R²=0.025",
                 "Total R²=0.05",
                 "Total R²=0.10",
                 "Total R²=0.15")
    ),
    ind_R2_lab = factor(
      target_R2_indirect,
      levels = c(0.15, 0.10, 0.05, 0.025, 0.015, 0.01, 0.005),
      labels = c("Ind. R²=0.15",
                 "Ind. R²=0.10",
                 "Ind. R²=0.05",
                 "Ind. R²=0.025",
                 "Ind. R²=0.015",
                 "Ind. R²=0.01",
                 "Ind. R²=0.005")
    )
  )



# Plot 1: focus on the individual coefficients

library(ggplot2)

results_no_lrt <- results_long %>%
  dplyr::filter(metric != "Model 2 vs Model 1: LRT (χ², df=1)")

pd <- position_dodge(width = 40)  # adjust depending on your x-scale spacing

p_all <- ggplot(
  results_no_lrt,
  aes(x = n_families, y = power, group = metric, color = metric,
      linetype = metric, shape = metric)
) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(size = 2, position = pd) +
  facet_grid(ind_R2_lab ~ child_R2_lab) +
  scale_x_continuous(breaks = n_families_grid) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  labs(
    title = "Power for shared-parent indirect strategy (Wald tests only)",
    subtitle = paste0(
      "Dashed line = 80% power | sims/point=", n_sims,
      " | prop_MZ=", prop_MZ, ", MZ cor=", mz_twin_cor, ", DZ cor=", dz_twin_cor
    ),
    x = "Number of families (each family = 1 twin pair)",
    y = "Estimated power",
    color = "Test",
    linetype = "Test",
    shape = "Test"
  )

print(p_all)



# to save the plot
ggsave(
  filename = "Plot_poweranalysis_alpha0.008.pdf",
  plot = p_all,
  device = cairo_pdf,
  width = 16, height = 10, units = "in",
  dpi = 300
)



# Plot 2: focus on LRT power (often your primary indirect-effect test)
lrt_only <- results_long %>% filter(metric == "Model 2 vs Model 1: LRT (χ², df=1)")

p_lrt <- ggplot(lrt_only, aes(x = n_families, y = power, group = 1)) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(ind_R2_lab ~ child_R2_lab) +
  scale_x_continuous(breaks = n_families_grid) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Power for indirect effect (Model 2 vs Model 1 likelihood-ratio test)",
    x = "Number of families",
    y = "Estimated power"
  )

print(p_lrt)



