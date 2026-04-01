## ============================
## RUN POWER GRID:
## child R2 scenarios  x  indirect R2 scenarios  x  N scenarios
## ============================

# Child direct-effect R2 scenarios
target_R2_child_vec <- c(0.15, 0.10, 0.05, 0.025)

# Indirect-effect R2 scenarios (same for mother and father)
target_R2_indirect_vec <- c(0.025, 0.05, 0.10)

# Family sample sizes (each family = 1 twin pair = 2 children)
n_families_grid <- c(200, 500, 1000, 2000)

# Sim settings
n_sims <- 500
alpha  <- 0.05

# Twin settings (as before)
prop_MZ     <- 0.50
mz_twin_cor <- 1.00
dz_twin_cor <- 0.50

results_all <- bind_rows(lapply(target_R2_child_vec, function(R2_child) {
  
  beta_child_true <- beta_from_R2(R2_child)
  
  bind_rows(lapply(target_R2_indirect_vec, function(R2_ind) {
    
    beta_ind_true <- beta_from_R2(R2_ind)  # same for mother and father
    
    bind_rows(lapply(n_families_grid, function(N) {
      
      out <- power_lmm_models_twins(
        n_families  = N,
        n_sims      = n_sims,
        alpha       = alpha,
        beta_child  = beta_child_true,
        beta_mother = beta_ind_true,
        beta_father = beta_ind_true,
        prop_MZ     = prop_MZ,
        mz_twin_cor = mz_twin_cor,
        dz_twin_cor = dz_twin_cor
      )
      
      out %>%
        mutate(
          target_R2_child    = R2_child,
          target_R2_indirect = R2_ind,
          beta_child_true    = beta_child_true,
          beta_ind_true      = beta_ind_true,
          prop_MZ            = prop_MZ,
          mz_twin_cor        = mz_twin_cor,
          dz_twin_cor        = dz_twin_cor
        )
    }))
  }))
}))

print(results_all)


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
    child_R2_lab = factor(target_R2_child,
                          levels = c(0.15, 0.10, 0.05, 0.025),
                          labels = c("Child R²=0.15", "Child R²=0.10", "Child R²=0.05", "Child R²=0.025")),
    ind_R2_lab = factor(target_R2_indirect,
                        levels = c(0.10, 0.05, 0.025),
                        labels = c("Indirect R²=0.10", "Indirect R²=0.05", "Indirect R²=0.025"))
  )

p_all <- ggplot(results_long, aes(x = n_families, y = power, color = effect, group = effect)) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(ind_R2_lab ~ child_R2_lab) +
  scale_x_continuous(breaks = n_families_grid) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  labs(
    title = "Power curves for Model 1 and Model 2 across child and indirect R² scenarios",
    subtitle = paste0("Dashed line = 80% power | prop_MZ=", prop_MZ,
                      ", MZ cor=", mz_twin_cor, ", DZ cor=", dz_twin_cor,
                      " | sims/point=", n_sims),
    x = "Number of families (each = 1 twin pair)",
    y = "Estimated power",
    color = "Effect"
  )

print(p_all)


indirect_only <- results_long %>%
  filter(effect %in% c("Model 2: Indirect (mother PGS)", "Model 2: Indirect (father PGS)"))

p_indirect <- ggplot(indirect_only, aes(x = n_families, y = power, color = effect, group = effect)) +
  geom_hline(yintercept = 0.80, linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_grid(ind_R2_lab ~ child_R2_lab) +
  scale_x_continuous(breaks = n_families_grid) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  labs(
    title = "Power for indirect effects (mother = father) in Model 2",
    x = "Number of families",
    y = "Estimated power",
    color = "Effect"
  )

print(p_indirect)


#note: In these simulations, you’re setting the data-generating β’s via R² using the Option 3 mapping (residual SD = 1). Because Model 2 includes correlated predictors (child/mother/father PGS), the observed partial R² in the fitted model won’t necessarily equal the target R² exactly; but this is still the right way to do scenario-based power checks.