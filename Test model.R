# Dossier contenant les fichiers
folder <- "C:/Users/Aurele/Documents/etude modèles crozet/Results/RepeatedCV/scores"

files <- list.files(folder, pattern = "\\.rds$", full.names = TRUE)

# Lecture de tous les fichiers
all_results <- map_dfr(files, readRDS)

# Vérification
glimpse(all_results)

all_results <- all_results %>%
  filter(!mesh %in% c("ScanWeightKg", "cutoff_11"))

summary_results_Crozet <- all_results %>%
  group_by(mesh, family) %>%
  summarise(
    mean_rmse = mean(rmse, na.rm = TRUE),
    mean_mae = mean(mae, na.rm = TRUE),
    mean_sum_loglik = mean(sum_loglik, na.rm = TRUE),
    mean_n = mean(n, na.rm = TRUE),
    mean_n_folds_tot = mean(n_folds_tot, na.rm = TRUE),
    mean_n_converged = mean(n_converged, na.rm = TRUE),
    mean_percent_converged = mean(percent_converged, na.rm = TRUE),
    n_repeats = n(),
    .groups = "drop"
  )

diff_results <- summary_results_Crozet %>%
  select(mesh, family, mean_rmse, mean_mae) %>%
  group_by(mesh) %>%
  mutate(
    rmse_ref = mean_rmse[family == "gamma"],
    mae_ref  = mean_mae[family == "gamma"],
    rmse_diff_pct = 100 * (mean_rmse - rmse_ref) / rmse_ref,
    mae_diff_pct  = 100 * (mean_mae - mae_ref) / mae_ref
  ) %>%
  ungroup()

diff_results

summary_results_Crozet <- all_results %>%
  select(mesh, family, rmse, mae) %>%
  pivot_longer(
    cols = c(rmse, mae),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(mesh, family, metric) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = sum(!is.na(value)),
    se = sd / sqrt(n),
    IC95 = qt(0.975, df = n - 1) * se,
    .groups = "drop"
  ) %>%
  mutate(
    model = paste(mesh, family, sep = " - "),
    metric = toupper(metric)
  )

ggplot(summary_results_Crozet,
       aes(x = reorder(model, mean), y = mean, colour = metric)) +
  geom_point(position = position_dodge(width = 0.5),
             size = 2.8) +
  geom_errorbar(aes(ymin = mean - IC95,
                    ymax = mean + IC95),
                width = 0.15,
                position = position_dodge(width = 0.5),
                linewidth = 0.6) +
  coord_flip() +
  labs(
    title = "Crozet",
    x = "Modèle",
    y = "Erreur moyenne",
    colour = "Métrique"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )




# Dossier contenant les fichiers
folder <- "C:/Users/Aurele/Documents/Results_Kerguelen/RepeatedCV/scores"

files <- list.files(folder, pattern = "\\.rds$", full.names = TRUE)

# Lecture de tous les fichiers
all_results <- map_dfr(files, readRDS)

# Vérification
glimpse(all_results)


summary_results_Kerguelen <- all_results %>%
  group_by(mesh, family) %>%
  summarise(
    mean_rmse = mean(rmse, na.rm = TRUE),
    mean_mae = mean(mae, na.rm = TRUE),
    mean_n = mean(n, na.rm = TRUE),
    mean_sum_loglik = mean(sum_loglik, na.rm = TRUE),
    mean_n_folds_tot = mean(n_folds_tot, na.rm = TRUE),
    mean_n_converged = mean(n_converged, na.rm = TRUE),
    mean_percent_converged = mean(percent_converged, na.rm = TRUE),
    n_repeats = n(),
    .groups = "drop"
  )



diff_results <- summary_results_Kerguelen %>%
  select(mesh, family, mean_rmse, mean_mae) %>%
  group_by(mesh) %>%
  mutate(
    rmse_ref = mean_rmse[family == "gamma"],
    mae_ref  = mean_mae[family == "gamma"],
    rmse_diff_pct = 100 * (mean_rmse - rmse_ref) / rmse_ref,
    mae_diff_pct  = 100 * (mean_mae - mae_ref) / mae_ref
  ) %>%
  ungroup()

diff_results


summary_results_Kerguelen <- all_results %>%
  select(mesh, family, rmse, mae) %>%
  pivot_longer(
    cols = c(rmse, mae),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(mesh, family, metric) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = sum(!is.na(value)),
    se = sd / sqrt(n),
    IC95 = qt(0.975, df = n - 1) * se,
    .groups = "drop"
  ) %>%
  mutate(
    model = paste(mesh, family, sep = " - "),
    metric = toupper(metric)
  )

ggplot(summary_results_Kerguelen,
       aes(x = reorder(model, mean), y = mean, colour = metric)) +
  geom_point(position = position_dodge(width = 0.5),
             size = 2.8) +
  geom_errorbar(aes(ymin = mean - IC95,
                    ymax = mean + IC95),
                width = 0.15,
                position = position_dodge(width = 0.5),
                linewidth = 0.6) +
  coord_flip() +
  labs(
    title = "Kerguelen",
    x = "Modèle",
    y = "Erreur moyenne",
    colour = "Métrique"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )
