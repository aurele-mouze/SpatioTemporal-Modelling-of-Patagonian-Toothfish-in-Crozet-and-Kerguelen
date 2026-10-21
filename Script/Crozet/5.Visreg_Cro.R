vars_std <- c(
  "TAAFseason",
  "prof_bathy",
  "slope",
  "HookLength",
  "NB_MEAN_CACHALOT",
  "NB_MEAN_ORQUE",
  "bottomT"
)

# Pour chaque variable standardisé on :
#   Fait une prédiction via visreg
#   Remet la variable sur l échelle réelle
#   crée le graphique de l effet

stats_std <- dat_pos %>%
  summarise(
    across(
      all_of(vars_std),
      list(
        mean = ~mean(.x, na.rm = TRUE),
        sd   = ~sd(.x, na.rm = TRUE)
      )
    )
  )


nb_ham_med <- median(dat_pos$nb_ham, na.rm = TRUE)

offset_med <- log(
  0.56 * 0.002 * (nb_ham_med - 1) + 0.246
)

#Hameçon
H <- visreg(
  m_1,
  "HookLength_sc",
  plot = FALSE
)

y_ref <- H$fit$visregFit[1]

H$fit$indice_relatif <-
  H$fit$visregFit / y_ref - 1

H$fit$indice_lwr <-
  H$fit$visregLwr / y_ref - 1

H$fit$indice_upr <-
  H$fit$visregUpr / y_ref - 1

H$fit$indice_pct <-
  100 * (H$fit$visregFit / y_ref - 1)

H$fit$indice_lwr_pct <-
  100 * (H$fit$visregLwr / y_ref - 1)

H$fit$indice_upr_pct <-
  100 * (H$fit$visregUpr / y_ref - 1)

H$fit$HookLength_sc_real <-
  H$fit$HookLength_sc *
  stats_std$HookLength_sd +
  stats_std$HookLength_mean


p_hook <- ggplot(
  H$fit,
  aes(x = HookLength_sc_real, y = indice_pct)
) +
  geom_ribbon(
    aes(
      ymin = indice_lwr_pct,
      ymax = indice_upr_pct
    ),
    alpha = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line(linewidth = 1) +
  geom_rug(
    data = dat_pos_std,
    aes(x = HookLength),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.2
  ) +
  labs(
    x = "Longueur d'hameçon",
    y = "Variation relative de la CPUE (%)"
  ) +
  theme_bw()


#bathymétrie
B <- visreg(
  m_1,
  "prof_bathy_sc",
  plot = FALSE,
  predict = list(offset = offset_med)
)

# Retour à l'échelle réelle
B$fit$prof_bathy_real <-
  B$fit$prof_bathy_sc *
  stats_std$prof_bathy_sd +
  stats_std$prof_bathy_mean

# Valeur de référence = première prédiction
y_ref <- B$fit$visregFit[1]

# Indice relatif (%)
B$fit$indice_pct <-
  100 * (B$fit$visregFit / y_ref - 1)

B$fit$indice_lwr_pct <-
  100 * (B$fit$visregLwr / y_ref - 1)

B$fit$indice_upr_pct <-
  100 * (B$fit$visregUpr / y_ref - 1)

# Graphique
p_bathy <- ggplot(
  B$fit,
  aes(x = prof_bathy_real, y = indice_pct)
) +
  geom_ribbon(
    aes(
      ymin = indice_lwr_pct,
      ymax = indice_upr_pct
    ),
    alpha = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line(linewidth = 1) +
  geom_rug(
    data = dat_pos_std,
    aes(x = prof_bathy),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.2
  ) +
  labs(
    x = "Bathymétrie (m)",
    y = "Variation relative de la CPUE (%)"
  ) +
  theme_bw()


#slope
S <- visreg(
  m_1,
  "slope_sc",
  plot = FALSE,
  predict = list(offset = offset_med)
)

# Retour à l'échelle réelle
S$fit$slope_real <-
  S$fit$slope_sc *
  stats_std$slope_sd +
  stats_std$slope_mean

# Valeur de référence = première prédiction
y_ref <- S$fit$visregFit[1]

# Indice relatif (%)
S$fit$indice_pct <-
  100 * (S$fit$visregFit / y_ref - 1)

S$fit$indice_lwr_pct <-
  100 * (S$fit$visregLwr / y_ref - 1)

S$fit$indice_upr_pct <-
  100 * (S$fit$visregUpr / y_ref - 1)

# Graphique
p_slope <- ggplot(
  S$fit,
  aes(x = slope_real, y = indice_pct)
) +
  geom_ribbon(
    aes(
      ymin = indice_lwr_pct,
      ymax = indice_upr_pct
    ),
    alpha = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line(linewidth = 1) +
  geom_rug(
    data = dat_pos_std,
    aes(x = slope),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.2
  ) +
  labs(
    x = "Pente",
    y = "Variation relative de la CPUE (%)"
  ) +
  theme_bw()


#Nombre orque
O <- visreg(
  m_1,
  "NB_MEAN_ORQUE_sc",
  plot = FALSE,
  predict = list(offset = offset_med)
)

# Retour à l'échelle réelle
O$fit$NB_MEAN_ORQUE_real <-
  O$fit$NB_MEAN_ORQUE_sc *
  stats_std$NB_MEAN_ORQUE_sd +
  stats_std$NB_MEAN_ORQUE_mean

# Valeur de référence = première prédiction
y_ref <- O$fit$visregFit[1]

# Indice relatif (%)
O$fit$indice_pct <-
  100 * (O$fit$visregFit / y_ref - 1)

O$fit$indice_lwr_pct <-
  100 * (O$fit$visregLwr / y_ref - 1)

O$fit$indice_upr_pct <-
  100 * (O$fit$visregUpr / y_ref - 1)

# Graphique
p_orque <- ggplot(
  O$fit,
  aes(x = NB_MEAN_ORQUE_real, y = indice_pct)
) +
  geom_ribbon(
    aes(
      ymin = indice_lwr_pct,
      ymax = indice_upr_pct
    ),
    alpha = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line(linewidth = 1) +
  geom_rug(
    data = dat_pos_std,
    aes(x = NB_MEAN_ORQUE),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.2
  ) +
  labs(
    x = "Nombre d'orques",
    y = "Variation relative de la CPUE (%)"
  ) +
  theme_bw()

#Cachalot
C <- visreg(
  m_1,
  "NB_MEAN_CACHALOT_sc",
  plot = FALSE,
  predict = list(offset = offset_med)
)

# Retour à l'échelle réelle
C$fit$NB_MEAN_CACHALOT_real <-
  C$fit$NB_MEAN_CACHALOT_sc *
  stats_std$NB_MEAN_CACHALOT_sd +
  stats_std$NB_MEAN_CACHALOT_mean

# Valeur de référence = première prédiction
y_ref <- C$fit$visregFit[1]

# Indice relatif (%)
C$fit$indice_pct <-
  100 * (C$fit$visregFit / y_ref - 1)

C$fit$indice_lwr_pct <-
  100 * (C$fit$visregLwr / y_ref - 1)

C$fit$indice_upr_pct <-
  100 * (C$fit$visregUpr / y_ref - 1)

# Graphique
p_cachalot <- ggplot(
  C$fit,
  aes(x = NB_MEAN_CACHALOT_real, y = indice_pct)
) +
  geom_ribbon(
    aes(
      ymin = indice_lwr_pct,
      ymax = indice_upr_pct
    ),
    alpha = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line(linewidth = 1) +
  geom_rug(
    data = dat_pos_std,
    aes(x = NB_MEAN_CACHALOT),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.2
  ) +
  labs(
    x = "Nombre de cachalots",
    y = "Variation relative de la CPUE (%)"
  ) +
  theme_bw()

#Année
A <- visreg(
  m_1,
  "TAAFseason_sc",
  plot = FALSE,
  predict = list(offset = offset_med)
)

# Retour à l'échelle réelle
A$fit$TAAFseason_real <-
  A$fit$TAAFseason_sc *
  stats_std$TAAFseason_sd +
  stats_std$TAAFseason_mean

# Valeur de référence = première prédiction
y_ref <- A$fit$visregFit[1]

# Indice relatif (%)
A$fit$indice_pct <-
  100 * (A$fit$visregFit / y_ref - 1)

A$fit$indice_lwr_pct <-
  100 * (A$fit$visregLwr / y_ref - 1)

A$fit$indice_upr_pct <-
  100 * (A$fit$visregUpr / y_ref - 1)

# Graphique
p_annee <- ggplot(
  A$fit,
  aes(x = TAAFseason_real, y = indice_pct)
) +
  geom_ribbon(
    aes(
      ymin = indice_lwr_pct,
      ymax = indice_upr_pct
    ),
    alpha = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line(linewidth = 1) +
  geom_rug(
    data = dat_pos_std,
    aes(x = TAAFseason),
    inherit.aes = FALSE,
    sides = "b",
    alpha = 0.2
  ) +
  labs(
    x = "Année",
    y = "Variation relative de la CPUE (%)"
  ) +
  theme_bw()


#Mois
M <- visreg(
  m_1,
  "Month",
  plot = FALSE,
  predict = list(offset = offset_med)
)

# Valeur de référence = première prédiction
y_ref <- M$fit$visregFit[1]

# Indice relatif (%)
M$fit$indice_pct <-
  100 * (M$fit$visregFit / y_ref - 1)

M$fit$indice_lwr_pct <-
  100 * (M$fit$visregLwr / y_ref - 1)

M$fit$indice_upr_pct <-
  100 * (M$fit$visregUpr / y_ref - 1)

# Graphique
p_month <- ggplot(
  M$fit,
  aes(x = Month, y = indice_pct)
) +
  geom_ribbon(
    aes(
      ymin = indice_lwr_pct,
      ymax = indice_upr_pct
    ),
    alpha = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line(linewidth = 1) +
  scale_x_continuous(
    breaks = 1:12
  ) +
  labs(
    x = "Mois",
    y = "Variation relative de la CPUE (%)"
  ) +
  theme_bw()

#Vessel
# Vessel
V <- visreg(
  m_1,
  "Vessel",
  plot = FALSE,
  predict = list(offset = offset_med)
)

# Modalité de référence = première modalité
y_ref <- V$fit$visregFit[1]

# Indice relatif (%)
V$fit$indice_pct <-
  100 * (V$fit$visregFit / y_ref - 1)

V$fit$indice_lwr_pct <-
  100 * (V$fit$visregLwr / y_ref - 1)

V$fit$indice_upr_pct <-
  100 * (V$fit$visregUpr / y_ref - 1)

# Anonymisation des navires
V$fit$Vessel_anonyme <- LETTERS[
  seq_along(unique(V$fit$Vessel))
]

# Conserver l'ordre des modalités
V$fit$Vessel_anonyme <- factor(
  V$fit$Vessel_anonyme,
  levels = unique(V$fit$Vessel_anonyme)
)

p_vessel <- ggplot(
  V$fit,
  aes(x = Vessel_anonyme, y = indice_pct)
) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = indice_lwr_pct,
      ymax = indice_upr_pct
    ),
    width = 0.2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Navire",
    y = "Variation relative de la CPUE (%)"
  ) +
  theme_bw()

#Cacholot fiabilité
C_obs <- visreg(
  m_1,
  "cachalot_obs",
  plot = FALSE,
  predict = list(offset = offset_med)
)

y_ref <- C_obs$fit$visregFit[1]

C_obs$fit$indice_pct <-
  100 * (C_obs$fit$visregFit / y_ref - 1)

C_obs$fit$indice_lwr_pct <-
  100 * (C_obs$fit$visregLwr / y_ref - 1)

C_obs$fit$indice_upr_pct <-
  100 * (C_obs$fit$visregUpr / y_ref - 1)

C_obs$fit$absence_obs <- factor(
  C_obs$fit$cachalot_obs,
  levels = c("Absent", "Présent"),
  labels = c("Oui", "Non")
)

p_cachalot_obs <- ggplot(
  C_obs$fit,
  aes(x = absence_obs, y = indice_pct)
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = indice_lwr_pct,
      ymax = indice_upr_pct
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Correction de l'effort d'observation des cachalots",
    y = "Variation relative de la CPUE (%)"
  ) +
  theme_bw()

#obs orque
O_obs <- visreg(
  m_1,
  "orque_obs",
  plot = FALSE,
  predict = list(offset = offset_med)
)

y_ref <- O_obs$fit$visregFit[1]

O_obs$fit$indice_pct <-
  100 * (O_obs$fit$visregFit / y_ref - 1)

O_obs$fit$indice_lwr_pct <-
  100 * (O_obs$fit$visregLwr / y_ref - 1)

O_obs$fit$indice_upr_pct <-
  100 * (O_obs$fit$visregUpr / y_ref - 1)

O_obs$fit$absence_obs <- factor(
  O_obs$fit$orque_obs,
  levels = c("Absent", "Présent"),
  labels = c("Oui", "Non")
)



p_orque_obs <- ggplot(
  O_obs$fit,
  aes(x = absence_obs, y = indice_pct)
) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = indice_lwr_pct,
      ymax = indice_upr_pct
    ),
    width = 0.15
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Correction de l'effort d'observation d'orques",
    y = "Variation relative de la CPUE (%)"
  ) +
  theme_bw()

#Cachalot x orque
CO <- visreg2d(
  m_1,
  xvar = "NB_MEAN_CACHALOT_sc",
  yvar = "NB_MEAN_ORQUE_sc",
  plot = FALSE
)

# ------------------------------------------------------------
# 1. Mettre la matrice des prédictions sous forme de data.frame
# ------------------------------------------------------------

df_CO <- expand.grid(
  NB_MEAN_CACHALOT_sc = CO$x,
  NB_MEAN_ORQUE_sc = CO$y
)

df_CO$visregFit <- as.vector(CO$z)


# ------------------------------------------------------------
# 2. Retour à l'échelle réelle
# ------------------------------------------------------------

df_CO$NB_MEAN_CACHALOT_real <-
  df_CO$NB_MEAN_CACHALOT_sc *
  stats_std$NB_MEAN_CACHALOT_sd +
  stats_std$NB_MEAN_CACHALOT_mean

df_CO$NB_MEAN_ORQUE_real <-
  df_CO$NB_MEAN_ORQUE_sc *
  stats_std$NB_MEAN_ORQUE_sd +
  stats_std$NB_MEAN_ORQUE_mean


# ------------------------------------------------------------
# 3. Indice relatif (%) par rapport à la première prédiction
# ------------------------------------------------------------

y_ref <- df_CO$visregFit[1]

df_CO$indice_pct <-
  100 * (df_CO$visregFit / y_ref - 1)


# ------------------------------------------------------------
# 4. Graphique
# ------------------------------------------------------------
min_pct <- min(df_CO$indice_pct, na.rm = TRUE)

p_cachalot_orque <- ggplot(
  df_CO,
  aes(
    x = NB_MEAN_CACHALOT_real,
    y = NB_MEAN_ORQUE_real
  )
) +
  geom_tile(
    aes(fill = indice_pct)
  ) +
  geom_contour(
    aes(z = indice_pct),
    breaks = 0,
    colour = "black",
    linewidth = 0.4
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "red",
    midpoint = min_pct/2,
    limits =  c(min_pct, 0),
    name = "Variation relative\nCPUE (%)"
  ) +
  labs(
    x = "Nombre de cachalots",
    y = "Nombre d'orques"
  ) +
  theme_bw()



#Création des groupement de figure en 2 catégories (variable de capturabilité ou de densité)

figure_effets_catch <- (
    p_cachalot +
    p_cachalot_obs +
    p_month +
    p_orque +
    p_orque_obs +
    p_hook +
    p_vessel +
    p_cachalot_orque
) +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A") &
  theme(
    axis.title.y = element_text(size = 8),
    axis.title.x = element_text(size = 9)
  )



figure_effets_dens <- (
  p_bathy +
    p_slope +
    p_annee )+
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A")

figure_effets_dens

#Sauvegarde

ggsave(
  here("data", output_plot, "figure_effets_catch_Cro.jpg"),
  figure_effets_catch,
  width = 12,
  height = 8,
  units = "in",
  dpi = 300
)

ggsave(
  here("data", output_plot, "figure_effets_dens_Cro.jpg"),
  figure_effets_dens,
  width = 12,
  height = 4.8,
  units = "in",
  dpi = 300
)
