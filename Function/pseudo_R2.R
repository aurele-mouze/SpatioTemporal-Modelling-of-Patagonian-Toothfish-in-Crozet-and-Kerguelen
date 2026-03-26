pseudo_R2 <- function(obs, pred_full, pred_reduced) {
  1 - var(obs - pred_full) / var(obs - pred_reduced)
}