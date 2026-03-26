# 4) Fonction utilitaire pour prédire par lots (évite OOM)
predict_in_chunks <- function(model, newdata, chunk_size = 20000) {
  n <- nrow(newdata)
  i_start <- seq(1, n, by = chunk_size)
  out_list <- vector("list", length(i_start))
  for (k in seq_along(i_start)) {
    i1 <- i_start[k]
    i2 <- min(n, i1 + chunk_size - 1)
    message("Predict chunk ", k, " : rows ", i1, " to ", i2)
    out_list[[k]] <- predict(model, newdata = newdata[i1:i2, , drop = FALSE])
  }
  bind_rows(out_list)
}