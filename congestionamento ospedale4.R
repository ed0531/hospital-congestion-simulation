# ============================================================
# PS simulation: Hawkes-like arrivals + triage + congestion
# ============================================================

# ── Parametri ─────────────────────────────────────────
lambda <- 7.5      # intensità base media
eta    <- 0.12     # auto-eccitazione Hawkes
beta   <- 1.2      # decadimento memoria

c      <- 25       # soglia congestione soft
C      <- 60       # soglia sovraffollamento hard
alpha  <- 0.08     # rallentamento dimissioni

T      <- 24
dt     <- 1/12     # 5 minuti
Nsim   <- 5000

steps <- round(T / dt)
times <- seq(0, T - dt, by = dt)

# ── Baseline giornaliera NHPP ─────────────────────────
lambda0_t <- function(t) {
  lambda * (1 + 0.5 * sin(2 * pi * (t - 10) / 24))
}
# ============================================================
# Calibrazione triage da psaccessi.csv
# ============================================================

read.csv(db))

clean_num <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub(",", ".", x)
  as.numeric(x)
}

data$`TOTALE ACCESSI` <- as.numeric(data$`TOTALE ACCESSI`)

triage_cols <- c(
  "TRIAGE BIANCO (%)",
  "TRIAGE VERDE (%)",
  "TRIAGE GIALLO (%)",
  "TRIAGE ROSSO (%)"
)

data[triage_cols] <- lapply(data[triage_cols], clean_num)

tot_accessi <- sum(data$`TOTALE ACCESSI`, na.rm = TRUE)

p_bianco <- sum(data$`TOTALE ACCESSI` * data$`TRIAGE BIANCO (%)` / 100, na.rm = TRUE) / tot_accessi
p_verde  <- sum(data$`TOTALE ACCESSI` * data$`TRIAGE VERDE (%)`  / 100, na.rm = TRUE) / tot_accessi
p_giallo <- sum(data$`TOTALE ACCESSI` * data$`TRIAGE GIALLO (%)` / 100, na.rm = TRUE) / tot_accessi
p_rosso  <- sum(data$`TOTALE ACCESSI` * data$`TRIAGE ROSSO (%)`  / 100, na.rm = TRUE) / tot_accessi

p_triage <- c(
  bianco = p_bianco,
  verde  = p_verde,
  giallo = p_giallo,
  rosso  = p_rosso
)

round(p_triage, 4)
sum(p_triage)









# ── Triage ────────────────────────────────────────────
p_triage <- c(
  bianco = 0.0543,
  verde  = 0.6884,
  giallo = 0.2304,
  rosso  = 0.0223
)

p_triage <- p_triage / sum(p_triage)

theta <- c(
  bianco = 1/1.2,
  verde  = 1/3,
  giallo = 1/5,
  rosso  = 1/8
)

# ── Matrici risultati ─────────────────────────────────
set.seed(42)

N_tot_mat <- matrix(0L, nrow = Nsim, ncol = steps)
N_B_mat   <- matrix(0L, nrow = Nsim, ncol = steps)
N_V_mat   <- matrix(0L, nrow = Nsim, ncol = steps)
N_G_mat   <- matrix(0L, nrow = Nsim, ncol = steps)
N_R_mat   <- matrix(0L, nrow = Nsim, ncol = steps)

lambda_mat <- matrix(NA_real_, nrow = Nsim, ncol = steps)

# ── Simulazione ───────────────────────────────────────
for (s in seq_len(Nsim)) {
  
  N_cls <- round(lambda * 3 * p_triage)
  names(N_cls) <- names(p_triage)
  
  H <- 0   # memoria Hawkes
  
  for (t in seq_len(steps)) {
    
    # intensità con auto-eccitazione
    lambda_eff <- lambda0_t(times[t]) + eta * H
    lambda_eff <- max(lambda_eff, 0)
    
    # arrivi totali
    arr_tot <- rpois(1, lambda_eff * dt)
    
    # aggiornamento memoria Hawkes
    H <- exp(-beta * dt) * H + arr_tot
    
    # assegnazione triage
    arr_cls <- as.vector(rmultinom(1, arr_tot, prob = p_triage))
    names(arr_cls) <- names(p_triage)
    
    # congestione sulle dimissioni
    N_tot <- sum(N_cls)
    slowdown <- 1 + alpha * pmax(0, N_tot - c) / c
    
    theta_eff <- theta / slowdown
    
    dep_cls <- rbinom(
      n    = length(N_cls),
      size = N_cls,
      prob = pmin(1, theta_eff * dt)
    )
    names(dep_cls) <- names(p_triage)
    
    # aggiornamento stato
    N_cls <- N_cls + arr_cls - dep_cls
    N_cls <- pmax(0L, N_cls)
    names(N_cls) <- names(p_triage)
    
    # salvataggio
    N_B_mat[s, t]   <- N_cls["bianco"]
    N_V_mat[s, t]   <- N_cls["verde"]
    N_G_mat[s, t]   <- N_cls["giallo"]
    N_R_mat[s, t]   <- N_cls["rosso"]
    N_tot_mat[s, t] <- sum(N_cls)
    
    lambda_mat[s, t] <- lambda_eff
  }
}

# ── Statistiche ───────────────────────────────────────
mean_N <- colMeans(N_tot_mat)
q05    <- apply(N_tot_mat, 2, quantile, 0.05)
q95    <- apply(N_tot_mat, 2, quantile, 0.95)
risk   <- colMeans(N_tot_mat > C)

mean_B <- colMeans(N_B_mat)
mean_V <- colMeans(N_V_mat)
mean_G <- colMeans(N_G_mat)
mean_R <- colMeans(N_R_mat)

mean_lambda <- colMeans(lambda_mat)

cat("Rischio massimo: ",
    round(max(risk) * 100, 1), "%\n")

cat("Occupazione media alle 14h: ",
    round(mean_N[14 * 12]), "paz.\n")

cat("Ore con rischio >=20%: ",
    round(mean(risk >= 0.20) * T, 1), "ore su 24\n")

# ── Plot ──────────────────────────────────────────────
par(mfrow = c(4, 1), mar = c(3, 4, 2, 1))

# 1. Intensità Hawkes
plot(times, mean_lambda, type = "l", lwd = 2,
     xlab = "", ylab = "Arrivi/ora",
     main = "Intensità arrivi con memoria Hawkes")

# 2. Occupazione totale
plot(times, mean_N, type = "l", col = "#3266ad", lwd = 2,
     ylim = c(0, max(q95, na.rm = TRUE) * 1.1),
     xlab = "", ylab = "Pazienti",
     main = "Occupazione simulata totale")

polygon(c(times, rev(times)), c(q95, rev(q05)),
        col = adjustcolor("#3266ad", 0.15), border = NA)

lines(times, mean_N, col = "#3266ad", lwd = 2)
abline(h = C, col = "#e24b4a", lty = 2, lwd = 1.5)

legend("topleft",
       c("Media", "IC 90%", "Soglia C"),
       col = c("#3266ad", adjustcolor("#3266ad", 0.4), "#e24b4a"),
       lty = c(1, 1, 2),
       lwd = c(2, 8, 1.5),
       bty = "n",
       cex = 0.85)

# 3. Composizione triage
plot(times, mean_B, type = "l", col = "grey40", lwd = 2,
     ylim = c(0, max(mean_B + mean_V + mean_G + mean_R, na.rm = TRUE) * 1.1),
     xlab = "", ylab = "Pazienti",
     main = "Occupazione media per triage")

lines(times, mean_V, col = "darkgreen", lwd = 2)
lines(times, mean_G, col = "orange", lwd = 2)
lines(times, mean_R, col = "red", lwd = 2)

legend("topleft",
       c("Bianco", "Verde", "Giallo", "Rosso"),
       col = c("grey40", "darkgreen", "orange", "red"),
       lty = 1,
       lwd = 2,
       bty = "n",
       cex = 0.85)

# 4. Rischio sovraffollamento
plot(times, risk * 100, type = "l", col = "#e24b4a", lwd = 2,
     ylim = c(0, 100),
     xlab = "Ora",
     ylab = "Probabilità (%)",
     main = "Rischio sovraffollamento P(N(t) > C)")

polygon(c(times, rev(times)),
        c(risk * 100, rep(0, steps)),
        col = adjustcolor("#e24b4a", 0.12),
        border = NA)

lines(times, risk * 100, col = "#e24b4a", lwd = 2)
abline(h = 20, col = "#ef9f27", lty = 2, lwd = 1.5)