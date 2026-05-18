rm(list = ls())

library(data.table)
library(R2jags)
library(zoo)

args <- commandArgs(trailingOnly = TRUE)
tt <- unlist(strsplit(args[1], split = "_"))[1]
vp <- unlist(strsplit(args[1], split = "_"))[2]
cat(sprintf("Ecoregion: %s | Cell ID: %s\n", tt, vp))

tt <- '2.2.2'
vp <- '619929'

# ── Load GDN input data ───────────────────────────────────────────────────────
path <- paste('/projectnb/modislc/users/seamorez/HLS_Pheno/GDN_climate_sensitivity/jags_input/', tt, sep='')
sstr <- paste(tt, '*', vp, '.rda', sep='')
file <- list.files(path, pattern = glob2rx(sstr), full.names = TRUE)
load(file)

# ── Gap-fill over full 320-DOY blocks ────────────────────────────────────────
gap_fill_block <- function(x) {
  x <- na.approx(x, na.rm = FALSE)
  x <- na.locf(x, na.rm = FALSE)
  x <- na.locf(x, fromLast = TRUE, na.rm = FALSE)
  return(x)
}

cat("NAs before gap-filling:", sum(is.na(dat$Tavg)), "\n")

ndoys_model       <- length(unique(dat$doy))   # 171 for DOY 151-320
n_blocks_raw <- nrow(dat) / ndoys_model        # 171 DOYs per block
for (i in 1:n_blocks_raw) {
  idx <- (ndoys_model*(i-1)+1):(ndoys_model*i)
  dat$Tavg[idx] <- gap_fill_block(dat$Tavg[idx])
  dat$SWin[idx] <- gap_fill_block(dat$SWin[idx])
  dat$Pprd[idx] <- gap_fill_block(dat$Pprd[idx])
}

cat("NAs after gap-filling:", sum(is.na(dat$Tavg)), "\n")

# ── Block structure and GUP onset ─────────────────────────────────────────────
dat <- as.data.table(dat)
setorder(dat, px_location, year, doy)
dat[, block := .GRP, by = .(px_location, year)]
n_blocks <- max(dat$block)

# Derive GUP DOY from full data BEFORE trimming — GUP occurs pre-solstice
# gup_doy pre-computed in input script — read directly from constant column
gup_onset <- unique(dat[, .(block, gup_doy)])
setorder(gup_onset, block)
GUP_sc <- as.vector(scale(gup_onset$gup_doy))

# ── Trim to post-solstice accumulation window (DOY 151–320) ──────────────────
# Pre-solstice T is collinear with GUP timing. GUP predictor absorbs those
# pre-season thermal effects plus genetic/GSL constraints on senescence.
# Post-solstice, all climate drivers have unambiguous senescence directionality.
gdn_start   <- 151L
dat_model   <- dat[doy >= gdn_start]
n_per_block <- 320L - gdn_start + 1L   # 171 DOYs per block

stopifnot(nrow(dat_model) == n_blocks * n_per_block)

# ── Design matrix — identical structure to GUP model ─────────────────────────
# Expected signs (all negative): cooler T, lower SW, shorter Pprd → more forcing
X <- as.matrix(data.table(
  intercept = 1,
  Tavg      = as.vector(scale(dat_model$Tavg)),
  SWin      = as.vector(scale(dat_model$SWin)),
  Pprd      = as.vector(scale(dat_model$Pprd))
))

Y <- dat_model$PS   # GDN binary: 0 before onset, 1 from onset DOY onwards

# ── JAGS data ─────────────────────────────────────────────────────────────────
data <- list(
  X          = X,
  Y          = Y,
  GUP        = GUP_sc,
  block      = dat_model$block,
  head_nodes = which(dat_model$doy == gdn_start),
  main_nodes = which(dat_model$doy != gdn_start),
  hmax       = 100,
  lambda     = 1
)
data$n  <- nrow(X)
data$np <- ncol(X)   # 4: intercept, Tavg, SWin, Pprd

# ── JAGS model ────────────────────────────────────────────────────────────────
# Structurally identical to GUP model — latent state h accumulates post-solstice
# forcing. GUP enters logit link only, testing legacy effect independent of
# post-solstice climate. Expected: beta_gup < 0 (earlier GUP → earlier GDN).
FitModel_GDN <- "
model {

  for (i in head_nodes) { h[i] <- 0 }

  for (i in main_nodes) {
    h[i] <- h[i-1] + max(0, X[i-1,] %*% beta) * (1 - h[i-1] / hmax)
  }

  for (i in 1:n) {
    Y[i]  ~ dbern(p[i])
    yp[i] ~ dbern(p[i])
    logit(p[i]) <- kappa + lambda * h[i] + beta_gup * GUP[block[i]]
  }

  kappa    ~ dnorm(0, 0.001) T(, 0)
  beta_gup ~ dnorm(0, 0.001)

  for (j in 1:np) {
    beta[j] ~ dnorm(0, 0.001)
  }
}
"

model <- jags(
  model.file         = textConnection(FitModel_GDN),
  data               = data,
  parameters.to.save = c('beta', 'beta_gup', 'kappa', 'yp', 'h'),
  n.chains  = 1,
  n.iter    = 10000,
  n.burnin  = 3000
)

gibbs <- model$BUGSoutput$sims.list

# ── Save ──────────────────────────────────────────────────────────────────────
outDir <- paste('/projectnb/modislc/users/seamorez/HLS_Pheno/GDN_climate_sensitivity/jags_output/', tt, sep='')
if (!dir.exists(outDir)) dir.create(outDir, recursive = TRUE)
setwd(outDir)
save(dat, dat_model, data, gibbs,
     file = paste('gdn_', tt, '_', vp, '.rda', sep=''))

# =============================================================================
# DIAGNOSTICS
# =============================================================================
library(ggplot2)

pred_names <- colnames(X)   # intercept, Tavg, SWin, Pprd
idx_tavg <- 2; idx_swin <- 3; idx_pprd <- 4

# ── 1. Beta posteriors ────────────────────────────────────────────────────────
beta_summary <- data.frame(
  predictor = c(pred_names[-1], "GUP (legacy)"),
  median    = c(apply(gibbs$beta[, -1], 2, median), median(gibbs$beta_gup)),
  lo        = c(apply(gibbs$beta[, -1], 2, quantile, 0.025), quantile(gibbs$beta_gup, 0.025)),
  hi        = c(apply(gibbs$beta[, -1], 2, quantile, 0.975), quantile(gibbs$beta_gup, 0.975))
)

ggplot(beta_summary, aes(x = predictor, y = median)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "GDN climate sensitivity (95% credible intervals)",
       subtitle = "Tavg, SWin, Pprd expected negative | GUP expected negative",
       x = NULL, y = "Beta (posterior median + 95% CI)") +
  theme_bw()

# ── 2. h trajectory and P(green-down) ────────────────────────────────────────
h_med      <- apply(gibbs$h, 2, median)
gup_by_obs <- GUP_sc[dat_model$block]
p_med      <- 1 / (1 + exp(-(median(gibbs$kappa) + data$lambda * h_med +
                               median(gibbs$beta_gup) * gup_by_obs)))
yp_med     <- apply(gibbs$yp, 2, median)

compare <- data.table(
  px_locs = dat_model$px_location,
  year    = dat_model$year,
  doy     = dat_model$doy,
  y       = data$Y,
  yp      = yp_med,
  p       = p_med,
  h       = h_med
)

h_doy <- compare[, .(h_med = median(h), p_med = median(p)), by = doy]

par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
plot(h_doy$doy, h_doy$h_med, type = "l", lwd = 2,
     xlab = "DOY", ylab = "Senescence state h",
     main = "Accumulated senescence state (avg across blocks)")
abline(h = data$hmax, lty = 2, col = "grey50")

plot(h_doy$doy, h_doy$p_med, type = "l", lwd = 2,
     xlab = "DOY", ylab = "P(green-down)",
     main = "Probability of green-down over season", ylim = c(0, 1))
abline(h = 0.5, lty = 2, col = "red")

# ── 3. Predicted vs. observed GDN DOY ────────────────────────────────────────
onset <- compare[, .(
  onset_obs  = min(c(Inf, doy[y  == 1]), na.rm = TRUE),
  onset_pred = min(c(Inf, doy[yp == 1]), na.rm = TRUE)
), by = .(px_locs, year)]

onset_clean <- onset[is.finite(onset_obs) & is.finite(onset_pred)]
rmse   <- sqrt(mean((onset_clean$onset_obs - onset_clean$onset_pred)^2))
ss_res <- sum((onset_clean$onset_obs - onset_clean$onset_pred)^2)
ss_tot <- sum((onset_clean$onset_obs - mean(onset_clean$onset_obs))^2)
r2     <- 1 - ss_res / ss_tot

ggplot(onset_clean, aes(x = onset_obs, y = onset_pred, color = factor(year))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = min(onset_clean$onset_obs), y = max(onset_clean$onset_pred),
           label = paste0("RMSE = ", round(rmse, 1), " days\nR² = ", round(r2, 3)),
           hjust = 0, fontface = "bold") +
  labs(title = "Predicted vs. observed green-down DOY",
       x = "Observed GDN (DOY)", y = "Predicted GDN (DOY)", color = "Year") +
  theme_bw()

cat(sprintf("R² (GDN DOY): %.3f\n", r2))
cat(sprintf("RMSE (days):  %.1f\n", rmse))

# ── 4. Single-pixel probability trajectories ──────────────────────────────────
px_focus <- unique(dat_model$px_location)[1]
traj <- compare[px_locs == px_focus]

ggplot(traj, aes(x = doy, y = p, group = year, color = factor(year))) +
  geom_line(alpha = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(title = paste("Green-down probability trajectories — pixel", px_focus),
       x = "DOY", y = "P(green-down)", color = "Year") +
  theme_bw()