rm(list = ls())

library(data.table)
library(R2jags)

## C6 one tile multi-year
args <- commandArgs(trailingOnly = TRUE)
# Load the input rda
tt <- unlist(strsplit(args[1], split = "_"))[1]
vp <- unlist(strsplit(args[1], split = "_"))[2]

cat(sprintf("Ecoregion: %s | Cell ID: %s\n", tt, vp))

# tt <- '2.2.2'
# vp <- '668891'

path <- paste('/projectnb/modislc/users/seamorez/HLS_Pheno/GUP_climate_sensitivity/jags_input/',tt,sep='')
sstr <- paste(tt,'*',vp,'.rda',sep='')
file <- list.files(path,pattern=glob2rx(sstr),full.names=T)
load(file)


############################################
library(zoo)

# Apply before the smoothing loop, after loading and ordering dat
gap_fill_block <- function(x){
  # linear interpolation for interior NAs
  x <- na.approx(x, na.rm = FALSE)
  # if NAs remain at edges (start/end of block), fill with nearest valid value
  x <- na.locf(x, na.rm = FALSE)           # forward fill
  x <- na.locf(x, fromLast = TRUE, na.rm = FALSE)  # backward fill for leading NAs
  return(x)
}

cat("NAs before gap-filling:", sum(is.na(dat$Tavg)), "\n")

for(i in 1:100){
  idx <- (250*(i-1)+1):(250*i)
  dat$Tavg[idx] <- gap_fill_block(dat$Tavg[idx])
  dat$SWin[idx] <- gap_fill_block(dat$SWin[idx])
  dat$Pprd[idx] <- gap_fill_block(dat$Pprd[idx])
}

cat("NAs after gap-filling:", sum(is.na(dat$Tavg)), "\n")
# should be 0 now
############################################


############################################
# Replace Tavg with its 7-day moving average
# ma <- function(x, n = 7){stats::filter(x, rep(1 / n, n), sides = 2)}
# stavg <- rep(NA, nrow(dat))
# mint <- 7
# 
# for(i in 1:100){
#   ta <- dat$Tavg[(250*(i-1)+1):(250*i)]
#   mta <- ma(ta, mint)
#   mta[1:(mint%/%2)] <- mta[(mint%/%2+1)]
#   mta[(250+1-mint%/%2):250] <- mta[(250-mint%/%2)]
#   stavg[(250*(i-1)+1):(250*i)] <- mta
# }
# 
# print(sum(is.na(stavg)))
# 
# dat$Tavg <- stavg
############################################
# Build design matrix X and binary outcome Y
# X <- as.matrix(data.table(
#   intercept = 1,
#   dat[,.(tmean = scale((tmin+tmax)/2),
#          dayl = scale(dayl),
#          cu = cu,
#          ittd = scale((tmin+tmax)/2)*scale(dayl),
#          ittc = scale((tmin+tmax)/2)*cu,
#          itdc = scale(dayl)*cu
#          )])
# )
# X <- as.matrix(data.table(
#   intercept = 1,
#   dat[,.(cu = cu,
#          ittd = scale((tmin+tmax)/2)*scale(dayl)
#   )])
# )
# X <- as.matrix(data.table(
#   intercept = 1,
#   dat[,.(tmean = scale((tmin+tmax)/2),
#          cu = cu
#          )])
# )
X <- as.matrix(data.table(
  intercept = 1,
  dat[, .(Tavg = scale(Tavg),     # replaces tmean = scale((tmin+tmax)/2)
          Pprd = scale(Pprd),     # replaces dayl  = scale(dayl)
          SWin = scale(SWin)      # new predictor
  )])
)
# X <- as.matrix(data.table(
#   intercept = 1,
#   dat[,.(tmean = scale((tmin+tmax)/2),
#          dayl = scale(dayl),
#          cu = cu
#   )])
# )
# X <- as.matrix(data.table(
#   intercept = 1,
#   dat[,.(tmax = scale(tmax),
#          tmin = scale(tmin),
#          dayl = scale(dayl),
#          cu = cu
#          )])
# )

Y <- dat$PS

data <- list(X=X,Y=Y,
             head_nodes = which(dat$doy==1),
             main_nodes = which(dat$doy!=1))

data$n <- nrow(data$X); data$np <- ncol(data$X)
data$hmax <- 100; data$lambda <- 1
############################################
# JAGS modeling
FitModel = "
model{
  # development state for head nodes (i.e. start of the driving forces)
  for (i in head_nodes){
    h[i] <- 0
  }
  # development state for main nodes (i.e. all nodes except for the head nodes)
  for (i in main_nodes){
    h[i] <- h[i-1]+ max(0,X[i-1,]%*%beta)*(1-h[i-1]/hmax)
  }
  # observation level
  for(i in 1:n){
    Y[i] ~ dbern(p[i])                   # binary outcome
    yp[i] ~ dbern(p[i])                  # predictions to validate
    logit(p[i]) <- kappa + lambda * h[i] # logit link
  }
  kappa ~ dnorm(0, 0.001)T(, 0)
  
  # priors for tb and beta's
  for(i in 1:np){
    beta[i] ~ dnorm(0, 0.001) # vague prior for beta  
  }
}
"  

model <- jags(model.file = textConnection(FitModel),
              data = data,
              parameters.to.save = c('beta','yp','h','kappa'),
              n.chains = 1,
              n.iter = 10000,
              n.burnin=3000)

gibbs <- model$BUGSoutput$sims.list


outDir <- paste('/projectnb/modislc/users/seamorez/HLS_Pheno/GUP_climate_sensitivity/jags_output/',tt,sep='')
if (!dir.exists(outDir)) {dir.create(outDir)}

setwd(paste('/projectnb/modislc/users/seamorez/HLS_Pheno/GUP_climate_sensitivity/jags_output/',tt,sep=''))
save(dat,data,gibbs,file=paste('ma_',tt,'_',vp,'.rda',sep=''))


#############################################################################################################################################################################################
library(data.table)
library(ggplot2)

n_sims   <- nrow(gibbs$beta)
np       <- ncol(gibbs$beta)
pred_names <- colnames(X)   # "intercept", "Tavg", "Pprd", "SWin"

# ── 1. Median development state h across all MCMC samples ─────────────────
h_med  <- apply(gibbs$h, 2, median)           # length n (all obs)
h_lo   <- apply(gibbs$h, 2, quantile, 0.025)
h_hi   <- apply(gibbs$h, 2, quantile, 0.975)

# ── 2. Probability of green-up p  ─────────────────────────────────────────
# kappa and h are already sampled; compute p from logit link
kappa_med <- median(gibbs$kappa)
p_med     <- 1 / (1 + exp(-(kappa_med + data$lambda * h_med)))

# ── 3. Predicted onset (first DOY where yp == 1) per block ────────────────
yp_med <- apply(gibbs$yp, 2, median)   # 0/1 per observation (median vote)

compare <- data.table(
  px_locs = dat$px_location,
  year    = dat$year,
  doy     = dat$doy,
  y       = data$Y,
  yp      = yp_med,
  p       = p_med,
  h       = h_med
)

onset <- compare[, .(
  onset_obs  = min(c(Inf, doy[y  == 1]), na.rm = TRUE),
  onset_pred = min(c(Inf, doy[yp == 1]), na.rm = TRUE)
), by = .(px_locs, year)]

##
beta_df <- as.data.frame(gibbs$beta)
colnames(beta_df) <- pred_names

beta_long <- reshape(beta_df,
                     varying   = pred_names,
                     v.names   = "value",
                     timevar   = "predictor",
                     times     = pred_names,
                     direction = "long")

ggplot(beta_long, aes(x = value, fill = predictor)) +
  geom_density(alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~predictor, scales = "free") +
  labs(title = "Posterior distributions of climate sensitivity coefficients",
       x = "Beta value", y = "Density") +
  theme_bw() + theme(legend.position = "none")
##
beta_summary <- data.frame(
  predictor = pred_names,
  median    = apply(gibbs$beta, 2, median),
  lo        = apply(gibbs$beta, 2, quantile, 0.025),
  hi        = apply(gibbs$beta, 2, quantile, 0.975)
)

ggplot(beta_summary, aes(x = predictor, y = median)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title  = "Climate sensitivity: 95% credible intervals",
       x = NULL, y = "Beta (posterior median + 95% CI)") +
  theme_bw()
##
# Average across all blocks for a "typical year" view
h_doy <- compare[, .(h_med = median(h),
                     p_med = median(p)),
                 by = doy]

par(mfrow = c(2,1), mar = c(4,4,2,1))

# Development state
plot(h_doy$doy, h_doy$h_med, type = "l", lwd = 2,
     xlab = "DOY", ylab = "Development state h",
     main = "Accumulated development state (avg across blocks)")
abline(h = data$hmax, lty = 2, col = "grey50")

# Probability of green-up
plot(h_doy$doy, h_doy$p_med, type = "l", lwd = 2,
     xlab = "DOY", ylab = "P(green-up)",
     main = "Probability of green-up over season",
     ylim = c(0,1))
abline(h = 0.5, lty = 2, col = "red")
##
# Remove Inf (blocks where event never occurred in obs or pred)
onset_clean <- onset[is.finite(onset_obs) & is.finite(onset_pred)]

rmse <- sqrt(mean((onset_clean$onset_obs - onset_clean$onset_pred)^2))

ggplot(onset_clean, aes(x = onset_obs, y = onset_pred, color = factor(year))) +
  geom_point(size = 2, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  annotate("text", x = min(onset_clean$onset_obs),
           y = max(onset_clean$onset_pred),
           label = paste0("RMSE = ", round(rmse, 1), " days"),
           hjust = 0, fontface = "bold") +
  labs(title = "Predicted vs. observed green-up DOY",
       x = "Observed onset (DOY)", y = "Predicted onset (DOY)",
       color = "Year") +
  theme_bw()
##
# Pick a single px_locs to examine
px_focus <- unique(dat$px_location)[1]

traj <- compare[px_locs == px_focus]

ggplot(traj, aes(x = doy, y = p, group = year, color = factor(year))) +
  geom_line(alpha = 0.7) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(title = paste("Green-up probability trajectories — pixel", px_focus),
       x = "DOY", y = "P(green-up)", color = "Year") +
  theme_bw()

## Rsq
onset_clean <- onset[is.finite(onset_obs) & is.finite(onset_pred)]

ss_res <- sum((onset_clean$onset_obs - onset_clean$onset_pred)^2)
ss_tot <- sum((onset_clean$onset_obs - mean(onset_clean$onset_obs))^2)
r2_onset <- 1 - ss_res / ss_tot

cat(sprintf("R² (onset DOY):  %.3f\n", r2_onset))
cat(sprintf("RMSE (days):     %.1f\n", sqrt(ss_res / nrow(onset_clean))))

# compare <- data.table(year = dat$year,
#                       id = dat$id,
#                       yp = apply(model$BUGSoutput$sims.list$yp, 2, median),
#                       y = data$Y)
# 
# onset = compare[,.(onset_pred = min(c(Inf, which(yp==1))),
#                    onset = min(c(Inf, which(y==1)))), .(year, id)]
# 
# ## plot
# layout(matrix(c(1:8,rep(9,8)), nrow = 4))
# for(i in 1:ncol(gibbs$beta)){
#   hist(gibbs$beta[,i],main=colnames(X)[i],50)
# }
# #
# hh <- apply(gibbs$h,2,median)
# hh <- matrix(hh,250,100)
# hh <- apply(hh,1,median)
# medip <- median(onset$onset_pred[seq(2,200,2)],na.rm=T)
# 
# plot(hh,type='l',lwd=2)
# abline(v=medip+30,lwd=2,lty=5)
# # 1:1
# plot(onset$onset, onset$onset_pred,xlim=c(0,250),ylim=c(0,250),
#      col = rainbow(50)[as.factor(onset$id)],pch =19)
# rmse <- sqrt(mean(((onset$onset-onset$onset_pred)^2),na.rm=T))
# mtext(paste('RMSE = ', round(rmse,2)), line = -2, adj = .1, font = 2)
# abline(0,1)
