rm(list = ls())

library(data.table)
library(R2jags)

## C6 one tile multi-year
args <- commandArgs()
print(args)

# Load the input rda

tt <- '2.2.2'
vp <- '668891'

tt <- unlist(strsplit(args[3], split = "_"))[1] 
vp <- unlist(strsplit(args[3], split = "_"))[2] 

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
ma <- function(x, n = 7){stats::filter(x, rep(1 / n, n), sides = 2)}
stavg <- rep(NA, nrow(dat))
mint <- 7

for(i in 1:100){
  ta <- dat$Tavg[(250*(i-1)+1):(250*i)]
  mta <- ma(ta, mint)
  mta[1:(mint%/%2)] <- mta[(mint%/%2+1)]
  mta[(250+1-mint%/%2):250] <- mta[(250-mint%/%2)]
  stavg[(250*(i-1)+1):(250*i)] <- mta
}

print(sum(is.na(stavg)))

dat$Tavg <- stavg
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
