library(INLA)
library(INLAspacetime)
library(raster)
library(ncdf4)
library(evd)
require(sf)

# -------------------------------------------------
# 1. Data upload and preparation

path_in <- ## where you store the data (already uploaded here in GitHub)
path_out <- ## where you want to put fitted parameters and plots if you save them

infile <- paste0(path_in,"AnnualMaximaSeries_box.nc")
r <- stack(infile)

#Basic statistics on the raw input
r_median <- calc(r, fun=median, na.rm=TRUE)
plot(r_median)

#Dimensions of the dataset
n_years <- dim(r)[3]
years <- 1986:2020
n_lon <- dim(r)[2]
n_lat <- dim(r)[1]

#Put netcdf values into a table
mat_val <- rasterToPoints(r)
n_pts <- dim(mat_val)[1]

#Standardized series
q95 <- quantile(mat_val, 0.95)
q05 <- quantile(mat_val, 0.05)

# Calculate the difference between the 95th and 5th quantiles
q_diff <- as.numeric(q95 - q05)
mat_std <- as.matrix(mat_val/q_diff)
colnames(mat_std) <- c("lon","lat",years)


# -------------------------------------------------
# 2. INLA parameters setting

# Hyperparameters prior distribution definitions
spread = 1    
tail = 0.1  

# Those probabilities define q_alpha (location) ans s_beta (spread) parameters
# as the p.alpha quantile and distance between quantiles -p.beta/2 and p.beta/2
p.alpha= 0.5    
p.beta= 0.8       

# Set hyperparameters
hyper.spread = list(initial = 1,
                    fixed=FALSE,
                    prior = "loggamma",
                    param = c(3, 3))

tail.interval = c(0, 0.5)
##! remember to run the function "map.tail" before this 
##(find it at the end of the code)

tail.intern = map.tail(tail, tail.interval, inverse=TRUE)
hyper.tail = list(initial = tail.intern,
                  prior = "pc.gevtail",
                  param = c(7, tail.interval),
                  fixed= FALSE)

# Bind hyperparameter together
hyper.bgev = list(spread = hyper.spread, tail = hyper.tail)

# Define covariates (null)
null.matrix = matrix(nrow = n_years, ncol= 0)
spread.x = null.matrix
tail.x = null.matrix

# Set control parameters 
control.bgev = list(q.location = p.alpha,
                    q.spread = p.beta,
                    # quantile levels for the mixing part
                    q.mix= c(0.05, 0.20),
                    # the Beta(s1, s2) mixing distribution parameters.
                    # Hard-coded for the moment: s1=s2=5
                    beta.ab = 5)

# Define fitting formula without covariates
formula = inla.mdata(y, spread.x, tail.x) ~ -1 + intercept


# -------------------------------------------------
# 3. Fitting the distribution (single point)

# choose and index from 1 to n_pts (625)
i <- 100
time_series_standardized <- mat_std[i,3:(n_years+2)]  
time_series_lon <- mat_val[i,1]  
time_series_lat <- mat_val[i,2]

# plot the series and its location
plot(years,time_series_standardized,type='l')
writeLines(sprintf("lon: %.2f, lat: %.2f",time_series_lon,time_series_lat))

# Prepare data for the fit
data.bgev = data.frame(y = time_series_standardized, intercept=1, spread.x = spread.x, tail.x = tail.x)
  
# INLA fit
fit1 = inla(formula,
            family = "bgev",
            data = data.bgev,
            control.family = list(hyper = hyper.bgev,
                                  control.bgev = control.bgev),
            control.predictor = list(compute = TRUE),
            control.fixed = list(prec=1/var(time_series_standardized), mean=median(time_series_standardized)),
            control.compute = list(cpo = TRUE, dic= TRUE, waic = TRUE),
            control.inla = list(int.strategy = "eb"),
            verbose=FALSE, safe=TRUE)

# Summary of the results
fit1$summary.fixed
fit1$summary.hyperpar

# Extract posterior distributions mean for each parameter
location_mean = fit1$summary.fixed[1,1]
spread_mean = fit1$summary.hyperpar[1,1]
tail_mean = fit1$summary.hyperpar[2,1]

# Traditional GEV fit

fit2 <- fgev(time_series_standardized)
location_GEV <- fit2$estimate[1]
scale_GEV <- fit2$estimate[2]
shape_GEV <- fit2$estimate[3]

# 5. Fit assessment

# Convert into traditional GEV parameter to compare with other works
# the inputs should be: giveme.gev.par(q, sbeta, alpha, beta, xi)
gev.par <- giveme.gev.par(q = location_mean, sbeta = spread_mean, alpha = p.alpha, beta = p.beta,
                          xi = tail_mean)

# Compare with parameters from traditional GEV fit

writeLines(sprintf("Location: %.3f (bGEV with INLA), %.3f (GEV)",gev.par$mu,location_GEV))
writeLines(sprintf("Scale: %.3f (bGEV with INLA), %.3f (GEV)",gev.par$sigma, scale_GEV))
writeLines(sprintf("Shape: %.3f (bGEV with INLA), %.3f (GEV)",gev.par$xi,shape_GEV))

# Quality control plots for standardized series
# Generate the GEV distribution
x_obs <- seq(min(time_series_standardized),max(time_series_standardized),length.out =10)
x_gev <- seq(min(time_series_standardized),max(time_series_standardized),length.out =1000)
y_bgev <- dgev(x_gev, loc=gev.par$mu, scale=gev.par$sigma, shape=gev.par$xi) 
y_gev <- dgev(x_gev, loc=location_GEV, scale=scale_GEV, shape=shape_GEV) 

# Plot density function over observed normalized data density barplot
hist(time_series_standardized,breaks = x_obs ,prob = TRUE,#ylim=c(0,3),
     main="density plot",xlab = "level", ylab = "probability")
lines(x_gev,y_bgev,col='red',lw=3)
lines(x_gev,y_gev,col='blue',lw=3)
legend(0.8, 2.5, legend=c("bGEV (INLA)", "GEV"),
       col=c("red", "blue"), lw=3, cex=0.6)

# PP-Plot 
hist_obj <- hist(time_series_standardized,breaks = x_obs,plot = FALSE)
p_obs <- hist_obj$density/sum(hist_obj$density)
p_bmod <- dgev(hist_obj$mids, loc=gev.par$mu, scale=gev.par$sigma, shape=gev.par$xi) 
p_bmod <- p_bmod/sum(p_bmod)
p_mod <- dgev(hist_obj$mids, loc=location_GEV, scale=scale_GEV, shape=shape_GEV) 
p_mod <- p_mod/sum(p_mod)

plot(p_obs,p_bmod,type='p',pch=19,col='red',xlim=c(0,max(c(p_obs,p_mod))),ylim=c(0,max(c(p_obs,p_mod))),
     main="probability plot",xlab = "observed probability (frequency)", ylab = "fitted probability")
points(p_obs,p_mod,type='p',pch=19,col='blue')
lines(c(0,max(c(p_obs,p_mod))),c(0,max(c(p_obs,p_mod))),type='l',lty='dashed')
legend(0.15, 0.05, legend=c("bGEV (INLA)", "GEV"),
       col=c("red", "blue"),pch=19)

# QQ-Plot
q_obs <- quantile(time_series_standardized,seq(0.1,0.9,0.1))
q_bmod <- qgev(seq(0.1,0.9,0.1), loc=gev.par$mu, scale=gev.par$sigma, shape=gev.par$xi)
q_mod <- qgev(seq(0.1,0.9,0.1),  loc=location_GEV, scale=scale_GEV, shape=shape_GEV) 

plot(q_bmod,q_obs,col='red',type='p',pch=19,xlim=c(0,max(c(q_obs,q_mod))),ylim=c(0,max(c(q_obs,q_mod))),
     main="qq plot",xlab = "fitted quantile", ylab = "observed quantile")
points(q_mod,q_obs,type='p',pch=19,col='blue')
lines(c(0,max(c(q_obs,q_mod))),c(0,max(c(q_obs,q_mod))),type='l',lty='dashed')
legend(0.4, 0.1, legend=c("bGEV (INLA)", "GEV"),
       col=c("red", "blue"),pch=19)


