library(INLA)
library(INLAspacetime)
library(raster)
library(ncdf4)
library(evd)
require(sf)
library(viridis)
library(extRemes)

# -------------------------------------------------
# 1. Data upload and preparation

path_in <- ## where you store the file AnnualMaximaSeries_box.nc
path_out <- ## where you want to save fitting output and plots

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
coords <- mat_val[,1:2]

#Standardized series
q95 <- quantile(mat_val, 0.95)
q05 <- quantile(mat_val, 0.05)

# Calculate the difference between the 95th and 5th quantiles
q_diff <- as.numeric(q95 - q05)
mat_std <- as.matrix(mat_val/q_diff)
colnames(mat_std) <- c("lon","lat",years)

#global variance:
hist(var(mat_std[,3:(n_years+2)]))
var_global <- var(c(mat_val[,3:(n_years+2)]))
var_global_std <- var(c(mat_std[,3:(n_years+2)]))

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

#-----------------------------------------------------
# 3 Cycle over all the points and get parameters

# for each point store the mean and the difference between 
# 0.025quant 0.975quant of the fitted posteriod distribution

# for bGEV INLA method
bGEV_location_post <- cbind(rep(0,dim(mat_val)[1]),rep(0,dim(mat_val)[1]))
bGEV_spread_post <- cbind(rep(0,dim(mat_val)[1]),rep(0,dim(mat_val)[1]))
bGEV_tail_post <- cbind(rep(0,dim(mat_val)[1]),rep(0,dim(mat_val)[1]))

# for traditional GEV method
GEV_location_post <- cbind(rep(0,dim(mat_val)[1]),rep(0,dim(mat_val)[1]))
GEV_scale_post <- cbind(rep(0,dim(mat_val)[1]),rep(0,dim(mat_val)[1]))
GEV_shape_post <- cbind(rep(0,dim(mat_val)[1]),rep(0,dim(mat_val)[1]))

for(i in 1:n_pts){
  
  writeLines(sprintf("%d / %d time series elaborated",i,n_pts))

  # select data from 
  # - mat_val for original series
  # - mat_std for standardized series
  ts <- mat_val[i,3:(n_years+2)]

  # A: INLA bGEV fit
  data.bgev = data.frame(y = ts, intercept=1, spread.x = spread.x, tail.x = tail.x)
  
  fit1 = inla(formula,
              family = "bgev",
              data = data.bgev,
              control.family = list(hyper = hyper.bgev,
                                    control.bgev = control.bgev),
              control.predictor = list(compute = TRUE),
              control.fixed = list(prec=1/var_global, mean=median(ts)),
              control.compute = list(cpo = TRUE),
              control.inla = list(int.strategy = "eb"),
              verbose=FALSE, safe=TRUE)
  
  fitpar = c(fit1$summary.fixed$mean,fit1$summary.hyperpar$mean[1],fit1$summary.hyperpar$mean[2])
  
  bGEV_location_post[i,1] <- fit1$summary.fixed$mean
  bGEV_location_post[i,2] <- fit1$summary.fixed[,5]-fit1$summary.fixed[,3]
  bGEV_spread_post[i,1] <- fit1$summary.hyperpar$mean[1]
  bGEV_spread_post[i,2] <- fit1$summary.hyperpar[1,5]-fit1$summary.hyperpar[1,3]
  bGEV_tail_post[i,1] <- fit1$summary.hyperpar$mean[2]
  bGEV_tail_post[i,2] <- fit1$summary.hyperpar[2,5]-fit1$summary.hyperpar[2,3]
  
  # B: Common GEV fit
  fit2 <- fgev(ts)
    
  GEV_location_post[i,1] <- fit2$estimate[1]
  GEV_location_post[i,2] <- fit2$std.err[1]
  GEV_scale_post[i,1] <- fit2$estimate[2]
  GEV_scale_post[i,2] <- fit2$std.err[2]
  GEV_shape_post[i,1] <- fit2$estimate[3]
  GEV_shape_post[i,2] <- fit2$std.err[3]

}


# -------------------------------------------------
# 4. (optional) Save output parameters

bGEV_location_info <- cbind(coords,bGEV_location_post)
colnames(location_info) <- c("lon","lat","mean","025to975quant")
bGEV_spread_info <- cbind(coords,bGEV_spread_post)
colnames(spread_info) <- c("lon","lat","mean","025to975quant")
bGEV_tail_info <- cbind(coords,bGEV_tail_post)
colnames(tail_info) <- c("lon","lat","mean","025to975quant")
write.table(location_info, file = paste0(path_out,"bGEV_loc_pd.txt"), row.names = FALSE, append = FALSE, quote = TRUE, sep = " ")
write.table(spread_info, file = paste0(path_out,"bGEV_spread_pd.txt"), row.names = FALSE, append = FALSE, quote = TRUE, sep = " ")
write.table(tail_info, file = paste0(path_out,"bGEV_tail_pd.txt"), row.names = FALSE, append = FALSE, quote = TRUE, sep = " ")

GEV_location_info <- cbind(coords,GEV_location_post)
colnames(location_info) <- c("lon","lat","mean","025to975quant")
GEV_scale_info <- cbind(coords,GEV_scale_post)
colnames(scale_info) <- c("lon","lat","mean","025to975quant")
GEV_shape_info <- cbind(coords,GEV_shape_post)
colnames(shape_info) <- c("lon","lat","mean","025to975quant")
write.table(location_info, file = paste0(path_out,"GEV_loc.txt"), row.names = FALSE, append = FALSE, quote = TRUE, sep = " ")
write.table(spread_info, file = paste0(path_out,"GEV_scale.txt"), row.names = FALSE, append = FALSE, quote = TRUE, sep = " ")
write.table(tail_info, file = paste0(path_out,"GEV_shape.txt"), row.names = FALSE, append = FALSE, quote = TRUE, sep = " ")

# -------------------------------------------------
# 5. Compare results

# Convert bGEV into usual parameters
GEVpar = giveme.gev.par(q = bGEV_location_info[,3], sbeta = bGEV_spread_info[,3], alpha = p.alpha, beta = p.beta, xi = bGEV_tail_info[,3])

# Get average values of parameters from bGEV and GEV
bGEV_loc_mean = GEVpar$mu
bGEV_scale_mean = GEVpar$sigma
bGEV_shape_mean = GEVpar$xi

GEV_loc_mean = GEV_location_post[,1]
GEV_scale_mean = GEV_scale_post[,1]
GEV_shape_mean = GEV_shape_post[,1]

# put into rasters for easily plot
r_bGEV_loc_mean = rasterFromXYZ(as.data.frame(cbind(coords,bGEV_loc_mean)))
r_bGEV_scale_mean = rasterFromXYZ(as.data.frame(cbind(coords,bGEV_scale_mean)))
r_bGEV_shape_mean = rasterFromXYZ(as.data.frame(cbind(coords,bGEV_shape_mean)))
r_GEV_loc_mean = rasterFromXYZ(as.data.frame(cbind(coords,GEV_loc_mean)))
r_GEV_scale_mean = rasterFromXYZ(as.data.frame(cbind(coords,GEV_scale_mean)))
r_GEV_shape_mean = rasterFromXYZ(as.data.frame(cbind(coords,GEV_shape_mean)))

#location
par(mfrow=c(1,2),m=c(2,2,2,2))
plot(r_bGEV_loc_mean, main="bGEV pd mean", zlim=c(35,165), col=rainbow(100), axes=FALSE)
plot(r_GEV_loc_mean, main="GEV", zlim=c(35,165), col=rainbow(100), axes=FALSE)
# main title location
mtext("Location parameter", outer=TRUE, cex=1.5)

#scale
par(mfrow=c(1,2),m=c(2,2,2,2))
plot(r_bGEV_scale_mean, main="bGEV pd mean", zlim=c(0,75), col=rainbow(100), axes=FALSE)
plot(r_GEV_scale_mean, main="GEV", zlim=c(0,75), col=rainbow(100), axes=FALSE)
# main title scale
mtext("Scale parameter", outer=TRUE, cex=1.5)

#shape
par(mfrow=c(1,2),m=c(2,2,2,2))
plot(r_bGEV_shape_mean, main="bGEV pd mean", zlim=c(-0.7,0.7),col=rainbow(100), axes=FALSE)
plot(r_GEV_shape_mean, main="GEV",  zlim=c(-0.7,0.7), col=rainbow(100), axes=FALSE)
# main title shape
mtext("Shape parameter", outer=TRUE, cex=1.5)

#return levels (give return time as input)

year_rl <- 100
p    <- 1/year_rl
yp   <- -log(1-p)
bGEV_rl <- (bGEV_loc_mean - (bGEV_scale_mean/0.2)*(1-yp^(-bGEV_shape_mean)) )
GEV_rl <- (GEV_loc_mean - (GEV_scale_mean/0.2)*(1-yp^(-GEV_shape_mean)) )

r_bGEV_rl = rasterFromXYZ(as.data.frame(cbind(coords,bGEV_rl)))
r_GEV_rl = rasterFromXYZ(as.data.frame(cbind(coords,GEV_rl)))

#return levels
par(mfrow=c(1,2),m=c(2,2,2,2))
plot(r_bGEV_rl, main="bGEV pd mean", col=rev(viridis(100)), axes=FALSE)
plot(r_GEV_rl, main="GEV", col=rev(viridis(100)), axes=FALSE)
# main title return levels
mtext(paste0("Return levels for ",year_rl," years"), outer=TRUE, cex=1.5)

