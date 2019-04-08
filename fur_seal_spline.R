library(crawl)
library(sf)
library(ggplot2)
library(mgcv)
library(mvtnorm)
library(foreach)
library(ggspatial)
library(tidyverse)

data("northernFurSeal")
northernFurSeal <- st_as_sf(northernFurSeal, coords=c("long","lat"), crs=4326) 
northernFurSeal <- st_transform(northernFurSeal, 3832)

# pause for a pic
np = ptolemy::npac() # install with-- devtools::install_github("jmlondon/ptolemy")
plt_data <- ggplot2::ggplot() + layer_spatial(northernFurSeal, color="blue") + 
  annotation_spatial(np, fill=1)
###

northernFurSeal <- cbind(northernFurSeal, st_coordinates(northernFurSeal))
northernFurSeal <- st_drop_geometry(northernFurSeal)
northernFurSeal$hour <- as.numeric(northernFurSeal$GMT)/3600

## Data for each coordinate needs to be stacked
nfs_stack <- data.frame(hour=northernFurSeal$hour, loc=northernFurSeal$X, 
                        quality = northernFurSeal$loc_class, coord="X")
nfs_stack <- rbind(
  nfs_stack, 
  data.frame(hour=northernFurSeal$hour, loc=northernFurSeal$Y, 
             quality=northernFurSeal$loc_class, coord="Y")
  )
head(nfs_stack)

## define Argos error variance form
nfs_stack$ones <- 1 # necessary to trick the variance function
fix <- c('3'=log(150), '1'=log(500), '2'=log(250))
var_func <- varExp(
  fixed = fix,
  form = formula(~ones|quality)
)

fit = mgcv::gamm(
  loc ~ 0 + coord + s(hour,by=coord, k=200, bs=c("ps","ps")), 
  weights = var_func,
  data=nfs_stack, method="REML")

### Fitted path
newdata <- data.frame(
  hour = rep(seq(ceiling(min(nfs_stack$hour)), floor(max(nfs_stack$hour)), 1),2)
) 
newdata$coord <- c(rep("X",nrow(newdata)/2), rep("Y",nrow(newdata)/2))
pred <- matrix(predict(fit$gam, newdata), nrow=nrow(newdata)/2, ncol=2) %>% 
  `colnames<-`(c("X","Y")) %>% as_tibble()
plt_pred <- plt_data + geom_path(aes(x=X, y=Y), data=pred, color='red')


### Make some simulated paths
B = predict(fit[["gam"]],newdata=newdata, type="lpmatrix")
beta = coef(fit$gam)
V = vcov(fit$gam, unconditional = T)
beta_smp = t(rmvnorm(20, beta, V)) 
paths <- foreach(i=1:20)%do%{
  p <- B %*% beta_smp[,i]
  p <- matrix(p, nrow=nrow(newdata)/2, ncol=2)
  p
}

paths <- do.call(rbind, paths)
paths <- data.frame(rep = rep(1:20, each=nrow(newdata)/2), paths)
colnames(paths)[2:3] <- c("X","Y")

## Plot them up
bb = c(range(paths$X), range(paths$Y)) + c(-1,1,-1,1)*100000
plt_tracks <- ggplot(data=paths) + geom_path(aes(x=X, y=Y, group=rep), color="red",alpha=0.2) +
  xlab("Longitude") + ylab("Latitude") + geom_sf(data=np, col=1, fill=1) + 
  coord_sf(xlim=bb[1:2], y=bb[3:4])


