library(tidyverse)
library(MASS)
library(mgcv)
library(cowplot)

if(!require(gpe)) devtools::install_github("goldingn/gpe")

### Picture of 3 GP realizations

set.seed(123)
k1 <- gpe::rbf("x")
suppressWarnings(gpe::demoKernel(k1, ndraw = 3))



### Fit a model to some data using process convolution appeoximation

### Plot the data
data(mcycle)
plt_mcycle <- ggplot(data=mcycle) + geom_point(aes(x=times, y=accel)) +
  ylab("") + xlab("") + ggtitle("The infamous motorcycle data")

### Discrete convolution
u <- seq(-10, 75, length=50)
sd_fac <- 1
# Basis matrix
B_gau <- outer(
  mcycle$times,
  u,
  FUN=function(x,y,sd){dnorm(x,y,sd)}, sd=sd_fac*diff(u[1:2])
)

### Fit the model using mgcv (I added a linear effect too.)
fit_gau <- mgcv::gam(
  accel ~ times + B_gau,
  data=mcycle,
  paraPen=list(B_gau=list(S=diag(length(u)))),
  method="REML"
)

### Make some predictions with the fitted model
pred_times <- seq(min(mcycle$times), 60, 0.5)
# Define the basis matrix for the predictions
B_pred <- outer(
  pred_times,
  u, 
  FUN=function(x,y,sd){dnorm(x,y,sd)}, sd=sd_fac*diff(u[1:2])
)

# Use mgcv to make the predictions
pred_data <- as.data.frame(
  list(
    times=pred_times,
    predict(fit_gau, newdata = list(times=pred_times, B_gau=B_pred), se.fit=T)
  )
)

### Make a nifty pic of the predictions
ggplot(data=mcycle) + geom_point(aes(x=times, y=accel)) +
  geom_path(aes(x=times, y=fit), data=pred_data) +
  geom_ribbon(
    aes(x=times, ymin=fit-2*se.fit, ymax=fit+2*se.fit), 
    data=pred_data, alpha=0.2
  ) + 
  geom_ribbon(
    aes(x=times, 
        ymin=fit-2*(se.fit+sqrt(fit_gau$sig2)), 
        ymax=fit+2*(se.fit+sqrt(fit_gau$sig2))
    ), 
    data=pred_data, alpha=0.1
  )

### Optimize the kernel width
sd_fac <- seq(1,5,0.1)
REML <- rep(NA, length(sd_fac))
for(i in seq_along(sd_fac)){
  sd <- sd_fac[i]*diff(u[1:2])
  B_gau <- outer(
    mcycle$times, 
    u, 
    FUN=function(x,y,sd){dnorm(x,y,sd)}, sd=sd
  )
  fit_gau <- mgcv::gam(
    accel ~ times + B_gau, 
    data=mcycle, 
    paraPen=list(B_gau=list(S=diag(length(u)))),
    method="REML"
  )
  REML[i] <- fit_gau$gcv.ubre
}

# plot the REML score
sd_fac_opt = sd_fac[which.min(REML)]
ggplot()+geom_path(aes(x=sd_fac, y=REML), lwd=2) + 
  geom_vline(xintercept = sd_fac_opt, color="blue", lwd=2)


### Refit model with optimum kernel width
sd <- sd_fac_opt*diff(u[1:2])
B_gau <- outer(
  mcycle$times, 
  u, 
  FUN=function(x,y,sd){dnorm(x,y,sd)}, sd=sd
)
fit_gau <- mgcv::gam(
  accel ~ times + B_gau, 
  data=mcycle, 
  paraPen=list(B_gau=list(S=diag(length(u)))),
  method="REML"
)

pred_times <- seq(min(mcycle$times), 60, 0.5)
B_pred = B_gau <- outer(
  pred_times,
  u, 
  FUN=function(x,y,sd){dnorm(x,y,sd)}, sd=sd
)

pred_data <- as.data.frame(
  list(
    times=pred_times,
    predict(fit_gau, newdata = list(times=pred_times, B_gau=B_pred), se.fit=T)
  )
)

plt_opt <- ggplot(data=mcycle) + geom_point(aes(x=times, y=accel)) +
  geom_path(aes(x=times, y=fit), data=pred_data) +
  geom_ribbon(
    aes(x=times, ymin=fit-2*se.fit, ymax=fit+2*se.fit), 
    data=pred_data, alpha=0.2
  ) + 
  geom_ribbon(
    aes(x=times, 
        ymin=fit-2*(se.fit+sqrt(fit_gau$sig2)), 
        ymax=fit+2*(se.fit+sqrt(fit_gau$sig2))
    ), 
    data=pred_data, alpha=0.1
  )
print(plt_opt)



