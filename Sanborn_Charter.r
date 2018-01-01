# Supplemental Information
# Bayesian Brains without Probabilities
# Adam N. Sanborn and Nick Chater
# University of Warwick, Gibbet Hill Rd., Coventry CV4 7AL, UK
# Correspondence: a.n.sanborn@warwick.ac.uk

# FIGURE 1A

# small probability of uniform distribution plus a narrow multivariate Gaussian
iso.density <- function(xy, n.params=n.params){
  sigma <- 0.01*diag(n.params)
  mu <- rep(0.75, n.params)
  if(n.params == 1){
    out <- 0.5 + 0.5*dnorm(xy, mean=mu, sd=sigma)
  }else{
    out <- 0.5 + 0.5*dmvnorm(xy, mean=mu, sigma=sigma)
  }
  if(is.vector(xy)){
    if(any(xy > 1) | any(xy < 0)){
      out <- 0
    }
  }else{
    out[apply(xy > 1, 1, any)] <- 0
    out[apply(xy < 0, 1, any)] <- 0
  }
  return(out)
}


# two univariate Gaussians, one very narrow 
bimodal.density <- function(x){
  sigma.1 <- 0.1
  sigma.2 <- 0.02
  mu.1 <- 0.2
  mu.2 <- 0.8
  out <- 0.000001 + 0.999999*((0.001/0.2)*dunif(x, min=0.1, max=0.3) + (0.999/0.2)*dunif(x, min=0.7, max=0.9))
  
  return(out)
}

# plot a bimodal distribution and the output of an MCMC sampler started in each mode
library(mvtnorm)

par(mfrow=c(1,1))

# plot the 1D marginalization
x <- seq(from=0, to=1, by=0.001)
x.density <- bimodal.density(x)
plot(x, log(x.density), type='l', axes=FALSE, xlab='', ylab='')
box()
#abline(v=0.5, lty=3)

# Metropolis-Hastings sampler for all dimensions
n.reps <- 10
n.samples <- 100000
start.states <- c(0.3,0.7)
samples <- c()
proposal.sigma <- 0.05
for(i in 1:length(start.states)){
  samples[[i]] <- rep(NaN,n.samples)
  state <- start.states[i]
  state.prob <- bimodal.density(state)
  for(j in 1:n.samples){
    proposal <- rnorm(1, mean=state, sd=proposal.sigma)
    proposal.prob <- bimodal.density(proposal)
    # M-H acceptance function
    if((proposal.prob > state.prob) | (runif(1) < (proposal.prob/state.prob))){
      state <- proposal
      state.prob <- proposal.prob
    }
    samples[[i]][j] <- state
  }
  samples[[i]] <- samples[[i]][(n.samples/2):n.samples]
}

hist(samples[[1]], 50, xlim=c(0,1), main='Sampler started in left mode', freq=F)
hist(samples[[2]], 50, xlim=c(0,1), main='Sampler started in right mode', freq=F)



#############################
# FIGURE 1B
# R Code for demonstrating autocorrelations of samplers for different problems using trace plots
library(rjags)
library(R2jags)
library(rstan)
library(gplots)

# values to use for different distributions
to.include <- 1:4
dist.names <- c('round','correlated','bimodal','close_bimodal')[to.include]
scatter.lims <- list(c(-1,5),c(-1,5),c(0,4),c(0,4))[to.include]

# sampling values
n.chains <- 1
n.iter <- 10000

jags.model <- c()

jags.model[1] <- paste('model{
    tau[1,1] <- 2
    tau[1,2] <- 0
    tau[2,1] <- 0
    tau[2,2] <- 2
    mu[1,1] <- 2
    mu[1,2] <- 2
    x ~ dmnorm(mu,tau)
    obs.tau[1,1] <- 1e-100
    obs.tau[1,2] <- 0
    obs.tau[2,1] <- 0
    obs.tau[2,2] <- 1e-100
    d ~ dmnorm(x,obs.tau)
}')

jags.model[2] <- paste('model{
    tau1 <- 2
    tau2 <- 100
    mu <- 2
    x[1] ~ dnorm(mu,tau1)
    x[2] ~ dnorm(x[1],tau2)
    obs.tau[1,1] <- 1e-100
    obs.tau[1,2] <- 0
    obs.tau[2,1] <- 0
    obs.tau[2,2] <- 1e-100
    d ~ dmnorm(x,obs.tau)
}')

jags.model[3] <- paste('model{
  tau[1,1] <- 10
  tau[1,2] <- 0
  tau[2,1] <- 0
  tau[2,2] <- 10
  mu[1,1] <- 1
  mu[1,2] <- 1
  mu[2,1] <- 3
  mu[2,2] <- 3
  w[1] <- 0.5
  w[2] <- 0.5
  z ~ dcat(w)
  x ~ dmnorm(mu[z,],tau)
  obs.tau[1,1] <- 1e-100
  obs.tau[1,2] <- 0
  obs.tau[2,1] <- 0
  obs.tau[2,2] <- 1e-100
  d ~ dmnorm(x,obs.tau)
}')

jags.model[4] <- paste('model{
                       tau[1,1] <- 10
                       tau[1,2] <- 0
                       tau[2,1] <- 0
                       tau[2,2] <- 10
                       mu[1,1] <- 1.25
                       mu[1,2] <- 1.25
                       mu[2,1] <- 2.75
                       mu[2,2] <- 2.75
                       w[1] <- 0.5
                       w[2] <- 0.5
                       z ~ dcat(w)
                       x ~ dmnorm(mu[z,],tau)
                       obs.tau[1,1] <- 1e-100
                       obs.tau[1,2] <- 0
                       obs.tau[2,1] <- 0
                       obs.tau[2,2] <- 1e-100
                       d ~ dmnorm(x,obs.tau)
                       }')

mh.density <- list()
mh.density[[1]] <- function(x){
  sigma <- solve(diag(c(2,2)))
  mu.1 <- c(2,2)
  out <- dmvnorm(x, mean=mu.1, sigma=sigma)
  return(out)
}

mh.density[[2]] <- function(x){
  sigma1 <- sqrt(1/2)
  sigma2 <- sqrt(1/100)
  mu <- 2
  out <- dnorm(x[1],mean=mu,sd=sigma1) * dnorm(x[2],mean=x[1],sd=sigma2)
  return(out)
}

mh.density[[3]] <- function(x){
  sigma <- solve(diag(c(10,10)))
  mu.1 <- c(1,1)
  mu.2 <- c(3,3)
  w <- 0.5
  out <- w*dmvnorm(x, mean=mu.1, sigma=sigma) + (1-w)*dmvnorm(x, mean=mu.2, sigma=sigma)
  return(out)
}

mh.density[[4]] <- function(x){
  sigma <- solve(diag(c(10,10)))
  mu.1 <- c(1.25,1.25)
  mu.2 <- c(2.75,2.75)
  w <- 0.5
  out <- w*dmvnorm(x, mean=mu.1, sigma=sigma) + (1-w)*dmvnorm(x, mean=mu.2, sigma=sigma)
  return(out)
}


# Metropolis-Hastings sampler for all dimensions
mh.sampler <- function(density.function,init,n.iter=1000,n.reps=1){
  samples <- c()
  #sigma <- ((2.38^2)/n.params)*(sigma+(1e-10)*diag(n.params))
  n.params <- length(init)
  sigma <- 0.1*diag(n.params) # starting value
  for(i in 1:n.reps){
    samples[[i]] <- matrix(NaN,nrow=n.iter,ncol=n.params)
    state <- init
    state.prob <- density.function(state)
    for(j in 1:n.iter){
      proposal <- rmvnorm(1, mean=state, sigma=sigma)
      proposal.prob <- density.function(proposal)
      # M-H acceptance function
      if((proposal.prob > state.prob) | (runif(1) < (proposal.prob/state.prob))){
        state <- proposal
        state.prob <- proposal.prob
      }
      samples[[i]][j,] <- state
    }
    samples[[i]] <- samples[[i]][(n.iter/2):n.iter,] # keeps last half of iterations
  }
  return(samples)
}

# plot the traces and scatterplots
plot.helper <- function(sample.matrix,file.prefix,file.suffix,scatter.lims){
  
  pdf(paste(file.prefix,'_trace_',file.suffix,'.pdf',sep=''))
  par(mfrow=c(2,1))
  plot(mcmc(sample.matrix),type='p',main='',density=T,pch='.')
  dev.off()
  
  png(paste(file.prefix,'_scatter_',file.suffix,'.png',sep=''))
  extended.sample.matrix <- rbind(sample.matrix, 
                                  c(scatter.lims[1], scatter.lims[1]),
                                  c(scatter.lims[1], scatter.lims[2]),
                                  c(scatter.lims[2], scatter.lims[1]),
                                  c(scatter.lims[2], scatter.lims[2]))
  smoothScatter(extended.sample.matrix[,1],extended.sample.matrix[,2],
                colramp=colorRampPalette(c("white", "black")), 
                nrpoints=0, bandwidth=0.025, xlab='', ylab='', 
                xlim=scatter.lims, ylim=scatter.lims, axes=FALSE, nbin=1000)
  dev.off()
}


# copied from http://stackoverflow.com/questions/16774928/removing-part-of-a-graphic-in-r
my.filled.contour <-
  function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1,
                                                         length.out = ncol(z)), z, xlim = range(x, finite = TRUE),
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE),
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors,
            col = color.palette(length(levels) - 1), plot.title, plot.axes,
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1,
            axes = TRUE, frame.plot = axes, ...)
  {
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0))
      stop("increasing 'x' and 'y' values expected")
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
    par(las = las)
    mar <- mar.orig
    mar[4L] <- mar[2L]
    mar[2L] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i",
                yaxs = "i")
    #    rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
    #    if (missing(key.axes)) {
    #        if (axes)
    #            axis(4)
    #    }
    #    else key.axes
    #    box()
    if (!missing(key.title))
      key.title
    mar <- mar.orig
    mar[4L] <- 1
    par(mar = mar)
    plot.new()
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1L || ncol(z) <= 1L)
      stop("no proper 'z' matrix specified")
    if (!is.double(z))
      storage.mode(z) <- "double"
    .filled.contour(as.double(x), as.double(y), z, as.double(levels),
                            col = col)
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(x, side = 1)
        Axis(y, side = 2)
      }
    }
    else plot.axes
    if (frame.plot)
      box()
    if (missing(plot.title))
      title(...)
    else plot.title
    invisible()
  }

# plot the traces and scatterplots
contour.plot.helper <- function(density.function,file.suffix,scatter.lims){

  
  par(mfrow=c(1,1))
  png(paste('contour_',file.suffix,'.png',sep=''))
  vals <- seq(from=scatter.lims[1],to=scatter.lims[2],length.out=201)
  heights <- matrix(0,nrow=length(vals),ncol=length(vals))
  for(i in 1:length(vals)){
    for(j in 1:length(vals)){
      heights[i,j] <- density.function(c(vals[i],vals[j]))
    }
  }
  my.filled.contour(x=vals,y=vals,z=heights, col=colorpanel(12, "white", "black"), 
                 nlevels=10, xlab='', ylab='', 
                xlim=scatter.lims, ylim=scatter.lims, axes=FALSE, frame.plot=TRUE)
  dev.off()
}



for(m.idx in 1:length(to.include)){
  writeLines(jags.model[to.include[m.idx]], con='jags_model.txt')
  na.data <- list(d=c(NA,NA))
  data <- list(d=c(2,2))
  
  # Plot the contours
  contour.plot.helper(mh.density[[to.include[m.idx]]],dist.names[m.idx],scatter.lims[[m.idx]])
  
  # Direct sampling (JAGS without any data samples directly from the prior)
  model.direct <- jags.model('jags_model.txt', na.data, inits = NULL, n.chains = n.chains)
  sample <- jags.samples(model.direct, 'x', n.iter)
  direct.sample.matrix <- matrix(sample$x, ncol=2, dimnames=list(NULL,c('x1','x2')), byrow=T)
  direct.sample.matrix <- direct.sample.matrix[(n.chains*n.iter/2+1):(n.chains*n.iter),]
  plot.helper(direct.sample.matrix,'direct',dist.names[m.idx],scatter.lims[[m.idx]])
  
  # JAGS sampling (adding non-informative data forces non-ancestral sampling)
  model.jags <- jags.model('jags_model.txt', data, inits = NULL, n.chains = n.chains)
  sample <- jags.samples(model.jags, 'x', n.iter)
  jags.sample.matrix <- matrix(sample$x, ncol=2, dimnames=list(NULL,c('x1','x2')), byrow=T)
  jags.sample.matrix <- jags.sample.matrix[(n.chains*n.iter/2+1):(n.chains*n.iter),]
  plot.helper(jags.sample.matrix,'jags',dist.names[m.idx],scatter.lims[[m.idx]])
  
  # Metropolis-Hastings sampling
  mh.sample.matrix <- mh.sampler(mh.density[[to.include[m.idx]]], c(0,0), n.iter=n.iter)[[1]]
  dimnames(mh.sample.matrix) <- list(NULL,c('x1','x2'))
  plot.helper(mh.sample.matrix,'mh',dist.names[m.idx],scatter.lims[[m.idx]])  
}


