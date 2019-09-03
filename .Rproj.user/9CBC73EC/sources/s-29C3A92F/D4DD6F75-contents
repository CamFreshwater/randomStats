## Modeling integers 
# Practice modeling with binomial, poisson and multinomial distributions
# using rethinking examples
# NOTE: 64 bit R necessary for rstan to function

install.packages(c("devtools", "mvtnorm", "loo", "coda"), dependencies = TRUE)

library(rethinking)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
library(ggplot2)

### BINOMIAL MODEL ###

## Model chimpanzees prosocial tendencies in experimental setting
data(chimpanzees)
d <- chimpanzees

d$treatment <- 1 + d$prosoc_left + 2*d$condition
xtabs( ~ treatment + prosoc_left + condition , d )

# prior trimmed data list 11.10
dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  treatment = as.integer(d$treatment) )
# particles in 11-dimensional space
m11.4 <- ulam( #replacement for map2stan
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + b[treatment] ,
    a[actor] ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 0.5 )
  ) ,
  data=dat_list)
precis( m11.4 , depth=2 )

#note link function is logit so inv_logit backtransformation necessary
post <- extract.samples(m11.4)
p_left <- inv_logit( post$a )
plot( precis( as.data.frame(p_left) ) , xlim=c(0,1) )

diffs <- list( db13 = post$b[,1] - post$b[,3],
               db24 = post$b[,2] - post$b[,4] )
plot( precis(diffs) )

d$side <- d$prosoc_left + 1 # right 1, left 2
d$cond <- d$condition + 1 # no partner 1, partner 2
d_aggregated <- aggregate(
  d$pulled_left ,
  list( treatment=d$treatment , actor=d$actor ,
        side=d$side , cond=d$cond ) ,
  sum )
colnames(d_aggregated)[5] <- "left_pulls"


### AGGREGATED BINOMIAL MODEL ###

## UCB admissions data
# Perhaps similar to aggregated multinomial data?
data(UCBadmit)
d <- UCBadmit
d$gid <- ifelse( d$applicant.gender=="male" , 1 , 2 )
d$dept_id <- rep(1:6,each=2) 

m11.8 <- quap(
  alist(
    admit ~ dbinom( applications , p ) ,
    logit(p) <- a[gid] + delta[dept_id] ,
    a[gid] ~ dnorm( 0 , 1.5 ) ,
    delta[dept_id] ~ dnorm( 0 , 1.5 )
  ) , data=d )
precis( m11.8 , depth=2 )

post <- extract.samples(m11.8)
#absolute contrast
diff_a <- post$a[,1] - post$a[,2] 
#relative constrast (i.e. outcome)
diff_p <- inv_logit(post$a[,1]) - inv_logit(post$a[,2])
precis( list( diff_a=diff_a , diff_p=diff_p ) )
