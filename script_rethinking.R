library(rethinking)
library(tidyverse)

# CHAPTER 3 ---------------------------------------------------------------

# Easy
p_grid <- seq(from = 0, to = 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(6, size = 9, prob = p_grid)
posterior <- likelihood*prior
posterior <- posterior/sum(posterior)

set.seed(100)
samples <- sample(p_grid, 1e4, replace = T, prob = posterior)

sum(samples < 0.2)/1e4
sum(samples > 0.8)/1e4
sum(samples > 0.2 & samples < 0.8)/1e4
quantile(samples, 0.2)
quantile(samples, 0.8)
HPDI(samples, 0.66)
PI(samples, 0.66)

# Medium
p_grid <- seq(from = 0, to = 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(8, size = 15, prob = p_grid)
posterior <- likelihood*prior
posterior <- posterior/sum(posterior)
set.seed(100)
samples <- sample(p_grid, 1e4, replace = T, prob = posterior)

HPDI(samples, 0.9)
w <- rbinom(1e4, size = 15, prob = samples)
sum(w == 8)/1e4
v <- rbinom(1e4, size = 9, prob = samples)
sum(v == 6)/1e4

prior_h <- ifelse(p_grid < 0.5, 0, 2)
posterior_h <- likelihood*prior_h
posterior_h <- posterior_h/sum(posterior_h)
set.seed(100)
samples_h <- sample(p_grid, 1e4, replace = T, prob = posterior_h)

HPDI(samples_h, 0.9)
w_h <- rbinom(1e4, size = 15, prob = samples_h)
sum(w_h == 8)/1e4
v_h <- rbinom(1e4, size = 9, prob = samples_h)
sum(v_h == 6)/1e4

dbinom(6, size = 9, prob = 0.7)
dbinom(8, size = 15, prob = 0.7)

# Hard
data(homeworkch3)
# 100 families with 2 children
Nm <- sum(birth1) + sum(birth2) # n of males
Nf <- sum(1 - birth1) + sum(1 - birth2) # n of females

p_grid <- seq(from = 0, to = 1, length.out = 1000)
prior <- rep(1, 1000)
likelihood <- dbinom(Nm, size = 200, prob = p_grid)
posterior <- likelihood*prior
posterior <- posterior/sum(posterior)
set.seed(100)
samples <- sample(p_grid, 1e4, replace = T, prob = posterior)

p_grid[which.max(posterior)]
# chainmode(samples)
HPDI(samples, 0.5); HPDI(samples, 0.89); HPDI(samples, 0.97)

# compare 111 boys in 200 births
w <- rbinom(1e4, size = 200, prob = samples)
simplehist(w)
chainmode(w)

# compare 51 boys in 100 births (first borns)
v <- rbinom(1e4, size = 100, prob = samples)
simplehist(v)
chainmode(v)

tab_born <- data.frame(first = birth1, second = birth2)
tab_born %>% filter(birth1 == 0) %>% table() # 39 males out of 49 second borns
z <- rbinom(1e4, size = 49, prob = samples)
simplehist(z)
chainmode(z)


# CHAPTER 4 ---------------------------------------------------------------

pos <- replicate(1000, sum(runif(16, -1, 1)))
hist(pos)
plot(density(pos))

growth <- replicate(1e4, prod(1 + runif(12, 0, 0.1)))
dens(growth, norm.comp = T)

# Howell data
data("Howell1")
d <- Howell1
d2 <- d[d$age >= 18,]

# prior predictive check
curve(dnorm(x, 178, 20), from = 100, to = 250)
curve(dunif(x, 0, 50), from = -10, to = 60)

sample_mu <- rnorm(1e4, 178, 20)
sample_sigma <- runif(1e4, 0, 50)
prior_h <- rnorm(1e4, sample_mu, sample_sigma)
dens(prior_h)

# grid
mu.list <- seq(from = 150, to  = 160, length.out = 100)
sigma.list <- seq(from = 7, to = 9, length.out = 100)
post <- expand.grid(mu = mu.list, sigma = sigma.list)
post$LL <- sapply(1:nrow(post), function(i) {
  sum(dnorm(d2$height, post$mu[i], post$sigma[i], log = T))
})
post$prod <- post$LL + dnorm(post$mu, 178, 20, T) + dunif(post$sigma, 0, 50, T)
post$prob <- exp(post$prod - max(post$prod))

contour_xyz(post$mu, post$sigma, post$prob)
image_xyz(post$mu, post$sigma, post$prob)

# sampling from the posterior
sample.rows <- sample(1:nrow(post), size = 1e4, replace = T, prob = post$prob)
sample.mu <- post$mu[sample.rows]
sample.sigma <- post$sigma[sample.rows]
plot(sample.mu, sample.sigma, cex = 0.5, pch = 16, col = col.alpha(rangi2, 0.1))
dens(sample.mu)
PI(sample.mu)

d3 <- sample(d2$height, size = 20)

# quap
flist <- alist(
  height ~ dnorm(mu, sigma),
  mu ~ dnorm(178, 20),
  sigma ~ dunif(0, 50)
)

m4.1 <- quap(flist, data = d2)
precis(m4.1)

m4.2 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu ~ dnorm(178, 0.1),
    sigma ~ dunif(0, 50)
  ), data = d2)
precis(m4.2)

vcov(m4.1)
diag(vcov(m4.1))    # mu and sigma variances
cov2cor(vcov(m4.1)) # correlation matrix

post <- extract.samples(m4.1, n = 1e4)
plot(post, cex = 0.5, pch = 16, col = col.alpha(rangi2, 0.2))

# First linear model
library(rethinking)
data("Howell1")
d <- Howell1; d2 <- d[d$age >= 18,]

plot(d2$height ~ d2$weight)

b <- rlnorm(1e4, 0, 1)
dens(b, xlim = c(0,5), adj = 0.1)

xbar <- mean(d2$weight)
m4.3 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*(weight - xbar),
    a ~ dnorm(178, 20),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ), data = d2)


# EX. CHAPTER 4 -----------------------------------------------------------

# 4M1
sample_mu <- rnorm(1e4, 0, 10)
sample_sigma <- rexp(1e4, 1)
prior_sim <- rnorm(1e4, sample_mu, sample_sigma)
dens(prior_sim)
PI(prior_sim, 0.89)

# 4M2
flist <- alist(
  y ~ dnorm(mu, sigma),
  mu ~ dnorm(0, 10),
  sigma ~ dexp(1)
)

# 4M4
data_2018 <- data.frame(height = rnorm(30, 170, 5),
                        year = "I",
                        student = 1:30)
data_2019 <- data.frame(height = data_2018$height + sample(1:4, 30, replace = T),
                        year = "II",
                        student = 1:30)
data_2020 <- data.frame(height = data_2019$height + sample(1:4, 30, replace = T),
                        year = "III",
                        student = 1:30)
data_4M4 <- rbind(data_2018, data_2019, data_2020)

data_4M4 %>%
  ggplot(aes(year, height)) +
  geom_point() +
  geom_line(aes(group = student)) +
  theme_bw() +
  labs(y = "height (cm)")

flist <- alist(
  height ~ dnorm(mu, sigma),
  mu <- a + b*year,
  a ~ dnorm(170, 10),
  b ~ dlnorm(0, 1),
  sigma ~ dexp(1)
)

# 4M7
data("Howell1"); d <- Howell1; d2 <- d[d$age >= 18,]
xbar <- mean(d2$weight)

m4.3 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*(weight - xbar),
    a ~ dnorm(178, 20),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ), data = d2
)

m4.3_mod <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*weight,
    a ~ dnorm(178, 20),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ), data = d2
)

precis(m4.3)
precis(m4.3_mod)
round(vcov(m4.3), 3)
round(vcov(m4.3_mod), 3)
pairs(m4.3)
pairs(m4.3_mod)

post_4.3 <- extract.samples(m4.3, 1e4)
post_4.3_mod <- extract.samples(m4.3_mod, 1e4)

d2 %>%
  ggplot(aes(weight, height)) +
  geom_point(alpha = 0.6, size = 2, shape = 1, color = "blue") +
  geom_function(fun = function(x) mean(post_4.3$a) + mean(post_4.3$b)*(x - mean(x))) +
  geom_function(fun = function(x) mean(post_4.3_mod$a) + mean(post_4.3_mod$b)*x, color = "red", lty = 2) +
  # geom_abline(intercept = mean(post_4.3_mod$a), slope = mean(post_4.3_mod$b)) +
  theme_bw()

d2 %>%
  ggplot(aes(weight, height)) +
  geom_point(alpha = 0.6, size = 2, shape = 1, color = "blue") +
  geom_abline(aes(intercept = a - b*mean(d2$weight), slope = b), data = post_4.3[sample(1:1e4, 10),]) +
  geom_abline(aes(intercept = a, slope = b), data = post_4.3_mod[sample(1:1e4, 10),], color = "red", lty = 2) +
  theme_bw()

mu4.3 <- link(m4.3, n = 1e4)
mu4.3_mod <- link(m4.3_mod, n = 1e4)

data_ribbon_4.3 <- rbind(apply(link(m4.3, n = 1e4), 2, mean), apply(link(m4.3, n = 1e4), 2, PI, 0.89),
                         apply(link(m4.3_mod, n = 1e4), 2, mean), apply(link(m4.3_mod, n = 1e4), 2, PI, 0.89),
                         apply(sim(m4.3, n = 1e4), 2, PI, 0.67), apply(sim(m4.3_mod, n = 1e4), 2, PI, 0.67)) %>%
  t() %>% as.data.frame() %>%
  rename(mean = 1, PI5 = 2, PI94 = 3, mean_mod = 4, PI5_mod = 5, PI94_mod = 6,
         PI5_h_sim = 7, PI94_h_sim = 8, PI5_h_mod_sim = 9, PI94_h_mod_sim = 10)

d2 %>% bind_cols(data_ribbon_4.3) %>%
  ggplot(aes(weight, height)) +
  geom_ribbon(aes(ymin = PI5_h_sim, ymax = PI94_h_sim), alpha = 0.3, fill = "cyan") +
  geom_ribbon(aes(ymin = PI5, ymax = PI94), alpha = 0.5) +
  geom_point(alpha = 0.6, size = 2, shape = 1, color = "blue") +
  geom_line(aes(weight, mean), color = "red") +
  theme_bw() +
  labs(x = "weight (kg)", y = "height (cm)", title = "Posterior predictive simulations of m4.3")

d2 %>% bind_cols(data_ribbon_4.3) %>%
  ggplot(aes(weight, height)) +
  geom_ribbon(aes(ymin = PI5_h_mod_sim, ymax = PI94_h_mod_sim), alpha = 0.3, fill = "cyan") +
  geom_ribbon(aes(ymin = PI5_mod, ymax = PI94_mod), alpha = 0.5) +
  geom_point(alpha = 0.6, size = 2, shape = 1, color = "blue") +
  geom_line(aes(weight, mean_mod), color = "red") +
  theme_bw() +
  labs(x = "weight (kg)", y = "height (cm)", title = "Posterior predictive simulations of m4.3_mod")

# 4M8
data("cherry_blossoms"); d <- cherry_blossoms
precis(d)

ggplot(d) +
  geom_point(aes(year, doy), alpha = 0.5, color = "blue") +
  theme_bw() +
  labs(y = "day in year")

d2 <- d[complete.cases(d$doy),]
n_knots <- 15
knot_list <- quantile(d2$year, probs = seq(0, 1, length.out = n_knots))
# mgcv::place.knots(d2$year, n_knots)

B <- splines::bs(d2$year, knots = knot_list[-c(1,n_knots)], degree = 3, intercept = T)

cbind(d2, as.data.frame(B) %>% rename_with(~ str_c("basis", .))) %>%
  pivot_longer(starts_with("basis"), names_to = "BASIS_NAME", values_to = "BASIS_VALUE") %>%
  ggplot() +
  geom_line(aes(year, BASIS_VALUE, group = BASIS_NAME)) +
  theme_bw() +
  labs(y = "basis value")

m4.7 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + B %*% w,
    a ~ dnorm(100, 10),
    w ~ dnorm(0, 10),
    sigma ~ dexp(1)
  ), data = list(D = d2$doy, B = B), start = list(w = rep(0, ncol(B)))
)

post <- extract.samples(m4.7)
w <- apply(post$w, 2, mean)

cbind(d2, as.data.frame(B %*% diag(w)) %>% rename_with(~ str_c("basis", .))) %>%
  pivot_longer(starts_with("basis"), names_to = "BASIS_NAME", values_to = "BASIS_VALUE") %>%
  ggplot() +
  geom_line(aes(year, BASIS_VALUE, group = BASIS_NAME)) +
  theme_bw() +
  labs(y = "basis by weight value")

data_ribbon_4.7 <- rbind(apply(link(m4.7), 2, mean), apply(link(m4.7), 2, PI, 0.97)) %>%
  t() %>% as.data.frame() %>%
  rename(mean = 1, PI2 = 2, PI98 = 3)

d2 %>% bind_cols(data_ribbon_4.7) %>%
  ggplot(aes(year, doy)) +
  # geom_ribbon(aes(ymin = PI5_h_sim, ymax = PI94_h_sim), alpha = 0.3, fill = "cyan") +
  geom_ribbon(aes(ymin = PI2, ymax = PI98), alpha = 0.5) +
  geom_point(alpha = 0.6, size = 2, shape = 1, color = "blue") +
  geom_line(aes(year, mean), color = "red") +
  theme_bw() +
  labs(y = "day of year", title = "Posterior predictive simulations of m4.7")

# 4H1
data("Howell1"); d <- Howell1; d2 <- d[d$age >= 18,]

xbar <- mean(d2$weight)

m.4h1 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*(weight - xbar),
    a ~ dnorm(178, 20),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ), data = d
)

pred_data <- list(weight = c(46.95,43.72,64.78,32.59,54.63))
sim.height <- sim(m.4h1, n = 1e4, data = pred_data)
height.PI <- apply(sim.height, 2, PI, 0.89)

post <- extract.samples(m.4h1)
mu <- link(m.4h1, n = 1e4)
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI, 0.97)

# 4H2
data("Howell1"); d <- Howell1; d2 <- d[d$age < 18,]

xbar <- mean(d2$weight)

m.4h2 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*(weight - xbar),
    a ~ dnorm(110, 25),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ), data = d2
)

precis(m.4h2)

data_ribbon_4h2 <- rbind(apply(link(m.4h2), 2, mean), apply(link(m.4h2), 2, PI, 0.89),
                         apply(sim(m.4h2, n = 1e4), 2, PI, 0.89)) %>%
  t() %>% as.data.frame() %>%
  rename(mean = 1, mu_PI5 = 2, mu_PI94 = 3, h_PI5 = 4, h_PI94 = 5)

d2 %>% bind_cols(data_ribbon_4h2) %>%
  ggplot(aes(weight, height)) +
  geom_ribbon(aes(ymin = h_PI5, ymax = h_PI94), alpha = 0.3, fill = "cyan") +
  geom_ribbon(aes(ymin = mu_PI5, ymax = mu_PI94), alpha = 0.5, color = "lightgrey") +
  geom_point(alpha = 0.6, size = 2, color = "blue") +
  geom_line(aes(weight, mean), color = "red") +
  theme_bw() +
  labs(x = "weight (kg)", y = "height (cm)", title = "Posterior predictive simulations of m.4h2")

# m.4h2 overestimates heights in the low- and high-end of the weight range.
# Maybe it would be better to add a quadratic term (weight^2) to the linear model.

# 4H3
data("Howell1"); d <- Howell1

d$log_weight <- log(d$weight)
xbar <- mean(d$log_weight)

m4h3 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*(log_weight - xbar),
    a ~ dnorm(110, 25),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ), data = d
)

data_ribbon_4h3 <- rbind(apply(link(m4h3), 2, mean), apply(link(m4h3), 2, PI, 0.97),
                         apply(sim(m4h3, n = 1e4), 2, PI, 0.97)) %>%
  t() %>% as.data.frame() %>%
  rename(mean = 1, mu_PI2 = 2, mu_PI98 = 3, h_PI2 = 4, h_PI98 = 5)

d %>% bind_cols(data_ribbon_4h3) %>%
  ggplot(aes(weight, height)) +
  geom_ribbon(aes(ymin = h_PI2, ymax = h_PI98), alpha = 0.3, fill = "cyan") +
  geom_ribbon(aes(ymin = mu_PI2, ymax = mu_PI98), alpha = 0.5, color = "lightgrey") +
  geom_point(alpha = 0.6, size = 2, color = "blue") +
  geom_line(aes(weight, mean), color = "red") +
  theme_bw() +
  labs(x = "weight [kg] (log scale)", y = "height [cm]", title = "Posterior predictive simulations of m4h3")

# 4H4
data("Howell1"); d <- Howell1

d$weight_s <- standardize(d$weight)
d$weight_s2 <- d$weight_s^2

m4h4 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b1*weight_s + b2*weight_s2,
    a ~ dnorm(178, 20),
    b1 ~ dlnorm(0, 1),
    b2 ~ dlnorm(0, 1),
    sigma ~ dunif(0, 50)
  ), data = d
)

set.seed(13100)
prior <- extract.prior(m4h4, n = 1e4)
mu <- link(m4h4, post = prior,
           data = list(weight_s = seq(from = -2, to = 2, by = 0.1),
                       weight_s2 = seq(from = -2, to = 2, by = 0.1)^2)) %>%
  magrittr::set_rownames(str_c("sim", 1:nrow(.))) %>%
  t() %>% as.data.frame() %>%
  mutate(weight_s = seq(from = -2, to = 2, by = 0.1)) %>%
  pivot_longer(-weight_s, names_to = "nsim", values_to = "mu")

mu %>%
  filter(nsim %in% sample(unique(nsim), size = 20)) %>%
  ggplot(aes(weight_s, mu)) +
  geom_line(aes(group = nsim)) +
  theme_bw()

# a ~ dnorm(178, 20)
# b1 ~ dlnorm(0, 1) più spostato a destra o a sinistra a seconda di b2
# b2 ~ dlnorm(0, 1) solo negativo o solo positivo

# 4H5
data("cherry_blossoms"); d <- cherry_blossoms

d2 <- d %>% select(year, doy, temp) %>% filter(complete.cases(doy, temp))
d2$t <- standardize(d2$temp)

# linear
m4.5 <- quap(
  alist(
    doy ~ dnorm(mu, sigma),
    mu <- a + b*t,
    a ~ dnorm(75, 8),
    b ~ dnorm(0, 1),
    sigma ~ dunif(0, 50)
  ), data = d2
)

data_ribbon_4.5 <- rbind(apply(link(m4.5), 2, mean), apply(link(m4.5), 2, PI, 0.97),
                         apply(sim(m4.5, n = 1e4), 2, PI, 0.97)) %>%
  t() %>% as.data.frame() %>%
  rename(mean = 1, mu_PI2 = 2, mu_PI98 = 3, h_PI2 = 4, h_PI98 = 5)

d2 %>% bind_cols(data_ribbon_4.5) %>%
  ggplot(aes(temp, doy)) +
  geom_ribbon(aes(ymin = h_PI2, ymax = h_PI98), alpha = 0.3, fill = "cyan") +
  geom_ribbon(aes(ymin = mu_PI2, ymax = mu_PI98), alpha = 0.5, color = "lightgrey") +
  geom_point(alpha = 0.6, size = 2, color = "blue") +
  geom_line(aes(temp, mean), color = "red") +
  theme_bw() +
  labs(x = "temperature", y = "day of year", title = "Posterior predictive simulations of m4.5")


# CHAPTER 5 ---------------------------------------------------------------

library(tidyverse)
library(rethinking)

data("WaffleDivorce")
d <- WaffleDivorce

d$D <- standardize(d$Divorce)
d$M <- standardize(d$Marriage)
d$A <- standardize(d$MedianAgeMarriage)

# divorce rate versus median age at marriage
m5.1 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bA*A,
    a ~ dnorm(0, 0.2),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)

set.seed(10)
prior <- extract.prior(m5.1)
mu <- link(m5.1, post = prior, data = list(A = c(-2,2)))

plot(NULL, xlim = c(-2,2), ylim = c(-2,2))
for (i in 1:50) lines(c(-2,2), mu[i,], col = col.alpha("black", 0.4))

A_seq <- seq(from = -3, to = 3.2, length.out = 30)
mu <- link(m5.1, data = list(A = A_seq))
mu.mean <- apply(mu , 2, mean)
mu.PI <- apply(mu, 2, PI)

# plot it all
plot( D ~ A , data=d , col=rangi2 )
lines( A_seq , mu.mean , lwd=2 )
shade( mu.PI , A_seq )

# divorce rate versus marriage rate
m5.2 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM*M,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)

# marriage rate versus median age at marriage
m5.3 <- quap(
  alist(
    M ~ dnorm(mu, sigma),
    mu <- a + bA*A,
    a ~ dnorm(0, 0.2),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)

# DAGs
library(dagitty)
dag5.1 <- dagitty("dag{A -> M; A -> D; M -> D}")
coordinates(dag5.1) <- list(x = c(A = 0, D = 1, M = 2), y = c(A = 0, D = 1, M = 0))
drawdag(dag5.1)

DMA_dag1 <- dagitty('dag{ D <- A -> M -> D }')
impliedConditionalIndependencies( DMA_dag1 )
DMA_dag2 <- dagitty('dag{ D <- A -> M }')
impliedConditionalIndependencies( DMA_dag2 )

# correlations
d %>% select(D, A, M) %>% cor()

m5.4 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bM*M + bA*A,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)

plot(coeftab(m5.1, m5.2, m5.4), par = c("bM", "bA"))

# Simulating the divorce example
N <- 50
age <- rnorm(N)       # sim A
mar <- rnorm(N, -age) # sim A -> M
div <- rnorm(N, -age)  # sim A -> D
div2 <- rnorm(N, age + mar)

# predictor residual plots
mu <- link(m5.3)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)
mu.resid <- d$M - mu.mean

# posterior prediction plots
mu <- link(m5.4)
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)

D_sim <- sim(m5.3, n = 1e4)
D_PI <- apply(D_sim, 2, PI)

plot( mu_mean ~ d$D , col=rangi2 , ylim=range(mu_PI) , xlab="Observed divorce" , ylab="Predicted divorce" )
abline( a=0 , b=1 , lty=2 )
for ( i in 1:nrow(d) ) lines( rep(d$D[i],2) , mu_PI[,i] , col=rangi2 )
identify( x=d$D , y=mu_mean , labels=d$Loc )

# simulating spurious association
N <- 100 # number of cases
x_real <- rnorm( N ) # x_real as Gaussian with mean 0 and stddev 1
x_spur <- rnorm( N , x_real ) # x_spur as Gaussian with mean=x_real
y <- rnorm( N , x_real ) # y as Gaussian with mean=x_real
d <- data.frame(y,x_real,x_spur) # bind all together in data frame

m5.5 <- quap(
  alist(
    y ~ dnorm(mu, sigma),
    mu <- a + bR*x_real + bS*x_spur,
    a ~ dnorm(0, 0.2),
    bR ~ dnorm(0, 0.5),
    bS ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)

pairs(d)
precis(d)

# counterfactual plots
data(WaffleDivorce)
d <- list()
d$A <- standardize( WaffleDivorce$MedianAgeMarriage )
d$D <- standardize( WaffleDivorce$Divorce )
d$M <- standardize( WaffleDivorce$Marriage )

m5.3_A <- quap(
  alist(
    ## A -> D <- M
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM*M + bA*A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 ),
    ## A -> M
    M ~ dnorm( mu_M , sigma_M ),
    mu_M <- aM + bAM*A,
    aM ~ dnorm( 0 , 0.2 ),
    bAM ~ dnorm( 0 , 0.5 ),
    sigma_M ~ dexp( 1 )
  ) , data = d )

A_seq <- seq( from=-2 , to=2 , length.out=30 )
sim_dat <- data.frame( A=A_seq )
s <- sim( m5.3_A , data=sim_dat , vars=c("M","D") )

graphics::plot( sim_dat$A , colMeans(s$D) , ylim=c(-2,2) , type="l" ,
                xlab="manipulated A" , ylab="counterfactual D" )
shade( apply(s$D,2,PI) , sim_dat$A )
mtext( "Total counterfactual effect of A on D" )
