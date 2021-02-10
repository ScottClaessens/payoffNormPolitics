# custom functions

# power analysis
runPowerAnalysis <- function(n = 450, nsims = 100) {
  # get previous data for simulation
  d <- 
    read.csv("data/previous/panizzaDG.csv") %>%
    left_join(read.csv("data/previous/panizzaKW.csv") %>% filter(right == 0), 
              by = c("part", "TL", "BL", "TR", "BR")) %>%
    left_join(read.csv("data/previous/panizzaKW.csv") %>% filter(right == 1), 
              by = c("part", "TL", "BL", "TR", "BR")) %>%
    drop_na() %>%
    transmute(
      part        = factor(part),                          # participant id
      type        = factor(paste0(TL, BL, TR, BR)),        # miniDG type
      non_selfish = 1 - selfish,                           # did participant choose unselfish option?
      self_adv    = self_adv,                              # payoffs to dictator - selfish option
      self_dis    = self_dis,                              # payoffs to dictator - non-selfish option
      diffPayoff  = (self_dis - self_adv) / 60,            # difference in payoffs to dictator
      norm_adv    = ifelse(TL > TR, rating.x, rating.y),   # normative rating - selfish option
      norm_dis    = ifelse(TL > TR, rating.y, rating.x),   # normative rating - non-selfish option
      diffNorm    = ordered(round(norm_dis - norm_adv, 2)) # difference in normative rating
    )
  # fit initial model to the complete dataset
  m1 <- brm(non_selfish ~ diffPayoff + mo(diffNorm) + (diffPayoff + mo(diffNorm) | part),
            data = d, family = bernoulli,
            prior = c(prior(normal(0, 1), class = Intercept),
                      prior(normal(0, 1), class = b),
                      prior(exponential(1), class = sd)), 
            iter = 3000, cores = 4, seed = 2113)
  # add random effects parameters (posterior means) to dataset
  getRanMean <- function(par) {
    m1 %>%
      spread_draws(r_part[part,parameter]) %>% 
      group_by(part, parameter) %>% 
      summarise(mean = mean(r_part)) %>% 
      filter(parameter == par) %>%
      mutate(part = factor(part)) %>%
      select(-parameter)
  }
  d <- left_join(d, getRanMean("diffPayoff") %>% rename(beta = mean), by = "part")
  d <- left_join(d, getRanMean("modiffNorm") %>% rename(gamma = mean), by = "part")
  # create function to randomly simulate vector with correlation coefficient
  # from https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
  complement <- function(y, rho, x) {
    if (missing(x)) x <- rnorm(length(y)) # Optional: supply a default if `x` is not given
    y.perp <- residuals(lm(x ~ y))
    rho * sd(y.perp) * y + y.perp * sd(y) * sqrt(1 - rho^2)
  }
  # create data simulation function
  simData <- function(seed) {
    # set seed
    set.seed(seed)
    # initialise empty data frame
    out <- data.frame()
    # sample n participants (repeats allowed)
    parts <- sample(1:166, n, replace = TRUE)
    # loop over these participants and add to "out" data.frame
    for (i in 1:n) out <- rbind(out, d %>% filter(part == parts[i]) %>% mutate(part = i))
    # add sdo and rwa, weakly positively correlated with beta and gamma parameters (based on previous research)
    beta <- out %>% group_by(part) %>% summarise(beta = mean(beta)) %>% pull(beta)
    gamma <- out %>% group_by(part) %>% summarise(gamma = mean(gamma)) %>% pull(gamma)
    out <- left_join(out, data.frame(part = 1:n, sdo = complement(beta, 0.24)), by = "part")
    out <- left_join(out, data.frame(part = 1:n, rwa = complement(gamma, 0.2)), by = "part")
    return(out)
  }
  # fit interaction model to first simulated dataset
  m2 <- brm(non_selfish ~ diffPayoff*sdo + mo(diffNorm)*rwa + (diffPayoff + mo(diffNorm) | part),
            data = simData(1), family = bernoulli,
            prior = c(prior(normal(0, 1), class = Intercept),
                      prior(normal(0, 1), class = b),
                      prior(exponential(1), class = sd)), 
            iter = 3000, cores = 4, seed = 2113)
  # summarise model fixed effects function
  summaryModel <- function(model) {summary(model)$fixed}
  # pass_through function for progress printing
  pass_through <- function(data, fun) {fun(data); data}
  # simulate data and fit model nsims times, extract parameters of interest
  s <-
    tibble(seed = 1:nsims) %>%
    mutate(d = map(seed, simData)) %>%
    mutate(pars = map2(d, seed, ~update(m2, newdata = .x, seed = .y, iter = 3000, cores = 4) %>% 
                         summaryModel() %>% 
                         pass_through(function(x) print(paste0("Simulation ", .y))))) %>%
    select(-d) %>% unnest(pars)
  return(s)
}

# dummy data for pre-registration
createDummyData <- function(n = 450, seed = 2113) {
  set.seed(seed)
  diffPayoff <- c(-60,-40,-36,-30,-9 ,-6 ,
                  -9 ,-4 ,-9 ,-30,-15,-6 ,
                  -3 ,-25,-1 ,-12,-21,-15)
  tibble(
    id      = rep(1:n, each = 18),
    dgType  = rep(paste0("dg", 1:18), times = n),
    nonself = rbinom(n*18, 1, 0.3),
    diffP   = rep(diffPayoff/60, times = n),
    diffN   = ordered(sample(-3:3, size = n*18, replace = T)),
    sdo     = rep(rnorm(n, 0, 1), each = 18),
    rwa     = rep(rnorm(n, 0, 1), each = 18),
    male    = rep(rbinom(n, 1, 0.5), each = 18),
    age     = rep(sample(18:70, size = n, replace = T), each = 18)
    )
}

# pre-registered model 0
fitModel0 <- function(d) {
  brm(nonself ~ 0 + Intercept,
      data = d, family = bernoulli,
      prior = c(prior(normal(0, 1), class = b)),
      cores = 4, seed = 2113)
}

# pre-registered model 1
fitModel1 <- function(d) {
  brm(nonself ~ 0 + Intercept + diffP + mo(diffN) + (1 + diffP + mo(diffN) | id),
      data = d, family = bernoulli,
      prior = c(prior(normal(0, 1), class = b),
                prior(exponential(1), class = sd)),
      cores = 4, seed = 2113)
}

# pre-registered model 2a (no controls)
fitModel2a <- function(d) {
  brm(nonself ~ 0 + Intercept + diffP*sdo + mo(diffN)*rwa + (1 + diffP + mo(diffN) | id),
      data = d, family = bernoulli,
      prior = c(prior(normal(0, 1), class = b),
                prior(exponential(1), class = sd)), 
      cores = 4, seed = 2113)
}

# pre-registered model 2b (with controls)
fitModel2b <- function(d) {
  brm(nonself ~ 0 + Intercept + diffP*sdo + mo(diffN)*rwa + male + age + (1 + diffP + mo(diffN) | id),
      data = d, family = bernoulli,
      prior = c(prior(normal(0, 1), class = b),
                prior(exponential(1), class = sd)), 
      cores = 4, seed = 2113)
}