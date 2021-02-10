# drake plan
plan <- drake_plan(
  # power analysis with previous data
  power = runPowerAnalysis(),
  # preregistration
  preregistration = rmarkdown::render(knitr_in("preregistration.Rmd"), quiet = TRUE),
  # create dummy data
  d = createDummyData(),
  # fit dummy models
  m0  = fitModel0(d),
  m1  = fitModel1(d),
  m2a = fitModel2a(d),
  m2b = fitModel2b(d),
  # loo
  loo0  = loo(m0),
  loo1  = loo(m1),
  loo2a = loo(m2a),
  loo2b = loo(m2b),
  # summaries
  s0  = summary(m0),
  s1  = summary(m1),
  s2a = summary(m2a),
  s2b = summary(m2b),
  # loo comparisons
  compare1 = loo_compare(loo0, loo1),
  compare2 = loo_compare(loo1, loo2a)
)