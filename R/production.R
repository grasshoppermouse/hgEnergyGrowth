
inv_logit <- plogis

# source(here("R/extract-digitized-plot-data.R"))
# source(here("R", "koster2020.R"))

# factor to multiply production functions so
# production = TEE at age_independent
# returns a new, normalized function
normalize_ind <- function(fprod, fwt, sex = 1, age_independent = 18, meat_kcal = 1500){
  wt <- fwt(age_independent, sex)
  a <- TEE(wt, age_independent, sex)/fprod(age_independent)
  function(age) a * fprod(age)
}

# factor to multiply production function so
# mean adult production = specified value
# returns a new, normalized function
normalize_adult <- function(fprod, adult_avg, min_age = 18, max_age = 60){
  a <- adult_avg / mean(fprod(min_age:max_age))
  function(age) a * fprod(age)
}

# From Hawkes et al. 1982
# The Ache eat every edible bit of an animal. We have estimated this to be 65 percent of the live
# weight for mammals and birds, 70 percent for reptiles and fish. Caloric values for most mammals are
# estimated at 300 Ca1/100 g edible portion (Meehan 1977; Lee 1979). Deer and monkey are estimated at
# 125 Ca1/100 g and 200 Ca1/100 g edible portion, respectively; birds at 190 Ca1/100 g; reptiles and fish at
# 150 Ca1/100 g and 137 Ca1/100 g. All of these estimates are derived from the USDA Agricultural Handbook
# No. 456, or Meehan (1977); the former is also the source of our caloric figure for honey. The
# caloric content of palm fiber may be inaccurately estimated-we have used the result from an
# analysis of the nutritional constituents of the liquid squeezed by hand
# (sucking may extract more from the fiber).


# Gurven and Kaplan 2006 -------------------------------------------
# Determinants of Time Allocation across the Lifespan

# Defaults for theoretical model
hg_strengthzGK2006 <- function(age, a1=1, b1=1.87, b2=-0.06){
  a1 * age^b1 * exp(b2 * age)
}

hg_skillGK2006 <- function(age, a2=1, b3=1.5, b4=-0.026){
  a2 * age^b3 * exp(b4 * age)
}

hg_productivityGK2006 <- function(age, alpha=1, beta=1){
  hg_strengthzGK2006(age)^alpha * hg_skillGK2006(age)^beta
}

# Empirical model for Machiguenga & Piro

# Male hunting params
piro_rate_m <- list(a = 4.8e-5, m = 6.837, n = -0.197)
piro_time_m <- list(a = 0.049, m = 3.157, n = -0.116)

# Female gathering params
piro_rate_f <- list(a = 0.006, m = 4.53, n = -0.132)
piro_time_f <- list(a = 3.354, m = 0.622, n = -0.003)

GK2006_machiguenga_piro <- function(age, p){
  p$a * age^p$m * exp(p$n * age)
}

GK2006_mp <- function(age, rate_params, time_params){
  GK2006_machiguenga_piro(age, rate_params) * GK2006_machiguenga_piro(age, time_params)
}

# McElreath and Koster 2014 ----------------------------------------------------

# Third degree polynomial model of men's hunting (kg meat)
# meat_kcal = 1700 based on Kaplan et al. 2000 (need better meat value)

ache_productivity_mcelreath <- function(age, meat_kcal = 1500){
  # convert age to sd units based on hunting.csv from McElreath & Koster 2014 SI
  age <- (age - 46.44)/13.55
  zeros <- inv_logit(0.04 + 0.14*age + 0.17*age^2 + -0.11*age^3)
  kg <- exp(1.83 + -0.06*age + -0.02*age^2 + 0.05*age^3)
  meat_kcal * (1 - zeros) * kg
}

# Predictions of Hagen GAM fit to Ache hunting data
ache_hunt_zi <- function(age, m = m_ache_zi, meat_kcal = 1500, type = "response", re.form = NULL){
    meat_kcal * glmmTMB:::predict.glmmTMB(m_ache_zi, newdata = data.frame(age = age, id = NA), re.form = re.form, type = 'response')
  }

# Koster et al. 2020 -----------------------------------------------

# Theoretical

# Declining component:
# M(age) = exp(-m * age) # m is rate of decline

M <- function(age, m){
  exp(-m * age)
}

# Growth component
# K(age) = 1 - exp(-k * age) # k is rate of increase

K <- function(age, k){
  1 - exp(-k * age)
}

# skill(age) = M(age) * K(age)^b

skill_koster2020 <- function(age, k, m, b, a = 1){
  age <- age / 80
  a * M(age, m) * K(age, k)^b
}

# harvest = skill^eta * labor^beta * alpha
# eta = skill elasticity
# beta = labor elasticity
# α is a linear model for covariates such as group size,
# the number of assistants, the use of dogs, and the use of firearms.

harvest <- function(age, m, k, b, eta, beta, alpha, labor){
  skill_koster2020(age, m, k, b)^eta * labor^beta * alpha
}

# Productivity = harvest * probability of success (age)

success <- function(age){
  0.5 #
}

productivity <- function(age, m, k, b, eta, beta, alpha, labor){
  success(age) * harvest(age, m, k, b, eta, beta, alpha, labor)
}

# Empirical
# units are kg meat per trip

# get_koster_model <- function(){
#   load(here("data", "Koster osfstorage-archive/model_fix_17092018.Rdata"))
#   rstan::extract(mfit)
# }

# Code adapted from cchunts

#' @title Skill Koster et al 2020
#' @description Relative skill as a function of age from Koster et al. (2020)
#' @param age Age in years
#' @return Relative skill level
#' @details See Koster et al. (2020)
#' @examples
#' \dontrun{
#' if(interactive()){
#'  koster2020skill(0:80)
#'  }
#' }
#' @rdname koster2020skill
#' @export
koster2020skill <- function(age){
  if (min(age) < 0 || max(age) > 80) stop("Ages must be between 0 and 80 years")

  k <- exp(koster2020posterior$lifehistmeans[, 1])
  m <- exp(koster2020posterior$lifehistmeans[, 2])
  b <- exp(koster2020posterior$lifehistmeans[, 3])

  # x_seq = seq(from = 0, to = 1, length.out = 81)
  skill_at_x <- sapply(age, function(x) skill_koster2020(x, k, m, b))
  skill <- apply(skill_at_x, 2, mean)
  skill / max(skill)
}

# According to the code, age 0.5 = 40, so 1 = 80 years old
# By leaving soc_id = j, hours = h_use, hunter_id as NULL
# (i.e., not in data.frame) the function will compute the average (kg per hour?)

get_koster2020production <- function(post){
  p <- cch_predict(post, data = data.frame(age = seq(0, 1, 1/80)), raw = F)
  Ep <- (1 - p$failure) * p$harvest
  apply(Ep, 2, mean)
}

# j <- 16 # ache
# cch_plot_curve(
#   map_id=j,
#   col="#D4ACE8",
#   alpha=0.5 ,
#   pts=TRUE,
#   sample_post=FALSE,
#   lwd=1,
#   skill_only=FALSE ,
#   main="Ache - production"
#   )

koster2020hagen <- function(age, m, re.form = NULL){
  predict(
    m,
    newdata = data.frame(age = age, independent = T, society = NA, uniqueid = NA),
    re.form = re.form,
    type = 'response',
    allow.new.levels=TRUE
  )
}

hadza_production_func <- function(age, sex){
    Sex <- factor(ifelse(sex == 0, "Female", "Male"))
    p <- mgcv::predict.gam(m_hadza_kcals, newdata = data.frame(Age = age, Sex = Sex), type = 'response')
    as.numeric(ifelse(age < 3, 0, p))
}

# units: kcals
# kaplan2000_hadza_male
# kaplan2000_hadza_female

# units: kcals
# hadza_production # my gam of marlowe's plot


# kaplan2000_ache_female

# units: kg of meat
# ache_productivity_mcelreath # 2014

# units: kg of meat
# ache_hunt_zi # my model of ache hunting data

# units: kcals? arbitrary
# GK2006_hg_productivity

# units: ? Need to look this up
# GK2006_mp(age, rate_params, time_params)

# units: kg of meat/hr
# koster2020production # their model of cchunts data

# units: kg meat/hr
# koster2020hagen # my model of cchunts data


# Generic production functions --------------------------------------------

production_male <- function(age, group){
  if (group == "ache"){
    # return(ache_productivity_mcelreath(age))
    return(ache_hunt_zi(age))
  } else if (group == "hadza"){
    # return(hadza_production(age, 1))
    return(kaplan2000_hadza_male(age))
  } else if (group == "kung"){
    return(normalize_adult(koster2020production, adult_avg = 3203)(age))
  } else if (group == "avg"){
    return(normalize_adult(koster2020production, adult_avg = 3840)(age)) # Avg HG_mean
  }
}

production_female <- function(age, group){
  if (group == "ache"){
    return(kaplan2000_ache_female(age))
  } else if (group == "hadza"){
    # return(hadza_production(age, 0)) # Brought into camp
    return(kaplan2000_hadza_female(age))
  } else if (group == "kung"){
    return(normalize_adult(\(age) hadza_production(age, 0), adult_avg = 3660)(age))
  } else if (group == "avg"){
    return(normalize_adult(\(age) hadza_production(age, 0), adult_avg = 2082)(age)) # Avg HG_mean
  }
}

production_child <- function(age, group){
  if (group == "ache"){
    return(kaplan2000_ache_female(age))
  } else if (group == "hadza"){
    # return((hadza_production(age, 0) + hadza_production(age, 1))/2) # Brought into camp
    return((kaplan2000_hadza_female(age) + kaplan2000_hadza_male(age))/2)
  } else if (group == "kung"){
    return(normalize_adult(\(age) hadza_production(age, 0), adult_avg = 3660)(age))
  } else if (group == "avg"){
    return((hadza_production(age, 0) + hadza_production(age, 1))/2)
  }
}

hg_strength <- function(age, sex, pregnant = 0, group = 'avg', max_age = 85, b=-0.15){
  wt <- hg_weight(age, sex, pregnant, group)
  max_wt <- max(hg_weight(20:60, 1, 0, group))
  wt * (1 - exp(b * (max_age-age))) / max_wt
}

# OLD
# b1 = 0.25: slow skill acquisition (midpoint ~ 28; asymptotes around 40-45)
# b1 = 0.50: fast skill acquisition (midpoint ~ 14-15; asymptotes around 25)
# hg_skill <- function(age, b0 = -7, b1, b2 = -0.15, max_age = 85){
#   a <- 0.001
#   y <- b1 * age
#   ifelse(age < 3, 0, ( a*exp(y)/(1 + a*exp(y))  ) * (1 - exp(b2 * (max_age-age))))
# }

# Fast: age50 = 10, age95 = 20, b1 = 0.3
# Medium: age50 = 15, age95 = 30, b1 = 0.2
# Slow: age50 = 20, age95 = 40, b1 = 0.15

#' @title Hunter-gatherer skill
#' @description Skill given age, rate of acquisition, and age of 50% level
#' @param age Age in years
#' @param b1 Rate of skill acquisition, e.g., 0.15 is slow, 0.40 is fast
#' @param age50 Age when skill is 50% of maximum, Default: 15
#' @param b2 Rate of skill decline in old age, Default: -0.15
#' @param max_age Maximum age, Default: 85
#' @return Skill level between 0 and 1
#' @details Skill is computed as \code{ifelse(age < 3, 0, ( 1/(1 + exp(-b1 * (age - age50)))  ) * (1 - exp(b2 * (max_age-age))))}
#' @examples
#' \dontrun{
#' if(interactive()){
#'  hg_skill(c(5, 10, 15), b1 = 0.25)
#'  }
#' }
#' @rdname hg_skill
#' @export
hg_skill <- function(age, b1, age50 = 15, b2 = -0.15, max_age = 85){
  ifelse(age < 3, 0, ( 1/(1 + exp(-b1 * (age - age50)))  ) * (1 - exp(b2 * (max_age-age))))
}

# alpha = 0: all skill; alpha = 1: all strength
# TEE_prop = TEE of this many adult men (TEE ~ 2800 kcals)

#' @title Hunter-gatherer productivity
#' @description Productivity as a function of age, sex, and group
#' @param age Age in years
#' @param sex 0: female; 1: male
#' @param TEE_prop Proportion of adult TEE, Default: 2
#' @param alpha Relative importance of strength vs skill. 0: all skill; 1: all strength, Default: 1
#' @param b1 Rate of skill acquisition, Default: 0.25
#' @param age50 Age at which skill is 50% of adult
#' @param group character 'ache', 'hadza', 'kung', 'avg'. Default: 'avg'
#' @return Productivity as a proportion of adult TEE
#' @details \code{TEE_prop * hg_strength^alpha * hg_skill^(1-alpha)}
#' @examples
#' \dontrun{
#' if(interactive()){
#'  hg_productivity(c(5, 15, 25), 0, TEE_prop = 1.5, alpha = 0.5, age50 = 15)
#'  }
#' }
#' @rdname hg_productivity
#' @export
hg_productivity <- function(age, sex, TEE_prop = 2, alpha=1, b1 = 0.25, age50, group = "avg"){
  ifelse(
    age > 80, # age of death
    0,
    TEE_prop * hg_strength(age, sex, group = group)^alpha * hg_skill(age, b1=b1, age50=age50)^(1 - alpha)
  )
}

#' @title Average hunter-gatherer productivity
#' @description Productivity as a function of age and group averaged across sex
#' @param age Age in years
#' @param sex 0: female; 1: male
#' @param TEE_prop Proportion of adult TEE, Default: 2
#' @param alpha Relative importance of strength vs skill. 0: all skill; 1: all strength, Default: 1
#' @param b1 Rate of skill acquisition, Default: 0.25
#' @param age50 Age at which skill is 50% of adult
#' @param group character 'ache', 'hadza', 'kung', 'avg'. Default: 'avg'
#' @return Productivity as a proportion of adult TEE
#' @details \code{TEE_prop * hg_strength^alpha * hg_skill^(1-alpha)}
#' @examples
#' \dontrun{
#' if(interactive()){
#'  hg_productivity(c(5, 15, 25), TEE_prop_f = 1.5, TEE_prop_m = 2.5, alpha = 0.5, age50 = 15)
#'  }
#' }
#' @rdname hg_productivity_avg
#' @export
hg_productivity_avg <- function(age, TEE_prop_f, alpha_f, b1_f, age50_f, TEE_prop_m, alpha_m, b1_m, age50_m){
  pf <- hg_productivity(age, 0, TEE_prop_f, alpha_f, b1_f, age50_f)
  pm <- hg_productivity(age, 1, TEE_prop_m, alpha_m, b1_m, age50_m)
  return((pf + pm)/2)
}

