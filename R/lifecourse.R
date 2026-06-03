

# Lifecourse functions -----------------------------------------------------

# child_values can be a vector of child_survival, childTEE, or child_production
sum_resident_children <- function(wife_age, afb, births, child_values){

  # Window: maximum child residence is afb - 1 years, i.e., age 0:(afb-1), or index 1:afb
  wife_age_then <- ifelse(wife_age - afb < afb, afb, wife_age - afb + 1)

  wife_now_index <- wife_age - afb + 1
  wife_then_index <- wife_age_then - afb + 1

  # births backwards from current wife age for afb years or to afb
  birth_window <- births[wife_now_index:wife_then_index]

  return(
    sum(
      birth_window * child_values[1:length(birth_window)]
    )
  )
}

#' @title Hunter-gatherer lifecourse
#' @description Simulation of a nuclear family from the wife's age at first birth to the maximum human lifespan
#' @param SRB Sex ratio at birth, Default: 1.05
#' @param afb Age at first birth, Default: 18
#' @param max_age Maximum human lifespan in years, Default: 80
#' @param IBI Interbirth interval (currently cannot be changed), Default: 3
#' @param preg_cost The cost of pregnancy (currently not used), Default: 0.1
#' @param alb The age at last birth in years (80: no menopause), Default: 80
#' @param group character, one of 'ache', 'hadza', 'kung', 'avg', Default: 'avg'
#' @param e0_f Female expected lifespan at birth, Default: 35
#' @param e0_m Male expected lifespan at birth , Default: 30
#' @param TEE_prop_m Male energy production as a proportion of adult male TEE, Default: 2
#' @param alpha_m The relative contribution of skill and strength to productivity (0: all skill; 1: all strength), Default: 0.5
#' @param b1_m The rate of skill acquisition (0.2: slow; 0.4: fast), Default: 0.25
#' @param age50_m The age in years at which male skill reaches 50% of maximum, Default: 15
#' @param TEE_prop_f Female energy production as a proportion of adult male (yes, male) TEE, Default: 1
#' @param alpha_f The relative contribution of skill and strength to productivity (0: all skill; 1: all strength), Default: 0.5
#' @param b1_f The rate of skill acquisition (0.2: slow; 0.4: fast), Default: 0.25
#' @param age50_f The age in years at which female skill reaches 50% of maximum, Default: 15
#' @return A tibble (data frame)
#' @details
#' The function returns a tibble (data frame) with these columns and one row
#' for each year of the wife's age:
#'
#' \code{wife_age, husband_age, wife_num, husband_num,
#' wife_survival, husband_survival, husband_wife_ratio, wife_ex, husband_ex,
#' pregnancy, births, fertility, child_age, child_index, girl_survival, boy_survival,
#' child_survival, wifeTEE, husbandTEE, childTEE, wife_production, husband_production,
#' child_production, resident_children, total_child_consumption, total_child_production,
#' family_size, family_consumption, family_production, energy_balance.}
#'
#' The age at first birth (afb) is also the age of child independence.
#'
#' The resident_children, total_child_consumption, total_child_production values are for
#' resident children only, and therefore are computed for the wife's latest birth
#' to her oldest resident child using the internal function \code{sum_resident_children}.
#' Additionally, whereas the other child values are conditional on the child being born,
#' these values are conditional on the wife surviving to give birth.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  hg_lifecourse()
#'  }
#' }
#' @rdname hg_lifecourse
#' @export
hg_lifecourse <- function(
  SRB = 1.05,
  afb = 18, # Age of first birth; also age of independence
  max_age = 80,
  IBI = 3,
  preg_cost = 0.1, # Not used
  alb = 80,
  group = "avg",
  e0_f = 35,
  e0_m = 30,

  TEE_prop_m = 2,
  alpha_m = 0.5,
  b1_m = 0.25,
  age50_m = 15,

  TEE_prop_f = 1,
  alpha_f = 0.5,
  b1_f = 0.25,
  age50_f = 15,
  ... # To absorb unused args
){
  if (alb <= afb) stop("alb is less than or equal to afb (age at first birth)")

  # Simulation starts at wife_age = afb and ends at wife_age = max_age
  num_ages <- max_age - afb + 1

  lt_f <- lifetable(group, sex = 0, e0 = e0_f, SRB = SRB)
  lt_m <- lifetable(group, sex = 1, e0 = e0_m, SRB = SRB)

  normalize <- 1 + SRB

  TEEadult_f <- mean(TEE2(20:60, 0, group))
  TEEadult_m <- mean(TEE2(20:60, 1, group))
  TEEadult_avg <- (TEEadult_f + TEEadult_m)/2

  preg_schedule <- c(T, rep(F, IBI - 1)) # Integer IBI
  pregnancies <- rep(preg_schedule, ceiling((max_age-afb+2)/IBI)+1)
  alb_index <- alb - afb + 2
  pregnancies[alb_index:length(pregnancies)] <- F
  births <- c(F, pregnancies) # Birth occurs in the year after pregnancy
  # Trim to total years, if necessary
  births <- births[2:(num_ages+1)]
  pregnancies <- pregnancies[2:(num_ages + 1)]

  tibble::tibble(
    wife_age = afb:max_age,
    husband_age = wife_age,
    wife_num = lt_f$lx[wife_age + 1],
    husband_num = lt_m$lx[husband_age + 1],
    wife_survival = wife_num/wife_num[1], # Simulation starts conditional on female survival to afb
    husband_survival = husband_num / wife_num[1], # Wives start out with less than 1 husband
    parents = wife_survival + husband_survival,
    # husband_survival = husband_num/husband_num[1],
    # husband_wife_ratio = husband_num / wife_num,
    wife_ex = lt_f$ex[wife_age + 1], # lifespan remaining
    husband_ex = lt_m$ex[husband_age + 1], # lifespan remaining
    pregnancy = pregnancies, #rep(c(F, F, T), 28)[1:num_ages],
    # pregnancy_cost = 0.05*cumsum(preg_cost*pregnancy)^2, # Needs work
    births = births * wife_survival, # births accounting for wife mortality
    fertility = cumsum(births > 0), # Total births for wives who survive to age x
    fertility2 = cumsum(births), # Taking wife mortality into account
    child_age = wife_age - afb, # Starts at 0
    child_index = child_age + 1, # Age + 1
    girl_survival = lt_f$lx[child_index] / normalize, # Presumably, mother survival is baked in to these values
    boy_survival = lt_m$lx[child_index] / normalize,
    child_survival = boy_survival + girl_survival,

    wifeTEE = wife_survival * TEE2(wife_age, 0, group, pregnancy),
    husbandTEE = husband_survival * TEE2(husband_age, 1, group),
    # childTEE = child_survival * TEE2(child_age, NA, group),
    girlTEE = girl_survival * TEE2(child_age, 0, group),
    boyTEE = boy_survival * TEE2(child_age, 1, group),

    # Production of both sexes as a proportion of average adult male TEE
    wife_production = wife_survival * TEEadult_avg * hg_productivity(wife_age, 0, TEE_prop = TEE_prop_f, alpha = alpha_f, b1 = b1_f, age50 = age50_f),
    husband_production = husband_survival * TEEadult_avg * hg_productivity(husband_age, 1, TEE_prop = TEE_prop_m, alpha = alpha_m, b1 = b1_m, age50 = age50_m),
    girl_production = girl_survival * TEEadult_avg * hg_productivity(child_age, 0, TEE_prop = TEE_prop_f, alpha = alpha_f, b1 = b1_f, age50 = age50_f),
    boy_production = boy_survival * TEEadult_avg * hg_productivity(child_age, 1, TEE_prop = TEE_prop_m, alpha = alpha_m, b1 = b1_m, age50 = age50_m),

    # The following values are for resident children only, and therefore must
    # be computed for the wife's latest birth to her oldest resident child.
    # Additionally, the following values are conditional on the wife
    # surviving to give birth.
    resident_girls = purrr::map_dbl(wife_age, \(wifeage) sum_resident_children(wifeage, afb, births, girl_survival)),
    resident_boys = purrr::map_dbl(wife_age, \(wifeage) sum_resident_children(wifeage, afb, births, boy_survival)),
    resident_children = resident_girls + resident_boys,

    total_girl_consumption = purrr::map_dbl(wife_age, \(wifeage) sum_resident_children(wifeage, afb, births, girlTEE)),
    total_boy_consumption = purrr::map_dbl(wife_age, \(wifeage) sum_resident_children(wifeage, afb, births, boyTEE)),
    total_child_consumption = total_girl_consumption + total_boy_consumption,

    total_girl_production = purrr::map_dbl(wife_age, \(wifeage) sum_resident_children(wifeage, afb, births, girl_production)),
    total_boy_production = purrr::map_dbl(wife_age, \(wifeage) sum_resident_children(wifeage, afb, births, boy_production)),
    total_child_production = total_girl_production + total_boy_production,

    family_size = resident_girls + resident_boys + wife_survival + husband_survival,
    family_consumption = total_girl_consumption + total_boy_consumption + wifeTEE + husbandTEE,
    family_production = total_girl_production + total_boy_production + wife_production + husband_production,
    energy_balance = family_production - family_consumption,

    energy_balance2 = wife_production + husband_production - family_consumption,
    cumulativeEB = cumsum(energy_balance)
  )
}

#' @title hg_lifecourse2
#' @description Similar to hg_lifecourse but without creating families
#' @param SRB Sex ratio at birth, Default: 1.05
#' @param max_age Maximum age, Default: 80
#' @param group one of 'ache', 'hadza', 'kung', 'avg', Default: 'avg'
#' @param e0_f Expected female lifespan at birth, Default: 35
#' @param e0_m Expected male lifespan at birth, Default: 30
#' @param TEE_prop_m Male productivity as a proportion of average adult TEE, Default: 2
#' @param alpha_m The relative contribution of skill and strength to productivity (0: all skill; 1: all strength), Default: 0.5
#' @param b1_m Male rate of skill acquisition, Default: 0.25
#' @param age50_m Age at which male skill is 50% of maximum, Default: 15
#' @param TEE_prop_f Female productivity as a proportion of average adult TEE, Default: 1
#' @param alpha_f The relative contribution of skill and strength to productivity (0: all skill; 1: all strength), Default: 0.5
#' @param b1_f Female rate of skill acquisition, Default: 0.25
#' @param age50_f Age at which female skill is 50% of maximum, Default: 15
#' @return A tibble (data frame)
#' @details One row for each age from birth to death
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #hg_lifecourse2()
#'  }
#' }
#' @rdname hg_lifecourse2
#' @export
hg_lifecourse2 <- function(
    SRB = 1.05,
    afb = 20,
    max_age = 80,
    group = "avg",
    e0_f = 35,
    e0_m = 30,

    TEE_prop_m = 2,
    alpha_m = 0.5,
    b1_m = 0.25,
    age50_m = 15,

    TEE_prop_f = 1,
    alpha_f = 0.5,
    b1_f = 0.25,
    age50_f = 15,
    ... # To absorb unused params
  ){

  lt_f <- lifetable(group, sex = 0, e0 = e0_f, SRB = SRB)
  lt_m <- lifetable(group, sex = 1, e0 = e0_m, SRB = SRB)

  normalize <- 1 + SRB

  TEEadult_f <- mean(TEE2(20:60, 0, group))
  TEEadult_m <- mean(TEE2(20:60, 1, group))
  TEEadult_avg <- (TEEadult_f + TEEadult_m)/2

  tibble::tibble(
    age = 0:max_age,
    age_index = age + 1,
    survival_f = lt_f$lx[age_index] / normalize,
    survival_m = lt_m$lx[age_index] / normalize,
    TEE_f = survival_f * TEE2(age, 0, group),
    TEE_m = survival_m * TEE2(age, 1, group),
    production_f = survival_f * TEEadult_avg * hg_productivity(age, 0, TEE_prop = TEE_prop_f, alpha = alpha_f, b1 = b1_f, age50 = age50_f),
    production_m = survival_m * TEEadult_avg * hg_productivity(age, 1, TEE_prop = TEE_prop_m, alpha = alpha_m, b1 = b1_m, age50 = age50_m),
    energy_balance = production_f + production_m - TEE_f - TEE_m,
    production_adult_f = ifelse(age < afb, 0, production_f),
    production_adult_m = ifelse(age < afb, 0, production_m),
    energy_balance_adult_prod = production_adult_m + production_adult_f - TEE_m - TEE_f,
  )
}

