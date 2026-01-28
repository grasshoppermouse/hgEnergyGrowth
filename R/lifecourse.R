

# Lifecourse functions -----------------------------------------------------

# child_values can be vector of child_survival, childTEE, or child_production
sum_resident_children <- function(wife_age, afb, menopause_age, births, child_values){

  # Window: maximum child residence is afb years, i.e., age 0:(afb-1), or index 1:afb
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
#' @param age_gap The age gap at marriage (e.g., 5: husband is 5 years older than the wife)
#' @param SRB Sex ratio at birth, Default: 1.05
#' @param afb Age at first birth, Default: 18
#' @param max_age Maximum human lifespan in years, Default: 80
#' @param IBI Interbirth interval (currently cannot be changed), Default: 3
#' @param preg_cost The cost of pregnancy (currently not used), Default: 0.1
#' @param menopause_age The age at menopause in years (80: no menopause), Default: 80
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
  age_gap = 5, # Marriage age gap
  SRB = 1.05,
  afb = 18, # Age of first birth; also age of independence
  max_age = 80,
  IBI = 3,
  preg_cost = 0.1,
  menopause_age = 80,
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
  age50_f = 15
){
  if (menopause_age <= afb) stop("menopause_age is less than or equal to afb (age at first birth)")

  # Simulation starts at wife_age = afb and ends at wife_age = max_age
  num_ages <- max_age - afb + 1

  lt_f <- lifetable(group, sex = 0, e0 = e0_f)
  lt_m <- lifetable(group, sex = 1, e0 = e0_m)

  TEEadult_f <- mean(TEE2(20:60, 0, group))
  TEEadult_m <- mean(TEE2(20:60, 1, group))
  TEEadult_avg <- (TEEadult_f + TEEadult_m)/2

  preg_schedule <- c(T, rep(F, IBI - 1)) # Integer IBI
  pregnancies <- rep(preg_schedule, ceiling((max_age-afb+2)/IBI)+1)
  menopause_index <- menopause_age - afb + 2
  pregnancies[menopause_index:length(pregnancies)] <- F
  births <- c(F, pregnancies) # Birth occurs in the year after pregnancy
  # Trim to total years, if necessary
  births <- births[2:(num_ages+1)]
  pregnancies <- pregnancies[2:(num_ages + 1)]

  tibble::tibble(
    wife_age = afb:max_age,
    husband_age = wife_age + age_gap,
    wife_num = lt_f$lx[wife_age + 1],
    husband_num = SRB * lt_m$lx[husband_age + 1],
    wife_survival = wife_num/wife_num[1],
    husband_survival = husband_num/husband_num[1],
    husband_wife_ratio = husband_num / wife_num,
    wife_ex = lt_f$ex[wife_age + 1], # lifespan remaining
    husband_ex = lt_m$ex[husband_age + 1], # lifespan remaining
    pregnancy = pregnancies, #rep(c(F, F, T), 28)[1:num_ages],
    # pregnancy_cost = 0.05*cumsum(preg_cost*pregnancy)^2, # Needs work
    births = births * wife_survival, # births accounting for wife mortality
    fertility = cumsum(births > 0), # Total births for wives who survive to age x
    fertility2 = cumsum(births), # Taking wife mortality into account
    child_age = wife_age - afb, # Starts at 0
    child_index = child_age + 1, # Age + 1
    girl_survival = lt_f$lx[child_index], # Presumably, mother survival is baked in to these values
    boy_survival = SRB * lt_m$lx[child_index],
    child_survival = (boy_survival + girl_survival)/2,

    wifeTEE = wife_survival * TEE2(wife_age, 0, group, pregnancy),
    husbandTEE = husband_survival * TEE2(husband_age, 1, group),
    childTEE = child_survival * TEE2(child_age, NA, group),
    wife_production = wife_survival * TEEadult_f * hg_productivity(wife_age, 0, TEE_prop = TEE_prop_f, alpha = alpha_f, b1 = b1_f, age50 = age50_f),
    husband_production = husband_survival * TEEadult_m * hg_productivity(husband_age, 1, TEE_prop = TEE_prop_m, alpha = alpha_m, b1 = b1_m, age50 = age50_m),
    child_production = child_survival * TEEadult_avg * hg_productivity_avg(child_age, TEE_prop_f, alpha_f, b1_f, age50_f, TEE_prop_m, alpha_m, b1_m, age50_m),

    # The following 3 values are for resident children only, and therefore must
    # be computed for the wife's latest birth to her oldest resident child.
    # Additionally, whereas the preceding child values are conditional on the
    # child being born, the following values are conditional on the wife
    # surviving to give birth
    resident_children = purrr::map_dbl(wife_age, \(wifeage) sum_resident_children(wifeage, afb, menopause_age, births, child_survival)),
    total_child_consumption = purrr::map_dbl(wife_age, \(wifeage) sum_resident_children(wifeage, afb, menopause_age, births, childTEE)),
    total_child_production = purrr::map_dbl(wife_age, \(wifeage) sum_resident_children(wifeage, afb, menopause_age, births, child_production)),

    family_size = resident_children + wife_survival + husband_survival,
    family_consumption = total_child_consumption + wifeTEE + husbandTEE,
    family_production = total_child_production + wife_production + husband_production,
    energy_balance = family_production - family_consumption
  )
}
