# Code from Gaddy et al. 2025, DOI: 10.1073/pnas.2508091122
# https://www.pnas.org/doi/10.1073/pnas.2508091122

# Definitions
# https://www.youtube.com/@biodemography4195/videos
# https://ihmeuw-demographics.github.io/demCore/articles/introduction_to_life_tables.html
#
# e0: life expectancy at birth
# mx: age-specific mortality RATE (age x to x+n); ~ observed deaths / mid-interval population
# ax: number of years lived between x and x + n for those who died in the interval;
#     ax will usually close to 0.5, except for infants and the very elderly
# qx: PROBABILITY of death between age x and x+n for those who live to age x

# Set demographic functions

#function to build life table

#' @title Life table from age-specific mortality rate
#' @description Build a life table from mx, ax, sex, and SRB
#' @param mx age-specific mortality RATE (age x to x+n); ~ observed deaths / mid-interval population
#' @param ax number of years lived between x and x + n for those who died in the interval; ax will usually close to 0.5, except for infants and the very elderly
#' @param sex character "Female" or "Male"
#' @param SRB Sex ratio at birth (e.g., 1.05)
#' @return A data frame with x, ax, mx, qx, px, lx, dx, Lx, Tx, ex
#' @details Code from Gaddy et al. 2025, DOI: 10.1073/pnas.2508091122; https://www.pnas.org/doi/10.1073/pnas.2508091122
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname LT_from_mx
#' @export
LT_from_mx <- function(mx, ax, sex, SRB) {

  x = 0:130
  qx = mx / (1 + (1 - ax) * mx) # For constant mortality rate in the interval: qx = 1 - exp(-n*mx)
  qx[131] = 1
  px = 1 - qx

  # if (sex == "Female") {l0 = 100000}
  # if (sex == "Male") {l0 = 100000 * SRB}

  if (sex == "Female") {l0 = 1}
  if (sex == "Male") {l0 = SRB}

  lx = cumprod(c(l0, px[1:(length(px) - 1)]))
  dx = c(-diff(lx), lx[length(lx)])
  Lx = lx[-1] + ax[1:(length(ax) - 1)] * dx[1:(length(dx) - 1)]
  Lx = c(Lx, lx[length(lx)] / mx[length(mx)])
  Tx = rev(cumsum(rev(Lx)))
  ex = Tx / lx
  ax[131] = ex[131]
  LT = cbind.data.frame(x, ax, mx, qx, px, lx, dx, Lx, Tx, ex)

  return(LT)
}

#function to interpolate between life table quantities

#' @title Interpolate life table quantities
#' @description Interpolates lifetable quantities for fractional life expectancies, e.g, e0 = 30.7
#' @param col The life table quantity to interpolate
#' @param e0 Life expectancy at birth
#' @param sex character "Female" or "Male"
#' @return The interpolated vector
#' @details Code from Gaddy et al. 2025, DOI: 10.1073/pnas.2508091122; https://www.pnas.org/doi/10.1073/pnas.2508091122
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname lt_col_interpolate
#' @export
lt_col_interpolate <- function(col, e0, sex) {

  if (e0 %% 1 == 0) {

    desired <- modelLT |>
      dplyr::filter(Sex == sex, E0 == e0) |>
      dplyr::pull(!!col)

  } else {

    lower <- modelLT |>
      dplyr::filter(Sex == sex, E0 == floor(e0)) |>
      dplyr::pull(!!col)

    upper <- modelLT |>
      dplyr::filter(Sex == sex, E0 == floor(e0) + 1) |>
      dplyr::pull(!!col)

    diff <- upper - lower

    prop_through <- e0 - floor(e0)

    desired <- lower + (diff * prop_through)

  }
  return(desired)
}

#' @title UN Model Life Table
#' @description Returns a UN Model Life Table for a given life expectancy, sex, and SRB
#' @param e0 Life expectancy at birth
#' @param sex character "Female" or "Male"
#' @param SRB Sex ratio at birth, Default: 1
#' @return a UN Model Life Table
#' @details Code from Gaddy et al. 2025, DOI: 10.1073/pnas.2508091122; https://www.pnas.org/doi/10.1073/pnas.2508091122
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @rdname lifetable_UN
#' @export
lifetable_UN <- function(e0, sex, SRB = 1.0){
  mx <- lt_col_interpolate('mx', e0, sex)
  ax <- lt_col_interpolate('ax', e0, sex)
  LT_from_mx(mx = mx, ax = ax, sex = sex, SRB = SRB)
}

# lifetable_f <- function(e0_f, SRB = 1.0){
#   mx_f <- lt_col_interpolate('mx', e0_f, 'Female')
#   ax_f <- lt_col_interpolate('ax', e0_f, 'Female')
#   LT_from_mx(mx = mx_f, ax = ax_f, sex = "Female", SRB = SRB)
# }


# Siler for Hadza and Ache ------------------------------------------------

# From Gurven and Kaplan 2007
# Not sex-specific

# hadza <- list(age = 10:75, a1=0.351, b1=0.895, a2=0.011, a3=6.70e-06, b3=0.125)
# ache <- list(age = 10:75, a1=0.157, b1=0.721, a2=0.013, a3=4.80e-05, b3=0.103)
# avgHG <- list(age = 10:75, a1=0.422, b1=1.131, a2=0.013, a3=1.47e-04, b3=0.086)
#
# siler0 <- expression(a1 * exp(-b1 * age) + a2 + a3 * exp(b3 * age))
# siler <- function(age, a1, b1, a2, a3, b3) siler0
# dxsiler <- D(siler0, "age")
# d2xsiler <- D(dxsiler, "age")
#
# lx0 <- expression(exp(-(a1/b1) * (1 - exp(-b1*age))) * exp(-a2 * age) * exp((a3/b3) * (1 - exp(b3*age))))
# lx <- function(age, a1, b1, a2, a3, b3) lx0
# dxlx <- D(lx0, "age")
# d2xlx <- D(dxlx, "age")


# Hunter-gatherer life tables -------------------------------------------------------------

# Using UN life table for !Kung and avg

#' @title Life tables for hunter-gatherer groups
#' @description Returns a life table for Ache, Hadza, !Kung, or average
#' @param group character "ache", "hadza", "kung", "avg"
#' @param sex 0: female; 1: male
#' @param e0 Life expectancy at birth (years), Default: NULL
#' @return life table
#' @details Ache life table from Ache Life History; Hadza life table from Blurton-Jones (2016) Demography and evolutionary ecology of Hadza hunter-gatherers; !Kung and average life tables use code from Code from Gaddy et al. 2025, DOI: 10.1073/pnas.2508091122; https://www.pnas.org/doi/10.1073/pnas.2508091122
#' @examples
#' \dontrun{
#' if(interactive()){
#'  lifetable("hadza", 0, e0 = 35) # Female life table for e0 = 35 years
#'  }
#' }
#' @rdname lifetable
#' @export
lifetable <- function(group, sex, e0 = NULL){
  if (group == "ache"){
    if (sex == 0) return(ache_mortality_female)
    return(ache_mortality_male)
  } else if (group == "hadza"){
    if (sex == 0) return(hadza_mortality_female)
    return(hadza_mortality_male)
  } else if (group == "kung"){
    if (sex == 0) return(lifetable_UN(e0 = e0, sex = "Female"))
    return(lifetable_UN(e0 = e0, sex = "Male"))
  } else if (group == "avg"){
    if (sex == 0) return(lifetable_UN(e0 = e0, sex = "Female"))
    return(lifetable_UN(e0 = e0, sex = "Male"))
  } else stop(paste("Unknown group:", group))
}
