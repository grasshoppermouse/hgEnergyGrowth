library(gamm4)
library(mgcv)

# JPPS 7 parameter model of human growth curves ---------------------------

#' @title JPPS human growth model
#' @description JPPS 7 parameter model of human growth curves
#' @param age Age in years
#' @param params A numeric vector c(B, C1, C2, C3, D1, D2, D3)
#' @return Weight in kilograms
#' @details Jolicoeur, Pontier, Pernin, and  Sempe (1988) A Lifetime Asymptotic Growth Curve for Human Height
#' @examples
#' \dontrun{
#' if(interactive()){
#'  Ache_female <- c(B=55.63, C1=15.60, C2=3.58, C3=0.67, D1=13.60, D2=14.24, D3=26.95)
#'  jpps(25, Ache_female)
#'  }
#' }
#' @rdname jpps
#' @export
jpps <- function(age, params){
  age_t <- age + 0.75 # Adding gestation
  params[1] * (1 - 1/(1 + sum((age_t / params[5:7])^params[2:4])))
}

# Kung --------------------------------------------------------------------

#' @title !Kung weights
#' @description Weight in kilograms of !Kung by age, sex, and pregnancy
#' @param age Age in years
#' @param sex 0: female, 1: male
#' @param pregnant Pregnant = 1, Default: 0
#' @return Weight in kilograms
#' @details Predictions of GAM model fit to Howell data
#' @examples
#' \dontrun{
#' if(interactive()){
#'  kung_wt(c(5, 10, 15), 0, 0) # Female
#'  }
#' }
#' @rdname kung_wt
#' @export
kung_wt <- function(age, sex, pregnant = 0){
  pregnant <- factor(ifelse(pregnant == 0, "No", "Yes"))
  sex <- ifelse(sex == 0, "Female", "Male")
  as.numeric(mgcv::predict.gam(m_kung_weight$gam, newdata = data.frame(age=age, sex2=factor(sex), pregnant = pregnant), re.form = NA))
}

#' @title !Kung heights
#' @description Height in cm of !Kung by age and sex
#' @param age Age in years
#' @param sex 0: female, 1: male
#' @return Height in centimeters
#' @details Predictions of GAM model fit to Howell data
#' @examples
#' \dontrun{
#' if(interactive()){
#'  kung_wt(c(5, 10, 15), 0) # Female
#'  }
#' }
#' @rdname kung_ht
#' @export
kung_ht <- function(age, sex){
  sex <- ifelse(sex == 0, "Female", "Male")
  as.numeric(mgcv::predict.gam(m_kung_height$gam, newdata = data.frame(age=age, sex2=factor(sex)), re.form = NA))
}

# Ache --------------------------------------------------------------------

# https://onlinelibrary.wiley.com/doi/epdf/10.1002/ajpa.20306

#' @title Ache weights
#' @description Ache weight by age, sex, and. pregnancy status
#' @param age Age in years
#' @param sex 0: female, 1: male
#' @param pregnant Pregnant: 1, Default: 0
#' @return Weight in kilograms
#' @details Uses parameters of the JPPS model fit to Ache data in https://onlinelibrary.wiley.com/doi/epdf/10.1002/ajpa.20306
#' @examples
#' \dontrun{
#' if(interactive()){
#'  ache_wt(c(5, 10, 15), 1) # Male
#'  }
#' }
#' @rdname ache_wt
#' @export
ache_wt <- function(age, sex, pregnant=0){
  # Female parameters
  ache_f <- c(B=55.63, C1=15.60, C2=3.58, C3=0.67, D1=13.60, D2=14.24, D3=26.95)

  # Male parameters
  ache_m <- c(B=59.47, C1=15.70, C2=3.00, C3=0.56, D1=15.24, D2=15.74, D3=50.18)

  wt <- purrr::map2_dbl(age, sex, \(age, sex) {if (sex == 0) params <- ache_f else params <- ache_m; jpps(age, params)})
  if (identical(0, pregnant)) return(wt)
  wt <- ifelse(pregnant == 1, wt + 2.4, wt) # Value from !Kung
  return(wt)
}

#' @title Ache heights
#' @description Ache height by age and sex
#' @param age Age in years
#' @param sex 0: female, 1: male
#' @return Height in centimeters
#' @details Based on a GAM fit to data points extracted from a plot in Ache Life History
#' @examples
#' \dontrun{
#' if(interactive()){
#'  ache_ht(c(5, 10, 15), 1) # Male
#'  }
#' }
#' @rdname ache_ht
#' @export
ache_ht <- function(age, sex){
  sex <- ifelse (sex == 0, "Female", "Male")
  as.numeric(mgcv::predict.gam(m_ache_height, newdata = data.frame(age=age, Sex=factor(sex))))
}

# Hadza -------------------------------------------------------------------

#' @title Hadza weights
#' @description Weight of Hadza by age and sex
#' @param age Age in years
#' @param sex 0: female, 1: male
#' @param pregnant Pregnant: 1, Default: 0
#' @return Weight in kilograms
#' @details  Predictions from a GAM model fit to data digitized from a scatterplot in the growth SI to Demography And Evolutionary Ecology Of Hadza (Blurton-Jones).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  hadza_wt(c(25, 30, 35), 0, 1) # Pregnant women
#'  }
#' }
#' @rdname hadza_wt
#' @export
hadza_wt <- function(age, sex, pregnant = 0){
    sex <- ifelse (sex == 0, "Female", "Male")
    wt <- as.numeric(mgcv::predict.gam(m_hadza_weight, newdata = data.frame(Age=age, Sex=factor(sex))))
    if (identical(0, pregnant)) return(wt)
    wt <- ifelse(pregnant == 1, wt + 2.4, wt) # value from !kung
    return(wt)
  }

#' @title Hadza heights
#' @description Height of Hadza by age and sex
#' @param age Age in years
#' @param sex 0: female, 1: male
#' @return Height in centimeters
#' @details  Predictions from a GAM model fit to data digitized from a scatterplot in the growth SI to Demography And Evolutionary Ecology Of Hadza (Blurton-Jones).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  hadza_ht(c(5, 10, 15), 0) # Female
#'  }
#' }
#' @rdname hadza_ht
#' @export
hadza_ht <- function(age, sex){
  sex <- ifelse (sex == 0, "female", "male")
  as.numeric(mgcv::predict.gam(m_hadza_height, newdata = data.frame(age=age, sex=factor(sex))))
}

# HG weight ---------------------------------------------------------------

hg_weight0 <- function(age, sex, pregnant, group){
  dplyr::case_when(
    group == "ache" ~ ache_wt(age, sex, pregnant),
    group == "hadza" ~ hadza_wt(age, sex, pregnant),
    group == "kung" ~ kung_wt(age, sex, pregnant),
    group == "avg" ~ (ache_wt(age, sex, pregnant) + hadza_wt(age, sex, pregnant) + kung_wt(age, sex, pregnant))/3
  )
}

#' @title Hunter-gatherer weights
#' @description Hunter gather weight by age and sex
#' @param age Age in years
#' @param sex 0: female, 1: male
#' @param pregnant Pregnant: 1, Default: 0
#' @param group One of 'ache', 'kung', 'hadza', 'avg'
#' @return Weight in kilograms
#' @details 'avg': weight averaged across the three populations
#' @examples
#' \dontrun{
#' if(interactive()){
#'  hg_weight(c(5, 10, 15), 1, 0, "avg") # The weights of 3 "average" hunter-gatherer boys
#'  }
#' }
#' @rdname hg_weight
#' @export
hg_weight <- function(age, sex, pregnant, group){
  if(length(sex) == 1 && is.na(sex)){
    wt_female <- hg_weight0(age, 0, pregnant, group)
    wt_male <- hg_weight0(age, 1, pregnant, group)
    return((wt_male + wt_female)/2)
  }
  hg_weight0(age, sex, pregnant, group)
}

# HG height ---------------------------------------------------------------

hg_height0 <- function(age, sex, group){
  dplyr::case_when(
    group == "ache" ~ ache_ht(age, sex),
    group == "hadza" ~ hadza_ht(age, sex),
    group == "kung" ~ kung_ht(age, sex),
    group == "avg" ~ (ache_ht(age, sex) + hadza_ht(age, sex) + kung_ht(age, sex))/3
  )
}

#' @title Hunter-gatherer heights
#' @description Hunter gather height by age and sex
#' @param age Age in years
#' @param sex 0: female, 1: male
#' @param group One of 'ache', 'kung', 'hadza', 'avg'
#' @return Height in centimeters
#' @details 'avg': height averaged across the three populations
#' @examples
#' \dontrun{
#' if(interactive()){
#'  hg_weight(c(5, 10, 15), 1, 0, "avg") # The heights of 3 "average" hunter-gatherer boys
#'  }
#' }
#' @rdname hg_height
#' @export
hg_height <- function(age, sex, group){
  if(length(sex) == 1 && is.na(sex)){
    ht_female <- hg_height0(age, 0, group)
    ht_male <- hg_height0(age, 1, group)
    return((ht_male + ht_female)/2)
  }
  hg_height0(age, sex, group)
}
