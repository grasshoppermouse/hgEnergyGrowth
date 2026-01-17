
# Equations from Pontzer et al. (2021). Science 373, 808–812, Table S2

#' @title TEE from body mass, age, sex, and percent fat free mass
#' @description TEE from body mass, age, sex, and percent fat free mass
#' @param bodymass Body mass in kilograms
#' @param age Age in years
#' @param sex 0: female; 1: male
#' @param FFMpct Percent fat free mass, Default: 0.85
#' @return Total energy expenditure (TEE) in kcals/day
#' @details Formula from Pontzer et al. (2021). Science 373, 808–812, Table S2
#' @examples
#' \dontrun{
#' if(interactive()){
#'  TEElog(50, 25, 0)
#'  }
#' }
#' @rdname TEElog
#' @export
TEElog <- function(bodymass, age, sex, FFMpct = 0.85){
  # TEE (MJ/day); bodymass (kg); sex (F=0, M=1); age (years)

  # Don't have FFM or FM, so estimate it
  FFM <- bodymass * FFMpct
  FM <- bodymass - FFM

  logtee <- dplyr::case_when(
    age >= 60 ~  0.092 + 0.736 * log(FFM) + -0.030 * log(FM) +  0.011 * sex + -0.008  * age,
    age >= 20 ~ -1.118 + 0.920 * log(FFM) + -0.032 * log(FM) + -0.002 * sex +  0.000 * age,
    age >=  1 ~ -0.348 + 0.784 * log(FFM) + -0.019 * log(FM) +  0.067 * sex + -0.012 * age,
    age >=  0 ~ -1.122 + 1.025 * log(FFM) +  0.034 * log(FM) + -0.014 * sex +  0.254 * age
  )

  # Convert megajoules/day to kcals/day
  exp(logtee) * 239.006 #* 365
}

TEElog2 <- function(bodymass, age, sex){
  # TEE (MJ/day); bodymass (kg); sex (F=0, M=1); age (years)

  if (sex == 0) FFMpct <- 0.79 else FFMpct <- 0.87

  # Don't have FFM or FM for most, so estimate it
  FFM <- bodymass * FFMpct
  FM <- bodymass - FFM

  logtee <- dplyr::case_when(
    age >= 60 ~ -0.773 + 0.797 * log(FFM) + -0.016 * log(FM),
    age >= 20 ~ -1.102 + 0.916 * log(FFM) + -0.030 * log(FM),
    age >=  1 ~ -0.121 + 0.696 * log(FFM) + -0.041 * log(FM),
    age >=  0 ~ -1.270 + 1.163 * log(FFM) +  0.053 * log(FM)
  )

  # Convert megajoules/day to kcals/day
  exp(logtee) * 239.006 #* 365
}

#' @title TEE from body mass, age, and sex
#' @description Total energy expenditure (TEE) from body mass, age, and sex
#' @param bodymass Body mass in kilograms
#' @param age Age in years
#' @param sex 0: female; 1: male
#' @return otal energy expenditure (TEE) in kcals/day
#' @details Formula from Pontzer et al. (2021). Science 373, 808–812, Table S2
#' @examples
#' \dontrun{
#' if(interactive()){
#'  TEE_Ponzter(50, 25, 1) # A 25 yr old male weighing 50 kg
#'  }
#' }
#' @rdname TEE_Ponzter
#' @export
TEE_Ponzter <- function(bodymass, age, sex){
  # TEE (MJ/day); bodymass (kg); sex (F=0, M=1); age (years)

  tee <- dplyr::case_when(
    # sex == 0 & age >= 75 ~ 10.917 + 0.048 * bodymass + 1.659 * sex + -0.080 * age,
    # sex == 0 & age >= 16 ~  5.984 + 0.065 * bodymass + 2.669 * sex + -0.025 * age,
    # sex == 1 & age >= 60 ~ 10.917 + 0.048 * bodymass + 1.659 * sex + -0.080 * age,
    # sex == 1 & age >= 20 ~  5.984 + 0.065 * bodymass + 2.669 * sex + -0.025 * age,
    age >= 60 ~ 10.917 + 0.048 * bodymass + 1.659 * sex + -0.080 * age,
    age >= 20 ~  5.984 + 0.065 * bodymass + 2.669 * sex + -0.025 * age,
    age >=  1 ~  2.592 + 0.080 * bodymass + 1.436 * sex +  0.183 * age,
    age >=  0 ~  0.255 + 0.205 * bodymass + 0.090 * sex +  0.951 * age
  )

  # Convert megajoules/day to kcals/day
  tee * 239.006
}

# https://pmc.ncbi.nlm.nih.gov/articles/PMC11772230/
# https://www.nature.com/articles/s43016-024-01089-5

#' @title TEE from body mass, height, age, sex, elevation
#' @description Total energy expenditure (TEE) from body mass, height, age, sex, and elevation
#' @param bodymass body mass in kilograms
#' @param height Height in centimeters
#' @param age Age in years
#' @param sex 0: female; 1: male
#' @param elevation elevation of residence in meters, Default: 1000
#' @return Total energy expenditure (TEE) in kcals/day
#' @details Formula from Bajunaid et al. (2025) https://www.nature.com/articles/s43016-024-01089-5
#' @examples
#' \dontrun{
#' if(interactive()){
#'  TEE_Bajunaid(50, 150, 25, 0)
#'  }
#' }
#' @rdname TEE_Bajunaid
#' @export
TEE_Bajunaid <- function(bodymass, height, age, sex, elevation = 1000){
  sex <- ifelse(sex == 0, 1, -1)
  lnTEE <-
    -0.21723930921 +
    0.01939639976 + # African ethnicity for now
    # -0.01554772302 + # Hispanic
    # 0.00623768257 + # Asian
    0.41666419569 * log(bodymass) +
    0.00656496388 * height +
    -0.02054339322 * age +
    0.00033079019 * age^2 +
    0.09126350903 * log(elevation) +
    -0.04091711710 * sex +
    -0.00067594646 * height * log(elevation) +
    -0.00000185178 * age^3 +
    0.00201815477 * age * log(elevation) +
    -0.00002262281 * age^2 * log(elevation) +
    -0.00694699228 * sex * log(elevation)
  239.006 * exp(lnTEE)
}

#' @title TEE from body mass, height, age, sex, elevation
#' @description Total energy expenditure (TEE) from body mass, height, age, sex, and elevation
#' @param bodymass body mass in kilograms
#' @param height Height in centimeters
#' @param age Age in years
#' @param sex 0: female; 1: male
#' @param elevation elevation of residence in meters, Default: 1000
#' @return Total energy expenditure (TEE) in kcals/day
#' @details Formula for infants (age 0-1) from Pontzer et al. (2021); otherwise from Bajunaid et al. (2025) https://www.nature.com/articles/s43016-024-01089-5
#' @examples
#' \dontrun{
#' if(interactive()){
#'  TEE(50, 150, 25, 0)
#'  }
#' }
#' @rdname TEE
#' @export
TEE <- function(bodymass, height, age, sex, elevation = 1000){
  # For infants, use values from Pontzer et al. 2021
  ifelse(age < 1, TEE_Ponzter(bodymass, age, sex), TEE_Bajunaid(bodymass, height, age, sex, elevation))
}

#' @title TEE from age, sex, group, elevation
#' @description Total energy expenditure (TEE) from hunter-gatherer age, sex, group, pregnancy status, and elevation
#' @param age Age in years
#' @param sex 0: female; 1: male
#' @param group character "ache", "hadza", "kung", "avg"
#' @param pregnant 0: No; 1: Yes, Default: 0
#' @param elevation Elevation in meters, Default: 1000
#' @return Total energy expenditure in kcals/day
#' @details Formula for infants (age 0-1) from Pontzer et al. (2021); otherwise from Bajunaid et al. (2025) https://www.nature.com/articles/s43016-024-01089-5; heights and weights from the specified group.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  TEE2(25, 1, "hadza") # TEE of 25 year old Hadza male
#'  }
#' }
#' @rdname TEE2
#' @export
TEE2 <- function(age, sex, group, pregnant = 0, elevation = 1000){
  if (length(sex) == 1 && is.na(sex)){
    wt_f <- hg_weight(age, 0, pregnant, group)
    wt_m <- hg_weight(age, 1, pregnant, group)

    ht_f <- hg_height(age, 0, group)
    ht_m <- hg_height(age, 1, group)
    return( (TEE(wt_f, ht_f, age, 0, elevation = elevation) + TEE(wt_m, ht_m, age, 1, elevation = elevation))/2 )
  }
  wt <- hg_weight(age, sex, pregnant, group)
  ht <- hg_height(age, sex, group)
  TEE(wt, ht, age, sex, elevation = elevation)
}
