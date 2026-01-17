
#' @title Hadza male productivity
#' @description Hadza male productivity (kcals/day) digitized from plot in Kaplan et al. 2000
#' @param age Age in years
#' @return kcals
#' @details Interpolation of digitized from plot in Kaplan et al. 2000
#' @examples
#' \dontrun{
#' if(interactive()){
#'  kaplan2000_hadza_male(5:60)
#'  }
#' }
#' @rdname kaplan2000_hadza_male
#' @export
kaplan2000_hadza_male <- function(age){
  approxfun(kaplan2000male$age, kaplan2000male$kcals)(age)
}

#' @title Hadza female productivity
#' @description Hadza female productivity (kcals/day) digitized from plot in Kaplan et al. 2000
#' @param age Age in years
#' @return kcals
#' @details Interpolation of digitized from plot in Kaplan et al. 2000
#' @examples
#' \dontrun{
#' if(interactive()){
#'  kaplan2000_hadza_female(5:60)
#'  }
#' }
#' @rdname kaplan2000_hadza_female
#' @export
kaplan2000_hadza_female <- function(age){
  approxfun(kaplan2000female$age[kaplan2000female$group == 'hadza'], kaplan2000female$kcals[kaplan2000female$group == 'hadza'])(age)
}

#' @title Ache female productivity
#' @description Ache female productivity (kcals/day) digitized from plot in Kaplan et al. 2000
#' @param age Age in years
#' @return kcals
#' @details Interpolation of digitized from plot in Kaplan et al. 2000
#' @examples
#' \dontrun{
#' if(interactive()){
#'  kaplan2000_ache_female(5:60)
#'  }
#' }
#' @rdname kaplan2000_ache_female
#' @export
kaplan2000_ache_female <- function(age){
  approxfun(kaplan2000female$age[kaplan2000female$group == 'ache'], kaplan2000female$kcals[kaplan2000female$group == 'ache'])(age)
}
