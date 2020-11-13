#' @title Server side function to calculates the survival object
#' @description This function is similar to R function \code{Surv}.
#' @details This function calculates the survival object
#' usually used as a response variable in a model formula..
#' @param x a numeric vector indicating the follow up time.
#' @param y a numeric vector with compatible dimensions to x, providing the name
#' of the status indicator(or event).
#' @return the results of the survival object stored on the server.
#' @author Sofack, Ghislain. (based on corTestDS by Demetris Avraam, for DataSHIELD Development Team)
#' @export
#'

SurvDS <- function (x,y) {

  x.var <- eval(parse(text=x), envir = parent.frame())
  y.var <- eval(parse(text=y), envir = parent.frame())

  # Calling the survival function

  out <- survival::Surv(x.var, y.var)

  # return the Surv object
  return(out)

}

#ASSIGN FUNCTION
# SurvDS
