#' Create a list of options to be used in optim()
#' @export
optimctrl <- function() {
    list(maxit = 5000, reltol = 1e-15, fnscale = -1)
}
