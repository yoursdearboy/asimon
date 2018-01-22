.onLoad <- function(lib, pkg) {
  library.dynam('asimon', pkg, lib)
}

g <- function(s1, r1, n1, s, m, r, n, p) {
    .Call('g_facade', s1, r1, n1, s, m, r, n, p, PACKAGE = packageName())
}

en <- function(s1, r1, n1, m2, n2, p) {
    .Call('en_facade', s1, r1, n1, m2, n2, p, PACKAGE = packageName())
}

asimon_calculate <- function(n1, n2, a_max, b1_max, b2_max, p0, p1, p2) {
    .Call('asimon_facade', n1, n2, a_max, b1_max, b2_max, p0, p1, p2, PACKAGE = packageName())
}

#' Calculates parameters for adaptive two-stage Simon's design for clinical trials.
#'
#' @param a_max type I error constraint.
#' @param b1_max, b2_max type II error constraints for stages 1 and 2.
#' @param p0 null hypothesis response rate.
#' @param p1, p2 alternative hypothesis rates.
#' @return List of 4 elements for each type of optimality.
#' @examples
#' asimon(0.05, 0.2, 0.1, 0.05, 0.2, 0.25)
#' @export
asimon <- function(a_max, b1_max, b2_max, p0, p1, p2) {
    simon.n1 <- clinfun::ph2simon(p0, p1, a_max, b1_max)
    simon.n2 <- clinfun::ph2simon(p0, p2, a_max, b2_max)
    n1 <- simon.n1$out[which.min(simon.n1$out[,5]),4]
    n2 <- simon.n2$out[which.min(simon.n2$out[,5]),4]
    x <- asimon_calculate(n1, n2, a_max, b1_max, b2_max, p0, p1, p2)
    for (i in 1:4) {
        colnames(x[[i]]) <- c('s1', 'r1', 'n1', 's', 'm', 'r', 'n', 'a', 'b1', 'b2', 'EN0', 'EN1', 'EN2')
    }
    class(x) <- 'asimon'
    x
}

print.asimon <- function(x, ...) {
    out <- do.call(rbind, lapply(x, function(r) r[1,]))
    rownames(out) <- paste('OT', 1:nrow(out), sep = '')
    print(out, ...)
}

summary.asimon <- function(x, ...) {
    for (i in 1:4) {
        cat(sprintf("Optimality Type %d\n", i))
        prmatrix(x[[i]], rowlab = rep('',nrow(x[[i]])), ...)
    }
}
