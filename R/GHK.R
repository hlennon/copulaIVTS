#' @export
GHK <-
function(m, Sigma, a, b){
        N      <- length(a)
        ut     <- matrix(0, N, m)
        u      <- matrix(0, N, m)

        pc     <- pnorm(a[1])
        pd     <- pnorm(b[1])
        ut[1,] <- runif(m)*(pd - pc) + pc
        u[1, ] <- qnorm(ut[1, ])

        L  <- t(chol(Sigma))
        for(i in 2:N){
                dum      <- L[i, ] %*% u
                pc       <- pnorm((a[i] - dum)/L[i, i])
                pd       <- pnorm((b[i] - dum)/L[i, i])
                ut[i, ]  <- runif(m)*(pd-pc) + pc
                u[i, ]   <- qnorm(ut[i, ])
        }
        return(L%*%u)
}
