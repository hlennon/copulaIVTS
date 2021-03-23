#' @export
Q_function <-
function(para, para.star, m, low, up, p, q){


            para.star[1:p] <- ensure_causality_invertibility(c(para.star[1:p]))         # Ensures parameters and current pa
            a.star     <- c(para.star[1:p])                      # ar current parameters
            b.star     <- c(para.star[(p+1):(p+q)])              # ma current parameters

            para[1:p]  <- ensure_causality_invertibility(para[1:p])
            a          <- para[1:p]
            b          <- para[(p+1):(p+q)]


            N          <- length(low)

            if(p==0){
                        r.star     <- as.numeric(ARMAacf(ma=b.star, lag.max=N, pacf=FALSE))   #r.star <- arma_cov(a.star, b.star, N); r.star <- r.star/r.star[1]
                        r          <- as.numeric(ARMAacf(ma=b,      lag.max=N, pacf=FALSE))                     # Theoretical autocovariances for parameters
            }else if(q==0){
                        r.star     <- as.numeric(ARMAacf(ar=a.star,  lag.max=N, pacf=FALSE))   #r.star <- arma_cov(a.star, b.star, N); r.star <- r.star/r.star[1]
                        r          <- as.numeric(ARMAacf(ar=a,       lag.max=N, pacf=FALSE))                     # Theoretical autocovariances for parameters
            }
            else{
            r.star     <- as.numeric(ARMAacf(ar=a.star, ma=b.star, lag.max=N, pacf=FALSE))   #r.star <- arma_cov(a.star, b.star, N); r.star <- r.star/r.star[1]
            r          <- as.numeric(ARMAacf(ar=a,      ma=b, lag.max=N, pacf=FALSE))                     # Theoretical autocovariances for parameters
            }



            output     <- LD(r,N-1)
            PHI        <- t(output[, 1:N])
            tau        <- output[,(N+1)]
            P          <- toeplitz(r.star[1:N])                # Covariance matrix for Mult Truncated Normal
            sim.data   <- GHK(m, P, low, up)
            xxx       <<- sim.data
            YY         <- sim.data%*%t(sim.data)/m
            tau.mat    <- diag(1/tau)                           # Matix with tau^-1 in diagonal elements
            A          <- PHI%*%tau.mat%*%t(PHI)                # Matrix form of equations
            like       <- sum(A*YY) + sum(log(tau))
            return(like)
}
