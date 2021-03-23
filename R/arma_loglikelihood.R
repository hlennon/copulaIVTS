#' @export
arma_loglikelihood <-
function(para, data, ord){

            p <- ord[1]
            q <- ord[2]
            #           para[1:p] <- Convert(para[1:p])                  # Ensure new parameter values give a stationary Yt time series
            #           para[(p+1):(p+q)]<- -Convert(-para[(p+1):(p+q)]) # Ensure new parameter values give a causal Yt time series

            a <- para[1:p]
            b <- para[(p+1):(p+q)]
            x <- data
            n <- length(x)
            if (p < q + 1){                                              # Make sure p = q + 1
                        a <- c(a, array(0,c(q + 1 - p)))                 # pad a with zeros
                        p <- q + 1
            }                                     	             # increase p to q+1
            if (p > q + 1){
                        b <- c(b, array(0,c(p - 1 - q))) 	             # pad b with zeros
                        q = p - 1;                	    	 # increase q to p-1
            }
            r      <- ARMAacf(ar = a, ma = b, pacf = FALSE)
            sigma2 <- 1/r[1]


            # Set up Kalman Filter to compute e_t's and tau_t's (aka S)
            F  <- rbind(a, cbind(diag(1,(p-1)), c(rep(0,(p-1)))))
            G  <- rbind(c(1, b))
            Q  <- diag(c(1, array(0,c(p-1, 1))))
            R  <- 0
            X0 <- array(0, c(p,1))
            P0 <- toeplitz(r[1:p])
            kal<- kalman_filter(F, G, Q, R, x, X0, P0)


            # Evaluates Reduced Gaussian log-likelihood
            g   <- n*log(2*pi) + n*log(sigma2)	+ sum(log(kal$S)) + sum(kal$E^2/kal$S)/sigma2
            g   <- as.numeric(-g/2)
            return(g)
}
