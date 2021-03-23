# if(marginals=="nonparametric"){
#                     Fn   <- ecdf(obs.data)                                  
#                     up   <- qnorm(Fn(obs.data))
#                     low  <- qnorm(Fn(obs.data-1))
# }
# 
# if(marginals=="nbinom"){
#                     up  <- qnorm(pnbinom(obs.data, 10, 0.5))
#                     low <- qnorm(pnbinom((obs.data-1), 10, 0.5))
# }
# 

compute_std_errors <-
function(theta, obs.data, m, p, q, low, up){
            n    <- length(obs.data)
            order<- c(p,q)
            a    <- theta[1:p]
            b    <- theta[(p+1):(p+q)]   
            r    <- as.numeric(ARMAacf(ar=a, ma=b, lag.max=n, pacf=FALSE))   
            P    <- toeplitz(r[1:n]) 
            
            g    <- matrix(0, nrow=(p+q), ncol=m)
            gg   <- NULL
            h    <- NULL       
            for(j in 1:m){
                        sim.data    <- GHK(1, P, low, up)
                        g[ ,j]      <- grad(func=arma_loglikelihood, theta, data=sim.data, ord=order)
                        h[[j]]      <- -hessian(func=arma_loglikelihood, theta, data=sim.data, ord=order)
            }
            
            for (i in 1:m){     gg[[i]] <- g[ ,i]%*%t(g[ ,i])        }
            
            
            D1  <- Reduce('+', h )/m
            D2  <- Reduce('+', gg)/m
            
            I   <- D1 - D2
            Iinv<- solve(I)
            
            se  <- sqrt(diag(Iinv))
            return(se)
}
# if(marginals=="nonparametric"){
#                     Fn   <- ecdf(obs.data)                                  
#                     up   <- qnorm(Fn(obs.data))
#                     low  <- qnorm(Fn(obs.data-1))
# }
# 
# if(marginals=="nbinom"){
#                     up  <- qnorm(pnbinom(obs.data, 10, 0.5))
#                     low <- qnorm(pnbinom((obs.data-1), 10, 0.5))
# }
# 

compute_std_errors <-
function(theta, obs.data, m, p, q, low, up){
            n    <- length(obs.data)
            order<- c(p,q)
            a    <- theta[1:p]
            b    <- theta[(p+1):(p+q)]   
            r    <- as.numeric(ARMAacf(ar=a, ma=b, lag.max=n, pacf=FALSE))   
            P    <- toeplitz(r[1:n]) 
            
            g    <- matrix(0, nrow=(p+q), ncol=m)
            gg   <- NULL
            h    <- NULL       
            for(j in 1:m){
                        sim.data    <- GHK(1, P, low, up)
                        g[ ,j]      <- grad(func=arma_loglikelihood, theta, data=sim.data, ord=order)
                        h[[j]]      <- -hessian(func=arma_loglikelihood, theta, data=sim.data, ord=order)
            }
            
            for (i in 1:m){     gg[[i]] <- g[ ,i]%*%t(g[ ,i])        }
            
            
            D1  <- Reduce('+', h )/m
            D2  <- Reduce('+', gg)/m
            
            I   <- D1 - D2
            Iinv<- solve(I)
            
            se  <- sqrt(diag(Iinv))
            return(se)
}
