generate_data <-
function (n, a, b, marginal_dist, marginal_parameters){
            
            if(marginal_dist=="poisson"){        lambda  <- marginal_parameters[1]}
            if(marginal_dist=="geometric"){      prob    <- marginal_parameters[1]}
            if(marginal_dist=="zeroinf_poisson"){lambda  <- marginal_parameters[1]}
            if(marginal_dist=="nbinom"){          nb.s   <- marginal_parameters[1]
                                                  nb.p   <- marginal_parameters[2]}
            if(marginal_dist=="zeroinf_nbinom"){  nb.s   <- marginal_parameters[1]
                                                  nb.p   <- marginal_parameters[2]}
            
            
            ar.sigma <- arma_cov(a, b, 0)
            ar.sigma <- sqrt(1/ar.sigma)
            
            if(a[1]==0){  y <- arima.sim(list(ma = b), n, sd = ar.sigma)}
            else if(b[1]==0){  y <- arima.sim(list(ar = a), n, sd = ar.sigma)}
            else{ y <- arima.sim(list(ar = a, ma = b), n, sd = ar.sigma)}
            
            if(marginal_dist=="poisson")  x <- qpois(pnorm(y), lambda)
            if(marginal_dist=="nbinom")   x <- qnbinom(pnorm(y), size = nb.s, prob = nb.p)
            if(marginal_dist=="geometric")  x <- qgeom(pnorm(y), prob)
            if(marginal_dist=="zeroinf_poisson") {
                        x <- qzipois(pnorm(y), lambda, pstr0 = 0)
            }
            if(marginal_dist=="zeroinf_nbinom") {
                        x <- qzinegbin(pnorm(y), size=nb.s, prob=nb.p, pstr0 = 0)
            }
            
            list(x = x, y = y)
}
