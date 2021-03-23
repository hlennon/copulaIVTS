#' @export
direct_GCNB_AR1_likelihood <-
function(rho, nb.p, nb.s, data, method="mnormt") {
                    n <- length(data)
                    u <- matrix(0, n, 2)
                    for (i in 1:n){
                                        u[i,1] <- qnorm(pnbinom(data[i],   size=nb.s, prob=nb.p))
                                        u[i,2] <- qnorm(pnbinom(data[i]-1, size=nb.s, prob=nb.p))
                    }
                    sigma <- toeplitz(rho^(0:(n-1)))
                    ell <- 0
                    j <- rep(0, n)


                    if(method=="mnormt"){

                                        while (sum(j) < n) {
                                                            #cat(j, '\n')
                                                            ell <- ell + (-1)^(sum(j))*mnormt::pmnorm(u[(1:n) + j*n], mean=rep(0, n), varcov=sigma)[1]
                                                            j <- add_1(j)
                                        }
                                        ell <- ell + (-1)^(sum(j))*mnormt::pmnorm(u[(1:n) + j*n], mean=rep(0, n), varcov=sigma)[1]
                    }

                    if(method=="mvtnorm"){

                                        while (sum(j) < n) {
                                                            #cat(j, '\n')
                                                            ell <- ell + (-1)^(sum(j))*mvtnorm::pmvnorm(upper=u[(1:n) + j*n], mean=rep(0, n), sigma=sigma, abseps=1e-26, algorithm="MIWA")[1]
                                                            j <- add_1(j)
                                        }
                                        ell <- ell + (-1)^(sum(j))*mvtnorm::pmvnorm(upper=u[(1:n) + j*n], mean=rep(0, n), sigma=sigma, abseps=1e-26, algorithm="MIWA")[1]

                    }




                    return(ell)
}
