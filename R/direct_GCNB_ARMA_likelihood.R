#' @export
direct_GCNB_ARMA_likelihood <-
function(covariances, nb.p, nb.s, data, method="mvtnorm") {
                    n <- length(data)
                    u <- matrix(0, n, 2)
                    u <- qnorm(pnbinom(cbind(data, data-1), size=nb.s, prob=nb.p))

                    sigma <- toeplitz(covariances)
                    ell   <- numeric(2^n)
                    j <- rep(0, n)
                    index <- 1

                    if(method=="mnormt"){
                                        while (sum(j) < n) {
                                                            #cat(j, '\n')
                                                            ell[index] <- (-1)^(sum(j))*mnormt::pmnorm(u[(1:n) + j*n], mean=rep(0, n), varcov=sigma)[1]
                                                            index <- index + 1
                                                            j <- add_1(j)


                                        }
                                        ell[index] <-  (-1)^(sum(j))*mnormt::pmnorm(u[(1:n) + j*n], mean=rep(0, n), varcov=sigma)[1]
                    }

                    if(method=="mvtnorm"){
                                        while (sum(j) < n) {
                                                            #cat(j, '\n')
                                                            ell[index] <- (-1)^(sum(j))*mvtnorm::pmvnorm(upper=u[(1:n) + j*n], mean=rep(0, n), sigma=sigma,  abseps=1e-26, algorithm="MIWA")[1]
                                                            index <- index + 1
                                                            j <- add_1(j)
                                        }

                                        ell [index] <-  (-1)^(sum(j))*mvtnorm::pmvnorm(upper=u[(1:n) + j*n], mean=rep(0, n), sigma=sigma,  abseps=1e-26,algorithm="MIWA")[1]
                    }

                    return(ell)
}
