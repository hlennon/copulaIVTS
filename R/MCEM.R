#' @export
MCEM <- function(obs.data, Monte_Carlo_samples, p, q, initial,  num_iterations,
                prec=0.01,  marginal_dist="negbin",
                optim_method="nmkb", compute_stderrors=TRUE){



            # A MCEM algorithm to computes the MLEs of the ARMA(p,q) parameters
            #              in the Gaussian copula model.
            #
            # Args:
            #   obs.data: Observed integer-valued time series.
            #   Monte_Carlo_samples: The number of Monte Carlo samples to be used in E-Step
            #         with no missing values.
            #   p, q: The order of the underlying ARMA model.
            #   initial: The initial values to begin the MCEM algorithm.
            #   num_iterations: The length of the output of the MCEM sequence.
            #   prec: The numeric value to determine the stopping criteria,
            #         default to 0.01.
            #   marginal_dist: The specification of the marginal distribution,
            #         default="negbin".
            #   optim_method: The numerical optimisation proceedure to be used,
            #         default="nmkb" in dfoptim package.
            #   compute_stderrors: Option for the output to include the standard
            #         errors of the MLEs, defualt=TRUE.
            #
            # Returns:
            #   The maximum likelihood estimates for the observed data.


                    # The ARMA order
                    p <<- p
                    q <- q

                    # Allocate space/Initial setup
                    mm        <- as.vector(Monte_Carlo_samples)
                    ml        <- length(mm)
                    n         <- length(obs.data)
                    t         <- 1
                    save      <- matrix(0,   nrow=(num_iterations+1)*ml, ncol=(p+q))
                    ll        <- numeric((num_iterations+1)*ml)
                    saved     <- numeric(p+q)
                    para      <- numeric(p+q)
                    diff      <- rep(10, (num_iterations*ml))
                    finish    <- FALSE
                    para.star <- initial
                    MLEs      <- 0


                    # Timestamp for file names
                    da        <- format(Sys.time(), "%a_%b_%d_%X_%Y")
                    file.name <- sprintf("MCEM %s%s results_%s.txt", p,q,da)
                    file.create(file.name)


                    # Specify estimation methods for distribution functions
                    if(marginal_dist=="non-parametric"){
                                        Fn        <- ecdf(obs.data)
                                        up        <- qnorm(Fn(obs.data))
                                        low       <- qnorm(Fn(obs.data-1))
                    }
                    if(marginal_dist=="poisson"){
                                        lambdahat  <- coef(fitdistr(obs.data, "poisson"))[1]
                                        up  <- c(qnorm(ppois(obs.data,      lambda=lambdahat)))
                                        low <- c(qnorm(ppois((obs.data-1),   lambda=lambdahat)))
                                        MLEs <- lambdahat
                    }
                    if(marginal_dist=="negbin"){
                                        nbs  <- coef(fitdistr(obs.data, "negative binomial"))[1]
                                        nbmu <- coef(fitdistr(obs.data, "negative binomial"))[2]
                                        nbp  <- nbs/(nbs+nbmu)
                                        up  <- c(qnorm(pnbinom(obs.data,      size=nbs, prob=nbp)))
                                        low <- c(qnorm(pnbinom((obs.data-1),  size=nbs, prob=nbp)))
                                        MLEs <- c(nbs, nbp)
                    }
                    if(marginal_dist=="zeroinf_poisson"){
                                        me <- mean(obs.data)
                                        s2 <- sum((obs.data-me)^2)/(n-1)
                                        if(me>s2){
                                                            lambdahat  <- me
                                                            pstr       <- 0
                                        }else{      sm         <- s2/me
                                                    lambdahat  <- me + sm - 1
                                                    pstr       <- (sm - 1)/lambdahat
                                        }
                                        up  <- c(qnorm(pzipois(obs.data,       lambda=lambdahat, pstr0 = pstr)))
                                        low <- c(qnorm(pzipois((obs.data-1),   lambda=lambdahat, pstr0 = pstr)))
                                        MLEs <- c(lambdahat, pstr)
                    }

                    # Check the results
                    print(MLEs);


                    # Begin MCEM algorithm iterations
                    for(k in 1:ml){
                                        # Set number of samples for Monte Carlo
                                        m  <- mm[k]
                                        cat('m=', m, '\n')

                                        # The E and the M-step
                                        # Each loop is one iteration
                                        # Different numerical optimisers can be chosen
                                        for(j in 1:num_iterations){
                                                            current <- para.star
                                                            xxx  <- matrix(0, ncol=length(obs.data), nrow=m)

                                                            if(optim_method=="L-BFGS-B"){
                                                                                para <- optim(current, Q_function, para.star=para.star, m=m, low=low, up=up, p=p, q=q, method="L-BFGS-B", control=list(maxit=3000))$par
                                                            }


                                                            if(optim_method=="BFGS"){
                                                                                para <- optim(current, Q_function, para.star=para.star, m=m, low=low, up=up, p=p, q=q, method="BFGS", control=list(maxit=3000))$par
                                                            }


                                                            if(optim_method=="CG"){
                                                                                para <- optim(current, Q_function, para.star=para.star, m=m, low=low, up=up, p=p, q=q, method="CG", control=list(maxit=3000))$par
                                                            }

                                                            if(optim_method=="nlm"){
                                                                                para <- nlm(f=Q_function, p=current, para.star=para.star, m=m, low=low, up=up, p=p, q=q)$estimate
                                                            }


                                                            if(optim_method=="nmkb"){
                                                                                lower = c(rep(-0.95, (p+q)))
                                                                                upper = c(rep( 0.95, (p+q)))
                                                                                para <- nmkb(fn = Q_function, par = current, para.star=para.star, m=m, low=low, up=up, p=p, q=q, lower = lower, upper = upper)$par
                                                            }

                                                            # Save the samples from the Monte Carlo Step
                                                            xdata <- xxx

                                                            # Check: Ensure the stationarity of the updated parameters
                                                            para[1:p]        <-  ensure_causality_invertibility(para[1:p])
                                                            para[(p+1):(p+q)]<- -ensure_causality_invertibility(-para[(p+1):(p+q)])

                                                            # Estimate the change in log-likelihood
                                                            ll[t]            <- Delta_l_x(xdata, para.star, para, c(p,q))
                                                            diff[t]          <- ll[t]

                                                            # Check: Has the stopping criteria been satisfied for the past parameter values?
                                                            if(t>6) finish   <- all(diff[t:(t-5)]<prec)


                                                            # Save and print updated parameter values
                                                            t                <- t + 1
                                                            save[t, ]        <- para.star  <- para
                                                              cat('index = ', t, ' parameters = ', para, ' Change_in_lik = ', diff[t-1], ' finish = ', finish, "\n");
                                                            line =  c("Iteration: ", t, " gives ", para, ' Change_in_lik = ', diff[t-1], '\n')

                                                            # Append the updated values to the results file
                                                            write(line, file=file.name, append=TRUE)
                                        }
                    }

                    # Computing the Std errs via numerical optimisation
                    if(compute_stderrors==TRUE) se <- compute_std_errors(para, obs.data,  m, p, q, low, up)

                    # Collect results for output
                    ret_list <- list(MCEM_Iterations=save[1:t,], diff=diff[1:t], dall=dall[1:t], Monte_Carlo_samples=Monte_Carlo_samples, t=t, MLEs=MLEs)

                    return(ret_list)
}
