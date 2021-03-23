MCEM_POLIO <-
function(obs.data, Monte_Carlo_samples, p, q, initial, marginal_dist, num_iterations,
                 prec=0.00001, optim_method="nmkb", compute_stderrors=TRUE){
                    if(compute_stderrors==TRUE) require(numDeriv)
                    p <<- p
                    q <- q
                    mm        <- as.vector(Monte_Carlo_samples)
                    ml        <- length(mm)
                    n         <- length(obs.data)
                    t         <- 1
                    pq        <- p+q
                    save      <- matrix(0,   nrow=(num_iterations+1)*ml, ncol=(p+q))
                    saved     <- numeric(p+q)
                    para      <- numeric(p+q)
                    ll        <- numeric((num_iterations+1)*ml)
                    diff      <- c(rep(10, (num_iterations*ml)))
                    finish    <- FALSE
                    para.star <- initial

                    da        <- format(Sys.time(), "%a_%b_%d_%X_%Y")
                    file.name <- sprintf("Polio_MCEM %s%s results_%s %s.txt", p,q,marginal_dist,da)
                    file.create(file.name)

## ---------------------------------------------------------------------
## This Section is specific to the Polio Dataset
                    time      <- (c(1:n)/1000) - 0.073
                    x         <- 2*pi*(0:(n-1))/6
                    sina      <- sin(x)
                    sinb      <- sin(x/2)
                    cosa      <- cos(x)
                    cosb      <- cos(x/2)

                    sim_data  <- data.frame(obs.data, time, cosb, sinb, cosa, sina)
                    names(sim_data)= c("y", "time", "cos12", "sin12",  "cos6", "sin6")

     if(marginal_dist=="semi"){
                    Fn        <- ecdf(obs.data)
                    up        <- qnorm(Fn(obs.data))
                    low       <- qnorm(Fn(obs.data-1))
                    MLEs      <- NULL
                    }

     if(marginal_dist=="poisson"){
                    # Estimate the beta parameters of xT beta
                    beta_hat <- coef(glm(y~., family=poisson, data=sim_data))


                    # Create design matrix X
                    x          <- data.frame(rep(1,n), sim_data[,(-1)])
                    x          <- as.matrix(x)
                    x_beta     <- x%*%beta_hat

                    # Estimate lambda for poisson marginal regression model (link function=log)
                    lambda_hat <-  exp(x_beta)

                    # Estimate the upper and lower truncatiojn points of the underlying arma
                    up  <- qnorm(ppois(obs.data, lambda=lambda_hat))
                    low <- qnorm(ppois(obs.data-1, lambda=lambda_hat))

                    MLEs <- list(beta_hat=beta_hat)
                    print(MLEs);

}
if(marginal_dist=="negbin"){
                    # Estimate the beta parameters of xT beta
                    fit <- MASS::glm.nb(y~., data=sim_data, link=log)

                    beta_hat <- fit$coefficients
                    theta    <- fit$theta

                    # Create design matrix X
                    x          <- data.frame(rep(1,n), sim_data[,(-1)])
                    x          <- as.matrix(x)
                    x_beta     <- x%*%beta_hat
                    kappa      <- 1/theta
                    # Estimate lambda for poisson marginal regression model (link function=log)
                    lambda_hat <-  exp(x_beta)

                    # Estimate the upper and lower truncatiojn points of the underlying arma

                    up   <- c(qnorm(pnbinom(obs.data,      size=theta, mu=lambda_hat)))
                    low  <- c(qnorm(pnbinom((obs.data-1),  size=theta, mu=lambda_hat)))
                    MLEs <- list(beta_hat=beta_hat, theta=theta, kappa=kappa)
                    print(MLEs);
}


## ---------------------------------------------------------------------
for(k in 1:ml){
                    m  <- mm[k]
                    cat('m=', m, '\n')


for(j in 1:num_iterations){
                    current <- para.star
                    xxx  <- matrix(0, ncol=length(obs.data), nrow=m)

                    if(optim_method=="L-BFGS-B"){
                                        para <- optim(current, fn=Q_function, para.star=para.star, m=m, low=low, up=up, p=p, q=q, method="L-BFGS-B", control=list(maxit=3000))$par
                                        }


                    if(optim_method=="BFGS"){
                                        para <- optim(current, fn=Q_function, para.star=para.star, m=m, low=low, up=up, p=p, q=q, method="BFGS", control=list(maxit=3000))$par
                                        }


                    if(optim_method=="CG"){
                                        para <- optim(current, fn=Q_function, para.star=para.star, m=m, low=low, up=up, p=p, q=q, method="CG", control=list(maxit=3000))$par
                                        }

                    if(optim_method=="nlm"){
                                        para <- nlm(f=Q_function, p=current, para.star=para.star, m=m, low=low, up=up, p=p, q=q)$estimate
                                        }


                    if(optim_method=="nmkb"){
                                        if(pq!=1){

                                        lower_bound <- c(rep(-0.95, (p+q)))
                                        upper_bound <- c(rep( 0.95, (p+q)))
                                        para <- dfoptim::nmkb(fn = Q_function, par = current, para.star=para.star, m=m, low=low, up=up, p=p, q=q, lower =lower_bound, upper = upper_bound)$par
                                        }
                                        if(pq==1){
                                                  para <- optimize(f=Q_function,  lower =0.95, upper=0.95, para.star=para.star, m=m, low=low, up=up)$par

                                        }
                    }


                    xdata <- xxx
                    para[1:p]        <-  ensure_causality_invertibility(para[1:p])
                    para[(p+1):(p+q)]<- -ensure_causality_invertibility(-para[(p+1):(p+q)])


                    ## ---------------------------------------------------------------------
                    ##  Compute the CHANGE in likelihood using the data simulated using new.para
                    if(t>2) diff[t]  <- Delta_l_x(xdata, para.star, para, c(p,q))

                    ## ---------------------------------------------------------------------
                    ## Stopping Criterion:
                    ##  Have the change in lik of the last 5 iterations not exceeded the precision?
                    if(t>6) finish   <- all(diff[t:(t-5)]<prec)


                    ## ---------------------------------------------------------------------
                    ##  Parameter Updates
                    t                <- t + 1
                    save[t, ]        <- para.star  <- para
                    cat('index = ', t, ' parameters = ', para,  ' Change_in_lik = ', diff[t-1], ' finish = ', finish, "\n");
                    line =  c("Iteration: ", t, " gives ", para, ' Change_in_lik = ', diff[t-1], '\n')
                    write(line, file=file.name, append=TRUE)


                    }
                    }

                    ## ---------------------------------------------------------------------

                    if(compute_stderrors==TRUE) { print("Computing Standard Errors for parameter estimates: Please be patient.");
                                        std_err <- compute_std_errors(para, obs.data,  m, p, q, low, up)
                    }
                    ret_list <- list(MCEM_Iterations=save[1:t,], change_in_liklhd=diff[1:(t-1)],  Monte_Carlo_samples=Monte_Carlo_samples, t=t, MLEs=MLEs, Std_Errors=std_err)
return(ret_list)
}
