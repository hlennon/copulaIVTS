#' @export
log_lik_arma <-function(xxx, para, order){
                    m <- ncol(xxx)
                    Ij   <- numeric(m)

                    for(i in 1:m){
                    Ij[i]   <-  exp(arma_loglikelihood(para,      xxx[, i], ord=order))
                    }

return(Ij)
}
