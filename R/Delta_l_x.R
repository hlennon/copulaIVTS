Delta_l_x <- function(xxx, para, para_old, order){
            m <- ncol(xxx)
            
            Ij   <- numeric(m)
            Ij_1 <- numeric(m)
            
            for(i in 1:m){
                        Ij_1[i] <-  arma_loglikelihood(para_old,  xxx[, i], ord=order)
                        Ij[i]   <-  arma_loglikelihood(para,      xxx[, i], ord=order)
            }
            
            change <- -log(mean(Ij_1/Ij))
            return(change)
}
# Delta_l_x <-
# function(xxx, para, para_old, order){
#             m <- ncol(xxx)
#             
#             Ij   <- numeric(m)
#             Ij_1 <- numeric(m)
#             
#             for(i in 1:m){
#                         Ij_1[i] <-  -arma_loglikelihood(para_old,  xxx[, i], ord=order)
#                         Ij[i]   <-  -arma_loglikelihood(para,      xxx[, i], ord=order)
#             }
#             ll     <- mean(Ij)
#             change <- ll - mean(Ij_1)
#             
#             return(list=list(change=change, loglik=ll))
# }