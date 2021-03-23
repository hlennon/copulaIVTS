LD <-
function(r, n){                      
            s         <- vector(length=n)
            phi.mat   <- diag(1, nrow=n+1, ncol=n+1);
            phi       <- r[2]/r[1]
            s[1]      <- r[1]*(1-phi^2) 
            phi.mat[1,2]<- -phi
            for(k in 2:n){
                        phi_k <- (r[k + 1] - sum(phi*r[k:2]))/s[k-1] 
                        phi   <- c(phi - phi_k*phi[(k-1):1], phi_k)
                        phi.mat[k:1, (k+1)]<- -phi
                        s[k]  <- s[k-1]*(1-phi_k^2)
            }
            s  <- c(r[1], s)
            return(cbind(phi.mat, s))
}
