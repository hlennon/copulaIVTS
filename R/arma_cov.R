arma_cov <-
function(a, b, n){
            p <- length(a)
            q <- length(b)
            c <- ARMAtoMA(a, b, q)
            
            #find coeff matrix of R(0), R(1),...,R(p)
            aa<- c(1, -a)
            A <- matrix(0, p+1, p+1)
            for(i in 0:p){
                        for(j in 0:p){
                                    k <- abs(i-j);
                                    A[i+1, k+1] <- A[i+1, k+1] + aa[j+1] }}              
            
            #find the column vector d: Ar=d
            b <- c(1, b)
            c <- c(1, c)
            d <- NULL
            for(i in 0:p){
                        d <- rbind(d, sum(b*c))
                        b    <- c(b[-1], 0)
            }
            r <- solve(A,d)                             #solve Ar=d for r
            
            #find r(p+1),...,r(n) if n > p
            for(j in (p+1):n){
                        dum <- 0
                        for(i in 1:p){
                                    dum  <- dum + a[i]*r[j-i+1]
                        }
                        r <- c(r, dum + sum(b*c))
                        b <- c(b[-1], 0)
            }
            return(r[1:(n+1)])
}
