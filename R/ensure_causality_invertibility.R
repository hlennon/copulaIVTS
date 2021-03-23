ensure_causality_invertibility <- function(parameters){
                    L   <- length(parameters)
                    ps  <- polynomial(c(-parameters[L:1],1)) 			
                    ps  <- solve(ps)
                    
                    if (all(abs(ps)<1)) return(parameters)
                    else{
                    j   <- numeric(L) 
                    k   <- 0
                    vec <- numeric(L)
                    # We want roots less than 1
                    ps[abs(ps)>=1] <- 1/ps[abs(ps)>=1]			    	        
                    while(k < L){		
                                        j 	 <- add_1(j)	 
                                        k 	 <- sum(j)		
                                        vec[k] <- vec[k] + (-1)^(k+1)*prod(ps^j)
                    }
                    #cat("Stationarity Test: SWITCH", '\n');
                    return(as.numeric(vec))	
}
}
