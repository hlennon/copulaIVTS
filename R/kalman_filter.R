kalman_filter <-
function(F, G, Q, R, Y, X0, P0){
            
            n  <- length(Y)
            X  <- NULL; 
            E  <- NULL; 
            S  <- NULL
            X1 <- F%*%X0                           	# E[X_t|Y_1,...,Y_t-1]
            P1 <- F%*%P0%*%t(F) +Q	 		# Var{X_t|Y_1,...,Y_t-1}
            E1 <- Y[1] - (G%*%X1)	 		# Y_t - E[Y_t|Y_1,...,Y_t-1]
            
            P1G<- P1%*%t(G)
            SE <- G%*%P1G + R                               # Var{Y_t - E[Y_t|Y_1,...,Y_t-1]}
            K  <- P1G%*%solve(SE)                           # Gain
            

            X0 <- X1 + K%*%E1			# E[X_t|Y_1,...,Y_t]
            P0 <- (P1 - K%*%G%*%P1)			# Var{X_t|Y_1...Y_t]
            X  <- cbind(X, X0)
            E  <- cbind(E, E1)
            S  <- cbind(S, SE)
            for (t in 2:n){
                        X1 <- F%*%X0		# E[X_t|Y_1,...,Y_t-1]
                        P1 <- F%*%P0%*%t(F) +Q	 	# Var{X_t|Y_1,...,Y_t-1}
                        E1 <- Y[t] - (G%*%X1)	 	# Y_t - E[Y_t|Y_1,...,Y_t-1]
                        SE <- G%*%P1%*%t(G) + R	 	# Var{Y_t - E[Y_t|Y_1,...,Y_t-1]}
                        K  <- P1%*%t(G)%*%solve(SE)	# Gain
                        X0 <- X1 + K%*%E1		# E[X_t|Y_1,...,Y_t]
                        P0 <- (P1 - K%*%G%*%P1)		# Var{X_t|Y_1...Y_t]
                        X  <- cbind(X, X0)
                        E  <- cbind(E, E1)
                        S  <- cbind(S, SE)
            }
            return(list(X=X, E=E, S=S))
}
kalman_filter <-
function(F, G, Q, R, Y, X0, P0){
            
            n  <- length(Y)
            X  <- NULL; 
            E  <- NULL; 
            S  <- NULL
            X1 <- F%*%X0                           	# E[X_t|Y_1,...,Y_t-1]
            P1 <- F%*%P0%*%t(F) +Q	 		# Var{X_t|Y_1,...,Y_t-1}
            E1 <- Y[1] - (G%*%X1)	 		# Y_t - E[Y_t|Y_1,...,Y_t-1]
            
            P1G<- P1%*%t(G)
            SE <- G%*%P1G + R                               # Var{Y_t - E[Y_t|Y_1,...,Y_t-1]}
            K  <- P1G%*%solve(SE)                           # Gain
            

            X0 <- X1 + K%*%E1			# E[X_t|Y_1,...,Y_t]
            P0 <- (P1 - K%*%G%*%P1)			# Var{X_t|Y_1...Y_t]
            X  <- cbind(X, X0)
            E  <- cbind(E, E1)
            S  <- cbind(S, SE)
            for (t in 2:n){
                        X1 <- F%*%X0		# E[X_t|Y_1,...,Y_t-1]
                        P1 <- F%*%P0%*%t(F) +Q	 	# Var{X_t|Y_1,...,Y_t-1}
                        E1 <- Y[t] - (G%*%X1)	 	# Y_t - E[Y_t|Y_1,...,Y_t-1]
                        SE <- G%*%P1%*%t(G) + R	 	# Var{Y_t - E[Y_t|Y_1,...,Y_t-1]}
                        K  <- P1%*%t(G)%*%solve(SE)	# Gain
                        X0 <- X1 + K%*%E1		# E[X_t|Y_1,...,Y_t]
                        P0 <- (P1 - K%*%G%*%P1)		# Var{X_t|Y_1...Y_t]
                        X  <- cbind(X, X0)
                        E  <- cbind(E, E1)
                        S  <- cbind(S, SE)
            }
            return(list(X=X, E=E, S=S))
}
