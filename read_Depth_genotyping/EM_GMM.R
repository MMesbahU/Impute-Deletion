
# ----------------------------------- Expectation Maximization Steps--------------------------
underFlow <- 1e-323
min_v <- 1e-20
low_lim_variance <- 1e-30
# Step 1: Expectation
expectation <- function(Data, Weight,  MU, VAR){
  
  # three equal weights of 1/3 each
  # assigning inputs to variables
  d <- Data  
  W <- Weight 
  M <- MU 
  V <- VAR
  
  # Estimating the expected copy-number class weight 
  E_Z_CN0 <-  (  W[1] * dnorm(d, M[1], sqrt(V[1]))) / 
              (( W[1] * dnorm(d, M[1], sqrt(V[1]))) + 
              (  W[2] * dnorm(d, M[2], sqrt(V[2]))) + 
              (  W[3] * dnorm(d, M[3], sqrt(V[3]))))
  
  E_Z_CN1 <-  (  W[2] * dnorm(d, M[2], sqrt(V[2]))) / 
              (( W[1] * dnorm(d, M[1], sqrt(V[1]))) + 
              (  W[2] * dnorm(d, M[2], sqrt(V[2]))) + 
              (  W[3] * dnorm(d, M[3], sqrt(V[3]))))
  
  E_Z_CN2 <-  (  W[3] * dnorm(d, M[3], sqrt(V[3]))) / 
              (( W[1] * dnorm(d, M[1], sqrt(V[1]))) + 
              (  W[2] * dnorm(d, M[2], sqrt(V[2]))) + 
              (  W[3] * dnorm(d, M[3], sqrt(V[3]))))
  
  # Estimating log likelihood using current parameters
  ll_cur <-  sum(log 
              ( pmax( ( W[1] * dnorm(d, M[1], sqrt(V[1]))) + 
              (  W[2] * dnorm(d, M[2], sqrt(V[2]))) +  
              (  W[3] * dnorm(d, M[3], sqrt(V[3]))), underFlow) 
               ), na.rm = TRUE)
  
  #Output of Expectation step  
  return( list(E_Z_CN0, 
               E_Z_CN1, 
               E_Z_CN2, 
               M, 
               ll_cur))

  
}


# Step 2: Maximization 

maximization <- function(Data, W_exp){
  
  min_v <- min_v
  d      <- Data 
  W_exp  <- W_exp 
  M_cur  <- W_exp[[4]]
  ll_cur <- W_exp[[5]]
  # Update Weights
  W_cur  <- c( mean(W_exp[[1]]), 
               mean(W_exp[[2]]), 
               mean(W_exp[[3]])
			  )
  
  ## No update for Mean, using same mean
  # For deletion genotyping, we kept expected read depth fixed and 
  # used genome-wide RD average matched with GC-content (as follows): 
	# (i) Calculate GC% from Reference Genome for the deletion locus 
	# (ii) Calculate Read Depth in bins of 100bp with similar GC-content after removing copy-number variable regions, 
		# such as, known CNVs and SVs from DGVa, Sex-chromosomes, unplaced contigs, repeats, segmental duplicates, etc. 
	# (iii) Used this average RD as expected RD (instead of global average).  
  # Update Means
  # M_cur  <- c( sum(W_exp[[1]]*d)/sum(W_exp[[1]]), 
   #          sum(W_exp[[2]]*d)/sum(W_exp[[2]]), 
    #          sum(W_exp[[3]]*d)/sum(W_exp[[3]]))
  
  
  # Update Variances 
  V_cur <- c( sum(W_exp[[1]]*(d - M_cur[1])^2, na.rm = TRUE)/sum(W_exp[[1]], na.rm = TRUE), 
              sum(W_exp[[2]]*(d - M_cur[2])^2, na.rm = TRUE)/sum(W_exp[[2]], na.rm = TRUE), 
              sum(W_exp[[3]]*(d - M_cur[3])^2, na.rm = TRUE)/sum(W_exp[[3]], na.rm = TRUE) 
			)
  # Checking variances not reaching to Zero, due to under flow 
  V_cur[1] <- if ( is.na(V_cur[1]) || (V_cur[1] < min_v) ) min_v else V_cur[1]
  V_cur[2] <- if ( is.na(V_cur[2]) || (V_cur[2] < min_v) ) min_v else V_cur[2]
  V_cur[3] <- if ( is.na(V_cur[3]) || (V_cur[3] < min_v) ) min_v else V_cur[3]

  # Log likelihood of data with new parameters
  
  ll_new <- sum( log( pmax( ( (W_cur[1]*dnorm(d,M_cur[1],sqrt(V_cur[1]))) +
                      (W_cur[2]*dnorm(d,M_cur[2],sqrt(V_cur[2]))) + 
                      (W_cur[3]*dnorm(d,M_cur[3],sqrt(V_cur[3])))
                      ), underFlow)), na.rm = TRUE) 
  
  
  return( list(W_cur, 
               M_cur, 
               V_cur,
               ll_cur, 
               ll_new))
  
}

# Step 3: Optimization

EM_geno <- function(Data, W_init, MU_init, VAR_init, maxiter = 1500, tol = 1e-6){
  Data  <- Data 
  W_cur <- W_init 
  M_cur <- MU_init
  V_cur <- VAR_init
  
  # Iterate between E and M-step
  i <- 1
  while(i <= maxiter){
    Data  <- Data 
    W_cur <- W_cur 
    M_cur <- M_cur
    V_cur <- V_cur
    
    #### New parameters
    new_Param <- maximization(Data=Data, W_exp=expectation(Data = Data,Weight = W_cur, MU = M_cur, VAR = V_cur))
    
    W_new  <- new_Param[[1]] 
    M_new  <- new_Param[[2]]
    V_new  <- new_Param[[3]]
    ll_cur <- new_Param[[4]] 
    ll_new <- new_Param[[5]]
    
    ## Checking convergence
    if( (abs(ll_new - ll_cur) < tol) || (V_new[1] <= low_lim_variance) || (V_new[2] <= low_lim_variance) || (V_new[3] <= low_lim_variance)) break
    
    # otherwise continue updating 
    cat(paste0(i," ",new_Param[[4]],"\n"))

    W_cur <- W_new 
    M_cur <- M_new 
    V_cur <- V_new
    i =  i+1
    
  }
    
  return(c(new_Param, i))
  
}
################################ EM ENDs ######################
