fitCubicAndStorage<-function(data_in, ldata_out, start_ind, temp3, R){ 
	#################################
	#      Fit a cubic model		#
	#################################
	# Initialise the quadratic of the time data
	data_in2 = data_in^2
	# Initialise the cubic of the time data
	data_in3 = data_in^3
	# Fit cubic model
	glm.cubic<-glm(ldata_out ~ data_in + data_in2 + data_in3)
	# Find the AIC value associated with the cubic model 
	AIC3 = AIC(glm.cubic)
	mean1 = mean(ldata_out)
	# Storage for the cubic model
	R[start_ind,6] = glm.cubic$coefficients[[2]]            # Slope coefficient (associated with data_in)
	R[start_ind,8] = glm.cubic$coefficients[[3]]            # Slope coefficient (associated with data_in^2)
	R[start_ind,10] = glm.cubic$coefficients[[4]]           # Slope coefficient (associated with data_in^3)
	# Isolate the standard errors for the cubic model 
 	SEs = sqrt(diag( vcov(glm.cubic)))
	R[start_ind,7] = SEs[[2]]	                              # SE of slope (associated with data_in)  
	R[start_ind,9] = SEs[[3]]                               # SE of slope (associated with data_in^2)  
	R[start_ind,11] = SEs[[4]]                              # SE of slope (associated with data_in^3)  
	R[start_ind,12] = 1 - sum(residuals(glm.cubic)^2)/sum((ldata_out-mean1)^2)   # Store R2 for cubic model
		
	R[temp3, 13] <-fitted.values(glm.cubic)        # Store fitted values
	R[temp3, 14] <- residuals(glm.cubic)           # Store residual values 
    if (length(data_in)>5){
        R[temp3, 26] <-ls.diag(glm.cubic)$stud.res      # Store student residual values
        R[temp3,27] <-ls.diag(glm.cubic)$std.res        # Store standardised residual values
    }
    
	R[start_ind, 17] <-summary.glm(glm.cubic)$coefficients[[2,4]]    # Store p-value associated with slope in cubic model (associated with data_in)
	R[start_ind, 18] <-summary.glm(glm.cubic)$coefficients[[3,4]]    # Store p-value associated with slope in cubic model (associated with data_in^2)
	R[start_ind, 19] <-summary.glm(glm.cubic)$coefficients[[4,4]]    # Store p-value associated with slope in cubic model (associated with data_in^3)
	R[start_ind,20] = AIC3                                   # Store AIC for cubic model
	
	return(list(glm.cubic = glm.cubic, R=R, AIC3 = AIC3))
}
