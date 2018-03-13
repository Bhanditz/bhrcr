fitQuadraticAndStorage<-function(data_in, ldata_out, start_ind, temp3, R){ 
    ###########################################
	#      Fit a quadratic model	          #
	###########################################
	# Initialise the square of the time data
	data_in2 = data_in^2
	# Fit quadratic model
	glm.quad<-glm(ldata_out ~ data_in + data_in2)
	# Find the AIC value associated with the quadratic model 
	AIC2 = AIC(glm.quad)
	mean1 = mean(ldata_out)
	# Storage for the quadratic model
	R[start_ind,30] = glm.quad$coefficients[[2]]            # Slope coefficient (associated with data_in)
	R[start_ind,32] = glm.quad$coefficients[[3]]            # Slope coefficient  (associated with data_in^2)
	# Isolate the standard errors for the quadratic model 
    SEs = sqrt(diag( vcov(glm.quad)))
	R[start_ind,31] = SEs[[2]]                               # SE of slope (associated with data_in)   
	R[start_ind,33] = SEs[[3]]                               # SE of slope (associated with data_in^2)
	R[start_ind,34] = 1 - sum(residuals(glm.quad)^2)/sum((ldata_out-mean1)^2) # Store R2 for quadratic model
	
	R[temp3, 35] <-fitted.values(glm.quad)          # Store fitted values
	R[temp3, 36] <- residuals(glm.quad)             # Store residual values
	if (length(data_in)>4){
        R[temp3, 26] <-ls.diag(glm.quad)$stud.res      # Store student residual values
        R[temp3,27] <-ls.diag(glm.quad)$std.res        # Store standardised residual values
    }
    
	R[start_ind, 39] <-summary.glm(glm.quad)$coefficients[[2,4]]   # Store p-value associated with slope in quadratic model (associated with data_in)
	R[start_ind, 40] <-summary.glm(glm.quad)$coefficients[[3,4]]   # Store p-value associated with slope in quadratic model (associated with data_in^2)
	R[start_ind,41] = AIC2                                   # Store AIC for quadratic model
	
	return(list(glm.quad = glm.quad, R=R, AIC2 = AIC2))
}
