fitLinearAndStorage<-function(data_in, ldata_out, start_ind, temp3, R){ 
	###########################################
	#      Fit a linear model               	#
	###########################################
	glm.linear<-glm(ldata_out ~ data_in) 
    # Find the AIC value associated with the linear model    
	AIC1 = AIC(glm.linear)

	# Storage for the linear model
	R[start_ind,21] = glm.linear$coefficients[[2]]     # Slope coefficient 
	R[start_ind,71] = glm.linear$coefficients[[1]]     # intercept coefficient - for output to user later on

    # Isolate the standard errors for the linear model           
	SEs = sqrt(diag( vcov(glm.linear)))                         
	R[start_ind,22] = SEs[[2]]                         # SE of slope    
    # Find the mean of the data for calculation of the R2 value of linear model              
	mean1 = mean(ldata_out)
	R[start_ind,23] = 1 - sum(residuals(glm.linear)^2)/sum((ldata_out-mean1)^2)    # R2 value for linear model
	
	R[temp3, 24] <-fitted.values(glm.linear)         # Store fitted values
    R[temp3, 25] <- residuals(glm.linear)            # Store residual values
    if (length(data_in)>3){
        R[temp3, 26] <-ls.diag(glm.linear)$stud.res      # Store student residual values
        R[temp3,27] <-ls.diag(glm.linear)$std.res        # Store standardised residual values
    }
	R[start_ind, 28] <-summary.glm(glm.linear)$coefficients[[2,4]]    # Store p-value associated with slope in linear model
	R[start_ind,29] = AIC1                                    # Store AIC for linear model

	
	return(list(glm.linear = glm.linear, R=R, AIC1 = AIC1))
}
