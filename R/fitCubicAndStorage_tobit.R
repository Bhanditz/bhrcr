fitCubicAndStorage_tobit<-function(data_in, ldata_out, start_ind, temp3, R, DT, MaxIters){ 
	#################################
	#      Fit a cubic model		#
	#################################
	# Initialise the quadratic of the time data
	data_in2 = data_in^2
	# Initialise the cubic of the time data
	data_in3 = data_in^3
	# Fit cubic model
	
	
	# transform data for fitting purposes:
	mu = mean(data_in)
	sd = sd(data_in)
	data_in_new = (data_in - mu)/sd
	data_in_new2 = data_in_new^2
	data_in_new3 = data_in_new^3
	scaled = 1
	
	ind1 = which(ldata_out == log(DT))
		if (length(ind1) != 0){
	m2 <- lm(ldata_out[-ind1] ~ poly(data_in_new[-ind1], 3, raw = TRUE)) 
		}
	else{	
		m2 <- lm(ldata_out[-length(data_in)] ~ poly(data_in_new[-length(data_in)], 3, raw = TRUE)) 
	}

	glm.cubic <- try(tobit(ldata_out ~ data_in_new + data_in_new2 + data_in_new3, left = log(DT), iter.max = MaxIters*10))
	if (class(glm.cubic)[1]=="try-error"){
		yy<-sprintf("Trying to fix error in tobit fitting (1) - cubic")
		print(yy)
		glm.cubic <-try(tobit(ldata_out ~ data_in_new + data_in_new2 + data_in_new3, left = log(DT), iter.max = MaxIters, init = coef(m2)) ) 
		
		}
	if (class(glm.cubic)[1]=="try-error"){
		yy<-sprintf("Trying to fix error in tobit fitting (2) - cubic")
		print(yy)
		m2 <- lm(ldata_out ~ poly(data_in_new, 3, raw = TRUE))
		glm.cubic <- try(tobit(ldata_out ~ data_in_new + data_in_new2 + data_in_new3, left = log(DT), iter.max = MaxIters, init = coef(m2)))
		}
	if (class(glm.cubic)[1]=="try-error"){
		yy<-sprintf("Trying to fix error in tobit fitting (3a) - cubic")
		print(yy)
		ldata_out[length(ldata_out)] = 0
		m2 <- lm(ldata_out ~ poly(data_in, 3, raw = TRUE))
		scaled = 0
		glm.cubic <- try(tobit(ldata_out ~ poly(data_in, 3, raw = TRUE), left = log(DT), iter.max = MaxIters, init = coef(m2)))
		}
	else if (coef(glm.cubic)[[1]]<0){
		yy<-sprintf("Trying to fix error in tobit fitting (3b) - cubic")
		print(yy)
		ldata_out[length(ldata_out)] = 0
		m2 <- lm(ldata_out ~ poly(data_in, 3, raw = TRUE))
		scaled = 0
		glm.cubic <- try(tobit(ldata_out ~ poly(data_in, 3, raw = TRUE), left = log(DT), iter.max = MaxIters, init = coef(m2)))
		
		}
	if (class(glm.cubic)[1]=="try-error"){
		yy<-sprintf("Trying to fix error in tobit fitting (4) - cubic")
		print(yy)
		glm.cubic <- try(tobit(ldata_out ~ poly(data_in, 3, raw = TRUE), left = log(DT), iter.max = MaxIters*10))
		scaled = 0

		}
	if (class(glm.cubic)[1]=="try-error"){
		yy<-sprintf("NOOOOOO -->Did not fix error in tobit fitting!! - cubic")
		print(yy)
		glm.cubic <- tobit(ldata_out ~ data_in_new + data_in_new2 + data_in_new3, left = log(DT), iter.max = MaxIters)
		}

	# Extract coefficents	
	a = coef(glm.cubic)[[1]]
	b = coef(glm.cubic)[[2]]
	c = coef(glm.cubic)[[3]]
	d = coef(glm.cubic)[[4]]

	# Transform coefficients back to our orignial model (if necessary)
	if (scaled == 1){
	a1 = a - b*mu/sd + c*mu^2/sd^2 - d*mu^3/sd^3
	b1 = b/sd - 2*mu*c/sd^2 + 3*mu^2*d/sd^3
	c1 = c/sd^2 - 3*mu*d/sd^3
	d1 = d/sd^3
	}
	else {
		a1 = a
		b1 = b
		c1 = c
		d1 = d
		}
	
	# Find the AIC value associated with the cubic model 
	AIC3 = AIC(glm.cubic)
	# Storage for the cubic model
	R[start_ind,6] = b1            # Slope coefficient (associated with data_in)
	R[start_ind,8] = c1           # Slope coefficient (associated with data_in^2)
	R[start_ind,10] = d1           # Slope coefficient (associated with data_in^3)
	
	# Isolate the standard errors for the cubic model 
 	SEs = sqrt(diag( vcov(glm.cubic)))
	R[start_ind,7] = SEs[[2]]	                              # SE of slope (associated with data_in)  
	R[start_ind,9] = SEs[[3]]                               # SE of slope (associated with data_in^2)  
	R[start_ind,11] = SEs[[4]]                              # SE of slope (associated with data_in^3)  

	inds_R2 = which(ldata_out != log(DT))
	mean1 = mean(ldata_out[inds_R2])
	R[start_ind,12] = 1 - sum(residuals(glm.cubic)[inds_R2]^2)/sum((ldata_out[inds_R2]-mean1)^2)   # Store R2 for cubic model
	
	
	R[temp3, 13] <-fitted.values(glm.cubic)        # Store fitted values
	R[temp3, 14] <- residuals(glm.cubic)           # Store residual values 
	
   res = residuals(glm.cubic)
   R[temp3, 14] <- res            # Store residual values  
   R[start_ind,20] = AIC3                                   # Store AIC for cubic model
	
	return(list(glm.cubic = glm.cubic, R=R, AIC3 = AIC3, scaled = scaled))
}
