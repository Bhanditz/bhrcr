fitQuadraticAndStorage_tobit<-function(data_in, ldata_out, start_ind, temp3, R, DT, MaxIters){ 
    ###########################################
	#      Fit a quadratic model	          #
	###########################################
	# Initialise the square of the time data
	data_in2 = data_in^2
	# Fit quadratic model
	scaled = 1
	
	# transform data for fitting purposes:
	mu = mean(data_in)
	sd = sd(data_in) 
	data_in_new = (data_in - mu)/sd
	data_in_new2 = data_in_new^2
	
	ind1 = which(ldata_out == log(DT))
	
	if (length(ind1) != 0){
		m2 <- lm(ldata_out[-ind1] ~ poly(data_in_new[-ind1], 2, raw = TRUE)) 
		}
	else{	
		m2 <- lm(ldata_out[-length(data_in)] ~ poly(data_in_new[-length(data_in)], 2, raw = TRUE)) 
	}
		
	glm.quad <- try(tobit(ldata_out ~ data_in_new + data_in_new2, left = log(DT), iter.max = MaxIters))
	if (class(glm.quad)[1]=="try-error"){
		yy<-sprintf("Trying to fix error in tobit fitting (1) - quad")
		print(yy)
		glm.quad <-try(tobit(ldata_out ~ data_in_new + data_in_new2, left = log(DT), iter.max = MaxIters, init = coef(m2))) 
		
		}
	if (class(glm.quad)[1]=="try-error"){
		yy<-sprintf("Trying to fix error in tobit fitting (2) - quad")
		print(yy)
		m2 <- lm(ldata_out ~ poly(data_in_new, 2, raw = TRUE))
		glm.quad <- try(tobit(ldata_out ~ data_in_new + data_in_new2, left = log(DT), iter.max = MaxIters, init = coef(m2)))
		}
	if (class(glm.quad)[1]=="try-error"){
		yy<-sprintf("Trying to fix error in tobit fitting (3) - quad")
		print(yy)
		ldata_out[length(ldata_out)] = 0
		m2 <- lm(ldata_out ~ poly(data_in, 2, raw = TRUE))
		scaled = 0
		glm.quad <- try(tobit(ldata_out ~ poly(data_in, 2, raw = TRUE), left = log(DT), iter.max = MaxIters, init = coef(m2)))
		}	
	if (class(glm.quad)[1]=="try-error"){
		yy<-sprintf("Trying to fix error in tobit fitting (4) - quad")
		print(yy)
		scaled = 0
		glm.quad <- try(tobit(ldata_out ~ poly(data_in, 2, raw = TRUE), left = log(DT), iter.max = MaxIters*1000))
		}
	if (class(glm.quad)[1]=="try-error"){
		yy<-sprintf("Trying to fix error in tobit fitting (5) - quad")
		print(yy)
		scaled = 0
		m2 <- lm(ldata_out[-ind1] ~ poly(data_in[-ind1], 2, raw = TRUE)) 
		glm.quad <- try(tobit(ldata_out ~ data_in + data_in2, left = log(DT), init = coef(m2)))
		}
	if (class(glm.quad)[1]=="try-error"){
		yy<-sprintf("Trying to fix error in tobit fitting (6) - quad")
		print(yy)
		scaled = 0
		m2 <- lm(ldata_out[-1] ~ poly(data_in[-1], 2, raw = TRUE)) 
		glm.quad <- try(tobit(ldata_out ~ data_in + data_in2, left = log(DT), init = coef(m2)))
		}
	if (class(glm.quad)[1]=="try-error"){
		yy<-sprintf("Trying to fix error in tobit fitting (7) - quad")
		print(yy)
		scaled = 0
		m2 <- lm(ldata_out[-1] ~ poly(data_in[-1], 2, raw = TRUE)) 
		glm.quad <- try(tobit(ldata_out ~ data_in + data_in2, left = log(DT), iter.max = MaxIters, init = c(ldata_out[1],coef(m2)[2:3])))
	}
	if (class(glm.quad)[1]=="try-error"){
		yy<-sprintf("NOOOOOO --> Did not fix error in tobit fitting!! - quad")
		print(yy)
		glm.quad <- tobit(ldata_out ~ data_in_new + data_in_new2, left = log(DT), iter.max = MaxIters)
		}
	
	# Extract coefficents	
	a = coef(glm.quad)[[1]]
	b = coef(glm.quad)[[2]]
	c = coef(glm.quad)[[3]]

	# Transform coefficients back to our orignial model (if necessary)
	if (scaled == 1){
	a1 = a - b*mu/sd + c*mu^2/sd^2
	b1 = b/sd - 2*mu*c/sd^2 
	c1 = c/sd^2 
	}
	else {
		a1 = a
		b1 = b
		c1 = c
		}
	
	# Find the AIC value associated with the quad model 
	AIC2 = AIC(glm.quad)

	# Storage for the quad model
	R[start_ind,30] = b1            # Slope coefficient (associated with data_in)
	R[start_ind,32] = c1           # Slope coefficient (associated with data_in^2)

   # Isolate the standard errors for the quadratic model 
    SEs = sqrt(diag( vcov(glm.quad)))
	R[start_ind,31] = SEs[[2]]                               # SE of slope (associated with data_in)   
	R[start_ind,33] = SEs[[3]]                               # SE of slope (associated with data_in^2)
	inds_R2 = which(ldata_out != log(DT))
	mean1 = mean(ldata_out[inds_R2])
	R[start_ind,34] = 1 - sum(residuals(glm.quad)[inds_R2]^2)/sum((ldata_out[inds_R2]-mean1)^2) # Store R2 for quadratic model
	
	

	
	R[temp3, 35] <-fitted.values(glm.quad)          # Store fitted values
	R[temp3, 36] <- residuals(glm.quad)             # Store residual values

   res = residuals(glm.quad)
   R[temp3, 25] <- res            # Store residual values
	R[start_ind,41] = AIC2                                   # Store AIC for quadratic model
	
	return(list(glm.quad = glm.quad, R=R, AIC2 = AIC2, scaled = scaled))
}
