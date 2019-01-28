fitLinearAndStorage_tobit<-function(data_in, ldata_out, start_ind, temp3, R, DT, MaxIters){ 
	###########################################
	#      Fit a linear model               	#
	###########################################
	
	# transform data for fitting purposes:
	mu = mean(data_in) 
	sd = sd(data_in)
	data_in_new = (data_in - mu)/sd
	scaled = 1
	
	ind1 = which(ldata_out == log(DT))
	if (length(ind1) != 0){
		m2 <- lm(ldata_out[-ind1] ~ poly(data_in_new[-ind1], 1, raw = TRUE)) 
		}
	else{
		m2 <- lm(ldata_out[-length(data_in)] ~ poly(data_in_new[-length(data_in)], 1, raw = TRUE)) 
		}
	glm.linear <- try(tobit(ldata_out ~ data_in_new, left = log(DT), iter.max = MaxIters, init = coef(m2)))
	if (class(glm.linear)[1]=="try-error"){
		yy<-sprintf("Trying to fix error in tobit fitting (1) - linear")
		print(yy)
		m2 <- lm(ldata_out ~ poly(data_in_new, 1, raw = TRUE))
		glm.linear <-try(tobit(ldata_out ~ data_in_new, left = log(DT), iter.max = MaxIters, init = coef(m2)) ) 
		
		}
	if (class(glm.linear)[1]=="try-error"){
		yy<-sprintf("Trying to fix error in tobit fitting (2) - linear")
		print(yy)
		glm.linear <- try(tobit(ldata_out ~ data_in_new, left = log(DT), iter.max = MaxIters))
		}
	if (class(glm.linear)[1]=="try-error"){
		yy<-sprintf("Trying to fix error in tobit fitting (3) - linear")
		print(yy)
		scaled = 0
		m2 <- lm(ldata_out ~ poly(data_in, 1, raw = TRUE)) 
		glm.linear <- try(tobit(ldata_out ~ poly(data_in, 1, raw = TRUE), left = log(DT), iter.max = MaxIters, init = coef(m2)))
		
		}
	if (class(glm.linear)[1]=="try-error"){
		yy<-sprintf("Trying to fix error in tobit fitting (4) - linear")
		print(yy)
		scaled = 0
		ldata_out[length(ldata_out)] = 0
		m2 <- lm(ldata_out ~ poly(data_in, 1, raw = TRUE)) 
		glm.linear <- try(tobit(ldata_out ~ poly(data_in, 1, raw = TRUE), left = log(DT), iter.max = MaxIters, init = coef(m2)))
		}
	if (class(glm.linear)[1]=="try-error"){
		yy<-sprintf("Trying to fix error in tobit fitting (5) - linear")
		print(yy)
		scaled = 0
		m2 <- lm(ldata_out[-ind1] ~ poly(data_in[-ind1], 1, raw = TRUE)) 
		glm.linear <- try(tobit(ldata_out ~ data_in, left = log(DT), init = coef(m2)))
		}
	if (class(glm.linear)[1]=="try-error"){
		yy<-sprintf("Trying to fix error in tobit fitting (6) - linear")
		print(yy)
		scaled = 0
		m2 <- lm(ldata_out[-1] ~ poly(data_in[-1], 1, raw = TRUE)) 
		glm.linear<- try(tobit(ldata_out ~ data_in, left = log(DT), init = coef(m2)))
		}
	if (class(glm.linear)[1]=="try-error"){
	  yy<-sprintf("Trying to fix error in tobit fitting (7) - linear")
	  print(yy)
	  scaled = 0
	  m2 <- lm(ldata_out[-length(data_in)] ~ poly(data_in_new[-length(data_in)], 1, raw = TRUE))
	  glm.linear<- try(tobit(ldata_out ~ data_in, left = log(DT), init = coef(m2)))
	}	
	if (class(glm.linear)[1]=="try-error"){
	  yy<-sprintf("Trying to fix error in tobit fitting (8) - linear")
	  print(yy)
	  scaled = 0
	  glm.linear<- try(tobit(ldata_out ~ data_in, left = log(DT),iter.max = MaxIters,  init = c(ldata_out[1], 3)))
	}	
	if (class(glm.linear)[1]=="try-error"){
	  yy<-sprintf("Trying to fix error in tobit fitting (9) - linear")
	  print(yy)
	  scaled = 0
	  glm.linear<- try(tobit(ldata_out ~ data_in, left = log(DT),iter.max = MaxIters,  init = c(ldata_out[1], coef(m2)[2])))
	}		
		
	if (class(glm.linear)[1]=="try-error"){
		yy<-sprintf("NOOOOOO --> Did not fix error in tobit fitting!! - linear")
		print(yy)
		glm.linear <- tobit(ldata_out ~ data_in_new, left = log(DT), iter.max = MaxIters)
		}

	# Extract coefficents	
	a = coef(glm.linear)[[1]]
	b = coef(glm.linear)[[2]]
	
	# Transform coefficients back to our orignial model (if necessary)	
	if (scaled == 1){
		a1 = a - b*mu/sd
		b1 = b/sd
		slope_se_factor = sd
	}
    else {
		a1 = a
		b1 = b
		slope_se_factor = 1
		}

	 	
	# Find the AIC value associated with the linear model 
	AIC1 = AIC(glm.linear)
	mean1 = mean(ldata_out)
	# Storage for the linear model
	R[start_ind,21] = b1            # Slope coefficient (associated with data_in)
	R[start_ind,71] = a1     # intercept coefficient - for output to user later on


    # Isolate the standard errors for the linear model           
	SEs = sqrt(diag( vcov(glm.linear)))                         
	R[start_ind,22] = 1/slope_se_factor*SEs[[2]]                         # SE of slope    
	
	# Find the mean of the data for calculation of the R2 value of linear model              
	inds_R2 = which(ldata_out != log(DT))
	mean1 = mean(ldata_out[inds_R2])
	R[start_ind,23] = 1 - sum(residuals(glm.linear)[inds_R2]^2)/sum((ldata_out[inds_R2]-mean1)^2)    # R2 value for linear model
	
	R[temp3, 24] <-fitted.values(glm.linear)         # Store fitted values
	res = residuals(glm.linear)
   R[temp3, 25] <- res            # Store residual values
  
 
	R[start_ind, 28] <-summary(glm.linear)$coefficients[[2,4]]    # Store p-value associated with slope in linear model
	R[start_ind,29] = AIC1                                    # Store AIC for linear model


	return(list(glm.linear = glm.linear, R=R, AIC1 = AIC1, a1 = a1))
}
