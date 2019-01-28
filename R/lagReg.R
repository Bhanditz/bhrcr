###################################################################
#    Function to fit a lag regression model as per specifiation   #
###################################################################
lagReg = function(R, data_in, data_out, i, start_ind, end_ind, code, Condition_on_24_hours, fact, Threshold2, Threshold3, plotting, Name, FirstHours, MaxDiffInFirstHours, NeededInFirstHours, imageType=c("pdf", "png")){
  par(ask = FALSE)
  imageType <- match.arg(imageType)
	# Initialise indexes for storage of results
	indexes = start_ind:end_ind
	# Initialise where to store results - only where outliers have not been detected
    temp2 = which(R[start_ind:end_ind,4]==0)
	# Store the locations in current data set (start_ind:end_ind) where outliers were not identified
    temp3 = indexes[temp2]          		
	# Pull out data to fit the model to - only use data that has not been identified as outliers
	data_out_store = data_out
	data_in_store = data_in
	data_in = data_in[temp2]
	data_out = data_out[temp2]	
	# Take the logarithm of the output data
	ldata_out = log(data_out)	
	# Find how many data points we have to fit the model to 
	n2 = length(ldata_out)	
  	# Store the log(para) data.
	R[temp3,5] = ldata_out	
	
	# assume there is no model estimation
	lag_phase_attempted = 2

	 NUM = data_in[which(data_in <= FirstHours)]
 	if (length(NUM)<=1){
 		NE = 1
 		}
 	else if ((max(diff(NUM))> MaxDiffInFirstHours) | (length(NUM)<NeededInFirstHours)){
 		NE = 1
    }
 	else{
 		NE = 0
 		}

   	filename2 <-sprintf("%s_%s.%s", Name, code, imageType)
   	
   	
    prof_type = -9999

 	# initialise k and tlag
	k = NULL; tlag = data_in[1]
	# Do not fit a model becuase of two few data points to fit even a linear model to 
	if (n2 < 3){
		xx<-sprintf("There are not enough time points - no model fitted: data set %s", code)
		R[start_ind,56] = 1
		}
	# Do not fit a model becauuse the baseline level is too low
	else if (data_out[1]< Threshold2){
		xx<-sprintf("Baseline level is too low - no model fitted: data set %s",code)
		R[start_ind,56] = 2

		}
	else if (data_out[length(data_out)]> Threshold3) {
		xx<-sprintf("Parasite is OBVIOUSLY not cleared - no model is fitted: data set %s",code)
		R[start_ind,56] = 3.2

		}
	else if (NE==1){
	###########################################
	#      fit a linear model (only option)   #
	###########################################
	Linear = fitLinearAndStorage(data_in, ldata_out, start_ind, temp3, R)
	glm.linear = Linear$glm.linear; R= Linear$R; AIC1 = Linear$AIC1
	
 	mod = glm.linear
 	LinPart = seq(1,n2)
    # Identify the rate of clearance (as the slope of the linear model)
    k = mod$coefficients[[2]]
    SEs2 = sqrt(diag(vcov(mod)))
	R[start_ind,53] = SEs2[[2]]  
    # Output message about which type of model is being fitted and why
    xx<-sprintf("Fitting a linear model since we don't have enough data during the first 24 hours so we can't measure a lag phase: data set %s",code)
    # Set the profile type 
    prof_type = 1
    # Store the R2 value for the linear model
    R[start_ind,44] = R[start_ind, 23]		
 	} 
    # If we have only 3 data points - only allow option to fit a linear model
	else if (n2==3){		
	###########################################
	#      fit a linear model (only option)   #
	###########################################
	Linear = fitLinearAndStorage(data_in, ldata_out, start_ind, temp3, R)
	glm.linear = Linear$glm.linear; R= Linear$R; AIC1 = Linear$AIC1
	
 	mod = glm.linear
 	LinPart = seq(1,n2)
    # Identify the rate of clearance (as the slope of the linear model)
    k = mod$coefficients[[2]]
    SEs2 = sqrt(diag(vcov(mod)))
	R[start_ind,53] = SEs2[[2]]  
    # Output message about which type of model is being fitted and why
    xx<-sprintf("Fitting a linear model since it has the lowest AIC n = 3: data set %s",code)
    # Set the profile type 
    prof_type = 1
    # Store the R2 value for the linear model
    R[start_ind,44] = R[start_ind, 23]		
 	}
   	# We don't have at least three data points in the first 14 hours so we can't measure a lag phase - fit a linear model 

  # If we have only 4 data points - only allow option to fit a linear or quadratic model
  else if (n2==4){ 
  	###########################################
	#      Fit a linear model               	#
	###########################################
	Linear = fitLinearAndStorage(data_in, ldata_out, start_ind, temp3, R)
	glm.linear = Linear$glm.linear; R= Linear$R; AIC1 = Linear$AIC1

    ###########################################
	#      Fit a quadratic model	          #
	###########################################
	Quad = fitQuadraticAndStorage(data_in, ldata_out, start_ind, temp3, R)
	glm.quad = Quad$glm.quad; R= Quad$R; AIC2 = Quad$AIC2

    # Find the model with the minimum AIC value
	ind = which(c(AIC1, AIC2) == min(c(AIC1, AIC2)))

	##################################################################
	#  Use model with minimum AIC	to find clearance rate etc         #
	##################################################################
 if (data_out[2] > (1+fact)*data_out[1]){
    # Find the indices from the maximum value of the parasites level to the end of the data set
    #inds =  which(max(ldata_out)==ldata_out)
    # Isolate the part of the input data set past the maximum value
    data_in_mod = data_in[2:n2]
    # Isolate the part of the data set past the maximum value
    ldata_out_mod = ldata_out[2:n2]
    tlag = data_in_mod[1]
   	# Fit a linear model to the data past the maximum value
    MaxReg<-glm(ldata_out_mod ~ data_in_mod) 
    # Find the AIC value associated with the linear model    
	AICMaxReg = AIC(MaxReg)
	#sum of squared residuals for all datapoints starting from the second.
	ss_1 = sum(residuals(glm.linear)[2:n2]^2)
	ss_2 = sum(residuals(glm.quad)[2:n2]^2)
	ss_3 = sum(residuals(MaxReg)[1:3]^2)
	ind = which(c(ss_1, ss_2, ss_3) == min(c(ss_1, ss_2, ss_3)))
  }	
 #Fit a linear model, since it has the lowest AIC  (ind ==1)
 if (ind == 1){ 
 	 mod = glm.linear
 	 LinPart = seq(1,n2)
    # Identify the rate of clearance (as the slope of the linear model)
    k = mod$coefficients[[2]]
    tlag = data_in[1]
    SEs2 = sqrt(diag(vcov(mod)))
	 R[start_ind,53] = SEs2[[2]]  
    # Output message about which type of model is being fitted and why
    xx<-sprintf("Fitting a linear model since it has the lowest SS n = 4: data set %s",code)
    # Set the profile type 
    prof_type = 1
    # Store the R2 value for the linear model
    R[start_ind,44] = R[start_ind, 23]		    
	} 
  	# Fit a quadratic model, since it has the lowest AIC
 else if (ind == 2){
    	# Set the model fitted to be the quadratic one 
    	mod = glm.quad
		
		# Find the a,b,c values in the quadratic fit (ln(para) = a*time^2 + b*time + c)
		a = mod$coefficients[[3]]
		b = mod$coefficients[[2]]
		c = mod$coefficients[[1]]
		
    	# Calculate the fitted values for the model
		mod.fitted <-fitted.values(mod)
		# Calculate the slopes between consecutive points (using the model fitted values)
		slopes = diff(mod.fitted)/diff(data_in)
		
		# Isolate the candidate points for the slowest and fastest rates: the left most point, the right most points
		c_points = c(data_in[1], max(data_in)) 
		# Isolate the slopes of the cubic at the candidate points (slope = first derivative of cubic function)
		c_devs = 2*a*c_points + b
	    
     # Find the fastest absolute rate out of the three candidate points
		ind_max = 2
		# Isolate the fastest rate
		fastest_rate = c_devs[ind_max]
		ratios = fastest_rate/slopes
		LinPart = c(which((ratios <= 5)&(ratios>=0)),n2)
		
		# Isolate the linear parts of the input data
		ldattemp = mod.fitted[LinPart] 
		
		# Convex: if the 'a' value is negative then we have a lag phase to identify -> Slow initial clearance followed by faster clearance later
    	if (a <0){
			# Output message about which type of model is being fitted and why
			xx<-sprintf("Fitting a quadratic (convex) model since it has the lowest SS and n = 4: data set %s",code)
			# Set the profile type 
			prof_type = 3			
		}
		# Concave: if the 'a' value is postive then we have a fast initial clearance followed by slower clearance later. There is no lag phase identified and 		we fit a linear model
		else if (a>=0) { 
        	 # Output message about which type of model is being fitted and why
   	     	xx<-sprintf("Fitting a linear model since we have identified a concave quadratic (quadratic has the lowest SS) and n = 4: data set %s",code)
         	# Set the profile type 
         	prof_type = 2
			LinPart = seq(1,n2)
			# Isolate the linear parts of the output data
			ldattemp = ldata_out[LinPart] 
		
    	}	   
    	
		tlag = data_in[LinPart[1]]	
		
		if (tlag == data_in[1]){
    		LinPart = seq(1,n2)
    		ldattemp = ldata_out[LinPart] 
    		# Store the clearance rate -> equal to the slope of the identified 'linear part'
    		to_include_R2 = LinPart
    	}else{
    		to_include_R2 = seq(1,length(LinPart))
    	}
    	
    		
		# Isolate the linear parts of the input data
		dattemp = data_in[LinPart]
		# Fit a linear model the identified 'linear part'
		mod<-glm(ldattemp ~ dattemp)
        # Store the locations in current data set (start_ind:end_ind) where outliers were not identified
        temp3 = indexes[temp2]     
         
		# Store the R2 value for the 'linear part' model
      #R[start_ind,44] = 1 - sum(residuals(mod)^2)/sum((ldattemp-mean(ldattemp))^2) # modified 29/5/2013, for the quadratic model
   		R[start_ind,44] = 1 - sum((fitted(mod)[to_include_R2]-ldata_out[LinPart])^2)/sum((ldata_out[LinPart]-mean(ldata_out[LinPart]))^2)

        # Store the clearance rate -> equal to the slope of the identified 'linear part'
		k = mod$coeff[[2]]			
	    SEs2 = sqrt(diag(vcov(mod)))
		R[start_ind,53] = SEs2[[2]]  

		# Store the difference between the fitted value from the 'linear part' at the time lag value and the first value in the 'linear part' model
		R[start_ind, 50] <- mod.fitted[which(data_in == tlag)[1]] - mod.fitted[1]
		R[start_ind, 50] <- mod.fitted[which(data_in == tlag)[1]] - mod.fitted[1]
		# Store the quadratic model fitted values (the original model fitted)    
        glm.quad.fitted <-fitted.values(glm.quad)
		# Find the difference between the maximum quadratic fitted value and the fitted value (from the quadratic model) at the time lag
        R[start_ind, 51] <- max(glm.quad.fitted) -  glm.quad.fitted[which(data_in == tlag)[1]]
	 
  	}
	 
 else if (ind == 3){
 	mod<-MaxReg 
    # Find the AIC value associated with the linear model    
	AIC1 = AIC(mod)
	LinPart = 2:n2
	
    # Identify the rate of clearance (as the slope of the linear model)
    k = mod$coefficients[[2]]
    SEs2 = sqrt(diag(vcov(mod)))
	R[start_ind,53] = SEs2[[2]]  
    # Output message about which type of model is being fitted and why
    xx<-sprintf("Fitting a linear model from max value since we only have 4 data points and an initial increase in parasite level (also - lowest SS): data set %s",code)
    # Set the profile type 
    prof_type = 5
         
    # Store the R2 value for the linear model
    R[start_ind,44] = 1 - sum(residuals(mod)^2)/sum((ldata_out_mod-mean(ldata_out_mod))^2)
    
 	}
 	}
  # If we have more than 4 data points - allow option to fit linear, quadratic or cubic model
  else {
	#############################################
	#      Fit a linear model               	#
	#############################################
	Linear = fitLinearAndStorage(data_in, ldata_out, start_ind, temp3, R)
	glm.linear = Linear$glm.linear; R= Linear$R; AIC1 = Linear$AIC1

    ###########################################
	#      Fit a quadratic model	          #
	###########################################
	Quad = fitQuadraticAndStorage(data_in, ldata_out, start_ind, temp3, R)
	glm.quad = Quad$glm.quad; R= Quad$R; AIC2 = Quad$AIC2

	#################################
	#      Fit a cubic model		#
	#################################
	Cubic = fitCubicAndStorage(data_in, ldata_out, start_ind, temp3, R)
	glm.cubic = Cubic$glm.cubic; R=Cubic$R; AIC3 = Cubic$AIC3
	
    # Find the model with the minimum AIC value
	ind = which(c(AIC1, AIC2, AIC3) == min(c(AIC1, AIC2, AIC3)))
	
    ##################################################################
	#  Use model with minimum AIC	to find clearance rate etc       #
	##################################################################
	# Fit a linear model, since it has the lowest AIC  (ind ==1)
  if (ind == 1){ 
	mod = glm.linear
	
	LinPart = seq(1,n2)
    # Identify the rate of clearance (as the slope of the linear model)
    k = mod$coefficients[[2]]
    SEs2 = sqrt(diag(vcov(mod)))
	R[start_ind,53] = SEs2[[2]]  
    # Output message about which type of model is being fitted and why
    xx<-sprintf("Fitting a linear model since it has the lowest AIC: data set %s",code)
    # Set the profile type 
    prof_type = 1
    # Store the R2 value for the linear model
    R[start_ind,44] = R[start_ind, 23]		
       
	}
  	# Fit a quadratic model, since it has the lowest AIC
  	else if (ind == 2){
		# Set the model fitted to be the quadratic one 
    	mod = glm.quad
		
		# Find the a,b,c values in the quadratic fit (ln(para) = a*time^2 + b*time + c)
		a = mod$coefficients[[3]]
		b = mod$coefficients[[2]]
		c = mod$coefficients[[1]]
		
    	# Calculate the fitted values for the model
		mod.fitted <-fitted.values(mod)
		# Calculate the slopes between consecutive points (using the model fitted values)
		slopes = diff(mod.fitted)/diff(data_in)
		
		# Isolate the candidate points for the slowest and fastest rates: the left most point, the right most points
		c_points = c(data_in[1], max(data_in)) 
		# Isolate the slopes of the cubic at the candidate points (slope = first derivative of cubic function)
		c_devs = 2*a*c_points + b
	
	    # Find the fastest absolute rate out of the candidate points
		ind_max = 2
		# Isolate the fastest rate
		fastest_rate = c_devs[ind_max]
		ratios = fastest_rate/slopes
		LinPart = c(which((ratios <= 5)&(ratios>=0)),n2)
		
		# Isolate the linear parts of the input data
		ldattemp = mod.fitted[LinPart] 
		
		# Convex: if the 'a' value is negative then we have a lag phase to identify -> Slow initial clearance followed by faster clearance later
    	if (a <0){
			# Output message about which type of model is being fitted and why
			xx<-sprintf("Fitting a quadratic (convex) model since it has the lowest AIC: data set %s",code)
			# Set the profile type 
			prof_type = 3			
		}
		# Concave: if the 'a' value is postive then we have a fast initial clearance followed by slower clearance later. There is no lag phase identified and 		we fit a linear model
		else if (a>=0) { 
        	 # Output message about which type of model is being fitted and why
   	     	xx<-sprintf("Fitting a linear model since we have identified a concave quadratic (quadratic has the lowest AIC): data set %s",code)
         	# Set the profile type 
         	prof_type = 2
			LinPart = seq(1,n2)
			# Isolate the linear parts of the output data
			ldattemp = ldata_out[LinPart] 
		
    	}	   
    	
		tlag = data_in[LinPart[1]]	
		
		if (tlag ==  data_in[1]){
    		LinPart = seq(1,n2)
    		ldattemp = ldata_out[LinPart] 
    		# Store the clearance rate -> equal to the slope of the identified 'linear part'
    		to_include_R2 = LinPart
    	}else{
    		to_include_R2 = seq(1,length(LinPart))
    	}
    		
		# Isolate the linear parts of the input data
		dattemp = data_in[LinPart]
		# Fit a linear model the identified 'linear part'
		mod<-glm(ldattemp ~ dattemp)
        # Store the locations in current data set (start_ind:end_ind) where outliers were not identified
        temp3 = indexes[temp2]     
          
		# Store the R2 value for the 'linear part' model
      #R[start_ind,44] = 1 - sum(residuals(mod)^2)/sum((ldattemp-mean(ldattemp))^2) # modified 29/5/2013, for the quadratic model
		R[start_ind,44] = 1 - sum((fitted(mod)[to_include_R2]-ldata_out[LinPart])^2)/sum((ldata_out[LinPart]-mean(ldata_out[LinPart]))^2)

        # Store the clearance rate -> equal to the slope of the identified 'linear part'
		k = mod$coeff[[2]]
		SEs2 = sqrt(diag(vcov(mod)))
		R[start_ind,53] = SEs2[[2]]  			

		# Store the difference between the fitted value from the 'linear part' at the time lag value and the first value in the 'linear part' model
		R[start_ind, 50] <- mod.fitted[which(data_in == tlag)[1]] - mod.fitted[1]
		R[start_ind, 50] <- mod.fitted[which(data_in == tlag)[1]] - mod.fitted[1]
		# Store the quadratic model fitted values (the original model fitted)    
        glm.quad.fitted <-fitted.values(glm.quad)
		# Find the difference between the maximum quadratic fitted value and the fitted value (from the quadratic model) at the time lag
        R[start_ind, 51] <- max(glm.quad.fitted) -  glm.quad.fitted[which(data_in == tlag)[1]]
	 }
 
   # Fit a cubic model, since it has the lowest AIC
   else if (ind == 3){
	   	# Set the model fitted to be the cubic one 
    	mod = glm.cubic
		
		# Find the 'a', 'b' and 'c' values in the cibic fit (ln(para) = a*time^3 + b*time^2 + c*time + d)
		a = mod$coefficients[[4]] # a value
		b = mod$coefficients[[3]] # b value
		c = mod$coefficients[[2]] # c value
		d = mod$coefficients[[1]] # b value
		
    	# Calculate the point of inflection of the cubic model (-b/3a) which is found by putting the second derivative to 0
		PoI = -b/(3*a)
		# Isolate the candidate points for the slowest and fastest rates: the left most point, the PoI, the right most points
		c_points = c(data_in[1], PoI,  max(data_in)) 
		# Isolate the slopes of the cubic at the candidate points (slope = first derivative of cubic function)
		c_devs = 3*a*c_points^2 + 2*b*c_points + c
		
		# Find the slowest absolute rate out of the three candidate points
		ind_min = which(abs(c_devs) == min(abs(c_devs)))
		# Isolate the slowest rate
		slowest_rate = c_devs[ind_min]
		# Find the fastest absolute rate out of the three candidate points
		ind_max = which(abs(c_devs) == max(abs(c_devs)))
		# Isolate the fastest rate
		fastest_rate = c_devs[which(c_devs == min(c_devs))]
		ind_max_neg = which(c_devs == min(c_devs))
		fastest_rate = c_devs[ind_max_neg]			
		# Calculate the fitted values for the model
		mod.fitted <-fitted.values(mod)
		# Calculate the slopes between consecutive points (using the model fitted values)
		slopes = diff(mod.fitted)/diff(data_in)
		ratios = fastest_rate/slopes
		FastParts = which((ratios <= 5)&(ratios>=0))
			
		if (((PoI < data_in[1]) | (PoI > max(data_in))) & (b<0)){ #  PoI is outside range and convex
 			# Output message about which type of model is being fitted and why
			xx<-sprintf("Fitting a cubic model for cubic case (d) where b <0: data set %s",code)
			
			fastest_rate = min(c(c_devs[1], c_devs[3]))
			ratios = fastest_rate/slopes
			# Set the profile type 
			prof_type = 7
			LinPart = which((ratios <= 5)&(ratios>=0))
			temp = seq(1,n2)
			LinPart = c(LinPart, temp[LinPart[length(LinPart)]+1])
			ldattemp = mod.fitted[LinPart] 
			
			
			}
		else if (((PoI < data_in[1]) | (PoI > max(data_in))) & (b>0)){ # PoI is outside range and concave
			# Output message about which type of model is being fitted and why
    	    xx<-sprintf("Fitting a linear model for cubic case (d) where b >0: data set %s",code)
        	# Set the profile type 
        	prof_type = 6
			LinPart = seq(1,n2)
			# Isolate the linear parts of the input data
			ldattemp = ldata_out[LinPart] 
			}
		else if (ind_max_neg==2) {# max negative rate of change occurs at the PoI
			# Output message about which type of model is being fitted and why
        	xx<-sprintf("Fitting a cubic model for cubic case (a): data set %s",code)
			# Set the profile type 
       		prof_type = 4
       		
          	LinPart = which((ratios <= 5)&(ratios>=0)&(slopes<0))
			temp = seq(1,n2)
			LinPart = c(LinPart, temp[LinPart[length(LinPart)]+1])
			ldattemp = mod.fitted[LinPart] 
			}
		else if ((ind_max_neg!=2) & (FastParts[1]==1)){  # max negative rate of change doesn't occur at the PoI and there is a fast change at the beginning of the profile
			## Output message about which type of model is being fitted and why
            xx<-sprintf("Fitting a linear model for cubic case (biii) where there is no lag phase: data set %s",code)
	       	# Set the profile type 
           	prof_type = 9
           	LinPart = seq(1,n2)
           	ldattemp = ldata_out[LinPart] 
			}
		else if  ((ind_max_neg!=2) & (FastParts[1]!=1)){ # max negative rate of change doesn't occur at the PoI and there is NOT a fast change at the beginning of the profile
			xx<-sprintf("Fitting a linear model for cubic case (bii) where there is a lag phase: data set %s",code)
			# Set the profile type 
          	prof_type = 8					
          	LinPart = which((ratios <= 5)&(ratios>=0)&(data_in[2:n2]>PoI))
			temp = seq(1,n2)
			LinPart = c(LinPart, temp[LinPart[length(LinPart)]+1])
			ldattemp = mod.fitted[LinPart]           	
			}
   		
		tlag = data_in[LinPart[1]]		
		
		if (tlag ==  data_in[1]){
    		LinPart = seq(1,n2)
    		ldattemp = ldata_out[LinPart] 
    		# Store the clearance rate -> equal to the slope of the identified 'linear part'
			to_include_R2 = LinPart
    	}else{
    		to_include_R2 = seq(1,length(LinPart))
    		}


		# Isolate the linear parts of the output data
		dattemp = data_in[LinPart]
		# Fit a linear model the identified 'linear part'
		mod<-glm(ldattemp ~ dattemp)
        # Store the locations in current data set (start_ind:end_ind) where outliers were not identified
        temp3 = indexes[temp2]     
          
		# Store the R2 value for the 'linear part' model
      # R[start_ind,44] = 1 - sum(residuals(mod)^2)/sum((ldattemp-mean(ldattemp))^2) # modified 29/5/2013, for the cubic model
		R[start_ind,44] = 1 - sum((fitted(mod)[to_include_R2]-ldata_out[LinPart])^2)/sum((ldata_out[LinPart]-mean(ldata_out[LinPart]))^2)

        # Store the clearance rate -> equal to the slope of the identified 'linear part'
		k = mod$coeff[[2]]		
		
	    SEs2 = sqrt(diag(vcov(mod)))
		R[start_ind,53] = SEs2[[2]]   
		# Store the difference between the fitted value from the 'linear part' at the time lag value and the first value in the 'linear part' model
		R[start_ind, 50] <- mod.fitted[which(data_in == tlag)[1]] - mod.fitted[1]
		R[start_ind, 50] <- mod.fitted[which(data_in == tlag)[1]] - mod.fitted[1]
		# Store the quadratic model fitted values (the original model fitted)    
        glm.cubic.fitted <-fitted.values(glm.cubic)
		# Find the difference between the maximum quadratic fitted value and the fitted value (from the quadratic model) at the time lag
        R[start_ind, 51] <- max(glm.cubic.fitted) -  glm.cubic.fitted[which(data_in == tlag)[1]]
   	}
  }  
  
  
  
    if ((plotting==TRUE)){     
	if (imageType == "png") {
    # MPF 2015/02/25 Cairo attached by Depends
    #library(Cairo)
    Cairo(file=filename2, type="png")
  } else if (imageType == "pdf") {
    pdf(file=filename2)
  }
   			filename <-sprintf("patient id: %s", code)
		 	par( mar = c( 5.1, 5.1, 4.1, 2.1 ) )
		 	x0 = data_in[1]
		 	xend = data_in[length(data_in)]*1.8
		 	y0 = 0
		 	yend = max(ldata_out)*1.05
		 	if (yend==-Inf){yend = 10}
		 	xlim=c(x0,xend)
		 	ylim=c(y0,yend)
		 	plot(data_in, ldata_out, main = filename, xlab="Time (hours)", ylab = expression(log[e](parasitaemia)), cex = 2,cex.main=2, cex.lab=2, cex.axis=2, cex.main=2,xlim=xlim, ylim=ylim, type = 'p', pch=19)
		 	
		 	if (R[start_ind,56]==0){
			if (prof_type != 5){
				lines(data_in[LinPart], fitted(mod), lty=1, col = "blue",lwd=3)
			}
			else {
				lines(data_in[LinPart+1], fitted(mod), lty=1, col = "blue",lwd=3)
				}
			}
				
			out_temp = floor(as.numeric(R[start_ind:end_ind,4]))
			outliers = which(out_temp ==2.0)
    	if (length(outliers)>0){
  			
  			points(data_in_store[outliers], log(data_out_store[outliers]), col = "yellow", pch = 19, cex = 2) # put outliers in pink
    		}
    	
    		out_temp = floor(as.numeric(R[start_ind:end_ind,4]))
			outliers = which(out_temp ==3.0)
    	
    	if (length(outliers)>0){
  			
  			points(data_in_store[outliers], log(data_out_store[outliers]), col = "red", pch = 19, cex = 2)
    		}
    		
    		
    	
    	if ((tlag >0)&(R[start_ind,56]==0)){
    		t_temp = which(data_in <tlag)
    		points(data_in[t_temp], log(data_out[t_temp]), col = "gray", pch = 19, cex = 2)
    		}
    		
		dev.off()
	}
	
	
	print(xx)
	if (is.null(k) == FALSE){
		   	# store R2 value for 'linear model' where measuremetns below the Detect level are ignored
			mod_lin <- lm(ldata_out~data_in)
			R[start_ind,70]= 1-sum(residuals(mod_lin)^2)/sum((ldata_out-mean(ldata_out))^2)
		if (tlag == data_in[1]){
			MSRE = (mean(residuals(mod)^2))^(1/2)
			}
		else {
			MSRE = sqrt(mean((fitted(mod) - ldata_out[LinPart])^2))
			}
		R[start_ind, 55] = MSRE
		R[temp3[LinPart], 52] = 1                             # Store an indicator value of 1 for those points used in the 'linear part'
		R[temp3[LinPart], 45] <-fitted.values(mod)            # Store fitted values
		R[temp3[LinPart], 46] <- residuals(mod)               # Store residual values
		R[temp3[LinPart], 47] <-ls.diag(mod)$stud.res         # Store student residual values
		R[temp3[LinPart], 48] <-ls.diag(mod)$std.res          # Store standardised residual values


		# Store the time lag for the model
    	R[start_ind, 49] <- tlag -  data_in[1]
    	# Store the profile type for the model
    	R[start_ind, 42] <- prof_type
    	# Store clearance rate for the model
   	    R[start_ind, 43] <- k   	    
   	    alin = mod$coeff[[1]]	 	
      	knum =as.numeric(k)
      	int = alin + knum*as.numeric(R[start_ind, 49]) 
       R[start_ind,67] = int # intercept at tlag

       
   	    PC50 = (log(data_out[1]*(1-0.5)) - alin)/knum # PC50
       PC90 = (log(data_out[1]*(1-0.9)) - alin)/knum  # PC90
       PC95 = (log(data_out[1]*(1-0.95)) - alin)/knum # PC95
       PC99 = (log(data_out[1]*(1-0.99)) - alin)/knum  # PC99
   	    
   	    if (PC50>0){R[start_ind,59] =PC50} # PC50
   	    if (PC90>0){R[start_ind,60] =PC90} # PC50
   	    if (PC95>0){R[start_ind,61] =PC95} # PC50
   	    if (PC99>0){R[start_ind,62] =PC99} # PC50
       

	   	lag_phase_attempted = 1 # assume no lag estimation
  	    if ((NE==0)|(n2==3)){
   	    		lag_phase_attempted = 0
   	    	}
	    	
  	   
   	    # store the change_max and change_tlag
		R[start_ind, 50] <- data_out[which(data_in == tlag)[1]] -  max(data_out[which(data_in <=tlag)])
	    R[start_ind, 51] <- data_out[which(data_in == tlag)[1]] -  data_out[1]

	}
 R[start_ind, 57] <- lag_phase_attempted
  
	# Return the array R
 	return(list(R=R))
}
