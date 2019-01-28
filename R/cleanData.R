###################################################################
#    Function to clean data for lag regression fitting procedure  #
##################################################################
cleanData<-function(data_in, data_out, i, start_ind, end_ind, R, code, Detec_limit, DT, Threshold2, DaysRecurrence, TimeDiffRecurrence,testing_tobit){
  par(ask = FALSE)
  # initialise indexes for storage of results
	indexes = start_ind:end_ind
	# initialise where to store results - only where outliers have not been detected
	temp2 = which(R[start_ind:end_ind,4]==0)
	data_in_store = data_in
	data_out_store = data_out
	# assume that there are no comments to be made about this data set
	xx <-sprintf("No comment")
	# length of the currnet data set
	n2 = length(data_in)
	# initialise that there are no 0's replaced with 10's
  	ind_DL =c()
	#################################
	#    Profile evaluation         #
	#################################
	# Part (a): Baseline parasitaemia or not enough points at all (ie 1 or 2 points only)
	if (n2 < 3){
		xx <-sprintf("In data set %g, there are not enought time points - no model fitted", i)
		# do not fit a model becuase of two few data points to fit even a linear model to 
		}
	else if (data_out[1]< Threshold2){
		xx <-sprintf("In data set %g, baseline level is too low - no model fitted",i)
		# do not fit a model becauuse the baseline level is too low
		}
	else {
		
	
	# remove data after DaysRecrud days	
	DD = which(data_in>24*DaysRecurrence)
	n2 = length(data_in)
	if (length(DD)>0){
		To_remove6 = seq(DD[1], n2)
		data_in = data_in[-To_remove6]
		data_out = data_out[-To_remove6]
		l = length(To_remove6)
		temp2 = which(R[start_ind:end_ind,4]==0)
		temp3 = indexes[temp2]
		R[temp3[To_remove6],4] = 4.1
		# if there are indices to remove, make a comment as such
		xx <-sprintf("In data set %g due to recrudesence, removing data points - %g ", i,To_remove6[1], To_remove6[l])	
		}	
		
		# remove data after DaysRecrud days	
	diff_datain = diff(data_in)
	ind5 = which(diff_datain> TimeDiffRecurrence)
	n2 = length(data_in)
	if (length(ind5)>1){
		if (ind5[1]>3){
		To_remove7 = seq(ind5[1]+1, n2)
		data_in = data_in[-To_remove7]
		data_out = data_out[-To_remove7]
		l = length(To_remove7)
		temp2 = which(R[start_ind:end_ind,4]==0)
		temp3 = indexes[temp2]
		R[temp3[To_remove7],4] = 4.2
		# if there are indices to remove, make a comment as such
		xx <-sprintf("In data set %g due to recrudesence, removing data points - %g ", i,To_remove7[1], To_remove7[l])		}
	}	
		
			
		
	 # remove data after the last zero if the time between surrounding positive is more than TimeDiffRecrud2 hours. 
	DD_zero = which(data_out==0)
	DD_nonzero = which(data_out!=0)
	n3 = length(DD_nonzero)
	
	n2 = length(data_in)
	if (length(DD_nonzero)>1){
			diff_times = data_in[DD_nonzero[n3]] -data_in[DD_nonzero[n3-1]]
		if ((length(DD_zero)>0) & (DD_nonzero[n3]==n2) & (diff_times > TimeDiffRecurrence)){ # if the last value of para is positive AND there # was a 0 somewhere, then we need to check distance between last two positives
			To_remove6 = seq(DD_zero[length(DD_zero)]+1, n2)
			data_in = data_in[-To_remove6]
			data_out = data_out[-To_remove6]
	 		l = length(To_remove6)
	 		temp2 = which(R[start_ind:end_ind,4]==0)
	 		temp3 = indexes[temp2]
	 		R[temp3[To_remove6],4] = 4.3
	 		# if there are indices to remove, make a comment as such
	 		xx <-sprintf("In data set %g due to recrudesence, removing data points - %g ", i,To_remove6[1], To_remove6[l])	 	}
	 }	
	
		
		#n2 = length(data_out)			
		## find the ratio of last two para values
		#ratio = data_out[n2]/data_out[n2-1]
		## if there are any points equal to 16
		#if (ratio > ratio_constant) & (){
		#	# remove the last data point ...
		#	data_in = data_in[-n2]
		#	data_out = data_out[-n2]
		#	temp2 = which(R[start_ind:end_ind,4]==0)
		#	temp3 = indexes[temp2]
		#	R[temp3[n2],4] = 4
		#	# if there are indices to remove, make a comment as such
		#	xx <-sprintf("In data set %g, removing data point number: %g", i, n2)
		#	}
			

		
		
		# Remove zeros at the end of the profile
		indt = which(data_out!=0) # find the last nonzero value
		n2 = length(data_in)
		if ((length(indt)>0) &((indt[length(indt)])<(n2-1))){
			# keep first 0 and remove the others after it
			To_remove7 = seq(indt[length(indt)]+2,n2)
			data_in = data_in[-To_remove7]
			data_out = data_out[-To_remove7]
			l = length(To_remove7)
			temp2 = which(R[start_ind:end_ind,4]==0)
			temp3 = indexes[temp2]
			R[temp3[To_remove7],4] = 5
			# if there are indices to remove, make a comment as such
			xx <-sprintf("In data set %g due to recrudesence, removing data points - %g ", i,To_remove7[1], To_remove7[l])	
			
			}

		
		# part (f): Identify tail 
		# find the data points which fall below Threshold2
		To_remove5 = c()

		indt = which(data_out>100)
		if (length(indt)>0){indt = seq(indt[length(indt)]+1, n2)}
		ptemp = data_out[indt]

		# find the unique values that are repeated that are less than 100
		ptemp2 = unique(ptemp)
		n2 = length(data_in)

		# if there are 2 values less than 100, and have a repeated value under 100 and if the values go under 100 and stay under 100-- there is a tail   & (length(which(exp(ldata_out)[indt[1]:n2]<100)) == (n2-indt[1]+1))
		if ((length(indt) > 1) &(length(ptemp)>length(ptemp2))){
			# initialisation 
			counts = array(0,c(1,length(ptemp2)))
			keep_excluding = 1
			# loop over the repeated values under 100, going in backwards order 
			for (m in c(length(ptemp2):1)){
				# set para_temp to be the m-th repeated value below 100
				para_temp = ptemp2[m]
				# find the values that are in the data set equal to the para_temp value
				indt2 = which(data_out[indt]== para_temp)
				counts[1,m] = length(indt2)
				# if there are more than one 
				if ((counts[1,m] >1) &(keep_excluding == 1) & (para_temp != 0)){
					To_remove5 = c(indt[indt2], To_remove5)
				}
				# we have come across an isolated value of a repeated value (ie not in order of given data), so we stop excluding data
				else if (counts[1,m]==1){
					#keep_excluding =0
				}
			}
			if (length(To_remove5)>0){
			# remove all the data points (except the first one)
			To_remove5 = seq(To_remove5[2], n2)
			# update where to store results - only where outliers have not been detected
			temp2 = which(R[start_ind:end_ind,4]==0)
			temp3 = indexes[temp2]
			# update the cleaning points indication vector
			R[temp3[To_remove5],4] = 3
			}
		}
	

		# if there are points to remove from the data set, do so now
		if (length(To_remove5)>0){
			data_in = data_in[-To_remove5]
			data_out = data_out[-To_remove5]
			l = length(To_remove5)
			# if there are indices to remove, make a comment as such
			xx <-sprintf("In data set %g due to tail, removing data points - %g ", i,To_remove5[1], To_remove5[l])
			}
		
		# Identify if the last point is 0 which should be replaced with the detection limit
		if (testing_tobit==FALSE){
    ind = which(data_out>0)	
		LastPos = ind[length(ind)]
		n2= length(data_out)
		if (LastPos < n2){
    	  	if (data_out[LastPos+1]==0){
    	  	 	# identify index where parasite value = 0 (and )to replace with detection limit
      			ind_DL = LastPos+1
            data_out[ind_DL] = DT
      				yy <-sprintf("Replacing a 0 parasite level (at data point %g) with detection level: data set %s", LastPos+1,code)      

       		}
       	}
    }
	   # Part (b) Repeated time/para points (with both the same time and parasitemia reported)
		# Repeated time (but different parasite levels reported) points 

		# initialisation of elements that need to be removed
		trepeats = c()
		To_remove = c()
		
		# find which times have been repeated and store them to be removed
		#for (j in c(1:n2)){
		#	ind = which(data_in == data_in[j])
		#	# only need to worry about repeats if there are two or more of the same entry
		#	if (length(ind)==2){
		#		# store the time which has been repeated
		#		trepeats = c(trepeats, j)
		#		if (data_out[ind[1]] == data_out[ind[2]]){
		#			# store the index to be removed because of the repeated time
		#			To_remove = c(To_remove, ind[1])
		#			}
		#		}	
		#	}
		# only remove unique indices
		To_remove = unique(To_remove)
					
		# if there are indices to remove because of repeated time points, make a comment as such
		if (length(To_remove)>0){
			xx <-sprintf("Warning! In data set %g, there are repeated time/parasite entries - we have removed one", i)
			}

		# We will only if the values are considered "weird" according to the outlier methodology - this will be picked up later when we run the outlier 				program (see part d below).
		# Part (c) unusual values - set threshold values
		PLowerThr = 0   # lower threshold value for parasitemia (in per microlitre)
		PUpperThr = 3e6   # upper threshold value for parasitemia (in per microlitre)
		TLowerThr = 0   # lower threshold value for time (in hours)
		#TUpperThr = 175   # upper threshold value for time (in hours)
	
		# remove values that don't abide by the threshold values we have set 

		To_remove = c(To_remove,c(which(data_out>PUpperThr), which(data_out<=PLowerThr), which(data_in<TLowerThr)))
			
		# Remove points that have been picked up during this process
		if (length(To_remove)>0){
			data_in = data_in[-To_remove]
			data_out = data_out[-To_remove]
			# update where to store results - only where outliers have not been detected
			temp2 = which(R[start_ind:end_ind,4]==0)
			temp3 = indexes[temp2]
			# update the cleaning points indication vector
			R[temp3[To_remove],4] = 1.0
			}
	
		# Part (d): outliers
		# take the logarithm of the data
		ldata_out = log(data_out)
		# set some tolerance levels for use in determinging outliers
		rel_tol1 = 10 
		rel_tol2 =5
		# calculate the slopes between data points for use in determining outliers
		slopes = diff(ldata_out)/diff(data_in)
		n2 = length(ldata_out) # this may be different to before if points were removed from the cleaning process above
		# calculate the normalised slopes between points, relative to the average rate of change between the last and first points
		norm_slopes = slopes/(-(ldata_out[n2]-ldata_out[1])/(data_in[n2]-data_in[1]))

		# update where to store results - only where outliers have not been detected
		temp2 = which(R[start_ind:end_ind,4]==0)
		temp3 = indexes[temp2]
					
		# If there are any changes which are substantial decreases in parasitaemia followed by increases, then there is an outlier
		To_remove2 = c()
		# for all the slopes in the data set
		if ((n2 >= 3) & ( (ldata_out[n2]-ldata_out[1]) !=0)){
		for (k in c(2:(n2-1))){
			# for time points after 12 hours
			if (data_in[k]>12){
				# if the criterion for removal is met, remove the kth point
				if (((norm_slopes[k-1]<(-0.75*rel_tol1)) & (norm_slopes[k]>rel_tol2)) | ((norm_slopes[k-1]<(-4*rel_tol1)) & (norm_slopes[k]>0.75*rel_tol2)) & (data_in[k-1]!=data_in[k])){
					# remove point k
					To_remove2  = c(To_remove2, k)
					# update the cleaning points indication vector
					R[temp3[k],4] = 2

					}
				}
			# for time points before 12 hours
			else if (data_in[k]<=12){
				# if the criterion for removal is met, remove the kth point
				if ((norm_slopes[k-1]<(-2*rel_tol1)) & (norm_slopes[k]>2*rel_tol1)& (data_in[k-1]!=data_in[k])){
					# remove point k
					To_remove2 = c(To_remove2, k)
					# update the cleaning points indication vector
					R[temp3[k],4] = 2

					}
				}			
			}
		}	
		# if there are points to remove from the data set, do so now
		if (length(To_remove2)>0){
			data_in = data_in[-To_remove2]
			data_out = data_out[-To_remove2]
			ldata_out = ldata_out[-To_remove2]
			l = length(To_remove2)
			#  if there are indices to remove, make a comment as such
			xx <-sprintf("In data set %g, removing data point number: %g", i, To_remove2[1:l])
			}

		# calculate the slopes between data points for use in determining outliers
		slopes = diff(ldata_out)/diff(data_in)
		n2 = length(ldata_out) # this may be different to before if points were removed from the cleaning process above
		# calculate the normalised slopes between points, relative to the average rate of change between the last and first points
		norm_slopes = slopes/(-(ldata_out[n2]-ldata_out[1])/(data_in[n2]-data_in[1]))

		# If after 12 hours, there are any changes which are substantial increases in parasitaemia followed by decreases, then there is a possible outlier
		ind2 = which(data_in>12)
		To_remove3 = c()
		n3 = length(ind2)
		# set some tolerance levels for use in determinging outliers
		rel_tol3 =2

	   	# update where to store results - only where outliers have not been detected
	   temp2 = which(R[start_ind:end_ind,4]==0)
    	temp3 = indexes[temp2]

		I1 = max(c(2, ind2[1]))
		# for all the slopes in the data set
		if (n3 > 1){
		for (k in (I1:ind2[n3-1])){
			# if the criterion for removal is met, remove the kth point
			if (((norm_slopes[k-1]>(rel_tol3)) & (norm_slopes[k]< (-rel_tol1))) | ((norm_slopes[k-1]>(rel_tol1)) & (norm_slopes[k]< (-rel_tol3)))  |((norm_slopes[k-1]>(0.5*rel_tol3)) & (norm_slopes[k]<(-2*rel_tol1))) | ((norm_slopes[k-1]>(5*rel_tol1)) & (norm_slopes[k]<(-0.2*rel_tol3)))& (data_in[k-1]!=data_in[k])){
				# remove point k
				To_remove3 = c(To_remove3, k)
				# update the cleaning points indication vector
				R[temp3[k],4] = 2
				}
			}			
		}	
			
			
		# if there are points to remove from the data set, do so now
		if (length(To_remove3)>0){
			data_in = data_in[-To_remove3]
			data_out = data_out[-To_remove3]
			ldata_out = ldata_out[-To_remove3]
			l = length(To_remove3)	
			#  if there are indices to remove, make a comment as such
			xx <-sprintf("In data set %g, removing data point number: %g", i, To_remove3[1:l])
			}
			

	# remove last point if its 'large' compared to second last one
	n2 = length(data_in)
  if (n2 > 1){
	if ((ldata_out[n2-1]<200) & (ldata_out[n2]>3*ldata_out[n2-1])& (data_out[n2]>100)){
		To_remove8 = n2
		data_in = data_in[-To_remove8]
		data_out = data_out[-To_remove8]
		temp2 = which(R[start_ind:end_ind,4]==0)
		temp3 = indexes[temp2]
		R[temp3[To_remove8],4] = 2.4
		# if there are indices to remove, make a comment as such
		xx <-sprintf("In data set %g due to last point being too large, removing data points - %g ", i,To_remove8[1])	
		}	
  }		
		# find the data points which are equal to 16	
		#indt = which(exp(ldata_out)<=16)	
		#To_remove4 = c()
		## if there are any points equal to 16
		#if (length(indt) >0){
		#	# loop over the data points equal to 16
		#	for (k in seq(1,length(indt))){
		#		temp_ind = indt[k]
		#		# if the next point after being equal to 16 is 'large' then remove the following data point
		#		if ((temp_ind+1) <= length(data_in)) {
		#			# if the criterion for removal is met, remove the temp_ind^th point
		#			if ((exp(ldata_out[temp_ind+1]) > 16*5) & ((data_in[temp_ind+1] - data_in[temp_ind]) < 10)){
		#				# remove point (temp_ind+1)
		#				To_remove4 = c(To_remove4, temp_ind+1)
		#				# update where to store results - only where outliers have not been detected
		#       temp2 = which(R[start_ind:end_ind,4]==0)
		#				temp3 = indexes[temp2]
		#				# update the cleaning points indication vector
		#				R[temp3[temp_ind+1],4] = 2
		#
		#				}
		#			}	
		#		}
		#	}
		
		
		# if there are points to remove from the data set, do so now
		#if (length(To_remove4)>0){
		#	data_in = data_in[-To_remove4]
		#	data_out = data_out[-To_remove4]
		#	ldata_out = ldata_out[-To_remove4]
		#	l = length(To_remove4)
		#	# if there are indices to remove, make a comment as such
		#	xx <-sprintf("In data set %g, removing data point number: %g", i, To_remove4[1:l])
		#	}


		# Identify if the last point is 0 which should be replaced with the detection limit
		if (length(ind_DL)>0){
    	  	if (R[start_ind+ind_DL-1,4]==0){
    	  	 	print(yy) 
       		}
       	else {
       		ind_DL = c()
       		}
       	}

	

		}	
	   if (all(R[start_ind:end_ind,4] == 0)){
	   		R[start_ind,58]=0 # all observations included
	   	}			
	   	else if ((any(floor(as.numeric(R[start_ind:end_ind,4])) == 2)) & (any(R[start_ind:end_ind,4] == 3))){
	   		R[start_ind,58]=3 # tail and outliers detected
	   		}
	   	else if (any(floor(as.numeric(R[start_ind:end_ind,4])) == 2)){
	   		R[start_ind,58]=1 # outliers detected
	   		}
	   	else if (any(R[start_ind:end_ind,4] == 3)){
	   		R[start_ind,58]=2 # tail detected
	   		}
	   	else {
	   		R[start_ind,58]=0 
	   		}
	   		
		# return the results (R) and other interesting information
		return(list(data_in = data_in, data_out = data_out, xx=xx, R=R, ind_DL = ind_DL))
	}