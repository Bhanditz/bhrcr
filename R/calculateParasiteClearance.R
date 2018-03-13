###################################################################
#    Script code to fit a lag regression model to PCT data        #
###################################################################

calculateParasiteClearance <- function(Name, minpara, data1, imageType=c("pdf", "png"), appRun=TRUE)
{
# MPF 2015/02/25
# No need to require packages, attached by Depends
#require(MASS)
#require(AER)
par(ask = FALSE)
signif_figures = 3

if(appRun==FALSE){
	Name_mod = sprintf("%s.csv",Name)
	data1 <- read.csv(Name_mod, header = TRUE)
}else if(appRun==TRUE){
      #----------in PCE web, use this------------
	signif_figures = 3
    imageType <- match.arg(imageType)
	#get the filename without path
	temp = unlist(strsplit(Name,"/"))
	csvName = temp[length(temp)]
	headers = names(data1)
	outputFiles = c()
}

# Constants used
#*********maximum regression********************
fact = 0.25 # the factor by which the second para level must exceed the first (when we have only 4 data points) in order that a straight line model be fitted from the max para level. 0.25 = 25% increase

#********initial parasitaemia******************
Threshold2 = 1000 # The threshold for the first para level. If its below this then no model is fitted.

plotting = TRUE  # don't save a plot for each patient dataset

#*****Condition for lag regression*************
FirstHours = 24  # The time period that we want a certain number of observations in

NeededInFirstHours = 3 # The number of observations we need in FirstHour

MaxDiffInFirstHours = 14  # The max difference in the FirstHours, before which we don't fit a lag model

#***Conditions for last positive observation********************
Threshold3 = 1000 # The threshold for the last para level. If its above this then no model is fitted (since para is not cleared).

TimeRes1 = 30 # the time restriction for not fitting a model if the LastPos between (LowerThresh, UpperThresh)
TimeRes2 = 14 #t he time restriction for not fitting a model if the LastPos > UpperThresh
Threshold4 = 500 # the upper threshold for not fitting model

#*************Conditions for recurrence**********************
DaysRecurrence = 7 # we demove data after DaysRecurrence days as Recurrence
TimeDiffRecurrence = 36 # if a zero is in the dataset AND the last para is positive AND the difference between the last two postive values exceeds this value, we remove the data after the last 0

testing_tobit = FALSE # testing the use of tobit regression -whether its effect is important or not

Condition_on_24_hours = 1

# MPF 2015/02/25: Could be problem with not 
# sourcing these files if the constants 
# defined here (inside this function)
# are used inside the sourced functions without 
# calling them in the argument list
# Check with Colin about that.
#source("cleanData.R")
#source("fitCubicAndStorage.R")
#source("fitQuadraticAndStorage.R")
#source("fitLinearAndStorage.R")
#source("lagReg.R")
#source("fitCubicAndStorage_tobit.R")
#source("fitQuadraticAndStorage_tobit.R")
#source("fitLinearAndStorage_tobit.R")
#source("lagReg_tobit.R")

# Isolate the patient id,  time and para for the entire data set
unicode = as.character(data1[,1])
time = data1[,2]
para = data1[,3]
# Isolate the number of individual data sets we have in total
n = length(unique(unicode))
# Initialise the storate array
R = array(NaN, c(length(time),71))
# Initialise the outlier detection column to be zeros (assume all points are not outliers)
R[,4]<-0
U=unique(unicode)
start_indexes = array(0,n)
if(appRun==TRUE){
  #--- only in PCE web ----
  detectionLimitUsed = c()
}
# Loop over all the data sets
for (i in 1:n) {

	code = U[i]
	inds = which(unicode == code)
	# Set the start_index
    start_ind = inds[1]
    start_indexes[i] = start_ind
    # Set the end index
	end_ind = inds[length(inds)]
	# Initialise a temperary array. This is to ensure that the times are strictly increasing in values
	temp_data = array(NA, c(end_ind+1-start_ind,2))
	# Set the first column to be the time data
	temp_data[,1]  = time[start_ind:end_ind]

	# Set the second column to be the para data
	temp_data[,2] = para[start_ind:end_ind]
	# Reorder the data so that the time column is in assending order
	temp_data2 <- temp_data[order(temp_data[,1]), ]

	if (start_ind == end_ind){# Initialise data_in and data_out vectors
		data_in = temp_data2[1]
		data_out = temp_data2[2]}
	else if (start_ind != end_ind)	{# Initialise data_in and data_out vectors
		data_in = temp_data2[,1]
		data_out = temp_data2[,2]
		}

	# Store the patient id, time data, para data
 	R[start_ind:end_ind,1] = U[i]
	R[start_ind:end_ind,2] = data_in
	R[start_ind:end_ind,3] = data_out

	# find appropriate detection limit and DT value
	mincount = try(min(data_out[data_out!=0]))
	if (class(mincount)[1]=="try-error"){
		mincount = 0
		}

	if (minpara>mincount){
		Detec_limit = mincount }
	else {
		Detec_limit = minpara
		}
	DT= Detec_limit-0.01 #take the logarithm later
	if(appRun==TRUE){
  		#--- only in PCE web ----
  		detectionLimitUsed = c(detectionLimitUsed, Detec_limit)
	}

	R[start_ind,56]= 0
	if (Detec_limit!=0){R[start_ind,54] = log(DT)}
	# Clean the current data set
	Out = cleanData(data_in, data_out,i, start_ind, end_ind, R, code, Detec_limit, DT, Threshold2, DaysRecurrence, TimeDiffRecurrence,testing_tobit)
	# Re-set R to be the modified version from the cleaning process
	R = Out$R
	ind_DL = Out$ind_DL
	if (length(ind_DL) ==0 ){
	   # Fit a lag regression model
	   Out = lagReg(R, data_in, data_out, i, start_ind, end_ind, code, Condition_on_24_hours, fact, Threshold2, Threshold3, plotting, Name, FirstHours, MaxDiffInFirstHours, NeededInFirstHours, imageType)
       # Re-set R to be the modified version from the cleaning process
       R = Out$R
    }
   else {
  	# Fit a lag regression model (with tobit)
  	data_out[ind_DL] = DT
    MaxIters = 25000
	if(appRun==FALSE){
		Out = lagReg_tobit(R, data_in, data_out, i, start_ind, end_ind, code, Condition_on_24_hours, fact, Threshold2, Threshold3, Detec_limit, DT, MaxIters, plotting, Name, ind_DL,FirstHours, MaxDiffInFirstHours, NeededInFirstHours, TimeRes1, TimeRes2, Threshold4)
	}else if(appRun==TRUE){
		Out = lagReg_tobit(R, data_in, data_out, i, start_ind, end_ind, code, Condition_on_24_hours, fact, Threshold2, Threshold3, Detec_limit, DT, MaxIters, plotting, Name, ind_DL,FirstHours, MaxDiffInFirstHours, NeededInFirstHours, TimeRes1, TimeRes2, Threshold4, imageType)
	}
	   # Re-set R to be the modified version from the cleaning process}

    R = Out$R
   }
   if(appRun==TRUE){
	#---only in PCE web-----
	outputFiles = c(outputFiles, sprintf("%s_%s.%s", Name, code, imageType))
   }
   R[start_ind,68] = log(0.5)/as.numeric(R[start_ind,43])  # half life constant
}
if(appRun==TRUE){
  #---only in PCE web-----

  # Sean 2012 May 25 - write out user minpara and computed Detec_limit
  # https://stat.ethz.ch/pipermail/r-help/2007-October/142023.html
  detectionLimitFilename = "detection_limit.csv";
  write.table(data.frame(minpara, minDetectionLimitUsed = min(detectionLimitUsed), maxDetectionLimitUsed = max(detectionLimitUsed)),file=detectionLimitFilename,sep=",",quote=FALSE,row.names=FALSE)
  outputFiles = c(outputFiles, detectionLimitFilename)
}
Name_mod1 = sprintf("%s_lagR.csv",Name)
# Convert the R array to a .csv file
write.matrix(R, Name_mod1, sep =", ")
if(appRun==TRUE){
	#---only in PCE web-----
	outputFiles=c(outputFiles, Name_mod1)
}
# Import .csv file from previous line as a data file
data1 <- read.csv(Name_mod1, header = FALSE)

# Convert back to a csv file with the appropriate headings and entries separated by a tab
write.table(data1, Name_mod1, sep =", ", col.names = c("patientid", "time", "para", "outlier detection", "lpara", "cslope1","cseslope1", "cslope2","cseslope2","cslope3", "cseslope3","cr2","cpred", "cres","crstud", "crsta", "cp1", "cp2", "cp3", "caic", "lslope1", "lseslope1", "lr2", "lpred", "lres", "lrstud", "lrsta", "lp1", "laic", "qslope1", "qseslope1", "qslope2", "qseslope2", "qr2", "qpred", "qres", "qrstud", "qrsta", "qp1", "qp2", "qaic", "profile_type", "final_slope", "fr2", "fpred", "fres", "frstud", "frsta", "tlag", "lparachange", "max_change","used","SE_clearance_rate", "detection limit ", "MSRE", "Model not fitted","estimation summary","excluded observations", "PC50","PC90","PC95","PC99","PC50 adjusted","PC90 adjusted","PC95 adjusted","PC99 adjusted", "Intercept at tlag","half life", "ratio_slopes", "R2 - full linear", "linear model intercept,"),quote = FALSE, row.names=FALSE, na = "NaN")

R2 = array(NaN, c(length(time),24))

R2[,1] = R[,1]
# col 1 = patient id (R col 1)

R2[,2] = R[,2]
# col 2 = times (R col 2)

R2[,3] = R[,3]
# col 3 = para (R col 3)

R2[,4] = R[,5]
# col 4 =log(para) (R col 5)

R2[,5] = R[,54]
# col 5 = detection limit (R col 54)

R2[,6] = floor(as.numeric(R[,4]))
# col 6 = outliers (R col 4)

R2[,7] = R[,57]
# col 7 = estimation summary (R col 57)

R2[,8] = floor(as.numeric(R[,56]))
# col 8 = if no fit is made, indication of why (R col 56)

R2[,9] = R[,49]
# col 9 = lag time (R col 49)

R2[,10] = as.character(-1*(as.numeric(R[,43])))
# col 10 = clearance rate (R col 43). x by -1 to get 'clearance rate'

R2[,11] = R[,53]
# col 11 = standard error of clearance rate estimate (R col 53)

R2[,12] = R[,67]
# col 12 = intercept at tlag (R col 67)

R2[,13] = R[,68]
# col 13 = half life (R col 68)

R2[,14]= R[,44]
# col 14 = final R2 value for final fitted model (R col 44)

R2[,15] = R[,59]
# col 15 = PC50 (R col 59)

R2[,16] = R[,60]
# col 16 = PC90 (R col 60)

R2[,17] = R[,61]
# col 17 = PC95 (R col 61)

R2[,18] = R[,62]
# col 18 = PC99 (R col 62)

R2[,19] = R[,45]
# col 19 = final model predicted values (R col 45)

R2[,20] = as.character(-1*(as.numeric(R[,21])))
# col 20 = linear model rate (R col 21)

R2[,21] = R[,22]
# col 21 = standard error of linear model rate (R col 22)

R2[,22] = R[,71]
# col 22 = intercept for linear model (R col 71)

R2[,23] = R[,70]
# col 23 = R2 for 'full' linear model (R col 70) with zeros excluded

R2[,24] = R[,24]
# col 24 = predicted values for linear model (R col 24)

Name_mod1 = sprintf("%s_Estimates.csv",Name)
write.matrix(R2, Name_mod1, sep =", ")
if(appRun==TRUE){
	#---only in PCE web---
	outputFiles=c(outputFiles, Name_mod1)
}
# Import .csv file from previous line as a data file
data1 <- read.csv(Name_mod1, header = FALSE)
# Convert back to a csv file with the appropriate headings and entries separated by a tab
if(appRun==FALSE){
	write.table(data1, Name_mod1, sep =", ", col.names = c("Id", "Time", "Para", "Lpara", "Detect", "Outlier", "Estimation_summary","No_fit", "Tlag", "Clearance_rate_constant", "SE_clearance", "Intercept_tlag", "Slope_half_life", "R2", "PC50", "PC90", "PC95", "PC99","Predicted", "Linear_slope","SE_linear_slope", "Intercept_linear", "R2_linear", "Predicted_linear,"),quote = FALSE, row.names=FALSE, na = "NaN")
}else if(appRun==TRUE){
	#---- only in PCE web ----
	write.table(data1, Name_mod1, sep =", ", col.names = c(headers[1], headers[2], headers[3], "Lpara", "Detect", "Outlier", "Estimation_summary","No_fit", "Tlag", "Clearance_rate_constant", "SE_clearance", "Intercept_tlag", "Slope_half_life", "R2", "PC50", "PC90", "PC95", "PC99","Predicted", "Linear_slope","SE_linear_slope", "Intercept_linear", "R2_linear", "Predicted_linear,"),quote = FALSE, row.names=FALSE, na = "NaN")
}
R3 = array(NaN, c(n,14))

R3[,1] = R[start_indexes,1]
# col 1 = patient id (R col 1)

R3[,2]= R[start_indexes,57]
# col 2 = estimation summary (R col 57)

R3[,3]= R[start_indexes,58]
# col 3 = excluded observations (R col 58)

R3[,4]= R[start_indexes,56]
# col 4 = no fit (R col 56)

R3[,5]= R[start_indexes,49]
# col 5 = Tlag (R col 49)

R3[,6]= as.character(-1*(as.numeric(R[start_indexes,43])))
# col 6 = clearance rate constant (R col 43)

R3[,7]= R[start_indexes,53]
# col 7 = standard error of clearance rate estimate (R col 53)

R3[,8]= R[start_indexes,67]
# col 8 = Intercept at tlag (R col 67)

R3[,9]= R[start_indexes,68]
# col 9 = slope half life (R col 68)

R3[,10]= R[start_indexes,44]
# col 10 = final R2 value for final fitted model (R col 44)

R3[,11]= R[start_indexes,59]
# col 11 = PC50 (R col 59)

R3[,12]= R[start_indexes,60]
# col 12 = PC90 (R col 60)

R3[,13]= R[start_indexes,61]
# col 13 = PC95 (R col 61)

R3[,14]= R[start_indexes,62]
# col 14 = PC99 (R col 62)




Name_mod1 = sprintf("%s_Estimates_short.csv",Name)
write.matrix(R3, Name_mod1, sep =", ")
if(appRun==TRUE){
	#---only in PCE web----
	outputFiles=c(outputFiles, Name_mod1)
}
# Import .csv file from previous line as a data file
data1 <- read.csv(Name_mod1, header = FALSE)
# Convert back to a csv file with the appropriate headings and entries separated by a tab
if(appRun==FALSE){
	write.table(data1, Name_mod1, sep =", ", col.names = c("Id", "Estimation_summary","Excluded_observations", "No_fit", "Tlag", "Clearance_rate_constant", "SE_clearance", "Intercept_tlag ", "Slope_half_life","R2",  "PC50","PC90","PC95", "PC99,"),quote = FALSE, row.names=FALSE, na = "NaN")
}else if(appRun==TRUE){
	#---only in PCE web----
	write.table(data1, Name_mod1, sep =", ", col.names = c(headers[1], "Estimation_summary","Excluded_observations", "No_fit", "Tlag", "Clearance_rate_constant", "SE_clearance", "Intercept_tlag ", "Slope_half_life","R2",  "PC50","PC90","PC95", "PC99,"),quote = FALSE, row.names=FALSE, na = "NaN")
}

CR = as.numeric(R3[,6])  # 6th column is clearance rate constant
CR2 = CR[which(CR != "NaN")]
S = summary(CR2)
maxCR = S[6][[1]]

breaks = 0.05*seq(floor(min(CR2/0.05)),ceiling(max(CR2/0.05)))
if(appRun==FALSE){
	filename <-sprintf("%s_Histogram.png", Name)
	png(filename)
}else if(appRun==TRUE){
	#---only in PCE web---
	filename <-sprintf("%s_Histogram.%s", Name, imageType)
	if (imageType == "png") {
     # MPF: 2015/02/25
     # no need for library, Cairo attached by Depends 
     #library(Cairo)
      Cairo(file=filename, type="png")
  	} else if (imageType == "pdf") {
    		pdf(file=filename)
  	}
}
par(mar = c(4.1, 5.1, 4.1, 2.1))
r <- hist(CR2, breaks = breaks,  col='blue1',main = "Distribution of clearance rate constant", xlab = "Clearance rate constant (/hour)", ylab = "Frequency (N)", cex.lab=1.5, cex.axis=1.5, cex = 1.5, cex.main=1.5)

text(r$mids, r$counts, r$counts, adj=c(.5, -.5), col='blue3')
lines(r, lty = 3, border = "purple", cex = 2) # -> lines.histogram(*)
dev.off()
if(appRun==TRUE){
	#---only in PCE web----
	outputFiles=c(outputFiles, filename)
}


breaks = 0.05*seq(floor(min(CR2/0.05)),ceiling(max(CR2/0.05)))
if(appRun==FALSE){
	filename <-sprintf("%s_Histogram2.png", Name)
	png(filename)
}else if(appRun==TRUE){
	#------only in PCE web----------
	filename <-sprintf("%s_Histogram2.%s", Name, imageType)
	if (imageType == "png") {
      #MPF 2015/02/25 do not need to attach Cairo, 
      # included by Depends
    	#library(Cairo)
    	Cairo(file=filename, type="png")
  	} else if (imageType == "pdf") {
    		pdf(file=filename)
 	}
}
par(mar = c(4.1, 5.1, 4.1, 2.1))
r <- hist(CR2, breaks = breaks,freq = FALSE, col='blue1',main = "Distribution of clearance rate constant", xlab = "Clearance rate constant (/hour)", ylab = "Density", cex.lab=1.5, cex.axis=1.5, cex = 1, cex.main=1.5)
dev.off()
if(appRun==TRUE){
 outputFiles=c(outputFiles, filename)
}

InDataSet = length(R3[,6]) # 6th column is clearance rate constant
Analysed = InDataSet - length(which(as.numeric(R3[,6])=="NaN")) # 6th column is clearance rate constant
WithLagPhase = Analysed - length(which(as.numeric(R3[,5])==0)) # 5th column is tlag
WithTail = length(which(as.numeric(R3[,3])==2)) + length(which(as.numeric(R3[,3])==3))# 3rd column is excluded observations
WithOutliers = length(which(as.numeric(R3[,3])==1)) + length(which(as.numeric(R3[,3])==3))# 3rd column is excluded observations

Report_output0 = c(InDataSet, Analysed, WithLagPhase, WithTail,WithOutliers )

write.table(t(Report_output0),file="Report_output0.csv",sep=",",
quote=FALSE,row.names=FALSE,col.names=FALSE)
if(appRun==TRUE){
	#---------in PCE web-------0
	outputFiles=c(outputFiles, "Report_output0.csv")
}

CR = as.numeric(R3[,6])# 6th column is clearance rate constant
CR2 = CR[which(CR != "NaN")]
S = summary(CR2)
maxCR = S[6][[1]]

Report_output = array(0,5 + length(breaks)-1)
Report_output[1:5] = c(S[3][[1]], S[1][[1]],S[6][[1]], S[2][[1]], S[5][[1]])
Report_output[6:length(Report_output)] = r$counts

write.table(t(Report_output),file="Report_output.csv",sep=",",
quote=FALSE,row.names=FALSE,col.names=FALSE)
if(appRun==TRUE){
	outputFiles=c(outputFiles, "Report_output.csv")
}

Percents = r$counts/length(CR2)*100
CumalSum = cumsum(Percents)

write.table(t(Percents),append=TRUE,file="Report_output.csv",
sep=",",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(t(CumalSum),append=TRUE,file="Report_output.csv"
,sep=",",quote=FALSE,row.names=FALSE,col.names=FALSE)

LT = as.numeric(R3[,5])  # 5th column is tlag
LT2 = LT[which(LT > 0)]
S_LT = summary(LT2)

Report_output2 = array(0,15)
if(length(LT2)>0){
	Report_output2[1:5] = c(S_LT[3][[1]], S_LT[1][[1]],S_LT[6][[1]], S_LT[2][[1]], S_LT[5][[1]])
	}
#Report_output2[6:10] = c(S_MC[3][[1]], S_MC[1][[1]],S_MC[6][[1]], S_MC[2][[1]], S_MC[5][[1]])
#Report_output2[11:15] = c(S_CP[3][[1]], S_CP[1][[1]],S_CP[6][[1]], S_CP[2][[1]], S_CP[5][[1]])


write.table(t(Report_output2),file="Report_output2.csv",sep=",",
quote=FALSE,row.names=FALSE,col.names=FALSE)
if(appRun==TRUE){
	#---ONLY IN pce WEB------
	outputFiles=c(outputFiles, "Report_output2.csv")
	#length of LT2
	write.table(t(length(LT2)), append=TRUE, file="Report_output2.csv",sep=",",
							quote=FALSE,row.names=FALSE,col.names=FALSE)
}


Report_output3 = c(fact, Threshold2, FirstHours, NeededInFirstHours, MaxDiffInFirstHours, Threshold3, TimeRes1, TimeRes2, Threshold4)
write.table(t(Report_output3),file="Report_output3.csv",sep=",",
quote=FALSE,row.names=FALSE,col.names=FALSE)
if(appRun==TRUE){
	#---ONLY IN pce WEB------
	outputFiles=c(outputFiles, "Report_output3.csv")
}

T1 = summary(as.numeric(R[start_indexes,59])[R[start_indexes,59]!="NaN"])  # PC50
T2 = summary(as.numeric(R[start_indexes,60])[R[start_indexes,60]!="NaN"])  # PC90
T3 = summary(as.numeric(R[start_indexes,61])[R[start_indexes,61]!="NaN"])  # PC95
T4 = summary(as.numeric(R[start_indexes,62])[R[start_indexes,62]!="NaN"])  # PC99
T5 = summary(as.numeric(R[start_indexes,63])[R[start_indexes,63]!="NaN"])  # PC50 adjusted
T6 = summary(as.numeric(R[start_indexes,64])[R[start_indexes,64]!="NaN"])  # PC90 adjusted
T7 = summary(as.numeric(R[start_indexes,65])[R[start_indexes,65]!="NaN"])  # PC95 adjusted
T8 = summary(as.numeric(R[start_indexes,66])[R[start_indexes,66]!="NaN"])  # PC99 adjusted
T9 = summary(as.numeric(R[start_indexes,68])[R[start_indexes,68]!="NaN"])  # half life


Report_output4 = array(0,c(9,5))
Report_output4[1,1:5] = c(T1[3][[1]], T1[1][[1]],T1[6][[1]], T1[2][[1]], T1[5][[1]])
Report_output4[2,1:5] = c(T2[3][[1]], T2[1][[1]],T2[6][[1]], T2[2][[1]], T2[5][[1]])
Report_output4[3,1:5] = c(T3[3][[1]], T3[1][[1]],T3[6][[1]], T3[2][[1]], T3[5][[1]])
Report_output4[4,1:5] = c(T4[3][[1]], T4[1][[1]],T4[6][[1]], T4[2][[1]], T4[5][[1]])
Report_output4[5,1:5] = c(T5[3][[1]], T5[1][[1]],T5[6][[1]], T5[2][[1]], T5[5][[1]])
Report_output4[6,1:5] = c(T6[3][[1]], T6[1][[1]],T6[6][[1]], T6[2][[1]], T6[5][[1]])
Report_output4[7,1:5] = c(T7[3][[1]], T7[1][[1]],T7[6][[1]], T7[2][[1]], T7[5][[1]])
Report_output4[8,1:5] = c(T8[3][[1]], T8[1][[1]],T8[6][[1]], T8[2][[1]], T8[5][[1]])
Report_output4[9,1:5] = c(T9[3][[1]], T9[1][[1]],T9[6][[1]], T9[2][[1]], T9[5][[1]])


write.table(t(Report_output4),file="Report_output4.csv",sep=",",
quote=FALSE,row.names=FALSE,col.names=FALSE)
if(appRun==TRUE){
	#---ONLY IN pce WEB------
	outputFiles=c(outputFiles, "Report_output4.csv")
}


T10 = as.numeric(R[start_indexes,56])[R[start_indexes,56]!="NaN"]

Y1 = which(T10 == 1)

write.table(t(U[Y1]),file="Report_output5.csv",sep=",",
quote=FALSE,row.names=FALSE,col.names=FALSE)
if(appRun==TRUE){
	#---ONLY IN pce WEB------
	outputFiles=c(outputFiles, "Report_output5.csv")
}

Y2 = which(T10 == 2)

write.table(t(U[Y2]),file="Report_output6.csv",sep=",",
quote=FALSE,row.names=FALSE,col.names=FALSE)
if(appRun==TRUE){
	#---ONLY IN pce WEB------
	outputFiles=c(outputFiles, "Report_output6.csv")
}

Y3 = which(T10 == 3.1)

write.table(t(U[Y3]),file="Report_output7.csv",sep=",",
quote=FALSE,row.names=FALSE,col.names=FALSE)
if(appRun==TRUE){
	#---ONLY IN pce WEB------
	outputFiles=c(outputFiles, "Report_output7.csv")
}

Y4 = which(T10 == 3.2)

write.table(t(U[Y4]),file="Report_output8.csv",sep=",",
quote=FALSE,row.names=FALSE,col.names=FALSE)
if(appRun==TRUE){
	#---ONLY IN pce WEB------
	outputFiles=c(outputFiles, "Report_output8.csv")
}

Y5 = which(T10 == 3.3)

write.table(t(U[Y5]),file="Report_output9.csv",sep=",",
quote=FALSE,row.names=FALSE,col.names=FALSE)
if(appRun==TRUE){
	#---ONLY IN pce WEB------
	outputFiles=c(outputFiles, "Report_output9.csv")
}



CR = as.numeric(R[start_indexes,68])
CR2 = CR[which(CR != "NaN")]
S = summary(CR2)
maxCR = S[6][[1]]

breaks = seq(floor(min(CR2)),ceiling(max(CR2)))

if(appRun==FALSE){
	filename <-sprintf("%s_Histogram3.png", Name)
	png(filename)
}else if(appRun==TRUE){
	# ---in PCE web---
	filename <-sprintf("%s_Histogram3.%s", Name, imageType)
	if (imageType == "png") {
     # MPF 2015/02/25 Do not need to attach Cairo
     # included by Depends
     # library(Cairo)
      Cairo(file=filename, type="png")
      } else if (imageType == "pdf") {
        pdf(file=filename)
      }
}
par(mar = c(4.1, 5.1, 4.1, 2.1))
r <- hist(CR2, breaks = breaks,  col='blue1',main = "Distribution of slope half life", xlab = "Half life (hours)", ylab = "Frequency (N)", cex.lab=1.5, cex.axis=1.5, cex = 1.5, cex.main=1.5)

text(r$mids, r$counts, r$counts, adj=c(.5, -.5), col='blue3')
lines(r, lty = 3, border = "purple", cex = 2) # -> lines.histogram(*)
dev.off()
if(appRun==TRUE){
	# ---in PCE web---
	outputFiles=c(outputFiles, filename)
}

breaks = seq(floor(min(CR2)),ceiling(max(CR2)))
if(appRun==FALSE){
	filename <-sprintf("%s_Histogram4.png", Name)
	png(filename)
}else if(appRun==TRUE){
	# ---in PCE web---
	filename <-sprintf("%s_Histogram4.%s", Name, imageType)
	if (imageType == "png") {
     #MPF 2015/02/25 no need to attach Cairo
     # included by Depends
     #library(Cairo)
      Cairo(file=filename, type="png")
      } else if (imageType == "pdf") {
        pdf(file=filename)
      }
}
par(mar = c(4.1, 5.1, 4.1, 2.1))
r <- hist(CR2, breaks = breaks,freq = FALSE, col='blue1',main = "Distribution of slope half life", xlab = "Half life (hours)", ylab = "Density", cex.lab=1.5, cex.axis=1.5, cex = 1.5, cex.main=1.5)
dev.off()
if(appRun==TRUE){
	# ---in PCE web---
	outputFiles=c(outputFiles, filename)
}


Report_output10 = array(0,c(4,length(breaks)-1))
Report_output10[1,1:dim(Report_output10)[2]] = r$breaks[2:length(r$breaks)]
Report_output10[2,1:dim(Report_output10)[2]] = r$counts

Percents_new = r$counts/length(CR2)*100
CumalSum_new = cumsum(Percents_new)

Report_output10[3,1:dim(Report_output10)[2]] = Percents_new
Report_output10[4,1:dim(Report_output10)[2]] = CumalSum_new


write.table(Report_output10,file="Report_output10.csv",sep=",",
quote=FALSE,row.names=FALSE,col.names=FALSE)
if(appRun==TRUE){
	# ---in PCE web---
	outputFiles=c(outputFiles, "Report_output10.csv")
  	outputFiles
}


}


