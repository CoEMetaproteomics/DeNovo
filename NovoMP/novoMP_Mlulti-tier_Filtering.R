setwd("D:/Feng/DDA_PASEF/20240605_Ultra_mouse_fractionas_Novor_Fasta_mapping_Mgnify/Ultra_Mouse_Feces_Fractions_Novor/Re-Process_Rep1_20240625")
library(readr)
library(ggplot2)
library(dplyr)
library(viridis)
library(ggpointdensity)
library(readxl)
library(ggExtra)
library(dplyr)
library(mgcv)
library(ggpubr)
library(gridExtra)
library(cowplot)
data8 <- read_csv("DeNovo_Peptides_F8.csv")
data7 <- read_csv("DeNovo_Peptides_F7.csv")
data6 <- read_csv("DeNovo_Peptides_F6.csv")
data5 <- read_csv("DeNovo_Peptides_F5.csv")
data4 <- read_csv("DeNovo_Peptides_F4.csv")
data3 <- read_csv("DeNovo_Peptides_F3.csv")
data2 <- read_csv("DeNovo_Peptides_F2.csv")
data1 <- read_csv("DeNovo_Peptides_F1.csv")

data1$sample_name <- "F1"
data2$sample_name <- "F2"
data3$sample_name <- "F3"
data4$sample_name <- "F4"
data5$sample_name <- "F5"
data6$sample_name <- "F6"
data7$sample_name <- "F7"
data8$sample_name <- "F8"
PSM_all <- rbind(data1,data2,data3,data4,data5,data6,data7,data8)
library(dplyr)

PSM_all_Unique <- PSM_all %>%
  group_by(stripped_peptide) %>%
  filter(denovo_score == max(denovo_score)) %>%
  distinct(stripped_peptide, .keep_all = TRUE) %>%
  ungroup()

write.csv(PSM_all_Unique, "PSM_all_Raw_Unique.csv", row.names = F)
## De novo peptides overlapped with the MGnify gut fasta in-silico
library(ggplot2)
Novor_PSM_Precision_Calculation_Overlapped_Peptides_Info <- read_csv("Novor_PSM_Precision_Calculation_Overlapped_Peptides_Info.csv")
Overlaplist <- list(Novor_PSM_Precision_Calculation_Overlapped_Peptides_Info$stripped_peptide)

Non_overlapped_peptide <- PSM_all_Unique[!PSM_all_Unique$stripped_peptide %in% unlist(Overlaplist),]

f1 <- ggplot() +  geom_histogram(data = PSM_all,aes(x=denovo_score, y=..count..),position = "dodge",fill="lightblue",colour="blue",alpha=0.7, linewidth=0.5,bins = 60,size=0.8) +scale_x_continuous(breaks = c(0,25,50,75,100))+
  theme_bw()+theme(legend.position = "top",legend.title = element_blank())+ theme(axis.text = element_text(size=12), panel.grid = element_blank())+ylab("Number of mapped peptides")

f2 <- ggplot() +  geom_histogram(data = PSM_all,aes(x=denovo_score, y=..count..),position = "dodge",fill="darkseagreen",colour="darkseagreen1",alpha=0.8, linewidth=0.5,bins = 60,size=0.8) +scale_x_continuous(breaks = c(0,25,50,75,100))+
  theme_bw()+theme(legend.position = "top",legend.title = element_blank())+ theme(axis.text = element_text(size=12), panel.grid.minor = element_blank())+ylab("Number of PSMs") + geom_vline(xintercept = 65, color="black", linetype=2,size=1)


## Filter 1: De novo score >=65 and charge state 1
PSM_all_1 <- PSM_all[PSM_all$denovo_score >= 65,]
count_PSM_1 <- as.data.frame(table(PSM_all_1$sample_name))

f3 <- ggplot() + 
  geom_bin2d(data = subset(PSM_all_1, charge == 2), aes(x=precursor_mz, y=ook0),fill="green4",alpha = 0.5,bins=100) +
  geom_bin2d(data = subset(PSM_all_1, charge == 3), aes(x=precursor_mz, y=ook0),fill="tan2",alpha = 0.6,bins=100) +
  geom_bin2d(data = subset(PSM_all_1, charge == 4), aes(x=precursor_mz, y=ook0),fill="firebrick",alpha = 0.7,bins=100)+
  geom_bin2d(data = subset(PSM_all_1, charge == 5), aes(x=precursor_mz, y=ook0),fill="violet",alpha = 0.8,bins=100)+
  geom_bin2d(data = subset(PSM_all_1, charge == 1), aes(x=precursor_mz, y=ook0),fill="black",alpha = 0.9,bins=100)+
  ylab("Ion mobility (1/K0)")+
  theme_bw()+theme(legend.position = "top")+ theme(axis.text = element_text(size=12), panel.grid.minor = element_blank())


count_charge <- as.data.frame(table(PSM_all_1$charge))
fix(count_charge)
f3 <- ggplot(count_charge, aes(x = Charge, y = Count)) +
  geom_point(aes(size = log10(Count)),fill = "darkseagreen", alpha=0.8,shape = 21, colour = "darkseagreen1") +
  scale_size(range = c(3, 9)) + # Adjust the range of dot sizes
  theme_bw()+ theme(legend.position = "top")+ theme(axis.text = element_text(size=12), panel.grid.minor = element_blank())+ 
  geom_vline(xintercept = 1, color="black", linetype=2,size=1)+ ylab("Number of PSMs")+ xlab("Charge state")+ scale_x_continuous(breaks = c(1,2,3,4,5))


PSM_all_2 <- PSM_all_1[PSM_all_1$charge >1,]
count_PSM_2 <- as.data.frame(table(PSM_all_2$sample_name))

write.csv(PSM_all_2, "DeNovo_PSMs_After_Filter1.csv", sep = ",", row.names = F)

## filter2: peptide length >= 7 amino acids

PSM_all_2$Length <- nchar(PSM_all_2$stripped_peptide)
count_length <- as.data.frame(table(PSM_all_2$Length))
fix(count_length)
f4 <- ggplot(count_length, aes(x = Length, y = Count)) +
  geom_point(aes(size = log10(Count)),fill = "darkseagreen", alpha=0.8,shape = 21, colour = "darkseagreen1") +
  scale_size(range = c(2, 7)) + # Adjust the range of dot sizes
  theme_bw()+ theme(legend.position = "top")+ theme(axis.text = element_text(size=12), panel.grid.minor = element_blank())+ 
  geom_vline(xintercept = 7, color="black", linetype=2,size=1)+ scale_x_continuous(breaks = c(4,7,10,15,20,25,30,36))+ ylab("Number of Peptides")+ xlab("Peptide length")

PSM_all_3 <- PSM_all_2[PSM_all_2$Length >= 7,]
count_PSM_3 <- as.data.frame(table(PSM_all_3$sample_name))
write.csv(PSM_all_3, "DeNovo_PSMs_After_Filter2.csv", sep = ",", row.names = F)
## filter 3: mass error. checking the Mass error distribution first
## calculate the 95% data under the curve
## PSM_all_0$ppm_zscore <- ((PSM_all_0$ppm_error)-mean(PSM_all_0$ppm_error))/sd(PSM_all_0$ppm_error)
lower_ppm_cutoff <- qnorm(0.025, mean = mean(PSM_all_3$ppm_error), sd=sd(PSM_all_3$ppm_error))
upper_ppm_cutoff <- qnorm(0.975, mean = mean(PSM_all_3$ppm_error), sd=sd(PSM_all_3$ppm_error))

f5 <- ggplot(PSM_all_3) +  geom_histogram(aes(x=ppm_error, y=..count..),bins = 80,fill="darkseagreen",colour="darkseagreen1",alpha=0.7, linewidth=0.5,size=0.8) + xlab("Mass shift (ppm)")+
  theme_bw()+theme(legend.position = "top", legend.title = element_blank()) + theme(axis.text = element_text(size=12), panel.grid.minor = element_blank())+
  geom_vline(xintercept = c(lower_ppm_cutoff,upper_ppm_cutoff), color="black", linetype=2,size=1) + ylab("Number of PSMs")

ggplot(PSM_all_3) + geom_bin2d(aes(y=ppm_error, x=denovo_score), bins=200) + scale_fill_continuous(type = "viridis")+
  theme_bw()+theme(legend.position = "top")+ theme(axis.text = element_text(size=12), panel.grid = element_blank())+geom_abline(slope = 0, intercept = c(lower_ppm_cutoff,upper_ppm_cutoff), color="black", linetype=2,size=1)


PSM_all_4 <- PSM_all_3[PSM_all_3$ppm_error <= upper_ppm_cutoff & PSM_all_3$ppm_error >= lower_ppm_cutoff,]
count_PSM_4 <- as.data.frame(table(PSM_all_4$sample_name))
write.csv(PSM_all_4, "DeNovo_PSMs_After_Filter3.csv", sep = ",", row.names = F)
## filter4 RT measured and prediction
write.csv(PSM_all_4, "DeNovo_PSMs_After_Filter1-3_for_RT_Prediction.csv", sep = ",", row.names = F)

## RT prediction was done using DeppLC. need to extract and format the necessary info from the csv file (BPSNovor_output_to_IM2Deep.py, Fragpipe_PSM_to_IM2Deep_for_Calibration_IncludingRT.py)
## import the prediction results 

Predicted_RT <- read_csv("DeNovo_PSMs_After_Filter1-3_for_RT_Prediction_Formatted_Input_For_DeepLC_deeplc_predictions.csv")

PSM_all_4$Observed_RT <- PSM_all_4$rt / 60
PSM_all_4$Predicted_RT <- Predicted_RT$`predicted retention time`

PSM_all_4$RT_shift <- PSM_all_4$Predicted_RT - PSM_all_4$Observed_RT

lower_RT_cutoff <- qnorm(0.025, mean = mean(PSM_all_4$RT_shift), sd=sd(PSM_all_4$RT_shift))
upper_RT_cutoff <- qnorm(0.975, mean = mean(PSM_all_4$RT_shift), sd=sd(PSM_all_4$RT_shift))
model_RT <- lm(Predicted_RT ~ Observed_RT, PSM_all_4)
intercept <- coef(model_RT)[1]
slope <- coef(model_RT)[2]

f7 <- ggplot(PSM_all_4) + geom_pointdensity(aes(y=Predicted_RT, x=Observed_RT), size=0.8) + 
  scale_color_viridis(alpha = 1,begin = 0.8,end = 0, option = "turbo")+ geom_smooth(aes(y=Predicted_RT, x=Observed_RT),method = "lm",se=T, colour="blue", linetype=2, size=0.8, alpha=0.5)+
  theme_bw()+theme(legend.position = "top")+ theme(axis.text = element_text(size=12), panel.grid.minor = element_blank())+xlab('Observed retention time (min)')+ylab('Predicted retention time (min)')
  
#geom_abline(intercept = upper_RT_cutoff+intercept, slope = slope, color = "black", linetype = 2, size=0.8) + geom_abline(intercept = lower_RT_cutoff+intercept, slope = slope, color = "black", linetype = 2, size=0.8) 

f8 <- ggplot(PSM_all_4) +  geom_histogram(aes(x=RT_shift, y=..count..),linewidth=0.5,bins = 80,fill="darkseagreen",colour="darkseagreen1",alpha=0.8, size=0.8) + scale_colour_brewer(palette="Paired")+ xlab("Retention time shift (min)")+ylab("Number of PSMs")+
  theme_bw()+theme(legend.position = "top", legend.title = element_blank()) + theme(axis.text = element_text(size=12), panel.grid.minor = element_blank())+ 
  geom_vline(xintercept = c(lower_RT_cutoff,upper_RT_cutoff), color="black", linetype=2,size=1) + xlab("Retention time shift (min)")
f78 <- ggarrange(f7, f8, widths = c(0.5, 0.5), heights = c(1,0.85))
f78

PSM_all_5 <- PSM_all_4[PSM_all_4$RT_shift <= upper_RT_cutoff & PSM_all_4$RT_shift >= lower_RT_cutoff,]
count_PSM_5 <- as.data.frame(table(PSM_all_5$sample_name))
write.csv(PSM_all_5, "DeNovo_PSMs_After_Filter4.csv", sep = ",", row.names = F)

write.csv(PSM_all_5, "DeNovo_PSMs_After_Filter1-4_for_CCS_Prediction.csv", sep = ",", row.names = F)

ggplot(PSM_all_5) + geom_pointdensity(aes(y=Predicted_RT, x=Observed_RT), size=0.8) + 
  scale_color_viridis(alpha = 1,begin = 0.8,end = 0, option = "turbo")+ geom_smooth(aes(y=Predicted_RT, x=Observed_RT),method = "lm",se=T, colour="blue", linetype=2, size=0.8, alpha=0.5)+
  theme_bw()+theme(legend.position = "top")+ theme(axis.text = element_text(size=12), panel.grid.minor = element_blank())+
  geom_abline(intercept = upper_RT_cutoff+intercept, slope = slope, color = "black", linetype = 2, size=0.8) + 
  geom_abline(intercept = lower_RT_cutoff+intercept, slope = slope, color = "black", linetype = 2, size=0.8) + xlab('Observed retention time (min)')+ylab('Predicted retention time (min)')

## filter5 CCS measured and prediction
Predicted_CCS <- read_csv("DeNovo_PSMs_After_Filter1-4_for_CCS_Prediction_Formatted_Converted_IM2Deep-predictions.csv")
Observed_CCS <- read.csv("DeNovo_PSMs_After_Filter1-4_for_CCS_Prediction_Formatted_Converted.csv")

PSM_all_5$Observed_CCS <- Observed_CCS$CCS
PSM_all_5$Predicted_CCS <- Predicted_CCS$`predicted CCS`
PSM_all_5$CCS_shift <- PSM_all_5$Predicted_CCS - PSM_all_5$Observed_CCS

f9 <- ggplot(PSM_all_5) + geom_pointdensity(aes(x=precursor_mz, y=Observed_CCS), size=0.8) + scale_color_viridis(alpha = 1,begin = 0.8,end = 0, option = "turbo")+
  theme_bw()+theme(legend.position = "top")+ theme(axis.text = element_text(size=12), panel.grid.minor = element_blank())#+geom_abline(slope = 0, intercept = c(lower_cutoff,upper_cutoff), color="black", linetype=2,size=1)

f10 <- ggplot(PSM_all_5) + geom_pointdensity(aes(x=precursor_mz, y=Predicted_CCS), size=0.8) + scale_color_viridis(alpha = 1,begin = 0.8,end = 0, option = "turbo")+
  theme_bw()+theme(legend.position = "top")+ theme(axis.text = element_text(size=12), panel.grid.minor = element_blank())#+geom_abline(slope = 0, intercept = c(lower_cutoff,upper_cutoff), color="black", linetype=2,size=1)


##PSM_all_5_1$CCS_shift_Zscore <- ((PSM_all_5_1$CCS_shift)-mean(PSM_all_5_1$CCS_shift))/sd(PSM_all_5_1$CCS_shift)

lower_CCS_cutoff <- qnorm(0.025, mean = mean(PSM_all_5$CCS_shift), sd=sd(PSM_all_5$CCS_shift))
upper_CCS_cutoff <- qnorm(0.975, mean = mean(PSM_all_5$CCS_shift), sd=sd(PSM_all_5$CCS_shift))


f11 <- ggplot(PSM_all_5) + geom_pointdensity(aes(y=Predicted_CCS, x=Observed_CCS), size=0.8) + 
  scale_color_viridis(alpha = 1,begin = 0.8,end = 0, option = "turbo")+ geom_smooth(aes(y=Predicted_CCS, x=Observed_CCS),method = "lm",se=T, colour="blue", linetype=2, size=0.8, alpha=0.5)+
  theme_bw()+theme(legend.position = "top")+ theme(axis.text = element_text(size=12), panel.grid.minor = element_blank())##+geom_text(x =400, y = 900,label = corr_eqn(y=CCS_Predction_Measured$Predicted_CCS, x=CCS_Predction_Measured$Measured_CCS), parse = T)
f12 <- ggplot(PSM_all_5) +  geom_histogram(aes(x=CCS_shift, y=..count..),linewidth=0.5,bins = 100,fill="darkseagreen",colour="darkseagreen1",alpha=0.8, size=0.8) + scale_colour_brewer(palette="Paired")+ xlab("CCS shift (Å)")+ylab("Number of PSMs")+
  theme_bw()+theme(legend.position = "top", legend.title = element_blank()) + theme(axis.text = element_text(size=12), panel.grid.minor = element_blank())+ 
  geom_vline(xintercept = c(lower_CCS_cutoff,upper_CCS_cutoff), color="black", linetype=2,size=1)
f1112 <- ggarrange(f11, f12, widths = c(0.5, 0.5), heights = c(1,0.85))
f1112

PSM_all_6 <- PSM_all_5[PSM_all_5$CCS_shift <= upper_CCS_cutoff & PSM_all_5$CCS_shift >= lower_CCS_cutoff,]
count_PSM_6 <- as.data.frame(table(PSM_all_6$sample_name))

write.csv(PSM_all_6, "Denovo_PSM_after_all_filters_final.csv", row.names=F)


## correlation betweent score and intensity
f13 <- ggplot() +geom_histogram(data=psm_PDB48970_quant,aes(x=log2(Intensity), y=..count..),position = "dodge",linewidth=0.5,bins = 80,fill="darkgray",colour="black",alpha=1, size=0.8) +
  geom_histogram(data=PSM_all_6,aes(x=log2(precursor_intensity), y=..count..),position = "dodge",linewidth=0.5,bins = 80,fill="darkseagreen",colour="darkseagreen1",alpha=0.8, size=0.8) + 
 xlab("Log2 Intensity")+ylab("Number of PSMs")+
  theme_bw()+theme(legend.position = "top", legend.title = element_blank()) + theme(axis.text = element_text(size=12), panel.grid.minor = element_blank())
  
