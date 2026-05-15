# Intro -----------------------------------------------------------

# Multiple Group Analysis
# Separated V.S. Non-separated Family 

# Load packages -----------------------------------------------------------
install.packages("usethis")
usethis::use_git
usethis::use_github()
install.packages(c("MatchIt", "cobalt", "tableone"))
install.packages("MplusAutomation")
library(MplusAutomation)
library(MatchIt)
library(cobalt)
library(tableone)
library(tidyverse)
library(eeptools)
library(naniar)
library(PerformanceAnalytics)
library(psych)
library(readxl)
library(parameters)
library(car)
library(MASS)
library(effectsize)
library(lmtest)
library(sfsmisc)
library(lm.beta)
library(sandwich)
library(lme4)
library(lmerTest)

# Set Working Directory ---------------------------------------------------

setwd("E:\\FRCD_Lab_UMassMedicalSchool\\R00\\Propensity_Score_Matching_R&Data\\Raw_Data")

# Multiple Group Analysis -----------------------------------------------

parent_child_data <- read_csv("NEW Merged Parent-Child Scored Data.csv")

# replace missing data to NA
parent_child_data <-  
  parent_child_data %>%
  replace_with_na_all(~.x %in% c(-999, 99)) # -999 missing data; 99: prefer not to answer


# Mean of sum score for depression and anxiety
parent_child_data <- parent_child_data %>%
  mutate(p_mentalhealth = (p_anx_sum + p_dep_sum)/2) 


df_parent_child_Data <- parent_child_data |> 
  dplyr::select(ParticipantID,
                c_anx_mean, 
                c_dep_mean, 
                lifequality_mean
  )


write.table(df_parent_child_Data, "df_parent_child_Data.txt", sep = "\t", row.names = FALSE)

#pbi_childseparation Fcattach_mean Mcattach_mean Psupport_mean peersup_mean Schoolclim_mean c_dep_mean c_anx_mean
#pbi_childseparation Fcattach_mean Mcattach_mean c_dep_mean

parent_child_data <- parent_child_data |> 
  select(pbi_childseparation, 
         Fcattach_mean,
         Mcattach_mean, 
         Psupport_mean, 
         peersup_mean,
         Schoolclim_mean, 
         c_dep_mean,
         c_anx_mean) |> 
  filter(!is.na(Fcattach_mean) & !is.na(Mcattach_mean) & !is.na(c_dep_mean) & !is.na(Psupport_mean)
         & !is.na(peersup_mean) & !is.na(Schoolclim_mean) & !is.na(c_anx_mean))|> 
  mutate(Fcattach_c = Fcattach_mean - mean(Fcattach_mean), 
         Mcattach_c = Mcattach_mean - mean(Mcattach_mean), 
         Psupport_c = Psupport_mean - mean(Psupport_mean), 
         peersup_c = peersup_mean - mean(peersup_mean), 
         Schoolclim_c = Schoolclim_mean - mean(Schoolclim_mean))|> 
  rename (sep = pbi_childseparation, 
          fa = Fcattach_c, 
          ma = Mcattach_c, 
          dep = c_dep_mean, 
          anx = c_anx_mean,
          ps = Psupport_c, 
          pes = peersup_c, 
          ss = Schoolclim_c) |> 
  mutate(faps = fa * ps, fapes = fa * pes, fass = fa * ss) |>  
  mutate (maps = ma * ps, mapes = ma * pes, mass = ma * ss)

MplusAutomation::prepareMplusData(parent_child_data,filename = "parent_child_data.dat")


# Other Models ------------------------------------------------------------


