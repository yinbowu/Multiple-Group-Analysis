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
install.packages("Hmisc") 
library("Hmisc")

# Set Working Directory ---------------------------------------------------

setwd("E:\\FRCD_Lab_UMassMedicalSchool\\R00\\Propensity_Score_Matching_R&Data\\Raw_Data")

# Multiple Group Analysis -----------------------------------------------

parent_child_data <- read_csv("NEW Merged Parent-Child Scored Data.csv")

matched_data <- read_csv("matched_full_R00Data.csv")

# replace missing data to NA
parent_child_data <-  
  parent_child_data %>%
  replace_with_na_all(~.x %in% c(-999, 99)) # -999 missing data; 99: prefer not to answer


# Mean of sum score for depression and anxiety
parent_child_data <- parent_child_data %>%
  mutate(p_mentalhealth = (p_anx_sum + p_dep_sum)/2) 

#pbi_childseparation Fcattach_mean Mcattach_mean Psupport_mean peersup_mean Schoolclim_mean c_dep_mean c_anx_mean
#pbi_childseparation Fcattach_mean Mcattach_mean c_dep_mean

# Flora's Model
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

# Separation Length

parent_child_data <- parent_child_data |>
  mutate(separation_length = psre_returnage - psre_separationage)

# Parental Negativity - New Variable
neg_new_parent_child_data <- parent_child_data |> 
  pivot_longer(cols = c (
    'pcnegative_01',
    'pcnegative_02',
    'pcnegative_03',
    'pcnegative_07',
    'pcnegative_08',
    'pcnegative_09',
    'pcnegative_10',
    'pcnegative_11',
    'pcnegative_12'),
    names_to = 'pc_negative',
    values_to = 'm_pcneg') |> 
  group_by(ParticipantID) |> 
  summarise(mean_neg = mean(m_pcneg, na.rm = TRUE)
  )


parent_child_data_neg <- left_join(parent_child_data, neg_new_parent_child_data, by = "ParticipantID")


# Used variables all included - full variables in a dataset

parent_child_data_negn <- parent_child_data_neg |> 
  rowwise()  |> 
  mutate(rel_s = # community religion-related support
           sum ( religion_01,
                 religion_02,
                 religion_03,
                 religion_04,
                 religion_05,
                 na.rm = TRUE
           )
  )|> 
  dplyr::select(P_ID = ParticipantID,
                F_ID = FamilyID,
                c_id = crecord_id,
                pbi_childseparation,
                separation_length,
                bi_sex,
                pChildage,
                pbi_education,
                p_mentalhealth,
                Psych_control_mean,
                Pcontrol_mean,
                mean_neg,
                bully_mean,
                Racism_mean,
                lifequality_mean,
                Mcattach_mean,
                Fcattach_mean,
                Psupport_mean,
                peersup_mean,
                Schoolclim_mean,
                Teachersup_mean,
                c_dep_mean,
                c_anx_mean,
                sleep_02,
                suicide_12,
                suicide_13,
                ER_cog_mean,
                ER_exp_mean,
                fare = pfr33_mean, # family resilience
                rel_s, # community religion support
                com_s = mhsupport_07, # community support
                coping_pe_mean,
                coping_pd_mean,
                coping_ee_mean,
                coping_ed_mean) |> 
  filter(!is.na(pbi_education) &!is.na(Psych_control_mean) & !is.na(Pcontrol_mean) & !is.na(bully_mean) & !is.na(Racism_mean)
         & !is.na(lifequality_mean) & !is.na(Mcattach_mean) & !is.na(Fcattach_mean) &
           !is.na(Psupport_mean) & !is.na(peersup_mean) & !is.na(Schoolclim_mean) & !is.na(Teachersup_mean)
         & !is.na(c_dep_mean) & !is.na(c_anx_mean) & !is.na(sleep_02)
         & !is.na(suicide_12) & !is.na(suicide_13) & !is.na(ER_cog_mean) &
           !is.na(ER_exp_mean) & !is.na(coping_pe_mean) & !is.na(coping_pd_mean) & !is.na(coping_ee_mean)
         & !is.na(coping_ed_mean)) |> 
  dplyr:: select(P_ID,
                 F_ID,
                 c_id,
                 c_sex =  bi_sex,
                 cage = pChildage,
                 edu = pbi_education,
                 p_mh = p_mentalhealth,
                 sep_s = pbi_childseparation,
                 s_len = separation_length,
                 psych_c = Psych_control_mean,
                 pcontrol = Pcontrol_mean,
                 pc_neg = mean_neg,
                 bully = bully_mean,
                 racism = Racism_mean,
                 life_qua = lifequality_mean,
                 mcatt = Mcattach_mean,
                 fcatt = Fcattach_mean,
                 p_supp = Psupport_mean,
                 peer_supp = peersup_mean,
                 sch_supp = Schoolclim_mean,
                 tea_supp = Teachersup_mean,
                 c_dep = c_dep_mean,
                 c_anx = c_anx_mean,
                 sleep = sleep_02,
                 sui_12 = suicide_12,
                 sui_13 = suicide_13,
                 er_cog = ER_cog_mean,
                 er_exp = ER_exp_mean,
                 cope_pe = coping_pe_mean,
                 cope_pd = coping_pd_mean,
                 cope_ee = coping_ee_mean,
                 cope_ed = coping_ed_mean) |> 
  mutate(bully_c = bully - mean(bully), 
         racism_c = racism - mean(racism), 
         fcatt_c = fcatt - mean(fcatt), 
         mcatt_c = mcatt - mean(mcatt), 
         p_supp_c = p_supp - mean(p_supp), 
         pe_sup_c = peer_supp - mean(peer_supp), 
         s_supp_c = sch_supp - mean(sch_supp),
         t_sup_c = tea_supp - mean(tea_supp))|> 
  mutate(bullyps = bully_c * p_supp_c, 
         bullypes = bully_c * pe_sup_c, 
         bullysc = bully_c * s_supp_c, 
         bullyts = bully_c * t_sup_c,
         bullyma = bully_c * mcatt_c, 
         bullyfa = bully_c * fcatt_c) |>  
  mutate (rac_ps = racism_c * p_supp_c, 
          rac_pes = racism_c * pe_sup_c, 
          rac_sc = racism_c * s_supp_c,
          rac_ts = racism_c * t_sup_c,
          rac_ma = racism_c * mcatt_c, 
          rac_fa = racism_c * fcatt_c) |> 
  mutate(psych_cc = psych_c - mean(psych_c), 
         pc_neg_c = pc_neg - mean(pc_neg), 
         er_cog_c = er_cog - mean(er_cog), 
         er_exp_c = er_exp - mean(er_exp))|> 
  mutate(conercog = psych_c * er_cog_c, 
         conerexp = psych_c * er_exp_c, 
         negercog = pc_neg_c  * er_cog_c, 
         negerexp = pc_neg_c  * er_exp_c) |> 
  mutate(#conercpe = psych_c * cope_pe, 
    #conerced = psych_c * cope_ed, 
    negercpe = pc_neg_c  * cope_pe, 
    negerced = pc_neg_c  * cope_ed)


parent_child_data_negn <- parent_child_data_negn |> # 1: boys; 2: girls
  mutate(female = ifelse(c_sex == 2, 1, 0))

MplusAutomation::prepareMplusData(parent_child_data_negn,filename = "parent_child_data_Mplus.dat")

# correlation

corr_data <- parent_child_data_negn |> 
  select(female,
         cage,
         p_mh,
         sep_s,
         s_len,
         psych_c,
         pc_neg, 
         er_cog, 
         er_exp,
         cope_pe,
         cope_ed,
         peer_supp, 
         sch_supp,
         tea_supp,
         c_dep,
         c_anx,
         life_qua) 


rcorr(as.matrix(corr_data))

## ParentalControl&SleepQuality ------------------------------------------------------------

parent_child_data <- parent_child_data |> 
  dplyr::select(pbi_childseparation,
                Psych_control_mean,
                Pcnegative_mean,
                c_dep_mean,
                c_anx_mean,
                sleep_02
  ) |> 
  rename (sep_s = pbi_childseparation, 
          psy_control = Psych_control_mean, 
          pcneg = Pcnegative_mean, 
          dep = c_dep_mean, 
          anx = c_anx_mean
  ) 


MplusAutomation::prepareMplusData(parent_child_data,filename = "ParentalControl_SleepQuality.dat")


## Bully&Depression, Anxiety ------------------------------------------------------------

parent_child_data <- parent_child_data |> 
  dplyr::select(pbi_childseparation,
                bully_mean,
                Racism_mean,
                lifequality_mean,
                Mcattach_mean,
                Fcattach_mean,
                Psupport_mean,
                peersup_mean,
                Schoolclim_mean,
                Teachersup_mean,
                c_dep_mean,
                c_anx_mean
  ) |> 
  mutate(bully_c = bully_mean - mean(bully_mean), 
         racism_c = Racism_mean - mean(Racism_mean), 
         Fattach_c = Fcattach_mean - mean(Fcattach_mean), 
         Mattach_c = Mcattach_mean - mean(Mcattach_mean), 
         parsup_c = Psupport_mean - mean(Psupport_mean), 
         peersup_c = peersup_mean - mean(peersup_mean), 
         schsup_c = Schoolclim_mean - mean(Schoolclim_mean),
         teasup_c = Teachersup_mean - mean(Teachersup_mean))|> 
  mutate(bullyps = bully_c * parsup_c, 
         bullypes = bully_c * peersup_c, 
         bullysc = bully_c * schsup_c, 
         bullyts = bully_c * teasup_c,
         bullyma = bully_c * Mattach_c, 
         bullyfa = bully_c * Fattach_c) |>  
  mutate (rac_ps = racism_c * parsup_c, 
          rac_pes = racism_c * peersup_c, 
          rac_sc = racism_c * schsup_c,
          rac_ts = racism_c * teasup_c,
          rac_ma = racism_c * Mattach_c, 
          rac_fa = racism_c * Fattach_c) |>
  select(sep_s = pbi_childseparation,
         bully = bully_mean,
         racism = Racism_mean,
         lifequa = lifequality_mean,
         # mattach = Mcattach_mean,
         # fattach = Fcattach_mean,
         # ps = Psupport_mean,
         pes = peersup_mean,
         ss = Schoolclim_mean,
         #  teasup = Teachersup_mean,
         dep = c_dep_mean,
         anx = c_anx_mean,
         # bullyps,
         bullypes,
         bullysc,
         #  bullyts,
         #  bullyma,
         #  bullyfa,
         #  rac_ps,
         rac_pes,
         rac_sc)
#  rac_ts,
#  rac_ma,
#  rac_fa)

MplusAutomation::prepareMplusData(parent_child_data,filename = "BullyDiscrimination_DepressionAnxiety.dat")

## Maltreatment & Suicide Thoughts ------------------------------------------------------------

parent_child_data <- parent_child_data |> 
  dplyr::select(pbi_childseparation,
                Psych_control_mean,
                Pcnegative_mean,
                suicide_12,
                suicide_13,
                ER_cog_mean,
                ER_exp_mean
  ) |> 
  mutate(Psycon_c = Psych_control_mean - mean(Psych_control_mean), 
         Pcneg_c = Pcnegative_mean - mean(Pcnegative_mean), 
         ER_cog_c = ER_cog_mean - mean(ER_cog_mean), 
         ER_exp_c = ER_exp_mean - mean(ER_exp_mean)
  )|> 
  mutate(conercog = Psycon_c * ER_cog_c, 
         conerexp = Psycon_c * ER_exp_c, 
         negercog = Pcneg_c  * ER_cog_c, 
         negerexp = Pcneg_c  * ER_exp_c
  ) |> 
  select(sep_s = pbi_childseparation,
         con = Psych_control_mean,
         Pcneg = Pcnegative_mean,
         sui12 = suicide_12,
         # sui13 = suicide_13,
         ERcog = ER_cog_mean,
         ERexp = ER_exp_mean,
         conercog,
         conerexp,
         negercog,
         negerexp)

MplusAutomation::prepareMplusData(parent_child_data,filename = "Maltreatment&SuicideThoughts.dat")


correlation <- parent_child_data %>% 
  dplyr::select(ER_cog_mean,
                ER_exp_mean
  )

pairs.panels(correlation[1:2], stars = TRUE) 


