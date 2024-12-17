library(readr)
library(tidyverse)

# Define your own location:
# dataloc <- XXXX

# Load all data from StudentLife dataset study 3
# Retrievable from Kaggle: https://www.kaggle.com/datasets/subigyanepal/college-experience-dataset
# Info about StudentLife: https://studentlife.cs.dartmouth.edu/dataset.html

data <- list()
setwd(dataloc)
for (folder in list.dirs(getwd())){
  setwd(folder)
  for(file in list.files()){
    print(file)
    if(str_sub(file, -3) == "csv"){
      filename <- gsub(".csv","",gsub(" ","_",gsub("\\)","",gsub("\\(","",tolower(file)))))
      data[[filename]] <- read_csv(file,show_col_types = F)
      print("file loaded")
    }else{
      print("no csv")
    }
  }
}

# Recode hased uid to participant id starting at 1, 2 ... n
uids <- tibble("uid" = unique(data$general_ema$uid)) %>% mutate(pid = row_number())

# EMA: ecological momentary assessment data
# Recode valence and arousal of PAM for optional use (not used in thesis)
data_ema <- data$general_ema %>%
  left_join(uids, by = "uid") %>%
  select(uid, pid, everything()) %>%
  mutate(mood_valence = ifelse(pam %in% c(1,2,5,6),1,
                               ifelse(pam %in% c(3,4,7,8),2,
                                      ifelse(pam %in% c(9,10,13,14),3,
                                             ifelse(pam %in% c(11,12,15,16),4, NA)))),
         mood_arousal = ifelse(pam %in% c(1,3,9,11),1,
                               ifelse(pam %in% c(2,4,10,12),2,
                                      ifelse(pam %in% c(5,7,13,15),3,
                                             ifelse(pam %in% c(6,8,14,16),4, NA)))))

# Further datasets of ema-data for data exploration of missingness between observations
data_eda1 <- data_ema %>%
  filter(!is.na(pam)) %>%
  group_by(pid) %>%
  mutate(date = as.Date(as.character(day), format = "%Y%m%d"),
         next_date = lead(date),
         difference = next_date - date) %>%
  select(uid:day, date:difference, everything())

data_eda2 <- data_eda1 %>%
  group_by(difference) %>%
  summarize(aantal = sum(n()))

# Recoding a feature dataset of steps
X_steps <- data$steps %>%
  mutate(feat_steps_nighttime = rowSums(select(.,step_hr_23,step_hr_0:step_hr_7)),
         feat_steps_daytime = rowSums(select(.,step_hr_8,step_hr_0:step_hr_22)))

# Recoding main feature dataset

#####
# 1 #
#####

# Recode hourly data to part-of-day values. Parts of day are:
# morning = 7 - 11h
# daytime = 11h - 20h
# evening = 20h-24h
# night = 24h-6h

X_features_hourly <- data$sensing %>% select(contains("hr"))
X_features_nothourly <- data$sensing %>% select(!names(X_features_hourly)) %>% select(!c(uid,day,is_ios))
basenames <- unique(sub("_hr_.*","",names(X_features_hourly)))

dayframes <- list()
dayframes$morning <- paste0("hr_",c(7:10))
dayframes$daytime <- paste0("hr_",c(11:19))
dayframes$evening <- paste0("hr_",c(20:23))
dayframes$night   <- paste0("hr_",c(0:6))

for(basename in basenames){
  for(k in c(1:length(dayframes))){
    X_features_hourly <- X_features_hourly %>% mutate(!!sym(paste(basename,names(dayframes)[k],sep="_")) := rowSums(select(.,all_of((paste(basename,dayframes[[k]],sep="_")))), na.rm=T))
  }
}
X_features_hourly_selected <- X_features_hourly %>% select(contains(names(dayframes)))

log_vars  <- c("act_on_bike","act_running","audio_amp_mean","audio_convo_num","audio_voice","loc_dist","loc_max_dis_from_campus")
nonlog_var <- c("act_still","act_walking","audio_amp_std","audio_convo_duration_evening","loc_visit_num","other_playing_duration","other_playing_num","unlock_duration","unlock_num","act_in_vehicle")
drop_vars <- c("act_on_foot","act_tilting","act_unknown")

log_vars_daily <- c("call_in_duration_ep_0_log","call_out_duration_ep_0_log","light_mean_ep_1_log","light_mean_ep_2_log","light_mean_ep_3_log","light_std_ep_1_log","light_std_ep_2_log",
                    "light_std_ep_3_log","loc_food_audio_amp_log","loc_food_audio_amp_log","loc_home_audio_amp_log","loc_other_dorm_audio_amp_log","loc_self_dorm_audio_amp_log",
                    "loc_social_audio_amp_log","loc_study_audio_amp_log")
nonlog_vars_daily <- c("call_in_num_ep_0","call_miss_num_ep_0","call_out_num_ep_0","loc_food_convo_duration","loc_food_convo_num","loc_food_dur","loc_food_still","loc_health_dur",
                       "loc_home_convo_duration","loc_home_dur","loc_home_still","loc_leisure_dur","loc_other_dorm_convo_duration","loc_other_dorm_convo_num","loc_other_dorm_dur",
                       "loc_other_dorm_still","loc_self_dorm_convo_duration","loc_self_dorm_convo_num","loc_self_dorm_dur","loc_self_dorm_still","loc_social_convo_duration","loc_social_convo_num",
                       "loc_social_dur","loc_social_still","loc_study_convo_duration","loc_study_dur","loc_study_still","loc_workout_dur","loc_worship_dur","sleep_duration","sleep_end","sleep_start",
                       "sleep_heathkit_dur","sms_in_num_ep_1","sms_in_num_ep_2","sms_in_num_ep_3","sms_out_num_ep_1","sms_out_num_ep_2","sms_out_num_ep_3")

##############################
#### All feature datasets ####
##############################

x_features_log <- X_features_hourly_selected %>% 
  mutate(across(all_of(grep(paste(log_vars, collapse = "|"), names(.), value = TRUE)), ~ log(. + 1), .names = "log_{.col}")) %>%
  select(starts_with("log"))

x_features_nonlog <- X_features_hourly_selected %>% 
  mutate(across(all_of(grep(paste(nonlog_var, collapse = "|"), names(.), value = TRUE)), ~ log(. + 1), .names = "nonlog_{.col}")) %>%
  select(starts_with("nonlog"))

X_features_nothourly <- X_features_nothourly %>% mutate(across(everything(), ~ log(. + 1), .names = "{.col}_log")) %>% select(order(names(.)))
x_features_daily <- X_features_nothourly %>% select(all_of(c(log_vars_daily,nonlog_vars_daily)))

x_features_rest <- data$demographics

features_all <- bind_cols(x_features_log,x_features_nonlog,x_features_daily)

df_all <- data$sensing %>% select(uid, day, is_ios) %>%
  left_join(data$demographics, by = "uid") %>%
  bind_cols(features_all) %>%
  filter(gender %in% c("F","M"))



###########################################
#### Cleaned IOS and ANDROID dataframes ###
###########################################

data_ema_selected <- data_ema %>% select(uid, pid, day, pam, phq4_score)

na_cutoff <- 0.25
min_count <- 25

df_android <- df_all %>% filter(is_ios != 1) %>%
  select(where(~ mean(is.na(.)) <= na_cutoff)) %>%
  select(where(~ any(. != 0))) %>%
  drop_na()

df_ios <- df_all %>% filter(is_ios == 1) %>%
  select(where(~ mean(is.na(.)) <= na_cutoff)) %>%
  select(where(~ any(. != 0))) %>%
  drop_na()

df_android <- df_android %>%
  left_join(data_ema_selected, by = c("uid","day")) %>%
  drop_na() %>%
  group_by(pid) %>%
  filter(n() > min_count) %>%
  ungroup()

df_ios <- df_ios %>%
  left_join(data_ema_selected, by = c("uid","day")) %>%
  drop_na() %>%
  group_by(pid) %>%
  filter(n() > min_count) %>%
  ungroup()

########################