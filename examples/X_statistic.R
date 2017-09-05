library(exSTRa)
library(truncnorm)
source("../R/Ttest.R")

# Data munging


str_score <- read_score (
  file = system.file("extdata", "HiSeqXTen_WGS_PCR_2.txt", package = "exSTRa"), 
  database = "./data/repeat_disorders_2017_06_10.xlsx", # for more control, use object from exstra_db_read() instead
  #database = "../disease_repeats/repeat_disorders_2017_06_10.xlsx",
  groups.regex = c(control = "^WGSrpt_0[24]$", case = ""), # here, matches on successive patterns override previous matches # (TODO: maybe should be reversed?)
  filter.low.counts = TRUE
)


# 1) Ricks method to calculate the Statistics 

str_score_p <- Statistics(str_score,20,hthresh=.5)

Classical_stats<- str_score_p$data %>% group_by(sample,locus) %>% summarise(pvlaue = mean(pvalue))
View(Classical_stats)

# 2) Terry's Bayes method

# Uncomment if using Epi25 data 
#Epi_output_files <- Sys.glob("/wehisan/bioinf/Bioinformatics/SNPchipdata/MPS_samples/Epilepsy_WEHI/exSTRa/test_Epi25/output/AUSAUS1*.txt")

alt_test <- Bayes_test(str_score)

Bayes_stats<- alt_test$data %>% group_by(sample,locus) %>% summarise(likelihoodR = mean(LR))

View(Bayes_test)


# 3) Quick plot 








# 4) Do same thing with Epi25 K

Epi_output_files <- Sys.glob("/wehisan/bioinf/Bioinformatics/SNPchipdata/MPS_samples/Epilepsy_WEHI/exSTRa/test_Epi25/output/AUSAUS1*.txt")

# Load all Epi25 scores
load_Epi25_scores <- function(filename){
  read_score(file=filename, database="./data/repeat_disorders_2017_06_10.xlsx",
             groups.regex=c(case="^[AN]", control=""), filter.low.counts=FALSE)
}
str_score<- rbind_score_list(lapply(Epi_output_files, load_Epi25_scores), fill=TRUE)









## -- Standard Bootstrapping way to have P value --##
