# An example of exSTRa usage, for known STR expansion disorder loci


rm(list=ls())
## ---- strexpansion_prepare
library(exSTRa)
setwd("~/github/piotrpython/exSTRa_forked")

str_score <- read_score (
  #file = system.file("extdata", "HiSeqXTen_WGS_PCR_2.txt", package = "exSTRa"), 
  #file = "/wehisan/bioinf/Bioinformatics/SNPchipdata/MPS_samples/Epilepsy_WEHI/exSTRa/test_Epi25/output/AUSAUS1*.txt",
  file = "/wehisan/bioinf/Bioinformatics/SNPchipdata/MPS_samples/Epilepsy_WEHI/exSTRa/test_Epi25/output/AUSAUS10042_exSTRa_scores.txt",
  database = "./data/repeat_disorders_2017_06_10.xlsx", # for more control, use object from exstra_db_read() instead
  #database = "../disease_repeats/repeat_disorders_2017_06_10.xlsx",
  groups.regex = c(control = "^AUSAUS$", case = ""), # here, matches on successive patterns override previous matches # (TODO: maybe should be reversed?)
  filter.low.counts = TRUE
)



# --- # --- #

Epi_output_files <- Sys.glob("/wehisan/bioinf/Bioinformatics/SNPchipdata/MPS_samples/Epilepsy_WEHI/exSTRa/test_Epi25/output/AUSAUS1*.txt")

# Load all Epi25 scores
load_Epi25_scores <- function(filename){
  read_score(file=filename, database="./data/repeat_disorders_2017_06_10.xlsx",
             groups.regex=c(case="^[AN]", control=""), filter.low.counts=FALSE)
}
str_score<- rbind_score_list(lapply(Epi_output_files, load_Epi25_scores), fill=TRUE)



# --- # --- #



test <- str_score$data  %>% select(sample,rep,locus) %>% group_by(sample,locus) %>% summarise(new = max(rep))

loci = "HD"
fudge_number <- 2.0
pi_E <- 0.01
pi_N <- 1- pi_E

normal    <- (test %>% filter(locus == loci))$new
expanded  <- (test %>% filter(locus == loci))$new + fudge_number 


#KDE_n <- approxfun(density(new_normal))
#KDE_e1 <- approxfun(density(new_normal + 1))

KDE_n <- approxfun(density(normal),yright=.00001)
KDE_e1 <- approxfun(density(expanded) ,yleft=0)


# --------

sn <- sample(expanded,1)


probE_d <- KDE_e1(sn)
probN_d <- KDE_n(sn)

post_odds <- (pi_E/pi_N) * (probE_d/probN_d)


# Quick Plot
df <- data.frame(normal,expanded)
df_melt <- melt(df)
g<-ggplot(df_melt)
g+geom_density(aes(x=value,color= variable)) +geom_point(aes(x=sn,y=probN_d)) + geom_point(aes(x=sn,y=probE_d))
