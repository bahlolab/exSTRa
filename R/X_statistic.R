library(exSTRa)
library(truncnorm)
source("R/Ttest.R")

# Data munging
#htresh <- .5

str_score <- read_score (
  file = system.file("extdata", "HiSeqXTen_WGS_PCR_2.txt", package = "exSTRa"), 
  database = "./data/repeat_disorders_2017_06_10.xlsx", # for more control, use object from exstra_db_read() instead
  #database = "../disease_repeats/repeat_disorders_2017_06_10.xlsx",
  groups.regex = c(control = "^WGSrpt_0[24]$", case = ""), # here, matches on successive patterns override previous matches # (TODO: maybe should be reversed?)
  filter.low.counts = TRUE
)





# Reads in each quantiles
df_real <- df_real_f(str_score,htresh)
# T test added to data frame
str_score_real <- T_test(df_real,str_score)
                    

#P_value <- data.frame(locus=character(),sample=character(),quantile=character(),rep=numeric(),stringsAsFactors = FALSE)
#I=seq(1,levels(factor(str_score$data$locus)))


I=0
B=20

str_score_p <- str_score
pvalue <- seq(1,100,length.out = dim(str_score_p$data)[1])
str_score_p$data <- cbind(str_score_p$data,pvalue)

df_p <- data.frame(locus=character(),sample=character(),Bnumber =character(), pvalue=numeric(),stringsAsFactors = FALSE)
for (b in seq(1,B)){
  
  # Reads in each quantile (simulate)
  df_simulate <- df_simulate_f(df_real)
  # T test added to data frame (simulate)
  str_score_simulate <- T_test(df_simulate,str_score) 
  
  print(b)
  for (loci in levels(factor(str_score_real$data$locus))){
    for (person in levels(factor(str_score_real$data$sample))){
      
      # Store 1,0 in df_p
      i <- mean((str_score_real$data %>% filter(locus == loci & sample==person) %>% select(xSTAT))[[1]]) > mean((str_score_simulate$data %>% filter(locus == loci & sample==person) %>% select(xSTAT))[[1]])
      df_p[nrow(df_p)+1,]  <- list(loci, person, sprintf("B%03d",b),i) 
    }
  }
}

# Move p value from df_p into str_score_p and sum over sample runs.
for (loci in levels(factor(str_score_real$data$locus))){
  for (person in levels(factor(str_score_real$data$sample))){
    a <- (sum(df_p %>% filter(locus ==loci & sample == person) %>% select(pvalue))+1)/(B+1)
    str_score_p$data[which(str_score_p$data$locus == loci & str_score_p$data$sample == person),]$pvalue <- a
  }
}


df_p %>% group_by(locus, sample) %>% summarise(pvaluesum = (sum(pvalue)+1)/(B+1))
# Quick check for indicator function

str_score_p$data %>% group_by(locus sample) %>% summarise(mean_p = mean(pvalue)) filter(sample == "WGSrpt_10") 









str_score_simulate$data[which(str_score_simulate$data$locus =="HD" & str_score_simulate$data$sample == "WGSrpt_09" ) ,]$xSTAT > str_score_simulate$data[which() ,]$xSTAT
names() <- levels(factor(str_score_real$data$sample))






### Plotting function  ##### 
library(grid)
library(gridBase)


Loci_plot <- levels(factor(str_score_p$data$locus))
#Loci_plot <- Loci_plot[-c(2,4)]


#####
DD1 <- str_score_p$data %>% group_by(locus,sample) %>% summarise(pmean=mean(pvalue))
for (loci_plot in Loci_plot){
  #jpeg(file = sprintf("/home/users/allstaff/degorski.p/exSTRa/plot_output/pvalues/pvalue%s.jpeg",loci_plot))
  #print(loci_plot)
 
  #DD <- unique(   str_score1$data[which(str_score1$data$locus == loci_plot),]$xSTAT)
  #DD <- unique(   str_score_p$data[which(str_score_p$data$locus == loci_plot),]$pvalue)
  
  DD <- (DD1 %>% filter(locus==loci_plot))$pmean 
  names(DD) <- (DD1 %>% filter(locus==loci_plot))$sample
  #names(DD) <- levels(factor(str_score_p$data$sample))
  ## Plot, but suppress the labels
  midpts <- barplot(DD, col=rainbow(20), names.arg="", main = sprintf("P_values for locus %s",loci_plot) , ylab = "T_value")
  
  ## Use grid to add the labels    
  vps <- baseViewports()
  pushViewport(vps$inner, vps$figure, vps$plot)
  
  grid.text(names(DD),
            x = unit(midpts, "native"), y=unit(-1, "lines"),
            just="right", rot=50)
  
  popViewport(3)  
  #dev.off()
}


