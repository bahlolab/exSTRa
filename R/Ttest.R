library(truncnorm)

remove_outliers <- function(x,alpha=.15, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(alpha,1-alpha), na.rm = na.rm, ...)
  y <- x
  y[x < qnt[1]] <- NA
  y[x > qnt[2]] <- NA
  y <- y[!is.na(y)]
  y
}
stester  <- function(df,locus,person)
{ 
  #Tij <- data.frame(Quantile=character(),rep=numeric(),stringsAsFactors = FALSE)
  Tij <- as.numeric(c())
  for (qnt in levels(factor(df$quantile))){
    
    #Tij[nrow(Tij)+1,] <- list( qnt ,(df[which(df$locus == locus & df$sample ==person & df$quantile == qnt)]  -  mean(df[which(df$locus == locus & df$quantile == qnt),4]))/(Error))
    
    # trimming data
    
    temp_quantile <- df[which(df$locus == locus & df$quantile == qnt),]$rep
    
    q70 <- remove_outliers(temp_quantile,alpha=.15)
    
    # Maths statistics 
    nj=length( q70   )
    sj = sd(q70)
    Error = sj*sqrt(1+1/nj)
    
    Tij[length(Tij)+1] <- (df[which(df$locus == locus & df$sample == person & df$quantile == qnt),]$rep  - mean(q70) )/(Error)    
    
  }
  Tij
}

T_test <- function(df,str_score){
  # Appending column onto data frame
  str_score_dummy <- str_score
  xSTAT <- seq(1,100,length.out = dim(str_score_dummy$data)[1])
  str_score_dummy$data <- cbind(str_score_dummy$data,xSTAT)
  
  # Calculate Ti 
  for (loci in levels(factor(str_score_dummy$data$locus))){
    #print(loci)
    for (person in levels(factor(str_score_dummy$data$sample))){
      Tij <- na.omit(stester(df,loci,person))
      
      D <- length(Tij)
      Ti <- (1/D)*sum(Tij)
   
      str_score_dummy$data[which(str_score_dummy$data$locus == loci & str_score_dummy$data$sample == person),]$xSTAT <- Ti 
      
      
      #str_score_dummy$data[which(str_score_dummy$data$locus == loci & str_score_dummy$data$sample == person),]$rep <- df_simulate[which(df_simulate$locus==loci & df_simulate$sample ==person),]
    }
  }

  str_score_dummy
}

# mreads for each quantile in df object
df_real_f <- function(str_score,htresh){
  # df is used to store quantiles.
  df <- data.frame(locus=character(),sample=character(),quantile=character(),rep=numeric(),stringsAsFactors = FALSE)
  for (loci in levels(str_score$data$locus)){
    
    #q_num <- max(summary(str_score$data$sample[which(str_score$data$locus==loci)]))
    q_num <- 10             
    for ( i in levels(str_score$data$sample)){
      
      i_df <- str_score$data[which(str_score$data$sample==i & str_score$data$locus==loci)]
      q_tile <- quantile(i_df$rep, probs = seq(htresh,1,length.out = q_num))
      # Make dataframe locus sample rep
      mint <- 0
      for (qnt in q_tile){
        mint <- mint + 1
    
        df[nrow(df)+1,] <- list(loci,i,sprintf("Q%03d",mint),qnt)
      }
    }
  }
  df
}

df_simulate_f <- function(df){
  df_simulate <- data.frame(locus=character(),sample=character(),quantile=character(),rep=numeric(),stringsAsFactors = FALSE)
  # Looping over quantile and locus
  for (loci in levels(factor(df$locus))){

    for (quant in levels(factor(df$quantile))){
      
      med_quan  <-  median(df[which(df$locus == loci & df$quantile == quant) ,]$rep)
      d_mad     <-  mad(df[which(df$locus == loci & df$quantile == quant) ,]$rep)/qnorm(0.75)
      m <- length(levels(factor(df$sample)))
      
      # CONDITIONAL SIMULATE dummy simulate after next loop called Prev_simulate
      simulates <- rnorm(m,med_quan,d_mad)
      # simulates <- ordered_statistics()
      count = 1
      for (sims in levels(factor(df$sample))){
        df_simulate[nrow(df_simulate)+1,] <- list(loci,sims,quant,simulates[count])
        count = count + 1 
      }
    }
    for (sim in levels(factor(df_simulate$sample))){
      
      df_simulate[which(df_simulate$locus == loci & df_simulate$sample == sim),]$rep <- sort(df_simulate[which(df_simulate$locus == loci & df_simulate$sample == sim),]$rep)
    }
  }
  df_simulate
}

vector_breeder <- function(df){

   df_breed <- data.frame(locus=character(),sample=character(),quantile=character(),vector=numeric(),stringsAsFactors = FALSE)
   #df_mean <-  data.frame(locus=character(),quantile=character(),mean=numeric(), stdev=numeric(),stringsAsFactors = FALSE)

   qnt <- 0
   for ( loci in levels(factor(df$locus))){
     qnt <- 0

     for ( q in levels(factor(df$quantile))){ #ignore last question the [-1]
       qnt <- qnt +1
       
       for ( person in levels(factor(df$sample))){

         if (q == levels(factor(df$quantile))[1]){
           vec <- (df %>% filter(locus==loci & sample==person))$rep[qnt]
           df_breed[nrow(df_breed)+1,]  <- list(loci,person,q, vec)
         }
         else{
           
           vec <- (df %>% filter(locus==loci & sample==person))$rep[qnt] - (df %>% filter(locus==loci & sample==person))$rep[qnt-1]
           df_breed[nrow(df_breed)+1,]  <- list(loci,person,q, vec)
          }
        }
               #df_mean[nrow(df_mean)+1,] <- list(loci,Q,mean(df_breed %>% filter(locus==loci & qunatile==Q))  ,std(df_breed %>% filter(locus==loci & qunatile==Q))
      }
      print(list(loci,qnt,person,(df %>% filter(locus==loci & sample==person))$rep[qnt]))
    }
    df_mean <- df_breed %>% group_by(locus, quantile) %>% summarise(mean= mean(vector), std= sd(vector))
    list(df_mean,df_breed) 
}




df_simulate_vector <- function(df){
  
  df_simulate <- data.frame(locus=character(),sample=character(),quantile=character(),rep=numeric(),stringsAsFactors = FALSE)
  # Looping over quantile and locus
  for (loci in levels(factor(df$locus))){
    print("top")
    print(loci)
    cnt_sim <- 1
    for (quant in levels(factor(df$quantile))){
      
      med_quan  <-  median(df[which(df$locus == loci & df$quantile == quant) ,]$vector)
      d_mad     <-  mad(df[which(df$locus == loci & df$quantile == quant) ,]$vector)/qnorm(0.75)
      m <- length(levels(factor(df$sample)))
      
      # CONDITIONAL SIMULATE dummy simulate after next loop called Prev_simulate
      simulates <- rtruncnorm(n=m,a=0,b=Inf,mean=med_quan,sd=d_mad)
      # simulates <- ordered_statistics()
      count = 1
      df.sample <- levels(factor(df$sample))
   
      for (sims in df.sample){
        if (quant == levels(factor(df$quantile))[1] ){
          df_simulate[nrow(df_simulate)+1,] <- list(loci,sims,quant,simulates[count])
          
        }
        else{ 
          print(loci)
          print(sims)
          print(cnt_sim)
          print("quantile")
          print(quant)
          print((df %>% filter(locus == loci & sample == sims  & quantile==levels(factor(df$quantile))[cnt_sim]))$vector )
          print(simulates[count])
          print("---------------")
          temp <- (df %>% filter(locus == loci & sample == sims  & quantile==levels(factor(df$quantile))[cnt_sim]))$vector  + simulates[count]
          #temp <- df_simulate[nrow(df_simulate),]$rep + simulates[count] # We all know it is safer if we introduce explicit past. 
          df_simulate[nrow(df_simulate)+1,] <- list(loci,sims,quant,temp)
          
        }
        
        count = count + 1 
      }
      cnt_sim <- cnt_sim + 1
    }
  }
  #print("hello")
  df_simulate
}


Statistics <- function(str_score,B,hthresh=.5){
  
  
  # Reads in each quantiles
  df_real <- df_real_f(str_score,hthresh)
  # T test added to data frame
  str_score_real <- T_test(df_real,str_score)
  
  
  #P_value <- data.frame(locus=character(),sample=character(),quantile=character(),rep=numeric(),stringsAsFactors = FALSE)
  #I=seq(1,levels(factor(str_score$data$locus)))
  
  
  I=0
  #B=20
  
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
  str_score_p
}
  
  
  
  
  
  #df_simulate_vector(df_mean,levels(factor(df_breed$sample)))
  
  
   
 
