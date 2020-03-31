library(fractal)
computQccm <- function(data,seq_P2P,seq_STOCK){
  # browser()
  down_value <- ((seq_P2P%*%seq_P2P)*(seq_STOCK%*%seq_STOCK))^(0.5)
  
  Qcc_m <- NULL
  for (m in 1:(nrow(data)-1)) {
    print(m)
    seq_P2P_t <- seq_P2P
    seq_STOCK_t <- seq_STOCK
    Qcc <- NULL
    for (i in 1:m) {
      seq_P2P_t <- seq_P2P_t[-1]
      seq_STOCK_t <- seq_STOCK_t[-length(seq_STOCK_t)]
      up_value <- seq_P2P_t %*% seq_STOCK_t
      c_i <- up_value/down_value
      Qcc <- c(Qcc,c_i^2/(nrow(data)-i))
      #Qcc <- Qcc + c_i^2/(nrow(data)-i)
    }
    #browser()
    Qcc_m <- c(Qcc_m,sum(Qcc)*(nrow(data))^2 )
  }
  Qcc_m
  
}

generHurstExp <- function(seq_P2P,seq_STOCK,q_segment,s_segment,fitP2p=3,fitStock=3) {
  
  result_square <- matrix(NA,ncol = length(s_segment),nrow = length(q_segment))
  rownames(result_square) <- q_segment
  colnames(result_square) <- s_segment
  
  generFluction <- function(seq_P2P,seq_STOCK,s,q,fitMethod="line",fitP2p,fitStock) {
    #if(s==20){browser()}
    seq_P2P_cent <- scale(seq_P2P, scale = F)
    seq_P2P_new <- cumsum(seq_P2P_cent)
    
    seq_STOCK_cent <- scale(seq_STOCK, scale = F)
    seq_STOCK_new <- cumsum(seq_STOCK_cent)
    
    mode <- length(seq_P2P_new) %% s
    if(mode != 0){
      seq_P2P_new_fw <- seq_P2P_new[1:(length(seq_P2P_new)-mode)]#向前
      seq_P2P_new_bk <- seq_P2P_new[(mode+1):length(seq_P2P_new)]#向后
      seq_P2P_new_tem <- c(seq_P2P_new_fw,seq_P2P_new_bk)
      
      seq_STOCK_new_fw <- seq_STOCK_new[1:(length(seq_STOCK_new)-mode)]#向前
      seq_STOCK_new_bk <- seq_STOCK_new[(mode+1):length(seq_STOCK_new)]#向后
      seq_STOCK_new_tem <- c(seq_STOCK_new_fw,seq_STOCK_new_bk)
      
      seq_P2P_STOCK <- rbind(seq_P2P_new_tem,seq_STOCK_new_tem)
      
    }else{
      seq_P2P_new_tem <- seq_P2P_new
      seq_STOCK_new_tem <- seq_STOCK_new
      seq_P2P_STOCK <- rbind(seq_P2P_new_tem,seq_STOCK_new_tem)
    }
    
    segment_n <- ncol(seq_P2P_STOCK)%/%s
    
    x <- 1:s
    #seq_P2P_STOCK[,1:s]
    
    star_flag <- 1
    
    
    if(q != 0){
      
      ever_segment_F <- sapply(1:segment_n, function(i){
        # browser()
        end_flag <- s*i
        if(fitMethod=="line"){
          
          P2P_fitResult <- lm(seq_P2P_STOCK[1,star_flag:end_flag]~x)
          STOCK_fitResult <- lm(seq_P2P_STOCK[2,star_flag:end_flag]~x)
        }else{
          #set.seed(1000) 
          P2P_fitResult <- lm(seq_P2P_STOCK[1,star_flag:end_flag]~poly(x,fitP2p))#取3最好
          STOCK_fitResult <- lm(seq_P2P_STOCK[2,star_flag:end_flag]~poly(x,fitStock))
        }  
        
        #P2P_fitResult <- lm(seq_P2P_STOCK[1,star_flag:end_flag]~x)
        P2P_residuals <- abs(P2P_fitResult$residuals)
        #STOCK_fitResult <- lm(seq_P2P_STOCK[2,star_flag:end_flag]~x)
        star_flag <<- end_flag +1
        STOCK_residuals <- abs(STOCK_fitResult$residuals)
        mean_residuals <- mean(P2P_residuals*STOCK_residuals)
        
        mean_residuals^(q/2)
        
      })  
      
      F_qs_value <- (mean(ever_segment_F))^(1/q)  
    }else{
      
      ever_segment_F <- sapply(1:segment_n, function(i){
        # browser()
        end_flag <- s*i
        if(fitMethod=="line"){
          
          P2P_fitResult <- lm(seq_P2P_STOCK[1,star_flag:end_flag]~x)
          STOCK_fitResult <- lm(seq_P2P_STOCK[2,star_flag:end_flag]~x)
        }else{
          #set.seed(1000)
          P2P_fitResult <- lm(seq_P2P_STOCK[1,star_flag:end_flag]~poly(x,2))#取3最好
          STOCK_fitResult <- lm(seq_P2P_STOCK[2,star_flag:end_flag]~poly(x,3))
        }  
        
        #P2P_fitResult <- lm(seq_P2P_STOCK[1,star_flag:end_flag]~x)
        P2P_residuals <- abs(P2P_fitResult$residuals)
        #STOCK_fitResult <- lm(seq_P2P_STOCK[2,star_flag:end_flag]~x)
        star_flag <<- end_flag +1
        STOCK_residuals <- abs(STOCK_fitResult$residuals)
        mean_residuals <- mean(P2P_residuals*STOCK_residuals)
        
        log(mean_residuals)
        
      })   
      F_qs_value <- exp(0.5 * mean(ever_segment_F)) 
    }
    
    F_qs_value
  }
  
  #square
  for (q in 1:length(q_segment)) {
    print(paste0("q--",q))
    #if(q==4){browser()}
    for (s in 1:length(s_segment)) {
      # print(paste0("s--",s))
      #browser()
      result_square[q,s] <-  generFluction(seq_P2P,seq_STOCK,s_segment[s],q_segment[q],"suqure",fitP2p,fitStock)
    }
    
  }
  
  #browser()
  lnValue <- log(result_square)
  result_final_square <- NULL
  for (i in 1:nrow(lnValue)) {
    #lmResult <- lm(lnValue[i,4:ncol(result_square)]~log(s_segment[4:ncol(result_square)]))#取s=40最好
    lmResult <- lm(lnValue[i,]~log(s_segment))
    temp <- c(lmResult$coefficients[1],lmResult$coefficients[2])
    result_final_square <- rbind(result_final_square,temp)
    
  }
  rownames(result_final_square) <- q_segment
  
  result_final_square <- cbind(q_segment,result_final_square)
  colnames(result_final_square) <- c("q","Intercept","Hxy")
  result_final_square
}

#------------------------------------

data <- read.csv("data.csv",header = T)
seq_P2P <- data$P2P
seq_STOCK <- data$STOCK
q_segment <- seq(-10,10,1)
s_segment <- seq(40,200,10)

#-----------------------------------
Qccm <- computQccm(data,seq_P2P,seq_STOCK)
#------------------------------------
result_daoshu <- generHurstExp(seq_P2P,seq_STOCK,q_segment,s_segment)

raw <- generHurstExp(seq_P2P,seq_STOCK,q_segment,s_segment)
one <- generHurstExp(seq_P2P,seq_STOCK,q_segment-0.0001,s_segment)
two <- generHurstExp(seq_P2P,seq_STOCK,q_segment+0.0001,s_segment)
#---------------------------------
#---------------------------------
alfa_value <- q_segment*(two[,3]- one[,3])/0.0002 + raw[,3]
f_alfa_value <- q_segment*(alfa_value - raw[,3]) + 1
plot(alfa_value,f_alfa_value)
#------------------------------------

raw_seq_P2P <- generHurstExp(seq_P2P,seq_P2P,q_segment,s_segment)
raw_seq_STOCK <- generHurstExp(seq_STOCK,seq_STOCK,q_segment,s_segment)
#---------------------------------

tau_q_raw <- raw[,1]*raw[,3] -1
tau_seq_P2P <- raw_seq_P2P[,1]*raw_seq_P2P[,3] -1
tau_seq_STOCK <-raw_seq_STOCK[,1]*raw_seq_STOCK[,3] -1
#save(Qccm,Hxy,alfa_value,f_alfa_value,tau_q,file = "allResult.RData")
