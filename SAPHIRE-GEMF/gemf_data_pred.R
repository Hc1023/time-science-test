if(T){
    rm(list = ls())
    library(data.table)
    library(dplyr)
    library(scales)
    realData_all <- read.csv("/share/home/jianglab/huangsisi/usr/time-science/Covid19CasesWH.csv", row.names = 1)  
    observe <- realData_all[25:47,1] # Jan 1 - Jan 23
    
}

#### read file
if(T){
    {
        path <- "/share/home/jianglab/huangsisi/usr/time-science/out.04.27/e6_150output/" 
        fileNames <- dir(path)
        filePath <- sapply(fileNames, function(x){ 
            paste(path,x,sep='/')})
    }
    #Col 1:Time_of_event 2:Total_rate 3:Infected_node
    #    4:Previous_state_of_the_node 5:New_state_of_node
    #    6:"S", 7"E", 8"P1", 9"P2", 10"I1", 11"I2", 12"A1", 13"A2", 14"H", 15"R"
}

## start at Time 0 E=1
gemfdata <- function(data){
    colnames(data) <- c("Time_of_event", "Total_rate", "Infected_node",
                        "Previous_state_of_the_node", "New_state_of_node",
                        # Number of nodes in each state
                        "S", "E", "P1", "P2", "I1", "I2", "A1", "A2", "H", "R",
                        "inducer_nodes")
    data_pred <- data %>% 
        mutate(Time = round(365 * Time_of_event)) %>% 
        group_by(Time) %>%
        summarize(delta_P = cumsum(P1)[length(cumsum(P1))],
                  Infected = 1000000-min(S),
                  P2 = P2[length(P2)],
                  delta_I = sum(I1),
                  delta_A = sum(A1))
    data_pred$P <- c(0, data_pred$P2[-1]) + data_pred$delta_P
    # r*P/Dp
    data_pred <- mutate(data_pred, Onset_expect = round(0.15 * P / 2.3))
        
    return(data_pred)
}

findretro <- function(data_pred, observe){
    mindist <- 1E11
    retro <- 0
    
    for (t in 0:(100-23)){
        
        ## fit with epidemic growth in Wuhan between 
        ## **01 January 2020 and 23 January 2020**
        idx <- data_pred$Time > t & data_pred$Time <= (t+23)
        time <- data_pred$Time[idx]-t
        if(length(time) == 0){next}
        dist_cur <- dist(rbind(observe[time], data_pred$Onset_expect[idx]))/length(time)
        if (dist_cur < mindist){
            retro <- t
            mindist <- dist_cur
        }
    }
    
    return(retro) # retro days
}

data_pred_list <- list()
retro_vec <- {}
count <- 0
for (i in 1:length(filePath)){
    
    data <- data.frame(fread(filePath[i], sep=" "))
    print(fileNames[i]) 
    # successful: ≥1000 people had become infected 
    ## meaningless!!! because of the saturation of nodes: and ≥1 person was still infectious
    ## forward retro!
    succ <- ((1e6 - data[nrow(data), 6]) >= 1000) & (sum(data[nrow(data),8:13]) >= 1)
    if(!succ){next}

    data_pred <- gemfdata(data)
    
    retro <- findretro(data_pred, observe)
    
    count = count+1
    
    retro_vec[count] <- retro
    data_pred_list[[count]] <- data_pred
    
}

# pass: 1 pass rejection sampling
#       0 successfully established dynamics
pass <- rep(0,count)
# data_pred_list
retro_first_cond_vec <- rep(0,count)

for (i in 1:count){
    data_pred <- data_pred_list[[i]]
    cond <- 1
    if(cond == 1){
        # I/A
        idx <- min(which(data_pred[,"delta_I"]==1), 
                   which(data_pred[,"delta_A"]==1))
    }else if(cond == 2){
        # I
        idx <- min(which(data_pred[,"delta_I"]==1))
    } 
    time <- data_pred$Time[idx]
    retro_first_cond <- retro_vec[i] + 1 - time
    retro_first_cond_vec[i] <- retro_first_cond
    
    date_first_cond <- as.Date("2020-01-01") - retro_first_cond
    
    doc_date <- c(as.Date("2019-11-17"), as.Date("2019-12-01"))
    if(date_first_cond < doc_date[2]){
        pass[i] <- 1
    }
}

xl <- max(retro_vec[pass==1]) + 23 
## start at Time 0 E=1
ptime <- 0:xl
anchor <- (xl-22) # anchor : retro + 1 : Jan 1

estN_mat0 <- {}

# anchor - Jan 1
# t = retro + 1 corresponds to Jan 1
# t = time corresponds to Jan 1 - (retro + 1 - time)
# t = 0 : Jan 1 - (retro + 1)
# x_start: anchor - (retro + 1 - time)
estN_mean <- rep(0,(xl+1))
estN_up <- rep(0,(xl+1))
estN_low <- rep(0,(xl+1))

for (xx in 0:xl) {
    
    pred_xx <- sapply(which(pass==1), function(i){
        s <- anchor-(retro_vec[i] + 1)
        ## xx: x coord
        data_pred_list[[i]]$Onset_expect[(data_pred_list[[i]]$Time+s)==xx]
    }) %>% unlist()
    
    if(!length(pred_xx)){next}
    estN_mean[xx+1] <- mean(pred_xx)
    estN_up[xx+1] <- quantile(pred_xx, 0.975)
    estN_low[xx+1] <- quantile(pred_xx, 0.025)
    
}

## plot canvas

# start: time=0
pdf(paste0("/share/home/jianglab/huangsisi/usr/time-science/out.04.27/output/", "gemf", ".pdf"), width = 9, height = 6)
par(mar = c(4, 5, 2.5, 1))
ptime <- 0:xl
{
    plot(ptime, ptime, ylim = c(0, max(observe) * 1.05), 
         xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
    mtext("Onset date (2019-2020)", side = 1, line  = 3, cex = 1.01)
    mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
    
    
    mydate <- as.Date("2020-01-01") - anchor + c(0:xl)
    Sys.setlocale("LC_TIME", "C")
    mydate <- format(mydate, "%b %d")
    axis(1, at = seq(0, xl, 10), labels = mydate[seq(0, xl, 10)+1])
    abline(v = c(1, 10, 23)-1+anchor, lty = 3, lwd = 2, col = "darkgrey")
    text(c(1, 10, 23)-1+anchor, par()$usr[4], 
         labels = mydate[c(1, 10, 23)+anchor-1+1], 
         col = "darkgrey", pos = 3, xpd = T)
    abline(v = anchor, lty = 3, lwd = 2)
    text(anchor, par()$usr[4], labels = mydate[anchor+1], pos = 3, xpd = T)
    text(anchor-4, max(observe), expression(2019))
    text(anchor+4, max(observe), expression(2020))

}

polygon(c(ptime, rev(ptime)), c(estN_up, rev(estN_low)), 
        col = "#F39B7FB2", border = NA)
points(ptime, estN_mean, col = "#BC3C29FF", pch = 17, cex = 0.8)

s <- anchor-retro_first_cond_vec[pass==1]

polygon(c(min(s):max(s), max(s):min(s)), 
        c(rep(200,max(s)-min(s)+1), rev(rep(250,max(s)-min(s)+1))), 
        col = alpha("#F8766D", 1), border = NA)
text((min(s)+max(s))/2, 140, 
     paste("Date of index case:", paste(mydate[min(s)+1],mydate[max(s)+1],sep = " - ")), cex = 1.1)
library(HDInterval)
low <- hdi(s, credMass = 0.95)[1]
up <- hdi(s, credMass = 0.95)[2]

polygon(c(low:up, rev(low:up)), c(rep(300,(up-low+1)), rev(rep(350,(up-low+1)))), 
        col = alpha("#F8766D", 0.5), border = NA) #7876B1FF
text((low+up)/2, 420, 
     paste("95% HPD:", paste(mydate[low+1],mydate[up+1],sep = " - ")), cex = 1.1)

text(median(s), 520, paste("median:", mydate[median(s)+1]), cex = 1.1)
# observed data
points((1:23)+(anchor-1), observe, col = "black", pch = 4, cex = 0.8)

legend("topleft", legend = c("Observed", "SAPHIRE-GEMF"), 
       col = c("black", "#BC3C29FF"), pch = c(4, 17), bty = "n")
dev.off()
save(data_pred_list, s, mydate,
     file = "/share/home/jianglab/huangsisi/usr/time-science/out.04.27/output/gemf.RData")

