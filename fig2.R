if(F){
    rm(list = ls())
    library(data.table)
    library(dplyr)
    library(scales)
    library(HDInterval)
    realData_all <- read.csv("D:/OneDrive - zju.edu.cn/lab/git/SAPHIRE/data/Covid19CasesWH.csv", row.names = 1)  
    observe <- realData_all[25:47,1] # Jan 1 - Jan 23
    load("D:/OneDrive - zju.edu.cn/lab/git/time-science-test/SAPHIRE-GEMF/gemf.RData")
    load("D:/OneDrive - zju.edu.cn/lab/git/time-science-test/SAPHIRE/retro_list.RData")
} 

#
outname="FIG2_back_0517"
pdf(paste0("D:/OneDrive - zju.edu.cn/lab/git/time-science-test/", 
           outname, ".pdf"), width = 15, height = 10)
par(mar = c(4, 5, 2.5, 1))
layout(matrix(c(1:4), byrow = T, nrow = 2))

#######################   Panel A  ############################

cond <- c(1,2)
te <- c(1,2)
# doc_date <- c(as.Date("2019-11-17"), as.Date("2019-12-01"))
outname <- c("1117", "1201", "1117I", "1201I")

cond <- cond[1]
te <- te[1]
outname <- outname[(cond-1)*2+te]

findretro <- function(data_pred, observe){
    mindist <- 1E11
    retro <- 0
    
    for (t in 1:(100-23)){
        
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


count <- length(data_pred_list)
retro_vec <- {}
for (i in 1:count){
    
    data_pred <- data_pred_list[[i]]
    retro <- findretro(data_pred, observe)
    
    retro_vec[i] <- retro
    
}

# pass: 1 pass rejection sampling
#       0 successfully established dynamics
pass <- rep(0,count)
# data_pred_list
retro_first_cond_vec <- rep(0,count)

for (i in 1:count){
    data_pred <- data_pred_list[[i]]

    if(cond == 1){
        # I/A
        idx <- min(which(data_pred[,"delta_I"]>=1), 
                   which(data_pred[,"delta_A"]>=1))
    }else if(cond == 2){
        # I
        idx <- min(which(data_pred[,"delta_I"]>=1))
    } 
    time <- data_pred$Time[idx]
    retro_first_cond <- retro_vec[i] + 1 - time
    
    
    date_first_cond <- as.Date("2020-01-01") - retro_first_cond
    
    
    doc_date <- c(as.Date("2019-11-17"), as.Date("2019-12-01"))
    if(date_first_cond < doc_date[te]){
        pass[i] <- 1
    }
    
    if(cond == 2){
        idx <- min(which(data_pred[,"delta_I"]>=1), 
                   which(data_pred[,"delta_A"]>=1))
        time <- data_pred$Time[idx]
        retro_first_cond <- retro_vec[i] + 1 - time
    }
    
    retro_first_cond_vec[i] <- retro_first_cond
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

# pdf(paste0("D:/OneDrive - zju.edu.cn/lab/git/time-science-test/0427/5000/", outname, ".pdf"), width = 9, height = 6)
# par(mar = c(4, 5, 2.5, 1))
ptime <- 0:xl
{
    plot(ptime, ptime, ylim = c(0, max(observe) * 1.05), 
         xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
    grid(10, 10, lty = 1, col = alpha("lightgray",0.4))
    mtext("Onset date (2019-2020)", side = 1, line  = 3, cex = 1.01)
    mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
    
    
    mydate <- as.Date("2020-01-01") - anchor + c(0:xl)
    Sys.setlocale("LC_TIME", "C")
    mydate <- format(mydate, "%b %d")
    axis(1, at = seq(0, xl, 10), labels = mydate[seq(0, xl, 10)+1])
    abline(v = c(10, 23)-1+anchor, lty = 3, lwd = 2, col = "darkgrey")
    text(c(10, 23)-1+anchor, par()$usr[4], 
         labels = mydate[c(10, 23)+anchor-1+1], 
         col = "darkgrey", pos = 3, xpd = T)
    abline(v = anchor, lty = 3, lwd = 2)
    REJ1 <- anchor-as.numeric(as.Date("2020-01-01") - as.Date("2019-11-17"))
    REJ2 <- anchor-as.numeric(as.Date("2020-01-01") - as.Date("2019-12-01"))
    text(c(REJ1, REJ2), par()$usr[4], 
         labels = mydate[c(REJ1+1, REJ2+1)], 
         col = "darkgrey", pos = 3, xpd = T)
    text(anchor, par()$usr[4], labels = mydate[anchor+1], pos = 3, xpd = T)
    text(anchor-4, max(observe), expression(2019))
    text(anchor+4, max(observe), expression(2020))
    
}

polygon(c(ptime, rev(ptime)), c(estN_up, rev(estN_low)), 
        col = "#F39B7FB2", border = NA)
points(ptime, estN_mean, col = "#BC3C29FF", pch = 16, cex = 0.8)

s <- anchor-retro_first_cond_vec[pass==1]

polygon(c(min(s):max(s), max(s):min(s)),
        c(rep(60,max(s)-min(s)+1), rev(rep(80,max(s)-min(s)+1))),
        col = alpha("#BC3C29FF", 1), border = NA)
text((min(s)+max(s))/2, 140,
     paste("Date of index case:", paste(mydate[min(s)+1],mydate[max(s)+1],sep = " - ")), cex = 1.1)

# text((min(s)+max(s))/2, 140,
#      "Date of the index case", cex = 1.1)

low <- hdi(s, credMass = 0.99)[1]
up <- hdi(s, credMass = 0.99)[2]
polygon(c(low:up, rev(low:up)), c(rep(200,(up-low+1)), rev(rep(250,(up-low+1)))), 
        col = alpha("#F8766D", 1), border = NA) #7876B1FF
text((low+up)/2, 300, 
     paste("99% HPD:", paste(mydate[low+1],mydate[up+1],sep = " - ")), cex = 1.1)


low <- hdi(s, credMass = 0.95)[1]
up <- hdi(s, credMass = 0.95)[2]

polygon(c(low:up, rev(low:up)), c(rep(350,(up-low+1)), rev(rep(400,(up-low+1)))), 
        col = alpha("#F8766D", 0.5), border = NA) #7876B1FF
text((low+up)/2, 450, 
     paste("95% HPD:", paste(mydate[low+1],mydate[up+1],sep = " - ")), cex = 1.1)

text(median(s), 540, paste("Median:", mydate[median(s)+1]), cex = 1.1)
# observed data
points((1:23)+(anchor-1), observe, col = "black", pch = 4, cex = 0.8)

legend("topleft", legend = c("Observed", "SAPHIRE-GEMF"), 
       col = c("black", "#BC3C29FF"), pch = c(4, 16), bty = "n")
# dev.off()
text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.1, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "A", xpd = T, cex = 2)

S.A <- s
REJ1.A <- REJ1
REJ2.A <- REJ2
#######################   Panel B  ############################

rejection_sampling <- function(pred_mat, doc_date, cond){
    date_first_cond <- min(which(pred_mat[,"I"]>=1), which(pred_mat[,"A"]>=1))
    if(cond == 1){
        # I/A
        date_first_cond2 <- min(which(pred_mat[,"I"]>=1), which(pred_mat[,"A"]>=1))
    }else if(cond == 2){
        # I
        date_first_cond2 <- min(which(pred_mat[,"I"]>=1))
    } 
    retro <- nrow(pred_mat)-68
    retro_index_case <- retro - (date_first_cond-1)
    retro_index_case2 <- retro - (date_first_cond2-1)
    date_first_cond2 <- as.Date("2020-01-01") - retro_index_case2
    if(date_first_cond2 < doc_date){
        return(c(1, retro_index_case)) 
    }
    return(c(0, retro_index_case))
}

nn <- 150

cond <- c(1,2)
te <- c(1,2)
# doc_date <- c(as.Date("2019-11-17"), as.Date("2019-12-01"))
outname <- c("1117", "1201", "1117I", "1201I")
cond <- cond[1]
te <- te[1]
outname <- outname[(cond-1)*2+te]

estN_mat0 <- {}
retro_vec <- {}
rvec <- {}
for (i in 1:1500) {
    pred_retro <- pred_retro_list[[i]]
    if (pred_retro$retro == 0){next}
    rvec <- c(rvec, pred_retro$retro)
    doc_date <- c(as.Date("2019-11-17"), as.Date("2019-12-01"))
    retro_index_case <- rejection_sampling(pred_retro$pred_mat, 
                                           doc_date[te], cond = cond)
    
    date_index_case <- as.Date("2020-01-01") - retro_index_case[2]
    date_index_case <- format(date_index_case, "%b %d")
    # print(paste("retro:", retro_index_case[2]))
    # print(paste("date of index case:", date_index_case))
    
    pred <- c(rep(0,nn-pred_retro$retro), pred_retro$pred)
    estN_mat0 <- cbind(estN_mat0, pred)
    retro_vec <- rbind(retro_vec, retro_index_case)
}
retro_vec1 <- retro_vec
retro_vec <- retro_vec[retro_vec[,1]==1,2]

if(nn == max(rvec)){
    print("warning")
    estN_mat <- estN_mat0
}else{
    estN_mat <- estN_mat0[-c(1:(nn-max(rvec))),]
}

estN_mean <- round(apply(estN_mat[-c((130-44):130),], 1, mean), 0)
estN_up <- round(apply(estN_mat[-c((130-44):130),], 1, function(x) quantile(x, 0.975)), 0)
estN_low <- round(apply(estN_mat[-c((130-44):130),], 1, function(x) quantile(x, 0.025)), 0)

# pdf(paste0("D:/OneDrive - zju.edu.cn/lab/git/time-science-test/0427/", outname, ".pdf"), width = 12, height = 6)
# par(mar = c(4, 5, 2.5, 1))
Sys.setlocale("LC_TIME", "C")
{
    {
        ptime <- 1:length(estN_mean)
        plot(ptime, estN_mean, ylim = c(0, max(estN_up, observe) * 1.05), 
             xlab = "", ylab = "", type = "p", col = "white", 
             pch = 16, xaxt="n", cex = 0.5)
        grid(10, 10, lty = 1, col = alpha("lightgray",0.4))
        mtext("Onset date (2019-2020)", side = 1, line  = 3, cex = 1.01)
        mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
        mydate <- as.Date("2020-01-01") - max(rvec) + c(0:(68+max(rvec)-1))
        mydate <- format(mydate, "%b %d")
        axis(1, at = seq(1, 68+max(rvec), 11), 
             labels = mydate[seq(1, 68+max(rvec), 11)])
        #
        abline(v = c(10, 23)+max(rvec), lty = 3, lwd = 2, col = "darkgrey")
        text(c(10, 23)+max(rvec), par()$usr[4], 
             labels = mydate[c(10, 23)+max(rvec)], 
             col = "darkgrey", pos = 3, xpd = T)
        abline(v = 1+max(rvec), lty = 3, lwd = 2)
        
        REJ1 <- 1+max(rvec)-as.numeric(as.Date("2020-01-01") - as.Date("2019-11-17"))
        REJ2 <- 1+max(rvec)-as.numeric(as.Date("2020-01-01") - as.Date("2019-12-01"))
        text(c(REJ1, REJ2), par()$usr[4], 
             labels = mydate[c(REJ1, REJ2)], 
             col = "darkgrey", pos = 3, xpd = T)
        
        text(1+max(rvec), par()$usr[4], labels = mydate[1+max(rvec)], pos = 3, xpd = T)
        text(max(rvec)-3.5, max(estN_up, observe), expression(2019))
        text(max(rvec)+5.5, max(estN_up, observe), expression(2020))
    }
    #
    polygon(c(ptime, rev(ptime)), c(estN_up, rev(estN_low)), 
            col = "#4DBBD5B2", border = NA)
    
    t <- 1+max(retro_vec)-min(retro_vec)
    dif <- max(rvec)-max(retro_vec)
    k <- max(estN_up)/max(observe)
    polygon(c(1:t + dif, rev(1:t + dif)), c(rep(60*k,t), rev(rep(80*k,t))),
            col = alpha("#0072B5FF", 1), border = NA)
    text(t/2+8+dif, 140*k,
         paste("Date of index case:", paste(mydate[1+dif],mydate[t+dif],sep = " - ")), cex = 1.1)
    
    # text(t/2+5+dif, 140, "Date of the index case", cex = 1.1)
    
    up <- max(retro_vec)-hdi(retro_vec, credMass = 0.99)[1]+1
    low <- max(retro_vec)-hdi(retro_vec, credMass = 0.99)[2]+1
    
    polygon(c(low:up+dif, rev(low:up+dif)), c(rep(200*k,(up-low+1)), rev(rep(250*k,(up-low+1)))), 
            col = alpha("#7876B1CC", 1), border = NA) #7876B1FF
    text((low+up)/2+5+dif, 300*k, 
         paste("99% HPD:", paste(mydate[low+dif],mydate[up+dif],sep = " - ")), cex = 1.1)
    
    up <- max(retro_vec)-hdi(retro_vec, credMass = 0.95)[1]+1
    low <- max(retro_vec)-hdi(retro_vec, credMass = 0.95)[2]+1
    
    polygon(c(low:up+dif, rev(low:up+dif)), c(rep(350*k,(up-low+1)), rev(rep(400*k,(up-low+1)))), 
            col = alpha("#7876B1CC", 0.5), border = NA) #7876B1FF
    text((low+up)/2+5+dif, 450*k, 
         paste("95% HPD:", paste(mydate[low+dif],mydate[up+dif],sep = " - ")), cex = 1.1)
    ME <- 1+max(retro_vec)-median(retro_vec)
    text(ME+dif, 540*k, paste("Median:", mydate[ME+dif]), cex = 1.1)
    
    points(ptime, estN_mean, col = "#0072B5FF", pch = 17, cex = 0.8)
    points((1:23)+max(rvec), observe, col = "black", pch = 4, cex = 0.8)
    #
    legend("topleft", legend = c("Observed", "SAPHIRE"), col = c("black",  "#0072B5FF"), pch = c(4, 17), bty = "n")
    
}
text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.1, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "B", xpd = T, cex = 2)
S.B <- 1+max(retro_vec)-retro_vec
REJ1.B <- REJ1-dif
REJ2.B <- REJ1-dif

#######################   Panel C  ############################

count <- length(data_pred_list)
retro_vec <- {}
for (i in 1:count){
    
    data_pred <- data_pred_list[[i]]
    retro <- findretro(data_pred, observe)
    
    retro_vec[i] <- retro
    
}

# no rejection sampling
# data_pred_list
retro_first_cond_vec <- rep(0,count)

for (i in 1:count){
    data_pred <- data_pred_list[[i]]
    
    # I/A
    idx <- min(which(data_pred[,"delta_I"]>=1), 
               which(data_pred[,"delta_A"]>=1))
    
    time <- data_pred$Time[idx]
    retro_first_cond <- max(1,retro_vec[i] + 1 - time)
    retro_first_cond_vec[i] <- retro_first_cond
}

xl <- max(retro_vec) + 23 
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
pass <- rep(1,count)
for (xx in 0:xl) {
    
    pred_xx <- sapply(1:1513, function(i){
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

# pdf(paste0("D:/OneDrive - zju.edu.cn/lab/git/time-science-test/0427/5000/", outname, ".pdf"), width = 9, height = 6)
# par(mar = c(4, 5, 2.5, 1))
ptime <- 0:xl
{
    plot(ptime, ptime, ylim = c(0, max(observe) * 1.05), 
         xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
    grid(10,10, lty = 1, col = alpha("lightgray",0.4))
    mtext("Onset date (2019-2020)", side = 1, line  = 3, cex = 1.01)
    mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
    
    
    mydate <- as.Date("2020-01-01") - anchor + c(0:xl)
    Sys.setlocale("LC_TIME", "C")
    mydate <- format(mydate, "%b %d")
    axis(1, at = seq(0, xl, 10), labels = mydate[seq(0, xl, 10)+1])
    abline(v = c(10, 23)-1+anchor, lty = 3, lwd = 2, col = "darkgrey")
    text(c(10, 23)-1+anchor, par()$usr[4], 
         labels = mydate[c(10, 23)+anchor-1+1], 
         col = "darkgrey", pos = 3, xpd = T)
    abline(v = anchor, lty = 3, lwd = 2)
    
    REJ1 <- anchor-as.numeric(as.Date("2020-01-01") - as.Date("2019-11-17"))
    REJ2 <- anchor-as.numeric(as.Date("2020-01-01") - as.Date("2019-12-01"))
    text(c(REJ1, REJ2), par()$usr[4], 
         labels = mydate[c(REJ1+1, REJ2+1)], 
         col = "darkgrey", pos = 3, xpd = T)
    text(anchor, par()$usr[4], labels = mydate[anchor+1], pos = 3, xpd = T)
    text(anchor-4, max(observe), expression(2019))
    text(anchor+4, max(observe), expression(2020))
    
}


polygon(c(ptime, rev(ptime)), c(estN_up, rev(estN_low)), 
        col = "#F39B7FB2", border = NA)
points(ptime, estN_mean, col = "#BC3C29FF", pch = 16, cex = 0.8)

s <- anchor-retro_first_cond_vec[pass==1]

polygon(c(min(s):max(s), max(s):min(s)),
        c(rep(60,max(s)-min(s)+1), rev(rep(80,max(s)-min(s)+1))),
        col = alpha("#BC3C29FF", 1), border = NA)
text((min(s)+max(s))/2, 140,
     paste("Date of index case:", paste(mydate[min(s)+1],mydate[max(s)+1],sep = " - ")), cex = 1.1)

# text((min(s)+max(s))/2, 140,
#      "Date of the index case", cex = 1.1)


low <- hdi(s, credMass = 0.99)[1]
up <- hdi(s, credMass = 0.99)[2]
polygon(c(low:up, rev(low:up)), c(rep(200,(up-low+1)), rev(rep(250,(up-low+1)))), 
        col = alpha("#F8766D", 1), border = NA) #7876B1FF
text((low+up)/2, 300, 
     paste("99% HPD:", paste(mydate[low+1],mydate[up+1],sep = " - ")), cex = 1.1)

low <- hdi(s, credMass = 0.95)[1]
up <- hdi(s, credMass = 0.95)[2]

polygon(c(low:up, rev(low:up)), c(rep(350,(up-low+1)), rev(rep(400,(up-low+1)))), 
        col = alpha("#F8766D", 0.5), border = NA) #7876B1FF
text((low+up)/2, 450, 
     paste("95% HPD:", paste(mydate[low+1],mydate[up+1],sep = " - ")), cex = 1.1)

text(median(s), 540, paste("Median:", mydate[median(s)+1]), cex = 1.1)
# observed data
points((1:23)+(anchor-1), observe, col = "black", pch = 4, cex = 0.8)

legend("topleft", legend = c("Observed", "SAPHIRE-GEMF"), 
       col = c("black", "#BC3C29FF"), pch = c(4, 16), bty = "n")

text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.1, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "C", xpd = T, cex = 2)

S.C <- s
REJ1.C <- REJ1
REJ2.C <- REJ2

#######################   Panel D  ############################
rejection_sampling <- function(pred_mat, doc_date, cond){
    date_first_cond <- min(which(pred_mat[,"I"]>=1), which(pred_mat[,"A"]>=1))
    if(cond == 1){
        # I/A
        date_first_cond2 <- min(which(pred_mat[,"I"]>=1), which(pred_mat[,"A"]>=1))
    }else if(cond == 2){
        # I
        date_first_cond2 <- min(which(pred_mat[,"I"]>=1))
    } 
    retro <- nrow(pred_mat)-68
    retro_index_case <- retro - (date_first_cond-1)
    retro_index_case2 <- retro - (date_first_cond2-1)
    date_first_cond2 <- as.Date("2020-01-01") - retro_index_case2
    if(date_first_cond2 < doc_date){
        return(c(1, retro_index_case))
    }
    return(c(1, retro_index_case))
}

nn <- 150
cond <- c(1,2)
te <- c(1,2)
# doc_date <- c(as.Date("2019-11-17"), as.Date("2019-12-01"))
outname <- c("1117", "1201", "1117I", "1201I")
cond <- cond[1]
te <- te[1]
outname <- outname[(cond-1)*2+te]

estN_mat0 <- {}
retro_vec <- {}
rvec <- {}
for (i in 1:1500) {
    pred_retro <- pred_retro_list[[i]]
    if (pred_retro$retro == 0){next}
    rvec <- c(rvec, pred_retro$retro)
    doc_date <- c(as.Date("2019-11-17"), as.Date("2019-12-01"))
    retro_index_case <- rejection_sampling(pred_retro$pred_mat, 
                                           doc_date[te], cond = cond)
    
    date_index_case <- as.Date("2020-01-01") - retro_index_case[2]
    date_index_case <- format(date_index_case, "%b %d")
    # print(paste("retro:", retro_index_case[2]))
    # print(paste("date of index case:", date_index_case))
    
    pred <- c(rep(0,nn-pred_retro$retro), pred_retro$pred)
    estN_mat0 <- cbind(estN_mat0, pred)
    retro_vec <- rbind(retro_vec, retro_index_case)
}
retro_vec1 <- retro_vec
retro_vec <- retro_vec[retro_vec[,1]==1,2]

if(nn == max(rvec)){
    print("warning")
    estN_mat <- estN_mat0
}else{
    estN_mat <- estN_mat0[-c(1:(nn-max(rvec))),]
}

estN_mean <- round(apply(estN_mat[-c((130-44):130),], 1, mean), 0)
estN_up <- round(apply(estN_mat[-c((130-44):130),], 1, function(x) quantile(x, 0.975)), 0)
estN_low <- round(apply(estN_mat[-c((130-44):130),], 1, function(x) quantile(x, 0.025)), 0)

# pdf(paste0("D:/OneDrive - zju.edu.cn/lab/git/time-science-test/0427/", outname, ".pdf"), width = 12, height = 6)
# par(mar = c(4, 5, 2.5, 1))
Sys.setlocale("LC_TIME", "C")
{
    {
        ptime <- 1:length(estN_mean)
        plot(ptime, estN_mean, ylim = c(0, max(estN_up, observe) * 1.05), 
             xlab = "", ylab = "", type = "p", col = "white", 
             pch = 16, xaxt="n", cex = 0.5)
        grid(10, 10, lty = 1, col = alpha("lightgray",0.4))
        mtext("Onset date (2019-2020)", side = 1, line  = 3, cex = 1.01)
        mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
        mydate <- as.Date("2020-01-01") - max(rvec) + c(0:(68+max(rvec)-1))
        mydate <- format(mydate, "%b %d")
        axis(1, at = seq(1, 68+max(rvec), 11), 
             labels = mydate[seq(1, 68+max(rvec), 11)])
        #
        abline(v = c(10, 23)+max(rvec), lty = 3, lwd = 2, col = "darkgrey")
        text(c(10, 23)+max(rvec), par()$usr[4], 
             labels = mydate[c(10, 23)+max(rvec)], 
             col = "darkgrey", pos = 3, xpd = T)
        abline(v = 1+max(rvec), lty = 3, lwd = 2)
        
        REJ1 <- 1+max(rvec)-as.numeric(as.Date("2020-01-01") - as.Date("2019-11-17"))
        REJ2 <- 1+max(rvec)-as.numeric(as.Date("2020-01-01") - as.Date("2019-12-01"))
        text(c(REJ1, REJ2), par()$usr[4], 
             labels = mydate[c(REJ1, REJ2)], 
             col = "darkgrey", pos = 3, xpd = T)
        
        text(1+max(rvec), par()$usr[4], labels = mydate[1+max(rvec)], pos = 3, xpd = T)
        text(max(rvec)-3.5, max(estN_up, observe), expression(2019))
        text(max(rvec)+5.5, max(estN_up, observe), expression(2020))
    }
    #
    polygon(c(ptime, rev(ptime)), c(estN_up, rev(estN_low)), 
            col = "#4DBBD5B2", border = NA)
    
    t <- 1+max(retro_vec)-min(retro_vec)
    dif <- max(rvec)-max(retro_vec)
    k <- max(estN_up)/max(observe)
    polygon(c(1:t + dif, rev(1:t + dif)), c(rep(60*k,t), rev(rep(80*k,t))),
            col = alpha("#0072B5FF", 1), border = NA)
    text(t/2+8+dif, 140*k,
         paste("Date of index case:", paste(mydate[1+dif],mydate[t+dif],sep = " - ")), cex = 1.1)

    # text(t/2+8+dif, 140, "Date of the index case", cex = 1.1)
    
    up <- max(retro_vec)-hdi(retro_vec, credMass = 0.99)[1]+1
    low <- max(retro_vec)-hdi(retro_vec, credMass = 0.99)[2]+1
    
    polygon(c(low:up+dif, rev(low:up+dif)), c(rep(200*k,(up-low+1)), rev(rep(250*k,(up-low+1)))), 
            col = alpha("#7876B1CC", 1), border = NA) #7876B1FF
    text((low+up)/2+5+dif, 300*k, 
         paste("99% HPD:", paste(mydate[low+dif],mydate[up+dif],sep = " - ")), cex = 1.1)
    
    up <- max(retro_vec)-hdi(retro_vec, credMass = 0.95)[1]+1
    low <- max(retro_vec)-hdi(retro_vec, credMass = 0.95)[2]+1
    
    polygon(c(low:up+dif, rev(low:up+dif)), c(rep(350*k,(up-low+1)), rev(rep(400*k,(up-low+1)))), 
            col = alpha("#7876B1CC", 0.5), border = NA) #7876B1FF
    text((low+up)/2+5+dif, 450*k, 
         paste("95% HPD:", paste(mydate[low+dif],mydate[up+dif],sep = " - ")), cex = 1.1)
    ME <- 1+max(retro_vec)-median(retro_vec)
    text(ME+dif, 540*k, paste("Median:", mydate[ME+dif]), cex = 1.1)

    points(ptime, estN_mean, col = "#0072B5FF", pch = 17, cex = 0.8)
    points((1:23)+max(rvec), observe, col = "black", pch = 4, cex = 0.8)
    #
    legend("topleft", legend = c("Observed", "SAPHIRE"), col = c("black",  "#0072B5FF"), 
           pch = c(4, 17), bty = "n")
    
}
text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.1, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "D", xpd = T, cex = 2)

S.D <- 1+max(retro_vec)-retro_vec
REJ1.D <- REJ1-dif
REJ2.D <- REJ2-dif

dev.off()


#
# layout(matrix(c(1:4), byrow = T, nrow = 2))

#######################   Panel A  ############################
outname="FIG2_d1"
pdf(paste0("D:/OneDrive - zju.edu.cn/lab/git/time-science-test/0427/0514/", 
           outname, ".pdf"), width = 6, height = 4)

ggdensity(data.frame(x=S.A, 
                     col=rep("SAPHIRE-GEMF-REJ",length(S.A))), 
          x = "x", add = "median", rug = TRUE,
          color = "black", fill = "col", alpha = 0.8,
          palette = c("#E18727FF"))
dev.off()
# text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.1, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "A", xpd = T, cex = 2)

#######################   Panel B  ############################
outname="FIG2_d2"
pdf(paste0("D:/OneDrive - zju.edu.cn/lab/git/time-science-test/0427/0514/", 
           outname, ".pdf"), width = 6, height = 4)
ggdensity(data.frame(x=S.B, col = rep("SAPHIRE-GEMF",length(S.B))), x = "x",
                  add = "median", rug = TRUE,
                  color = "black", fill = "col",
                  palette = c("#EFC000FF")) +
    geom_vline(xintercept = c(REJ1.B, REJ2.B), linetype="dotted",
               color = "darkgrey", size = 1.5)
# text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.1, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "D", xpd = T, cex = 2)
dev.off()
#######################   Panel C  ############################
outname="FIG3_d3"
pdf(paste0("D:/OneDrive - zju.edu.cn/lab/git/time-science-test/0427/0514/", 
           outname, ".pdf"), width = 6, height = 4)
ggdensity(data.frame(x=S.C, 
                     col=rep("SAPHIRE-REJ",length(S.C))), 
          x = "x", add = "median", rug = TRUE,
          color = "black", fill = "col", alpha = 0.8,
          palette = c("#E18727FF"))
# text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.1, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "D", xpd = T, cex = 2)
dev.off()
#######################   Panel D  ############################
outname="FIG3_d4"
pdf(paste0("D:/OneDrive - zju.edu.cn/lab/git/time-science-test/0427/0514/", 
           outname, ".pdf"), width = 6, height = 4)
ggdensity(data.frame(x=S.D,col=rep("SAPHIRE",length(S.D))), x = "x",
          add = "median", rug = TRUE,
          color = "black", fill = "col",
          palette = c("#EFC000FF")) +
    geom_vline(xintercept = c(REJ1.D, REJ2.D), linetype="dotted",
               color = "darkgrey", size = 1.5)
dev.off()
# text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.1, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "D", xpd = T, cex = 2)
# library(cowplot)
# plot_grid(den1, den2, den3, den4, 
#           labels=c("A", "B", "C", "D"),
#           ncol=2, nrow=2)
# 
# dev.off()

