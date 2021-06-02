if(F){
    rm(list = ls())
    library(data.table)
    library(dplyr)
    library(scales)
    realData_all <- read.csv("D:/OneDrive - zju.edu.cn/lab/git/SAPHIRE/data/Covid19CasesWH.csv", row.names = 1)  
    realData <- realData_all[25:46,]
    observe <- realData[,1]
    load("D:/OneDrive - zju.edu.cn/lab/git/time-science-test/SAPHIRE-GEMF/gemf.RData")
    load("D:/OneDrive - zju.edu.cn/lab/git/time-science-test/SAPHIRE-GEMF/e5.RData")
}

pdf(paste0("D:/OneDrive - zju.edu.cn/lab/git/time-science-test/", "fig1-06-02", ".pdf"), width = 9, height = 6)
par(mar = c(4, 5, 2.5, 1))

## plot canvas
ptime <- 1:100
{
    plot(ptime, 1:100, ylim = c(0, max(observe) * 1.05), 
         xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
    mtext("Onset date (2019-2020)", side = 1, line  = 3, cex = 1.01)
    mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
    
    retro <- 58
    mydate <- c(paste("Nov", 4:30), paste("Dec", 1:31), paste("Jan", 1:31), paste("Feb", 1:11))
    axis(1, at = seq(1, 100, 10), labels = mydate[seq(1, 100, 10)])
    abline(v = c(1, 10, 23, 33)+retro, lty = 3, lwd = 2, col = "darkgrey")
    text(c(1, 10, 23, 33)+retro, par()$usr[4], labels = mydate[c(1, 10, 23, 33)+retro], col = "darkgrey", pos = 3, xpd = T)
    abline(v = 1+retro, lty = 3, lwd = 2)
    text(1+retro, par()$usr[4], labels = mydate[1+retro], pos = 3, xpd = T)
    text(retro-3.5, max(observe), expression(2019))
    text(retro+5.5, max(observe), expression(2020))
    
}
glwd=2 # line
gcex=0.3 # points
#### plot
if(T){
    for (i in 1:20){
    data_pred <- data_pred_list[[i]]
    pred <- data_pred$Onset_expect
    idx <- min(which(data_pred[,"delta_I"]>=1), 
               which(data_pred[,"delta_A"]>=1))
    time <- data_pred$Time[idx]
    idx <- data_pred$Time >= time
    x <- data_pred$Time[idx]-time+1
    y <- pred[idx]
    smoothingSpline = smooth.spline(x, y, spar=0.35)
    lines(smoothingSpline, col=alpha("#f7a35c",0.5), lwd=glwd, pch=19)
    # lines(x, y, col=alpha("#BC3C29FF",0.5), lwd=glwd, pch=19)
    points(x, y, col = "#f7a35c", pch = 16, cex = gcex)
}}



for (i in 1:50){
    
    data_pred <- data_pred_list_e5[[i]]
    
    pred <- data_pred$Onset_expect
    idx <- min(which(data_pred[,"delta_I"]>=1), 
               which(data_pred[,"delta_A"]>=1))
    time <- data_pred$Time[idx]
    idx <- data_pred$Time >= time
    x <- data_pred$Time[idx]-time+1
    y <- pred[idx]
    smoothingSpline = smooth.spline(x, y, spar=0.35)
    lines(smoothingSpline, col=alpha("#8085e8",0.5), lwd=glwd, pch=19)
    # lines(x, y, col=alpha("#20854EFF",0.5), lwd=glwd, pch=19)
    points(x, y, col = "#8085e8", pch = 16, cex = gcex)
    
} 


## SAPHIRE
{
    {
        r = c(0.15, 0.15, 0.14, 0.1, 0.16)
        b = c(1.31, 1.31, 0.4, 0.17, 0.1)
        Di = 2.9
        Dp = 2.3
        De = 2.9
        Dq = c(21, 15, 10, 6, 3)
        alpha = 0.55
        Dh = 30
        N = 10000000
        
    }
    
    ## init states and real data
    #      S       E       P       I       A       H       R 
    # 9999021     478     326      34     114      27       0 
    
    {
        R0 <- 0
        H0 <- 27
        
        # realData_all <- read.csv("D:/OneDrive - zju.edu.cn/lab/git/SAPHIRE/data/Covid19CasesWH.csv", row.names = 1)  
        # realData <- realData_all[-c(1:24), ]
        jan1_idx =25
        r0=0.23
        ##
        E0 <- sum(realData_all[(jan1_idx+round(Dp)):(jan1_idx+round(Dp)+round(De)-1),1]) / r0 ## Jan 3-5 for De=2.9 and Dp=2.3
        # E0 <- (40 + 23 + 47) / r0                                              
        P0 <- sum(realData_all[jan1_idx:(jan1_idx+round(Dp)-1),1]) / r0                     ## Jan 1-2 for Dp=2.3
        # P0 <- (41 + 34) / r0
        I0 <- sum(realData_all[(jan1_idx-round(Di)):(jan1_idx-1),1])                             ## Dec 29-31 for Di=2.9
        # I0 <- 11 + 13 + 10                                     
        A0 <- I0 * (1 - r0) / r0
        S0 <- N - E0 - P0 - I0 - A0 - H0 - R0
        init_states <- round(c(S = S0, E = E0, P = P0, I = I0, A = A0, H = H0, R = R0), 0)
    }
    
    days_to_fit=1:68
    stage_intervals=list(
        c(start=1, end=9),
        c(start=10, end=22),
        c(start=23, end=32),
        c(start=33, end=47),
        c(start=48, end=68)
    )
    
    update_func <- function(stage_pars, states_old) {
        alpha = 0.55
        ## stage pars
        b <- stage_pars[1]
        r <- stage_pars[2]
        Dq <- stage_pars[3]
        n <- stage_pars[4]
        ## old states number: c(S, E, P, I, A, H, R)
        S <- states_old[1]
        E <- states_old[2]
        P <- states_old[3]
        I <- states_old[4]
        A <- states_old[5]
        H <- states_old[6]
        R <- states_old[7]
        
        ## S
        ## meaning S->E, S->, S->S
        pS_vec <- c(b * (alpha * P + I + alpha * A) / N, n / N, 1 - b * (alpha * P + I + alpha * A) / N - n / N)
        sample_S <- rmultinom(1, size = S, prob = pS_vec)
        ## E
        ## meaning E->P, E->, E->E
        pE_vec <- c(1 / De, n / N, 1 - 1 / De - n / N)
        sample_E <- rmultinom(1, size = E, prob = pE_vec)
        ## P
        ## meaning P->I, P->A, P->, P->P
        pP_vec <- c(r / Dp, (1 - r) / Dp, n/N, 1 - 1 / Dp - n/N)
        sample_P <- rmultinom(1, size = P, prob = pP_vec)
        ## I
        ## meaning I->H, I->R, I->I
        pI_vec <- c(1 / Dq, 1 / Di, 1 - 1 / Dq - 1 / Di)
        sample_I <- rmultinom(1, size = I, prob = pI_vec)
        ## A
        ## meaning A->R, A->, A->A
        pA_vec <- c(1 / Di, n / N, 1 - 1 / Di - n / N)
        sample_A <- rmultinom(1, size = A, prob = pA_vec)
        ## H
        ## meaning H->R, H->H
        pH_vec <- c(1 / Dh, 1 - 1 / Dh)
        sample_H <- rmultinom(1, size = H, prob = pH_vec)
        ## R
        ## meaning R->, R->R
        pR_vec <- c(n / N, 1 - n / N)
        sample_R <- rmultinom(1, size = R, prob = pR_vec)
        ## new values
        S_new <- sample_S[3] + n
        E_new <- sample_E[3] + sample_S[1]
        P_new <- sample_P[4] + sample_E[1]
        I_new <- sample_I[3] + sample_P[1]
        A_new <- sample_A[3] + sample_P[2]
        H_new <- sample_H[2] + sample_I[1]
        R_new <- sample_R[2] + sample_I[2] + sample_A[1] + sample_H[1]
        Onset_expect <- sample_P[1]
        
        return(c(S_new, E_new, P_new, I_new, A_new, H_new, R_new, Onset_expect))
    }
    
    
    simu <- function(){
        days_to_fit <- 1:22
        
        ## matrix for results
        states_mat <- matrix(0, length(days_to_fit), length(init_states) + 2)
        states_mat[, 1] <- days_to_fit
        colnames(states_mat) <- c("time", "S", "E", "P", "I", "A", "H", "R", "Onset_expect")
        
        myold_states <- init_states
        
        n_stage=length(stage_intervals)
        
        for (i_stage in 1:2) {
            stage_pars_setings <- c(b = b[i_stage], r = r[i_stage],
                                    Dq = Dq[i_stage], n = flowN[i_stage])
            for (d in stage_intervals[[i_stage]][["start"]]:stage_intervals[[i_stage]][["end"]]) {
                states_mat[d, -1] <- update_func(stage_pars = stage_pars_setings, states_old = myold_states)
                myold_states <- states_mat[d, -1]
            }
        }
        
        return(states_mat)
    }
    
}

### SAPHIRE(b=1.31,travel)
{
    flowN = c(500000, 800000, 0, 0, 0)
    estN_mat <- apply(matrix(0, 1000, 2), 1, function(x) 
        simu()[, "Onset_expect"])
    estN_mean <- round(apply(estN_mat, 1, mean), 0)
    estN_up <- round(apply(estN_mat, 1, function(x) quantile(x, 0.975)), 0)
    estN_low <- round(apply(estN_mat, 1, function(x) quantile(x, 0.025)), 0)
    
    polygon(c(ptime[1:22]+retro, rev(ptime[1:22]+retro)), 
            c(estN_up[1:22], rev(estN_low[1:22])), 
            col = "#4DBBD5B2", border = NA)
    points(ptime[1:22]+retro, estN_mean[1:22], 
           col = "#0072B5FF", pch = 17, cex = 0.8)
}
# observed data
points((1:22)+retro, observe, col = "black", pch = 4, cex = 0.8)

mydate <- c(paste("Nov", 4:30), paste("Dec", 1:31), paste("Jan", 1:31), paste("Feb", 1:11))
axis(1, at = seq(1, 100, 10), labels = mydate[seq(1, 100, 10)])
abline(v = c(1, 10, 23, 33)+retro, lty = 3, lwd = 2, col = "darkgrey")
text(c(1, 10, 23, 33)+retro, par()$usr[4], labels = mydate[c(1, 10, 23, 33)+retro], col = "darkgrey", pos = 3, xpd = T)
abline(v = 1+retro, lty = 3, lwd = 2)
text(1+retro, par()$usr[4], labels = mydate[1+retro], pos = 3, xpd = T)
text(retro-3.5, max(observe), expression(2019))
text(retro+5.5, max(observe), expression(2020))

legend("topleft", legend = c("Observed", "SAPHIRE", 
                             "SAPHIRE-GEMF \n(No. of nodes = 1E5)",
                             "SAPHIRE-GEMF \n(No. of nodes = 1E6)"), 
       col = c("black", "#0072B5FF", "#8085e8", "#f7a35c"), pch = c(4, 17, 16, 16), 
       lwd = c(NA, NA, 0.5, 0.5), bty = "n", y.intersp = 2)
dev.off()
