## parameters
rm(list = ls())
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
    flowN = c(500000, 800000, 0, 0, 0)
}

## init states and real data
#      S       E       P       I       A       H       R 
# 9999021     478     326      34     114      27       0 
## realData     : real data from the CDC

{
    
    ## N            : population size
    ## H0           : initial number of hospitalized cases based on the reports
    ## R0           : initial number of removed individuals
    ## De           : latent period
    ## r0           : initial ascertainment rate
    
    
    
    R0 <- 0
    H0 <- 27
    
    realData_all <- read.csv("D:/OneDrive - zju.edu.cn/lab/git/SAPHIRE/data/Covid19CasesWH.csv", row.names = 1)  
    realData <- realData_all[-c(1:24), ]
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
    ## new values
    S_new <- S - b * S * (alpha * P + I + alpha * A) / N + n - n * S / N
    E_new <- E + b * S * (alpha * P + I + alpha * A) / N - E / De - n * E / N
    P_new <- P +  E / De  - P / Dp - n * P / N
    I_new <- I + r * P / Dp - I / Di - I / Dq
    A_new <- A + (1 - r) * P / Dp - A / Di - n * A / N
    H_new <- H + I / Dq - H / Dh
    R_new <- R + H / Dh + (A + I) / Di - n * R / N
    Onset_expect <- r * P / Dp
    ##
    return(c(S_new, E_new, P_new, I_new, A_new, H_new, R_new, Onset_expect))
}

## matrix for results
states_mat <- matrix(0, length(days_to_fit), length(init_states) + 2)
states_mat[, 1] <- days_to_fit
colnames(states_mat) <- c("time", "S", "E", "P", "I", "A", "H", "R", "Onset_expect")

myold_states <- init_states

n_stage=length(stage_intervals)

for (i_stage in 1:5) {
    stage_pars_setings <- c(b = b[i_stage], r = r[i_stage], Dq = Dq[i_stage], n = flowN[i_stage])
    for (d in stage_intervals[[i_stage]][["start"]]:stage_intervals[[i_stage]][["end"]]) {
        states_mat[d, -1] <- update_func(stage_pars = stage_pars_setings, states_old = myold_states)
        myold_states <- states_mat[d, -1]
    }
}

observe <- realData[,1]

ptime <- 1:68
mydate <- c(paste("Jan", 1:31), paste("Feb", 1:29), paste("Mar", 1:8))
pred0 <- states_mat[,9]
# plot(ptime, pred, ylim = c(0, max(pred,observe) * 1.05), 
#      xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
# mtext("Onset date (2020)", side = 1, line  = 3, cex = 1.01)
# mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
# axis(1, at = seq(1, 68, 11), labels = mydate[seq(1, 68, 11)])
# #
# abline(v = c(10, 23, 33, 48, 68), lty = 3, lwd = 2, col = "darkgrey")
# text(c(10, 23, 33, 48, 68), par()$usr[4], labels = mydate[c(10, 23, 33, 48, 68)], col = "darkgrey", pos = 3, xpd = T)
# points(ptime[1:68], pred[1:68], col = "#BC3C29FF", pch = 16, cex = 0.8)
# points(ptime, observe, col = "black", pch = 4, cex = 0.8)
# 
# legend("topleft", legend = c("Observed", "Simulated"), col = c("black",  "#BC3C29FF", "#0072B5FF"), pch = c(4, 16, 17), bty = "n")
# text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "A", xpd = T, cex = 2)

##################################################################################
### No traveler component
{
    update_func <- function(stage_pars, states_old) {
        alpha = 0.55
        ## stage pars
        b <- stage_pars[1]
        r <- stage_pars[2]
        Dq <- stage_pars[3]
        n <- 0
        ## old states number: c(S, E, P, I, A, H, R)
        S <- states_old[1]
        E <- states_old[2]
        P <- states_old[3]
        I <- states_old[4]
        A <- states_old[5]
        H <- states_old[6]
        R <- states_old[7]
        ## new values
        S_new <- S - b * S * (alpha * P + I + alpha * A) / N + n - n * S / N
        E_new <- E + b * S * (alpha * P + I + alpha * A) / N - E / De - n * E / N
        P_new <- P +  E / De  - P / Dp - n * P / N
        I_new <- I + r * P / Dp - I / Di - I / Dq
        A_new <- A + (1 - r) * P / Dp - A / Di - n * A / N
        H_new <- H + I / Dq - H / Dh
        R_new <- R + H / Dh + (A + I) / Di - n * R / N
        Onset_expect <- r * P / Dp
        ##
        return(c(S_new, E_new, P_new, I_new, A_new, H_new, R_new, Onset_expect))
    }
    
    ## matrix for results
    states_mat <- matrix(0, length(days_to_fit), length(init_states) + 2)
    states_mat[, 1] <- days_to_fit
    colnames(states_mat) <- c("time", "S", "E", "P", "I", "A", "H", "R", "Onset_expect")
    
    myold_states <- init_states
    
    n_stage=length(stage_intervals)
    
    for (i_stage in 1:5) {
        stage_pars_setings <- c(b = b[i_stage], r = r[i_stage], Dq = Dq[i_stage], n = flowN[i_stage])
        for (d in stage_intervals[[i_stage]][["start"]]:stage_intervals[[i_stage]][["end"]]) {
            states_mat[d, -1] <- update_func(stage_pars = stage_pars_setings, states_old = myold_states)
            myold_states <- states_mat[d, -1]
        }
    }
    
    
 
    observe <- realData[,1]
    
    ptime <- 1:68
    mydate <- c(paste("Jan", 1:31), paste("Feb", 1:29), paste("Mar", 1:8))
    pred <- states_mat[,9]
    plot(ptime, pred, ylim = c(0, max(pred,observe) * 1.05), 
         xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
    mtext("Onset date (2020)", side = 1, line  = 3, cex = 1.01)
    mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
    axis(1, at = seq(1, 68, 11), labels = mydate[seq(1, 68, 11)])
    #
    abline(v = c(10, 23, 33, 48, 68), lty = 3, lwd = 2, col = "darkgrey")
    text(c(10, 23, 33, 48, 68), par()$usr[4], labels = mydate[c(10, 23, 33, 48, 68)], col = "darkgrey", pos = 3, xpd = T)
    points(ptime[1:68], pred0[1:68], col = "#BC3C29FF", pch = 16, cex = 0.8)
    points(ptime[1:68], pred[1:68], col = "#0072B5FF", pch = 17, cex = 0.8)
    points(ptime, observe, col = "black", pch = 4, cex = 0.8)
    
    legend("topleft", legend = c("Observed", "Simulated", "Simulated(No travel component)"), col = c("black",  "#BC3C29FF", "#0072B5FF"), pch = c(4, 16, 17), bty = "n")
    # text(par()$usr[1] - (par()$usr[2] -par()$usr[1]) * 0.12, par()$usr[4] + (par()$usr[4] - par()$usr[3]) * 0.06, labels = "A", xpd = T, cex = 2)
}








