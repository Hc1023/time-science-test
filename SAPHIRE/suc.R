# In Hao: b=1.31, In pekar: b=0.385 0.3 less infectious 0.55 more infectious
if (F){
    rm(list = ls())
    
    realData_all <- read.csv("D:/OneDrive - zju.edu.cn/lab/git/SAPHIRE/data/Covid19CasesWH.csv", row.names = 1)  
    realData <- realData_all[-c(1:24), ]
    observe <- realData[,1]
    N = 10000000
    Di = 2.9
    Dp = 2.3
    De = 2.9
    alpha = 0.55
    Dh = 30
    r = c(0.15, 0.15, 0.15, 0.14, 0.1, 0.16)
    Dq = c(21, 21, 15, 10, 6, 3)
}

update_func <- function(stage_pars, states_old) {
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

simu <- function(b, init_states, flowN){
    
    
    ## matrix for results
    states_mat <- matrix(0, 100, length(init_states) + 2)
    states_mat[, 1] <- 1:100
    colnames(states_mat) <- c("time", "S", "E", "P", "I", 
                              "A", "H", "R", "Onset_expect")
    
    myold_states <- init_states
    stage_pars_setings <- c(b = b, r = 0.15, Dq = 21, n = 50000)
    for (d in 1:100) {
        states_mat[d, -1] <- update_func(stage_pars = stage_pars_setings, states_old = myold_states)
        myold_states <- states_mat[d, -1]
    }
    
    return(states_mat)
}

forwardsim <- function(b0, pl = 0, ntry = 0){
    if (ntry!=0){
        print(paste("try:",ntry))
        
    }
    # b0 = 0.385
    # travel = 1
    # b0: transmission rate 1.31, 0.385, 0.3, 0.55
    b = b0
    flowN = 500000
    
    # initial value
    {
        R0 <- 0
        H0 <- 0
        
        E0 <- 1 
        P0 <- 0
        I0 <- 0                                    
        A0 <- 0
        S0 <- N - E0 - P0 - I0 - A0 - H0 - R0
        init_states <- round(c(S = S0, E = E0, P = P0, I = I0, A = A0, H = H0, R = R0), 0)
    }
    
    
    pred_mat <- simu(b, init_states, flowN)
    pred <- pred_mat[, "Onset_expect"]
    
    ## plot
    if(pl){
        
        ptime <- 1:100
        plot(ptime, ptime, ylim = c(0, max(pred) * 1.05), 
             xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
        mtext("Onset date (2019-2020)", side = 1, line  = 3, cex = 1.01)
        mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
        axis(ptime)
        points(ptime[1:(22+retro)], pred[1:(22+retro)], col = "#BC3C29FF", pch = 16, cex = 0.8)
    }
    return(list(pred = pred, retro = retro, pred_mat = pred_mat))
}

p <- forwardsim(b0 = 1.31, pl = 1)

suc <- 0
for (i in 1:5000){
    p <- forwardsim(b0 = 1.31, pl = 0)
    s <- p[["pred_mat"]][72,]
    # successfully established dynamics
    if(((1e7 - s["S"]) >= 1000) & (sum(s[c("P","I","A")]) >= 1)){
        suc = suc+1
    }
        
}
# 
# > suc
# [1] 4488

