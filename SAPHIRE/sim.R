# In Hao: b=1.31, In pekar: b=0.385 0.3 less infectious 0.55 more infectious
if (F){
    rm(list = ls())
    
    realData_all <- read.csv("/share/home/jianglab/huangsisi/usr/time-science/Covid19CasesWH.csv", row.names = 1)  
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

simu <- function(b, stage_intervals, init_states, flowN, nstage=6){
    
    days_to_fit <- 1:stage_intervals[[nstage]][["end"]]
    
    ## matrix for results
    states_mat <- matrix(0, length(days_to_fit), length(init_states) + 2)
    states_mat[, 1] <- days_to_fit
    colnames(states_mat) <- c("time", "S", "E", "P", "I", 
                              "A", "H", "R", "Onset_expect")
    
    myold_states <- init_states
    
    n_stage=length(stage_intervals)
    
    for (i_stage in 1:nstage) {
        stage_pars_setings <- c(b = b[i_stage], r = r[i_stage],
                                Dq = Dq[i_stage], n = flowN[i_stage])
        for (d in stage_intervals[[i_stage]][["start"]]:stage_intervals[[i_stage]][["end"]]) {
            states_mat[d, -1] <- update_func(stage_pars = stage_pars_setings, states_old = myold_states)
            myold_states <- states_mat[d, -1]
        }
    }
    
    return(states_mat)
}

indexcase <- function(b0, travel, nmax = 100, retro = 0, pl = 0, ntry = 0){
    if (ntry!=0){
        print(paste("try:",ntry))
    }
    # b0 = 0.385
    # travel = 1
    # b0: transmission rate 1.31, 0.385, 0.3, 0.55
    b = c(b0, 1.31, 1.31, 0.4, 0.17, 0.1)
    # travel: choice for travel component
    {
        if (travel == 1){
            flowN = c(500000, 500000, 800000, 0, 0, 0)
            # print(flowN)
        }else if (travel == 2){
            # no travel component for the whole course
            flowN = rep(0,6)
            # print(flowN)
        }else if (travel == 3){
            # no travel component in late 2019
            flowN = c(0, 500000, 800000, 0, 0, 0)
            # print(flowN)
        }
    }
    
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
    
    if(retro != 0){
        pred_mat <- simu(b, stage_intervals, init_states, flowN)
        pred <- pred_mat[, "Onset_expect"]
    }
    
    if (retro == 0){
        mindist <- 1E11
        minretro <- 0
        for (retro in 1:nmax){

            days_to_fit=1:(68+retro)
            
            stage_intervals=list(
                c(start=1, end=retro),
                c(start=retro+1, end=retro+9),
                c(start=retro+10, end=retro+22),
                c(start=retro+23, end=retro+32),
                c(start=retro+33, end=retro+47),
                c(start=retro+48, end=retro+68)
            )
            pred_mat <- simu(b, stage_intervals, init_states, flowN)
            
            ## successfully established epidemics
            ## On Jan 23, ≥1000 people had become infected and 
            cond1 <- 1e7-pred_mat[retro+23,"S"] >= 1000
            ## ≥1 person was still infectious (I, A, P) at the end of the simulation
            cond2 <- sum(pred_mat[retro+23,c("P","I","A")]) >= 1
            if(!(cond1 & cond2)){next} ## minretro=0 if next always
            
            ## for successfully established epidemics
            ## fit with epidemic growth in Wuhan between **01 January 2020 and 23 January 2020** 
            pred <- pred_mat[, "Onset_expect"]
            dist_cur <- dist(rbind(realData[1:23,1], pred[(retro+1):(retro+23)]))
            if (dist_cur < mindist){
                minretro <- retro
                mindist <- dist_cur
                minpred <- pred
                minmat <- pred_mat
            }
            
        }
        
        retro <- minretro
        pred <- minpred
        pred_mat <- minmat
        
    }
    
    ## plot
    if(pl){
        
        ptime <- 1:(68+retro)
        Sys.setlocale("LC_TIME", "C")
        
        mydate <- as.Date(date_index_case, format = "%b %d") + c(0:(68+retro-1))
        mydate <- format(mydate, "%b %d")
        
        plot(ptime, ptime, ylim = c(0, max(observe) * 1.05), 
             xlab = "", ylab = "", type = "p", col = "white", pch = 16, xaxt="n", cex = 0.5)
        mtext("Onset date (2019-2020)", side = 1, line  = 3, cex = 1.01)
        mtext("No. of ascertained cases", side = 2, line = 3, cex = 1.01)
        axis(1, at = seq(1, 68+retro, 11), labels = mydate[seq(1, 68+retro, 11)])
        #
        abline(v = c(10+retro, 23+retro, 33+retro, 48+retro, 68+retro), lty = 3, lwd = 2, col = "darkgrey")
        text(c(10+retro, 23+retro, 33+retro, 48+retro, 68+retro), par()$usr[4], labels = mydate[c(10+retro, 23+retro, 33+retro, 48+retro, 68+retro)], col = "darkgrey", pos = 3, xpd = T)
        abline(v = 1+retro, lty = 3, lwd = 2)
        text(retro-5, max(pred, observe), expression(2019))
        text(retro+5, max(pred, observe), expression(2020))
        text(retro, par()$usr[4], labels = "Jan 1", pos = 3, xpd = T)
        points(ptime[1:(68+retro)], pred[1:(68+retro)], col = "#BC3C29FF", pch = 16, cex = 0.8)
        points(ptime[1:68+retro], observe, col = "black", pch = 4, cex = 0.8)
        
        legend("topleft", legend = c("Observed", "Simulated"), 
               col = c("black",  "#BC3C29FF"), 
               pch = c(4, 16), bty = "n", y.intersp=0.2)
    }
    return(list(pred = pred, retro = retro, pred_mat = pred_mat))
}

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


pred_retro_list <- list()
nn <- 150
Sys.setlocale("LC_TIME", "C")
for (i in 1:1500){
    pred_retro <- indexcase(b0 = 1.31, travel = 1, nmax = nn, pl = 0, ntry = i)
    pred_retro_list[[i]] <- pred_retro
}

# save(pred_retro_list, 
#      file = "D:/OneDrive - zju.edu.cn/lab/git/time-science-test/SAPHIRE/retro_list.RData")
# load("D:/OneDrive - zju.edu.cn/lab/git/time-science-test/SAPHIRE/retro_list.RData")

