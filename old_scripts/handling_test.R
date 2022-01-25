cr_RM_pred <- function(time, y, p){
  N <- y[1]; P <- y[2]
  with(as.list(p),{
    dN <- r*N*(1-N/K) - a*N*P/(1+a*h*N)
    dP <- rp*P*(1-P/Kp)
  return(list(c(dN, dP)))
  })
}

t <- 0:100; y0 <-c(N=50, P=50)

handling_trials = c(0.005, 0.01, 0.05, 0.15, 0.25)
equilibrium = matrix(NA, nrow = length(handling_trials), ncol = 2)
trajectories = matrix(NA, nrow = length(handling_trials), ncol = 101)
colnames(equilibrium) = c("N", "P")
for (i in 1:length(handling_trials)) {
p <- list(r=1.2, K=100, a=.001, rp = 1.4, h=handling_trials[i], Kp = 150)
outdf <- as.data.frame( ode(y0, t, cr_RM_pred,  p) )
equilibrium[i,1] = outdf$N[100]
equilibrium[i,2] = outdf$P[100]
trajectories[i,] = outdf$N
}

jpeg('plots/handling_test_equilvals.jpg')
plot(equilibrium[,1])
dev.off()
jpeg('plots/handling_test.jpg')
matplot(t(trajectories), type = "l", ylim = c(84,90), xlim = c(0,120), ylab = "abundance", 
        xlab = "time", lty = 1, col = c("red", "orange", "green", "blue", "black"))
legend(100, 90, legend=c("0.005", "0.01", "0.05", "0.15", "0.25"),
       col=c("red", "orange", "green", "blue", "black"), lty=1, cex=0.8)
dev.off()

vary_handling_test_percdiff = (trajectories[length(handling_trials),length(t)] - trajectories[2,length(t)]) / trajectories[1,length(t)] * 100
 # 0.29 %

