##########################################
# Bayesian Analysis Project: Vinho Verde
# Assessment-3
# Submitted by - Rushmila Islam
##########################################

# Install libraries
library(mcmc)
library(MCMCpack)
library(MASS)
library(mvtnorm)

# Set the seed value
set.seed(1234)

###############################################################################
# Q2. Read the dataset into R and 
# check if there are missing values (NA) and, in case there are, remove them.
###############################################################################

# Read dataset
Wine <- read.csv("~/Library/CloudStorage/OneDrive-Personal/UNSW/2024-02-ZZSC5960-Bayesian Inference and Computation for Data Scientists/Assessment-3/winequality-red.csv", 
                 header = TRUE, sep = ",") 
#Wine 
View(Wine)
summary(Wine)

# Find missing values
print("Position of missing values ")
which(is.na(Wine))

# count total missing values 
print("Count of total missing values  ")
sum(is.na(Wine))

###############################################################################
# Q3. Consider "good" wine with quality above 6.5 (included)
###############################################################################

# Create a new variable good_wine for response variable with values 0, 1 
Wine$good_wine <- ifelse(Wine$quality >= 6.5, 1, 0)
View(Wine)
table(Wine$good_wine)


###############################################################################
#Q4. Run frequentist analysis on the logistic model, using glm() function.
###############################################################################

# Logistic model using glm()
model <- glm(Wine$good_wine ~ Wine$fixed.acidity + Wine$volatile.acidity + 
               Wine$citric.acid + Wine$residual.sugar + Wine$chlorides + 
               Wine$free.sulfur.dioxide + Wine$total.sulfur.dioxide + 
               Wine$density + Wine$pH + Wine$sulphates + Wine$alcohol, 
               data = Wine, family = binomial(link="logit"))

# Summary of the model
summary(model)

###############################################################################
# Q5. 
# Estimate the probabilities of having a "success": fix each covariate at its mean level, 
# and compute the probabilities for a wine to score "good" varying 
# total.sulfur.dioxide and plot the results.
###############################################################################
total_sulfur_dioxide_range <- Wine$total.sulfur.dioxide

total_sulfur_dioxide_range1 <- seq(from=min(Wine$total.sulfur.dioxide), 
                                   to=max(Wine$total.sulfur.dioxide), length.out = 1000)

total_sulfur_dioxide_range1

length(total_sulfur_dioxide_range)


fixed_acidity_mean <- mean(Wine$volatile.acidity) 
volatile_acidity_mean <- mean(Wine$volatile.acidity) 
citric_acid_mean <- mean(Wine$citric.acid)
residual_sugar_mean <- mean(Wine$residual.sugar)
chlorides_mean <-  mean(Wine$chlorides) 
free_sulfur_dioxide_mean <- mean(Wine$free.sulfur.dioxide)
total_sulfur_dioxide_mean <- mean(Wine$total.sulfur.dioxide)
density_mean <- mean(Wine$density)
pH_mean <- mean(Wine$pH)
sulphates_mean <- mean(Wine$sulphates)
alcohol_mean <- mean(Wine$alcohol)

b0 <- model$coef[1]
b1 <- model$coef[2]
b2 <- model$coef[3]
b3 <- model$coef[4]
b4 <- model$coef[5]
b5 <- model$coef[6]
b6 <- model$coef[7]
b7 <- model$coef[8]
b8 <- model$coef[9]
b9 <- model$coef[10]
b10 <- model$coef[11]
b11 <- model$coef[12]



total_sulfur_dioxide_model <- b0 + b1*fixed_acidity_mean + b2*volatile_acidity_mean +
  b3*citric_acid_mean + b4*residual_sugar_mean + b5*chlorides_mean +
  b6*free_sulfur_dioxide_mean + b7*total_sulfur_dioxide_range + 
  b8*density_mean + b9*pH_mean + b10*sulphates_mean + b11*alcohol_mean

length(total_sulfur_dioxide_model)

total_sulfur_dioxide_model1 <- b0 + b1*fixed_acidity_mean + b2*volatile_acidity_mean +
  b3*citric_acid_mean + b4*residual_sugar_mean + b5*chlorides_mean +
  b6*free_sulfur_dioxide_mean + b7*total_sulfur_dioxide_range1 + 
  b8*density_mean + b9*pH_mean + b10*sulphates_mean + b11*alcohol_mean

length(total_sulfur_dioxide_model1)


total_sulfur_dioxide_probs <- exp(total_sulfur_dioxide_model)/(1 + exp(total_sulfur_dioxide_model))
total_sulfur_dioxide_probs

total_sulfur_dioxide_probs1 <- exp(total_sulfur_dioxide_model1)/(1 + exp(total_sulfur_dioxide_model1))

length(total_sulfur_dioxide_probs1)
length(total_sulfur_dioxide_range1)

par(mfrow=c(1,1))

plot(total_sulfur_dioxide_range, total_sulfur_dioxide_probs, 
     #ylim=c(0,1),
     #type="l", 
     lwd=3, 
     lty=2, 
     col="cyan", 
     xlab="total.sulfur.dioxide", ylab="P(success)", main="Probabilities for a good wine varying total.sulfur.dioxide")

par(mfrow=c(1,1))

plot(total_sulfur_dioxide_range1, total_sulfur_dioxide_probs1, 
     #ylim=c(0,1),
     type="l", 
     #lwd=2, 
     lty=1, 
     col="blue", 
     xlab="total.sulfur.dioxide", ylab="P(success)", 
     main="Probabilities for a good wine varying total.sulfur.dioxide")

###############################################################################
# Q6. Perform Bayesian analysis of the logistic model 
# i.e. approximate the posterior distributions of the regression coefficients.
###############################################################################

# Write R function for the log posterior distribution
lpost.LR <- function(beta,x,y)
{
  eta <- as.numeric(x %*% beta)
  
  logp <- eta - log(1+exp(eta))
  logq <- log(1-exp(logp))
  
  #log-likelihood 
  logl <- sum(logp[y==1]) + sum(logq[y==0]) 
  
  #prior
  lprior <- sum(dnorm(beta,0,10,log=T)) 
  
  #posterior
  posterior <- logl + lprior
  
  return(posterior)
}

# Fix number of simulations
S <- 10^4

X1=cbind(rep(1,nrow(Wine)),Wine$fixed.acidity, Wine$volatile.acidity, Wine$citric.acid,
         Wine$residual.sugar, Wine$chlorides, Wine$free.sulfur.dioxide, Wine$total.sulfur.dioxide,
         Wine$density, Wine$pH, Wine$sulphates, Wine$alcohol)

X <- X1[1,]

y <- Wine$good_wine[1]


# Checking for NA
for(i in 1:nrow(X1))
{
  if(sum(is.na(X1[i,]))==0){
    X <- rbind(X,X1[i,])
    y <- c(y,Wine$good_wine[i])
  }  
}



# Create beta matrix
beta_mat <- matrix(NA,nrow=S,ncol=ncol(X))

# Initialization
beta_mat[1,] <- as.numeric(coefficients(model))


# For prediction
y_new <- c(1)
x_new <- c(1,7.5,0.6,0.0,1.70,0.085,5,45,0.9965,3.40,0.63,12)

Omega_prop <- solve(t(X) %*% X)
k <- ncol(beta_mat)
acc <- 0

for(iter in 2:S)
{
  # 1. Propose a new set of values
  beta_star <- rmvnorm(1,beta_mat[iter-1,], 0.7*Omega_prop) 
  
  # 2. Compute the posterior density on the proposed value and on the old value  
  newpost=lpost.LR(t(beta_star),X,y)
  oldpost=lpost.LR(matrix(beta_mat[iter-1,],ncol=1),X,y)
  
  # 3. Acceptance step
  if(runif(1,0,1)>exp(newpost-oldpost)){
    beta_mat[iter,]=beta_mat[iter-1,]
  } else{
    beta_mat[iter,]=beta_star
    acc=acc+1
  }
  # 4. Print the stage of the chain
  if(iter%%1000==0){print(c(iter,acc/iter))}
  
  # 5. Prediction 
  p_new <- exp(sum(beta_mat[iter,] * x_new) ) / (1 + exp(sum(beta_mat[iter,] * x_new) ))
  y_new[iter] <- rbinom(1,1,prob=p_new)

}

acc/iter

#Plotting Markov chains with MLE initialization
par(mfrow=c(4,3))
plot(beta_mat[,1],type="l", ylab=expression(beta[0])) 
abline(h=model$coefficients[1],col="red",lty=2)
plot(beta_mat[,2],type="l", ylab=expression(beta[1]))
abline(h=model$coefficients[2],col="red",lty=2)
plot(beta_mat[,3],type="l", ylab=expression(beta[2]))
abline(h=model$coefficients[3],col="red",lty=2)
plot(beta_mat[,4],type="l", ylab=expression(beta[3]))
abline(h=model$coefficients[4],col="red",lty=2)
plot(beta_mat[,5],type="l", ylab=expression(beta[4]))
abline(h=model$coefficients[5],col="red",lty=2)
plot(beta_mat[,6],type="l", ylab=expression(beta[5]))
abline(h=model$coefficients[6],col="red",lty=2)
plot(beta_mat[,7],type="l", ylab=expression(beta[6]))
abline(h=model$coefficients[7],col="red",lty=2)
plot(beta_mat[,8],type="l", ylab=expression(beta[7]))
abline(h=model$coefficients[8],col="red",lty=2)
plot(beta_mat[,9],type="l", ylab=expression(beta[8]))
abline(h=model$coefficients[9],col="red",lty=2)
plot(beta_mat[,10],type="l", ylab=expression(beta[9]))
abline(h=model$coefficients[10],col="red",lty=2)
plot(beta_mat[,11],type="l", ylab=expression(beta[10]))
abline(h=model$coefficients[11],col="red",lty=2)
plot(beta_mat[,12],type="l", ylab=expression(beta[11]))
abline(h=model$coefficients[12],col="red",lty=2)


###############################################################################
# Q6.3 Choose 4 different initialization coefficients.
###############################################################################
# Simulation with 4 initializations

# Initializing with 4 different set of coefficients
init1 <- as.numeric(coefficients(model)) 
init2 <- as.numeric(coefficients(model)*0.5) 
init3 <- as.numeric(coefficients(model)*0.33) 
init4 <- as.numeric(coefficients(model)*0) 

# For prediction
y_new <- c(1)
x_new <- c(1, 7.5, 0.6, 0.0, 1.70, 0.085, 5,45, 0.9965, 3.40, 0.63, 12)

# Simulation function
run_sim <- function(S, init, X, y, c){
  
  beta <- matrix(NA, nrow=S, ncol=ncol(X))
  beta[1,] <- as.numeric(init)
  
  Omega_prop <- solve(t(X) %*% X)

  y_new <- c()
  acc <- 0
  
  for(i in 2:S)
  {
    # 1. Propose a new set of values
    beta_star <- rmvnorm(1, beta[i-1,], c*Omega_prop)
    
    # 2. Compute the posterior density on the proposed value and on the old value  
    newpost=lpost.LR(t(beta_star), X, y)
    oldpost=lpost.LR(matrix(beta[i-1,], ncol=1), X, y)
    
    # 3. Acceptance step
    if(runif(1,0,1)>exp(newpost-oldpost)){
      beta[i,]=beta[i-1,]
    }else{
      beta[i,]=beta_star
      acc=acc+1
    }
    
    # 4. Print the stage of the chain
    if(i%%1000==0){print(c(i, acc/i))}
    acc_ratio <- acc/i
    
    # 5. Prediction
    p_new <- exp(sum(beta[i-1,] * x_new) ) / (1 + exp(sum(beta[i-1,] * x_new) ))
    y_new[i-1] <- rbinom(1, 1, prob=p_new)
  }
  
  
  return(list(beta=beta, y_new=y_new, acc_ratio=acc_ratio))
}



# Run four simulations with four different initialization for the coefficients and constant values

run1 <- run_sim(S, init1, X, y, 0.7) 
run1$acc_ratio
run2 <- run_sim(S, init2, X, y, 0.7) 
run2$acc_ratio
run3 <- run_sim(S, init3, X, y, 0.7) 
run3$acc_ratio
run4 <- run_sim(S, init4, X, y, 0.7) 
run4$acc_ratio


pred_run1 <- table(run1$y_new)
pred_run2 <- table(run2$y_new)
pred_run3 <- table(run3$y_new)
pred_run4 <- table(run4$y_new)


acc_ratios <- c(run1$acc_ratio, run2$acc_ratio, run3$acc_ratio, run4$acc_ratio)
acc_ratios

# Plotting the chains for each coefficients
#par(mfrow=c(4,3))
par(mfrow=c(1,1))
plotting <- function(run1, run2, run3, run4){
  
  for(i in 1:ncol(X1)){
    #par(mfrow=c(4,3))
    plot(run1$beta[, i], type = 'l', 
         ylab = expression(beta[i]), 
         main = paste('Plot of i=', i),
         ylim = c(min(min(run1$beta[, i]), min(run2$beta[, i]), 
                      min(run3$beta[, i]), min(run4$beta[, i])),
                      max(max(run1$beta[, i]), max(run2$beta[, i]), 
                      max(run3$beta[, i]), max(run4$beta[, i]))))
    lines(run2$beta[, i], type = "l", col = "green")
    lines(run3$beta[, i], type = "l", col = "orange")
    lines(run4$beta[, i], type = "l", col = "blue")
    abline(h=model$coefficients[i], col="red", lty=2)
    legend('topright', 
           legend = c('init1:MLE', 'init2', 'init3', 'init4:0s'),
           col = c('black', 'green', 'orange', 'blue'),
           lty=1,
           cex=0.8)
  }
}

plotting(run1, run2, run3, run4)

par(mfrow=c(1,1))
plot(beta_mat[2000:10000,4],type="l", ylab=expression(beta[3]))
     #xlim = c(2000, 8000), ylim = c(-4,-2))
abline(h=model$coefficients[4],col="red",lty=2)
###############################################################################
# Q7.
###############################################################################

# Probabilities of given data
prob_run1 <- table(run1$y_new)[2]/table(run1$y_new)[1]
prob_run1
prob_run2 <- table(run2$y_new)[2]/table(run2$y_new)[1]
prob_run2
prob_run3 <- table(run3$y_new)[2]/table(run3$y_new)[1]
prob_run3
prob_run4 <- table(run4$y_new)[2]/table(run4$y_new)[1]
prob_run4

probs <- c(prob_run1, prob_run2, prob_run3, prob_run4)
probs

#Barplots
par(mfrow=c(2,2))
barplot(pred_run1, main="MLE", ylab="count", xlab="y")
barplot(pred_run2, main="init2", ylab="count", xlab="y")
barplot(pred_run3, main="init3", ylab="count", xlab="y")
barplot(pred_run4, main="all 0", ylab="count", xlab="y")

#PDF
par(mfrow=c(2,2))
# Visualize the simulated data
hist(run1$y_new, freq = FALSE, 
     main = "Posterior Predictive Distribution run1 ", xlab = "Simulated Data")
lines(density(run1$y_new), col = "blue") 
hist(run2$y_new, freq = FALSE, 
     main = "Posterior Predictive Distribution run2 ", xlab = "Simulated Data")
lines(density(run2$y_new), col = "blue")  
hist(run3$y_new, freq = FALSE, 
     main = "Posterior Predictive Distribution run3 ", xlab = "Simulated Data")
lines(density(run3$y_new), col = "blue")  
hist(run4$y_new, freq = FALSE, 
     main = "Posterior Predictive Distribution run4 ", xlab = "Simulated Data")
lines(density(run4$y_new), col = "blue")  

par(mfrow=c(1,1))
plot(table(run1$y_new[2000:10000]))

par(mfrow=c(1,1))
# Visualize the simulated data
plot(density(run1$y_new), col = "magenta",
     main = "Posterior Predictive Distribution (initialized with MLE)")  # Add density plot for visualization


###############################################################################
# Q8. Metrop()
###############################################################################
x0 <- X
y0 <- Wine$good_wine[1]
Omega_prop0 <- solve(t(x0) %*% x0)


m_out1 <- metrop(obj=lpost.LR, initial=init1, nbatch=S, 
                 #x=x0, y=y0, scale=0.02*Omega_prop0)
                 x=x0, y=y0, scale=6)
m_out2 <- metrop(obj=lpost.LR, initial=init2, nbatch=S, 
                 #x=x0, y=y0, scale=0.02*Omega_prop0)
                 x=x0, y=y0, scale=6)

m_out3 <- metrop(obj=lpost.LR, initial=init3, nbatch=S, 
                 #x=x0, y=y0, scale=0.02*Omega_prop0) #0.023
                 x=x0, y=y0, scale=6)

m_out4 <- metrop(obj=lpost.LR, initial=init4, nbatch=S, 
                 #x=x0, y=y0, scale=0.02*Omega_prop0)
                 x=x0, y=y0, scale=6)


# Print summary of the MCMC output
(m_out1$accept) 
(m_out2$accept) 
(m_out3$accept) 
(m_out4$accept) 

# Plot ACF
plot_acf <- function(out){
  par(mfrow=c(4,3))
  for (i in 1:ncol(out$batch)) {
    acf(out$batch[ , i], lag.max = 5000)
  }
}


plot_acf(m_out1)
plot_acf(m_out2)
plot_acf(m_out3)
plot_acf(m_out4)
#ncol(out1$batch)

acf(run1$beta[, 1], lag.max = 5000)

#run2$beta[, i]
plot_acf1 <- function(out){
  par(mfrow=c(4,3))
  for (i in 1:ncol(out$beta)) {
    acf(out$beta[ , i], lag.max = 5000)
  }
}

plot_acf1(run1)
plot_acf1(run2)
plot_acf1(run3)
plot_acf1(run4)


# Plotting metrop() output for a single run
poltting_mcmc <- function(m_out1){
  par(mfrow=c(4,3))
  plot(m_out1$batch[ , 1], type="l", ylab=expression(beta[0]))
  abline(h=model$coefficients[1],col="red",lty=2)
  plot(m_out1$batch[ , 2], type="l", ylab=expression(beta[1]))
  abline(h=model$coefficients[2],col="red",lty=2)
  plot(m_out1$batch[ , 3], type="l", ylab=expression(beta[3]))
  abline(h=model$coefficients[3],col="red",lty=2)
  plot(m_out1$batch[ , 4], type="l", ylab=expression(beta[4]))
  abline(h=model$coefficients[4],col="red",lty=2)
  plot(m_out1$batch[ , 5], type="l", ylab=expression(beta[5]))
  abline(h=model$coefficients[5],col="red",lty=2)
  plot(m_out1$batch[ , 6], type="l", ylab=expression(beta[6]))
  abline(h=model$coefficients[6],col="red",lty=2)
  plot(m_out1$batch[ , 7], type="l", ylab=expression(beta[7]))
  abline(h=model$coefficients[7],col="red",lty=2)
  plot(m_out1$batch[ , 8], type="l", ylab=expression(beta[8]))
  abline(h=model$coefficients[8],col="red",lty=2)
  plot(m_out1$batch[ , 9], type="l", ylab=expression(beta[9]))
  abline(h=model$coefficients[9],col="red",lty=2)
  plot(m_out1$batch[ , 10], type="l", ylab=expression(beta[10]))
  abline(h=model$coefficients[10],col="red",lty=2)
  plot(m_out1$batch[ , 11], type="l", ylab=expression(beta[11]))
  abline(h=model$coefficients[11],col="red",lty=2)
  plot(m_out1$batch[ , 12], type="l", ylab=expression(beta[12]))
  abline(h=model$coefficients[12],col="red",lty=2)
}


poltting_mcmc(m_out1)
poltting_mcmc(m_out2)
poltting_mcmc(m_out3)
poltting_mcmc(m_out4)

# Plotting metrop() output for 4 initialization 
par(mfrow=c(4,3))
plotting_metrop <- function(m_out1, m_out2, m_out3, m_out4){
  
  for(i in 1:ncol(X1)){
    plot(m_out1$batch[, i], type = 'l', 
         ylab = expression(beta[i]), 
         main = paste('Metrop Plot of i=', i),
         ylim = c(min(min(m_out1$batch[, i]), min(m_out2$batch[, i]), 
                      min(m_out3$batch[, i]), min(m_out4$batch[, i])),
                      max(max(m_out1$batch[, i]), max(m_out2$batch[, i]), 
                      max(m_out3$batch[, i]), max(m_out4$batch[, i]))))
    lines(m_out2$batch[, i], type = "l", col = "green")
    lines(m_out3$batch[, i], type = "l", col = "orange")
    lines(m_out4$batch[, i], type = "l", col = "blue")
    abline(h=model$coefficients[i], col="red", lty=2)
    legend('topright', 
           legend = c('init1:MLE', 'init2', 'init3', 'init4:0s'),
           col = c('black', 'green', 'orange', 'blue'),
           lty=1,
           cex=0.2)
  }
}

plotting_metrop(m_out1, m_out2, m_out3, m_out4)


###############################################################################
# Visualizations of data
###############################################################################
par(mfrow=c(3,4))
plot(density(Wine$good_wine), col = "orange")  # Add density plot for visualization
plot(density(Wine$fixed.acidity), col = "grey")  # Add density plot for visualization
plot(density(Wine$volatile.acidity), col = "grey")  # Add density plot for visualization
plot(density(Wine$citric.acid), col = "grey")  # Add density plot for visualization
plot(density(Wine$residual.sugar), col = "grey")  # Add density plot for visualization
plot(density(Wine$chlorides), col = "grey")  # Add density plot for visualization
plot(density(Wine$free.sulfur.dioxide), col = "grey")  # Add density plot for visualization
plot(density(Wine$total.sulfur.dioxide), col = "grey")  
plot(density(Wine$density), col = "grey") 
plot(density(Wine$pH), col = "grey") 
plot(density(Wine$sulphates), col = "grey") 
plot(density(Wine$alcohol), col = "grey") 

par(mfrow=c(3,4))
plot(density(log(Wine$good_wine)), col = "orange")  
plot(density(log(Wine$fixed.acidity)), col = "grey")  
plot(density(log(Wine$volatile.acidity)), col = "grey")  
plot(density(log(Wine$citric.acid)), col = "grey")  
plot(density(log(Wine$residual.sugar)), col = "grey")  
plot(density(log(Wine$chlorides)), col = "grey")  
plot(density(log(Wine$free.sulfur.dioxide)), col = "grey")  
plot(density(log(Wine$total.sulfur.dioxide)), col = "grey")  
plot(density(log(Wine$density)), col = "grey") 
plot(density(log(Wine$pH)), col = "grey") 
plot(density(log(Wine$sulphates)), col = "grey") 
plot(density(log(Wine$alcohol)), col = "grey") 



