#librairies needed 
library(ggplot2)

library(car)
library(agricolae)

pros=read.table("http://statweb.stanford.edu/~tibs/ElemStatLearn/datasets/prostate.data")
pros
summary(pros)




lpsa <-pros$lpsa

ggplot(pros, aes(x = lpsa)) +
  geom_boxplot( ) +
  labs(x = "lpsa", y = "Frequency", title = "Distribution de logarithme PSA")



summary(pros$lpsa)
ggplot(pros, aes(x = lpsa)) +
  geom_boxplot(fill = "steelblue", color = "black") +
  labs(x = "Lpsa", title = "Boxplot de Lpsa") + scale_x_continuous(breaks = seq(-4, 15, by = 0.5))

summary(pros$age)

library(ggplot2)

# Assuming you have a dataframe called 'pros' with 'gleason' and 'lpsa' columns

ggplot(pros, aes(x = factor(gleason), y = lpsa, fill = factor(gleason))) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA) +
  xlab("Gleason") +
  ylab("lpsa") +
  ggtitle("Distribution of lpsa by Gleason") +
  theme_minimal()

# Correlation matrix
Matcor=cor(pros[,1:9])
Matcor

# Effect of gleason variable  on target 


pros$gleason <- as.factor(pros$gleason)



# Perform one-way ANOVA
anova_result <- aov(lpsa ~ gleason, data = pros)

# Summarize the ANOVA results

summary_anova<-summary(anova_result,alpha=0.05)

print(summary_anova)


# Perform Tukey's HSD test for multiple comparisons
tukey_result <- TukeyHSD(anova_result)

# Print the comparison results
print(tukey_result)



#Replace values in the 'gleason'  variable

pros$gleason<-as.numeric(pros$gleason)

pros$gleason <- replace(pros$gleason, pros$gleason == 6, 1)
pros$gleason <- replace(pros$gleason, pros$gleason == 7, 2)
pros$gleason <- replace(pros$gleason, pros$gleason == 8, 3)
pros$gleason <- replace(pros$gleason, pros$gleason == 9, 3)

#la nouvelle corrélation entre les deux variables
cor(pros$lpsa,pros$gleason)








# pairplot
ggpairs(pros[,1:9])


# Scatterplot lcavol vs lpsa

ggplot(pros, aes(x = lcavol, y = lpsa)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  xlab("lcavol") +
  ylab("lpsa") +
  ggtitle("lcavol vs. lpsa") +
  theme_minimal() +
  geom_text(aes(x = min(lcavol), y = max(lpsa), label = paste("Equation: ", round(coef(summary(lm(lpsa ~ lcavol, data = pros)))[1, 1], 2), " + ", round(coef(summary(lm(lpsa ~ lcavol, data = pros)))[2, 1], 2), " * lcavol", "\n", "R-squared: ", round(summary(lm(lpsa ~ lcavol, data = pros))$r.squared, 2))), 
            hjust = 0, vjust = 1, color = "black", parse = TRUE)
#scale the data


scaled_data <- scale(pros)

# View the standardized data
print(scaled_data)

# Drop the train variable from scaled_data
scaled_data <- subset(scaled_data, select = -train)

# Combine the scaled_data without train column with the train column
final_df<- cbind(scaled_data, pros$train)

# Split of dataset

train <- scaled_data[pros$train,]
test <- scaled_data[pros$train==FALSE,]

test<-data.frame(test)
train<-data.frame(train)



dim(train)

dim(test)



############################# OLS model ###########################################
# Ols model

model <- lm( lpsa~lcavol+lweight+svi, data = train)

# Print the model summary
summary(model)
#Evaluate model 
data = test[,c("lcavol","lweight","svi")]
predictions <- predict(model, newdata = data)

# Calculate the mean squared error
Rmse_ols <- sqrt(mean((test$lpsa - predictions)^2))

Rmse_ols



################################################Regression Lasso#######################################################
coord_descLasso<-function(X,Y,lambda){
  beta1<-rep(0,dim(X)[2])
  beta2<-rep(1,dim(X)[2])
  while(sqrt(t(beta1-beta2)%*%(beta1-beta2))>0.0000001){
    
    if(dim(X)[2]==2)
    { beta2<-beta1
    for(j in 1:2)
    {Rj<-t(X[,j])%*%(Y-X[,-j]*beta1[-j])
    betaj<-Rj*max(1/(t(X[,j])%*%X[,j])-lambda/(2*abs(Rj)*(t(X[,j])%*%X[,j])),0)
    beta1[j]<-betaj}
    }
    else
    {  beta2<-beta1
    for(j in 1:dim(X)[2])
    {Rj<-t(X[,j])%*%(Y-X[,-j]%*%(beta1[-j]))
    betaj<-Rj*max(1/(t(X[,j])%*%X[,j])-lambda/(2*abs(Rj)*(t(X[,j])%*%X[,j])),0)
    beta1[j]<-betaj}}
    
  }
  return(beta1)
}

lambda_lasso<-seq(4,6,0.01)
rmse_lasso_trials<-c()
min_index<-0
minimum<-100000



for (k in 1:length(lambda_lasso)){
  
  beta_tmp<-coord_descLasso(as.matrix(train[,1:8]),as.matrix(train[,9]),lambda_lasso[k])
  

  
  y_pred = t(beta_tmp)%*%t(as.matrix(test[,1:8]))
  
  rmse <- sqrt(mean((as.matrix(test[,9]) - t(y_pred))^2))
  
  rmse_lasso_trials<-c(rmse_lasso_trials,rmse)
  
  if( min(rmse_lasso_trials) < minimum){
    
    minimum <- min(rmse_lasso_trials)
    
    min_index<-k
  }
  
  
}

# View the mean squared error and trials
print("trials")
print(rmse_lasso_trials)

#best trial
print("best trial ")
sprintf("lowest RMSE %f",minimum)
sprintf("best lambda %f",lambda_lasso[min_index])

df_lasso_result<-data.frame(lambda = lambda_lasso, rmse =rmse_lasso_trials)


ggplot(data = df_lasso_result, aes(x = lambda, y = rmse)) +
  geom_point() +
  xlab("Lambda") +
  ylab("RMSE") +
  ggtitle("Lambda vs. RMSE")

#########Feature selection LASSO###########################################
best_combination= c("lcavol","lweight","svi")

lambda_lasso<-seq(1,30,0.01)
rmse_lasso_trials<-c()
min_index<-0
minimum<-100000



for (k in 1:length(lambda_lasso)){
  
  beta_tmp<-coord_descLasso(as.matrix(train[,best_combination]),as.matrix(train[,9]),lambda_lasso[k])
  
  
  
  y_pred = t(beta_tmp)%*%t(as.matrix(test[,best_combination]))
  
  rmse <- sqrt(mean((as.matrix(test[,9]) - t(y_pred))^2))
  
  rmse_lasso_trials<-c(rmse_lasso_trials,rmse)
  
  if( min(rmse_lasso_trials) < minimum){
    
    minimum <- min(rmse_lasso_trials)
    
    min_index<-k
  }
  
  
}

# View the mean squared error and trials
print("trials")
print(rmse_lasso_trials)

#best trial
print("best trial ")
sprintf("lowest RMSE %f",minimum)
sprintf("best lambda %f",lambda_lasso[min_index])

df_lasso_result<-data.frame(lambda = lambda_lasso, rmse =rmse_lasso_trials)


ggplot(data = df_lasso_result, aes(x = lambda, y = rmse)) +
  geom_point() +
  xlab("Lambda") +
  ylab("RMSE") +
  ggtitle("Lambda vs. RMSE")








################ Elastic net #####################################################################
coord_descElasticnet<-function(X,Y,lambda1,lambda2){
  beta1<-rep(0,dim(X)[2])
  beta2<-rep(1,dim(X)[2])
  while(sqrt(t(beta1-beta2)%*%(beta1-beta2))>0.0000001){
    
    if(dim(X)[2]==2)
    { beta2<-beta1
    for(j in 1:2)
    {Rj<-t(X[,j])%*%(Y-X[,-j]*beta1[-j])
    betaj<-(1/(1+(lambda2/(t(X[,j])%*%X[,j]))))*Rj*max(1/(t(X[,j])%*%X[,j])-lambda1/(2*abs(Rj)*(t(X[,j])%*%X[,j])),0)
    beta1[j]<-betaj}
    }
    else
    {  beta2<-beta1
    for(j in 1:dim(X)[2])
    {Rj<-t(X[,j])%*%(Y-X[,-j]%*%(beta1[-j]))
    betaj<-(1/(1+(lambda2/(t(X[,j])%*%X[,j]))))*Rj*max(1/(t(X[,j])%*%X[,j])-lambda1/(2*abs(Rj)*(t(X[,j])%*%X[,j])),0)
    beta1[j]<-betaj}}
    
  }
  return(beta1)
}

param1<- seq(0,0.06,0.001)
param2<- seq(0,0.01,0.01)
#param1<- seq(0,50,1)
#param2<- seq(0,50,1)

rmse_min_elastic<-1000
lambda1_opti<-NULL
lambda2_opti<-NULL

#lambda1_plot<-c()
#lambda2_plot<-c()
rmse_elastic_trial<-c()

for (i in 1:length(param1)){
  
  for (j in 1:length(param2)){
    #lambda1_plot<-c(lambda1_plot,i)
    #lambda2_plot<-c(lambda2_plot,j)
    
  beta_tmp<-coord_descElasticnet(as.matrix(train[,1:8]),as.matrix(train[,9]),i,j)

  
  y_pred = t(beta_tmp)%*%t(as.matrix(test[,1:8]))
  
  rmse <- sqrt(mean((as.matrix(test[,9]) - t(y_pred))^2))
  
  rmse_elastic_trial<-c(rmse_elastic_trial,rmse)
  print(rmse)
  
  if( rmse_min_elastic > min(rmse_elastic_trial)){
    
    rmse_min_elastic<- rmse
    
    lambda1_opti<-param1[i]
    lambda2_opti<-param2[j]
  }
  
  
}
}
sprintf("lambda1_opti :%f",lambda1_opti)
sprintf("lambda2_opti :%f",lambda2_opti)
sprintf("RMSE_opti :%f",rmse_min_elastic)

#library(scatterplot3d)


# Create a 3D scatterplot
#scatterplot3d(rmse_elastic_trial, lambda1_plot, lambda2_plot, highlight.3d = TRUE,
 #             type = "h", angle = 45, main = "RMSE against Lambda1 and Lambda2",
  #            xlab = "RMSE", ylab = "Lambda1", zlab = "Lambda2")




######################################### Regression ridge #####################################################################

lambda_ridge<-seq(10.4,12,0.02)
rmse_ridge_trials<-c()
min_index<-0
minimum<-100000

coord_ridge<-function(X,Y,lambda){
  solved_mat<-solve((t(X)%*%X+lambda*diag(dim(X)[2])))
  
  beta<-solved_mat%*%(t(X)%*%Y)
  
  return(beta)
}

for (k in 1:length(lambda_ridge)){
  
  beta_tmp<-coord_ridge(as.matrix(train[,1:8]),as.matrix(train[,9]),lambda_ridge[k])
  y_pred = t(beta_tmp)%*%t(as.matrix(test[,1:8]))
  rmse <- sqrt(mean((as.matrix(test[,9]) - t(y_pred))^2))
  
  rmse_ridge_trials<-c(rmse_ridge_trials,rmse)
  
  if( min(rmse_ridge_trials) < minimum){
    
    minimum <- min(rmse_ridge_trials)
    
    min_index<-k
  }
  
  
}





# View the mean squared error and trials
sprintf("trials")
sprintf(rmse_ridge_trials)

#best trial
sprintf("best trial ")
sprintf("lowest RMSE %f",minimum)
sprintf("best lambda %f",lambda_ridge[min_index])

df_ridge<-data.frame(lambda = lambda_ridge, rmse =rmse_ridge_trials)


ggplot(data = df_ridge, aes(x = lambda, y = rmse)) +
  geom_point() +
  xlab("Lambda") +
  ylab("RMSE") +
  ggtitle("Lambda vs. RMSE")


############################################# beta cp ############################################################################
library(Matrix)

ACP_parameters<-matrix(rep(0,64),8,8)

min_acp<-1000
ncomponents<-0


#Matrix to store our PCA estimators
r<-rankMatrix(A)[1]


 
#fonction qui calcule l'estimateur ACP avec k composantes
PCA<- function(X, Y, K){
  X_SVD <-svd(X)
  Vk<-X_SVD$v[,1:K]
  sing_val_k<-diag(X_SVD$d)[1:K,1:K]
  
  inverted_mat=solve(((sing_val_k)%*%(sing_val_k)))
  
  Beta_ACP = Vk%*%inverted_mat%*%(sing_val_k)%*%t(X_SVD$u[,1:K])%*%Y
  
  return (Beta_PCA)
}

for(i in 1:r ){
  ACP_parameters[,i]<-ACP(as.matrix(train[,1:8]),as.matrix(train[,9]), i)
}





rmse_acp<-c()

for (k in 1:r){
  
  beta_tmp<-as.matrix(ACP_parameters[,k])
  
  y_pred = as.matrix(test[,1:8])%*%beta_tmp
  
  
  
  rmse <- sqrt(mean((as.matrix(test[,9]) - (y_pred))^2))
  
  rmse_acp<-c(rmse_acp,rmse)
  if(rmse<min_acp){
    ncomponents<-k
    min_acp<-rmse
    
  }
}

sprintf("lowest acp rmse : %f",min_acp)
sprintf("meilleure combinaison : %d",ncomponents)

## plotting

df_acp_result<-data.frame(components = seq(1:8), rmse =rmse_acp)


ggplot(data = df_acp_result, aes(x = components, y = rmse)) +
  geom_point() +
  xlab("nb_components") +
  ylab("RMSE") +
  ggtitle("nb_components vs. RMSE")


















