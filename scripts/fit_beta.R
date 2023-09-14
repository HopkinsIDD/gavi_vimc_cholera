# Fit Kirsten's predictive distribution to beta distribution 

library(MASS)
library(dplyr)
library(ggplot2)

# inputpath <- "/Users/hammer/Downloads/Cholera/VIMC_project/params"
# load(paste0(inputpath, "/prop_pos.RData"))

# load posterior draws of cholera positivity from Kirsten's paper
load(paste0("/prop_pos.RData"))


# pre-calculation for the start values
calc_start <- function(x){
    
  m_x <- mean(x, na.rm = TRUE)
  s_x <- sd(x, na.rm = TRUE)
    
  alpha <- m_x*((m_x*(1 - m_x)/s_x^2) - 1)
  beta <- (1 - m_x)*((m_x*(1 - m_x)/s_x^2) - 1)
    
  return(list(alpha = alpha, beta = beta))
    
}
  
start <- calc_start(prop_pos)
shape1 <- start$alpha 
shape2 <- start$beta
  
fit <- fitdistr(prop_pos, densfun = "beta", start = list(shape1 = shape1, shape2 = shape2))
  
df_param <- cbind(as.data.frame(fit[1]), as.data.frame(fit[2]))
  
df_param_fn <- paste0("beta_params.csv")
write.csv(df_param, df_param_fn)
  
#}


# plot the fitted distribution over the original data
dis <- as.data.frame(prop_pos)
dis %>%
  ggplot(aes(x = prop_pos)) + geom_density(alpha = 0.3)

x <- seq(0,1,length=1000)
db <- dbeta(x, df_param$estimate[1], df_param$estimate[2])
# estimation: shape1 = 1.562326; shape2 = 1.638044 (2022)
# estimation after Kirsten's update on her results: shape1 = 5.6048124; shape2 =  5.1251011


ggplot() + 
  geom_histogram(data = dis, aes(x = prop_pos, y = ..density..),
                 colour = 1, fill = "lightgrey") +
  geom_density() +
  geom_line(aes(x,db)) +
  xlab("Predictive proportion positive") + 
  ylab("Density")


# make the new parameter file as model input
parameters <- data.frame(
  shape1 = df_param$estimate[1],
  shape2 = df_param$estimate[2]
)

write.csv(parameters, "parameters.csv")






