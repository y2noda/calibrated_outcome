library(simex)

x <- rnorm(200,0,100)
u <- rnorm(200,0,25)
w <- x+u
y <- x +rnorm(200,0,9)
true.model <- lm(y~x) # True model
naive.model <- lm(y~w, x=TRUE)
simex.model <- simex(model = naive.model
                     , SIMEXvariable = "w"
                     , measurement.error= 25)
plot(x,y)
abline(true.model,col="darkblue")
abline(simex.model,col ="red")
abline(naive.model,col = "green")
legend(min(x),max(y),legend=c("True Model","SIMEX model","Naive Model")
       , col = c("darkblue","red","green"),lty=1)

plot(simex.model, mfrow = c(2,2))
