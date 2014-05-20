#Sara Cannon, Homework 6

#Problem 1

library(deSolve)

LV_Pred_RM<-function(t,y,params){
  V<-y[1]
  P<-y[2]
  with(as.list(params),{
    dVdt<-b*V*(1-(1/K1)*V)-a*V*P/(1+a*h*V)
    dPdt<-e*a*V*P/(1+a*h*V)-d*P
    return(list(c(dVdt,dPdt)))
  })}

a = 0.0001;      # attack rate
e = 0.1;       # conversion efficiency rate
d = 0.1;       # predator death rate
b = 0.5;       # prey birth rate
h = 0.0015;       # prey handling time
K1 = 15000;       #Prey carrying capacity (Initial scenario)
K2 = 50000;       #Prey carrying capacity (2nd scenario)
K3 = 100000;     #Prey carrying capacity (3rd scenario)
K4 = 1000000;    #Prey carrying capacity (4th scenario)
K5 = 10000000;   #Prey carrying capacity (5th scenario)
y0 = c(V0=8000,P0=2000);  # Initial values for V and P populations

#Evaluate scenario for K initial

params<-c(b=b,a=a,e=e,d=d,K1=K1,h=h)
tspan<-seq(0,300,by=0.1)  # Timespan to evaluate (in months)

isoV1<-function(x){(b+b*x*(a*h-(1/K1)-a*(1/K1)*h*x))/a}
isoP<-d/(a*e-a*d*h)

out<-ode(y0,tspan,LV_Pred_RM,params)
out<-data.frame(out)

par(mfrow=c(1,2))
plot(out[,1],out[,2],type='l',col='green',xlab='Month',ylab='Popn size, V and P',ylim=c(0,12000))
lines(out[,1],out[,3],col='blue')

plot(out[,2],out[,3],type='l',col='red',xlab='V, Prey abundance',ylab='P, Predator abundance',
     ylim=c(0,2000),xlim=c(0,12000))
segments(isoP, 0, isoP, 100,lwd=2) #what is the 100 referring to here?
curve(isoV1, 0, 300, add = T) #what is the 300 referring to here?

#It looks like Graph 1 (Pop size, V&P vs Month) looks okay, but prey abundance vs. predator abundance is off somehow.
