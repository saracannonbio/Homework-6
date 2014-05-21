#Sara Cannon, Homework 6

#Problem 1

rm=(list())

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

#Evaluate scenario for K2

LV_Pred_RM2<-function(t,y,params2){
  V<-y[1]
  P<-y[2]
  with(as.list(params2),{
    dVdt2<-b*V*(1-(1/K2)*V)-a*V*P/(1+a*h*V)
    dPdt2<-e*a*V*P/(1+a*h*V)-d*P
    return(list(c(dVdt2,dPdt2)))
  })}

params2<-c(b=b,a=a,e=e,d=d,K2=K2,h=h)

isoV2<-function(x){(b+b*x*(a*h-(1-K2)-a*(1/K2)*h*x))/a}
out2<-ode(y0,tspan,LV_Pred_RM2,params2)
out2<-data.frame(out2)

#Evaluate scenario for K3

LV_Pred_RM3<-function(t,y,params3){
  V<-y[1]
  P<-y[2]
  with(as.list(params3),{
    dVdt3<-b*V*(1-(1/K3)*V)-a*V*P/(1+a*h*V)
    dPdt3<-e*a*V*P/(1+a*h*V)-d*P
    return(list(c(dVdt3,dPdt3)))
  })}

params3<-c(b=b,a=a,e=e,d=d,K3=K3,h=h)
isoV3<-function(x){(b+b*x*(a*h-(1-K3)-a*(1/K3)*h*x))/a}
out3<-ode(y0,tspan,LV_Pred_RM3,params3)
out3<-data.frame(out3)

#Evaluate scenario for K4

LV_Pred_RM4<-function(t,y,params4){
  V<-y[1]
  P<-y[2]
  with(as.list(params4),{
    dVdt4<-b*V*(1-(1/K4)*V)-a*V*P/(1+a*h*V)
    dPdt4<-e*a*V*P/(1+a*h*V)-d*P
    return(list(c(dVdt4,dPdt4)))
  })}

params4<-c(b=b,a=a,e=e,d=d,K4=K4,h=h)
isoV4<-function(x){(b+b*x*(a*h-(1-K4)-a*(1/K4)*h*x))/a}
out4<-ode(y0,tspan,LV_Pred_RM4,params4)
out4<-data.frame(out4)


#Evaluate scenario for K5

LV_Pred_RM5<-function(t,y,params5){
  V<-y[1]
  P<-y[2]
  with(as.list(params5),{
    dVdt5<-b*V*(1-(1/K5)*V)-a*V*P/(1+a*h*V)
    dPdt5<-e*a*V*P/(1+a*h*V)-d*P
    return(list(c(dVdt5,dPdt5)))
  })}

params5<-c(b=b,a=a,e=e,d=d,K5=K5,h=h)
isoV5<-function(x){(b+b*x*(a*h-(1-K5)-a*(1/K5)*h*x))/a}
out5<-ode(y0,tspan,LV_Pred_RM5,params5)
out5<-data.frame(out5)

#Plot projected population trajectory for each value of K

par(mfrow=c(2,3)) #show all plots on same screen

plot(out[,1],out[,2],type='l',col='green',xlab='Month',ylab='Popn size, V and P',ylim=c(0,12000))
lines(out[,1],out[,3],col='blue')
title(main="K = 15,000")

plot(out2[,1],out2[,2],type='l',col='green',xlab='Month',ylab='Popn size, V and P',ylim=c(0,12000))
lines(out2[,1],out2[,3],col='blue')
title(main="K = 50,000")

plot(out3[,1],out3[,2],type='l',col='green',xlab='Month',ylab='Popn size, V and P',ylim=c(0,12000))
lines(out3[,1],out3[,3],col='blue')
title(main="K = 100,000")

plot(out4[,1],out4[,2],type='l',col='green',xlab='Month',ylab='Popn size, V and P',ylim=c(0,12000))
lines(out4[,1],out4[,3],col='blue')
title(main="K = 1,000,000")

plot(out5[,1],out5[,2],type='l',col='green',xlab='Month',ylab='Popn size, V and P',ylim=c(0,12000))
lines(out5[,1],out5[,3],col='blue')
title(main="K = 10,000,000")

#Plot the P-V isoclines for each value of K
#Need to work on these graphs

par(mfrow=c(2,3)) #show all plots on same screen

plot(out[,2],out[,3],type='l',col='red',xlab='V, Prey abundance',ylab='P, Predator abundance',
     ylim=c(0,5000),xlim=c(0,15000))
segments(isoP, 0, isoP, 100,lwd=2) 
curve(isoV1, 0, 30000, add = T)
title(main="K = 15,000") 

plot(out2[,2],out2[,3],type='l',col='red',xlab='V, Prey abundance',ylab='P, Predator abundance',
     ylim=c(0,5000),xlim=c(0,20000))
segments(isoP, 0, isoP, 100,lwd=2)
curve(isoV2, 0, 30000, add = T)
title(main="K = 50,000")

plot(out3[,2],out3[,3],type='l',col='red',xlab='V, Prey abundance',ylab='P, Predator abundance',
     ylim=c(0,8000),xlim=c(0,30000))
segments(isoP, 0, isoP, 100,lwd=2)
curve(isoV3, 0, 30000, add = T)
title(main="K = 100,000")

plot(out4[,2],out4[,3],type='l',col='red',xlab='V, Prey abundance',ylab='P, Predator abundance',
     ylim=c(0,10000),xlim=c(0,50000))
segments(isoP, 0, isoP, 100,lwd=2)
curve(isoV4, 0, 30000, add = T)
title(main="K = 1,000,000")

plot(out5[,2],out5[,3],type='l',col='red',xlab='V, Prey abundance',ylab='P, Predator abundance',
     ylim=c(0,10000),xlim=c(0,50000))
segments(isoP, 0, isoP, 100,lwd=2)
curve(isoV5, 0, 30000, add = T)
title(main="K = 10,000,000")