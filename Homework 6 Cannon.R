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

title("Population projections for V and P",outer=TRUE)


#Plot the P-V isoclines for each value of K
#Need to work on these graphs
#They all show the population projection, but not the predator or prey isoclines.

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

title("P-V Isoclines", outer=T)

print("At the lowest K, the prey population increases and the predator population decreases slightly, before they stabilize.")
print("As K increases, oscillations in the populations of prey and predator increase, and the time it takes to reach equilibrium also increases.")
print("When K gets too high, the system appears to go into chaos.")


#Problem #2

#Evaluate Jacobian matrix for K = 15,000
V<-d/(a*e-a*d*h)
P1<-b/a+(b*V*(a*h-1/K1-a/K1*h*V)/a)
a11<-b-(P1*a)/(V*a*h+1)^2-(2*V*b)/K1
a12<-(V*a)/(V*a*h+1)
a21<-(P1*a*e)/(V*a*h+1)^2
a22<-(V*a*e)/(V*a*h+1)-d
  
matrx1<-function(a11,a12,a21,a22){
  J1<-matrix(c(a11,    a12,
              a21,    a22),byrow=T,ncol=2)
}

J1<-matrx1(a11,a12,a21,a22)
print(J1)
eigen(J1)

print("Because the diagnol of the Jacobian matrix is not equal to zero, this indicates there is density dependence.")
print("The eigenvalues are all real numbers, so indicates there are no oscillations.")
print("Because some of the real parts are <0, this indicates there is an attractor-repellor.")


#Evaluate Jacobian matrix for K = 50,000

P2<-b/a+(b*V*(a*h-1/K2-a/K2*h*V)/a)
b11<-b-(P2*a)/(V*a*h+1)^2-(2*V*b)/K2
b12<-(V*a)/(V*a*h+1)
b21<-(P2*a*e)/(V*a*h+1)^2
b22<-(V*a*e)/(V*a*h+1)-d

matrx2<-function(b11,b12,b21,b22){
  J2<-matrix(c(b11,    b12,
               b21,    b22),byrow=T,ncol=2)
}

J2<-matrx2(b11,b12,b21,b22)

print(J2)
eigen(J2)
print("The non-zero values in the diagonal of the Jacobian matrix show that there is density dependence.")
print("Because some of the values of the eigenvalues are <0, this indicates that there is an attractor-repellor.")
print("Because all of the eigenvalues are real numbers (imaginary parts are absent), there are no oscillations.")

#Evaluate Jacobian matrix for K = 100,000

P3<-b/a+(b*V*(a*h-1/K3-a/K3*h*V)/a)
c11<-b-(P3*a)/(V*a*h+1)^2-(2*V*b)/K3
c12<-(V*a)/(V*a*h+1)
c21<-(P3*a*e)/(V*a*h+1)^2
c22<-(V*a*e)/(V*a*h+1)-d

matrx3<-function(c11,c12,c21,c22){
  J3<-matrix(c(c11,    c12,
               c21,    c22),byrow=T,ncol=2)
}

J3<-matrx3(c11,c12,c21,c22)

print(J3)
eigen(J3)
print("Because the Jacobian matrix has non-zero numbers in the diagonal, this shows density dependence.")
print("The eigenvalues are all real numbers, indicating there are no oscillations.")
print("Because some of the eigenvalues are <0, this indicates there is an attractor-repellor.")

#Evaluate Jacobian matrix for K = 1,000,000

P4<-b/a+(b*V*(a*h-1/K4-a/K4*h*V)/a)
d11<-b-(P4*a)/(V*a*h+1)^2-(2*V*b)/K4
d12<-(V*a)/(V*a*h+1)
d21<-(P4*a*e)/(V*a*h+1)^2
d22<-(V*a*e)/(V*a*h+1)-d

matrx4<-function(d11,d12,d21,d22){
  J4<-matrix(c(d11,    d12,
               d21,    d22),byrow=T,ncol=2)
}

J4<-matrx4(d11,d12,d21,d22)

print(J4)
eigen(J4)
print("Because the Jacobian matrix has non-zero numbers in the diagonal, this shows density dependence.")
print("The eigenvalues are all real numbers, indicating there are no oscillations.")
print("Because some of the eigenvalues are <0, this indicates there is an attractor-repellor.")


#Evaluate Jacobian matrix for K = 10,000,000

P5<-b/a+(b*V*(a*h-1/K4-a/K4*h*V)/a)
e11<-b-(P5*a)/(V*a*h+1)^2-(2*V*b)/K5
e12<-(V*a)/(V*a*h+1)
e21<-(P5*a*e)/(V*a*h+1)^2
e22<-(V*a*e)/(V*a*h+1)-d

matrx5<-function(e11,e12,e21,e22){
  J5<-matrix(c(e11,    e12,
               e21,    e22),byrow=T,ncol=2)
}

J5<-matrx5(e11,e12,e21,e22)

print(J5)
eigen(J5)

print("Because the Jacobian matrix has non-zero numbers in the diagonal, this shows density dependence.")
print("The eigenvalues are all real numbers, indicating there are no oscillations.")
print("Because some of the eigenvalues are <0, this indicates there is an attractor-repellor.")
