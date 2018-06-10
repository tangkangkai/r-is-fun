#FUNCTIONS
Gm = function(t) {
  b = 80
  k = 5100
  p = -(t ^ 2) / (2 * b ^ 2)
  return(k * t / (b * b) * exp(p))
}

Gin_At_m = function(t, tm)
{
  u = 0
  if (t > tm)
  {
    u = 1
  }
  return(u * Gm(t - tm))
}

Gin = function(t) {
  
  first_stage = Gin_At_m(t, 50)
  second_stage = Gin_At_m(t, 150)
  third_stage = Gin_At_m(t, 900)
  
  return(first_stage + second_stage + third_stage)
}

G=function(t) {
  if(t <= 0){
    return (Init_G)
  }
  else if(Gmatrix[t]!=0){
    return (Gmatrix[t])
  } 
  else {
    res=G(t-1)+G_change(t-1)
    Gmatrix[t]=res
    assign('Gmatrix',Gmatrix,envir=.GlobalEnv)
    return (res)
  }
}

G_change=function(t) {
  if(t <=0){
    return (0)
  }
  else if (G_Change_Data[t] != 0)
  {
    return (G_Change_Data[t])
  }
  else {
    g = G(t)
    i = I(t)
    data = Gin(t)+f5(I(t-theta2))*f6(g)-(f2(g)+beta*f3(g)*f4(i))
    G_Change_Data[t] = data
    assign('G_Change_Data',G_Change_Data,envir=.GlobalEnv)
    return (data)
  }
}

I=function(t) {
  if (t <=0) {
    return (Init_I)
  }
  else if (Imatrix[t] != 0) {
    return (Imatrix[t])
  } 
  else {
    res = I(t-1) + I_change(t-1)
    Imatrix[t] = res
    assign('Imatrix',Imatrix,envir=.GlobalEnv)
    return (res)
  }
}

I_change= function(t) {
  if(t <= 0){
    return(0)
  }
  else if (I_Change_Data[t] != 0)
  {
    return (I_Change_Data[t])
  }
  else {
    data = alpha*f1(G(t-theta1))-d*I(t)
    I_Change_Data[t] = data
    assign('I_Change_Data',I_Change_Data,envir=.GlobalEnv)
    return (data)
  }
}

f1=function(g) {
  return(210/(1+exp((2000-g/10)/300)))
}
f2=function(g) {
  return(72*(1-exp(-g/144/10)))
}
f3=function(g) {
  return(g/10000)
}
f4=function(i) {
  return(40+(940-40)/(1+exp(-1.77*log(i/80*(1/11+1/0.2/100), exp(1)))))
}
f5=function(i) {
  return(180/(1+exp(0.29*(i/3))))
}
f6=function(g) {
  return(1/(1+exp(5*(g/10000-2))))
}
f7=function(g) {
  return(20+(140-20)/(1+exp(-2.4*(g/10000-2))))
}

alpha=0.896
beta=0.818
theta1=5.02
theta2=15.2
d=0.107

g_real<- c(1044,1044,1044,1026,1026,1026,1026,1044,
           1080,1134,1170,1206,1242,1260,1278,1260,1224,1170,
           1116,1062,1026,990,972,990,1008,1008,1008,1008,
           1008,990,972,972,954,954,954,954,954,972,
           972,972,972,954,936,918,918,918,918,918,
           936,954,990,1044,1098,1134,1152,1170,1188,1170,
           1170,1134,1116,1080,1044,1026,1008,1008,1008,1026,
           1044,1062,1098,1116,1134,1134,1134,1134,1116,1098,
           1062,1044,1026,1008,1008,1026,1044,1062,1116,1134,
           1152,1152,1152,1134,1116,1098,1062,1062,1044,1026,
           1008,972,972,954,954,954,972,972,990,1008,1026,1026,1026,
           1026,1026,1026,1026,1026,1026,1026,1026,1008,1008,1008,
           1008,1008,990,990,972,972,972,972,954,954,954,954,954,
           954,954,972,972,1008,1026,1044,1044,1062,1062,1080,1098,1098,
           1116,1116,1116,1098,1080,1044,1026,1026,1008,1008,1026,
           1026,1026,1044,1044,1044,1044,1044,1026,1026,1026,1026,
           1026,1026,1026,1026,1026,1026,1026,1026,
           1026,1026,1044,1044,1044,1044,1044,1062,
           1062,1062,1062,1044,1044,1044,1026,1026,1026,1008,
           1008,1008,1008,1008,1008,1008,1008,1008,1008,990,
           990,972,990,990,990,1008,1008,1008,1008,1008,
           1026,1026,1026,1026,1044,1044,1044,1044,1026,1026,
           1026,1026,1026,1026,1026,1026,1026,1026,1026,1026,
           1026,1026,1026,1026,1026,1026,1026,1026,1026,1026,
           1026,1026,1026,1026,1026,1026,1026,1026,1026,1026,
           1026,1026,1026,1026,1026,1026,1026,1026,1026,1026,
           1026,1026,1026,1008,1008,990,972,972,972,972,990,1008,1026,1026,1026,1044,1044,1044,1044,1044,
           1044,1044)


#INITIAL VALUES
Init_G = 90
Init_I = 10
g_total = length(g_real)
t_total = 5 * g_total

Gmatrix=matrix(0,t_total,1)
Imatrix=matrix(0,t_total,1)
Gmatrix[1,]=Init_G
Imatrix[1,]=Init_I

G_Change_Data = rep(0, t_total)
I_Change_Data = rep(0, t_total)

ABC=matrix(0,t_total,13)
ABC[1,]=c(Init_G,Init_I,G_change(1),I_change(1),0,0,0,0,0,0,0,0,1)

for (t in 2:t_total) {
  g = G(t)
  i = I(t)
  ABC[t,]=c(G(t),I(t),G_change(t),I_change(t),Gin(t),f1(G(t-theta1)),f2(g),f3(g),f4(i),f5(I(t-theta2)),f6(g),f7(g-330),t)
}

ABC[,1]->gt
ABC[,13]->t
plot(t,gt)

g_fit = rep(0, g_total)
for (t in 1:t_total) {
  g_fit_at_t = G(t)
  if (t %% 5 == 1) {
    g_fit[floor(t/5) + 1] = g_fit_at_t
  }
}

t = rep(0,g_total)
for (i in 1:g_total) {
  t[i] = i
}


plot(t, g_fit, 
     col="blue",
     cex=0.1,
     ylim=c(0,300),
     xlab = "Time Point",
     ylab = "Glucose Concentration (mg/L)",
     main="The fitting reslut for subject 01-011 on the second day of visit 2 ",
     type = 'l')

lines(t,g_fit,col="blue",lwd=2)




#expect g to between 900~2500