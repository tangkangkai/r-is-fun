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
  cat(t)
  writeLines("is the t")
  
  if(t < 1){
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
  if(t < 1){
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
  if (t < 1) {
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
  if(t < 1){
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

#INITIAL VALUES
Init_G = 900
Init_I = 30
g_total = 1000

Gmatrix=matrix(0,g_total,1)
Imatrix=matrix(0,g_total,1)
Gmatrix[1,]=Init_G
Imatrix[1,]=Init_I

G_Change_Data = rep(0, g_total)
I_Change_Data = rep(0, g_total)

ABC=matrix(0,g_total,13)
ABC[1,]=c(Init_G,Init_I,G_change(1),I_change(1),0,0,0,0,0,0,0,0,1)

for (t in 2:g_total) {
  g = G(t)
  i = I(t)
  ABC[t,]=c(G(t),I(t),G_change(t),I_change(t),Gin(t),f1(G(t-theta1)),f2(g),f3(g),f4(i),f5(I(t-theta2)),f6(g),f7(g-330),t)
}

ABC[,1]->gt
ABC[,13]->t
plot(t,gt)

t = rep(0,g_total)
for (i in 1:g_total) {
  t[i] = i
}


plot(t, gt, 
     col="blue",
     cex=0.1,
     ylim=c(200,5000),
     xlab = "Time Point",
     ylab = "Glucose Concentration (mg/L)",
     main="The fitting reslut for subject 01-011 on the second day of visit 2 ",
     type = 'l')

lines(t,gt,col="blue",lwd=2)




#expect g to between 900~2500