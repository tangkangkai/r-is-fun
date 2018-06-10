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

t = rep(0,600)
m_fit = rep(0, 600)
for (i in 1:600) {
  t[i] = i
  m_fit[i] = Gin(i)
}

plot(t, m_fit)