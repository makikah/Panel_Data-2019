#Modelo Salário e Produtividade

dados=ts(earns[,2:3],start=1947,freq=1)
View(dados)

layout(1:2)
plot(dados[,1],xlab='Ano',ylab='Produtividade')
plot(dados[,2],xlab='Ano',ylab='Salário')

outphr=log(dados[,1])
hrwage=log(dados[,2])

n=length(outphr)
trend=ts(1:n)

modelo1=lm(hrwage ~outphr + trend)
summary(modelo1)

aux1=lm(hrwage ~trend)
summary(aux1)

aux2=lm(outphr ~trend)
summary(aux2)

plot.ts(aux1$residuals)
plot.ts(aux2$residuals)

modelo2=lm(aux1$residuals ~aux2$residuals)
summary(modelo2)

cor(modelo2$residuals[1:(n-1)],modelo2$residuals[2:n])
layout(1:1)
acf(modelo2$residuals)

modelo3=lm(diff(hrwage)~diff(outphr))
summary(modelo3)


cor(modelo3$residuals[1:(n-2)],modelo3$residuals[2:(n-1)])
layout(1:1)
acf(modelo3$residuals)


#Análise gráfica dos resíduos
layout(1:2)
plot.ts(modelo2$residuals)
plot.ts(modelo3$residuals)

plot(modelo2$residuals[1:(n-1)],modelo2$residuals[2:n])
plot(modelo3$residuals[1:(n-2)],modelo3$residuals[2:(n-1)])














