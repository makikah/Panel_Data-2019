#Curva de Phillips

dados=ts(Phillips[,2:3],start=1948,freq=1)

layout(1:2)
plot(dados[,1],xlab='anos',ylab='inflação')
plot(dados[,2],xlab='anos',ylab='Desemprego')

install.packages('lmtest')
library(lmtest)

install.packages('dynlm')
library(dynlm)

reg1 = lm(INF~UNEM,data=dados)
summary(reg1)

#análise gráfica
res1=resid(reg1)

layout(1:2)
plot.ts(res1)
n=length(res1)
plot(res1[1:(n-1)],res1[2:n])

#teste de homocedasticidade
bptest(reg1)

#Teste AR(1)
mod=lm(res1[2:n]~res1[1:(n-1)])
summary(mod)

#teste de Durbin-Watson
dwtest(reg1)
cor(res1[1:(n-1)],res1[2:n])

#Modelo de Philpis expectativas aumentadas
reg.ea = dynlm(d(INF)~UNEM,data=dados)
summary(reg.ea)

#Taxa natural
tn=reg.ea$coef[1]/(-reg.ea$coef[2])
tn

residual.ea=resid(reg.ea)
coeftest(dynlm(residual.ea~L(residual.ea)))





























