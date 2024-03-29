#Aula 12 - Parte 2 - Correla��o serial

library(readxl)
Phillips <- read_excel("R:/HO-231-Rosangela/Dados/Phillips.xlsx")
View(Phillips)

dados = ts(Phillips[,2:3], start=1948, frequency = 1)

layout(1:2)
plot(dados[,1], xlab='anos', ylab='infla��o')
plot(dados[,2], xlab='anos', ylab='desemprego')

install.packages('lmtest')
library(lmtest)
install.packages('dynlm')
library(dynlm)

reg1 = lm (INF ~ UNEM, data = dados)
summary(reg1)

#An�lise gr�fica dos res�duos

res1 = resid(reg1)

layout(1:2)
plot.ts(res1)
n=length(res1)
plot(res1[1:(n-1)], res1[2:n])

#Teste de homocedasticidade

bptest(reg1)

#Teste AR(1) #Verificar se h� autocorrela��o serial

mod = dynlm(res1 ~ L(res1))
mod = lm(res1[2:n] ~ res1[1:(n-1)])
summary(mod)
#Isto �, h� correla��o residual.

#Teste de Durbin-Watson
dwtest(reg1)
#Isto �, rejeita a hip�tese nula de aus�ncia de correla��o serial
cor(res1[1:(n-1)], res1[2:n])

#Modelo de philips 
reg.ea = dynlm(d(INF)~UNEM, data = dados)
summary(reg.ea)

#Taxa natural
tn = reg.ea$coefficients[1]/(-reg.ea$coefficients[2])
tn

residual.ea = resid(reg.ea)
coeftest(dynlm(residual.ea ~ L(residual.ea))) #L cad lag
## p-values = 0.7798 on ne rejet pas lhypothese nulle. 
## travailler avec une serie non stationnaire, resulte de probleme de 
## correlation serial
