#Aula 12 - Parte 2 - Correlação serial

library(readxl)
Phillips <- read_excel("R:/HO-231-Rosangela/Dados/Phillips.xlsx")
View(Phillips)

dados = ts(Phillips[,2:3], start=1948, frequency = 1)

layout(1:2)
plot(dados[,1], xlab='anos', ylab='inflação')
plot(dados[,2], xlab='anos', ylab='desemprego')

install.packages('lmtest')
library(lmtest)
install.packages('dynlm')
library(dynlm)

reg1 = lm (INF ~ UNEM, data = dados)
summary(reg1)

#Análise gráfica dos resíduos

res1 = resid(reg1)

layout(1:2)
plot.ts(res1)
n=length(res1)
plot(res1[1:(n-1)], res1[2:n])

#Teste de homocedasticidade

bptest(reg1)

#Teste AR(1) #Verificar se há autocorrelação serial

mod = dynlm(res1 ~ L(res1))
mod = lm(res1[2:n] ~ res1[1:(n-1)])
summary(mod)
#Isto é, há correlação residual.

#Teste de Durbin-Watson
dwtest(reg1)
#Isto é, rejeita a hipótese nula de ausência de correlação serial
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
