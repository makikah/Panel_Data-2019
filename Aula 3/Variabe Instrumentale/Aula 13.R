
# Aula 13

## Exemplo Emprego e Salario em porto rico

dados = ts(PRMINWGE, start = 1950, frequency = 1)
trend = ts(1:nrow(dados))
reg = lm(log(PREPOP) ~ log(MINCOV) + log(USGNP) + log(PRGNP) + trend, data = dados)
summary(reg)

# Analise dos residuos
res1 = resid(reg)

n = nrow(dados)
layout(1:2)
plot.ts(res1)
plot(res1[1:(n-1)], res1[2:n])

#Test de homocedasticidade
bptest(reg)

library("car")
coeftest(reg, vcov. = hccm(reg, type ="hc0")) # test robusto:on verifie si la
                                              # variance est constante

library("stargazer")
stargazer(coeftest(reg), coeftest(reg, vcov. = hccm(reg, type = "hc0")),
          digits = 5, column.labels = c("Usual", "Robusto"), type = "text")

# Test de correlacao de Breuch-Godfrey
acf(res1)
pacf(res1)
bgtest(reg, order = 1)

#Estatistica robusto
install.packages("sandwich")
library(sandwich)

stargazer(coeftest(reg), coeftest(reg,vcovHAC), 
          digits = 5, column.labels = c("Usual", "Robusto"), type = "text")

#Correção de COchrane & Orcutt
install.packages("orcutt")
library("orcutt")

coch = cochrane.orcutt(reg)
coch
summary(coch)

#correção de Prais & Winsten
install.packages("prais")
library("prais")

PW = prais_winsten(reg)
summary(PW)
