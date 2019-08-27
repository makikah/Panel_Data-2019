#***************************************************************
#  CLASS 1 - MULTIPLE REGRESSION - REVIEW
#  AUTHOR:  ALEXANDRE GORI MAIA
#           UNIVERSITY OF CAMPINAS
#           CENTER FOR AGRICULTURAL AND ENVIRONMENTAL ECONOMICS
#***************************************************************

# import dataset with information for CO2, GDP and % manufacturing 
countries <- read.csv("Z:/HO-235/Data/Data_CO2.csv") 



# linear regression by ols
ols <- lm(co2 ~ gdp + ind, data=countries)
summary(ols)

