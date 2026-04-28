rm(list = ls())
OverdoseDeaths = read.csv("~/Box/P-IG/OverdoseDeaths.csv", header=T, stringsAsFactors = F, colClasses = c("character", "character", "integer", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
OverdoseDeaths = OverdoseDeaths[,2:10]
rates = OverdoseDeaths[,4:9]/apply(OverdoseDeaths[,4:9],1,sum)

df = data.frame(value = c(rates[,1], rates[,2], rates[,3], rates[,4], rates[,5], rates[,6]), Drug = rep(c("Cocaine", "Synthetic", "Heroin", "Methadone", "NatOpioids", "PsyStim"), each = 76) )
ggplot(df, aes(x = value, fill = Drug)) + 
  geom_histogram()+
  xlab("Proportion of overdose deaths")


