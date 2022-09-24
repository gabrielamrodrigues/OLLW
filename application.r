rm(list=ls(all=TRUE))

##'Packages

library(ggplot2)
library(readxl)
library(gamlss)
library(gamlss.cens)


##'Data
covid <- read.delim('https://raw.githubusercontent.com/gabrielamrodrigues/OLLW/main/data.covid.txt')

#Gamlss odd log logistic Weibull model
source('https://raw.githubusercontent.com/gabrielamrodrigues/OLLW/main/ollw.gamlss.r')

##'Marginal Model

fit1 <- gamlss(Surv(times,censor)~ 1,family=cens(WEI)(sigma.link="log"),data=covid,c.crit=0.001, n.cyc=200)

fit2 <-gamlss(Surv(times,censor)~ 1,family=cens(OLLW)(),data=covid, c.crit=0.001,sigma.start = fit1$sigma.fv, n.cyc=2000, mu.start = fit1$mu.fv)

BIC(fit1,fit2)
c(AIC(fit1),AIC(fit2))
LR.test(fit1,fit2)

km1 <- survfit(Surv(times,censor) ~1, data=covid)
plot(km1,conf.int=F,mark.time=T,xlab="",ylab="S(t)",lwd=3,lty=1,col="black")
curve(1-pWEI(x,mu=fit1$mu.fv[1],sigma=fit1$sigma.fv[1]),add=T,type="l",col='blue',lwd=2)
curve(1-pOLLW(x,mu=fit2$mu.fv[1],sigma=fit2$sigma.fv[1],nu=fit2$nu.fv[1]),add=T,type="l",col='red',lwd=2)
legend("topright",c("OLLW","Weibull"),bty='n',cex=1,
       col=c('red',"blue"), lty=c(1,1), lwd=c(2,2))


##'Complete model
 
fit2r <-gamlss(Surv(times,censor)~ sex+age+heart+asthma+diab+neuro+obesity, sigma.formula=~sex+age+heart+asthma+diab+neuro+obesity,family=cens(OLLW)(),data=covid, c.crit=0.001,sigma.start = fit1$sigma.fv, n.cyc=2000, mu.start = fit1$mu.fv)

summary(fit2r,type = 'qr')
BIC(fit2r)
AIC(fit2r)

##'Selection of covariates

m0 <- gamlss(Surv(times,censor)~  1,family=cens(OLLW),data=covid,
             c.crit=0.001,sigma.start = fit1$sigma.fv,
             n.cyc=2000, mu.start = fit1$mu.fv)

m11 <- stepGAICAll.A(m0, scope=list(lower=~1,upper=~sex+age+heart+asthma+diab+neuro+obesity), nu.try = FALSE)

#Final model 

model <- gamlss(formula = Surv(times,censor) ~ age +   
                  diab+neuro, sigma.formula = ~age + neuro + obesity,  
                family = cens(OLLW), data = covid, mu.start = fit1$mu.fv,  
                sigma.start = fit1$sigma.fv, c.crit = 0.001, n.cyc = 2000, trace = TRUE) 

summary(model, type='qr')


AIC(model)
BIC(model)


##'Residuals analysis

##'Modified deviance residual

fit_final <- model
mu<-     fit_final$mu.fv
sigma<-   fit_final$sigma.fv    
nu <- fit_final$nu.fv
n <- fit_final$N
rm= covid$censor+log(1-pOLLW(covid$times,mu=mu, sigma=sigma, nu=nu))
rd = sign(rm)*(-2*(rm+(covid$censor*log(covid$censor-rm))))^(0.5)
plot(rd, ylab = "Modified deviance residual", pch =20 , ylim =c(-4 ,4))
abline(h=c(-3,3))


index <- seq(1:length(rd))
d1r <- data.frame(index,rd)
ggplot(d1r, aes(x=index, y=rd)) + ylim(-4,4)+
  geom_point()+ylab('Modified deviance residual')+
  xlab('Index')+
  geom_hline(yintercept = 3,linetype = "dashed",size=0.8)+
  geom_hline(yintercept = -3,linetype = "dashed",size=0.8)+
  geom_hline(yintercept = 0,color='darkblue')+
  theme(legend.position="none",panel.background = element_rect(fill = "white", colour = "grey50"),axis.title = element_text( size = (12)),axis.text = element_text(size=11))
