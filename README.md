---
title: "Measuring father involvement across Europe with Bayesian hierarchical models"
# author: "Brett Ory"
thumbnailImagePosition: left
thumbnailImage: https://images.fun.com/products/26675/2-1-75976/darth-vader-1-dad-mens-t-shirt1.jpg
coverImage: https://images.fun.com/products/26675/2-1-75976/darth-vader-1-dad-mens-t-shirt1.jpg
metaAlignment: center
coverMeta: out
date: 2018-01-21T21:13:14-05:00
categories: ["Father involvement"]
tags: ["Bayes", "multilevel analysis", "jags", "R2Jags"]
---


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
devtools::install_github("hadley/emo")
library(dplyr)
library(Hmisc)
```


In ** breaking ** news surprising no one except sociologists: European countries are not all that dissimilar! 



Forgive me for starting with the conclusion. It took me four years to actually reach one, so I'm a little attached to it. For this week I will be showing some code from the first chapter of my dissertation on cross-national variation in the association between mother's work hours and father's share of childcare. Don't worry, I'll explain what that is as briefly and jargon-freely as I can.  

<br>

## Father involvement  across Europe

My dissertation asks what factors in men's families, social class, and country make it possible for them to spend time with their children. We hear a lot about how moms are systematically disadvantaged because they have no maternity leave, for example, but we hear less about what impact that has on dads. My dissertation consists of four related studies. The one I'm going to discuss today is titled:

_Mother's work hours and father's share of involvement in cross-national perspective_

Since the 1960s, women have been gradually increasing the amount of time they spend in paid employment and men have been increasing the amount of time they spend caring for children. Because these two trends occurred more or less together, people often assume that they are related. For example, that in families where mothers work more hours, men will be more involved in childcare. This assumption is even enshrined in European policy on parental leave. One [EU directive](http://eur-lex.europa.eu/legal-content/EN/TXT/?uri=celex%3A32010L0018) stipulates that governments should encourage fathers to take leave in order to help increase female employment. 

My research asks three things:   
1. whether this assumption is true,    
2. is it true for all<sup><a href="#fn1" id="ref1">1</a></sup> European countries and Australia, and     
3. if there are cross-national differences in the strength of the relationship between mothers' work hours and father involvement, can they be explained by paternity leave policy? 


In my paper I consider other measures as well, but my results are the same for all of them so I will use paternity leave as an example here.  

<br> 

## The data

I test my data using the first wave of the [Generations and Gender Survey](http://www.ggp-i.org/), downloaded in 2015 (yes, three years ago). Data cleaning for this project was brutal, and way too long to describe in this post. I will put it on GitHub soon. The data is open to the public, but you have to request access from the data collectors, so unfortunately I can't post it online. 

```{r load data, include=FALSE}
load("ggs.mf_011216.RData")
ggs.resp <- subset(ggs.resp, select=c(country, fistr, hrwkMom.c, pat.c,
                               maleage.c, maleedu.c, femedu.c, numkid, 
                               anyunder4, female))
```


Here is a preview of the data
```{r data table, warning = FALSE}
library(knitr)

Varname <- names(ggs.resp)
Description <- c("Country","Share of father involvement, 0 = always the mother, 4 = always the father", "Hours the mother works, mean centered", "Effective paternity leave, mean centered", "Age of the father, mean centered", "Education of the father, mean centered", "Education of the mother, mean centered", "Number of children living at home", "Any children under 4", "The respondent is female, 1 = female, 0 = male")
Range <- round(data.frame(Mean=sapply(ggs.resp[,-c(1)],mean,na.rm=T),Min=sapply(ggs.resp[,-c(1)],min,na.rm=T),Max=sapply(ggs.resp[,-c(1)],max,na.rm=T)),2) # I remove the first row (country) because I want to round the output to the hundreth, and you can't round strings
Range[nrow(Range) + 1,] = list(NA,"Australia","Russia") # now re-adding country value
Range <- Range[c(10,1:9),] # put in order again
rownames(Range) <- c()

table <- as.data.frame(cbind(Varname,Description, Range))
kable(table[1:10,])
```

```{r, echo=F}
# clean up global environment
rm(Range, table, Description, Varname)
```

The data consists of 10 variables. Continuous variables are mean centered, meaning the survey mean is subtracted from each respondent's value, such that a negative number means they are below the mean, a zero means they are average, and a positive number means they are greater than the mean. I did this to help with interpretation of the final results--it shouldn't have an effect on the estimation of regression parameters. Effective paternity leave is the number of days available to fathers nationwide, multiplied by the percentage they get paid. 5 days of paternity leave paid at 50% would be the same effect amount of leave as 2.5 days at 100%. 

<br> 

## The method

I use Bayesian multilevel analysis in R2Jags to model a) fixed effects for individual level variables and the effect of paternity leave; b) random effects for the effect of mothers' work hours, and c) cross-level interaction between random effects and paternity leave. Mathematically, this can be expressed as:

Level 1 (Fixed effects): 

$$Fatherinvovlement_{i,j} =  \alpha_{0,j} + \beta_{1,j}*mother work hours_{i,j} + \beta_{2,j}X_{i,j} + \sigma_{i,j}$$

Level 2 (Random effects):

$$\alpha_{0,j} = \gamma_{0,0} + \gamma_{0,j}*paternity leave_{j} + \sigma_{0,j}$$

Level 2 (Cross level interaction):

$$\beta_{1,j} = \gamma_{1,0} + \gamma_{1,1}*paternity leave_{j} + \sigma_{1,j}$$


I'm not very good at mathematical notation, so I've written out the names of the variables of interest to make this more easily readable. In short the first equation represents the individual-level regression, which has the parameters mother's work hours and matrix X. X stands for the other individual-level variables listed above which I include as control variables. If we were not accounting for the multilevel structure this would be the total equation. Both the constant as denoted by $\alpha_{0,j}$ and the slope as denoted by $\beta_{1,j}$ are allowed to vary across countries (that's what the subscript $_j$ indicates) 

<br>

## Jags

Ok, so now we know what we want our equation to look like, how do we program this in jags (and what is jags)? JAGS stands for Just Another Gibbs Sampler, and it's basically the mac version of WinBUGS from the perspective of the user (i.e. the syntax is mostly the same). If you've ever tried to run WinBUGS via a windows platform like WineBottler on a mac, you will be very, very grateful that there's a program native to macs. In fact, there are two, and a very nice exploration of the differences can be found [here](http://www.jkarreth.net/files/bayes-cph_Tutorial-JAGS.pdf). I used R2jags, so this example will be for that program. 


Install the package
```{r install R2jags, message=F}
library(R2jags)
```


R2jags uses two files (which can also be combined in one):   
&nbsp;&nbsp;&nbsp;&nbsp;   * Model file    
&nbsp;&nbsp;&nbsp;&nbsp;   * Call file    

The model file is the regression formula written out in computer friendly language. Mine looks like this. I've tried to use the same greek letters in the model and I did in the equation above. I saved it in a file called "model_pat.R"
```{r, eval=FALSE}
model{
  
  # Main model level 1
  
  for (i in 1:N){ # N represents each individual
    fistr[i] ~ dnorm(mu[i], tau[country[i]])    # normal distribution to dependent variable with mean mu and sd tau
      mu[i] <- alpha[country[i]] +              # mean mu consists of the regression formula level 1: alpha + 
               beta1[country[i]]*hrwkMom.c[i] +   # beta1*hrwkMom + beta2-7*X
               beta2*femedu.c[i] + beta3*maleage.c[i] + beta4*maleedu.c[i] + beta5*numkid[i] + beta6*anyunder4[i] + beta7*female[i] 
    
    # we can impute missings by giving them the same distribution observed in the non-missing data
    hrwkMom[i] ~ dnorm(-0.98, 0.01) 
    femedu[i] ~ dnorm(0.13, 0.84)
    age[i] ~ dnorm(-5.08, 0.02)
    edu[i] ~ dnorm(0.06, 0.84)
    under4[i] ~ dbern(0.31)
  }
  
  
  # Main model level 2
  
  for (j in 1:J){ # j represents each country
    beta1[j] ~ dnorm(mu.b1[j], tau.b1)          # coefficient for beta1 is given a normal distribution with mean mu.b1 and sd tau.b1
    mu.b1[j] <- gb0 + gb1*pat.c[j]              # mu.b1 consists of a constant gb0 and paternity leave
    tau[j] ~ dgamma(tau.a, tau.b)               # the sd tau from the level 1 equation also varies per country with a gamma distribution.
    alpha[j] ~ dnorm(mu.alpha[j], tau.alpha)    # alpha is the mean of father invovlement per country, which we model with a normal distribution with mean mu.alpha and sd tau.alpha
    mu.alpha[j] <- gamma0 + gamma1*pat.c[j]     # mu.alpha consists of a constant gamma0 and paternity leave
  }
  
  
  
  # Priors 
  # we need to specify prior distributions for each parameter which is not created from existing data. The distributions here represent conservative or "flat" priors. 
  tau.a ~ dgamma(1, 1)  
  tau.b ~ dgamma(1, 1)
  tau.b1 ~ dgamma(1, 1)
  tau.alpha ~ dgamma(1, 1)
  
  beta2 ~ dnorm(0, .01)
  beta3 ~ dnorm(0, .01)
  beta4 ~ dnorm(0, .01)
  beta5 ~ dnorm(0, .01)
  beta6 ~ dnorm(0, .01)
  beta7 ~ dnorm(0, .01)
  
  gb0 ~ dnorm(0, .01)
  gb1 ~ dnorm(0, .01)
  
  gamma0 ~ dnorm(0, .01)
  gamma1 ~ dnorm(0, .01)

  sigma.b <- 1/(tau.b)            # These are the error terms
  sigma.b1 <- 1/(tau.b1)
  sigma.alpha <- 1/(tau.alpha)
  
}

```

Phew! `r emo::ji("sweat_smile")`  Mathematic notation is much shorter! Now the hard part is done, we just need to write the call file. We start by defining some data
```{r, eval=FALSE}
# sample size
J <- length(unique(ggs.resp$country)) # number of countries
N <- length(ggs.resp$fistr) # number of people

# amount of iterations and burn-in period
itt <- 10000
bi <- 1000    # what this means is our model runs 10,000 iterations and drops the first 1000. This is because it can take some time to stabilize and we want to base our parameter estimates off of the stabilized coefficients. 


# data
data <- list(fistr = ggs.resp$fistr, country = ggs.resp$country, hrwkMom.c = ggs.resp$hrwkMom.c, 
             pat.c = ggs.resp$pat.c, maleage.c = ggs.resp$maleage.c, maleedu.c = ggs.resp$maleedu.c, 
             femedu.c = ggs.resp$femedu.c, numkid = ggs.resp$numkid, anyunder4 = ggs.resp$anyunder4, 
             female = ggs.resp$female, J = J, N = N)


# define initial values. The sampler works by starting with these initial values and then trying different values until it gets to the parameter estimates. We need to define an initial value for every parameter for which we defined a prior distribution in the model. It's good to specify multiple starting values in case the sampler gets stuck somewhere. This way it can start from two different places and then hopefully converge at similar values on each.
inits1 <- list(tau.a = 1, tau.b = 1, tau.b1 = 1, tau.alpha = 1, 
               beta2 = 0, beta3 = 0, beta4 = 0, beta5 = 0, beta6 = 0,
               beta7 = 0, 
               gb0 = 1, gb1 = 1, 
               gamma0 = 0, gamma1 = 0)
inits2 <- list(tau.a = .1, tau.b = .1, tau.b1 = .1, tau.alpha = .1, 
               beta2 = 1, beta3 = 1, beta4 = 1, beta5 = 1, beta6 = 1,
               beta7 = 1, 
               gb0= 0, gb1 = 0,
               gamma0 = 1, gamma1 = 1)
               
inits <- list(inits1, inits2)
```

Run the model. This takes a while for me (~20 minutes)
```{r, eval=FALSE}
fistr_pat <- jags(data, inits, model.file = "model_pat.R",
                     # all parameters specified here will be traced, meaning you will be able to see how well they converged, plus means and distributions
                     parameters  = c("alpha", 
                                     "beta1", "beta2", "beta3", "beta4", "beta5", 
                                     "beta6", "beta7",
                                     "gb0", "gb1",
                                     "gamma0", "gamma1",
                                     "sigma.b", "sigma.b1", "sigma.alpha"), 
                     n.chains = 2,    # number of chains (I chose 2)
                     n.iter = itt,    # number of iterations per chain (10,000)
                     n.burnin = bi,   # amount of iterations discarded (1,000)
                     n.thin = 10      # discards every 10th iteration. 
                     )
```

```{r, include=FALSE}
load("fistr_pat.RData")
```

<br>

## Did the model converge? 
fistr_pat is an rjags object consisting of 6 elements: 1) model; 2) BUGSoutput; 3) parameters.to.save; 4) model.file; 5) n.iter; 6) DIC
We can do a lot of cool things with this object, starting with visual diagnostics. Did the model converge?

```{r, results='hide', message=F, warning=F, error=F}
library(mcmcplots)
fistr_pat_mcmc <- as.mcmc(fistr_pat) # convert to type mcmc
```

```{r, results='hide', message=F, warning=F, error=F, eval=FALSE}
Plots <- mcmcplot(fistr_pat_mcmc, dir=".")
```

<a href="MCMCoutput.html" target="_blank">Plots</a>

The link above should open the trace plots in a new tab. If the pink and blue lines are evenly mixed, it means the parameter estimates converged. Deviance refers to how well the overall model converged. Visually, these look good!

<br>

## Results

The model seems trustworthy, now we can try to answer our research questions. Our first two questions were if fathers spend more time with children when their partners work more, and if this is true across all countries in our study. We'll produce a table eventually, but let's start by visualizing it. 


First extract the data we want to plot and create a list of countries
```{r}
# create summary data as list
sum <- summary(fistr_pat_mcmc)

# list of countries
countries <- c("Australia", "Austria", "Belgium", "Bulgaria", "CzechRepublic", 
               "Estonia", "France", "Georgia", "Germany", "Hungary", "Italy", 
               "Lithuania", "Netherlands", "Norway", "Poland", "Romania", "Russia")

```


Second, create a data frame of means, and lower and upper bounds for 95% credible intervals
```{r}
# create means for beta 1
meansb1 <- as.data.frame(sum$statistics[c(18:34), 1]) # the object fistr_pat_mcmc automatically averages results from both chains

# highest posterior density quantiles for beta 1
# the quantile output from fistr_pat_mcmc assumes normal distribution. As we saw in the trace plots, this is a reasonable assumption, but here I calculate the quantiles from the highest posterior density rather than assume normal distribution
quants12 <- HPDinterval(fistr_pat_mcmc) # we had two chains, so we have to do this twice 
quants1 <- quants12[[1]] # chain 1
quants2 <- quants12[[2]] # chain 2

# manipulate quantiles data to get it into a dataframe we can easily merge with means
quantsdf <- data.frame(quants1,quants2)
quantslower <- data.frame(rowMeans(quantsdf[,c(1,3)]))
quantsupper <- data.frame(rowMeans(quantsdf[c(2,4)]))
quants <- cbind.data.frame(quantslower,quantsupper)
quants$lower <- quants$rowMeans.quantsdf...c.1..3...
quants$rowMeans.quantsdf...c.1..3... <- NA
quants$upper <- quants$rowMeans.quantsdf.c.2..4...
quants$rowMeans.quantsdf.c.2..4... <- NA
quants <- quants[,c(3:4)]
quantsb1 <- quants[c(18:34),c(1,2)]

# combine meansb1 and quantsb1 in a dataframe
plot.datb1 <- as.data.frame(cbind(meansb1, quantsb1))
plot.datb1$means <- plot.datb1$`sum$statistics[c(18:34), 1]` # rename column
plot.datb1 <- plot.datb1[,c(2:4)] # erase extra columns

# using data labels from dataset
plot.datb1$countries <- countries
```

Order data and create range for x-axis
```{r}
# BETA1 ordering the values by size of effect
plot.datb1o <- plot.datb1[order(plot.datb1$means, decreasing=TRUE),]
plot.datb1o$countries <- reorder(plot.datb1o$countries, plot.datb1o$means)

# create range for x-axis in figure
rg <- diff(range(c(plot.datb1o$upper, plot.datb1o$lower)))
```

Graph beta 1
```{r}
# graph beta1
dotplot(countries ~ means, data=plot.datb1o ,scales=list(y=list(cex=.85)), xlim=c(min(plot.datb1o$lower)-.1*rg, max(plot.datb1o$upper)+.1*rg), xlab="effect", panel=function(x,y, subscripts){
  panel.abline(h = as.numeric(y), col = "gray80", lty = 10, v = 0)
  panel.segments(plot.datb1o$lower[subscripts], y, plot.datb1o$upper[subscripts], y, lty=1, col="gray40")
  panel.points(x,y, pch=16, col="black")})


```

To answer our questions:    
1. yes, fathers are more involved in childcare the more their partners work.    
2. no, this is not true for every country. In particular, in the Czech Republic and Estonia the credible intervals overlap 0, indicating that a significant proportion of fathers are not responsive to their wives' work hours. Furthermore, we can see clearly that the strength of this effect size varies across countries. In the Netherlands the mean effect size is .02 while in Estonia this is .002, and the confidence intervals don't overlap at all.     
3. This information allows us to answer our third question, whether paternity leave explains the cross-national difference in the effect size of mother's work hours on father involvement with a resounding NO. If it did, these lines would be aligned vertically.

The average effect size across all countries can be found by:
```{r}
# mean of beta1
mean(plot.datb1$means)
```

```{r}
# range of beta1
lower <- mean(plot.datb1$means) - 1.96*sd(plot.datb1$means)
upper <- mean(plot.datb1$means) + 1.96*sd(plot.datb1$means)
lower; upper
```

We can also get a table of results. Note that the parameter "hrwkMom * pat leave" is not significantly different from 0. 

```{r}
# table of results
coefficients <- (sum$statistics[c(18:40,42:48), 1])
quantiles <- (quants[c(18:40,42:48),c(1,2)])
analysisresults.fistr_pat <- as.data.frame(cbind(coefficients, quantiles))
analysisresults.fistr_pat$var <- c("hrwkMom.Austrlia", "hrwkMom.Austria", "hrwkMom.Belgium", "hrwkMom.Bulgaria", 
                                  "hrwkMom.CZ", "hrwkMom.Estonia", "hrwkMom.France", "hrwkMom.Georgia", 
                                  "hrwkMom.Germany", "hrwkMom.Hungary", "hrwkMom.Italy", 
                                  "hrwkMom.Lithuania", "hrwkMom.Netherlands", "hrwkMom.Norway", 
                                  "hrwkMom.Poland", "hrwkMom.Romania", "hrwkMom.Russia",
                                  "mom edu", "dad age", "dad edu", "number of kids", 
                                  "any kids under 4", "gender of respondent", 
                                  "intercept.country", "pat leave", "intercept beta1", 
                                  "hrwkMom * pat leave", "level2 var",
                                  "beta1 var", "level1 var")
analysisresults.fistr_pat
```

```{r, echo=FALSE}
# clean up workspace
rm(analysisresults.fistr_pat, fistr_pat, fistr_pat_mcmc, meansb1, plot.datb1, plot.datb1o, quantiles, quants, quants1, quants12, quants2, quantsb1, quantsdf, quantslower, quantsupper, sum, coefficients, countries, lower, rg, upper)
```

<br>

## Conclusions

For those policymakers and other ecological fallacy-makers who see higher averages of women's labor market participation and men's childcare and jump to the conclusion that this is caused by a reciprocal relationship at the household level: congratulations, you're right. 

### Mother's work hours and father involvement
However, this relationship is far from 1:1. The average effect size across countries was 0.01, meaning that for every hour a woman works, her husband is likely to do a .01 greater share of childcare. Because our measure of father involvement was on a scale of 0-4, this is higher than it sounds, but still fairly low. Let's compare households where one woman works 0 hours and another works 40 hours per week. All things being equal, the father whose wife works 40 hours will be 0.4 points on a 5-point scale more involved. 

### Paternity leave
Not shown here, but I tried about 20 different country level indicators and was never able to explain why men are more responsive to their wives' work hours in some countries than others. I tried measures of culture, policy, and economics, and nothing explained the difference between countries. In short, I had what is known as a null effect. After roughly four years of scratching my head about why there was no significant effect of not just paternity leave, but ANY macro-level variable, I finally came up with an answer:

European countries are pretty similar with regard to how responsive fathers are to their wives' work hours. Yes, there are significant cross-national differences, but the majority of countries don't differ from each other. If we exclude the Netherlands, Romania, Australia, and Estonia, we have a story of similarity rather than difference across countries.


<br>

## Map

And with that, I leave you with this map of total father involvement across Europe.

Load packages 
```{r, results='hide', message=F, warning=F, error=F}
library(ggplot2)
library(grid)
library(rworldmap)
library(mapproj)
library(xtable)
library(data.table)
```


Get world map and select countries
```{r}
# Get the world map
worldMap <- getMap()

# Countries we want in the map
MAPcountries <- c("Albania", "Andorra", "Austria", "Belarus", "Belgium", "Bosnia and Herz.", "Bulgaria", 
                  "Croatia", "Cyprus", "Czech Rep.", "Denmark", "Estonia", "Finland", "France",
                  "Georgia", "Germany", "Greece", "Hungary", "Iceland", "Ireland", "Italy", 
                  "Kosovo", "Latvia", "Liechtenstein", "Lithuania", "Luxembourg", 
                  "Macedonia", "Malta", "Moldova", "Monaco", "Montenegro", "Netherlands", "Norway",
                  "Poland", "Portugal", "Russia", "Romania", "San Marino", "Serbia", "Slovakia", "Slovenia", 
                  "Spain", "Sweden", "Switzerland", "Turkey", "Ukraine", "United Kingdom", "Vatican")

# Select only the index of countries in Europe
index <- which(worldMap$NAME%in%MAPcountries)

# Extract longitude and latitude border's coordinates of members states of E.U. 
Coords <- lapply(index, function(i){
  df <- data.frame(worldMap@polygons[[i]]@Polygons[[1]]@coords)
  df$region =as.character(worldMap$NAME[i])
  colnames(df) <- list("long", "lat", "region")
  return(df)
})

Coords <- do.call("rbind", Coords)
```

Create data frame of aggregated data
```{r}
# create a data frame aggdat with aggregated mean of fistr per country
tempdata <- data.table(ggs.resp$fistr, ggs.resp$country)
aggdat <- tempdata[, mean(V1, na.rm=T), by=V2] 
aggdat$fistr <- aggdat$V1
aggdat$V1 <- NULL
aggdat$country <- aggdat$V2
aggdat$V2 <- NULL

# create fistr column in Coords dataframe
Coords$fistr <- aggdat$fistr[match(Coords$region,aggdat$country)]
```

Plot the map
```{r}
P <- ggplot() + geom_polygon(data = Coords, aes(x = long, y = lat, group = region, fill = fistr),
                             colour = "black", size = 0.1) +
     coord_map(xlim = c(-13, 35),  ylim = c(32, 71))

P <- P + scale_fill_gradient(name = "Father's share of childcare", low = "#FF0000FF", high = "#FFFF00FF", na.value = "grey50")

P <- P + theme(panel.grid.minor = element_line(colour = NA), 
  panel.background = element_rect(fill = NA, colour = NA),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(), axis.ticks.x = element_blank(),
  axis.ticks.y = element_blank(), axis.title = element_blank(),
  rect = element_blank(),
  plot.margin = unit(0 * c(-1.5, -1.5, -1.5, -1.5), "lines"))
P
```



This blog post can be found on [GitHub](https://github.com/brettory/R2Jags)

<sup id="fn1">1. Countries in the analysis include: Australia, Austria, Belgium, Bulgaria, Czech Republic, Estonia, France, Georgia, Germany, Hungary, Italy, Lithuania, Netherlands, Norway, Poland, Romania, Russia  <a href="#ref1" title="Jump back to footnote 1 in the text.">↩</a></sup> 






