#Script debuged/refactoring by Bragatte.

#Preparation
#install.packages("laercio")
#install.packages(("MASS"))
#install.packages("sommer")
#install.packages("lme4")
#install.packages("ExpDes")
#install.packages("ExpDes.pt")
#install.packages("dplyr")

#Libraries
library(readxl)
library(laercio)
library(MASS)
library(sommer)
library(lme4)
library(ExpDes)
library(ExpDes.pt)
library(dplyr)

my_data <- read.delim(file.choose()) # or
setwd("/home/bragatte/Downloads/")

my_data<-read_xlsx(path = "my_data2.xlsx",col_names = TRUE,col_types = c("text","text","text","numeric","numeric","numeric","numeric"),na = c("","0"),skip = 1)
my_data$Block<-as.factor(my_data$Block)
is.factor(my_data$Block)
my_data$Genotype<-as.factor(my_data$Genotype)
is.factor(my_data$Genotype)
my_data$Treatment<-as.factor(my_data$Treatment)
is.factor(my_data$Treatment)
is.numeric(my_data$`10 days`)

# 10 days #
data10<-na.exclude(select(my_data,Block:`10 days`))
tapply(data10$`10 days`,data10$Genotype,sd) # 10 days
tapply(data10$`10 days`,data10$Block,sd)
tapply(data10$`10 days`, data10$Treatment, sd)
tapply(data10$`10 days`, data10$Treatment, var)
tapply(data10$`10 days`, data10$Treatment, sd)/ sqrt(nrow(data10))

m0 <- lmer(data10$`10 days` ~ (1|data10$Genotype) + data10$Treatment ,data=data10, REML = FALSE)
m1 <- lmer(data10$`10 days` ~ (1|data10$Genotype), data=data10, REML = FALSE)
m2<- lm(data10$`10 days`~ data10$Treatment,data = data10)
anova(m0,m1) #Test for treatment
anova(m0,m2) #Test for genotype

# 15 days #
data15<-na.exclude(select(my_data,Block:Genotype,`15 days`))
tapply(data15$`15 days`,data15$Genotype,sd) # 15 days
tapply(data15$`15 days`,data15$Block,sd)
tapply(data15$`15 days`, data15$Treatment, sd)
tapply(data15$`15 days`, data15$Treatment, var)
tapply(data15$`15 days`, data15$Treatment, sd)/ sqrt(nrow(data15))

m0 <- lmer(data15$`15 days` ~ (1|data15$Genotype) + data15$Treatment ,data=data15, REML = FALSE)
m1 <- lmer(data15$`15 days` ~ (1|data15$Genotype), data=data15, REML = FALSE)
m2<- lm(data15$`15 days`~ data15$Treatment,data = data15)
anova(m0, m1) #Test for treatment
anova(m0, m2) #Test for genotype

# 20 days #
data20<-na.exclude(select(my_data,Block:Genotype,`20 days`))
tapply(data20$`20 days`,data20$Genotype,sd) # 20 days
tapply(data20$`20 days`,data20$Block,sd)
tapply(data20$`20 days`, data20$Treatment, sd)
tapply(data20$`20 days`, data20$Treatment, var)
tapply(data20$`20 days`, data20$Treatment, sd)/ sqrt(nrow(data20))

m0 <- lmer(data20$`20 days` ~ (1|data20$Genotype) + data20$Treatment ,data=data20, REML = FALSE)
m1 <- lmer(data20$`20 days` ~ (1|data20$Genotype), data=data20, REML = FALSE)
m2 <- lm(data20$`20 days`~ data20$Treatment,data = data20)
anova(m0,m1) #Test for treatment
anova(m0,m2) #Test for genotype

# 30 days #
data30<-na.exclude(select(my_data,Block:Genotype,`30 days`))
tapply(data30$`30 days`,data30$Genotype,sd) # 30 days
tapply(data30$`30 days`,data30$Block,sd)
tapply(data30$`30 days`, data30$Treatment, sd)
tapply(data30$`30 days`, data30$Treatment, var)
tapply(data30$`30 days`, data30$Treatment, sd)/ sqrt(nrow(data30))

m0 <- lmer(data30$`30 days` ~ (1|data30$Genotype) + data30$Treatment ,data=data30, REML = FALSE)
m1 <- lmer(data30$`30 days` ~ (1|data30$Genotype), data=data30, REML = FALSE)
m2<- lm(data30$`30 days`~ data30$Treatment,data = data30)
anova(m0,m1) #Test for treatment
anova(m0,m2) #Test for genotype

# Heredability #
colnames(data10)<-c("Block","Treatment","Genotype","D10days")
colnames(data15)<-c("Block","Treatment","Genotype","D15days")
colnames(data20)<-c("Block","Treatment","Genotype","D20days")
colnames(data30)<-c("Block","Treatment","Genotype","D30days")
Heredability<-matrix(ncol = 2, nrow = 4)
colnames(Heredability)<-c("control","drought")
rownames(Heredability)<-c("10d","15d","20d","30d")
# Heredability 10 days #
control10 <- data.frame()
for (i in seq(1,nrow(data10))){
  for (j in seq(1,ncol(data10))){
    if (data10[i,2]=="control"){
      control10[i,j]<-data10[i,j]
    }}}
control10<-na.exclude(control10)
colnames(control10)<-c(colnames(data10))

drought10 <- data.frame()
for (i in seq(1,nrow(data10))){
  for (j in seq(1,ncol(data10))){
    if (data10[i,2]=="drought"){
      drought10[i,j]<-data10[i,j]
    }}}
drought10<-na.exclude(drought10)
colnames(drought10)<-c(colnames(data10))

an10control=aov(control10$D10days ~ control10$Genotype + control10$Block,data = control10)
n.env <- length(levels(control10$Block))
s10control<-summary(an10control)
Heredability[1,1]<-((s10control[[1]][1,3]-s10control[[1]][3,3])/n.env)/(s10control[[1]][1,3]/n.env)

an10drought=aov(drought10$D10days ~ drought10$Genotype + drought10$Block,data = drought10)
n.env <- length(levels(drought10$Block))
s10drought<-summary(an10drought)
Heredability[1,2]<-((s10drought[[1]][1,3]-s10drought[[1]][3,3])/n.env)/(s10drought[[1]][1,3]/n.env)

# Heredability 15 days
control15 <- data.frame()
for (i in seq(1,nrow(data15))){
  for (j in seq(1,ncol(data15))){
    if (data15[i,2]=="control"){
      control15[i,j]<-data15[i,j]
    }}}
control15<-na.exclude(control15)
colnames(control15)<-c(colnames(data15))

drought15 <- data.frame()
for (i in seq(1,nrow(data15))){
  for (j in seq(1,ncol(data15))){
    if (data15[i,2]=="drought"){
      drought15[i,j]<-data15[i,j]
    }}}
drought15<-na.exclude(drought15)
colnames(drought15)<-c(colnames(data15))

an15control=aov(control15$D15days ~ control15$Genotype + control15$Block,data = control15)
n.env <- length(levels(control15$Block))
s15control<-summary(an15control)
Heredability[2,1]<-((s15control[[1]][1,3]-s15control[[1]][3,3])/n.env)/(s15control[[1]][1,3]/n.env)

an15drought=aov(drought15$D15days ~ drought15$Genotype + drought15$Block,data = drought15)
n.env <- length(levels(drought15$Block))
s15drought<-summary(an15drought)
Heredability[2,2]<-((s15drought[[1]][1,3]-s15drought[[1]][3,3])/n.env)/(s15drought[[1]][1,3]/n.env)

# Heredability 20 days #
control20 <- data.frame()
for (i in seq(1,nrow(data20))){
  for (j in seq(1,ncol(data20))){
    if (data20[i,2]=="control"){
      control20[i,j]<-data20[i,j]
    }}}
control20<-na.exclude(control20)
colnames(control20)<-c(colnames(data20))

drought20 <- data.frame()
for (i in seq(1,nrow(data20))){
  for (j in seq(1,ncol(data20))){
    if (data20[i,2]=="drought"){
      drought20[i,j]<-data20[i,j]
    }}}
drought20<-na.exclude(drought20)
colnames(drought20)<-c(colnames(data20))

an20control=aov(control20$D20days ~ control20$Genotype + control20$Block,data = control20)
n.env <- length(levels(control20$Block))
s20control<-summary(an20control)
Heredability[3,1]<-((s20control[[1]][1,3]-s20control[[1]][3,3])/n.env)/(s20control[[1]][1,3]/n.env)

an20drought=aov(drought20$D20days ~ drought20$Genotype + drought20$Block,data = drought20)
n.env <- length(levels(drought20$Block))
s20drought<-summary(an20drought)
Heredability[3,2]<-((s20drought[[1]][1,3]-s20drought[[1]][3,3])/n.env)/(s20drought[[1]][1,3]/n.env)

# Heredability 30 days #
control30 <- data.frame()
for (i in seq(1,nrow(data30))){
  for (j in seq(1,ncol(data30))){
    if (data30[i,2]=="control"){
      control30[i,j]<-data30[i,j]
    }}}
control30<-na.exclude(control30)
colnames(control30)<-c(colnames(data30))

drought30 <- data.frame()
for (i in seq(1,nrow(data30))){
  for (j in seq(1,ncol(data30))){
    if (data30[i,2]=="drought"){
      drought30[i,j]<-data30[i,j]
    }}}
drought30<-na.exclude(drought30)
colnames(drought30)<-c(colnames(data30))

an30control=aov(control30$D30days ~ control30$Genotype + control30$Block,data = control30)
n.env <- length(levels(control30$Block))
s30control<-summary(an30control)
Heredability[4,1]<-((s30control[[1]][1,3]-s30control[[1]][3,3])/n.env)/(s30control[[1]][1,3]/n.env)

an30drought=aov(drought30$D30days ~ drought30$Genotype + drought30$Block,data = drought30)
n.env <- length(levels(drought30$Block))
s30drought<-summary(an30drought)
Heredability[4,2]<-((s30drought[[1]][1,3]-s30drought[[1]][3,3])/n.env)/(s30drought[[1]][1,3]/n.env)

readr::write_delim(as.data.frame(Heredability), path = "Heredability.txt", delim = "\t")

#<<<-------------------------------------QUESTION TWO--------------------------------------------->>>

setwd("~/Downloads/")
library(onemap)
mouse<-read_onemap(inputfile = "mouse_data.txt")
class(mouse)
class(mouse)

#Visualization 
plot(mouse)
plot_by_segreg_type(mouse)

#Testing segregation at molecular markers 
f2_test <- test_segregation(mouse)
class(f2_test)
f2_test
print(f2_test)

#Verifing distorcions at segregations 
Bonferroni_alpha(f2_test)

plot(f2_test)

#Which markers are distorted?
select_segreg(f2_test)
select_segreg(f2_test, distorted = TRUE)

#Recombining frequences
twopts_f2 <- rf_2pts(mouse)

(LOD_sug <- suggest_lod(mouse))

#Visualization of recombining frequences at molecular markers
print(twopts_f2, c("M5", "M10"))
write.table(twopts_f2$analysis,file = "Question2_rec_freq.txt")
