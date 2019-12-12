#Link RStudio (https://rstudio.cloud/project/792398)
#Link Github (https://github.com/BragatteMAS/R-Quantitative_Genetics)
# Aula 1 - Genética quantitativa - conceitos básicos e revisão ANOVA

#Install.packages
install.packages("installr"); #error as 'lib' is unspecified (is not available (for R version 3.6.0))
install.packages("ExpDes");
install.packages("MASS");
install.packages("laercio")
install.packages("corrplot");
install.packages("PerformanceAnalytics");

updateR() # updating R.(True) = auto

#Activate packages
library(installr) #is not available (for R version 3.6.0)
library(ExpDes)
library(MASS)
library(laercio)
library(corrplot)
library(PerformanceAnalytics)

#Upload file
my_data <- read.delim(file.choose()) 

exp_gen_seca <- read.delim(file="1_exp_gen_seca.txt")

dados <- exp_gen_seca

attach(dados)
names(dados)
is.factor(bloco)
as.factor(bloco)
is.factor(genotipo)
as.factor(genotipo)
is.factor(tratamento)
is.numeric(DPV3)

tapply(PR, tratamento, sd)
tapply(PR, tratamento, var)
tapply(PR, tratamento, sd)/ sqrt(36)

#TESTES DE NORMALIDADE
shapiro.test(DPV3)

shapiro.test(PR)

shapiro.test(H)

shapiro.test(MST)

shapiro.test(MSR)

bartlett.test(MST ~ tratamento)

bartlett.test(MST ~ genotipo)

bartlett.test(MST ~ bloco)

require(MASS) 
boxcox(MST ~ genotipo) 
locator()

MSTt = ((MST ^ 0.06246118)-1)/0.06246118

bartlett.test(MSTt ~ genotipo)

shapiro.test(MSTt)

an = aov(MSTt ~ bloco + tratamento + genotipo + tratamento*genotipo)
summary(an)

#TESTES DE MÉDIAS
require(laercio)
LTukey(an)
LDuncan(an)
LScottKnott(an, "genotipo")
LScottKnott(an, "tratamento")

#SE NÃO DER NORMALIDADE, TESTES NÃO-PARAMÉTRICOS
kruskal.test(PR ~ genotipo) #caso não tenha normalidade
kruskalmc(PR ~ genotipo) #caso não tenha normalidade ?

#PACOTE EXPDES
require(ExpDes)
fat2.rbd(genotipo, tratamento, bloco, DPV3)
fat2.rbd(genotipo, tratamento, bloco, DPV3)
fat2.rbd(genotipo, tratamento, bloco, PR)

fat2.crd(factor1, factor2, resp) #  para fatorial duplo em delineamento em blocos casualizados ?
fat2.rbd(factor1, factor2, block, resp) # para fatorial triplo em delineamento inteiramente casualizado ?
fat3.crd(factor1, factor2,factor3,resp) # para fatorial triplo em delineamento em blocos casualizados? 
fat3.rbd(factor1, factor2, factor3, block, resp) # ?

#<<<---------------------------------------------------------------------------------->>>

# Aula 2 - Componentes genéticos de médias e variâncias

#Aula teórica / prática ao final do script da aula 1

#Correlation

correlacao <- read.delim(file.choose()) or

correlacao <- read.delim(file="2_correlacao.txt")

res <- cor(correlacao)
round(res, 2)
library(corrplot)
corrplot(res, type="upper", order="hclust", tl.col="black",tl.srt=45)
library("PerformanceAnalytics")
chart.Correlation(correlacao, histogram=TRUE, pch=19)

#<<<---------------------------------------------------------------------------------->>>
 
# Aula 3 - Análise de experimentos e estimativa de parâmetros genéticos
#install.packages
install.packages("ExpDes.pt")
install.packages("ExpDes")  
install.packages("MASS")  
install.packages("corrplot")  
install.packages("PerformanceAnalytics")  

#Activet.packages
library("ExpDes.pt")
library("ExpDes")
library("MASS")  
library("corrplot")  
library("PerformanceAnalytics")

#PARA ANÁLISE DO EXPERIMENTO EM BLOCOS CASUALIZADOS
require(ExpDes.pt)
dados_e = read.table(file="3_experimento.txt", sep="", header=TRUE)
attach(dados_e)
names(dados_e)
dbc(genotype, block, plant_height, quali = TRUE, mcomp = "tukey", nl=FALSE,
    hvar='oneillmathews', sigT = 0.05, sigF = 0.05)

shapiro.test(plant_height)
bartlett.test(plant_height ~ genotype)
require(MASS)
boxcox(plant_height ~ genotype)
locator()

#PARA CORRELAÇÃO LINEAR
dados_m = read.table(file="3_molusco.txt", sep="", header=TRUE)
attach(dados_m)
cor(P,O)

#data3 <- read.delim(file.choose()) #pode chamar arquivo "molusco' ou novo arquivo(data3)
res <- cor(dados_m)
round(res, 2)
library(corrplot)
corrplot(res, type="upper", order="hclust", tl.col="black",tl.srt=45)
library("PerformanceAnalytics")
chart.Correlation(dados_m, histogram=TRUE, pch=19)

#PARA REGRESSÃO LINEAR
dado_m = read.table(file="3_molusco.txt", sep="", header=TRUE)
attach(dados_m)
names(dados_m)
plot(P, O,xlab="P(mm)",ylab="O(mm)",cex.lab=2,cex.axis=2,bty="l",cex=2,pch=16, col="skyblue3")
?pch
a = lm(O~P)
abline(a, col="red3, lw2=2")
text(locator(), "Equação de regressão", cex=2, col="red")
summary(a)

#<<<---------------------------------------------------------------------------------->>>

# Aula 4 - Análise de experimentos e estimativa de parâmetros genéticos
#install.packages
install.packages("lme4")
install.packages("sommer")
install.packages("laercio")
install.packages("onemap")

#Activet.packages
library("lme4")
library("sommer")
library("laercio")
library("onemap")

dados4 = read.table(file="4.1_DT_EXP.txt", sep="", header=TRUE)
attach(dados4)

#Para dados balanceados (mesmo número de repetições, mesmo número de unidades experimentais)
anova = aov(produtividade ~ genotipo + ambiente + bloco + genotipo:ambiente + ambiente:bloco)
summary(anova)

#Para dados desbalanceados (algumas repetições foram perdidas) 
modelo1 <- lm(produtividade ~ ambiente + ambiente/bloco + genotipo + genotipo:ambiente, data=dados4)
anova(modelo1)

#Testes complementares
require(laercio)
LTukey(anova)

#Modelo misto: quando há fatores fixos e aleatórios a serem analisados dentro de um mesmo modelo
require(lme4)
#Modelo em que a variável PRODUTIVIDADE é função do GENOTIPO (fator fixo), AMBIENTE (fator fixo) + componente de interação AMBIENTE:GENÓTIPO (fator aleatório) e erro experimental (resíduo). 
m0 <- lmer(produtividade ~ genotipo + ambiente + (1|ambiente:genotipo),data=dados4)

#Testando a significância do modelo: TESTE DA RAZÃO DE VEROSSIMILHANÇA para ver se há interação de genótipos por ambientes. 
m0 <- lmer(produtividade ~ genotipo + ambiente + (1|ambiente:genotipo),data=dados4, REML=FALSE)
anova(m0) #A ANOVA para o modelo m0, individualmente, não vai mostrar p-valor. Uma explicação para o porquê disso pode ser encontrada no seguinte site: https://stat.ethz.ch/pipermail/r-help/2006-May/094765.html

#Criamos o modelo m1, usando lm somente, pois deixa de ser um modelo misto e passa a ser apenas linear com efeitos fixos. 
m1 <- lm(produtividade ~ genotipo + ambiente,data=dados4)
anova(m1)

#Vamos fazer uma ANOVA comparando os dois modelos.
anova(m0,m1)

#No nosso conjunto de dados há interação de genótipos por ambientes. Podemos visualizar isso graficamente por: 
xyplot(produtividade~ambiente, groups=genotipo, data=dados4, type="o", cex=1)

#Computando a herdabilidade para este experimento.
ans1 <- mmer(produtividade~1, random= ~ genotipo + ambiente + genotipo:ambiente + ambiente:bloco, rcov= ~ units, data=dados4)

summary(ans1)$varcomp
(n.env <- length(levels(dados4$ambiente)))

#Computando a herdabilidade. 
pin(ans1, h2 ~ V1 / ( V1 + (V3/n.env) + (V5/(2*n.env)) ) )

#Cana de açucar...
dados_c = read.table(file="4.2_clones.txt", sep="", header=TRUE)
dados_ca = read.table(file="4.3_cana_acucar_2.txt", sep="", header=TRUE)
attach(dados_c)
attach(dados_ca)
# 1.lm or anova
anova = aov(rendimento_acucar~1, random= ~ clone + corte + clone:corte + corte:bloco)
summary(anova)
modelo1 <- lm(rendimento_acucar~1, random= ~ clone + corte + clone:corte + corte:bloco, data=dados_ca)
anova(modelo1)
# 2.Análises complementares de médias usando os testes de Tukey e Scott-Knott (pacote “laercio”).
require(laercio)
LTukey(anova)
LDuncan(anova)
LScottKnott(anova, "corte")
LScottKnott(anova, "clone")
# 3.o fator clone como aleatório, assim como a interação entre o corte e o clone (aleatório).
##O fator corte continuará sendo fixo. (pacote lme4) interação genótipos/ambientes é significativa?
require(lme4)
#Modelo em que a variável rendimento_acucar é função do CORTE (fator fixo), CLONE (fator aleatório) + componente de interação CORTE:CLONE (fator aleatório) e erro experimental (resíduo). 
m0 <- lmer(rendimento_acucar ~ corte + clone + (1|corte:clone),data=dados_ca)
anova(m0)
m1 <- lm(rendimento_acucar ~ corte + clone,data=dados_ca)
anova(m1)
##Compare os resultados obtidos nas tarefas (1) e (2). Qual mais adequado modelo fixo ou aleatório ?
anova(m0,m1)
##Esboce graficamente o rendimento de açúcar de todos os clones em relação a cada ambiente (corte).
xyplot(rendimento_acucar~clone, groups=clone, data=dados_ca, type="o", cex=1)
# 4.(pacote sommer) - compute os componentes de variância:
ans1 <- mmer(rendimento_acucar~1, random= ~ clone + corte + clone:corte + corte:bloco)
summary(ans1)$varcomp
(n.env <- length(levels(dados_ca$clone)))

#determine a herdabilidade ampla a partir da seguinte relação: 
pin(ans1, h2 ~ V1 / ( V1 + (V3/3) + (V4/(3*n.env)) ) )

#<<<---------------------------------------------------------------------------------->>>

# Aula 5 - Marcadores moleculares, seleção assistida e mapeamento de ligação

#Abrindo o onemap 
library(onemap) 

#Escolha o diretório  
data("5_onemap_example_f2") 

#Execute o arquivo para ver o seu conteúdo 
#Importanto um arquivo vcfR 
library(vcfR) 
vcfR.object <- read.vcfR("5_vcf_example_f2.vcf") 

#Convertendo a um objeto onemap 
vcf_example_f2 <- onemap_read_vcfR(vcfR.object = vcfR.object, 
                                   parent1 = "P1",  
                                   parent2 = "P2",  
                                   cross = "f2 intercross") 

#Verificando se a conversão foi realizada  
class(onemap_example_f2) 
class(vcf_example_f2) 

#Visualizando os dados graficamente  
plot(onemap_example_f2)
plot(vcf_example_f2) 
plot_by_segreg_type(onemap_example_f2) 
plot_by_segreg_type(vcf_example_f2) 

#Combinando os dois arquivos de genotipagem  
comb_example <- combine_onemap(onemap_example_f2, vcf_example_f2) 
comb_example 

plot(comb_example) 

#Exportando o novo arquivo – com a combinação 
write_onemap_raw(comb_example, file.name = "new_dataset_bragatte.raw", cross="f2 intercross") 

#Testando a segregação dos marcadores moleculares  
f2_test <- test_segregation(comb_example) 
class(f2_test) 
f2_test 
print(f2_test) 

#Para verificar distorções de segregação  
Bonferroni_alpha(f2_test) #se passar do marcador rejeita, não tem segre mendel
plot(f2_test) 

#Quais marcadores não têm distorção de segregação e quais não têm? 
select_segreg(f2_test) 
select_segreg(f2_test, distorted = TRUE) 

#Estimando frequências de recombinação entre dois marcadores  
twopts_f2 <- rf_2pts(comb_example) 
(LOD_sug <- suggest_lod(comb_example)) 

#Visualizando frequências de recombinação entre marcadores específicos  
print(twopts_f2, c("M5", "M10")) 
##abaixo do limite anterior de LOD, provável estar ligado

#Designando marcadores aos grupos de ligação. Quando você tem informações sobre a posição que certos marcadores ocupam no cromossomo, você pode indica-las. Caso contrário, o pacote tem outra estratégia para realizar o mapeamento.  
mark_all_f2 <- make_seq(twopts_f2, "all") 

#Constituindo os grupos de ligação  
LGs_f2 <- group(mark_all_f2) 
LGs_f2 
(LGs_f2 <- group(mark_all_f2, LOD = LOD_sug, max.rf = 0.5)) 
##ordenado os marcadores, começar a trabalhar pelos menores

#Ordenando os marcadores dentro dos grupos de ligação. Para isso, escolha a função de correção que você deseja.  
set_map_fun(type = "kosambi") 
##set_map_fun(type = "haldane") outra opção

#Vamos iniciar ordenando os marcadores do grupo 2  
LG2_f2 <- make_seq(LGs_f2, 2) 
LG2_f2 

#Vários métodos possíveis para ordenar. Vamos usar o Rapid Chain Delineation (Doerge, 1996) 
LG2_rcd_f2 <- rcd(LG2_f2) ##número do marcador e ordem no mapa
LG2_f2_ord <- order_seq(input.seq = LG2_f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3) 
LG2_f2_ord ##** mais asteríscos, mais chances
(LG2_f2_all <- make_seq(LG2_f2_ord, "force")) ## forçar marcador a ficar na posição mais significativa 
##não a evidência estatística de que o marcador estará presente!

#Para uma ordem mais “segura” dos marcadores  
LG2_f2_ord <- order_seq(input.seq = LG2_f2, n.init = 5, 
                        subset.search = "twopt", 
                        twopt.alg = "rcd", THRES = 3,
                        touchdown = TRUE) 
(LG2_f2_final <- make_seq(LG2_f2_ord, "force")) 

#Versão final do grupo de ligação  
LG2_f2_final 

#Linkage group 1 
LG1_f2 <- make_seq(LGs_f2, 1) 
LG1_f2_ord <- order_seq(input.seq = LG1_f2, n.init = 5, 
                        subset.search = "twopt", 
                        twopt.alg = "rcd", THRES = 3, 
                        touchdown = TRUE) 
(LG1_f2_final <- make_seq(LG1_f2_ord, "force")) 
LG1_f2_final 

#Linkage group 3 
LG3_f2 <- make_seq(LGs_f2, 3) 
LG3_f2_ord <- order_seq(input.seq = LG3_f2, n.init = 5, 
                        subset.search = "twopt", 
                        twopt.alg = "rcd", THRES = 3, 
                        touchdown = TRUE) 
(LG3_f2_final <- make_seq(LG3_f2_ord, "force")) 
LG3_f2_final 

#Heatmaps de matrizes de frequências de recombinação  
#Vamos colocar um marcador em uma posição errada  
temp_seq <- drop_marker(LG3_f2_final, 38) 
(temp_seq <- add_marker(temp_seq, 38)) 
(LG3_f2_wrong <- map(temp_seq)) 
rf_graph_table(LG3_f2_wrong) ##evidência locus efetivamente ligados

#Como temos evidência de que esse marcador está na posição errada, podemos tentar remapeá-lo. 
temp_seq <- drop_marker(LG3_f2_wrong, 38) 
temp_map <- map(temp_seq) 
(temp_try <- try_seq(temp_map, 38)) 
(LG3_f2_final <- make_seq(temp_try, 4)) 

#Checagem final  
rf_graph_table(LG1_f2_final) 

#Vamos desenhar o mapa final? 
maps_list <- list(LG1_f2_final, LG2_f2_final, LG3_f2_final) 
draw_map(maps_list, names = TRUE, grid = TRUE, cex.mrk = 0.7) 

#É possível desenhar mapas para um grupo de ligação apenas  
draw_map(LG1_f2_final, names = TRUE, grid = TRUE, cex.mrk = 0.7) 

#Output em PDF  
draw_map2(LG1_f2_final, LG2_f2_final, LG3_f2_final, main = "Only linkage information",  
          group.names = c("LG1", "LG2", "LG3"), output="map.pdf") 

#Visualizando os dados graficamente  
plot(onemap_example_f2) 
plot(vcf_example_f2) 
plot_by_segreg_type(onemap_example_f2) 
plot_by_segreg_type(vcf_example_f2) 

#<<<---------------------------------------------------------------------------------->>>#

# Aula 6 - MAPA DE LIGAÇÃO ASSOCIADO A MAPEAMENTO DE QTL

#Construa o mapa de ligação com o pacote onemap 
fake.f2.onemap <- read_mapmaker(file="fake.f2.onemap.raw")
twopts.f2 <- rf_2pts(fake.f2.onemap)
mark.all.f2 <- make_seq(twopts.f2, "all")
(LGs.f2 <- group(mark.all.f2, LOD=3, max.rf=0.5))
set_map_fun(type="kosambi")

LG2.f2 <- make_seq(LGs.f2, 2)
LG2.rcd.f2 <- rcd(LG2.f2)

LG2.f2.ord <- order_seq(input.seq=LG2.f2, n.init = 5,subset.search = "twopt",twopt.alg = "rcd", THRES = 3)
LG2.f2.ord
LG2.f2.safe <- make_seq(LG2.f2.ord,"safe")
(LG2.f2.all <- make_seq(LG2.f2.ord,"force"))
LG2.f2.ord <- order_seq(input.seq=LG2.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        touchdown=TRUE)
(LG2.f2.final<-make_seq(LG2.f2.ord, "force"))
LG2.f2.final
LG1.f2 <- make_seq(LGs.f2, 1)
LG1.f2.ord <- order_seq(input.seq=LG1.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        touchdown=TRUE)
(LG1.f2.final <- make_seq(LG1.f2.ord,"force"))
LG1.f2.final
LG3.f2 <- make_seq(LGs.f2, 3)
LG3.f2.ord <- order_seq(input.seq=LG3.f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        touchdown=TRUE)
(LG3.f2.final <- make_seq(LG3.f2.ord,"force"))
LG3.f2.final
maps.list<-list(LG1.f2.final, LG2.f2.final, LG3.f2.final)
draw_map(maps.list, names= TRUE, grid=TRUE, cex.mrk=0.7)
write_map(maps.list, "fake.f2.onemap.map")


#Agora, instale o pacote “qtl” e o abra. É hora de realizar o mapeamento de QTL 
library("qtl")
require(qtl)
fake.f2.qtl <- read.cross("mm", file="fake.f2.onemap.raw", mapfile="fake.f2.onemap.map")

#Vamos estimar um outro mapa no programa qtl
newmap <- est.map(fake.f2.qtl, tol=1e-6, map.function="kosambi")
newmap

#Visualize os dois mapas comparativamente. Houve muitas diferenças? 
plot.map(fake.f2.qtl, newmap)

#Agora, determine QTL com base nos métodos EM (EXPECTATION MAXIMIZATION) e HK (REGRESSÃO DE HALEY e KNOTT)
fake.f2.qtl <- calc.genoprob(fake.f2.qtl, step=2)
out.em <- scanone(fake.f2.qtl, method="em")
out.hk <- scanone(fake.f2.qtl, method="hk")
plot(out.em, out.hk, col=c("blue","red"))
write.cross(fake.f2.qtl, format="qtlcart", filestem="fake.f2.onemap")
max(out.hk)

#<<<---------------------------------------------------------------------------------->>>#

#7_Mapeamento de associação

#INSTALAÇÃO 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

install.packages("gplots")
install.packages("LDheatmap")
install.packages("genetics")
install.packages("ape")
install.packages("EMMREML")
install.packages("scatterplot3d")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("multtest")

source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("http://zzlab.net/GAPIT/emma.txt")

#ABRA OS SEGUINTES PACOTES
library(multtest)
library(gplots)
library(LDheatmap)
library(genetics)
library(compiler) #this library is already installed in R
library("scatterplot3d")

#ANÁLISE
#setwd("C:\\Users\\....")
myY <- read.table(file="mdp_traits.txt", head = TRUE) # y = pheontip
myG <- read.table(file="mdp_genotype.txt" , head = FALSE) #g = genotip
myGAPIT <- GAPIT(Y=myY, G=myG, PCA.total=3)
