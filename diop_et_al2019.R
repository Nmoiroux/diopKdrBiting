# This is code to replicate the analyses and figures from the article entitled 
# "Sub-lethal insecticide exposure affects host biting behaviour of Kdr-resistant 
# Anopheles gambiae" by DIOP Malal M, CHANDRE Fabrice, ROSSIGNOL Marie, 
# PORCIANI Angélique, CHATEAU Mathieu, MOIROUX Nicolas and PENNETIER Cédric (2019)
# Code developed by Nicolas Moiroux: nicolas.moiroux@ird.fr

library(tidyverse)
library(emmeans)
library(lme4)
library(glmmTMB)
library(coxme)
library(Hmisc)
library(multipanelfigure)
library(lubridate)
source("modifying-facet-scales-in-ggplot2.R")

###### load data
biting <- read.delim("data_diop_2019.txt")
### metadata: dataframe "biting" has 13 columns which are described below:
# 	Date : date in days after 1900-01-01
# 	Trials : a unique code for each mosquito (primary key, the first character correspond to the treatment (C: control, P: permanet or O: Olysetç the second character is an order number)
# 	genotype : genotype for the kdr mutation (RR, RS or SS) of the tested mosquito 
# 	ttmt : treatment of pre-exposure (Control = Untreated Net, Olyset = Permethrin teated net, Permanet = Deltamethrin Treated net)
# 	perf : biting succes (0= unfed, 1 = fed)
# 	T_KD : time to be recorded as KD (in seconds)
# 	T_prob : total duration of probing (in seconds)
# 	N_prob : number of probing events
# 	T_feed : duration of feeding (in seconds)
# 	N_feed : number of feeding events
# 	T_predi : duration of pre-diuresis (in seconds)
# 	V_bloodmeal : volume of blood (µL)
# 	Av_weight : Average weight of 5 mosquitoes from the same rearing cage (mg)

biting$Date <- as.factor(biting$Date)
biting$kd <- as.numeric(biting$T_KD > 0)
biting$kd[is.na(biting$kd)] <- 0
biting$ttmt <- fct_rev(biting$ttmt)

## calculate veighted volume of blood
biting$V_Weight <- biting$V_bloodmeal/biting$Av_weight

###### Table 1 ----
fed <- biting %>% filter(perf=="1") %>%     # count fed mosquitoes among treatment and genotypes
	group_by(genotype, ttmt) %>% 
	summarise(fed_kd=sum(kd),fed = n()) %>% 
	mutate(fed_unkd=fed-fed_kd)

unfed <- biting %>% filter(perf=="0") %>%   # count unfed mosquitoes among treatment and genotypes
	group_by(genotype, ttmt) %>% 
	summarise(unfed_kd=sum(kd), unfed = n())  %>% 
	mutate(unfed_unkd=unfed-unfed_kd)

table1 <- left_join(fed,unfed, by=c("genotype", "ttmt")) %>% # join fed and unfed tables
	mutate(total=fed+unfed) %>% select(-fed, -unfed)

table1 %>% mutate(total_kd=fed_kd+unfed_kd) %>% 							# compute knockdown rates and 95% confidence interval
	mutate(kdR=binconf(total_kd,total,alpha = 0.05)[,1]) %>% 
	mutate(kdR_lo=binconf(total_kd,total,alpha = 0.05)[,2]) %>%
	mutate(kdR_hi=binconf(total_kd,total,alpha = 0.05)[,3])

###### binomial mixed-effect model of feeding success  ----
glmm_perf <- glmer(perf~genotype*ttmt + (1|Date), data=biting, family=binomial)

# multiple comparisons
summary(emmeans(glmm_perf, pairwise~genotype | ttmt, type="response"), infer = TRUE) # among genotypes comparison
summary(emmeans(glmm_perf, pairwise~ttmt | genotype, type="response"), infer = TRUE) # among treatments comparison

###### Zero-trunctated Negative-Binomial (ZTNB) mixed-effect model of the number of probing attempts -----
glmm_Nprob <- glmmTMB(N_Prob~genotype*ttmt + (1|Date), data=biting, family=truncated_poisson(link = "log") ) # ZT poisson model
glmm_Nprob_nb <- glmmTMB(N_Prob~genotype*ttmt + (1|Date), data=biting, family=truncated_nbinom2(link = "log") ) # ZT NB
anova(glmm_Nprob, glmm_Nprob_nb) # model NB is the best

# multiple comparisons
summary(emmeans(glmm_Nprob_nb, pairwise~genotype | ttmt, type="response"), infer = TRUE) # among genotypes comparison
summary(emmeans(glmm_Nprob_nb, pairwise~ttmt | genotype, type="response"), infer = TRUE) # among treatments comparison

###### Cox proportional hazard mixed-effect model of the probing duration-----
coxm_prob<-coxme(Surv(T_Prob,perf)~genotype*ttmt + (1|Date), biting) 

# multiple comparisons
summary(emmeans(coxm_prob, pairwise~genotype | ttmt, type="response"), infer = TRUE) # among genotypes comparison
summary(emmeans(coxm_prob, pairwise~ttmt | genotype, type="response"), infer = TRUE) # among treatments comparison

###### Cox proportional hazard mixed-effect model of the feeding duration-----
coxm_feed<-coxme(Surv(T_feed,perf)~genotype*ttmt + (1|Date), biting) 

# multiple comparisons
summary(emmeans(coxm_feed, pairwise~genotype | ttmt, type="response"), infer = TRUE) # among genotypes comparison
summary(emmeans(coxm_feed, pairwise~ttmt | genotype, type="response"), infer = TRUE) # among treatments comparison

###### Cox proportional hazard mixed-effect model of the prediuresis duration----
coxm_predi<-coxme(Surv(T_Predi,perf)~genotype*ttmt + (1|Date), biting) 

# multiple comparisons
summary(emmeans(coxm_predi, pairwise~genotype | ttmt, type="response"), infer = TRUE) # among genotypes comparison
summary(emmeans(coxm_predi, pairwise~ttmt | genotype, type="response"), infer = TRUE) # among treatments comparison

###### linear mixed-effect model of the weighted volume of bloodmeal----
lmm_V <- lmer(V_Weight~genotype*ttmt + (1|Date), data=biting)

# multiple comparisons
summary(emmeans(lmm_V, pairwise~genotype | ttmt, type="response"), infer = TRUE) # among genotypes comparison
summary(emmeans(lmm_V, pairwise~ttmt | genotype, type="response"), infer = TRUE) # among treatments comparison

###### feeding duration vs blood-meal size/pre-diuresis duration correlation analysis (supplementary Table 1 and 2) ----

### feeding duration vs blood-meal size
# mixed-effect linear model
lmm_V_feed <- lmer(V_Weight~T_feed*genotype*ttmt + (1|Date), data=biting)
# is the slope different than zero ?
slopeV <- summary(emtrends(lmm_V_feed, ~genotype*ttmt, var = "T_feed"), infer = TRUE) 

### feeding duration vs pre-diuresis duration
# mixed-effect cox model
lmm_P_feed <- coxme(Surv(T_Predi,perf)~T_feed*genotype*ttmt + (1|Date), data=biting) 
# is the slope diffrent than zero ?
slopeP <- summary(emtrends(lmm_P_feed, ~genotype*ttmt, var = "T_feed", transform ="response"), infer = TRUE) 

##### supplementary Table 1 and 2 ---
supTable1 <- as.data.frame(slopeV)
supTable2 <- as.data.frame(slopeP)
supTable1
supTable2

##### binomial model of KD  
biting_ttmt <- biting %>% filter(ttmt != c("Control")) %>% mutate(perf=as.factor(perf))			# subset of mosquitoes exposed to insecticitide-treated nets
brglm_kd <- glm(kd~genotype*ttmt*perf, data=biting_ttmt, family=binomial, method = "brglmFit")

summary(emmeans(brglm_kd , pairwise~genotype | ttmt +perf, type="response"), infer = TRUE) # among genotypes comparison
summary(emmeans(brglm_kd , pairwise~ttmt | genotype, type="response"), infer = TRUE) # among treatments comparison


###### data for Figures-----
Perf_Biting <- biting
Perf_Biting$genotype <- fct_rev(Perf_Biting$genotype) # reverse factors level order
Perf_Biting$ttmt <- fct_rev(Perf_Biting$ttmt)			    # reverse factors level order

# table used for Figure 2A, sup Figure 1 and Figure 3
df.graph.bit<- Perf_Biting %>% group_by(ttmt, genotype) %>% summarise(P=sum(perf==1),Total=length(perf)) # count fed (data for figure 2A and 3)
df.graph.bit$perf<- binconf(df.graph.bit$P,df.graph.bit$Total,alpha = 0.05)[,1]                          # compute performanace rate
df.graph.bit$lower<- binconf(df.graph.bit$P,df.graph.bit$Total,alpha = 0.05)[,2]												 # and binomial CI
df.graph.bit$upper<- binconf(df.graph.bit$P,df.graph.bit$Total,alpha = 0.05)[,3]
df.graph.bit$genotype <-as.factor(df.graph.bit$genotype)
df.graph.bit$ttmt <- fct_recode(df.graph.bit$ttmt , Deltamethrin="Permanet",Untreated = "Control", Permethrin="Olyset")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")   # colorblindness friendly palette

###### Figure 2 -----

df.Fig2A <- filter(df.graph.bit, ttmt == "Untreated")                       # select treatment to plot
Fig2A <- ggplot(df.Fig2A , aes(genotype, y=perf)) +   										# panel A
	geom_bar(aes(fill=genotype), position = "dodge", stat= "identity" )+
	geom_errorbar(aes(ymin = lower, ymax = upper), width=0.2)+
	scale_fill_manual(values=cbPalette)+
	ylab("Feeding success rate")+
	theme(legend.position="none") +
	ylim(0,1)+
	xlab("")


df.fig2BCD <- biting %>% gather("var", "value", T_feed, T_Predi, V_Weight)        # data for panel B, C and D 
df.fig2BCD$genotype <- fct_rev(df.fig2BCD$genotype)
lab2 <- c(T_feed = "Feeding duration (s)", T_Predi = "Prediuresis duration (s)", V_Weight = "Blood-meal size (µL/µg)")

Fig2B <- ggplot(filter(df.fig2BCD, var=="T_feed"), aes(x=genotype, y=value)) +   # panel B
	geom_boxplot(aes(fill=genotype))+
	scale_fill_manual(values=cbPalette)+
	ylab(lab2[1])+
	xlab("")+
	theme(legend.position="none")

Fig2C <- ggplot(filter(df.fig2BCD, var=="T_Predi"), aes(x=genotype, y=value)) +   # panel C
	geom_boxplot(aes(fill=genotype))+
	scale_fill_manual(values=cbPalette)+
	ylab(lab2[2])+xlab("")+
	theme(legend.position="none")

Fig2D <- ggplot(filter(df.fig2BCD, var=="V_Weight"), aes(x=genotype, y=value)) +   # panel D
	geom_boxplot(aes(fill=genotype))+
	scale_fill_manual(values=cbPalette)+
	ylab(lab2[3])+xlab("")+
	theme(legend.position="none")

figure2 <- multi_panel_figure(columns = 3, rows = 2, width = 210, height = 160)   # create multipanel figure
figure2 %<>%
	fill_panel(Fig2A, column = 1, row = 1:2) %<>%
	fill_panel(Fig2B, column = 2, row = 1) %<>%
	fill_panel(Fig2C, column = 3, row = 1) %<>%
	fill_panel(Fig2D, column = 2, row = 2)

figure2


###### Supplementary Figure 1 ----

ggplot(df.graph.bit, aes(genotype, group = ttmt,  y=perf)) +   
	geom_bar(aes(fill=genotype), position = "dodge", stat= "identity" )+
	geom_errorbar(aes(ymin = lower, ymax = upper), width=0.2)+
	scale_fill_manual(values=cbPalette)+
	ylab("Feeding success rate")+
	theme(legend.position="none") +
	ylim(0,1)+
	facet_grid(~ttmt)+
	xlab("")+
	theme(strip.background = element_blank(),strip.placement = "outside", panel.spacing.x=unit(2, "lines"),legend.position="none")



###### Figure 3 -----
figure3 <- ggplot(df.graph.bit, aes(ttmt, group = genotype,  y=perf)) +   # data used generated in Figure 2 § (above)
	geom_bar(aes(fill=ttmt), position = "dodge", stat= "identity" )+
	geom_errorbar(aes(ymin = lower, ymax = upper), width=0.2)+
	scale_fill_manual(values=rev(cbPalette))+
	ylab("Feeding success rate")+
	theme(legend.position="none") +
	ylim(0,1)+
	facet_grid(~genotype)+
	xlab("")+
	theme(strip.background = element_blank(),strip.placement = "outside", panel.spacing.x=unit(2, "lines"),legend.position="none")
figure3


###### Figure 4 ----
df_fig4 <- biting %>% gather("var", "value", T_feed, T_Predi, V_Weight) # arrange data for figure 4
df_fig4$ttmt <- fct_rev(df_fig4$ttmt)
df_fig4$genotype <- fct_rev(df_fig4$genotype)

lab3 <- c(RR="", RS="", SS="")
lab2 <- c(T_feed = "Feeding duration (s)", T_Predi = "Prediuresis duration (s)", V_Weight = "Blood-meal size (µL/mg)")
df_fig4$ttmt <- fct_recode(df_fig4$ttmt, Deltamethrin = "Permanet", Permethrin = "Olyset", Untreated = "Control")

figure4 <- ggplot(df_fig4, aes(x=ttmt, y=value)) +
	geom_boxplot(aes(fill=ttmt))+
	scale_fill_manual(values=rev(cbPalette))+
	facet_wrap(genotype~var, scales="free", strip.position = "left", labeller=labeller(var = lab2, genotype=lab3))+
	theme(strip.background = element_blank(),strip.placement = "outside", panel.spacing.y=unit(3, "lines"),legend.position="none")+
	ylab("")+xlab("")

### change y_axis limits (need to load the functions 'scale_override', 'CustomFacetWrap', 'facet_wrap_custom' that are in "modifying-facet-scales-in-ggplot2.R")
figure4 + facet_wrap_custom(genotype~var, scales="free", strip.position = "left", labeller=labeller(var = lab2, genotype=lab3), 
														scale_overrides = list(
															scale_override(1, scale_y_continuous(limits = c(0, 550))),
															scale_override(2, scale_y_continuous(limits = c(0, 450))),
															scale_override(3, scale_y_continuous(limits = c(0, 8))),
															scale_override(4, scale_y_continuous(limits = c(0, 550))),
															scale_override(5, scale_y_continuous(limits = c(0, 450))),
															scale_override(6, scale_y_continuous(limits = c(0, 8))),
															scale_override(7, scale_y_continuous(limits = c(0, 550))),
															scale_override(8, scale_y_continuous(limits = c(0, 450))),
															scale_override(9, scale_y_continuous(limits = c(0, 8)))
														))
figure4


##### END -----
