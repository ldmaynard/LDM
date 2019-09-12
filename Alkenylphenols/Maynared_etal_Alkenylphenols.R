rm(list=ls())


# Load Libraries ----------------------------------------------------------

library(dplyr)
library(tidyr)
library(Rmisc)
library(reshape)
library(multcomp)
library(MASS)
library(RColorBrewer)
library(ggplot2)
library(multcompView)
library(emmeans)
library(lme4)
library(AICcmodavg)
library(betareg)
library(car)
library(viridis)
library(knitr)
library(kableExtra)
library(tinytex)
library(ggpubr)
library(pals)


# Objective 2: Quantitative chemistry -------------------------------------

df_all <- read.csv(file="Maynard_etal_AlkenylphenolQuantChem.csv",head=TRUE)
#removing unnec. cols
df_all <- subset(df_all, select = -c(1, 6:32,43:52))
df_all[df_all == "#VALUE!"] <- NA


df_all$A_pdw<-as.numeric(as.character(df_all$A_pdw))
df_all$B_pdw<-as.numeric(as.character(df_all$B_pdw))
df_all$C_pdw<-as.numeric(as.character(df_all$C_pdw))
df_all$D_pdw<-as.numeric(as.character(df_all$D_pdw))
df_all$E_pdw<-as.numeric(as.character(df_all$E_pdw))
df_all$F_pdw<-as.numeric(as.character(df_all$F_pdw))
df_all$G_pdw<-as.numeric(as.character(df_all$G_pdw))
df_all$H_pdw<-as.numeric(as.character(df_all$H_pdw))
df_all$I_pdw<-as.numeric(as.character(df_all$I_pdw))
df_all$J_pdw<-as.numeric(as.character(df_all$J_pdw))

df_all$tissue<-as.factor(df_all$tissue)
df_all$plant<-as.factor(df_all$plant)
df_all$stage<-as.factor(df_all$stage)

#gather rows into columns, more data cleaning
datall <- gather(df_all, "compound", "per_dry_wt", 5:14)
datall$per_dry_wt<-as.numeric(datall$per_dry_wt)
datall$compound<-as.character(datall$compound)
datall$stage<-as.character(datall$stage)

#dry weight % -> proportions
datall$props<-(datall$per_dry_wt/100)

##splitting between tissue types
##only leaves
dfl <- filter(datall, tissue == "L")

##only pulp
dat <- filter(datall, tissue == "P") 


#labeling
datall <- datall[order(datall$tissue),]
datall <- datall[order(datall$stage),]
datall$tissue<-as.character(datall$tissue)
datall$plant<-as.character(datall$plant)

#removing unripe seeds
datall<-datall[-c(441:460, 671:690),]

#sep stages within tissues 
datall$tissue[1:210]<-"R"
datall$tissue[231:440]<-"U2"
datall$tissue[441:650]<-"U3"
datall$tissue[651:860]<-"F4"
datall$tissue[861:1020]<-"F5"
datall$tissue[1021:1170]<-"F6"

datall$tissue[datall$tissue=="R"]="Ripe pulp (1)"
datall$tissue[datall$tissue=="U2"]="Unripe pulp (2)"
datall$tissue[datall$tissue=="U3"]="Unripe pulp (3)"
datall$tissue[datall$tissue=="F4"]="Flowers (4)"
datall$tissue[datall$tissue=="F5"]="Flowers (5)"
datall$tissue[datall$tissue=="F6"]="Dev. flowers (6)"
datall$tissue[datall$tissue=="S"]="Seeds"
datall$tissue[datall$stage=="E"]="Exp. leaves"
datall$tissue[datall$stage=="M"]="Mature leaves"

datall$compound[datall$compound=="A_pdw"]="A"
datall$compound[datall$compound=="B_pdw"]="B"
datall$compound[datall$compound=="C_pdw"]="C"
datall$compound[datall$compound=="D_pdw"]="D"
datall$compound[datall$compound=="E_pdw"]="E"
datall$compound[datall$compound=="F_pdw"]="F"
datall$compound[datall$compound=="G_pdw"]="G"
datall$compound[datall$compound=="H_pdw"]="H"
datall$compound[datall$compound=="I_pdw"]="I"
datall$compound[datall$compound=="J_pdw"]="J"

a5<-aggregate(props~tissue+plant, data=datall, FUN=sum) 

beta1<-betareg(props~tissue, data=a5)
summary(beta1)
Anova(beta1)

d1<-emmeans(beta1,pairwise~tissue, type="response")
CLD(d1$emmeans,  Letters ='ABCDEFGHIJKLMNOPQRS')


##plot
tissueplot2<-ggplot(a5, aes(x=tissue, y=props)) +
	geom_boxplot() + geom_point(aes(color=as.character(plant)))+
	labs(x=" ", y="Total alkenylphenols (% dry wt)")+
	theme_minimal()+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds","Exp. leaves", "Mature leaves"))+
	scale_color_manual(values=as.vector(kelly(21)), name = "Plant ID")
tissueplot2

#tiff('supp_tissue_plot', units="in", width=10, height=4, res=500)
#tissueplot2
#dev.off()

##too busy, selecting representative stages for each tissue

##selected stages for each tissue
dat_select<-datall%>%
	filter(tissue==c("Ripe pulp (1)","Unripe pulp (2)", "Flowers (4)","Seeds","Mature leaves"))

a10<-aggregate(props~tissue+plant, data=dat_select, FUN=sum)
a10<-a10[-c(35),]

beta10<-betareg(props ~ tissue, data=a10)
summary(beta10)

d10<-emmeans(beta10,pairwise~tissue, type="response")
d10
CLD(d10$emmeans,  Letters ='ABCDEFGHIJKLMNOPQRS')


##plot
a10$tissue[a10$tissue=="Unripe pulp (2)"]="Unripe pulp"
a10$tissue[a10$tissue=="Flowers (4)"]="Flowers"
a10$tissue[a10$tissue=="Ripe pulp (1)"]="Ripe pulp"

tissueplot10<-ggplot(a10, aes(x=tissue, y=props)) +
	geom_boxplot() + geom_point(aes(color=as.character(plant)))+
	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
	theme_minimal()+
	scale_x_discrete(limits=c("Mature leaves","Flowers", "Unripe pulp","Ripe pulp",
							  "Seeds"))+
	scale_color_manual(values=as.vector(kelly(21)), name = "Plant ID")
tissueplot10

#tiff('tissueplot_select.tiff', units="in", width=6, height=4, res=500)
#tissueplot10
#dev.off()

#plot w/o color
tissueplot_bwj<-ggplot(a10, aes(x=tissue, y=props)) +
	geom_boxplot(outlier.shape = NA) + geom_jitter(position=position_jitter(width = 0.04), alpha=0.4)+
	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
	theme_classic()+
	scale_x_discrete(limits=c("Mature leaves","Flowers", "Unripe pulp","Ripe pulp",
							  "Seeds"))+
	stat_summary(geom = 'text', label = c("AB","A","BC","A","C"),
				 fun.y = max, vjust = -0.8)+
	scale_y_continuous(limits = c(0, 0.055), breaks = c(0,0.025,0.05))+
	theme(text = element_text(size=15),
		  axis.text.x = element_text(angle=45, hjust=1))
tissueplot_bwj

tiff('tissueplot_bwj.tiff', units="in", width=6, height=5, res=500)
tissueplot_bwj
dev.off()

postscript('tissueplot_bwj.svg', width=6, height=5)
tissueplot_bwj
dev.off()

tissueplot_lett<-ggplot(a10, aes(x=tissue, y=props)) +
	geom_boxplot(outlier.shape = NA) + geom_jitter(position=position_jitter(width = 0.04), alpha=0.4)+
	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
	theme_classic()+
	scale_x_discrete(limits=c("Mature leaves","Flowers", "Unripe pulp","Ripe pulp",
							  "Seeds"))+
	stat_summary(geom = 'text', label = c("AB","A","BC","A","C"),
				 fun.y = max, vjust = -0.8)+
	scale_y_continuous(expand = c(0, 0), limits = c(0, 0.055))

tissueplot_lett


#alkenylphenols over fruit ripening

#creating columns for stage as continuous variable and stage^2 (quadratic term)
dat[is.na(dat)] <- 0
dat$stage2 = (as.numeric(dat$stage))^2
dat$stage.num = as.numeric(dat$stage)

#aggregating dataset so each point=one plant and isn't sep by compound
ag_dat<-aggregate(props~plant+stage.num+stage2,data=dat,FUN=sum)
ag_dat$plant<-as.character(ag_dat$plant)
##removing weird chromatograms
ag_dat<-ag_dat[-c(2, 51, 52, 55,94,107),]

#betaregression with quadratic term
beta12<-betareg(props~stage.num+stage2 , data=ag_dat)
summary(beta12)#stage and quad. term are significant
#since quad. term is significant, suggests the data are not linear

#prediction plot
ag_dat$yhat_12<-predict(beta12)
predplot12<-ggplot(ag_dat)+
	geom_point(aes(x=stage.num, y=props))+
	geom_point(aes(x=stage.num, y=yhat_12), color="red", size=2)+
	geom_line(aes(x=stage.num, y=yhat_12) ,color="red", size=1)
predplot12

#regular betaregression
beta13<-betareg(props~stage.num, data=ag_dat)
summary(beta13)#stage is significant

#prediction plot
ag_dat$yhat_13<-predict(beta13)
predplot13<-ggplot(ag_dat)+
	geom_point(aes(x=stage.num, y=props))+
	geom_point(aes(x=stage.num, y=yhat_13), color="red", size=2)+
	geom_line(aes(x=stage.num, y=yhat_13) ,color="red", size=1)
predplot13

#prediction plot with all three models
predplot_all<-ggplot(ag_dat)+
	geom_point(aes(x=stage.num,y=props,color=plant))+
	geom_line(aes(x=stage.num, y=yhat_12),color="red",size=1)+
	geom_line(aes(x=stage.num, y=yhat_13),color="blue",size=1)
predplot_all

predplot_data<-ggplot(ag_dat)+
	geom_point(aes(x=stage.num,y=props,color=plant))+
	geom_line(aes(x=stage.num, y=yhat_12),color="black",size=1,linetype="dashed")+
	geom_line(aes(x=stage.num, y=yhat_13),color="black",size=1)+
	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
	theme_minimal()+
	scale_color_manual(values=as.vector(kelly(21)), name = "Plant ID")+
	scale_x_reverse(breaks=c(6,5,4,3,2,1))+
	scale_y_continuous(limits = c(0,.2))
predplot_data

##export graph
#tiff('predplot_data.tiff', units="in", width=8, height=5, res=500)
#predplot_data
#dev.off()

null.mod<-betareg(props~1, data=ag_dat)

modcomp<-aictab(cand.set=list(beta12,beta13,null.mod),
	   modnames=c("nonlinear","linear","null"))#AIC table
#top model is betareg with quad term
modcomp

stage_line_plot<-ggplot(ag_dat)+
	geom_point(aes(x=stage.num,y=props,color=plant))+
	stat_smooth(aes(x=stage.num,y=props),method = "lm", formula = y ~ x + I(x^2), size = 1,
				linetype="dashed", color="black")+
	stat_smooth(aes(x=stage.num,y=props),method = "lm", formula = y ~ x, size = 1, color="black")+
	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
	theme_minimal()+
	scale_color_manual(values=as.vector(kelly(21)), name = "Plant ID")+
	scale_x_reverse(breaks=c(6,5,4,3,2,1))+
	scale_y_continuous(expand=c(0,0), limits=c(-.01,2))+
	coord_cartesian(xlim=c(1,6), ylim=c(0,.2))
stage_line_plot

##EXPORT GRAPH
#tiff('stage_line_plot.tiff', units="in", width=7, height=4, res=500)
#stage_line_plot
#dev.off()

library(gridExtra)
library(grid)
DFtxt <- textGrob("Developing\nflowers", gp=gpar(fontsize=10, fontface="bold"))
FLtxt <- textGrob("Flowers", gp=gpar(fontsize=10, fontface="bold"))
UPtxt <- textGrob("Unripe pulp", gp=gpar(fontsize=10, fontface="bold"))
RPtxt<- textGrob("Ripe pulp", gp=gpar(fontsize=10, fontface="bold"))


stage_line_plot2<-ggplot(ag_dat)+
	geom_line(aes(x=stage.num,y=props,color=plant),alpha=0.5,size=1,show.legend = F)+
	stat_smooth(aes(x=stage.num,y=props),method = "lm", formula = y ~ x + I(x^2), size = 1.5,
				linetype="solid", color="black")+
	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
	theme_classic()+
	scale_color_viridis(discrete = T, option = "D")+
	scale_x_reverse(breaks=c(6,5,4,3,2,1))+
	theme(text = element_text(size = 15))+
	scale_y_continuous(expand=c(0,0))+
	coord_cartesian(xlim=c(1,6), ylim=c(0,.1))
stage_line_plot2

stage_line_plot2+
	theme(plot.margin = unit(c(1,1,2,1), "lines")) +
	annotation_custom(DFtxt,xmin=1,xmax=1,ymin=-0.07,ymax=-0.01) + 
	annotation_custom(FLtxt,xmin=2,xmax=3,ymin=-0.07,ymax=-0.01)+ 
	annotation_custom(UPtxt,xmin=4,xmax=5,ymin=-0.07,ymax=-0.01)+ 
	annotation_custom(RPtxt,xmin=5,xmax=6,ymin=-0.07,ymax=-0.01)+
	coord_cartesian(clip = "off")

stage_line_plot2+
	theme(plot.margin = unit(c(1,1,2,1), "lines"))+
	annotate(geom="text", label="Ripe\npulp", 
			 x=1, y=-0.005, size=4, color="black")+
	annotate(geom="text", label="Unripe pulp", 
			 x=2.5, y=-0.01, size=5, color="black")+
	annotate(geom="text", label="Flowers", 
			 x=4.5, y=-0.01, size=5, color="black")+
	annotate(geom="text", label="Developing\nflowers", 
			 x=6, y=-0.005, size=4, color="black")



tiff('stage_line_plot2.tiff', units="in", width=8, height=5, res=500)
stage_line_plot2
dev.off()

#leaf analysis
a20<-aggregate(props~plant+stage,data=dfl,FUN=sum)
a20$stage<-as.factor(a20$stage)
levels(a20$stage)
a20$prop1<-(a20$props+0.01)

shapiro.test(a20$props)#not norm dist

beta20<-betareg(prop1~stage, data=a20)
summary(beta20)

d20<-emmeans(beta20,pairwise~stage, type="response")
d20


#Mann-Whiteney U test b/c data is not normally distributed
#Mann-Whitney U test is a non-parametric test that can be used in place of an unpaired t-test
wilcox.test(a20$props~a20$stage)#W=20, p=0.016

a20$stagen[a20$stage=="E"]="Expanding leaf"
a20$stagen[a20$stage=="M"]="Mature leaf"

leaf_plot_bw<-ggplot(a20, aes(x=stagen,y=props))+
				  	geom_boxplot() + 
	geom_jitter(position=position_jitter(width = 0.04), alpha=0.5, size=2)+
				  	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
				  	scale_x_discrete(limits=c("Expanding leaf", "Mature leaf"))+
				  	theme_classic()+
	scale_y_continuous(limits = c(0,0.008))+theme(text = element_text(size=15))
				 
leaf_plot_bw

tiff('leaf_plot_bw.tiff', units="in", width=5, height=4, res=500)
leaf_plot_bw
dev.off()

##summary stats
library(plyr)
cdata <- ddply(datall, c("tissue", "compound"), summarise,
			   N    = length(props),
			   mean = mean(props),
			   sd   = sd(props),
			   se   = sd / sqrt(N))

cdata <- cdata[order(cdata$tissue),]
cdata[is.na(cdata)] <- 'nd'

write.table(cdata, file = "stats.csv", sep = ",", quote = FALSE, row.names = F)

##SUMMARY PLOT
datall <- datall[order(datall$compound),]
datall[is.na(datall)] <- 0
a.dat<-slice(datall, 1:126)
b.dat<-slice(datall, 127:252)
c.dat<-slice(datall, 253:378)
d.dat<-slice(datall, 379:504)
e.dat<-slice(datall, 505:630)
f.dat<-slice(datall, 631:756)
g.dat<-slice(datall, 757:882)
h.dat<-slice(datall, 883:1008)
i.dat<-slice(datall, 1009:1134)
j.dat<-slice(datall, 1135:1260)

a.plot<-ggplot(a.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds","Exp. leaves", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15))

a.plot

b.plot<-ggplot(b.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds","Exp. leaves", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15))

c.plot<-ggplot(c.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds","Exp. leaves", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15))

d.plot<-ggplot(d.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds","Exp. leaves", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15))

e.plot<-ggplot(e.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds","Exp. leaves", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15))

f.plot<-ggplot(f.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds","Exp. leaves", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15))

g.plot<-ggplot(g.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds","Exp. leaves", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15))

h.plot<-ggplot(h.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds","Exp. leaves", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15))

i.plot<-ggplot(i.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds","Exp. leaves", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15))

j.plot<-ggplot(j.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds","Exp. leaves", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15))

tiff('supp_fig.tiff', units="in", width=25, height=25, res=500)
ggarrange(a.plot, b.plot,c.plot,d.plot,e.plot,f.plot,g.plot,h.plot,i.plot,j.plot, 
		  labels = c("A", "B","C","D","E","F","G","H","I","J"),heights = c(20, 20),
		  ncol = 2, nrow = 5, align = "v",
		  font.label = list(size = 22))
dev.off()

tiff('supp_fig_ab.tiff', units="in", width=12, height=10, res=500)
ggarrange(a.plot, b.plot, 
		  labels = c("A", "B"),heights = c(5, 5),
		  ncol = 1, nrow = 2, align = "v",
		  font.label = list(size = 22))
dev.off()

tiff('supp_fig_cd.tiff', units="in", width=12, height=10, res=500)
ggarrange(c.plot, d.plot, 
		  labels = c("C", "D"),heights = c(5, 5),
		  ncol = 1, nrow = 2, align = "v",
		  font.label = list(size = 22))
dev.off()

tiff('supp_fig_ef.tiff', units="in", width=12, height=10, res=500)
ggarrange(e.plot, f.plot, 
		  labels = c("E", "F"),heights = c(5, 5),
		  ncol = 1, nrow = 2, align = "v",
		  font.label = list(size = 22))
dev.off()

tiff('supp_fig_gh.tiff', units="in", width=12, height=10, res=500)
ggarrange(g.plot, h.plot, 
		  labels = c("G", "H"),heights = c(5, 5),
		  ncol = 1, nrow = 2, align = "v",
		  font.label = list(size = 22))
dev.off()

tiff('supp_fig_ij.tiff', units="in", width=12, height=10, res=500)
ggarrange(i.plot, j.plot, 
		  labels = c("I", "J"),heights = c(5, 5),
		  ncol = 1, nrow = 2, align = "v",
		  font.label = list(size = 22))
dev.off()

supp.plot.all<-ggplot(datall, aes(x=tissue, y=props,fill=compound)) +
	geom_boxplot() + #geom_point()+
	labs(x=" ", y="Total alkenylphenols (proportion dry wt)")+
	theme_minimal()+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds","Exp. leaves", "Mature leaves"))+
	scale_fill_manual(values=as.vector(tol(10)), name = "Compound")+ 
	theme(text = element_text(size=10))
supp.plot.all	

tiff('supp.plot.all.tiff', units="in", width=10, height=4, res=500)
supp.plot.all
dev.off()

####Fungal bioassays####
datf <- read.csv(file="Maynard_etal_FungalBioassays.csv",head=TRUE,fill=T)
datf <- datf[1:81,]

datf$Conc<-as.numeric(datf$Conc)
datf$fungi<-as.character(datf$fungi)
datf$nfungi[datf$fungi=="R3"]="Microdochium lycopodium"
datf$nfungi[datf$fungi=="R23"]="Fusarium verticillioides"
datf$nfungi[datf$fungi=="R26"]="Fusarium sp."

datf <- datf[order(datf$fungi),]
R23<-slice(datf, 1:27)
R26<-slice(datf, 28:54)
R3<-slice(datf, 55:81)

mod23<-lm(data = R23, abs_corr~Conc)
mod26<-lm(data = R26, abs_corr~Conc)
mod3<-lm(data = R3, abs_corr~Conc)

summary(mod23)#p<0.0001
summary(mod26)#p=0.275
summary(mod3)#p<0.0001

R23$yhat<-predict(mod23)
predplot23<-ggplot(R23)+
	geom_smooth(aes(x=Conc, y=yhat), color="red")+
	geom_point(aes(x=Conc, y=abs_corr))
predplot23

R26$yhat<-predict(mod26)
predplot26<-ggplot(R26)+
	geom_smooth(aes(x=Conc, y=yhat), color="red")+
	geom_point(aes(x=Conc, y=abs_corr))
predplot26

R3$yhat<-predict(mod3)
predplot3<-ggplot(R3)+
	geom_smooth(aes(x=Conc, y=yhat), color="red")+
	geom_point(aes(x=Conc, y=abs_corr))
predplot3

leg_fung <- c(expression(paste(italic("Fusarium "), "sp.")),
			 expression(paste(italic("Fusarium verticillioides"))), 
			 expression(paste(italic("Microdochium lycopodium"))))


plota<-ggplot(datf, aes(x=Conc, y=abs_corr, group=nfungi))+
	geom_smooth(aes(color=nfungi), method = "lm", se=T)+
	geom_point(aes(color=nfungi))+
	theme_classic()+
	scale_color_viridis(discrete = T, option = "D", labels=leg_fung)+
	labs(x="Concentration (Prop. in ripe infruct.)", y="Absorbance (OD)", color=" ")+
	theme(legend.text.align = 0, text = element_text(size=18),legend.position="top")+ 
	annotate("text", x = 2.20, y = 0.66,
			 label = "paste(italic(R) ^ 2, \" = 0.84\")", parse = TRUE, size =5)+ 
	annotate("text", x = 2.20, y = 0.54,
			 label = "paste(italic(R) ^ 2, \" = 0.04\")", parse = TRUE, size =5)+ 
	annotate("text", x = 2.20, y = 0.415,
			 label = "paste(italic(R) ^ 2, \" = 0.77\")", parse = TRUE, size =5)
plota

tiff('fungiplot.tiff', units="in", width=8, height=5, res=500)
plota
dev.off()

##bw plot
plot.bw<-ggplot(datf, aes(x=Conc, y=abs_corr, group=nfungi))+
	geom_smooth(aes(linetype=nfungi),color="black", method = "lm", se=T)+
	scale_linetype_manual(values=c("dotted", "solid", "twodash"),labels=leg_fung)+
	geom_point()+
	theme_minimal()+
	labs(x="Concentration (Prop. in ripe infruct.)", y="Absorbance (OD)", linetype="Fungi")+
	theme(legend.text.align = 0)+ 
	annotate("text", x = 2.20, y = 0.68,
			 label = "paste(italic(R) ^ 2, \" = 0.84\")", parse = TRUE)+ 
	annotate("text", x = 2.20, y = 0.54,
			 label = "paste(italic(R) ^ 2, \" = 0.04\")", parse = TRUE)+ 
	annotate("text", x = 2.20, y = 0.42,
			 label = "paste(italic(R) ^ 2, \" = 0.77\")", parse = TRUE)
plot.bw

tiff('fungiplot_bw.tiff', units="in", width=8, height=5, res=500)
plot.bw
dev.off()


####Animal feeding trials####
animal <- read.csv(file="Maynard_etal_AlkenylphenolAnimalTrials.csv",head=TRUE)

#only animals that participated
animal<-animal%>%
	filter(participate==1)

#separating bats and birds 
bat<-slice(animal, 1:32)
bird<-slice(animal, 33:58)


birdg <- gather(bird, "treatment", "amount_eaten", 17:18)
birdg$amount_eaten<-as.numeric(birdg$amount_eaten)
birdg$amount_eaten[birdg$amount_eaten>3] <- 3
birdg$treatment<-as.character(birdg$treatment)
bird_ag<-aggregate(amount_eaten~ID+treatment, data=birdg, FUN=mean) 
bird_1<-spread(data=bird_ag, treatment, amount_eaten)

shapiro.test(bird_1$c.eaten)#norm dist
shapiro.test(bird_1$t.eaten)#norm dist

batg <- gather(bat, "treatment", "amount_eaten", 17:18)
batg$amount_eaten<-as.numeric(batg$amount_eaten)
batg$treatment<-as.character(batg$treatment)
bat_ag<-aggregate(amount_eaten~ID+treatment, data=batg, FUN=mean) 
bat_1<-spread(data=bat_ag, treatment, amount_eaten)

shapiro.test(bat_1$c.eaten)#not norm dist?
shapiro.test(bat_1$t.eaten)#norm dist

t.test(bird_1$c.eaten, bird_1$t.eaten, paired = T)#t=0.30,df=9,p=0.775
t.test(bat_1$c.eaten, bat_1$t.eaten, paired = T)#t=2.73,df=16,p=0.015

#change names for graph
bat_ag$treatment[bat_ag$treatment=="t.eaten"]="Treatment"
bat_ag$treatment[bat_ag$treatment=="c.eaten"]="Control"
bird_ag$treatment[bird_ag$treatment=="t.eaten"]="Treatment"
bird_ag$treatment[bird_ag$treatment=="c.eaten"]="Control"

batpref<-ggplot(bat_ag,aes(x=treatment,y=amount_eaten))+geom_boxplot()+
	geom_jitter(position=position_jitter(width = 0.03), alpha=0.5, size=1.7)+
	theme_classic()+
	labs(x=" ",y="Average amount eaten (g)")+
	scale_y_continuous(limits =  c(0,3))+
	theme(text = element_text(size = 18))
batpref

birdpref<-ggplot(bird_ag,aes(x=treatment,y=amount_eaten))+geom_boxplot()+
	geom_jitter(position=position_jitter(width = 0.03), alpha=0.5,size=1.7)+
	theme_classic()+
	labs(x=" ",y="")+
	scale_y_continuous(limits =  c(0,3.00))+
	theme(text = element_text(size = 18))
birdpref

tiff('animal_pref.tiff', units="in", width=8, height=4, res=500)
ggarrange(batpref, birdpref, 
		  labels = c("a", "b"),heights = c(2, 2),
		  ncol = 2, nrow = 1, align = "v")
dev.off()

tiff('animal_pref1.tiff', units="in", width=4, height=8, res=500)
ggarrange(batpref, birdpref1, 
		  labels = c("A", "B"),heights = c(5, 5),
		  ncol = 1, nrow = 2, align = "v")
dev.off()

####Removal study####
remo <- read.csv(file="Maynard_etal_PiperRemovalStudy.csv",head=TRUE,fill=T)

remo<-remo[-c(71),]

remo$Ripeness <- NA
for(i in 1:length(remo$Ripeness)){
	if(remo$ripeness[i]=="1"){remo$Ripeness[i]="Ripe"}
	if(remo$ripeness[i]=="2"){remo$Ripeness[i]="Unripe"}
}

remo$Time <- NA
for(i in 1:length(remo$Time)){
	if(remo$time.of.day[i]=="1"){remo$Time[i]="Nocturnal"}
	if(remo$time.of.day[i]=="2"){remo$Time[i]="Diurnal"}
}

tbl = table(remo$Ripeness, remo$Time) 
tbl
chisq.test(tbl)

tbl.1 = table(remo$Time, remo$Ripeness)
tbl.1
chisq.test(tbl.1)

tbl2<-as.data.frame.matrix(tbl)
tbl2

tbl3<-as.data.frame(tbl)
tbl3

bp <- ggplot(tbl3, aes(x=Var2, y=Freq, fill=Var1))+
	geom_bar(aes(fill=Var1),stat = "identity", position = "stack")+
	theme_classic()+
	labs(x="", y="No. of fruits removed",fill="")+
	theme(text = element_text(size = 18))+
	scale_fill_manual(values=c('#006837','#addd8e'))
bp

