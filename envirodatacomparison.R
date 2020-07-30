### Comparison of environmental data ###
##### GITHUB VERSION ######


library("car")
library("ggplot2")
library("plyr")
library("reshape")
library("tidyr")
library("dplyr")
library("cowplot")
library("agricolae")

#keeping everything in this folder
setwd("/Users//tdolan/Documents//R-Github/WFCageStudy")

#inspect the data.
sev<-read.csv("vert2_2017.csv")
six<-read.csv("vert2_2016.csv")

sev <-mutate(sev, Year="2017")
six <-mutate(six, Year="2016")
vert <- bind_rows(six,sev)


#the geom_flat violin function#
##### The geom flat violin function ####

#https://github.com/tidyverse/ggplot2/issues/2459

library(tidyverse)
#devtools::install_github(repo = "IndrajeetPatil/ggstatsplot")
library(ggstatsplot)
library(ggplot2)
library(dplyr)


"%||%" <- function(a, b) {
  if (!is.null(a))
    a
  else
    b
}

#=========================== function definition ===========================

geom_flat_violin <-
  function(mapping = NULL,
           data = NULL,
           stat = "ydensity",
           position = "dodge",
           trim = TRUE,
           scale = "area",
           show.legend = NA,
           inherit.aes = TRUE,
           ...) {
    ggplot2::layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomFlatViolin,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(trim = trim,
                    scale = scale,
                    ...)
    )
  }
#@rdname ggplot2-ggproto
#@Format NULL
#@Usage NULL
#@export

GeomFlatViolin <-
  ggproto(
    "GeomFlatViolin",
    Geom,
    setup_data = function(data, params) {
      data$width <- data$width %||%
        params$width %||% (resolution(data$x, FALSE) * 0.9)
      
      # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
      data %>%
        dplyr::group_by(.data = ., group) %>%
        dplyr::mutate(
          .data = .,
          ymin = min(y),
          ymax = max(y),
          xmin = x,
          xmax = x + width / 2
        )
    },
    
    draw_group = function(data, panel_scales, coord)
    {
      # Find the points for the line to go all the way around
      data <- base::transform(data,
                              xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
      
      # Make sure it's sorted properly to draw the outline
      newdata <-
        base::rbind(
          dplyr::arrange(.data = base::transform(data, x = xminv), y),
          dplyr::arrange(.data = base::transform(data, x = xmaxv), -y)
        )
      
      # Close the polygon: set first and last point the same
      # Needed for coord_polar and such
      newdata <- rbind(newdata, newdata[1,])
      
      ggplot2:::ggname("geom_flat_violin",
                       GeomPolygon$draw_panel(newdata, panel_scales, coord))
    },
    
    draw_key = draw_key_polygon,
    
    default_aes = ggplot2::aes(
      weight = 1,
      colour = "grey20",
      fill = "white",
      size = 0.5,
      alpha = NA,
      linetype = "solid"
    ),
    
    required_aes = c("x", "y")
  )
######


#### Temperature #####
# raincloud boxplot Temperature
#+coord_flip() to change the orientation
#https://peerj.com/preprints/27137v1.pdf
#2016
p6 <- ggplot(six,aes(x=station,y=temp, fill = site, colour = site))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust
                   =2, trim = FALSE)+
  geom_point(position = position_jitter(width = .15), size = .50)+
  geom_boxplot(aes(x = station, y = temp),
               outlier.shape=NA, alpha = 0.2, width = .1, colour = "BLACK") +
  ylab('Temperature C')+xlab('station')+theme_cowplot()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_brewer(palette = "Greys")+
  scale_fill_brewer(palette = "Greys")
  #facet_wrap(~Year)+
  #ggtitle("Temperature (C)")
p6
#ggsave('temp_violins6.png', width = 7, height = 7, path = "/Users//tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()

#2017
p6 <- ggplot(sev,aes(x=station,y=temp, fill = site, colour = site))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust
                   =2, trim = FALSE)+
  geom_point(position = position_jitter(width = .15), size = .50)+
  geom_boxplot(aes(x = station, y = temp),
               outlier.shape=NA, alpha = 0.2, width = .1, colour = "BLACK") +
  ylab('Temperature C')+xlab('station')+theme_cowplot()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_brewer(palette = "Greys")+
  scale_fill_brewer(palette = "Greys")
#facet_wrap(~Year)+
#ggtitle("Temperature (C)")
p6
#ggsave('temp_violins7.png', width = 7, height = 7, path = "/Users//tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()


#Temperature 2016 anova. 
anovaT16 <-lm(temp~station, na.action=na.omit, data=six)
ano2 <-car::Anova(anovaT16)
ano2
df<-df.residual(anovaT16)
MSerror<-deviance(anovaT16)/df
comparison <- HSD.test(anovaT16,c("station"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05, group=TRUE)
comparison
## I think this means they are all different. 

#Temperature Site and Depth 
anovaT16 <-lm(temp~site*depth, na.action=na.omit, data=six)
ano2 <-car::Anova(anovaT16)
ano2
df<-df.residual(anovaT16)
MSerror<-deviance(anovaT16)/df
comparison <- HSD.test(anovaT16,c("site","depth"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05/15, group=TRUE)
comparison
## I think this means they are all different. 

#Temperature 2017 anova. 
anovaT17 <-lm(temp~station, na.action=na.omit, data=sev)
ano2 <-car::Anova(anovaT17)
ano2
df<-df.residual(anovaT17)
MSerror<-deviance(anovaT17)/df
comparison <- HSD.test(anovaT17,c("station"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05/15, group=TRUE)
comparison
## I think this means they are all different. 

#Temperature Site and Depth 
anovaT17 <-lm(temp~site*depth, na.action=na.omit, data=six)
ano2 <-car::Anova(anovaT17)
ano2
df<-df.residual(anovaT17)
MSerror<-deviance(anovaT17)/df
comparison <- HSD.test(anovaT17,c("site","depth"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05, group=TRUE)
comparison
## I think this means they are all different. 
#######

#### Dissolved Oxygen ####

#2016
p6 <- ggplot(six,aes(x=station,y=DO, fill = site, colour = site))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust
                   =2, trim = FALSE)+
  geom_point(position = position_jitter(width = .15), size = .50)+
  geom_boxplot(aes(x = station, y = DO, color=site),
               outlier.shape=NA, alpha = 0.2, width = .1, colour = "BLACK") +
  ylab('Dissolved Oxygen (mg/L)')+xlab('station')+theme_cowplot()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_brewer(palette = "Greys")+
  scale_fill_brewer(palette = "Greys")
#facet_wrap(~Year)+
#ggtitle("Temperature (C)")
p6
#ggsave('DO_violins6.png', width = 7, height = 7, path = "/Users//tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()

#2017
p6 <- ggplot(sev,aes(x=station,y=DO, fill = site, colour = site))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust
                   =2, trim = FALSE)+
  geom_point(position = position_jitter(width = .15), size = .50)+
  geom_boxplot(aes(x = station, y = DO),
               outlier.shape= NA, alpha = 0.3, width = .1, colour = "BLACK") +
  ylab('Dissolved Oxygen (mg/L)')+xlab('station')+theme_cowplot()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_brewer(palette = "Greys")+
  scale_fill_brewer(palette = "Greys")
#facet_wrap(~Year)+
#ggtitle("Temperature (C)")
p6
#ggsave('DO_violins7.png', width = 7, height = 7, path = "/Users//tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()


#DO 2016 anova. 
anovaT16 <-lm(DO~station, na.action=na.omit, data=six)
ano2 <-car::Anova(anovaT16)
ano2
df<-df.residual(anovaT16)
MSerror<-deviance(anovaT16)/df
comparison <- HSD.test(anovaT16,c("station"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05/15, group=TRUE)
comparison
## I think this means they are all different. 

#DO Site and Depth 
anovaT16 <-lm(DO~site*depth, na.action=na.omit, data=six)
ano2 <-car::Anova(anovaT16)
ano2
df<-df.residual(anovaT16)
MSerror<-deviance(anovaT16)/df
comparison <- HSD.test(anovaT16,c("site","depth"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05/15, group=TRUE)
comparison
## I think this means they are all different. 

#DO 2017 anova. 
anovaT17 <-lm(DO~station, na.action=na.omit, data=sev)
ano2 <-car::Anova(anovaT17)
ano2
df<-df.residual(anovaT17)
MSerror<-deviance(anovaT17)/df
comparison <- HSD.test(anovaT17,c("station"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05/15, group=TRUE)
comparison


#DO Site and Depth 
anovaT17 <-lm(DO~site*depth, na.action=na.omit, data=six)
ano2 <-car::Anova(anovaT17)
ano2
df<-df.residual(anovaT17)
MSerror<-deviance(anovaT17)/df
comparison <- HSD.test(anovaT17,c("site","depth"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05/15, group=TRUE)
comparison

###########

###### Salinity ########
#2016
p6 <- ggplot(six,aes(x=station,y=sal, fill = site, colour = site))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust
                   =2, trim = FALSE)+
  geom_point(position = position_jitter(width = .15), size = .50)+
  geom_boxplot(aes(x = station, y = sal, color=site),
               outlier.shape=NA, alpha = 0.2, width = .1, colour = "BLACK") +
  ylab('Salinity (ppt)')+xlab('station')+theme_cowplot()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_brewer(palette = "Greys")+
  scale_fill_brewer(palette = "Greys")
#facet_wrap(~Year)+
#ggtitle("salerature (C)")
p6
#ggsave('sal_violins6.png', width = 7, height = 7, path = "/Users//tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()

#2017
p6 <- ggplot(sev,aes(x=station,y=sal, fill = site, colour = site))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust
                   =2, trim = FALSE)+
  geom_point(position = position_jitter(width = .15), size = .50)+
  geom_boxplot(aes(x = station, y = sal),
               outlier.shape= NA, alpha = 0.3, width = .1, colour = "BLACK") +
  ylab('Salinity (ppt)')+xlab('station')+theme_cowplot()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_brewer(palette = "Greys")+
  scale_fill_brewer(palette = "Greys")
#facet_wrap(~Year)+
#ggtitle("salerature (C)")
p6
#ggsave('sal_violins7.png', width = 7, height = 7, path = "/Users//tdolan/Documents/WIP research/Caging paper/caging manuscript/cage_figs")
#dev.off()


#sal 2016 anova. 
anovaT16 <-lm(sal~station, na.action=na.omit, data=six)
ano2 <-car::Anova(anovaT16)
ano2
df<-df.residual(anovaT16)
MSerror<-deviance(anovaT16)/df
comparison <- HSD.test(anovaT16,c("station"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05/15, group=TRUE)
comparison
## I think this means they are all different. 

#sal Site and Depth 
anovaT16 <-lm(sal~site*depth, na.action=na.omit, data=six)
ano2 <-car::Anova(anovaT16)
ano2
df<-df.residual(anovaT16)
MSerror<-deviance(anovaT16)/df
comparison <- HSD.test(anovaT16,c("site","depth"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05/15, group=TRUE)
comparison

#sal 2017 anova. 
anovaT17 <-lm(sal~station, na.action=na.omit, data=sev)
ano2 <-car::Anova(anovaT17)
ano2
df<-df.residual(anovaT17)
MSerror<-deviance(anovaT17)/df
comparison <- HSD.test(anovaT17,c("station"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05/15, group=TRUE)
comparison
 

#sal Site and Depth 
anovaT17 <-lm(sal~site*depth, na.action=na.omit, data=six)
ano2 <-car::Anova(anovaT17)
ano2
df<-df.residual(anovaT17)
MSerror<-deviance(anovaT17)/df
comparison <- HSD.test(anovaT17,c("site","depth"),MSerror=MSerror, unbalanced=TRUE,alpha=0.05/15, group=TRUE)
comparison