rm(list = ls())
library(VennDiagram)  
library(UpSetR)
library(readxl)
library(ggsci)
library(ggplot2)
library(sf)
library(ggVennDiagram)


data = readxl::read_excel('./sampleoverlap.xlsx')
data = data.frame(data)

cohorts = unique(c(data$EducationalAttainment, data$CUD,
              data$LCU, data$Intelligence, 
              data$ADHD, data$cigarette.smoking))

data_new = data.frame(cohorts=cohorts)
data_new$EduAttain = ifelse(data_new$cohorts %in% data$EducationalAttainment,1,0)
data_new$CUD = ifelse(data_new$cohorts %in% data$CUD,1,0)
data_new$LCU = ifelse(data_new$cohorts %in% data$LCU,1,0)
data_new$Intelligence = ifelse(data_new$cohorts %in% data$Intelligence,1,0)
data_new$ADHD = ifelse(data_new$cohorts %in% data$ADHD,1,0)
data_new$CigSmk = ifelse(data_new$cohorts %in% data$cigarette.smoking,1,0)

UpSetR::upset(data_new, 
      sets = c("EduAttain", "CUD", "LCU", "Intelligence", "ADHD", "CigSmk"), 
      order.by = c("degree", "freq"),
      decreasing = c(TRUE, FALSE), 
      sets.bar.color = ggpubr::get_palette('npg',dim(data_new)[2]-1),                  
      matrix.color = 'black', 
      main.bar.color = 'black', 
      mainbar.y.label = 'Intersection Size',
      sets.x.label = 'Number of cohorts', 
      point.size = 1.5, 
      line.size = 0.5, 
      att.pos = 'bottom',
      set_size.show = T, 
      att.color = 'yellowgreen', 
      number.angles = 0, 
      # empty.intersections = 'on'  
      queries = list(
        list(
          query = intersects, 
          params = list("EduAttain", "CigSmk", "CUD", "LCU"), 
          color = "orange", 
          active = T), 
        list(
          query = intersects, 
          params = list("EduAttain", "CigSmk", "Intelligence"), 
          color = "#dda898", 
          active = T), 
        list(
          query = intersects,
          params = list("EduAttain", "CigSmk", "CUD"), 
          color = "#3fbdbc",
          active = T),
        list(
          query = intersects,
          params = list("EduAttain", "CigSmk"), 
          color = "brown",
          active = T),
        list(
          query = intersects,
          params = list("EduAttain", "LCU"), 
          color = "#a278d2",
          active = T),
        list(
          query = intersects,
          params = list("EduAttain", "ADHD"), 
          color = "#b53c92",
          active = T),
        list(
          query = intersects,
          params = list("EduAttain", "CigSmk", "LCU"), 
          color = "#297c5c",
          active = T)
        )
      )

savePlot(filename = 'Rplot', type = 'png')


venn.diagram(x = list(A = data$EducationalAttainment,B = data$CUD),
             height = 3000, 
             width = 3000, 
             resolution = 500,
             main = 'test', 
             main.cex = 4, 
             main.col = 'red', 
             fill = c('red','blue'), 
             col = 'black',
             lwd = 2, 
             alpha = .5, 
             cat.cex = 2, 
             filename = './venn.test.tiff'
)


color1 <- alpha("#f8766d",0.9)
color2 <- alpha("#FF99CC",0.7)
color3 <- alpha("#c77cff",0.5)
color4 <- alpha("#99CC00",0.5)
color5 <- alpha("#5b39ac",0.7)
ggVennDiagram(data[c(1,6)], label_alpha=0) +
scale_fill_gradient(low="white", high =color3, guide="none")

