library(ggplot2)
library(readxl)
library(cowplot)

d1 = read_excel('./GWAS/PowerCalculation/PowCal1.xlsx')
d1 = data.frame(d1)
d2 = read_excel('./GWAS/PowerCalculation/PowCal2.xlsx')
d2 = data.frame(d2)


g1 = ggplot(d1, mapping=aes(x=OR,y=Power,group=Outcome))+
    geom_line(aes(color=Outcome, size=Outcome)) + 
    geom_point(aes(color=Outcome)) +
    theme_classic() +  
    geom_hline(yintercept = 0.8, size=1, lty=2, color='#4836a2') +
    theme(legend.position="top") + 
    scale_color_manual(values=c('#c1466d','#349c4d'))+  
    scale_size_manual(values=c(0.5, 0.5)) + 
    scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 0.8, 1.00)) + 
    labs(title='Educational attainment on CUD or LCU',
         x="Odds ratio per standard deviation increase in educational attainment",
         y="Power")

g2 = ggplot(d2, mapping=aes(x=BETA,y=Power,group=Exposure))+
  geom_line(aes(color=Exposure, size=Exposure)) + 
  geom_point(aes(color=Exposure)) +
  theme_classic() +  
  geom_hline(yintercept = 0.8, size=1, lty=2, color='#4836a2') +
  theme(legend.position="top") + 
  scale_color_manual(values=c('#c1466d','#349c4d'))+  
  scale_size_manual(values=c(0.5, 0.5)) + 
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 0.8, 1.00)) + 
  labs(title='CUD or LCU on educational attainment',
       x="Effect size per log odds ratio change in CUD or LCU",
       y="Power")

# for 2016EA to LCU
new_plot = plot_grid(g1, g2,
                     labels = "AUTO", ncol = 2,
                     align = 'h', axis = "l",
                     label_size = 12)
save_plot("./GWAS/PowerCalculation/power_plot.jpg",
          new_plot,
          ncol = 2, 
          nrow = 2,
          base_height = 3.71,
          dpi = 800)










