library('survminer')
library('survival')
library(ggplot2)
library(gridExtra)

cal_KM <- function(item, legend_labels, titleitem){
  
  custom_theme <- function() {
    theme_survminer() %+replace%
      theme(
        plot.title=element_text(hjust=0.5)
      )
  }

  p <- ggsurvplot(fit, data=sub_colon,
                  risk.table = TRUE,
                  pval = TRUE,
                  legend.title = '',
                  legend.labs = c('greater than 3.3%', 'lower than 3.3%'),
                  title = paste('Progression Free Interval',title_item, sep='\n'),
                  font.title = c(20, 'bold'),
                  font.x = 17,
                  font.y = 17,
                  font.tickslab=15,
                  font.legend = 15,
                  xlab = 'Time(months)',
                  xlim = c(0,150),
                  break.x.by = 50,
                  lwd = 1,
                  ggtheme=custom_theme()

  )
  surv_plot <- p$plot
  table_plot <- p$table
  combined <- arrangeGrob(surv_plot, table_plot, ncol=1, heights=c(7,3))
  
  file_path <- paste('/Users/yuw/Desktop/pfi-',item,'.png', sep = "")
  ggsave(file_path, width=5, height=5, plot=combined,dpi=300)
}



plotlist <- list('til_status',
                 'hajcc_til', 'lajcc_til',
                 'lsize_til', 'hsize_til', 'lscore_til', 'hscore_til', 'ltumor_til', 'htumor_til')
titlelist <- list('UKentucky: TILs%',
                  'Stages III-IV: TILs%', 'Stages I-II: TILs%',
                  'Tumor size ≤ 5cm: TILs%', 'Tumor size > 5cm: TILs%',
                  'G1/G2: TILs%', 'G3: TILs%', 'Tumor% ≤ 30.0%: TILs%', 'Tumor% > 30.0%: TILs%')


plotlist <- list('til_status', 'hajcc_til', 'lajcc_til',
                 'lsize_til', 'hsize_til', 'lmsi_til', 'hmsi_til', 'htumor_til', 'ltumor_til')
titlelist <- list('TCGA READ: TILs%',
                  'Stages III-IV: TILs%', 'Stages I-II: TILs%',
                  'Tumor size ≤ 4.9cm: TILs%', 'Tumor size > 4.9cm: TILs%', 'MSS/MSI-L: TILs%', 'MSI-H: TILs%',
                  'Tumor% > 69.1%: TILs%', 'Tumor% ≤ 69.1%: TILs%')

plotlist <- list('til_status', 'hajcc_til', 'lajcc_til',
                 'Basal_til', 'LumA_til', 'LumB_til', 'Her_til',
                 'pER_til', 'nER_til', 'htumor_til', 'ltumor_til')
titlelist <- list('UNC: TILs%','Stages III-IV: TILs%', 'Stages I-II: TILs%',
                  'Basal: TILs%', 'LumA: TILs%', 'LumB: TILs%', 'Her2: TILs%',
                  'ER Positive: TILs%', 'ER Negative: TILs%', 'Tumor% >15.68%: TILs%', 'Tumor% ≤ 15.68%: TILs%')


for (i in seq_along(plotlist)){
  item = unlist(plotlist[i])
  title_item = unlist(titlelist[i])
  print(item)
  #colon <-  read.csv('/Users/yuw/Documents/Development/Colon_Data/data/coad+uk.csv')
  colon <- read.csv('/Users/yuw/Documents/Development/UNC/UNC-Breast/organized/incep.csv')
  #legend_labels <- unique(sub_colon[[item]])
  sub_colon <- subset(colon, colon[[item]] %in% c('high', 'low'))
  #sub_colon$OS.month <- as.numeric(sub_colon$OS.month)
  #sub_colon <- subset(colon, colon[[item]]%in% c('Q1', 'Q2', 'Q3', 'Q4'))
  fit <- survfit(Surv(OS.month, OS) ~ sub_colon[[item]], data=sub_colon)
  legend_labels <- unique(sub_colon[[item]])
  cal_KM(item, legend_labels, title_item)
}