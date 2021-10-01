library(UpSetR)
input <- Transfac_DC2_June2020

upset(input, 
      order.by = c("degree"),
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1
      )

