library(UpSetR)
input <- read.csv("https://raw.githubusercontent.com/SameOldSamOld/RNA-Seq_OG4953/master/data/REACTOME_COMPARISON_LIST.csv",
                  row.names = 1)

upset(input,
      order.by = c("degree"),
      mb.ratio = c(0.6, 0.4),
      number.angles = 0,
      text.scale = 1.1,
      point.size = 2.8,
      line.size = 1
      )

