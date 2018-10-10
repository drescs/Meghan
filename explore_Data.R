library(readr)
library(stringr)
library(reshape2)
library(ggpubr)
mypath="C:/Users/Sara/Documents/HHMI_projects/Meghan_Garrett/2018_10_05"
counts=read.table("C:/Users/Sara/Documents/HHMI_projects/Meghan_Garrett/2018_10_05/2018.10.05.annotatedCounts.txt", header=T, stringsAsFactors = F)

# "X4E10.2"              "Beads.only.Lib1.1"   
# [19] "Beads.only.Lib1.2"    "Beads.only.Lib2.1"    "Beads.only.Lib2.2"    "Beads.only.Lib2a"     "Beads.only.Lib2b"     "Input.Lib1"          
# [25] "Input.Lib2a"          "Input.Lib2b"          "Input.Lib2"          
# 


myxtitle="Input.Lib1"
myytitle="Input.Lib2"
print(ggscatter(counts, x = myxtitle, y = myytitle, 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = myxtitle, ylab = myytitle))

myxtitle="Beads.only.Lib1.1"
myytitle="Beads.only.Lib1.2"
print(ggscatter(counts, x = myxtitle, y = myytitle, 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = myxtitle, ylab = myytitle))




myxtitle="Beads.only.Lib1.1"
myytitle="Beads.only.Lib1.2"
temp_counts=counts[-which(counts[,myxtitle]==max(counts[,myxtitle])),]
print(ggscatter(temp_counts, x = myxtitle, y = myytitle, 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = myxtitle, ylab = myytitle))


myxtitle="Input.Lib2a"
myytitle="Input.Lib2b"
temp_counts=counts[-which(counts[,myxtitle]==max(counts[,myxtitle])),]
print(ggscatter(temp_counts, x = myxtitle, y = myytitle, 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = myxtitle, ylab = myytitle))


