#Step 1 is to import an annotated counts file:
#install.packages(c("readr", "stringr", "reshape2","ggpubr", "olsrr"))
library(readr)
library(stringr)
library(reshape2)
library(ggpubr)
library(olsrr)
mypath="C:/Users/Sara/Documents/HHMI_projects/Meghan_Garrett/2018_10_05"
counts=read.table("C:/Users/Sara/Documents/HHMI_projects/Meghan_Garrett/2018_10_05/2018.10.05.annotatedCounts.txt", header=T, stringsAsFactors = F)

#First, check that there is good correlation between replicates if two replicates exist.
#need to determine which reps are which; this is contained in the column names
#Multiple formats exist for column names, some don't include all parts

#pull out relevant parts of names here
sepnames=str_split_fixed(colnames(counts), fixed("."), 2)
colnames(sepnames)=c("Antibody", "rep")
rep1=which(sepnames[,2]%in%c("1"))
rep2=which(sepnames[,2]%in%c("2"))
replicate=rep(NA, length(colnames(counts)))
replicate[rep1]=1
replicate[rep2]=2

#Generate correlation plots for replicates--all in one PDF
correlations=rep(NA, length(rep1))
names(correlations)=sepnames[rep1,1]
pdf(file = "C:/Users/Sara/Documents/HHMI_projects/Meghan_Garrett/2018_10_05/corrplots_10_05.pdf", title="Replicate correlation plots")
for (i in 1:length(rep1)){
  myxtitle=colnames(counts)[rep1[i]]
  myytitle=colnames(counts)[rep2[i]]
  print(ggscatter(counts, x = myxtitle, y = myytitle, 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = myxtitle, ylab = myytitle))
  correlations[i]=cor(counts[,myxtitle], counts[,myytitle])
  
  model <- lm(myxtitle ~ myytitle, data = counts)
  ols_plotcooksd_bar(model)
}

#Do the same for input! Not included in forloop as name structure differs.
myxtitle="Input.Lib1"
myytitle="Input.Lib2"
print(ggscatter(counts, x = myxtitle, y = myytitle, 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = myxtitle, ylab = myytitle))



dev.off() #this command stops writing to PDF

#Now we want to remove the antibodies/plasma samples where the two replicates don't correlate well
#We chose a threshhold of 0.8 for the Pearson coefficient afer visual inspection


#Another troubleshooting step will involve checking if outliers unduly influence our correlation.
#To be efficient, here we'll trim off the max value. 
#This is quick and dirty. Would be better to remove, for example, all values greater than some # of SDs above mean, or some other metric
#Arbitrarily, I am removing the max value of replicate 1
#--we'll then inspect the plots to see if this misses any important outliers


#Generate correlation plots for replicates--all in one PDF
correlations2=rep(NA, length(rep1))
names(correlations2)=sepnames[rep1,1]
pdf(file = "C:/Users/Sara/Documents/HHMI_projects/Meghan_Garrett/2018_10_05/corrplots_no_outlier_10_05.pdf", title="Replicate correlation plots no outlier")
for (i in 1:length(rep1)){
  myxtitle=colnames(counts)[rep1[i]]
  myytitle=colnames(counts)[rep2[i]]
  temp_counts=counts[-which(counts[,myxtitle]==max(counts[,myxtitle])),]
  print(ggscatter(temp_counts, x = myxtitle, y = myytitle, 
                  add = "reg.line", conf.int = TRUE, 
                  cor.coef = TRUE, cor.method = "pearson",
                  xlab = myxtitle, ylab = myytitle))
  
  correlations2[i]=cor(counts[,myxtitle], counts[,myytitle])
}

#Do the same for input! Not included in forloop as name structure differs.
myxtitle="Input.Lib1"
myytitle="Input.Lib2"
print(ggscatter(counts, x = myxtitle, y = myytitle, 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = myxtitle, ylab = myytitle))



dev.off() #this command stops writing to PDF




samples_to_drop=c(colnames(counts)[rep1[which(correlations<0.8)]],colnames(counts)[rep2[which(correlations<0.8)]])

# samples_to_drop
# [1] "X14639.20ug.1" "X2988.10ug.1"  "X2988.20ug.1"  "X2C6.10ng.1"   "X6305.10ug.1"  "X6305.20ug.1"  "X6F5.10ng.1"   "X9577.10ug.1" 
# [9] "X9577.20ug.1"  "beads.only.1"  "PGT145.10ng.1" "QA255.10ug.1"  "QA255.20ug.1"  "X14639.20ug.2" "X2988.10ug.2"  "X2988.20ug.2" 
# [17] "X2C6.10ng.2"   "X6305.10ug.2"  "X6305.20ug.2"  "X6F5.10ng.2"   "X9577.10ug.2"  "X9577.20ug.2"  "beads.only.2"  "PGT145.10ng.2"
# [25] "QA255.10ug.2"  "QA255.20ug.2" 

#Don't actually want to drop both bead counts--we know which of these didn't work. Need to keep beads.only.2.

#samples_to_drop=samples_to_drop[-which(samples_to_drop=="beads.only.2")]
#
#counts[,samples_to_drop]=NULL

#I'll use this in the next file, "search for epitopes"
write.csv(counts, file="C:/Users/Sara/Documents/HHMI_projects/Meghan_Garrett/2018_10_05/only_clean_replicates.csv")

#Next we'll normalize the data--subtract the bead counts and input and divide by total reads



#First generate the sum of all the columns
#First several columns are metadata
head(colnames(counts))
# [1] "id"                   "Virus_Protein_Loc_AA" "Oligo"                "X006.10ng.1"          "X006.10ng.2"         
# [6] "X016.10ng.1"      

#sum all the columns except first three
col_sums=apply((counts[,-c(1:3)]),2, sum)

#divide each count by the correct column sum
temp_counts=cbind(counts[,1:3],(sweep((counts[,-c(1:3)]), 2, col_sums, "/")))


#some oligos had 0 input reads detected. To avoid dividing by zero, I'll set these equal to the mean of the input for all of the oligos
#How many had 0?

temp_counts$Input.Lib1[which(temp_counts$Input.Lib1==0)]=mean(temp_counts$Input.Lib1)
temp_counts$mean_beads=apply(cbind(temp_counts$Beadsonly.Lib1.1 ,temp_counts$Beadsonly.Lib1.2),1, mean)
write.csv(temp_counts,"C:/Users/Sara/Documents/HHMI_projects/Meghan_Garrett/2018_10_05/proportions.csv")

#divided the normalized output by the normalized input counts
temp_counts2=cbind(counts[,1:3],(sweep(temp_counts[,-c(1:3)], 1, temp_counts$Input.Lib1, "/")))

#now subtract the beads only controls
#counts$beads.only.1=NULL #this beads attempt failed
#leaving for now to not mess up our indices

temp_counts2$mean_beads=apply(cbind(temp_counts2$Beadsonly.Lib1.1 ,temp_counts$Beadsonly.Lib1.2),1, mean)
standardized_enrichment=cbind(counts[,1:3],(sweep(temp_counts2[,-c(1:3)], 1, (temp_counts2$mean_beads), FUN = "-", check.margin = TRUE)))


#Save this file:
write.csv(standardized_enrichment,"C:/Users/Sara/Documents/HHMI_projects/Meghan_Garrett/2018.10.05.standardized_enrichment.csv")

#time to clean the environment

rm(temp_counts,temp_counts2,col_sums, myxtitle, myytitle,i, samples_to_drop)



#Next, we'll make datasets for each strain
#However, some oligomers are commnon between two strains



#find the oligomers common to two strains
#if there are two strains the names are separated by "."

find_doubles=str_split(counts$Virus_Protein_Loc_AA, fixed("."))
double_positions=vector()
for (i in 1:length(find_doubles)){
  if (length(find_doubles[[i]])>1){
    double_positions=c(double_positions,i)
  }
}


#double_positions gives the rows where one oligomer is the same between two strains
#we want to make this portion of the dataframe twice, once for each strain with one title for each strain.
#a dataframe of only those rows that count for two strains:
doubledcounts1=standardized_enrichment[double_positions,]


#Splits the name into two separate names, one for each strain
doublenames=matrix(unlist(str_split(doubledcounts1$Virus_Protein_Loc_AA, fixed("."))), ncol=2, byrow=TRUE)

#Duplicate these rows and rename with 1 row as each strain.
doubledcounts2=doubledcounts1
doubledcounts1$Virus_Protein_Loc_AA=doublenames[,1]
doubledcounts2$Virus_Protein_Loc_AA=doublenames[,2]
#combine the new, correctly named copies of these with the rest of the counts dataframe
counts_no_doubles=rbind(standardized_enrichment[-double_positions,], doubledcounts1,doubledcounts2)

#Now we'll need to separate by strains. 
#split out the strain name for each oligo

seplocus=str_split_fixed(counts_no_doubles$Virus_Protein_Loc_AA, fixed("_"),4)
rownames(seplocus)=rownames(counts_no_doubles)
colnames(seplocus)=c("Virus", "Protein", "Locus", "AA")
counts_no_doubles=cbind(counts_no_doubles, seplocus)
counts_no_doubles$Locus=as.numeric(as.character(counts_no_doubles$Locus))
head(seplocus)
# Virus    Protein Locus AA 


# 1 "BG505"  "gp41"  "1"   "G"
# 2 "BG505"  "gp41"  "1"   "E"
# 3 "BG505"  "gp41"  "1"   "I"
# 4 "BG505"  "gp41"  "6"   "G"
# 5 "BG505"  "gp41"  "51"  "G"
# 6 "ZA1197" "gp41"  "155" "G"
table(seplocus[,"Virus"], seplocus[,"Protein"])
# gp41   V3
# BF520  3461  700
# BG505  3461  700
# ZA1197 3460  700
#we have six protein/strain combos represented, so we'll split into 6 data frames

vir=unique(seplocus[,"Virus"])
prot=unique(seplocus[,"Protein"])
# counts_no_doubles$Input.Lib1=NULL
# counts_no_doubles$Input.Lib2a=NULL
# counts_no_doubles$Input.Lib2b=NULL
# counts_no_doubles$Input.Lib2=NULL
# counts_no_doubles$mean_beads=NULL
# counts_no_doubles$Beadsonly.Lib1.1=NULL


samples=colnames(counts_no_doubles[-c(1:3,66:69)])

for (i in 1:length(vir)){
  for (j in 1:length(prot)){
    myfile=counts_no_doubles[intersect(which(counts_no_doubles$Virus==vir[i]), which(counts_no_doubles$Protein==prot[j])),]
    filename=paste(vir[i], prot[j], "SE.csv", sep="_")
    print(filename)
    filepath=paste(mypath, filename, sep="/")
    write.csv(myfile, filepath)
    #file for each virus/prot combo with all AA variants
    #Now we'll separate out just the wild-type for each
    #first load the wild type sequence
    wtseqfile=paste(vir[i], prot[j], "protein", "seq.txt", sep=" ")
    print(wtseqfile)
    wtseqpath=paste(mypath, wtseqfile, sep="/")
    wtseq=read_file(wtseqpath)
    wtseq=unlist(strsplit(wtseq, character(0)))
    is.wt=rep(NA, length(myfile$Locus))
    for (k in 1:length(myfile$Locus)){
      is.wt[k]=wtseq[myfile$Locus[k]]==myfile$AA[k]
    }
    wtseqfile=myfile[which(is.wt),]
    wtfilename=paste("wtonly",filename,sep="_")
    print(wtfilename)
    wtfilepath=paste(mypath,wtfilename,sep="/")
    write.csv(wtseqfile, wtfilepath)
    #for (m in 1:length(samples)){
    #  ycolumn=samples[m]
# 
#     print(ggline(wtseqfile, x="Locus", y = ycolumn,
#                     add = NULL, conf.int = FALSE,
#                     cor.coef = FALSE,
#                     xlab = "Locus", ylab = "Fold enrichment", title=samples[m]))
    }
    
    
  }







