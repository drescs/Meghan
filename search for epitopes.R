library(readr)
library(stringr)
library(reshape2)
library(ggpubr)
mypath="C:/Users/Sara/Documents/HHMI_projects/Meghan_Garrett"
pathdate="2018.07.25"
fullpath=paste(mypath, pathdate, sep="/")
proportions=read.csv("C:/Users/Sara/Documents/HHMI_projects/Meghan_Garrett/2018.07.25.proportions.csv", header=T, stringsAsFactors = F)

#Next, we'll make datasets for each strain
#However, some oligomers are commnon between two strains


#pretty sure I don't use these, will delete later
#tenng=which(sepnames[,2]%in%c("", "only", "1", "2", "10ng"))
#shortcolname=paste(sepnames[,1], sepnames[,2],sep=".")
#counts=rbind(counts, shortcolname)
#sepnames[which(sepnames[,2]==""),2]="ng"
#sepnames[which(sepnames[,3]==""),3]="rep"
#ng=rep(NA, length(colnames(counts))

#find the oligomers common to two strains
#if there are two strains the names are separated by "."

find_doubles=str_split(proportions$Virus_Protein_Loc_AA, fixed("."))
double_positions=vector()
for (i in 1:length(find_doubles)){
  if (length(find_doubles[[i]])>1){
    double_positions=c(double_positions,i)
  }
}


#double_positions gives the rows where one oligomer is the same between two strains
#we want to make this portion of the dataframe twice, once for each strain with one title for each strain.
#a dataframe of only those rows that count for two strains:
doubledcounts1=proportions[double_positions,]


#Splits the name into two separate names, one for each strain
doublenames=matrix(unlist(str_split(doubledcounts1$Virus_Protein_Loc_AA, fixed("."))), ncol=2, byrow=TRUE)

#Duplicate these rows and rename with 1 row as each strain.
doubledcounts2=doubledcounts1
doubledcounts1$Virus_Protein_Loc_AA=doublenames[,1]
doubledcounts2$Virus_Protein_Loc_AA=doublenames[,2]
#combine the new, correctly named copies of these with the rest of the counts dataframe
propnd=rbind(proportions[-double_positions,], doubledcounts1,doubledcounts2)

#Now we'll need to separate by strains. 
#split out the strain name for each oligo

seplocus=str_split_fixed(propnd$Virus_Protein_Loc_AA, fixed("_"),4)
rownames(seplocus)=rownames(propnd)
colnames(seplocus)=c("Virus", "Protein", "Locus", "AA")
propnd=cbind(propnd, seplocus)
propnd$Locus=as.numeric(as.character(propnd$Locus))
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
propnd$input.1=NULL
propnd$input.2=NULL
propnd$X=NULL
#this is just the actual samples, not input or metadata
samples=colnames(propnd[-c(1:3,46,67:71)])
for (i in 1:length(vir)){
  for (j in 1:length(prot)){
    myfile=propnd[intersect(which(propnd$Virus==vir[i]), which(propnd$Protein==prot[j])),]
    filename=paste(vir[i], prot[j], "prop.csv", sep="_")
    filepath=paste(fullpath, filename, sep="/")
    write.csv(myfile, filepath)
    #file for each virus/prot combo with all AA variants
    #Now we'll separate out just the wild-type for each
    #first load the wild type sequence
    wtseqfile=paste(vir[i], prot[j], "protein", "seq.txt", sep=" ")
    wtseqpath=paste(mypath, wtseqfile, sep="/")
    wtseq=read_file(wtseqpath)
    wtseq=unlist(strsplit(wtseq, character(0)))
    is.wt=rep(NA, length(myfile$Locus))
    for (k in 1:length(myfile$Locus)){
      is.wt[k]=wtseq[myfile$Locus[k]]==myfile$AA[k]
    }
    wtseqfile=myfile[which(is.wt),]
    wtfilename=paste("wtonly",filename,sep="_")
    wtfilepath=paste(fullpath,wtfilename,sep="/")
    write.csv(wtseqfile, wtfilepath)
    
  }
}

#same forloop structure, now to make the summed totals for each position

for (i in 1:length(vir)){
  for (j in 1:length(prot)){
    infilename=paste("wtonly", vir[i], prot[j], "prop.csv", sep="_")
    infile=read.csv(paste(fullpath, infilename, sep="/"))
    #We need our data in order by locus
    sorted=infile[order(infile$Locus),]
    newsorted=sorted
    pdffile=paste(vir[i], prot[j], "wt_foldenrich_totals.pdf", sep="_")
    pdfpath=paste(fullpath, pdffile, sep="/")
    pdf(file = pdfpath)
    
    for (m in 1:length(samples)){
      locustotal=rep(NA, length(sorted$Locus))
      for (k in 1:length(sorted$Locus)){
        overlap=intersect(c((k-15):(k+15)),c(1:length(sorted$Locus)))
        locustotal[k]=(sum(sorted[overlap,samples[m]])-sum(sorted[overlap, "beads.only.2"]))/sum(sorted[overlap,"meaninput"])
      }
      newsorted[,samples[m]]=locustotal
      ycolumn=samples[m]
      
      print(ggline(newsorted, x="Locus", y = ycolumn,
                   add = NULL, conf.int = FALSE,
                   cor.coef = FALSE,
                   xlab = "Locus", ylab = "Fold enrichment", title=samples[m]))
      
    } 
    dev.off()
    outfilename=paste("wtonly", vir[i], prot[j], "locustotal.csv", sep="_")
    outfilepath=paste(fullpath, outfilename, sep="/")
    write.csv(newsorted, outfilepath)
  }
}











