#Here I want to simpliy her dataset in order to to some practice computations/analyses.
#Step 1 is to import an annotated counts file:
library(stringr)
counts=read.table("C:/Users/Sara/Documents/HHMI_projects/Meghan_Garrett/2018.07.25.annotatedCounts.txt", header=T, stringsAsFactors = F)

#pull out relevant parts of names here
sepnames=str_split_fixed(colnames(counts), fixed("."), 3)
colnames(sepnames)=c("ID", "ng", "rep")
#going to trim this down to include only rep1, and samples with 10ng of input for the moment. 
#Later will do some quality control and decide where to average or drop relicates
rep1=which(sepnames[,3]%in%c("", "1"))
counts_rep1=counts[,rep1]
table(sepnames)
tenng=which(sepnames[,2]%in%c("", "only", "1", "2", "10ng"))

sepnames[which(sepnames[,2]==""),2]="ng"
sepnames[which(sepnames[,3]==""),3]="rep"

keep_samples=intersect(rep1, tenng)
rownames(sepnames)=colnames(counts)

#Next step will be to look first at the native AA sequence. I'll pick one for now.
#Keep oligos from BF520 with original amino acid

tempbf520="AVGIGAVFIGFLGAAGSTMGAASVTLTVQARQLLSGIVQQQSNLLRAIEAQQHLLKLTVWGIKQLQARVLAVERYLKDQQLLGIWGCSGKLICTTNVPWNSSWSNKSQDEIWGNMTWLQWDKEVSNYTQIIYTLIEESQNQQEKNEQDLLALDKWASLWNWFNISQWLWYIKI"
bf520=unlist(strsplit(tempbf520, character(0)))

seplocus=str_split_fixed(counts$Virus_Protein_Loc_AA, fixed("_"),4)

rownames(seplocus)=rownames(counts)
counts=cbind(seplocus, counts)

originalaa=vector()
for (i in 1:length(seplocus[,1])){
  originalaa[i]=(bf520[as.numeric(seplocus[i,3])]==seplocus[i,4])
}

keep_AA=which(originalaa )
keep_BF520=which(seplocus[,1]=="BF520")
keep_locus=intersect(keep_AA, keep_BF520)
toy_counts=counts[keep_locus, keep_samples]



new_toy=sweep(toy_counts[-c(1:5)], 1, (toy_counts$beads.only.1), FUN = "-", check.margin = TRUE)

toy_means=apply(new_toy, 2, mean)

final_toys=sweep(new_toy, 2, toy_means, "/")
#Now we can do some analysis of this
#we have input and beads only to control/normalize our outputs. Probably good to subtract beads only average from the data and then divide by input?


