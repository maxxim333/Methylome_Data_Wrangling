#DataWrangling Exercise

#vever

#This dataset is of a DNA CpG islands methylation profiling arra by Illumina of patients
#with 3 different Crohn´s diseases (and healthy).

#The goal of this exercise is more oriented toward Data Wrangling with R and training using 
#some of the libraries, rather than actually analyzing methylomes. I just picket the methylome
#dataset for fun. I don´t actually know if this is the correct pipeline of methylome analysis
#(it´s probably not). The scope and the time I dedicated to this exercise was too small to
#actually learn to analyze such datasets and this is based on only one paper I skim-read.

#I will pick a subset of data and work with it because the original data is too big

setwd("C:/Users/Maksym/Desktop/UPM/Analisis Estadistico/datawrang")
mydata <- read.table("GSE99788_methylated_and_unmethylated_signals.txt", header = TRUE, sep = "\t", dec = ".", row.names = 1) 

#Pick subset
mydata2 <-mydata[1:10000,]
newdata <- write.csv(mydata2,file="newdata.csv", sep = ",",row.names = TRUE, col.names = TRUE)
finaldata <- read.table("newdata.csv", header = TRUE, sep = ",", dec = ".", row.names = 1) 

#Packages
install.packages("dplyr")
library(dplyr)
library( stringr)
library(tibble)  # for `rownames_to_column` and `column_to_rownames`

#Check the names of the data
names(finaldata)


#I removed the column indicating the quality (the "confidence") of methylation assessment
#of each CpG
newdata_1 <- finaldata[1:10000,c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26,28,29,31,32,34,35,37,38,40,41,43,44,46,47,49,50,52,53)]
names(newdata_1)

#According to the paper I read, "At each CpG site, methylation is quantified by the beta 
#value b:=M/(M+U+a), where M>0 and U>0 denote the methylated and unmethylated signal 
#intensities, respectively, measured by the Illumina 450k array. The offset a???0 is usually 
#set equal to 100 and is added to M+U to stabilize beta values when both M and U are small.
#From https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5120494/

#Apply beta values
newdata_beta <- as.data.frame(c((newdata_1[1]/(newdata_1[1]+newdata_1[2]+100)),newdata_1[3]/(newdata_1[3]+newdata_1[4]+100),
                              newdata_1[5]/(newdata_1[5]+newdata_1[6]+100),newdata_1[7]/(newdata_1[7]+newdata_1[8]+100),
                              newdata_1[9]/(newdata_1[9]+newdata_1[10]+100),newdata_1[11]/(newdata_1[11]+newdata_1[12]+100),
                              newdata_1[13]/(newdata_1[13]+newdata_1[14]+100),newdata_1[15]/(newdata_1[15]+newdata_1[16]+100),
                              newdata_1[17]/(newdata_1[17]+newdata_1[18]+100),newdata_1[19]/(newdata_1[19]+newdata_1[20]+100),
                              newdata_1[21]/(newdata_1[21]+newdata_1[22]+100),newdata_1[23]/(newdata_1[23]+newdata_1[24]+100),
                              newdata_1[25]/(newdata_1[25]+newdata_1[26]+100),newdata_1[27]/(newdata_1[27]+newdata_1[28]+100),
                              newdata_1[29]/(newdata_1[29]+newdata_1[30]+100),newdata_1[31]/(newdata_1[31]+newdata_1[32]+100),
                              newdata_1[33]/(newdata_1[33]+newdata_1[34]+100),newdata_1[35]/(newdata_1[35]+newdata_1[36]+100)),row.names = rownames(newdata_1))
names(newdata_beta)
#This dataset has methylation data of healthy individuals and 
#individuals with different types of Crohn Diseases (inflamatory, non-inflamatory,
#and stenosis). I want to separate the dataset by types of disease
#Joining by types of Crohn Diseases
healthy <-subset(newdata_beta[,c(13,17,1,3,14)])
inf <-subset(newdata_beta[,4:5])
noninf <-subset(newdata_beta[,c(2,6,7,8,9,15,18)])
sten <-subset(newdata_beta[,c(10,11,16,12)])
names(newdata_beta)

#Now, I want to calculate the mean and standart deviations of methylations of patients
#grouped by disease type
rownames(newdata_beta) #CpG are annotated according to Illumina internal codes. In real applications, probably it would be a good idea to retreive the glocation on genome of these CpG islands
meanhealthy <- as.data.frame(rowMeans(healthy), row.names=rownames(newdata_beta))
meanhealthysd <- as.data.frame(apply(healthy, FUN=sd, MARGIN = 1, na.rm = FALSE),row.names=rownames(newdata_beta) )

meaninf <- as.data.frame(rowMeans(inf), row.names=rownames(newdata_beta))
meaninfsd <- as.data.frame(apply(inf, FUN=sd, MARGIN = 1, na.rm = FALSE),row.names=rownames(newdata_beta) )


meansten <- as.data.frame(rowMeans(sten), row.names=rownames(newdata_beta))
meanstensd  <- as.data.frame(apply(sten, FUN=sd, MARGIN = 1, na.rm = FALSE),row.names=rownames(newdata_beta) )


meannoninf <- as.data.frame(rowMeans(noninf), row.names=rownames(newdata_beta))
meannoninfsd <- as.data.frame(apply(noninf, FUN=sd, MARGIN = 1, na.rm = FALSE),row.names=rownames(newdata_beta) )




#Cbind to have the full dataset of processed data
fulldata <- cbind.data.frame(meanhealthy, meanhealthysd, meaninf, meaninfsd,meannoninf, meannoninfsd,meansten, meanstensd)

#How is the distribution of methylation?
hist(fulldata$`rowMeans(healthy)`)
?hist
require(ggplot2)
require(reshape2)

dfz <- data.frame(healthy = fulldata$`rowMeans(healthy)`,
                 inflamed = fulldata$`rowMeans(inf)`,
                 non_inflamed = fulldata$`rowMeans(noninf)`,
                 stenosis = fulldata$`rowMeans(sten)`)



ggplot(melt(dfz), aes(value, fill = variable)) + geom_histogram(position = "dodge")


#Summary doesnt really provide us any useful information in here because methylation data
#is very complex
summary(fulldata)

#Now I want to perform a t.test for each row (CpG island) between healthy and inflamatory
#Crohn disease patients. I understand that in real methylome study it would not be so straight
#forward as that and more control and data cleaning should be performed, probably including
#normalization by background methylation levels etc...
#I picked up an even smaller subset of CpG island because of limited computation power of my PC
help("t.test")

testsubset <- subset(healthy[1:100,])
testsubset2 <- subset(inf[1:100,])


healthyinf<-apply(testsubset,testsubset2, FUN=t.test, MARGIN=1)
summary(healthyinf)
class(healthyinf)

#Retreive only pvalues. It is returned in a very weird format. I did a loop to paste
#the pvalues in a new vector and then transformed this vector in a dataframe
a<-do.call(rbind,lapply(healthyinf,paste0))
data.frame(pos=rownames(a),word=a)
class(a)
pvalues <- c()
x <- 0
for (i in 1:100) {
  pvalues <- c(pvalues,(a[i,][3]))
  x <- x+1
}

pvaluesdf<-as.data.frame(pvalues[1:100])
row.names(pvaluesdf) <- rownames(a)
colnames(pvaluesdf) <-c("pvalues")


#Now this will show the cpg sites with significantly different methylation between healthy and inflamatory Chron Disease individuals
#If we apply Bonferroni correction for multiple hipothesis testing, conf level should be
#1-(0.05/100)=1-0.0005=0.9995. So, only <0.0005 will be significative.
#First, a conversion to numeric is needed.
pvaluesdf$pvalues <-as.numeric(as.character(pvaluesdf$pvalues))
pvaluesdffiltered <- pvaluesdf %>%
  rownames_to_column('gene') %>%
  filter_if(is.numeric, all_vars(. < 0.0005)) %>%
  column_to_rownames('gene')
pvaluesdffiltered

#Now we can plot some things. For this is better to work with full data.
summary(fulldata)
names(fulldata)

#Here I plot with ggplot2, mean methylation levels of each CpG for healthy
#(x axis) and inflamatory Crohn disease (y axis) and draw a correlation line
#(linear regression). Same could be done for other types of Crohn diseases and compared
#Or even compared all types of Crohn together VS healthy
library( ggplot2)
help(ggplot2)
ggplot(fulldata, aes(x=rowMeans(healthy), y=rowMeans(inf))) + geom_point(size=0.5)+ geom_smooth(method = "lm", se = TRUE)


#In the graph above, we see that in general, there will be a very high correlation
#of methylation levels between disease and healthy, which is not surprising, as
#we expect only few CpG islands to be differentially methylated in two conditions
#Perhaps, if we want to see the global levels of methylation, a more interesting
#approach would be to see boxplots

#We see in these boxplots that there seems to be no difference in global distribution of
#methylation levels between conditions.
minidata<-fulldata[,c(1,3,5,7)]
ggplot(stack(minidata), aes(x = ind, y = values),outline=FALSE) + geom_boxplot()

