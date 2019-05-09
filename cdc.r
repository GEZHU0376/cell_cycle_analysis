# cdc = cell division cycle
library(tidyverse)
#2
cdcdata = read.table(file ="combined.txt", header = TRUE, sep = "\t")
cdcdata 
rownames(cdcdata)
head(cdcdata)

#select rows starting with "cdc15_"
names(cdcdata)
row.names(cdcdata) <- cdcdata$X
nameselect <- cdcdata %>% select(starts_with("cdc15_")) 
nameselect
#created a logical vector with all the NA in an
Means<- apply(cdcdata, 1, mean)
Means

#delete row with 'NA' in it
newselect = na.omit(nameselect)
cdcdata[!is.na(Means)]
View(nameselect)

nameselect2=nameselect %>% filter(complete.cases(.))
head(nameselect)
na.omit(nameselect)
#identify the positions of 'NA'
nameselect %>% filter(!is.na(cdc15_10),)
colSums(is.na(nameselect))
#want all the columns

filter_data = nameselect[!is.na(rowSums(nameselect)),]
colSums(is.na(filter_data))
filter_data
#------------------------------------------------------------------

#4: use substring to get the time portion of cdc15_
times <- gsub("cdc15_","",nameselect %>% colnames() )
print(times)

#mean(nameselect$cdc15_190, na.rm = TRUE)
#5: extract first gene in the table then plot its expression as a function of time. 
gene1 <- nameselect[1,]
gene1 %>% gather(key = timepoint, value = value) %>% 
  ggplot(aes(x=timepoint,y=value, group = 1)) + geom_point() + geom_line()

#6 select the first five genes into a 5x24 array. 
# cor command to find correlation among each pairs. 

five_genes <- nameselect2[1:5,1:24]
five_genes

array_genes <- array(c(times,five_genes),dim = c(5,24,1))
array_genes
# x input has to be a matrix.
cor(t(five_genes), y = NULL, use = "pairwise.complete.obs")

#7 pairs command - create matrix of scatter plots
pairs(t(five_genes))

#8 apply command generate a vector that has the mean expression for each gene across all time points. 
mean_expression <- apply(nameselect, 1, mean, na.rm=TRUE) 
view(mean_expression)

#9 calculate mean and standard deviation of all the row means. 
sd(a, na.rm = TRUE)
mean(a, na.rm = TRUE)

#10
a <- rep(NA,22)
a
b <- rep(NA,22)
b

for (i in 1:22) {
  if(!is.na(nameselect[1,i])&!is.na(nameselect[1,i+1])&!is.na(nameselect[1,i+2])){
    a[i] = mean(c(nameselect[1,i+2],nameselect[1,i]))
    b[i] = nameselect[1,i+1]
  }
}
  
t.test(a,b)

#11a
rownumber = dim(nameselect)[1]
rownumber
#11.b
array_a <- array(data = NA, c(rownumber,22))
array_a

array_b <- array(data = NA, c(rownumber,22))
array_b
#11.c
for (j in 1:rownumber) {
  for (i in 1:22) {
    if(!is.na(nameselect[j,i])&!is.na(nameselect[j,i+1])&!is.na(nameselect[j,i+2])){
      array_a[j,i] = mean(c(nameselect[j,i+2],nameselect[j,i]))
      array_b[j,i] = nameselect[j,i+1]
    }
  }
}
#11.d
a2 <- as.vector(a, mode = "any")
a2
b2 <- as.vector(b, mode = "any")
b2
#11.e
t.test(a2,b2)

#12
Means<- apply(nameselect, 1, mean)
summary(Means)
nrow(nameselect)

#13
nameselect_nona = nameselect
for (j in 1:rownumber) {
  for (i in 2:(ncol(nameselect)-1)) {
    if(is.na(nameselect[j,i])){
      #before and after if na
      nameselect_nona[j,i] = mean(c(nameselect[j,i-1],nameselect[j,i+1]))
    }
  }
}


#14 left over should be the NA located at col1 and last column 
Means<- apply(nameselect_nona, 1, mean)
summary(Means)


#15
no_na <- na.omit(nameselect_nona)
no_na
#16
heatmap(as.matrix(no_na))

#17
pc.cdc <- princomp(no_na)
pc.cdc

#plot the cell cycle data
plot(pc.cdc)
summary(pc.cdc)

pc.cdc$loadings 
plot(1:ncol(no_na),pc.cdc$loadings[,1], type = "o")
plot(1:ncol(no_na),pc.cdc$loadings[,2], type = "o")
plot(1:ncol(no_na),pc.cdc$loadings[,3], type = "o")
plot(1:ncol(no_na),pc.cdc$loadings[,4], type = "o")

#20
pc.cdc$scores
pc.cdc.scores.sum <- rowSums(abs(pc.cdc$scores))
pc.cdc.scores.sum

per_1 = abs(pc.cdc$scores[,1]/pc.cdc.scores.sum)
per_1
per_2 = abs(pc.cdc$scores[,4]/pc.cdc.scores.sum)
per_2

sort_value_1 = sort(per_1, decreasing = TRUE)
sort_value_1
gene1 <- names(sort_value_1[1])
gene1

gene2 <- names(sort_value_1[length(sort_value_1)])
gene2

sort_value_2 = sort(per_2, decreasing = TRUE)
sort_value_2

gene3 <- names(sort_value_2[1])
gene3

gene4 <- names(sort_value_2[length(sort_value_2)])
gene4

#plot all the graph and use legend to distinguish them

plot(times,no_na[gene1,], type = "b")

plot(times,no_na[gene2,], type = "b", col = "red", lwd =1)
lines(times,no_na[gene2,], type = "b", col = "red", lwd =1)

points(times,no_na[gene3,], type = "b", col = "green", pch ="+")
lines(times,no_na[gene3,], type = "b", col = "green", pch ="+")

points(times,no_na[gene4,], type = "b", col = "blue" ,pch ="*")
lines(times,no_na[gene4,],  type = "b", col = "blue", pch ="*")

legend("topright", legend=c("gene2", "gene3", "gene4"),
       col=c("red", "green", "blue"), lty=1:2, cex=0.8)

