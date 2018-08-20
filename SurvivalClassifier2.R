#creating a function called rmcfsSchi that takes these parameters


options(java.parameters = "-Xmx16g")
#setting java parameter to 16 GB otherwise it will cause out of memory #error for the java heap space
#including library rmcfs 
library(rmcfs)

print("start :")
print(date())
#Reading the csv file in a data frame
# df <- as.data.frame(read.csv("/home/saryou/",header = TRUE))
df <- as.data.frame(read.csv("/home/saryou/SurvivalClassifierMCFS.csv",header = TRUE))

# calling mcfs with the required parameters 
#1:name of the decision coloumn
#2:The dataframe
#3:The number of projections
#4:The projections size
#5:The number of cutoff permutations
#6:setting two parameters the Final cross validation (FinalCV and #FinalRuleset to false) if you need them to be true just change the #script
#if you want to add other parameters check the rmcfs manual first
#7:Threads number are reminiscent to the number of cores you want the #rmcfs to work on
#result <- mcfs(decision~., df, threadsNumber=12)
result=Boruta(decision~., data = as.data.frame(Decisiontable), doTrace = 2)

#save the result to RDS output file

#saveRDS(result,"/home/saryou/Patricia/outexemple.rds")
#export.result(result, path='/home/saryou/Patricia/', label = "outputexemple", zip = TRUE)
export.result(result, path = "/home/saryou/outputclassifier”, label = "output01”, zip = TRUE)
saveRDS(result,"/home/saryou/out02.rds")

#print("end :")
#print(date())

#In the same script you call the function rmcfsSchi 
#For the schezophernia methylation sites filtered data I called the #functions with the below parameters
#Function call parameters :Input File path and name, #projections,projectionsize, permutations,threads,output file path and #na