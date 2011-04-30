# analyze a bayesallfile
# using coda
#
#
# User question: directory
cat("\n","Specify the working directory [return for current directory]","\n") 
workingdir<-scan(n=1)
ifelse(is.null(workingdir), getwd(), setwd(workingdir))
#
# User question: bayesallfile
cat("\n","Specify the filename that was used to record the Bayes mcmc output [return for bayesallfile]","\n") 
bfile<-scan(n=1)
ifelse(is.null(bfile), bayesallfile<-"bayesallfile", bayesallfile<-bfile) 

# Load CODA
library(coda)


# Load the bayesallfile
data<-read.table(bayesallfile, header=T)

# Split the whole list into the individual replicates, so that you will
# have a list of lists
data.list<-split(data,data$Replicate)

# Subset the parameters of interest, in other words, you don't need the
# step count for instance. In this case, chose columns 9 through 34
data.list.1<-lapply(data.list,subset,select=c(9:34)

# Transform each element of the list (or each replicate) into an mcmc
# object, which is required by coda
data.list.2<-lapply(data.list.1,mcmc) # Here, you can set: start, end,
# and thin in case you don't want all the data. If left as is, it will
# use all the data available

# Transform the list of mcmc objects above into an mcmc.list, the
# object required by gelman.diag(), for instance:
data.list.3<-mcmc.list(data.list.2)

# Now, the data should be ready for the CODA package.
# For instance, to plot all the traces and densities for each
# parameter - Replicates will be plotted in different colours for each
# parameter:
plot(data.list.3, ask=T) #The ask=T will force R to plot a page with
# parameters, then stop until you hit <return>, then it will start to
# print the next page.

# Find out what other commands are available on CODA:
library(help=coda)

#End of commands