
setwd('/home/divyae/divyae2/RFtrain/FinalSets/6')
library(randomForest)
N<-read.table(file='TrainNegatives.tsv',sep= '\t')
P<-read.table(file='TrainPositives.tsv',sep='\t')

total <- rbind(N,P)
total$V1<-as.factor(total$V1)



which(is.na(P), arr.ind=TRUE)
#row col
#[1,] 5466  29
which(is.na(N), arr.ind=TRUE)
#row col
#[1,] 3098  29

table(total$V1)
t1<-total
VB_pairs<-total$V2
total$V2<-NULL
total$V3<-NULL
total$V4<-NULL

total1<-total

set.seed(123)

colnames(total) <- c("Class","Blen","S_len","Bscore","S_score","eval","Bbscore","BbscoreBYBlen","S_bscore","S_BbscoreBYBlen",
                     "BC1", "BC2", "BC3", "BC4", "BC5", "BC6", 
                     "CB1", "CB2", "CB3", "CB4", "CB5", "CB6",
                     "Cos1", "Cos2", "Cos3", "Cos4", "Cos5", "Cos6",
                     "Crr1", "Crr2", "Crr3", "Crr4", "Crr5", "Crr6",
                     "ChebyS1", "ChebyS2", "ChebyS3", "ChebyS4", "ChebyS5", "ChebyS6",
                     "Euc1", "Euc2", "Euc3", "Euc4", "Euc5", "Euc6"
)

samp  <- sample(nrow(total), 0.75 * nrow(total))
train <- total[samp, ]
test  <- total[-samp, ]


model<-randomForest(Class ~ .,data=train,ntree=200,na.action = na.omit)
model
#oob.err=double(39)
#test.err=double(39)

plot(model)
importance(model)
varImpPlot(model)
Prediction <- predict(model, test)
table(Prediction, test$Class)
sum(diag(table(Prediction, test$Class)))/nrow(test)
#plot(margin(model, test$Class))


var.imp <- data.frame(importance(model,type=2))
var.imp$Variable <-row.names(var.imp)
var.imp[order(var.imp$MeanDecreaseGini,decreasing = T),]


#########################################################
#everything same as above, just 1/10th number of tress

t<-total
##t$V2<-NULL

#colnames(t) <- c("Class","qlen","DBlen","Blen","S_len","Bscore","S_score","eval","Bbscore","BbscoreBYBlen","S_bscore","S_BbscoreBYBlen",
#                     "Euc1", "Euc2", "Euc3", "Euc4", "Euc5",
#                     "BC1", "BC2", "BC3", "BC4", "BC5", 
#                     "CB1", "CB2", "CB3", "CB4", "CB5", 
#                     "Cos1", "Cos2", "Cos3", "Cos4", "Cos5",
#                     "Crr1", "Crr2", "Crr3", "Crr4", "Crr5",
#                     "ChebyS1", "ChebyS2", "ChebyS3", "ChebyS4", "ChebyS5")




samp4 <- sample(nrow(t), 0.75 * nrow(t))
train4 <- t[samp4, ]
test4 <- t[-samp4, ]


model4<-randomForest(Class ~ .,data=train4,ntree=200,na.action = na.omit,importance=T)
model4

plot(model4)
importance(model4)
varImpPlot(model4)
Prediction4 <- predict(model4, test4)
table(Prediction4, test4$Class)
sum(diag(table(Prediction4, test4$Class)))/nrow(test4)
#plot(margin(model4, test4$Class))
#MDSplot(model4, fac, k=2)

var.imp4 <- data.frame(importance(model4,type=2))
var.imp4$Variable <-row.names(var.imp4)
var.imp4[order(var.imp4$MeanDecreaseGini,decreasing = T),]


#par(mfrow=c(1,1))
#par(pty="s")
#varImpPlot(model4, type=1, pch=19, col=1, cex=.5, main="")
#varImpPlot(model4, type=2, pch=19, col=1, cex=.5, main="")
#dev.off()

#########################################################
#a RF after removing blast hits and scores, and JUST USING KMERS


total1$V5<-NULL
total1$V6<-NULL
total1$V7<-NULL
total1$V8<-NULL
total1$V9<-NULL
total1$V10<-NULL
total1$V11<-NULL
total1$V12<-NULL
total1$V13<-NULL

colnames(total1) <- c("Class", "BC1", "BC2", "BC3", "BC4", "BC5", "BC6", 
                      "CB1", "CB2", "CB3", "CB4", "CB5", "CB6",
                      "Cos1", "Cos2", "Cos3", "Cos4", "Cos5", "Cos6",
                      "Crr1", "Crr2", "Crr3", "Crr4", "Crr5", "Crr6",
                      "ChebyS1", "ChebyS2", "ChebyS3", "ChebyS4", "ChebyS5", "ChebyS6",
                      "Euc1", "Euc2", "Euc3", "Euc4", "Euc5", "Euc6")

# "Euc1", "Euc2", "Euc3", "Euc4", "Euc5", "Euc6",
# "BC1", "BC2", "BC3", "BC4", "BC5", "BC6",
# "CB1", "CB2", "CB3", "CB4", "CB5", "CB6",
# "Cos1", "Cos2", "Cos3", "Cos4", "Cos5", "Cos6",
# "Crr1", "Crr2", "Crr3", "Crr4", "Crr5", "Crr6",
# "ChebyS1", "ChebyS2", "ChebyS3", "ChebyS4", "ChebyS5", "ChebyS6",

samp <- sample(nrow(total1), 0.75 * nrow(total1))
train1 <- total1[samp, ]

test1 <- total1[-samp, ]

model2<-randomForest(Class ~ .,data=train1,ntree=1000,na.action = na.omit,importance=T)
model2
plot(model2)
plot(model2,log="x")
importance(model2)
varImpPlot(model2, sort =T, main= "Variable Importance", n.var=10)
varImpPlot(model2, sort =T, main= "Variable Importance", n.var=20)
Prediction2 <- predict(model2, test1)
table(Prediction2, test1$Class)
sum(diag(table(Prediction2, test1$Class)))/nrow(test1)
#plot(margin(model2, test1$Class))

#plot(margin(model2,sort=T))

#getTree(model2,k=2,labelVar = FALSE)

var.imp2 <- data.frame(importance(model2,type=2))
var.imp2$Variable <-row.names(var.imp2)
var.imp2[order(var.imp2$MeanDecreaseGini,decreasing = T),]

####################


samp3 <- sample(nrow(t), 0.75 * nrow(t))
train3 <- t[samp3, ]
test3 <- t[-samp3, ]


model3<-randomForest(Class ~ .,data=train4,ntree=5000,na.action = na.omit,importance=T)
model3

plot(model3)
importance(model3)
varImpPlot(model3)
Prediction3 <- predict(model3, test3)
table(Prediction4, test3$Class)
sum(diag(table(Prediction3, test3$Class)))/nrow(test3)
#plot(margin(model4, test4$Class))

var.imp3 <- data.frame(importance(model3,type=2))
var.imp3$Variable <-row.names(var.imp3)
var.imp3[order(var.imp3$MeanDecreaseGini,decreasing = T),]


########



#model2 = randomForest(Class ~ .,data=train1,ntree=1000) 
#print(model2)
#plot(model2,log="x",main="black default, red samplesize, green tree depth")
#plot(model2,log="x")

#reduced sample size
#model6 = randomForest(Class ~ .,data=train1,ntree=1000,sampsize=.4*length(train1)) 
#print(model6)
#plot(model6,log="x")

#limiting tree depth (not exact )
#rf3 = randomForest(Class ~ .,data=train1,ntree=1000,sampsize=.5*length(train1),maxnodes=6)
#print(rf3)
#plot(rf3,log="x")
