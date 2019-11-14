library(e1071)
library(RcppCNPy)

args=(commandArgs(trailingOnly = F))
len_args = length(args)
if(len_args==0){
    print("No arguments supplied. specify input_csv_file height")
  }
test_fn = args[len_args-1]; test_fn = sub("-","",test_fn)
height = args[len_args]; height = as.numeric(sub("-","",height))
print(test_fn);print(height)

x.test = npyLoad(test_fn)/255
print(dim(x.test))

load('/home/schnablelab/cmiao/MyRepo/schnablelab/CNN/R_models/SVMmodel.RData')
model = SVM_list$model_raw_fold_3

pred_prob <- predict(model, x.test, decision.values = TRUE, probability = TRUE)
pred.prob.matrix=attr(pred_prob, "probabilities")
labels = as.numeric(colnames(pred.prob.matrix))
y.predict = labels[apply(pred.prob.matrix,1,which.max)]
print(length(y.predict))

predict = matrix(y.predict, nrow=height, byrow=TRUE)
print(dim(predict))
real_pred_fn = gsub('.npy', '.svm.prd.csv', test_fn)
write.table(predict, file=real_pred_fn, row.names=FALSE, col.names=FALSE)
print('Done!')
