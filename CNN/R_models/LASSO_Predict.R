library(glmnet)
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

load('/home/schnablelab/cmiao/MyRepo/schnablelab/CNN/R_models/LASSOmodel.RData')
model = LASSO_list$model_raw_fold_0
y.predict = as.numeric(as.character(predict(model, newx=x.test, type='class', s=0.01)))
print(length(y.predict))

predict = matrix(y.predict, nrow=height, byrow=TRUE)
print(dim(predict))
real_pred_fn = gsub('.npy', '.lasso.prd.csv', test_fn)
write.table(predict, file=real_pred_fn, row.names=FALSE, col.names=FALSE)
print('Done!')
