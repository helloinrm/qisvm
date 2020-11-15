# using LIBSVM
# function svm3(X::Array{Float64,2}, y::Array{Float64}, x::Array{Float64})
#     model = svmtrain(X, y, kernel = Kernel.Linear)
#     (predicted_labels, decision_values) = svmpredict(model, x)
#     return predicted_labels
# end
# n=m=10000
# ratio=0.2
# train_set=rand(n,m-Int(m*ratio))
# test_set=rand(n,Int(m*ratio))
# train_label=rand([0.0,1.0],m-Int(m*ratio))
# test_label=rand([0.0,1.0],Int(m*ratio))
# svm3(train_set, train_label,[test_set train_set])

# train_set,test_set,train_label,test_label,a=generate(n,m,1,0.2,0)
# set=[test_set train_set]
# label3 = svm4(train_set, train_label,set)
# s3 = acc(label3[1:size(test_label)[1]], test_label)
# c3 = acc(label3[(size(test_label)[1]+1):end], train_label)
# (s3,c3)|>println
