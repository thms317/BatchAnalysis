import numpy as np

n_min=10
X = [1,2,3,4,5,6,7,8,9,10]

resample_i = np.floor(np.random.rand(n_min) * len(X)).astype(int)

X = np.array(X)
print(resample_i)
X_resample = np.array(X[resample_i])
print(X_resample)
# X_resample = np.array(X[resample_i])
#
#
# x = bootstrap_resample(X,n_min)
# print(x)
