import numpy as np
from numpy import empty
np.set_printoptions(threshold=np.inf)

dmatrix = np.zeros((69,69))
#a = np.load('mean1.npy')
#b = np.load('mean2.npy')

for i in range(1,70):
    for j in range(1,70):
        a = np.load('mean'+str(i)+'.npy')
        b = np.load('mean'+str(j)+'.npy')
        dmatrix[i-1,j-1]= np.linalg.norm(a-b)

for k in range(1,70):
    print('['+str(k)+']'+str(dmatrix[k-1]))