from scipy.io import FortranFile
import numpy as np
import h5py

xfile   = FortranFile('./images/vel_y0.05391', 'r')
train_dim = 5000
dev_dim = 1000
test_dim  = 1000

x_train = np.empty([train_dim,128,128,3], 'float64')
for i in range(train_dim):
    x_train[i] = np.moveaxis(xfile.read_reals('float64').reshape(3,128,128), 0, -1)
x_train -= np.mean(x_train, axis=(1,2)).reshape(-1,1,1,3)

x_dev = np.empty([dev_dim,128,128,3], 'float64')
for i in range(dev_dim):
    x_dev[i] = np.moveaxis(xfile.read_reals('float64').reshape(3,128,128), 0, -1)
x_dev -= np.mean(x_dev, axis=(1,2)).reshape(-1,1,1,3)

x_test = np.empty([test_dim,128,128,3], 'float64')
for i in range(test_dim):
    x_test[i] = np.moveaxis(xfile.read_reals('float64').reshape(3,128,128), 0, -1)
x_test -= np.mean(x_test, axis=(1,2)).reshape(-1,1,1,3)
xfile.close()

print("size of train dataset:", x_train.shape)

rx0,rx1 = np.amin(x_train[:,:,:,0]),np.amax(x_train[:,:,:,0])
ry0,ry1 = np.amin(x_train[:,:,:,1]),np.amax(x_train[:,:,:,1])
rz0,rz1 = np.amin(x_train[:,:,:,2]),np.amax(x_train[:,:,:,2])
rmaxx = max( abs(rx0), abs(rx1) )
rmaxy = max( abs(ry0), abs(ry1) )
rmaxz = max( abs(rz0), abs(rz1) )

x_train[:,:,:,0] = x_train[:,:,:,0] / rmaxx
x_train[:,:,:,1] = x_train[:,:,:,1] / rmaxy
x_train[:,:,:,2] = x_train[:,:,:,2] / rmaxz

x_dev[:,:,:,0] = x_dev[:,:,:,0] / rmaxx
x_dev[:,:,:,1] = x_dev[:,:,:,1] / rmaxy
x_dev[:,:,:,2] = x_dev[:,:,:,2] / rmaxz

x_test[:,:,:,0] = x_test[:,:,:,0] / rmaxx
x_test[:,:,:,1] = x_test[:,:,:,1] / rmaxy
x_test[:,:,:,2] = x_test[:,:,:,2] / rmaxz

print( "-----------------HACK ----------FORCING range vx,vy,vz:max")
print("rmaxx,rmaxy,rmaxz=",rmaxx,rmaxy,rmaxz)
minx,maxx = np.amin(x_train[:,:,:,0]),np.amax(x_train[:,:,:,0])
miny,maxy = np.amin(x_train[:,:,:,1]),np.amax(x_train[:,:,:,1])
minz,maxz = np.amin(x_train[:,:,:,2]),np.amax(x_train[:,:,:,2])
print( train_dim, "elements normalized into: x",minx,maxx,"y",miny,maxy,"z",minz,maxz) 

hf = h5py.File('./images/vel_y0.05391.h5', 'w')

hf.create_dataset('train', data=np.swapaxes(x_train,1,2))
hf.create_dataset('dev', data=np.swapaxes(x_dev,1,2))
hf.create_dataset('test', data=np.swapaxes(x_test,1,2))

hf.close()