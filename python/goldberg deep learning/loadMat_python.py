# -*- coding: utf-8 -*-
"""
Created on Thu Feb 18 00:26:04 2016

@author: Peter Rupprecht (c) 2016
"""

## load mat data into python

import scipy.io as sio
import os, glob
import numpy as np
import matplotlib.pyplot as plt

# this folder contains the subfolders with mat files, which in turn contain the downsampled sound files
os.chdir('/home/pierre/_goldberg')

folderlist = os.walk('/home/pierre/_goldberg').next()[1]

# size of dataset to be loaded (limited by my PC RAM >> no hard constraint)
# 20000 10-sec fragments for both training and test dataset
Total = np.zeros((63000, 20000*2))
Total_label = np.zeros((20000*2,1))
counter = 0
for i in range(0,7):
    print i
    os.chdir(folderlist[i])
    FileList = glob.glob("*.mat")
    folder_counter = 0
    for k in range(0,len(FileList)):
        if folder_counter < 19500:
            temp_data = sio.loadmat(FileList[k])
            temp_data2 = np.asarray(temp_data['chunked_song'])
            if i == 2: # in my case, only the second folder is Schubert ...
                Total_label[counter:counter+temp_data2.shape[1]] = 1 # Schubert
                Total[:,counter:counter+temp_data2.shape[1]] = temp_data2
            else: # ... everything else is Bach
                Total_label[counter:counter+temp_data2.shape[1]] = 2 # Bach
                Total[:,counter:counter+temp_data2.shape[1]] = temp_data2
            counter += temp_data2.shape[1]
            folder_counter += temp_data2.shape[1]
    os.chdir("../")
Total = Total[:,0:counter]
Total_label = Total_label[0:counter]

## generate test and training datasets

full_dataset = np.swapaxes(Total,0,1)
full_labels = Total_label

sum(full_labels == 2)
sum(full_labels == 1)

## random permutations

np.random.seed(1383)
def randomize(dataset, labels):
  permutation = np.random.permutation(labels.shape[0])
  shuffled_dataset = dataset[permutation,:]
  shuffled_labels = labels[permutation]
  return shuffled_dataset, shuffled_labels
full_dataset, full_labels = randomize(full_dataset, full_labels)

train_dataset = full_dataset[1:30000,:]
train_labels = full_labels[1:30000,:]
test_dataset = full_dataset[30001:,:]
test_labels = full_labels[30001:,:]

# delete unused variables
del full_dataset, full_labels

print('Training', train_dataset.shape, train_labels.shape)
print('Testing', test_dataset.shape, test_labels.shape)

# save data ...

np.save('goldberg_train_dataset.npy',train_dataset)
np.save('goldberg_train_labels.npy',train_labels)
np.save('goldberg_test_dataset.npy',test_dataset)
np.save('goldberg_test_labels.npy',test_labels)

# load data

train_dataset = np.load('goldberg_train_dataset.npy')
train_labels = np.load('goldberg_train_labels.npy')
test_dataset = np.load('goldberg_test_dataset.npy')
test_labels = np.load('goldberg_test_labels.npy')

train_labels = np.ravel(train_labels) - 1
test_labels = np.ravel(test_labels) - 1