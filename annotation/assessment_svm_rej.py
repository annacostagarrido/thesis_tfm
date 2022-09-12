

import os 
import numpy as np
import pandas as pd
import time as tm
import scanpy
from sklearn.svm import SVC
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn import preprocessing
from sklearn.calibration import CalibratedClassifierCV
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import StratifiedKFold
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.preprocessing import label_binarize
from sklearn.metrics import roc_curve, auc, f1_score, balanced_accuracy_score
from itertools import cycle
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV 



# Assessment of rejection option

os.chdir("C:/Users/Usuario/Desktop/master/tfm/HCL_data")

    # Load the test data
test_data = pd.read_csv("mca_processed_all.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("refImmScaled.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_refImm.csv", header=0,index_col=None, sep='|')	

# Tranpose genes (columns), cells (rows)
train_data = train_data.T
test_data = test_data.T


# Select the same features (genes) in both datasets
train_data = train_data[train_data.columns.intersection(test_data.columns)]
test_data = test_data[test_data.columns.intersection(train_data.columns)]	
    
print(train_data.shape)
print(test_data.shape)

    
train_data = train_data.sort_index(axis=1)
test_data = test_data.sort_index(axis=1)


# 0) Original prediction for SVM with linear kernel and rejection option

Classifier = LinearSVC()
   
cv = StratifiedKFold(n_splits= 5)
clf = CalibratedClassifierCV(Classifier, cv = cv)
        
clf.fit(train_data, np.ravel(train_labels))

     
# prediction
predicted = clf.predict(test_data)
prob = np.max(clf.predict_proba(test_data), axis = 1)
unlabeled = np.where(prob < 0.7)
predicted[unlabeled] = 'Unknown'

pred = pd.DataFrame(predicted)

np.unique(pred, return_counts=True) # All of them are T cells except 1 (Unknown)


# 1) Reduce the number of samples as if training one unique fold (664 samples)

cv = StratifiedKFold(n_splits=5)

index_x = None
index_y = None
for train_index, test_index in cv.split(train_data, train_labels):
    index_x = train_index
    index_y = test_index
    
X_train = train_data.iloc[index_x]
y_train = train_labels.iloc[index_x]
    

Classifier = LinearSVC()

start=tm.time()
Classifier.fit(X_train, np.ravel(y_train))
tr_time = tm.time()-start

start = tm.time()
predicted = Classifier.predict(test_data)
ts_time = tm.time() - start


pred = pd.DataFrame(predicted)
pred.to_csv("results/" + "SVM_ref_reduced_trainSample.csv", index = False)



# 2) CalibratedClassifierCV with only selected cells seen in both test and reference datasets

    # Load the test
test_data = pd.read_csv("mca_processed_select.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("refImmScaledSelect.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_refImm_select.csv", header=0,index_col=None, sep='|')	

# Tranpose genes (columns), cells (rows)
train_data = train_data.T
test_data = test_data.T


# Select the same features (genes) in both datasets
train_data = train_data[train_data.columns.intersection(test_data.columns)]
test_data = test_data[test_data.columns.intersection(train_data.columns)]	
    
print(train_data.shape)
print(test_data.shape)

    
train_data = train_data.sort_index(axis=1)
test_data = test_data.sort_index(axis=1)


Classifier = LinearSVC()
   
cv = StratifiedKFold(n_splits= 5)
clf = CalibratedClassifierCV(Classifier, cv = cv)
        
clf.fit(train_data, np.ravel(train_labels))

     
# prediction
predicted = clf.predict(test_data)
prob = np.max(clf.predict_proba(test_data), axis = 1)
unlabeled = np.where(prob < 0.7)
predicted[unlabeled] = 'Unknown'


pred = pd.DataFrame(predicted)
pred.to_csv("results/" + "SVM_ref_calibrated_select.csv", index = False)



# 3) CalibratedClassifierCV with prefit


test_data = pd.read_csv("mca_processed_all.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("refImmScaled.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_refImm.csv", header=0,index_col=None, sep='|')	

# Tranpose genes (columns), cells (rows)
train_data = train_data.T
test_data = test_data.T


# Select the same features (genes) in both datasets
train_data = train_data[train_data.columns.intersection(test_data.columns)]
test_data = test_data[test_data.columns.intersection(train_data.columns)]	
    
print(train_data.shape)
print(test_data.shape)

    
train_data = train_data.sort_index(axis=1)
test_data = test_data.sort_index(axis=1)


X_train, X_val, y_train, y_val = train_test_split(train_data, train_labels)

y_val["x"].values.tolist().index("B cells, pro")

X_val.drop(index=X_val.index[31], 
        axis=0, 
        inplace=True)

y_val.drop(index=y_val.index[31], 
        axis=0, 
        inplace=True)


y_train["x"].values.tolist().index("Basophils")


X_train.drop(index=X_train.index[508], 
        axis=0, 
        inplace=True)

y_train.drop(index=y_train.index[508], 
        axis=0, 
        inplace=True)


y_train["x"].values.tolist().index("Microglia")
X_train.drop(index=X_train.index[348], 
        axis=0, 
        inplace=True)

y_train.drop(index=y_train.index[348], 
        axis=0, 
        inplace=True)




Classifier = LinearSVC()

Classifier.fit(X_train, np.ravel(y_train))


clf = CalibratedClassifierCV(Classifier, cv = "prefit")

clf.fit(X_val, np.ravel(y_val))

predicted = clf.predict(test_data)
prob = np.max(clf.predict_proba(test_data), axis = 1)
unlabeled = np.where(prob < 0.7)
predicted[unlabeled] = 'Unknown'


pred = pd.DataFrame(predicted)
pred.to_csv("results/" + "SVM_rej07.csv", index = False)


predicted = clf.predict(test_data)
prob = np.max(clf.predict_proba(test_data), axis = 1)
unlabeled = np.where(prob < 0.3)
predicted[unlabeled] = 'Unknown'


pred = pd.DataFrame(predicted)
pred.to_csv("results/" + "SVM_rej03.csv", index = False)


predicted = clf.predict(test_data)
prob = np.max(clf.predict_proba(test_data), axis = 1)
unlabeled = np.where(prob < 0.1)
predicted[unlabeled] = 'Unknown'


pred = pd.DataFrame(predicted)
pred.to_csv("results/" + "SVM_rej01.csv", index = False)


# 4) SVC with probability = True 

clf_svm = SVC(random_state = 42, C = 1, kernel  = "linear",
              probability= True)

clf_svm.fit(train_data, np.ravel(train_labels))

predicted = clf_svm.predict(test_data)
prob = np.max(clf.predict_proba(test_data), axis = 1)
unlabeled = np.where(prob < 0.1)
predicted[unlabeled] = 'Unknown'




# 5) Train classifier with Logistic Regression

Classifier = LogisticRegression(multi_class="multinomial", max_iter = 10000)

Classifier.fit(train_data, np.ravel(train_labels))

predicted = Classifier.predict(test_data)
prob = np.max(Classifier.predict_proba(test_data), axis = 1)
unlabeled = np.where(prob < 0.7)
predicted[unlabeled] = 'Unknown'



# 6) Train a Random Forest

classifier = RandomForestClassifier(n_estimators = 100, 
                                    random_state = 42, bootstrap = True)
classifier.fit(train_data, np.ravel(train_labels))


predicted = classifier.predict(test_data)
prob = np.max(classifier.predict_proba(test_data), axis = 1)
unlabeled = np.where(prob < 0.7)
predicted[unlabeled] = 'Unknown'

pred = pd.DataFrame(predicted)
pred.to_csv("results/" + "RF_07.csv", index = False)


predicted = classifier.predict(test_data)
prob = np.max(classifier.predict_proba(test_data), axis = 1)
unlabeled = np.where(prob < 0.3)
predicted[unlabeled] = 'Unknown'

pred = pd.DataFrame(predicted)
pred.to_csv("results/" + "RF_03.csv", index = False)


predicted = classifier.predict(test_data)
prob = np.max(classifier.predict_proba(test_data), axis = 1)
unlabeled = np.where(prob < 0.1)
predicted[unlabeled] = 'Unknown'

pred = pd.DataFrame(predicted)
pred.to_csv("results/" + "RF_01.csv", index = False)



