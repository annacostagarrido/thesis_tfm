# Annotation - SVM with linear kernel #
#######################################

# Packages
import os
import numpy as np
import pandas as pd
from sklearn.svm import LinearSVC



# Function
def SVM_linear_kernel(train_data, train_labels, test_data):
    
    # PROCESSING
    
    # Tranpose genes (columns), cells (rows)
    train_data = train_data.T
    test_data = test_data.T
    
    
    if change_features == True:
        train_data.columns = train_data.columns.str.upper()
    
    
    # Select the same features (genes) in both datasets
    train_data = train_data[train_data.columns.intersection(test_data.columns)]
    test_data = test_data[test_data.columns.intersection(train_data.columns)]	
        
    print(train_data.shape)
    print(test_data.shape)
    
    # Features with the same order    
    train_data = train_data.sort_index(axis=1)
    test_data = test_data.sort_index(axis=1)
    
    
    # Train and prediction 
    Classifier = LinearSVC(max_iter=10000)
        
    Classifier.fit(train_data, np.ravel(train_labels))
    
    predicted = Classifier.predict(test_data)
        
    
    pred = pd.DataFrame(predicted)

    return(pred)
    


##############################
# RUN SVM with linear kernel # 
##############################
   
os.chdir("C:/Users/Usuario/Desktop/master/tfm/scripts_SingleR_merce")
   
# Setting 1a - 1b
##################

    # Load the test data
test_data = pd.read_csv("JArribas_processed_select.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("refImmScaled.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_specific_refImm.csv", header=0,index_col=None, sep='|')	

    # Load prediction and time
pred = SVM_linear_kernel(train_data, train_labels, test_data)

    # Export results
pred.to_csv("SVM_specific_SI_feature.csv", index = False)


# Setting 2a - 2b
##################

    # Load the test data
test_data = pd.read_csv("JArribas_processed_select.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("refImmScaled.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_main_refImm.csv", header=0,index_col=None, sep='|')	

    # Load prediction and time
pred = SVM_linear_kernel(train_data, train_labels, test_data)

    # Export results
pred.to_csv("SVM_main_SI_feature.csv", index = False)


# Setting 3a - 3b
##################

    # Load the test data
test_data = pd.read_csv("JArribas_processed_all.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("refImmScaled.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_main_refImm.csv", header=0,index_col=None, sep='|')	

    # Load prediction and time
pred = SVM_linear_kernel(train_data, train_labels, test_data)

pred.to_csv("SVM_main_NO_feature.csv", index = False)

