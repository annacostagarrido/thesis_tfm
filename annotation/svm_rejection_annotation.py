
# Annotation of SVM with linear kernel adding or not rejection option


# Packages
import os
import numpy as np
import pandas as pd
import time as tm
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.model_selection import StratifiedKFold


# Function
def SVM_linear_kernel(train_data, train_labels, test_data,
                      rejection_option, n_folds = None, probab = None,
                      seurat_processing = None, train_mouse = None):
    
    # Tranpose genes (columns), cells (rows)
    train_data = train_data.T
    test_data = test_data.T
    
    
    # If contains seurate processing
    if seurat_processing == True:
    
        if train_mouse == True:
            train_data.columns = train_data.columns.str.upper()
        
        if train_mouse == False:
            test_data.columns = test_data.columns.str.upper()
        
        # Select the same features (genes) in both datasets
        train_data = train_data[train_data.columns.intersection(test_data.columns)]
        test_data = test_data[test_data.columns.intersection(train_data.columns)]	
            
        print(train_data.shape)
        print(test_data.shape)
        
        # Features with the same order    
        train_data = train_data.sort_index(axis=1)
        test_data = test_data.sort_index(axis=1)
        
    
    # Train and prediction with lineal kernel and without rejection option
    if rejection_option == False:
        
        Classifier = LinearSVC(max_iter=10000)
        
        start=tm.time()
        Classifier.fit(train_data, np.ravel(train_labels))
        tr_time = tm.time()-start
        
        start = tm.time()
        predicted = Classifier.predict(test_data)
        ts_time = tm.time() - start
    
    # Train and prediction with lineal kernel and rejection option   
    if rejection_option == True:
        
        Classifier = LinearSVC(max_iter=10000)
   
        cv = StratifiedKFold(n_splits= n_folds)
        clf = CalibratedClassifierCV(Classifier, cv = cv)
        
        start = tm.time()
        clf.fit(train_data, np.ravel(train_labels))
        tr_time = tm.time()-start
        
        # prediction
        start = tm.time()
        predicted = clf.predict(test_data)
        prob = np.max(clf.predict_proba(test_data), axis = 1)
        unlabeled = np.where(prob < probab)
        predicted[unlabeled] = 'Unknown'
        ts_time = tm.time() - start
    
    prediction_time = dict(); 
    prediction_time['prediction'] = pd.DataFrame(predicted)
    prediction_time['time']   = pd.DataFrame({'Training time': [ tr_time],
                         'Test time': [ts_time]})
    
    return(prediction_time)



os.chdir("C:/Users/Usuario/Desktop/master/tfm/analysis_2sc")

# RAW COUNTS
##############

# Train SVM with rejection PBMC (train data) and MCA (test data)

    # Load the test data
test_data = pd.read_csv("mca_data.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("pbmc_data.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_pbmc.csv", header=0,index_col=None, sep='|')	

    # Load prediction and time
prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = True, n_folds = 5, probab = 0.7)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_rej70_PBMC_train_MCA_test_notMNN.csv", index = False)
time.to_csv("results/" + "SVM_rej70_Time_PBMC_train_MCA_test_notMNN.csv", index = False)


# Train SVM without rejection PBMC (train data) and MCA (test data)

prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = False)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_PBMC_train_MCA_test_notMNN.csv", index = False)
time.to_csv("results/" + "SVM_Time_PBMC_train_MCA_test_notMNN.csv", index = False)




# Train SVM with rejection MCA (train data) and PBMC (test data)


    # Load the test data
test_data = pd.read_csv("pbmc_data.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("mca_data.csv" ,index_col=0,sep=',')
train_data.rows = train_data.columns.str.upper()
train_labels = pd.read_csv("labels_mca.csv", header=0,index_col=None, sep='|')	

    # Load prediction and time
prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = True, n_folds = 5, 
                                    probab = 0.7)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_rej70_MCA_train_PBMC_test_notMNN.csv", index = False)
time.to_csv("results/" + "SVM_rej70_Time_MCA_train_PBMC_test_notMNN.csv", index = False)



# Train SVM without rejection MCA (train data) and PBMC (test data)

prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = False)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_MCA_train_PBMC_test_notMNN.csv", index = False)
time.to_csv("results/" + "SVM_Time_MCA_train_PBMC_test_notMNN.csv", index = False)



# Seurat processing
#####################

# Train SVM with rejection PBMC (train data) and MCA (test data)

    # Load the test data
test_data = pd.read_csv("mca_processed_all.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("pbmc_processed_all.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_pbmc.csv", header=0,index_col=None, sep='|')	

    # Load prediction and time
prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = True, n_folds = 5, probab = 0.7,
                                    seurat_processing = True, train_mouse = False)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_rej70_PBMC_train_MCA_test_seurat.csv", index = False)
time.to_csv("results/" + "SVM_rej70_Time_PBMC_train_MCA_test_seurat.csv", index = False)


# Train SVM without rejection PBMC (train data) and MCA (test data)

prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = False, 
                                    seurat_processing = True, train_mouse = False)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_PBMC_train_MCA_test_seurat.csv", index = False)
time.to_csv("results/" + "SVM_Time_PBMC_train_MCA_test_seurat.csv", index = False)




# Train SVM with rejection MCA (train data) and PBMC (test data)


    # Load the test data
test_data = pd.read_csv("pbmc_processed_all.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("mca_processed_all.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_mca.csv", header=0,index_col=None, sep='|')	

    # Load prediction and time
prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = True, n_folds = 5, 
                                    probab = 0.7, seurat_processing = True, 
                                    train_mouse = True)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_rej70_MCA_train_PBMC_test_seurat.csv", index = False)
time.to_csv("results/" + "SVM_rej70_Time_MCA_train_PBMC_test_seurat.csv", index = False)



# Train SVM without rejection MCA (train data) and PBMC (test data)

prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = False, 
                                    seurat_processing = True, 
                                    train_mouse = True)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_MCA_train_PBMC_test_seurat.csv", index = False)
time.to_csv("results/" + "SVM_Time_MCA_train_PBMC_test_seurat.csv", index = False)




# MNN COUNTS
##############

# Train SVM with rejection PBMC (train data) and MCA (test data)

    # Load the test data
test_data = pd.read_csv("mca_data_MNN.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("pbmc_data_MNN.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_pbmc.csv", header=0,index_col=None, sep='|')	

    # Load prediction and time
prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = True, n_folds = 5, probab = 0.7)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_rej70_PBMC_train_MCA_test_MNN.csv", index = False)
time.to_csv("results/" + "SVM_rej70_Time_PBMC_train_MCA_test_MNN.csv", index = False)



# Train SVM without rejection PBMC (train data) and MCA (test data)

prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = False)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_PBMC_train_MCA_test_MNN.csv", index = False)
time.to_csv("results/" + "SVM_Time_PBMC_train_MCA_test_MNN.csv", index = False)



# Train SVM with rejection MCA (train data) and PBMC (test data)


    # Load the test data
test_data = pd.read_csv("pbmc_data_MNN.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("mca_data_MNN.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_mca.csv", header=0,index_col=None, sep='|')	

    # Load prediction and time
prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = True, n_folds = 5, prob = 0.7)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_rej70_MCA_train_PBMC_test_MNN.csv", index = False)
time.to_csv("results/" + "SVM_rej70_Time_MCA_train_PBMC_test_MNN.csv", index = False)



# Train SVM without rejection MCA (train data) and PBMC (test data)

prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = False)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_MCA_train_PBMC_test_MNN.csv", index = False)
time.to_csv("results/" + "SVM_Time_MCA_train_PBMC_test_MNN.csv", index = False)

