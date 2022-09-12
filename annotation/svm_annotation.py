
# Annotation - SVM with linear kernel #
#######################################

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
                      change_features = False):
    
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
    
    
    # Train and prediction with lineal kernel without rejection option
    if rejection_option == False:
        
        Classifier = LinearSVC(max_iter=10000)
        
        start=tm.time()
        Classifier.fit(train_data, np.ravel(train_labels))
        tr_time = tm.time()-start
        
        start = tm.time()
        predicted = Classifier.predict(test_data)
        ts_time = tm.time() - start
        
    if rejection_option == True:
        
        Classifier = LinearSVC(max_iter=100000)
   
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
    prediction_time["prob"] = prob
    prediction_time['prediction'] = pd.DataFrame(predicted)
    prediction_time['time']   = pd.DataFrame({'Training time': [ tr_time],
                         'Test time': [ts_time]})
    
    return(prediction_time)
    


##############################
# RUN SVM with linear kernel # 
##############################
   
    
# MCA - ImmGen
###############
os.chdir("C:/Users/Usuario/Desktop/master/tfm/HCL_data")

# All cells:
    
    # Load the test data
test_data = pd.read_csv("mca_processed_all.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("refImmScaled.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_refImm.csv", header=0,index_col=None, sep='|')	

    # Load prediction and time
prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = False)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_linear_ImmGen_scaled_All.csv", index = False)
time.to_csv("results/" + "SVM_linear_Time_ImmGen_scaled_All.csv", index = False)



# Select cells

test_data = pd.read_csv("mca_processed_select.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("refImmScaledSelect.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_refImm_select.csv", header=0,index_col=None, sep='|')	

    # Load prediction and time
prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = False)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_linear_ImmGen_scaled_Select.csv", index = False)
time.to_csv("results/" + "SVM_linear_Time_ImmGen_scaled_Select.csv", index = False)





# PBMC - Monaco
################

os.chdir("C:/Users/Usuario/Desktop/master/tfm/PBMC_data")

# All cells

    # Load the test data
test_data = pd.read_csv("pbmc_processed_all.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("refMonacoScaled.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_refMonaco.csv", header=0,index_col=None, sep='|')	

    # Load prediction and time
prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = False)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_linear_Monaco_scaled_All.csv", index = False)
time.to_csv("results/" + "SVM_linear_Time_Monaco_scaled_All.csv", index = False)



# Select cells

    # Load the test
test_data = pd.read_csv("pbmc_processed_select.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("refMonacoScaledSelect.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_refMonaco_select.csv", header=0,index_col=None, sep='|')	

    # Load prediction and time
prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = False)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_linear_Monaco_scaled_Select.csv", index = False)
time.to_csv("results/" + "SVM_linear_Time_Monaco_scaled_Select.csv", index = False)



# PBMC - ImmGen
################

# All cells

    # Load the test data
test_data = pd.read_csv("pbmc_processed_all.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("refImmScaled.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_refImm.csv", header=0,index_col=None, sep='|')	

    # Load prediction and time
prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = False, change_features = True)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_linear_Imm_scaled_All.csv", index = False)
time.to_csv("results/" + "SVM_linear_Time_Imm_scaled_All.csv", index = False)



# Select cells

    # Load the test
test_data = pd.read_csv("pbmc_processed_select.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("refImmScaledSelect.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_refImm_select.csv", header=0,index_col=None, sep='|')	

    # Load prediction and time
prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = False, change_features = True)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_linear_Imm_scaled_Select.csv", index = False)
time.to_csv("results/" + "SVM_linear_Time_Imm_scaled_Select.csv", index = False)





####################################################
# RUN SVM with linear kernel with rejection option # 
####################################################


# MCA - ImmGen
###############
os.chdir("C:/Users/Usuario/Desktop/master/tfm/HCL_data")

# All cells

test_data = pd.read_csv("mca_processed_all.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("refImmScaled.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_refImm.csv", header=0,index_col=None, sep='|')	


train_labels["x"].values.tolist().index("B cells, pro")

train_data.drop(train_data.columns[428], 
        axis=1, 
        inplace=True)

train_labels.drop(index=train_labels.index[428], 
        axis=0, 
        inplace=True)


train_labels["x"].values.tolist().index("Microglia")


train_data.drop(train_data.columns[774], 
        axis=1, 
        inplace=True)

train_labels.drop(index=train_labels.index[774], 
        axis=0, 
        inplace=True)


train_data.drop(train_data.columns[774], 
        axis=1, 
        inplace=True)

train_labels.drop(index=train_labels.index[774], 
        axis=0, 
        inplace=True)


train_data.drop(train_data.columns[774], 
        axis=1, 
        inplace=True)

train_labels.drop(index=train_labels.index[774], 
        axis=0, 
        inplace=True)


train_labels["x"].values.tolist().index("Eosinophils")


train_data.drop(train_data.columns[25], 
        axis=1, 
        inplace=True)

train_labels.drop(index=train_labels.index[25], 
        axis=0, 
        inplace=True)


train_data.drop(train_data.columns[25], 
        axis=1, 
        inplace=True)

train_labels.drop(index=train_labels.index[25], 
        axis=0, 
        inplace=True)



train_data.drop(train_data.columns[25], 
        axis=1, 
        inplace=True)

train_labels.drop(index=train_labels.index[25], 
        axis=0, 
        inplace=True)


train_data.drop(train_data.columns[25], 
        axis=1, 
        inplace=True)

train_labels.drop(index=train_labels.index[25], 
        axis=0, 
        inplace=True)





    # Load prediction and time
prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = True, n_folds = 5, probab = 0.7)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_linear_rej_ImmGen_scaled_All.csv", index = False)
time.to_csv("results/" + "SVM_linear_rej_Time_ImmGen_scaled_All.csv", index = False)




os.chdir("C:/Users/Usuario/Desktop/master/tfm/HCL_data")

# All cells

test_data = pd.read_csv("mca_processed_all.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("ImmGen_Downsampled.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_ImmGen_Downsampled.csv", header=0,index_col=None, sep='|')	


    # Load prediction and time
prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = True, n_folds = 5, probab = 0.7)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_linear_rej_ImmGen_downsampled.csv", index = False)
time.to_csv("results/" + "SVM_linear_rej_Time_ImmGen_downsampled.csv", index = False)





# PBMC - Monaco
################

os.chdir("C:/Users/Usuario/Desktop/master/tfm/PBMC_data")


# All cells

test_data = pd.read_csv("pbmc_processed_all.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("refMonacoScaled.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_refMonaco.csv", header=0,index_col=None, sep='|')	


prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = True, n_folds = 5, 
                                    probab = 0.7, change_features = True)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_linear_rej_Monaco_scaled_All.csv", index = False)
time.to_csv("results/" + "SVM_linear_rej_Time_Monaco_scaled_All.csv", index = False)




# PBMC - ImmGen
################


# All cells

    # Load the test data
test_data = pd.read_csv("pbmc_processed_all.csv",index_col=0,sep=',')

    # Load the train data
train_data = pd.read_csv("refImmScaled.csv" ,index_col=0,sep=',')
train_labels = pd.read_csv("labels_refImm.csv", header=0,index_col=None, sep='|')	


train_labels["x"].values.tolist().index("B cells, pro")

train_data.drop(train_data.columns[428], 
        axis=1, 
        inplace=True)

train_labels.drop(index=train_labels.index[428], 
        axis=0, 
        inplace=True)


train_labels["x"].values.tolist().index("Microglia")


train_data.drop(train_data.columns[774], 
        axis=1, 
        inplace=True)

train_labels.drop(index=train_labels.index[774], 
        axis=0, 
        inplace=True)


train_data.drop(train_data.columns[774], 
        axis=1, 
        inplace=True)

train_labels.drop(index=train_labels.index[774], 
        axis=0, 
        inplace=True)


train_data.drop(train_data.columns[774], 
        axis=1, 
        inplace=True)

train_labels.drop(index=train_labels.index[774], 
        axis=0, 
        inplace=True)


train_labels["x"].values.tolist().index("Eosinophils")


train_data.drop(train_data.columns[25], 
        axis=1, 
        inplace=True)

train_labels.drop(index=train_labels.index[25], 
        axis=0, 
        inplace=True)


train_data.drop(train_data.columns[25], 
        axis=1, 
        inplace=True)

train_labels.drop(index=train_labels.index[25], 
        axis=0, 
        inplace=True)



train_data.drop(train_data.columns[25], 
        axis=1, 
        inplace=True)

train_labels.drop(index=train_labels.index[25], 
        axis=0, 
        inplace=True)


train_data.drop(train_data.columns[25], 
        axis=1, 
        inplace=True)

train_labels.drop(index=train_labels.index[25], 
        axis=0, 
        inplace=True)


    # Load prediction and time
prediction_time = SVM_linear_kernel(train_data, train_labels, test_data, 
                                    rejection_option = True, n_folds = 5, 
                                    probab = 0.7, change_features = True)

pred = prediction_time["prediction"]
time = prediction_time["time"]


    # Export results
pred.to_csv("results/" + "SVM_linear_rej_Imm_scaled_All.csv", index = False)
time.to_csv("results/" + "SVM_linear_rej_Time_Imm_scaled_All.csv", index = False)
