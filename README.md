# An evaluation of automated methods for cell type annotation in scRNA-seq data

This repository contains the code for the analyses developed in my thesis.

- **code_QC_mca.py** and **code_QC_pbmc.R** files contain the code to perform a quality control of MCA and PBMCs datasets, respectively.
- **code_processing.R** file contains the code for the processing of the test data (MCA and PBMCs) and reference data (ImmGen and Monaco).
- **code_correlation_matrix.R** file contains the code for obtaining the cell types correlation matrix of MCA and PBMCs dataset.
- **annotation** folder contains different files with the code for obtaining SingleR, SVM, SVMrejection and scType annotations. These annotations were analysed in Section 3.1 and 3.2.
- **evaluation.R** file contains the code for evaluating the SingleR, SVM, SVMrejection and scType annotations, giving the different tables and figures shown in Section 3.1 and 3.2. 
- **JArribasComparison** folder contains the processing of JArribas dataset, the annotations obtained from SingleR and SVM, and the comparison that was explained in Section 3.3.
