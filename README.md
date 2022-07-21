# These are the files used for this paper.
The ANOVA output is a mandatory step from Olink proteomics analysis (included in OlinkAnalyze R package)

In this case, the original data (Prep_CSV.csv) will produce one ANOVA analysis per top-10 protein, because it analyses 50 patients in day0 and day14 (where protein differntial expression NPX is exracted)

Our R routine (Olink_Population_Generator.R) creates 500 additional fictional patients (1000 rows, as 500 fictional patients will contain day0 and day14 proteomic measurements). It also creates the ANOVA data sets explained in sections II-A through II-E within the paper.

The ML model will be trained with anova_ML_Training_class.csv. Although this csv contains 1080 rows (see section II-E), each ML model will receive 360 rows, by representing 2 correlated categories each (Age>=50 and Age<50, BMI>=30 and BMI<30, and gender=male or female)

File names in this Git:
- Accuracy_LogLoss.chart.png: ML classifier results, per accuracy score and logarithmic loss
- Accuracy_Results.xlsx: Excel chart that generates Accuracy_LogLoss.chart.png graphics
- Accuracy_chart.png: ML classifier results, per accuracy score, only
- LogLoss_chart.png: ML classifier results, per logarithmic loss, only
- Onlink_ML_Age_KNN_All.ipynb: Jupyter notebook Python 3.9 ML tests for age
- Olink_ML_BMI_KNN_All.ipynb: Jupyter notebook Python 3.9 ML tests for BMI
- Olink_ML_Gender_KNN_All.ipynb: Jupyter notebook Python 3.9 ML tests for gender
- Olink_Population_Generator.R: R program to process OlinkAnalyze package, ANOVA and data sets to feed the ML models in Python Jupyter notebooks
- Prep_CSV.csv: CSV input data used by Olink_Population_Generator.R, based on Zhong et al. (2021)
- Prep_CSV.xlsx: Excel chart that produces Prep_CSV.csv
- anova_ML_Training_class.csv: input file that will feed the Python ML models in Jupyter notebooks

#Steps to reproduce this paper:
1. Install R language and R Studio (R language version used in this paper: 4.2.1)
