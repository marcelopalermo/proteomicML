# These are the files used for this paper.
The ANOVA output is a mandatory step from Olink proteomics analysis (included in OlinkAnalyze R package)
In this case, the original data (Prep_CSV.csv) will produce one ANOVA analysis per top-10 protein, because it analyses 50 patients in day0 and day14 (where protein differntial expression NPX is exracted)

Our R routine (Olink_Population_Generator.R) creates 500 additional fictional patients (1000 rows, as 500 fictional patients will contain day0 and day14 proteomic measurements). It also creates the ANOVA data sets explained in sections II-A through II-E within the paper.
The ML model will be trained with anova_ML_Training_class.csv. Although this csv contains 1080 rows (see section II-E), each ML model will receive 360 rows, by representing 2 correlated categories each (Age>=50 and Age<50, BMI>=30 and BMI<30, and gender=male or female)
File names in this Git:
- Accuracy_LogLoss.chart.png: 
