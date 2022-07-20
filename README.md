# These are the files used for this paper.
The ANOVA output is a mandatory step from Olink proteomics analysis (included in OlinkAnalyze R package)
In this case, the original data (Prep_CSV.csv) will produce one ANOVA analysis per top-10 protein, because it analyses 50 patients in day0 and day14 (where protein differntial expression NPX is exracted)
Our R routine creates 17 additional fictional patients to produce 17 ANOVA analysis per top-10 proteins, so that we will have enough data to enter the ML models
The ML model will be trained on the top-10 proteins to be able to predict potential patients that would present biomarker for such top-10 proteins as potential sequelae/long-COVID patients
