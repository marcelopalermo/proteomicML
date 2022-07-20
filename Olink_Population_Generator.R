#install.packages("remotes")
# remotes::install_github(repo ='Olink-Proteomics/OlinkRPackage/OlinkAnalyze', ref = "main", build_vignettes = TRUE)

#install.packages("OlinkAnalyze")

# open package
library(OlinkAnalyze)
#install.packages("dplyr")
library(dplyr)

getwd()
setwd("/Users/marcelo_palermo/Doutorado/R/Olink")
getwd()

# para manter os mesmos números aleatórios
set.seed(212)

# 4.2 /Library/Frameworks/R.framework/Versions/4.2/Resources/library/OlinkAnalyze/extdata
# COPIAR O ARQUIVO: exdata vem do caminho /Users/marcelo_palermo/Library/R/x86_64/4.1/library/OlinkAnalyze/extdata
#file <- system.file("extdata", "Example_NPX_Data.csv", package = "OlinkAnalyze")

file <- system.file("extdata", "Prep_CSV.csv", package = "OlinkAnalyze")

#file <- system.file("extdata", "Olink_NPX_COVID_Sample.csv", package = "OlinkAnalyze")


# Este é o diretório onde fica o arquivo
print(file)

NPX_csv <- read_NPX(file)

print(NPX_csv)

#typeof(NPX_csv)
#typeof(npx_data1)

#Tendencia correta crescente conforme cada proteina
#Significant -> Dados de treinamento de entrada
anova_results_time <-
  olink_anova(df = as.data.frame(NPX_csv), variable = 'Time')
anova_results_real_data <- anova_results_time
#--OR--
#nova_results_time <- olink_anova(df = NPX_csv,variable = 'Time')
write.csv(anova_results_time, "anova_article_ranking.csv", row.names = TRUE)


#Artigo
anova_results_article <- olink_anova(df = as.data.frame(NPX_csv),
                                     variable = 'Time',
                                     covariates = 'Gender:Age:BMI')

#Tendencia correta crescente conforme cada proteina
#Significant

#Non-Significant -> Um dos dados de teste de entrada
anova_results_subject <-
  olink_anova(df = NPX_csv, variable = 'Subject')

#Filtros se necessario em SCARB2
NPX_csv_df = as.data.frame(NPX_csv)
NPX_csv_df_SCARB2 <- NPX_csv_df[NPX_csv_df$Assay %in% c("SCARB2"), ]

NPX_csv_df_female <- NPX_csv_df[which(NPX_csv_df$Gender == 'f'),]
#NPX_csv_df_female <- NPX_csv_df[which(NPX_csv_df$Gender=='f' & NPX_csv_df$Assay=='SCARB2'), ]

NPX_csv_df_male <- NPX_csv_df[which(NPX_csv_df$Gender == 'm'),]
#NPX_csv_df_male <- NPX_csv_df[which(NPX_csv_df$Gender=='m' & NPX_csv_df$Assay=='SCARB2'), ]


NPX_csv_df_old <- NPX_csv_df[which(NPX_csv_df$Age >= 50),]
#NPX_csv_df_old <- NPX_csv_df[which(NPX_csv_df$Age>=60 & NPX_csv_df$Assay=='SCARB2'), ]

NPX_csv_df_young <- NPX_csv_df[which(NPX_csv_df$Age < 50),]
#NPX_csv_df_young <- NPX_csv_df[which(NPX_csv_df$Age<60 & NPX_csv_df$Assay=='SCARB2'), ]


NPX_csv_df_high_BMI <- NPX_csv_df[which(NPX_csv_df$BMI > 30),]
NPX_csv_df_good_BMI <- NPX_csv_df[which(NPX_csv_df$BMI <= 30),]
NPX_csv_df_test <-
  NPX_csv_df[which(NPX_csv_df$Age >= 50 & NPX_csv_df$BMI > 30),]
NPX_csv_df_test2 <-
  NPX_csv_df[which(NPX_csv_df$Gender == 'f' &
                     NPX_csv_df$BMI > 30),]
NPX_csv_df_test3 <-
  NPX_csv_df[which(NPX_csv_df$Gender == 'm' &
                     NPX_csv_df$BMI > 30),]

#FEMALE - Tendencia correta crescente conforme cada proteina
#Significant
anova_results_time_female <-
  olink_anova(df = NPX_csv_df_female, variable = 'Time')
#MALE - Tendencia correta crescente conforme cada proteina
#Significant
anova_results_time_male <-
  olink_anova(df = NPX_csv_df_male, variable = 'Time')
#OLD- Tendencia correta crescente conforme cada proteina
#Significant
anova_results_time_old <-
  olink_anova(df = NPX_csv_df_old, variable = 'Time')
#YOUNG - Tendencia correta crescente conforme cada proteina
#Significant
anova_results_time_young <-
  olink_anova(df = NPX_csv_df_young, variable = 'Time')
#High BMI - Tendencia correta crescente conforme cada proteina
#Significant
anova_results_time_high_BMI <-
  olink_anova(df = NPX_csv_df_high_BMI, variable = 'Time')
#Good BMI - Tendencia correta crescente conforme cada proteina
#Significant
anova_results_time_high_BMI <-
  olink_anova(df = NPX_csv_df_good_BMI, variable = 'Time')
#Gender=='f' & NPX_csv_df$BMI>30
anova_results_time_BMI_high <-
  olink_anova(df = NPX_csv_df_test2, variable = 'Time')
#Gender=='m' & NPX_csv_df$BMI>30
anova_results_time_BMI_low <-
  olink_anova(df = NPX_csv_df_test3, variable = 'Time')


#Ajuda em /Users/marcelo_palermo/Downloads/OlinkRPackage-main/OlinkAnalyze/vignettes/Vignett.Rmd

##Gerador de valores randômicos
Max_Min_Rand <- function(df, protein) {
  options(digits = 9)
  #randomicos day0
  df_filtered0 = df[which(NPX_csv_df$Assay == protein &
                            NPX_csv_df$Time == 'day0'),]
  # lim_day0 <- runif(50, min=min(df_filtered0$NPX), max=max(df_filtered0$NPX))
  lim_day0 <- df_filtered0$NPX * runif(1, 0.95, 1.05)
  
  #substituir no dataframe
  df_filtered0$NPX = c(lim_day0)
  df_filtered0$SampleID <- gsub('A', 'FICA', df_filtered0$SampleID)
  
  #randomicos day14
  df_filtered14 = df[which(NPX_csv_df$Assay == protein &
                             NPX_csv_df$Time == 'day14'),]
  #lim_day14 <- runif(50, min=min(df_filtered14$NPX), max=max(df_filtered14$NPX))
  lim_day14 <- df_filtered14$NPX * runif(1, 0.95, 1.05)
  
  #substituir no dataframe
  df_filtered14$NPX = c(lim_day14)
  df_filtered14$SampleID <- gsub('B', 'FICB', df_filtered14$SampleID)
  #retornar o dataframe
  
  
  return(rbind(df_filtered0, df_filtered14))
  
}

#Data Frame vazio
anova_results_test = anova_results_time[0,]
anova_results_random_old =  anova_results_time_old[0,]
anova_results_random_young =  anova_results_time_young[0, ]
anova_results_random_high_BMI =  anova_results_time_BMI_high[0, ]
anova_results_random_low_BMI =  anova_results_time_BMI_low[0, ]
anova_results_random_female =  anova_results_time_female[0, ]
anova_results_random_male =  anova_results_time_male[0, ]

for (j in 1:17) {
  for (i in 1:10) {
    #First 10 random anova rows for each category
    NPX_csv_random = Max_Min_Rand(NPX_csv_df, 'SCARB2')
    #anova_results_random <-
    # olink_anova(df = as.data.frame(NPX_csv_random), variable = 'Time')
    NPX_csv_random = merge(NPX_csv_random, Max_Min_Rand(NPX_csv_df, 'SIGLEC1'), all = TRUE)
    NPX_csv_random = merge(NPX_csv_random, Max_Min_Rand(NPX_csv_df, 'CTSO'), all = TRUE)
    NPX_csv_random = merge(NPX_csv_random, Max_Min_Rand(NPX_csv_df, 'CXCL10'), all = TRUE)
    NPX_csv_random = merge(NPX_csv_random, Max_Min_Rand(NPX_csv_df, 'LAG3'), all = TRUE)
    NPX_csv_random = merge(NPX_csv_random, Max_Min_Rand(NPX_csv_df, 'CCL8'), all = TRUE)
    NPX_csv_random = merge(NPX_csv_random, Max_Min_Rand(NPX_csv_df, 'GRN'), all = TRUE)
    NPX_csv_random = merge(NPX_csv_random, Max_Min_Rand(NPX_csv_df, 'IFNL1'), all = TRUE)
    NPX_csv_random = merge(NPX_csv_random, Max_Min_Rand(NPX_csv_df, 'LAMP3'), all = TRUE)
    NPX_csv_random = merge(NPX_csv_random, Max_Min_Rand(NPX_csv_df, 'CSF1'), all = TRUE)
    
    #anova_results_random <-
    # olink_anova(df = as.data.frame(NPX_csv_random), variable = 'Time')
    
    #anova_results_test = merge(anova_results_test, anova_results_random, all = TRUE)
    
    
    print(j)
  }
  
  #Gravar CSV para treinamento da rede ML Classificatória
  #write.csv(anova_results_time, "anova_ML_Training.csv", row.names = TRUE)
  
  filename = paste("NPX_csv_random",j,".csv",sep="")
  write.csv(NPX_csv_random, filename, row.names = TRUE)
  
  #Produce Test Data 2x (180 row in anova_results_class)
  NPX_csv_df_random_old <-
    NPX_csv_random[which(NPX_csv_random$Age >= 50),]
  NPX_csv_df_random_young <-
    NPX_csv_random[which(NPX_csv_random$Age < 50),]
  NPX_csv_df_random_high_BMI <-
    NPX_csv_random[which(NPX_csv_random$BMI > 30),]
  NPX_csv_df_random_good_BMI <-
    NPX_csv_random[which(NPX_csv_random$BMI <= 30),]
  NPX_csv_df_random_female <-
    NPX_csv_random[which(NPX_csv_random$Gender == 'f'),]
  NPX_csv_df_random_male <-
    NPX_csv_random[which(NPX_csv_random$Gender == 'm'),]
  
  anova_results_random_old =  merge(
    anova_results_random_old,
    olink_anova(df = NPX_csv_df_random_old, variable = 'Time'),
    all = TRUE
  )
  anova_results_random_young =  merge(
    anova_results_random_young,
    olink_anova(df = NPX_csv_df_random_young, variable = 'Time'),
    all = TRUE
  )
  anova_results_random_high_BMI =  merge(
    anova_results_random_high_BMI,
    olink_anova(df = NPX_csv_df_random_high_BMI, variable = 'Time'),
    all = TRUE
  )
  anova_results_random_low_BMI =  merge(
    anova_results_random_low_BMI,
    olink_anova(df = NPX_csv_df_random_good_BMI, variable = 'Time'),
    all = TRUE
  )
  anova_results_random_female =  merge(
    anova_results_random_female,
    olink_anova(df = NPX_csv_df_random_female, variable = 'Time'),
    all = TRUE
  )
  anova_results_random_male =  merge(
    anova_results_random_male,
    olink_anova(df = NPX_csv_df_random_male, variable = 'Time'),
    all = TRUE
  )
  
  
}


#Join Train + Test Data
#Adding Classifications for significant thresholds
anova_results_time_BMI_high = anova_results_time_BMI_high %>% mutate(class = 'High BMI',type = 'Real')
anova_results_time_BMI_low = anova_results_time_BMI_low %>% mutate(class = 'Low BMI',type = 'Real')
anova_results_time_old = anova_results_time_old %>% mutate(class = '50 above',type = 'Real')
anova_results_time_young = anova_results_time_young %>% mutate(class = '49 below',type = 'Real')
anova_results_time_male = anova_results_time_male %>% mutate(class = 'male',type = 'Real')
anova_results_time_female = anova_results_time_female %>% mutate(class = 'female',type = 'Real')

#Adding Classifications for significant thresholds random
anova_results_random_high_BMI = (as.data.frame(anova_results_random_high_BMI) %>% mutate(class = 'High BMI',type = 'Simulated'))
anova_results_random_low_BMI = (as.data.frame(anova_results_random_low_BMI) %>% mutate(class = 'Low BMI',type = 'Simulated'))
anova_results_random_old = (as.data.frame(anova_results_random_old) %>% mutate(class = '50 above',type = 'Simulated'))
anova_results_random_young = (as.data.frame(anova_results_random_young) %>% mutate(class = '49 below',type = 'Simulated'))
anova_results_random_male = (as.data.frame(anova_results_random_male) %>% mutate(class = 'male',type = 'Simulated'))
anova_results_random_female = (as.data.frame(anova_results_random_female) %>% mutate(class = 'female',type = 'Simulated'))




#put all data frames into list
anova_results_class <-
  list(
    anova_results_time_BMI_high,
    anova_results_time_BMI_low,
    anova_results_time_old,
    anova_results_time_young,
    anova_results_time_male,
    anova_results_time_female,
    anova_results_random_high_BMI,
    anova_results_random_low_BMI,
    anova_results_random_old,
    anova_results_random_young,
    anova_results_random_male,
    anova_results_random_female
  )
#merge all data frames in list
anova_results_class <-
  Reduce(function(x, y)
    merge(x, y, all = TRUE), anova_results_class)


#Gravar CSV para treinamento da rede ML Classificatória
write.csv(anova_results_class,
          "anova_ML_Training_class.csv",
          row.names = TRUE)









# UMAP Plot
NPX.labels = as.factor(NPX_csv_df_random_old[, "Time"])
NPX.data = as.data.frame(NPX_csv_df_random_old[, "NPX"])

#install.packages("umap")
library(umap)
NPX.umap = umap(NPX.data)
head(NPX.umap$layout, 3)



plot.NPX <- function(x,
                     labels,
                     text,
                     main = text,
                     colors = c("#ff7f00", "#e377c2", "#17becf"),
                     pad = 0.1,
                     cex = 0.65,
                     pch = 19,
                     add = FALSE,
                     legend.suffix = "",
                     cex.main = 1,
                     cex.legend = 1) {
  layout = x
  if (is(x, "umap")) {
    layout = x$layout
  }
  
  xylim = range(layout)
  xylim = xylim + ((xylim[2] - xylim[1]) * pad) * c(-0.5, 0.5)
  if (!add) {
    par(mar = c(0.2, 0.7, 1.2, 0.7), ps = 10)
    plot(xylim,
         xylim,
         type = "n",
         axes = F,
         frame = F)
    rect(xylim[1],
         xylim[1],
         xylim[2],
         xylim[2],
         border = "#aaaaaa",
         lwd = 0.25)
  }
  points(layout[, 1],
         layout[, 2],
         col = colors[as.integer(labels)],
         cex = cex,
         pch = pch)
  mtext(side = 3, main, cex = cex.main)
  
  labels.u = unique(labels)
  legend.pos = "topright"
  legend.text = as.character(labels.u)
  if (add) {
    legend.pos = "bottomright"
    legend.text = paste(as.character(labels.u), legend.suffix)
  }
  legend(
    legend.pos,
    legend = legend.text,
    col = colors[as.integer(labels.u)],
    bty = "n",
    pch = pch,
    cex = cex.legend
  )
}

plot.NPX(NPX.umap, NPX.labels, "UMAP Over 50 Random")

NPX.labels = as.factor(NPX_csv_df_old[, "Time"])
NPX.data = as.data.frame(NPX_csv_df_old[, "NPX"])
NPX.umap = umap(NPX.data)
head(NPX.umap$layout, 3)
plot.NPX(NPX.umap, NPX.labels, "UMAP Over 50 Real")



## Variar UMAPS
## Interessates: male e female
## Old e Young
# anova_results_class <- anova_results_class[which(anova_results_class$class == 'male' | anova_results_class$class == 'female'),]
#anova_results_class <- anova_results_class[which(anova_results_class$class != 'male' & anova_results_class$class != 'female'),]


#NPX.labels = as.factor(anova_results_class[,"class"])
#NPX.data = as.data.frame(anova_results_class[,"Adjusted_pval"])
#NPX.umap = umap(NPX.data)
#head(NPX.umap$layout, 3)
#plot.NPX(NPX.umap, NPX.labels, "ANOVA")



