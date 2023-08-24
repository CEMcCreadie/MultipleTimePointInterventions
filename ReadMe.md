The implementation of <a href="https://arxiv.org/pdf/1809.00734.pdf">Robust Estimation of Data-Dependent Causal Effects based on Observing a Single Time-Series</a> (Ivanka Malenica and Mark van der Laan) from my MSc Thesis.

Performs TMLE estimate through sequential regression (5.3) of ATE over multiple interventions in a time series.

Implementation of Single Time Point inspired from <a href="https://github.com/imalenica/tstmle3">tstmle3</a>, Multiple Time Point implementation original.

<em>Deprecated</em> notebooks included for exposition of code. 

<h2>Example Use</h2>
Runs TMLE estimate over multiple time series data points in a folder and saves all estimates to a .csv file.

```
library(glue)
source("multi_interventions.R")

scratch_dir <- "/users/john/documents/TMLE"
experiment_folder <- "experimentOne"
files <- list.files(path = glue("{scratch_dir}/{experiment_folder}"), recursive=TRUE)
print(files)
K <- 5 # Number of Time Point Interventions

tmle_fits <- matrix(nrow=0, ncol = 3)
for (file in files) {
    data <- preprocess_data(read.csv(glue("{scratch_dir}/{experiment_folder}/{file}")))

    node_list <- list(
        W = c("activity_type"),
        A = "pm2_5",
        Y = "target"
    ) # See tlverse tutorial Chapter 7

    result <- run_multi_interventions_tmle(data, node_list, K)
    subject <- strsplit(file, ".csv")[[1]][1]

    tmle_fits <- rbind(tmle_fits, c(subject, result))

    output_file <- glue("{scratch_dir}/{experiment_folder}/result{K}.csv")
    write.csv(tmle_fits, output_file, row.names=FALSE)
}
```

<h2>Packages Required</h2>

- data.table
- glue
- hal9001
- Hmisc
- origami
- ranger
- rje
- sl3
- tmle3
- xgboost
