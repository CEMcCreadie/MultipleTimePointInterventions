library(origami)
library(glue)
library(data.table)
source("tstmle.R")

ATE_SPEC <- tmle_ATE(treatment_level = 1, control_level=0)
TREATMENT_THRESHOLD <- 150

preprocess_data <- function(data){
    data$pm2_5 <- ifelse(data$pm2_5 >= TREATMENT_THRESHOLD, 1, 0)
    max_y <- max(data$breathing_rate)
    data$breathing_rate <- data$breathing_rate / max_y
    data$target <- data$breathing_rate
    data$target <- c(data$target[-1], mean(data$target)) # Shift time step so data is now i.i.d A(t) --> Y(t)
    return (data)
}

int_out_ak <- function(data, node_list, targeted_likelihood, g, k) {
    counterfactual_data <- copy(data)
    treatment <- node_list$"A"
    counterfactual_data$treatment <- 1

    folds <- folds_vfold(nrow(counterfactual_data), V=1)

    counterfactual_task <- make_tmle_task(ATE_SPEC, counterfactual_data, node_list, folds, lag=k)
    counterfactual_vector <- matrix((targeted_likelihood$factor_list[["Y"]]$get_likelihood(counterfactual_task)))

    counterfactual_data$treatment <- 0
    counterfactual_task <- make_tmle_task(ATE_SPEC, counterfactual_data, node_list, folds, lag=k)
    counterfactual_vector_2 <- matrix((targeted_likelihood$factor_list[["Y"]]$get_likelihood(counterfactual_task)))

    #weights <- counterfactual_vector 
    weights <- (counterfactual_vector *  matrix(g)) + (counterfactual_vector_2 * (1-matrix(g)))
    return (weights)
}

run_multi_interventions_tmle <- function(data, node_list, K){

    N <- nrow(data)
    learner_list = get_learner_list()
    se <- 0

    for (k in K:0){
        folds <- origami::make_folds(fold_fun = folds_rolling_window, n=N, window_size = N/5, validation_size=20, gap=10, batch=N/25)
        #folds <- folds_vfold(N, V=20)

        if (k == 0){
            tmle_task <- ATE_SPEC$make_tmle_task(data, node_list, folds)
        } else {
            tmle_task <- make_tmle_task(ATE_SPEC, data, node_list, folds, k)    
        }

        initial_likelihood <- ATE_SPEC$make_initial_likelihood(
            tmle_task,
            learner_list
        )

        g <- initial_likelihood$get_likelihood(tmle_task, node = "A")

        updater <- ATE_SPEC$make_updater(cvtmle=FALSE)
        targeted_likelihood <- ATE_SPEC$make_targeted_likelihood(initial_likelihood, updater)

        tmle_params <- ATE_SPEC$make_params(tmle_task, targeted_likelihood)
        updater$tmle_params <- tmle_params
        tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
        se <- se + tmle_fit$summary$se

        if (k != 0){
            q_nAk <- int_out_ak(data, node_list, targeted_likelihood, g, k) # $\bar Q_{A(t,K)}$
            data$target <- q_nAk
            node_list$Y <- "target"
        }
        
    }

    return(c(tmle_fit$summary$tmle_est, se))
}

