library(tmle3)
library(sl3)
library(origami)


lag_column = function(column, lag){
    lag_col <- lapply(seq_len(lag), function(x) {
        Hmisc::Lag(column[, 1], x)
    })
    lag_cols <- data.frame(matrix(unlist(lag_col), nrow = length(lag_col[[1]])),
        stringsAsFactors = FALSE)
    return(lag_cols)
}


generate_lag_table = function(W, A, Y, lag){

    lag_Ws <- data.frame(matrix(NA, nrow= nrow(W), ncol=0))
    for (i in 1:ncol(W)){
        colname <- colnames(W)[i]
        lag_col <- lag_column(W[colname], lag)
        names(lag_col) <- paste0(glue("{colname}_lag_"), seq(lag))
        lag_Ws <- cbind.data.frame(lag_Ws, lag_col)
    }

    lag_As = lag_column(A, lag) 
    names(lag_As) <- paste0("A_lag_", seq(lag))

    lag_Ys = lag_column(Y, lag) 
    names(lag_Ys) <- paste0("Y_lag_", seq(lag))

    combined_lag_columns <- cbind.data.frame(lag_Ws, lag_As, lag_Ys)
    combined_lag_columns[is.na(combined_lag_columns)] <- 0

    return(cbind.data.frame(W, A, Y, combined_lag_columns))
}


make_tmle_task = function(tmle_spec, data, node_list, folds, lag){

    data_lag <- generate_lag_table(data.frame(W=data[, node_list$W]),
        A = data.frame(A=data[, node_list$A]),
        Y = data.frame(Y=data[, node_list$Y]),
        lag
    )

    # I think this fixes a bug in the original tstmle3, where data_lag has column Y, 
    # but node_list_lag refers to it as "cnt" in this example
    node_list_lag <- list(
        W = colnames(data_lag)[-(which(names(data_lag)=="A" | names(data_lag)=="Y"))],
        A = "A",
        Y = "Y"
    )

    npsem <- list(
        define_node("W", node_list_lag$W),
        define_node("A", node_list_lag$A, c("W")),
        define_node("Y", node_list_lag$Y, c("A","W"))
    )

    #return(tmle_spec$make_tmle_task(data_lag, node_list_lag, folds=folds))
    return(tmle3_Task$new(data_lag, npsem=npsem, folds=folds))
}


run_tstmle <- function(tmle_spec, data, node_list, learner_list=NULL, markov_order = 5){

    N <- nrow(data)
    # Starting Window Size, Validation Window, Gap between end of Training Set and Validation Set in time units, The increment in Training set size
    folds <- origami::make_folds(fold_fun = folds_rolling_window, n=N, window_size = N/5, validation_size=20, gap=10, batch=N/25)

    tmle_task <- make_tmle_task(tmle_spec, data, node_list, folds, lag=markov_order)
    if(is.null(learner_list)){
        learner_list <- get_learner_list()
    }

    initial_likelihood <- tmle_spec$make_initial_likelihood(
        tmle_task,
        learner_list
    )

    updater <- ate_spec$make_updater(cvtmle=FALSE)
    targeted_likelihood <- ate_spec$make_targeted_likelihood(initial_likelihood, updater)
    tmle_params <- ate_spec$make_params(tmle_task, targeted_likelihood)

    updater$tmle_params <- tmle_params

    tmle_fit <- fit_tmle3(tmle_task, targeted_likelihood, tmle_params, updater)
    return(tmle_fit)
}

get_learner_list <- function(){
    lrn_hal <- make_learner(Lrnr_hal9001)
    lrn_glm <- make_learner(Lrnr_glm)
    lrn_xgboost <- make_learner(Lrnr_xgboost)
    lrn_rf <- make_learner(Lrnr_ranger)
    lrn_ridge <- Lrnr_glmnet$new(alpha = 0)
    lrn_lasso <- Lrnr_glmnet$new(alpha = 1)
    lrn_enet <- Lrnr_glmnet$new(alpha = 0.5)
    #lrn_pois <- Lrnr_glmnet$new(family="poisson")
    #lrn_link <- Lrnr_glmnet$new(family="link")
    lrn_gam <- make_learner(Lrnr_gam)
    #lrn_polspline <- Lrnr_polspline$new()
    #lrn_svm <- make_learner(Lrnr_svm)
    #lrnr_glm_interaction <- Lrnr_glm$new(formula = "~.^2")

    grid_params <- list(
        max_depth = c(3, 8),
        eta = c(0.001, 0.1),
        nrounds = 20
    )
    grid <- expand.grid(grid_params, KEEP.OUT.ATTRS = FALSE)

    xgb_learners <- apply(grid, MARGIN = 1, function(tuning_params) {
        do.call(Lrnr_xgboost$new, as.list(tuning_params))
    })

    learners <- c(list(lrn_glm, lrn_rf, lrn_ridge, lrn_lasso, lrn_hal, lrn_gam, lrn_enet), xgb_learners)

    # define metalearners appropriate to data types
    treatment_metalearner <- make_learner(
        Lrnr_solnp, 
        loss_function = loss_loglik_binomial,
        learner_function = metalearner_logistic_binomial, 
        eval_function = loss_squared_error
    )

    cv_selector <- Lrnr_cv_selector$new(eval_function = loss_squared_error)
    dSL <- Lrnr_sl$new(learners = stack, metalearner = cv_selector)

    outcome_metalearner <- make_learner(
        Lrnr_solnp, 
        # If continous loss_loglik_true_cat and metalearner_linear
        loss_function = loss_squared_error,   
        learner_function =  metalearner_logistic_binomial,
        eval_function = loss_squared_error
    )

    sl_A <- Lrnr_sl$new(
        learners = learners,
        metalearner = treatment_metalearner 
    )

    sl_Y <- Lrnr_sl$new(
        learners = learners,
        metalearner = outcome_metalearner 
    )

    learner_list <- list(A = sl_A, Y = sl_Y)
    return(learner_list)
}