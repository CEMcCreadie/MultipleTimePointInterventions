# Single Time Point Interventions Tests from original paper.
# Nota bene: Simulation 1b in our Thesis refers to Simulation 1c in original paper.


source("tstmle.R")
library(rje)
library(data.table)
library(glue)

N <- 100
monte_carlo_draws <- 500
folder <- "/home/simulation_data"

simulation_1a_true_ATE <- 0.218
simulation_1b_true_ATE <- 0.279

simulation_1a <- function(){
    # Using rbinom(_, 1, _) seems easiest way to get the behaviour of a bernoulli trial 
    W <- data.frame(
        "W1" = rbinom(N, 1, 0.5),
        "W2" = ceiling(runif(N, 0, 3)), # Equivalent to uniform 1 to 3
        "W3" = rbinom(N, 1, 0.5)
    )

    A <- numeric(N)
    Y <- numeric(N)
    A[1:4] <- rbinom(4, 1, 0.5)
    Y[1:4] <- rbinom(4, 1, 0.5)


    # W[i, j] <- i'th row, j'th column
    for(i in 5:N){
        p <- expit(0.25*W[(i-1), 1] - 0.2*W[(i-1), 2] + 0.3*Y[(i-1)] 
            - 0.2*A[(i-1)] + 0.2*W[(i-2), 3]
        )
        A[i] <- rbinom(1, 1, p)
        p <- expit(0.3 - 0.8*W[(i-1), 1] + 0.1*W[(i-1), 2] + 0.2*W[(i-1),3] 
            + A[i] - 0.5*W[(i-2), 1] + 0.2*W[(i-2), 3]
        )
        Y[i] <- rbinom(1, 1, p)
    }

    data <- cbind(W,A,Y)
    return(data)
}

file <- glue("{folder}/simulation1a_MC{monte_carlo_draws}_N{N}.rds")

node_list <- list(
    W = c("W1", "W2", "W3"),
    A = "A",
    Y = "Y"
)
ate_spec <- tmle_ATE(treatment_level = 1, control_level=0)

tmle_fits <- data.table()

# for (i in 1:monte_carlo_draws){
#     draw <- simulation_1a()
#     tmle_fit <- run_tstmle(ate_spec, draw, node_list, markov_order=3)
#     tmle_fits <- rbind(tmle_fits, tmle_fit$summary)
# }
# saveRDS(tmle_fits, file)


simulation_1b <- function(){
    # Using rbinom(_, 1, _) seems easiest way to get the behaviour of a bernoulli trial 
    W1 <- numeric(N)
    W2 <- numeric(N)

    W1[1:6] <- rbinom(6, 1, 0.5)
    W2[1:6] <- rnorm(6, 0, 1)

    A <- numeric(N)
    Y <- numeric(N)
    A[1:6] <- rbinom(6, 1, 0.5)
    Y[1:6] <- rbinom(6, 1, 0.5)

    diff <- numeric(N)
    # W[i, j] <- i'th row, j'th column
    for(i in 7:N){
        p <- expit(0.5*W1[(i-1)] - 0.5*Y[(i-1)] + 0.1*W2[i-1])
        W1[i] <- rbinom(1, 1, p)

        mu <- 0.6 * A[i-1] + Y[i-1] - W1[i-1]
        W2[i] <- rnorm(1, mu, 1)

        p <- expit( 0.7 * W1[i-2] - 0.3*A[i-1] + 0.2 * sin(W1[i-2] * A[i-3]))
        A[i] <- rbinom(1, 1, p)

        p <- expit(1.5 * A[i] - (W1[i-1] * A[i-2])^2 + 0.9 * sin(W2[i-4])*A[i-3]*cos(W2[i-6]) - abs(W2[i-5]))

        Y[i] <- rbinom(1, 1, p)
    }

    data <- data.frame(W1,W2,A,Y)
    return(data)
}

file <- glue("{folder}/simulation1c_MC{monte_carlo_draws}_N{N}.rds")

node_list <- list(
    W = c("W1", "W2"),
    A = "A",
    Y = "Y"
)
ate_spec <- tmle_ATE(treatment_level = 1, control_level=0)

tmle_fits <- data.table()

for (i in 1:monte_carlo_draws){
    draw <- simulation_1b()
    tmle_fit <- run_tstmle(ate_spec, draw, node_list, markov_order=6)
    tmle_fits <- rbind(tmle_fits, tmle_fit$summary)
}
saveRDS(tmle_fits, file)