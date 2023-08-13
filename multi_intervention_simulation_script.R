library(glue)
library(rje)
library(data.table)
source("tstmle.R")
source("multi_interventions.R")

N <- 500

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
            + A[i] - 0.5*W[(i-2), 1] + 0.2*W[(i-2), 3] - (0.6*A[i-1]) + (0.2 * A[i-2])
        )
        Y[i] <- rbinom(1, 1, p)
    }

    data <- cbind(W,A,Y)
    return(data)
}

node_list <- list(
    W = c("W1","W2","W3"),
    A = "A",
    Y = "Y"
)

monte_carlo_draws <- 50

tmle_fits <- matrix(nrow=0, ncol = 2)
for (i in 1:monte_carlo_draws){
    sim <- simulation_1a()
    result <- run_multi_interventions_tmle(sim, node_list, 3)
    tmle_fits <- rbind(tmle_fits, result)
}

folder <- "simulation_data"

file <- glue("{folder}/multi_intervention_{monte_carlo_draws}_N{N}.csv")
print(tmle_fits)

write.csv(tmle_fits, file, row.names = FALSE)