The implementation of <a href="https://arxiv.org/pdf/1809.00734.pdf">Robust Estimation of Data-Dependent Causal Effects based on Observing a Single Time-Series</a> (Ivanka Malenica and Mark van der Laan) from my MSc Thesis.

Performs TMLE estimate through sequential regression (5.3) of ATE over multiple interventions in a time series.

Implementation of Single Time Point inspired from <a href="https://github.com/imalenica/tstmle3">tstmle3</a>, Multiple Time Point implementation orginal.

<em>Deprecated</em> notebooks included for exposition of code. 

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