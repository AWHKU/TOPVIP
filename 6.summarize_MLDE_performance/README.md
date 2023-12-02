## Evaluate MLDE performance ##
The R scripts took in the predictedFitness.csv from MLDE runs and calculate the performance metrics NDCG, Spearman's *rho*, enrichment, selectivity, top 50 identified, top1% identified, etc. used in the paper. 


Use the MLDE output in 

data/execute_MLDE_benchmark

data/execute_MLDE_multi_round

data/MLDE_train_test_datasets

as input for 231201_evaluate_MLE_performance.R. 
