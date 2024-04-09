# contamination.and.boostedGS
Code for contamination experiments and boosted versions of Standard Gibbs Sampling 
This code goes with the paper by García-Donato and Castellanos (2024)

How to use:
Use the sentences in the file ExampleALLmethods-vc.R to obtain exact inclusion probabilities with different values of p for a (virtually) contaminated experiment (n=500, kt=4) and their estimation using PARNI, ASI, TGS, wTGS and SGS.

Use the sentences in the file ExampleALLmethods-sc.R to obtain estimation of inclusion probabilities with different values of p for a contaminated -via simulation- experiment (n=500, kt=4) and their estimation using PARNI, ASI, TGS, wTGS and SGS.

Since the virtually contaminated experiment "emulates" the simulated experiment, you should expect see similar results for the inclusion probabilities.

Acknowledgements:
In several places we are using functions coming from the github repositories (accessed in 2023) [1] XitongLiang/The-PARNI-scheme and [2] gZanella/TGS that, accompany the papers Liang, Livingstone, Griffin (2021) and Zanella and Roberts (2018) respectively. In particular ASI.R and PARNI.R are essentially the files with the same names in [1]; other_supportive_functions.R is a light version of the file with the same name in [1]; make_hyper_par.R is a modified version of the file with the same name in [1] but adapted to include our contributions; functions_for_BVS.R is a modified version of the file with the same name in [2] adapted to make it possible computation of TGS and wTGS handling g-Zellner prior in a fast and more reliable way; bf_g_L_4GS.R is an adaptation of the file with the name bf_g_L.R in [1] to be run with the vectorized version of Gibbs Sampling. All the other functions are built from scratch. 
