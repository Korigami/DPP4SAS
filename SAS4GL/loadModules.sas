%macro load_modules(dirPath);

filename exact "&dirPath.\DPP4SAS\SASMINER\DPP\ExactSampler_source.sas";
filename mcmc "&dirPath.\DPP4SAS\SASMINER\DPP\MCMCSampler_source.sas";
filename beta "&dirPath.\DPP4SAS\SASMINER\DPP\random_matrices.sas";
filename poison "&dirPath.\DPP4SAS\SASMINER\DPP\poisson_final.sas";
filename stat_1_d "&dirPath.\DPP4SAS\SASMINER\DPP\stationary_1_d.sas";
filename utils "&dirPath.\DPP4SAS\SASMINER\DPP\utils.sas";
%include exact;  
%include mcmc; 
%include beta; 
%include poison; 
%include stat_1_d; 
%include utils;
%mend;

%load_modules(dirPath);
