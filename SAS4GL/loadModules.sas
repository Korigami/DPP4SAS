%macro load_modules(dirPath);

filename exact "&dirPath.\SAS4GL\macros\ExactSampler_source.sas";
filename mcmc "&dirPath.\SAS4GL\macros\MCMCSampler_source.sas";
filename beta "&dirPath.\SAS4GL\macros\random_matrices.sas";
filename poison "&dirPath.\SAS4GL\macros\poisson_final.sas";
filename stat_1_d "&dirPath.\SAS4GL\macros\stationary_1_d.sas";
filename utils "&dirPath.\SAS4GL\macros\utils.sas";
%include exact;  
%include mcmc; 
%include beta; 
%include poison; 
%include stat_1_d; 
%include utils;
%mend;
