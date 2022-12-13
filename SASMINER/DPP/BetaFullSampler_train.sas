                                                                                                                                         
filename temp2 catalog 'sashelp.dpp.random_matrices.source';                                                                             
%include temp2;                                                                                                                  
filename temp2;                                                                                                                  
                                                                                                                                 
%macro boolean_decoder(encoded);                                                                                                         
    %local decoded;                                                                                                                  
    %if &encoded. = Y %then %let decoded = 1;                                                                                        
    %else %let decoded = 0;                                                                                                          
    &decoded.                                                                                                                        
%mend boolean_decoder;                                                                                                                   
                                                                                                                                 
%macro train;                                                                                                                            
                                                                                                                                    
                                                                                                                                    
    %em_getname(key=Sample, type=data);                                                                                                  
                                                                                                                                    
    %sample_from_beta_ensemble_full(                                                                                                     
        result_eigvals=&EM_USER_SAMPLE.,                                                                                                 
        ensemble_version=&EM_PROPERTY_ENSEMBLEVERSION.,                                                                                  
        M_1=&EM_PROPERTY_M_1., M_2=&EM_PROPERTY_M_2., /* M variables only for Laguerre (first) and Jacobi (both) */                      
        size=&EM_PROPERTY_SIZE.,                                                                                                         
        beta=&EM_PROPERTY_BETA.,                                                                                                         
        normalize=%boolean_decoder(&EM_PROPERTY_NORMALIZE.),                                                                             
        haar_mode=&EM_PROPERTY_HAAR_MODE., /* haar_mode only available for circular ensemble */                                          
        heurestic_fix=%boolean_decoder(&EM_PROPERTY_HEURESTIC_FIX.), /* heurestic_fix only available for circular ensemble */            
        random_state=&EM_PROPERTY_RANDOM_STATE.                                                                                          
    );                                                                                                                                   
                                                                                                                                    
    data &EM_EXPORT_TRAIN;                                                                                                               
    set &em_user_Sample.;                                                                                                          
    run;                                                                                                                                 
                                                                                                                                 
%mend train;
