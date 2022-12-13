                                                                                                                                         
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
          
    
    %sample_from_beta_ensemble_banded(
        result_eigvals=&EM_USER_SAMPLE.,     
        ensemble_version=&EM_PROPERTY_ENSEMBLEVERSION., 
        size=&EM_PROPERTY_SIZE.,
        beta=&EM_PROPERTY_BETA.,    
        /* distribution specific parameters */
        loc=&EM_PROPERTY_LOC.,  
        scale=&EM_PROPERTY_SCALE.,  
        shape=&EM_PROPERTY_SHAPE.,  
        a=&EM_PROPERTY_A.,  
        b=&EM_PROPERTY_B.,  
        normalize=%boolean_decoder(&EM_PROPERTY_NORMALIZE.),
        heurestic_fix=%boolean_decoder(&EM_PROPERTY_HEURESTIC_FIX.), 
        random_state=&EM_PROPERTY_RANDOM_STATE.     
        );

    data &EM_EXPORT_TRAIN;                                                                                                               
        set &em_user_Sample.;                                                                                                          
    run;                                                                                                                                 
                                                                                                                                 
%mend train;
