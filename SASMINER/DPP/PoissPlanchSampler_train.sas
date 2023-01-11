                                                                                                                                         
filename temp2 catalog 'sashelp.dpp.poisson_final.source';                                                                             
%include temp2;                                                                                                                  
filename temp2;                                                                                                                  
                                                                                                              
                                                                                                                                 
%macro train;                                                                                                                            
                                                                                                                                    
                                                                                                                                    
    %em_getname(key=Sample, type=data);                                                                                                  
          
    
    %sample_poissonized_plancherel(
        theta=&EM_PROPERTY_THETA,
        dest=&EM_USER_SAMPLE.,     
        random_state=&EM_PROPERTY_RANDOM_STATE.     
        );

    data &EM_EXPORT_TRAIN;                                                                                                               
        set &EM_USER_SAMPLE.;                                                                                                          
    run;                                                                                                                                 
                                                                                                                                 
%mend train;
