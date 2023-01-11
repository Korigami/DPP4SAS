                                                                                                                                       
filename temp2 catalog 'sashelp.dpp.stationary_1_d.source';                                                                             
%include temp2;                                                                                                                  
filename temp2;                                                                                                                  
                                                                                                              
                                                                                                                                 
%macro train;                                                                                                                            
                                                                                                                                    
                                                                                                                                    
    %em_getname(key=Sample, type=data);                                                                                                  
          
    %stat_1_dep_sampler(
        size=&EM_PROPERTY_SIZE.,
        dest=&EM_USER_SAMPLE.,
        mode=&EM_PROPERTY_MODE., 
        base=&EM_PROPERTY_BASE., 
        x0=&EM_PROPERTY_X0., 
        random_state=&EM_PROPERTY_RANDOM_STATE.);

    data &EM_EXPORT_TRAIN;                                                                                                               
        set &EM_USER_SAMPLE.;                                                                                                          
    run;                                                                                                                                 
                                                                                                                                 
%mend train;