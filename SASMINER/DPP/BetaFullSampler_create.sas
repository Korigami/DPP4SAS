%macro create;
    %em_register(key=Sample, type=data);  
    
    %em_property(name=EnsembleVersion,    value=HERMITE,    action=TRAIN);
    %em_property(name=M_1,    value=10.0,    action=TRAIN);
    %em_property(name=M_2,    value=10.0,     action=TRAIN);
    %em_property(name=Size,    value=10,    action=TRAIN);
    %em_property(name=Beta,    value=2,    action=TRAIN);
    %em_property(name=Normalize,    value=Y,    action=TRAIN);
    %em_property(name=Haar_Mode,    value=HERMITE,    action=TRAIN);
    %em_property(name=Heurestic_Fix,    value=Y,    action=TRAIN);
    %em_property(name=Random_State,    value=1618,    action=TRAIN);
%mend create;