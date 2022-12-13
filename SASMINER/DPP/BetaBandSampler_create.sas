%macro create;
    %em_register(key=Sample, type=data);  
    
    %em_property(name=EnsembleVersion,    value=HERMITE,    action=TRAIN);
    %em_property(name=Size,    value=10,    action=TRAIN);
    %em_property(name=Beta,    value=2,    action=TRAIN);

    %em_property(name=Loc,    value=0.0,    action=TRAIN);
    %em_property(name=Scale,    value=1.0,    action=TRAIN);
    %em_property(name=Shape,    value=1.0,    action=TRAIN);
    %em_property(name=A,    value=1.0,    action=TRAIN);
    %em_property(name=B,    value=1.0,    action=TRAIN);

    %em_property(name=Normalize,    value=Y,    action=TRAIN);
    %em_property(name=Heurestic_Fix,    value=Y,    action=TRAIN);
    %em_property(name=Random_State,    value=1618,    action=TRAIN);
%mend create;