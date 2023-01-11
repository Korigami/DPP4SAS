%macro create;
    %em_register(key=Sample, type=data);  
    %em_property(name=Mode,         value=descent,    action=TRAIN);
    %em_property(name=Size,         value=100,          action=TRAIN);
    %em_property(name=Base,         value=1,            action=TRAIN);
    %em_property(name=X0,           value=0.5,          action=TRAIN);
    %em_property(name=Random_State, value=1618,         action=TRAIN);
%mend create;