
%macro create;
    %em_register(key=Sample, type=data);  
    %em_property(name=Theta,    value=1,    action=TRAIN);
    %em_property(name=Random_State,    value=1618,    action=TRAIN);
%mend create;
