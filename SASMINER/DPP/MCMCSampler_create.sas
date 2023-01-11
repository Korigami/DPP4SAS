%macro create;
%em_register(key=Sample, type=data); 

%em_property(name=Method, value=AED);
%em_property(name=Kernel_Type, value=correlation);
%em_property(name=Seed, value=.);

%mend create;