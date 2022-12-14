%macro create;
/* Training Properties */
%em_property(name=Method, value=GS);
%em_property(name=Projection, value=N);
%em_property(name=Size, value=.);
%em_property(name=Kernel_Type, value=correlation);
%em_property(name=KDPP, value=N);
%em_property(name=Seed, value=.);

/* Scoring Properties */
/* Register Data Sets */
%EM_REGISTER(key=OUTEST, type=DATA);
%EM_REGISTER(key=EFFECTS, type=DATA);
%mend create;