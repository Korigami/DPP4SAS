%macro report();

proc print data=&EM_EXPORT_TRAIN.;
run;

%mend;