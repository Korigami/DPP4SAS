filename temp3 catalog 'sashelp.dpp.utils.source';                                                                             
%include temp3;                                                                                                                  
filename temp3;


filename temp2 catalog 'sashelp.dpp.MCMCSampler_source.source';                                                                             
%include temp2;                                                                                                                  
filename temp2;
%macro train;


%em_getname(key=Sample, type=data); 


%if %sysfunc(index(&EM_DEBUG, SOURCE))>0 or
%sysfunc(index(&EM_DEBUG, ALL))>0 %then %do;
options mprint;
%end;
/**/
/**/
/**/
/*%if (^%sysfunc(exist(&EM_IMPORT_DATA)) and*/
/*	^%sysfunc(exist(&EM_IMPORT_DATA, VIEW)))*/
/*	or "&EM_IMPORT_DATA" eq "" %then %do;*/
/*		%let EMEXCEPTIONSTRING = exception.server.IMPORT.NOTRAIN,1;*/
/*		%goto doenda;*/
/*%end;*/
/**/
/**/
%if (^%sysfunc(exist(&EM_IMPORT_VALIDATE)) and
	^%sysfunc(exist(&EM_IMPORT_VALIDATE, VIEW)))
	or "&EM_IMPORT_VALIDATE" eq "" %then %do;
		%let s_init = .;
%end;
%else %do;
	%let s_init = &EM_IMPORT_VALIDATE.;
%end;

%let x = %size(&em_import_data.);
%put &=x;

/*data _null_;*/
/*set &em_import_data. end=end;*/
/*put _all_;*/
/*run;*/
/*%let k =&sysnobs;*/
/*%put &=k;*/


/*%sample_mcmc(K=&em_import_data., dest =&EM_USER_SAMPLE., sampling_mode = AD ,kernel_type = likelihood ,projection = Y, random_state=547);*/

%sample_mcmc(K=&em_import_data., dest = &EM_USER_SAMPLE., sampling_mode = &EM_PROPERTY_Method. ,
kernel_type = &EM_PROPERTY_KERNEL_TYPE. ,projection = Y, random_state=&EM_PROPERTY_SEED.,s_init = &s_init);

data &EM_EXPORT_TRAIN;                                                                                                               
    set &em_user_Sample.;                                                                                                          
run; 

%mend train;
