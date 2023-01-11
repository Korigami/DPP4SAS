filename temp3 catalog 'sashelp.dpp.utils.source';                                                                             
%include temp3;                                                                                                                  
filename temp3;


filename temp2 catalog 'sashelp.dpp.ExactSampler_source.source';                                                                             
%include temp2;                                                                                                                  
filename temp2;

%macro train;
                                                                                                                                                                                                                       
%if (^%sysfunc(exist(&EM_IMPORT_VALIDATE)) and
	^%sysfunc(exist(&EM_IMPORT_VALIDATE, VIEW)) and ^%sysfunc(exist(&EM_IMPORT_TEST)) and
	^%sysfunc(exist(&EM_IMPORT_TEST, VIEW)))
	or (("&EM_IMPORT_VALIDATE" eq "") and ("&EM_IMPORT_TEST" eq ""))%then %do;
	%let e_vec = .;
	%let e_val = .;
	%if (^%sysfunc(exist(&EM_IMPORT_DATA)) and
	^%sysfunc(exist(&EM_IMPORT_DATA, VIEW)))
	or "&EM_IMPORT_DATA" eq "" %then %do;
		%let EMEXCEPTIONSTRING = exception.server.IMPORT.NOTRAIN,1;
		%goto doenda;
	%end;
%end;
%else %do;
	%let e_vec = &EM_IMPORT_VALIDATE.;
	%let e_val = &EM_IMPORT_TEST.;
%end;



%if &EM_PROPERTY_KDPP. = Y %then %do;
	%let s = &em_property_Size.;
%end;
%else %do;
	%let s=.;
%end;

/*%put &=EM_PROPERTY_KERNEL_TYPE;*/
/*%put &=em_import_data;*/
/*%put &=EM_PROPERTY_Method;*/
/*%put &=em_export_train;*/
/*%put &=p;*/
/*%put &=EM_PROPERTY_SEED;*/
/*%put &=s;*/
/*%put &=e_val;*/
/*%put &=e_vec;*/

%if &EM_PROPERTY_KERNEL_TYPE. = correlation %then %do;
%sample_exact(&EM_PROPERTY_KERNEL_TYPE., K_kernel=&em_import_data. ,
					mode = &EM_PROPERTY_Method. ,dest=&em_export_train. ,
					projection=&EM_PROPERTY_Projection., random_state=&EM_PROPERTY_SEED., size=&s.,
					K_e_vals = &e_val., K_e_vecs = &e_vec.);
%end;
%else %do;
%sample_exact(&EM_PROPERTY_KERNEL_TYPE., L_kernel=&em_import_data. ,
					mode = &EM_PROPERTY_Method. ,dest=&em_export_train. ,
					projection=&EM_PROPERTY_Projection., random_state=&EM_PROPERTY_SEED., size=&s.,
					L_e_vals = &e_val., L_e_vecs = &e_vec.);
%end;



%doenda:

/*proc print data=sashelp.eig_vals;*/
/*run;*/
%mend train;