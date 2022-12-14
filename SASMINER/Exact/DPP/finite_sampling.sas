options minoperator mindelimiter=',';
/*%macro init(kernel_type ,K_kernel = ., L_kernel = .,K_e_vals=., K_e_vecs=.,L_e_vals=., L_e_vecs=., dest = sample, mode=GS, random_state=., projection=0, size=.);*/
/*%if ^(&kernel_type. in (correlation, likelihood)) %then %do;*/
/*	%put ERROR: Select kernel from: correlation, likelihood;*/
/*	%goto exit;*/
/*%end;*/
/**/
/*%if ((^%sysfunc(exist(&K_kernel.))) and &K_kernel. ne .) %then %do;*/
/*	%put ERROR: Data set doesnt exist;*/
/*	%goto exit;*/
/*%end;*/

/*check input*/

/*check if K symetric	*/
/*check if projection*/
/*get eig vals and eig vecs*/
/*is_equal_to_0_or_1 for e_vals*/
/*is_in_01 for e_vals)*/
/*is_orthonormal_columns for e_vecs*/

/*the same for L*/

/*something with L Phi and gram_factor*/
/**/
/*%exit:*/
/*%mend();*/

%macro proj_dpp_sampler_kernel(kernel,
								dest ,
								mode=GS ,
								size=.,
								random_state=.);

%if &random_state. = . %then %let random_state = %random_int();

%if &size. ne . %then %do;
	%let rank = %trace( &kernel. );
	%if &size. > &rank. %then %do;
		%put ERROR: size k=&size. > rank=&rank.; 
		%goto exit;

	%end;
%end;

%if &mode.=GS %then %do;
	%proj_dpp_sampler_kernel_GS(&kernel.,&dest., size = &size., random_state=&random_state.);
%end; 
%else %if &mode. = Chol %then %do;
	%proj_dpp_sampler_kernel_Chol(&kernel.,&dest., size = &size., random_state=&random_state.);
%end;

%else %if &mode. = Schur %then %do;
	%proj_dpp_sampler_kernel_Schur(&kernel.,&dest., size = &size., random_state=&random_state.);
%end;

%else %do;
	%let EMEXCEPTIONSTRING = ERROR;
	%put &em_codebar;
	%put ERROR: Invalid sampling mode, choose from: GS, Chol, Schur. Given: &mode.; 
	%put &em_codebar;
	%goto doendm;
%end;
%doendm:
%mend;

%macro proj_dpp_sampler_eig(eig_vecs,
								dest ,
								mode=GS ,
								size=.,
								random_state=.);

%if &random_state. = . %then %let random_state = %random_int();


%if &mode.=GS %then %do;
	%proj_dpp_sampler_eig_GS(&eig_vecs., &dest., size = &size., random_state = &random_state.);
%end; 
%else %if &mode. = GS_bis %then %do;
	%proj_dpp_sampler_eig_GS_bis(&eig_vecs., &dest., size = &size., random_state = &random_state.);
%end;

%else %if &mode. = KuTa12 %then %do;
	%proj_dpp_sampler_eig_KuTa12(&eig_vecs., &dest., size = &size., random_state = &random_state.);
%end;

%else %do;
	%let EMEXCEPTIONSTRING = ERROR;
	%put &em_codebar;
	%put ERROR: Invalid sampling mode, choose from: GS, GS_bis, KuTa12. Given: &mode.; 
	%put &em_codebar;
	%goto doendm;
%end;
%doendm:
%mend;

%macro K_eig_vecs(eig_vals,eig_vecs);
proc iml;
	use &eig_vals.; read all into e_vals; close;
	use &eig_vecs.; read all into e_vecs; close;
	e_vals = e_vals`;
	V = {};
	do i=1 to ncol(e_vals);
		if e_vals[i] > 0.5 then V = V || i;
	end;
	V = e_vecs[1:nrow(e_vecs),V];
	create proj_dpp_V from V; append from V;
quit;
%mend;

%macro sample_exact(kernel_type ,K_kernel = ., L_kernel = .,K_e_vals=., K_e_vecs=.,L_e_vals=., L_e_vecs=., dest = sample, mode=GS, random_state=., projection=0, size=.);
/*%init(&kernel_type. ,K_kernel = &K_kernel., L_kernel = &L_kernel.,K_e_vals=&K_e_vals., K_e_vecs=&K_e_vecs.,L_e_vals=&L_e_vals., L_e_vecs=&L_e_vecs., dest = &dest., mode=&mode., random_state=&random_state., projection=&projection., size=&size.);*/

%if &random_state. = . %then %let random_state = %random_int();
%if &mode. = Schur %then %do;
	%put NOTE: Schur;
	%if &kernel_type. = correlation & &projection. = Y %then %do;
		%computeK(K_kernel=K_kernel., L_kernel = &L_kernel. ,L_e_vals=&L_e_vals., L_e_vecs=&L_e_vecs.,K_e_vals=&K_e_vals., K_e_vecs=&K_e_vecs.);
		%if &K_kernel. = . %then %let K_kernel = K;
		%proj_dpp_sampler_kernel(&K_kernel., &dest.,mode=&mode., random_state = &random_state., size = &size.);
	%end;
	%else %do;
		%let EMEXCEPTIONSTRING = ERROR;
		%put &em_codebar;
		%put ERROR: Schur sampling mode is only available for kernel_type=correlation and projection=1. Given kernel_type=&kernel_type. and projection=&projection.;
		%put &em_codebar;
		%goto doendm;
	%end;

%end;

%else %if &mode. = Chol %then %do;
	
	%computeK(K_kernel=&K_kernel., L_kernel = &L_kernel. ,L_e_vals=&L_e_vals., L_e_vecs=&L_e_vecs.,K_e_vals=&K_e_vals., K_e_vecs=&K_e_vecs.);
	%if &K_kernel. = . %then %let K_kernel = K;
	%if &kernel_type. = correlation & &projection. = Y %then %do;
		%put NOTE: Chol proj;
		%proj_dpp_sampler_kernel(&K_kernel., &dest.,mode=&mode., random_state = &random_state.,  size = &size.);
	%end;
	%else %do;
		%put NOTE: Chol;
		%dpp_sampler_generic_kernel(&K_kernel., &dest., random_state = &random_state.,  size = .);
	%end;
%end;



%else %if &K_e_vals. ne . %then %do;
	%put NOTE: K_e_vals ne ., eig sampler;
	%if &kernel_type. = correlation & &projection. = Y %then %do;
		%K_eig_vecs(&K_e_vals.,&K_e_vecs.);
	%end;
	%else %do;
		%dpp_eig_vecs_selector(&K_e_vals.,&K_e_vecs., random_state = &random_state.);
	%end;
	%proj_dpp_sampler_eig(proj_dpp_V, &dest., mode = &mode.,random_state = &random_state.,  size = &size.);
/*	proc datasets nolist lib=work;*/
/*		delete proj_dpp_V;*/
/*	run;*/
	
%end;

%else %if &L_e_vals. ne . %then %do;
	%put NOTE: L_e_vals ne ., eig sampler;
	%compute_K_evals_from_L_evals(&L_e_vals.);
	%let K_e_vals = K_e_vals;
	%let K_e_vecs = &L_e_vecs.;
	%if &kernel_type. = correlation & &projection. = Y %then %do;
		%K_eig_vecs(&K_e_vals.,&K_e_vecs.);
	%end;
	%else %do;
		%dpp_eig_vecs_selector(&K_e_vals.,&K_e_vecs., random_state = &random_state.);
	%end;
	%proj_dpp_sampler_eig(proj_dpp_V, &dest., mode = &mode.,random_state = &random_state.,  size = &size.);
	proc datasets nolist lib=work;
		delete proj_dpp_V K_e_vals;
	run;
%end;

%else %if (&K_kernel. ne .) and (&projection. = Y) %then %do;
	%put NOTE: K_kernel ne . and  p=1;
	%proj_dpp_sampler_kernel(&K_kernel., &dest.,mode=&mode., random_state = &random_state.,  size = &size.);
%end;
%else %if &K_kernel. ne .  %then %do;
	%put NOTE: K_kernel ne . ;
	%get_eigendecomposition(&K_kernel., eig_vals, eig_vecs);
	%let  K_e_vals = eig_vals;
	%let  K_e_vecs = eig_vecs;
	%dpp_eig_vecs_selector(&K_e_vals.,&K_e_vecs., random_state = &random_state.);
	%proj_dpp_sampler_eig(proj_dpp_V, &dest., mode = &mode.,random_state = &random_state.,  size = &size.);
	proc datasets nolist lib=work;
		delete proj_dpp_V eig_vals eig_vecs;
	run;
%end;
%else %if &L_kernel. ne . %then %do;
	%put NOTE: L_kernel ne .;
	%get_eigendecomposition(&L_kernel., eig_vals, eig_vecs);
	%let  L_e_vals = eig_vals;
	%let  L_e_vecs = eig_vecs;

	%compute_K_evals_from_L_evals(&L_e_vals.);
	%let K_e_vals = K_e_vals;
	%let K_e_vecs = &L_e_vecs.;
	%if &kernel_type. = correlation & &projection.=Y %then %do;
		%K_eig_vecs(&K_e_vals.,&K_e_vecs.);
	%end;
	%else %do;
		%dpp_eig_vecs_selector(&K_e_vals.,&K_e_vecs., random_state = &random_state.);
	%end;
	%proj_dpp_sampler_eig(proj_dpp_V, &dest., mode = &mode.,random_state = &random_state.,  size = &size.);
	proc datasets nolist lib=work;
		delete proj_dpp_V K_e_vals eig_vals eig_vecs;
	run;
%end;

%else %do;
	%let EMEXCEPTIONSTRING = ERROR;
	%put &em_codebar;
	%put ERROR: None of the available samplers could be used based on the current DPP representation.;
	%put &em_codebar;
	%goto doendm;
%end;




proc datasets nolist lib=work;
	delete computed_K K;
run;

%doendm:
%mend;

%macro compute_K_from_eig(K_e_vals, K_e_vecs);
	proc iml;
		use &K_e_vals.; 	read all  into K_e_vals; 		close;
		use &K_e_vecs.; 	read all  into K_e_vecs; 		close;
		if dimension(K_e_vals)[1] ^= 1 then K_e_vals = K_e_vals`;
		K = (K_e_vecs # K_e_vals)*K_e_vecs`;
		create K from K;
				append from K;
		close K;
	quit;
%mend;
%macro compute_K_evals_from_L_evals(L_e_vals);
	proc iml;
		use &L_e_vals.; 	read all  into L_e_vals; 		close;
		if dimension(L_e_vals)[1] ^= 1 then L_e_vals = L_e_vals`;
		K_e_vals = L_e_vals / (1 + L_e_vals);
		create K_e_vals from K_e_vals;
				append from K_e_vals;
		close K_e_vals;
	quit;
%mend;

%macro computeK(K_kernel = ., L_kernel = .,K_e_vals=., K_e_vecs=.,L_e_vals=., L_e_vecs=.);
	%if &K_kernel. ne . %then %do;
		%put &=K_kernel;
		%put NOTE: K already computed;
		%goto exit;
	%end;
	%else %if (&K_e_vals. ne .) & (&K_e_vecs. ne .) %then %do;
		%compute_K_from_eig(&K_e_vals., &K_e_vecs.);
		%let K_kernel = K;
	%end;
	%else %if &L_e_vals. ne . %then %do;
		%put liczy K z L;
		%compute_K_evals_from_L_evals(&L_e_vals.);
		%compute_K_from_eig(K_e_vals, &L_e_vecs.);
		proc datasets nolist lib=work;
		delete K_e_vals;
		run;
		%let K_kernel = K;
	%end;
	%else %if &L_kernel. ne . %then %do;
		%compute_K(&L_kernel., K);
		%let K_kernel = K;
	%end;
	%exit:
%mend;

%macro dpp_eig_vecs_selector(eig_vals, eig_vecs,random_state=.);
	proc iml;
		use &eig_vals.; 	read all  into eig_vals; 		close;
		CALL STREAMINIT(&random_state.);
		if dimension(eig_vals)[1] ^= 1 then eig_vals = eig_vals`;
		y = {};
		do i = 1 to dimension(eig_vals)[2];
			x = rand("uniform");
			if x < eig_vals[i] then do;
				y = y || i;
/*				print(x);*/
/*				print(eig_vals[i]);*/
/*				print(y);*/
			end;
		end;
		
		use &eig_vecs.; 	read all  into eig_vecs; 	close;
		eig_vecs = eig_vecs[,y];
		create proj_dpp_V from eig_vecs; append from eig_vecs;
	quit;
%mend;

/*%sample_exact(correlation, K_e_vals = SASHELP.eig_vals, K_e_vecs = SASHELP.eig_vecs,dest=eig_GS_dpp, projection=0);*/

/*%sample_exact(correlation, K_kernel=SASHELP.K ,dest=GS, projection=1);*/

/*%matrix_K(10, 8);*/

/*Working methods - K kernel and K eig*/
/*DPP*/
/*%sample_exact(correlation, L_kernel=K ,dest=GS, projection=1);*/
/*%sample_exact(correlation, K_kernel=K ,mode = Schur ,dest=Schur, projection=1);*/
/*%sample_exact(correlation, K_kernel=K ,mode = Chol ,dest=Chol, projection=1);*/
/*%sample_exact(correlation, K_kernel=K ,dest=generic, mode=Chol, projection=0);/*generic sampler*/*/
/*%sample_exact(likelyhood, K_kernel= K ,dest=generic2, mode=Chol, projection=1);/*generic sampler*/*/
/**/
/*%get_eigendecomposition(K, eig_vals, eig_vecs);*/
/*%sample_exact(correlation, K_e_vals = eig_vals, K_e_vecs = eig_vecs,dest=eig_GS, projection=1);*/
/*%sample_exact(correlation, K_e_vals = eig_vals, K_e_vecs = eig_vecs,dest=KuTa12, mode = KuTa12, projection=1);*/

/*K-DPP*/
/*%sample_exact(correlation, K_kernel=K ,dest=GS_kdpp, projection=1, size = 4);*/
/*%sample_exact(correlation, K_kernel=K ,mode = Schur ,dest=Schur_kdpp, projection=1 , size = 4);*/
/*%sample_exact(correlation, K_kernel=K ,mode = Chol ,dest=Chol_kdpp, projection=1 , size = 4);*/
/**/
/*%get_eigendecomposition(K, eig_vals, eig_vecs);*/
/*%sample_exact(correlation, K_e_vals = eig_vals, K_e_vecs = eig_vecs,dest=eig_GS_dpp, projection=1 , size = 4);*/
/*%sample_exact(correlation, K_e_vals = eig_vals, K_e_vecs = eig_vecs,dest=KuTa12_dpp, mode = KuTa12, projection=1 , size = 4);*/



/*Working methods - L kernel */
/*DPP*/
/**/
/*%sample_exact(correlation, L_kernel=K,  dest = test1, mode=Schur, projection=1);*/
/*%sample_exact(correlation, L_kernel=K ,dest = test2, mode=Chol, projection=1);*/
/*%sample_exact(correlation, L_kernel=K,  dest = test3, projection=1);*/


/*Error check*/


/*%sample_exact(correlation, K_kernel=K ,dest=test1, mode=Schur, projection=0);/*error check*/*/
/*%sample_exact(likelyhood, K_kernel= K ,dest=test1, mode=Schur, projection=1);/*error check*/*/

/*%sample_exact(correlation, K_kernel= K,dest=eig_GS, projection=0);*/
/*%sample_exact(correlation, K_kernel= K ,dest=KuTa12, mode = KuTa12, projection=1);*/
/**/
/**/
/**/
/*%sample_exact(correlation, K_e_vals = eig_vals, K_e_vecs = eig_vecs, projection=0);*/
/*%sample_exact(correlation, K_e_vals = eig_vals, K_e_vecs = eig_vecs, mode = KuTa12, projection=0);*/

*/*/;








