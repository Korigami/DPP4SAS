/*%macro init();*/
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
	%put ERROR: Invalid sampling mode, choose from: GS, Chol, Schur. Given: &mode.; 
	%goto exit;
%end;
%exit:
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
	%put ERROR: Invalid sampling mode, choose from: GS, GS_bis, KuTa12. Given: &mode.; 
	%goto exit;
%end;
%exit:
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
/*%init()*/

%if &random_state = . %then %let random_state = %random_int();

%if &mode. = Schur %then %do;
	%if &kernel_type. = correlation & &projection. %then %do;
		%let K = &K_kernel.;
		%if &K_kernel. = . & &L_kernel. ne . %then %do;
			%compute_K(&L_kernel., computed_K);
			%let K = computed_K; 
		%end;
		%else %if &K_kernel. = . %then %do;
			%put ERROR: Only compute from L To K;
			%goto exit;
		%end;
		%proj_dpp_sampler_kernel(&K., &dest.,mode=&mode., random_state = &random_state., size = &size.);

	%end;
	%else %do;
		%put ERROR: Schur sampling mode is only available for kernel_type=correlation and projection=1. Given kernel_type=&kernel_type. and projection=&projection.;
		%goto exit;
	%end;

%end;

%else %if &mode. = Chol %then %do;
	%let K = &K_kernel.;
	%if &K_kernel. = . & &L_kernel. ne . %then %do;
		%compute_K(&L_kernel., computed_K);
		%let K = computed_K; 
	%end;
	%else %if &K_kernel. = . %then %do;
		%put ERROR: Only compute from L To K;
		%goto exit;
	%end;
	%if &kernel_type. = correlation & &projection. %then %do;
		%proj_dpp_sampler_kernel(&K., &dest.,mode=&mode., random_state = &random_state.,  size = &size.);
	%end;
	%else %do;
		%dpp_sampler_generic_kernel(&K., &dest., random_state = &random_state.,  size = &size.);
	%end;
%end;



%else %if &K_e_vals. ne . %then %do;
	%if &kernel_type. = correlation & &projection. %then %do;
		%K_eig_vecs(&K_e_vals.,&K_e_vecs.);
		%proj_dpp_sampler_eig(proj_dpp_V, &dest., mode = &mode.,random_state = &random_state.,  size = &size.);
		proc datasets nolist lib=work;
			delete proj_dpp_V;
		run;
	%end;
	%else %do;
		%put ERROR: Sampling stoped. Need for eig generator;
		%goto exit;
	%end;
	
%end;
%else %if &K_kernel. ne . & &projection. = 1 %then %proj_dpp_sampler_kernel(&K_kernel., &dest.,mode=&mode., random_state = &random_state.,  size = &size.);
%else %if &K_kernel. = . & &projection. = 1 %then %do;
	%let K = &K_kernel.;
	%if &K_kernel. = . & &L_kernel. ne . %then %do;
		%compute_K(&L_kernel., computed_K);
		%let K = computed_K; 
	%end;
	%else %if &K_kernel. = . %then %do;
		%put ERROR: Only compute from L To K;
		%goto exit;
	%end;
	%proj_dpp_sampler_kernel(&K., &dest.,mode=&mode., random_state = &random_state.,  size = &size.);
%end;
%else %if &K_kernel. ne . %then %do;
	%put ERROR: If you want to generate without projection pass K eigen vector and K eigen values or select different mode;
%end;




proc datasets nolist lib=work;
	delete computed_K;
run;

%exit:
%mend;


%matrix_K(10, 8);

/*Working methods - K kernel and K eig*/
/*DPP*/
%sample_exact(correlation, K_kernel=K ,dest=GS, projection=1);
%sample_exact(correlation, K_kernel=K ,mode = Schur ,dest=Schur, projection=1);
%sample_exact(correlation, K_kernel=K ,mode = Chol ,dest=Chol, projection=1);
%sample_exact(correlation, K_kernel=K ,dest=generic, mode=Chol, projection=0);/*generic sampler*/
%sample_exact(likelyhood, K_kernel= K ,dest=generic2, mode=Chol, projection=1);/*generic sampler*/

%get_eigendecomposition(K, eig_vals, eig_vecs);
%sample_exact(correlation, K_e_vals = eig_vals, K_e_vecs = eig_vecs,dest=eig_GS, projection=1);
%sample_exact(correlation, K_e_vals = eig_vals, K_e_vecs = eig_vecs,dest=KuTa12, mode = KuTa12, projection=1);

/*K-DPP*/
%sample_exact(correlation, K_kernel=K ,dest=GS_kdpp, projection=1, size = 4);
%sample_exact(correlation, K_kernel=K ,mode = Schur ,dest=Schur_kdpp, projection=1 , size = 4);
%sample_exact(correlation, K_kernel=K ,mode = Chol ,dest=Chol_kdpp, projection=1 , size = 4);

%get_eigendecomposition(K, eig_vals, eig_vecs);
%sample_exact(correlation, K_e_vals = eig_vals, K_e_vecs = eig_vecs,dest=eig_GS_dpp, projection=1 , size = 4);
%sample_exact(correlation, K_e_vals = eig_vals, K_e_vecs = eig_vecs,dest=KuTa12_dpp, mode = KuTa12, projection=1 , size = 4);



/*Working methods - L kernel */
/*DPP*/

%sample_exact(correlation, L_kernel=K,  dest = test1, mode=Schur, projection=1);
%sample_exact(correlation, L_kernel=K ,dest = test2, mode=Chol, projection=1);
%sample_exact(correlation, L_kernel=K,  dest = test3, projection=1);


/*Error check*/


%sample_exact(correlation, K_kernel=K ,dest=test1, mode=Schur, projection=0);/*error check*/
%sample_exact(likelyhood, K_kernel= K ,dest=test1, mode=Schur, projection=1);/*error check*/

%sample_exact(correlation, K_kernel= K,dest=eig_GS, projection=0);
%sample_exact(correlation, K_kernel= K ,dest=KuTa12, mode = KuTa12, projection=1);



%sample_exact(correlation, K_e_vals = eig_vals, K_e_vecs = eig_vecs, projection=0);
%sample_exact(correlation, K_e_vals = eig_vals, K_e_vecs = eig_vecs, mode = KuTa12, projection=0);










