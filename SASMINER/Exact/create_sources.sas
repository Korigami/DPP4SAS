%macro CreateSources(dirPath, libraryname, catname);

	libname &libraryname "&dirPath";
	filename src_cat catalog "&libraryname..&catname..exact.source";
	filename mydata "&dirPath.\main.sas";
	data _null_;
		file src_cat;
		infile mydata;
		input;
		put _infile_;
	run;
	filename src_cat catalog "&libraryname..&catname..exact_create.source";
	filename mydata "&dirPath.\create.sas";
	data _null_;
		file src_cat;
		infile mydata;
		input;
		put _infile_;
	run;
	filename src_cat catalog "&libraryname..&catname..exact_train.source";
	filename mydata9 "&dirPath.\proj_chol_sampl.sas";
	filename mydata2 "&dirPath.\proj_dpp_GS_cor.sas";
	filename mydata3 "&dirPath.\utils.sas";
	filename mydata4 "&dirPath.\proj_dpp_sampler_kernel_eig_GS_corre.sas";
	filename mydata5 "&dirPath.\proj_dpp_sampler_kernel_eig_KuTa12.sas";
	filename mydata6 "&dirPath.\proj_dpp_sampler_kernel_schur.sas";
	filename mydata7 "&dirPath.\sampler_generic.sas";
	filename mydata8 "&dirPath.\finite_sampling.sas";
	filename mydata1 "&dirPath.\train.sas";
	options mprint;
	%do i=1 %to 9;
	
	data _null_;
		%if &i = 1 %then %do;
		file src_cat;
		%end;
		%else %do;
		file src_cat mod;
		%end;
		infile mydata&i.;
		input;
		put _infile_;
	run;
	%end;

	filename src_cat catalog "&libraryname..&catname..exact_score.source";
	filename mydata "&dirPath.\score.sas";
	data _null_;
		file src_cat;
		infile mydata;
		input;
		put _infile_;
	run;


	filename src_cat catalog "&libraryname..&catname..exact_report.source";
	filename mydata "&dirPath.\report.sas";
	data _null_;
		file src_cat;
		infile mydata;
		input;
		put _infile_;
	run;

	proc catalog cat=&libraryname..&catname.;
		copy out=sashelp.&catname.;
	run;
%mend CreateSources;

%macro createCatalog(userPath);
	%CreateSources(&userPath, exact, exact)
	LIBNAME exact CLEAR;
%mend createCatalog;




/*proc catalog cat=sashelp.exact;*/
/*   delete exact_train.source;*/
/*run;*/

/*%get_eigendecomposition(SASHELP.K, SASHELP.eig_vals, SASHELP.eig_vecs);*/
/*%matrix_K(10, 10);*/
/*%compute_L(K, SASHELP.L);*/
/*%get_eigendecomposition(SASHELP.L, SASHELP.L_eig_vals, SASHELP.L_eig_vecs);*/




%createCatalog(C:\Users\kasia\OneDrive\Dokumenty\GitHub\DPP4SAS\SAS4GL\macros\finite_sampling);



