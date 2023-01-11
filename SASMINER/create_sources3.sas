%macro CreateSources(dirPath);
	
	%let list = ExactSampler MCMCSampler BetaFullSampler BetaBandSampler PoissPlanchSampler StatOneDSampler;
	%let list2 = Exact MCMC BTFLL BTBND PSSPL STTOD;
	%let libraryname = DPP;
	%do i = 1 %to 2;
		%let catname = %scan(&list2.,&i.);
		%let modulename = %scan(&list.,&i.);

		libname &libraryname "&dirPath";
		filename src_cat catalog "&libraryname..&catname..&modulename..source";
		filename mydata "&dirPath.\&modulename..sas";
		data _null_;
			file src_cat;
			infile mydata;
			input;
			put _infile_;
		run;
		filename src_cat catalog "&libraryname..&catname..&modulename._create.source";
		filename mydata "&dirPath.\&modulename._create.sas";
		data _null_;
			file src_cat;
			infile mydata;
			input;
			put _infile_;
		run;
		filename src_cat catalog "&libraryname..&catname..&modulename._train.source";
		filename mydata "&dirPath.\&modulename._train.sas";
		data _null_;
			file src_cat;
			infile mydata;
			input;
			put _infile_;
		run;
		filename src_cat catalog "&libraryname..&catname..&modulename._score.source";
		filename mydata "&dirPath.\&modulename._score.sas";
		data _null_;
			file src_cat;
			infile mydata;
			input;
			put _infile_;
		run;
		filename src_cat catalog "&libraryname..&catname..&modulename._report.source";
		filename mydata "&dirPath.\&modulename._report.sas";
		data _null_;
			file src_cat;
			infile mydata;
			input;
			put _infile_;
		run;

	%end;
	filename src_cat catalog "&libraryname..&catname..ExactSampler_source.source";
	filename mydata "&dirPath.\ExactSampler_source.sas";
	data _null_;
		file src_cat;
		infile mydata;
		input;
		put _infile_;
	run;

	filename src_cat catalog "&libraryname..&catname..MCMCSampler_source.source";
	filename mydata "&dirPath.\MCMCSampler_source.sas";
	data _null_;
		file src_cat;
		infile mydata;
		input;
		put _infile_;
	run;

	filename src_cat catalog "&libraryname..&catname..utils.source";
	filename mydata "&dirPath.\utils.sas";
	data _null_;
		file src_cat;
		infile mydata;
		input;
		put _infile_;
	run;

	filename src_cat catalog "&libraryname..&catname..random_matrices.source";
	filename mydata "&dirPath.\random_matrices.sas";
	data _null_;
		file src_cat;
		infile mydata;
		input;
		put _infile_;
	run;

	filename src_cat catalog "&libraryname..&catname..poisson_final.source";
	filename mydata "&dirPath.\poisson_final.sas";
	data _null_;
		file src_cat;
		infile mydata;
		input;
		put _infile_;
	run;

	filename src_cat catalog "&libraryname..&catname..poisson_final.source";
	filename mydata "&dirPath.\stationary_1_d.sas";
	data _null_;
		file src_cat;
		infile mydata;
		input;
		put _infile_;
	run;


	proc catalog cat=&libraryname..&catname.;
		copy out=sashelp.&libraryname.;
	run;


%mend CreateSources;





%CreateSources(path to ... \DPP4SAS\SASMINER\DPP);
LIBNAME DPP CLEAR;


