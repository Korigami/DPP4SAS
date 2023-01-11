%macro CreateSources(dirPath, libraryname, catname, modulename);

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

	%if &catname. = BTFLL | &catname. = BTBND %then
	%do;
	    filename src_cat catalog "&libraryname..&catname..random_matrices.source";
		filename mydata "&dirPath.\random_matrices.sas";
		data _null_;
			file src_cat;
			infile mydata;
			input;
			put _infile_;
		run;
	%end;

	%if &catname. = PSSPL %then
	%do;
		filename src_cat catalog "&libraryname..&catname..poisson_final.source";
		filename mydata "&dirPath.\poisson_final.sas";
		data _null_;
			file src_cat;
			infile mydata;
			input;
			put _infile_;
		run;
	%end;

	%if &catname. = STTOD %then
	%do;
		filename src_cat catalog "&libraryname..&catname..stationary_1_d.source";
		filename mydata "&dirPath.\stationary_1_d.sas";
		data _null_;
			file src_cat;
			infile mydata;
			input;
			put _infile_;
		run;
	%end;

	proc catalog cat=&libraryname..&catname.;
		copy out=sashelp.&libraryname.;
	run;

%mend CreateSources;

%CreateSources(C:\Users\sas\Desktop\DPPMINER-20221213T231047Z-001\DPPMINER\BetaFull\DPP, DPP, BTFLL, BetaFullSampler);
LIBNAME DPP CLEAR;
%CreateSources(C:\Users\sas\Desktop\DPPMINER-20221213T231047Z-001\DPPMINER\BetaFull\DPP, DPP, BTBND, BetaBandSampler);
LIBNAME DPP CLEAR;
%CreateSources(C:\Users\sas\Desktop\asd, DPP, PSSPL, PoissPlanchSampler);
LIBNAME DPP CLEAR;
%CreateSources(C:\Users\sas\Desktop\as\minerowe_makra-20230110T165658Z-001\minerowe_makra, DPP, STTOD, StatOneDSampler);
LIBNAME DPP CLEAR;
