%macro main;
	%if %upcase(&EM_ACTION) = CREATE %then %do;
		filename temp catalog 'sashelp.DPP.MCMCSampler_create.source';
		%include temp;
		filename temp;
		%create;
	%end;

	%else
	%if %upcase(&EM_ACTION) = TRAIN %then %do;
		filename temp catalog 'sashelp.DPP.MCMCSampler_train.source';
		%include temp;
		filename temp;
		%train;
	%end;

	%if %upcase(&EM_ACTION) = SCORE %then %do;
		filename temp catalog 'sashelp.DPP.MCMCSampler_score.source';
		%include temp;
		filename temp;
		%score;
	%end;
	
	%if %upcase(&EM_ACTION) = REPORT %then %do;
		filename temp catalog 'sashelp.DPP.MCMCSampler_report.source';
		%include temp;
		filename temp;
		%report;
	%end;

%mend main;
%main;