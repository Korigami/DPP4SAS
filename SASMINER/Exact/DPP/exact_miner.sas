%macro main;
%if %upcase(&EM_ACTION) = CREATE %then %do;
/*add CREATE code */
%else;
%if %upcase(&EM_ACTION) = TRAIN %then %do;
/*add TRAIN code */
%else;
%if %upcase(&EM_ACTION) = SCORE %then %do;
/*add SCORE code */
%else;
%if %upcase(&EM_ACTION) = REPORT %then %do;
/*add REPORT code */
%mend main;

%main;




%macro createSource(dirPath, libraryname, catname);
	libname &libraryname "&dirPath";
	filename source_cat catalog "&libraryname..&catname..dpp.source";
%mend;