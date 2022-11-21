
%macro trace(K);
	/*œlad macierzy*/
	%local tr trace;
	%let tr=%sysfunc(
			dosubl('data _null_; set &K. end=end; 
					array columns[*] _numeric_; 
					trace + columns[_N_];
					if end then call symputx("trace", round(trace), "G");
					run;'));
	&trace
%mend;

%macro size(set);
		/*liczba wierszy*/
		%local rc k;
        %let rc = %sysfunc(
			dosubl('data _null_; set &set end=end;run; %let k =&sysnobs;')
		);
        &k
%mend;

%macro alter_avail( avail, sel);
/*zamieñ w avail 1 na 0 o ineksie sel*/
	%let avail2 = ;

	%let i = 1;
	%let el = %scan(&avail., &i.);
	%do %while( &el. = 1 or &el. = 0);

		%if &i. = &sel. %then
			%let avail2 = &avail2. 0;
		%else
			%let avail2 = &avail2. %scan(&avail., &i.);

		%let i = %eval( &i.+1);
		%let el = %scan(&avail., &i.);
	%end;

	&avail2.
%mend;


%macro get_avail_indices(avail);
/*poka¿ indeksy avail gdzie s¹ jedynki*/
	%let ind = ;
	%let i = 1;
	%let el = %scan(&avail., &i.);
	%do %while( &el. = 1 or &el. = 0);

		%if %scan(&avail., &i.) = 1 %then
				%let ind = &ind. &i.;

		%let i = %eval( &i.+1);
		%let el = %scan(&avail., &i.);
	%end;
	&ind.
%mend;





