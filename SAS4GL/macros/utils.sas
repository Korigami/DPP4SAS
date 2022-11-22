
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



%macro get_emp_kernel_from_sample(sample);

	%local no_od_samples max_index;
	%let no_of_samples = %size(&sample.);%put &no_of_samples.;
	%let max_index = %numeric_var_count(&sample.);
	


	data tmp(keep=%do i=1 %to %eval(&max_index.*&max_index.); occurence_counter&i. %end; );
		do i=1 to &no_of_samples.;

			set &sample. end=end;
			array cols[%eval(&max_index.+1)] _numeric_; 

			array occurence_counter[%eval(&max_index.*&max_index.)] (%eval(&max_index.*&max_index.) * 0);

			do k=1 to &max_index.;
				if cols[k+1] ne . then
				do;
					do l=1 to &max_index.;
						if cols[l+1] ne . then
						 	occurence_counter[ cols[k+1] + &max_index. * (cols[l+1]-1)] + 1;

					end;
				end;
			end;
		end;
	run;

	data emp_ker_from_&sample.(keep=col1-col&max_index.);
		set tmp;
		array count[*] _numeric_;

		array col[&max_index.] (&max_index.*0);

		do i=1 to &max_index.;
			do j=1 to &max_index.;
		
				if i = j then
					col[j] = count[j + (i-1)*&max_index.]/ &no_of_samples.;
				else
					do;
						tmp = count[j + (j-1)*&max_index.]/ &no_of_samples.
								* count[i + (i-1)*&max_index.]/ &no_of_samples.
								- count[j + (i-1)*&max_index.]/ &no_of_samples.;
						if tmp < 0 then	
							tmp = 0;
						col[j] = sqrt(tmp);
					end;
			end;
			output;
		end;

	run;

	proc delete data = tmp; run;
	
%mend;


%macro plot_matrix(matrix);
	data to_plot(keep=i j x);
		set &matrix.;
		array k[*] _numeric_;
		i=_N_;
		do j=1 to %size(&matrix.);
			x = k[%size(&matrix.)+1-j];
			output;
		end;
	run;

	proc sgrender data=to_plot template=Heatmap; 
	   dynamic _X='i' _Y='j' _Z='x' _T="Basic Heat Map";
	run;

	proc delete data = to_plot;
	run;
%mend;
