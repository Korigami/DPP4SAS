
%macro select_index_Chol(rank, it);
	%local sel;
	%let sell = %sysfunc(dosubl('
							data _null_; 
								a=&it.;
								set R point=a;
								cum_sum = 0;
								k=&it-1;
								d_point = &it;
								do while(x > cum_sum);
									set d point=d_point;
									cum_sum + abs( col1  / (&rank. - &it. + 1) );
									k + 1;
									d_point + 1;
								end;
								call symputx("sel", k, "L");
								stop;
							run;
						')
				);
		&sel.
%mend;

%macro adjust_probabilities_chol(sel, it, size);
	%let tmp = %sysfunc(dosubl('
proc iml;
	use A; 	read all  into A; 	close;
	if &sel.+1 <= &size. then do;
		A[&sel.+1:&size.,{&it.,&sel.}] = A[&sel.+1:&size.,{&sel., &it.}];
	end;
	A[{&it., &sel.},{&it., &sel.}] = A[{&sel., &it.},{&sel., &it.}];
	if &it.> 1 then do;
		A[{&it., &sel.}, 1:&it.-1] = A[{&sel., &it.}, 1:&it.-1];

	end;
	use orig_indices; 	read all  into orig_indices; 	close;
	orig_indices[{&it., &sel.},] = orig_indices[{&sel.,&it.},];
	use d; 	read all  into d; 	close;
	d[{&it., &sel.}] = d[{&sel., &it.}];
	A[&it., &it.] = sqrt(d[&it.]);

	if &it > 1 then do;
		A[&it. + 1:&size., &it.] = A[&it. + 1:&size., &it.] - A[&it. + 1:&size., 1:&It.-1]*A[&it., 1:&it.-1]`;
	end;
	A[&it. + 1:&size., &it.] = A[&it. + 1:&size., &it.]/A[&it., &it.];
	d[&it + 1:&size.] = d[&it + 1:&size.] - A[&it. + 1:&size., &it.]##2;

	create A from A;
					append from A;
	close A;
	create orig_indices from orig_indices;
					append from orig_indices;
	close orig_indices;
	create d from d;
					append from d;
	close d;

quit;
'));
%mend;


%macro proj_dpp_sampler_kernel_Chol(K_set, dest, size=., random_state=.);
	%let N = %size( &K_set. );
	%let rank = %trace( &K_set. );
	%put &=rank;
	%if &size = . %then %let size = &rank.;
	%if &random_state = . %then %let random_state = %random_int();
	%let tmp = %sysfunc(dosubl('
		data A;
			set &K_set.;
		run;
	'));

	%let tmp = %sysfunc( dosubl('
		data d( keep=COL1 ); 
			set &K_set.; 
			array K_col[*] _numeric_; 
			COL1 = K_col[_N_];
		run;
	'));

/*	%let tmp = %sysfunc( dosubl('*/
/*	data r(keep=x);*/
/*		call streaminit(&random_state.);*/
/*		do i=1 to &size.;*/
/*			x = rand("UNIFORM");*/
/*			output;*/
/*		end;*/
/*	run;*/
/*	'));*/
	%generate_random(&size., random_state = &random_state.);

	%let tmp = %sysfunc(dosubl('
		data orig_indices;
			do i=1 to &N.;
				output;
			end;
		run;
	'));

	%let sampl = ;
	%do it=1 %to &size.;
	
		%if &it <= &size. %then	
			%do;
				%let sel = %select_index_Chol(&rank., &it.);
				
				%if &it ne &size. %then %do;
					%adjust_probabilities_chol( &sel., &it., &size.);
				%end;
		%end;
	%end;

/*	%let tmp = %sysfunc( dosubl('*/
/*	data &dest.;*/
/*		set orig_indices;*/
/*		output;*/
/*		if _N_ = &rank. then stop;*/
/*	run;*/

	data &dest. (keep=Col1-Col&N.);
		array COL[&N.] (&N.*.);
		%do i=1 %to &size.;
			j = &i;
			set orig_indices(rename =(col1 = x)) point=j;
			put x=;
			Col[&i.] = x;
			put Col&i.=;
		%end;
		output;
		stop;
	run;
/*	'));*/
	proc datasets nolist lib=work;
	 delete d r orig_indices A;
	quit;


%mend;

