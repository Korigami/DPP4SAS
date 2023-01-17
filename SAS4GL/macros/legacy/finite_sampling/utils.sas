

%macro numeric_var_count(set);
	%local no_of_vars rc;
	%let rc = %sysfunc(dosubl('
		data _null_; 
			set &set; 
			array vars_num _numeric_;
			call symputx("no_of_vars", dim(vars_num), "L"); 
			stop; 
		run;
		'));
	&no_of_vars.
%mend;

%macro trace(K);
	%local tr trace;
	%let tr=%sysfunc(
			dosubl('data _null_; set &K. end=end; 
					array columns[*] _numeric_; 
					trace + columns[_N_];
					if end then call symputx("trace", round(trace), "L");
					run;'));
	&trace
%mend;

%macro size(set);
		%local rc k;
        %let rc = %sysfunc(
			dosubl('data _null_; set &set end=end;run; %let k =&sysnobs;')
		);
        &k
%mend;

%macro alter_avail( avail, sel);
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



/* template for the matrix graphing */
proc template;              
define statgraph heatmap;
dynamic _X _Y _Z _T;       
 begingraph;
 entrytitle _T;             
  layout overlay;
    heatmapparm x=_X y=_Y colorresponse=_Z /  
       name="heatmap" primary=true
       xbinaxis=false ybinaxis=false;  
    continuouslegend "heatmap";
  endlayout;
endgraph;
end;
run;

%macro plot_matrix(matrix, title);
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
	   dynamic _X='i' _Y='j' _Z='x' _T="&title.";
	run;

	proc delete data = to_plot;
	run;
%mend;

%macro alter_kernel_for_comparison(kernel, com_kernel);
	data &com_kernel.(drop = i j);
		set &kernel.;
		array Col[*] _numeric_;
		do i=1 to dim(col);
			if col[i] < 0 then
				col[i] = - col[i];
		end;
	run;
%mend;

%macro compute_K(L, K);

	proc iml;
		use &L.; read all into L; close;
		ones = I(dimension(L)[1]);
		K =  ones - inv(L+ones);
		create &K. from K; append from K; close K;
	quit;	

%mend;

%macro compute_L(K, L);
	proc iml;
		
		use &K.; read all into K; close;
		no_of_rows = dimension(K)[1];
		if rank(K) = no_of_rows then 
			L = K*inv( I(no_of_rows) - K);
		else
			L = .;
		create &L. from L; append from L; close L;
	quit;	

%mend;

%macro get_eigendecomposition(matrix, eig_vals, eig_vecs);
	proc iml;
		use &matrix.; read all into A; close;

		call eigen(val, rvec, A) vecl="lvec";
		
		create &eig_vals. from val; append from val;
		create &eig_vecs. from rvec; append from rvec;
	quit;	

%mend;

%macro generate_random(size, random_state=1);
	/* generating set with size observables of random numbers from uniform distribution on (0,1) */
	data r(keep=x);
		call streaminit(&random_state.);
		do i=1 to &size.;
			x = rand("UNIFORM");
			output;
		end;
	run;
%mend;

%macro create_identity_matrix( size);
	data unit&size.;
		%do i = 1 %to &size.;
			%do j=1 %to &size.;
				%if &i. = &j. %then
					col&j. = 1;
				%else
					col&j. = 0;
				;
			%end;
			output;
		%end;
	run;
%mend;

%macro matrix_K(size, rank);
proc iml;
	call streaminit(0);
	call streamrewind;
	x = j(&rank.*&size.,1);
	do i=1 to &rank.*&size.;
		x[i] = rand('uniform');
	end;
	x = shape(x,&size.,&rank.);
	call qr(q, r, piv, lindep, x);
	I = diag(j(&size.,1))[,1:&rank.];
	q = q*I;
	K = q*q`;
	create K from K;
			append from K;
	close K; 
quit;
%mend;


%macro random_int(n=.);
%local r;
%if &n = . %then %let n=100000;
%let r = %sysfunc(round(%sysevalf(%sysfunc(ranuni(0))*&n)));
&r.
%mend;
