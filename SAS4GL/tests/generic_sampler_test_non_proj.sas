%macro size(set);
	%local k;

	%let rc = %sysfunc(dosubl('
			data _null_; 
				set &set. end=end;
				if _N_ = 1 then stop;
			run; 
			%let k =&sysnobs;
		'));
    &k
%mend;

%put %size(full_res);

%macro numeric_var_count(set);

	%local no_of_vars;

	%let rc = %sysfunc(dosubl('
		data _null_; 
			set &set; 
			array vars_num _numeric_;
			call symputx('no_of_vars', dim(vars_num), "L"); 
			stop; 
		run;
		'));
	&no_of_vars.
%mend;

%put %numeric_var_count(full_res);


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

%create_identity_matrix( 5);

%macro generate_random(size);
	%let tmp = %sysfunc( dosubl('
	data r(keep=x);
		call streaminit(&random_state.);
		do i=1 to &size.;
			x = rand("UNIFORM");
			output;
		end;
	run;
	'));
%mend;



%macro sampler_generic(K_set, dest, size);
	proc iml;
		use &K_set.; read all into A; close;
		use R; read all into R; close;
		
		res = shape(., 1, &size.);

		do i=1 to &size.;
			
			if R[i] < A[i,i] then
				res[1,i] = i;
			else
				A[i, i] = A[i, i] -1;

			if i < &size. then	
				do;

					A[i+1:&size., i] = A[i+1:&size., i] / A[i, i];
					A[i+1:&size., i+1:&size.] 
							= A[i+1:&size., i+1:&size.] 
								- A[i+1:&size., i]*A[i, i+1:&size.];
				end;
		end;
		
		create &dest. from res; append from res;
	quit;
%mend;


%macro dpp_sampler_generic_kernel(K_set, dest, random_state);

	%let size = %size(&K_set.);
	%generate_random(&size.);	
	%sampler_generic(&K_set., &dest., &size.);
%mend;




%macro matrix_K(size, rank);


%macro compute_K(L, K);

	proc iml;
		use &L.; read all into L; close;
		ones = I(nrow(L)[1]);
		K =  ones - inv(L+ones);
		create &K. from K; append from K; close K;
	quit;	

%mend;

proc iml;
	call streaminit(0);
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


%matrix_K(10, 8);

%compute_K(K, K2);

%macro test();
	%do i=1 %to 1000;
		%dpp_sampler_generic_kernel(K2, res&i., &i.);
	%end;

	data full_res;
		set %do i=1 %to 1000; res&i. %end;;
	run;

	proc delete data = %do i=1 %to 1000; res&i. %end;; quit;
%mend;

%test;




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
	/* rysuje wykres macierzy */
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


%get_emp_kernel_from_sample(full_res);

%macro alter_kernel_for_comparison(kernel, com_kernel);
	/* zamienia ujemne wartosci kernela na dodatnie */
	data &com_kernel.(drop = i j);
		set &kernel.;
		array Col[*] _numeric_;
		do i=1 to dim(col);
			if col[i] < 0 then
				col[i] = - col[i];
		end;
	run;
%mend;

%alter_kernel_for_comparison(K2, K2_com);
%plot_matrix(K2_com, Original Kernel);
%plot_matrix(emp_ker_from_full_res, Empirical Kernel);
