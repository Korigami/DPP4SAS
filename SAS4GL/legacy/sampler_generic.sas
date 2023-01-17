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
