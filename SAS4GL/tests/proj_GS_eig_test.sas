options spool;



%macro get_eigendecomposition(matrix, eig_vals, eig_vecs);

	proc iml;
		use &matrix.; read all into A; close;

		call eigen(val, rvec, A) vecl="lvec";
		
		create &eig_vals. from val; append from val;
		create &eig_vecs. from rvec; append from rvec;
	quit;	

%mend;

%get_eigendecomposition(K, eval, evec);
data evec2; 
	set evec(keep=Col1 Col2 Col3 Col4 Col5 Col6 Col7 Col8);
run;

/* we assume that K has trace = 8 */

%macro generate_random(size, random_state);
	/* generating set with size observables of random numbers from uniform distribution on (0,1) */
	data r(keep=x);
		call streaminit(&random_state.);
		do i=1 to &size.;
			x = rand("UNIFORM");
			output;
		end;
	run;
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
	%local avail2 i el;
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
	%local ind i el;
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

%macro select_index_GS(rank, it);
	%local sell sel;
	%let sell = %sysfunc(dosubl('
							data _null_; 
								a=&it.;
								set R point=a;
								cum_sum = 0;
								k=0;
								do while(x > cum_sum);
									set norm2;
									cum_sum + abs( col1);
									k + 1;
								end;
								call symputx("sel", k, "L");
								stop;
							run;
						')
				);
		&sel.
%mend;


%macro adjust_probabilities_eig( eig_vecs, avail_indices, sel, it, rank);
  %local tmp;
	%let tmp = %sysfunc( dosubl('
	proc iml;
			use &eig_vecs.; read all  into V; 		close;
			use C;     		read all  into C; 		close;
			use norm2;		read all  into norm2;	close;
			use norm22;		read all  into norm22;	close;
			if &it. > 1 then 
				do;
					C[{&avail_indices.}, &it.] 
						= (
							V[{&avail_indices.},] * V[&sel.,]`
							- C[{&avail_indices.}, 1:&it.-1]*C[&sel., 1:&it.-1]`
						   )/ sqrt( norm22[&sel.]);
				end;
			else  
				do;
					C[{&avail_indices.}, &it.]  
						= V[{&avail_indices.},] * V[&sel.,]` /sqrt( norm22[&sel.]);
				end;
			norm22[{&avail_indices.}] 
					= norm22[{&avail_indices.}] 
						- C[{&avail_indices.}, &it.]##2;
			norm2[&sel.] = 0;
			norm2[{&avail_indices.}] = abs(norm22[{&avail_indices.}]) / (&rank. - &it.);
			create C from C;
					append from C;
			close C; 
			create norm2 from norm2;
					append from norm2;
			close norm2; 
			create norm22 from norm22;
					append from norm22;
			close norm22; 
	quit;
	'));
%mend;



%macro proj_dpp_sampler_eig_GS(eig_vecs, dest, size, random_state);
  	%local N rank avail sampl avail_indices;

	%let N = %size(&eig_vecs.);
	%let rank = %numeric_var_count(&eig_vecs.);

	%if &size. = . %then
		%let size = &rank.;

	%let avail = 1; 
	%do i=2 %to &N.;
		%let avail = &avail. 1;
	%end;

	proc iml;
		use &eig_vecs.; read all into V; close;
		norm2 = shape(0, &N., 1);
		C = shape(0, &N., &size.);
		do i=1 to &N.;
			norm2[i] = V[i,]*V[i,]`;
		end;

		norm22 = norm2;

		create norm2 from norm2; append from norm2; close norm2;
		create norm22 from norm22; append from norm22; close norm22;
		create C from C; append from C; close C;
	quit;


	%generate_random(&size., &random_state.);	

	
	%let sampl = ;
	%do it=1 %to &size.; 
	
			%let sel = %select_index_GS(&rank., &it.);
			%let sampl = &sampl. &sel. ;
			%let avail = %alter_avail(&avail., &sel.);
			%let avail_indices = %get_avail_indices( &avail.);
			
			%if &it ne &size. %then 
				%do;
					%adjust_probabilities_eig( &eig_vecs., &avail_indices., &sel., &it., &rank.);
				%end;
	%end;
	

	data &dest. (keep=Col1-Col&N.);
		array COL[&N.] (&N.*.);
		%do i=1 %to &size.;
			Col[&i.] = %scan( &sampl., &i. );
		%end;
	run;
	

	
%mend;


%proj_dpp_sampler_eig_GS(evec2, dest = res1, size = ., random_state = 113)
%proj_dpp_sampler_eig_GS(evec2, dest = res2, size = ., random_state = 11342)
%proj_dpp_sampler_eig_GS(evec2, dest = res3, size = ., random_state = 1512)
%proj_dpp_sampler_eig_GS(evec2, dest = res4, size = ., random_state = 4651)
%proj_dpp_sampler_eig_GS(evec2, dest = res5, size = ., random_state = 1543)
%proj_dpp_sampler_eig_GS(evec2, dest = res6, size = ., random_state = 134)
%proj_dpp_sampler_eig_GS(evec2, dest = res7, size = ., random_state = 1132)
%proj_dpp_sampler_eig_GS(evec2, dest = res8, size = ., random_state = 1543)
%proj_dpp_sampler_eig_GS(evec2, dest = res9, size = ., random_state = 1123)
%proj_dpp_sampler_eig_GS(evec2, dest = res10, size = ., random_state = 1123)
%proj_dpp_sampler_eig_GS(evec2, dest = res11, size = ., random_state = 154)

data full_res_eig_GS;
		set res1 res2 res3 res4 res5 res6 res7 res8 res9 res10 res11;
	run;

%get_emp_kernel_from_sample(full_res_eig_GS);



%alter_kernel_for_comparison(K, K_com);

%plot_matrix(K_com, Original Kernel);
%plot_matrix(emp_ker_from_full_res_eig_GS, Empirical Kernel);

