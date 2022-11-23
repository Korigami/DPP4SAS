
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
									cum_sum + abs( col1 /* / (&rank. - &it. + 1) */);
									k + 1;
								end;
								call symputx("sel", k, "L");
								stop;
							run;
						')
				);
		&sel.
%mend;



%macro adjust_probabilities( K_set, avail_indices, sel, it, rank);
	%local tmp;
	%let tmp = %sysfunc( dosubl('
	proc iml;
			use &K_set.; 	read all  into K; 		close;
			use C;     		read all  into C; 		close;
			use norm2;		read all into norm2;	close;
			use norm22;		read all into norm22;	close;

			if &it. > 1 then do;
			C[{&avail_indices.}, &it.] 
					= (
						K[{&avail_indices.}, &sel.] 
						- C[{&avail_indices.}, 1:&it.-1]*C[&sel., 1:&it.-1]`
					   )/ sqrt( norm22[&sel.]);
			end;
			else do;
				C[{&avail_indices.}, &it.]  = K[{&avail_indices.}, &sel.] /sqrt( norm22[&sel.]);
			end;
/*			norm2[{&avail_indices.}] */
/*					= norm2[{&avail_indices.}] */
/*						- C[{&avail_indices.}, &it.]##2;*/
			norm22[{&avail_indices.}] 
					= norm22[{&avail_indices.}] 
						- C[{&avail_indices.}, &it.]##2;
			norm2[&sel.] = 0;
/*			norm22[&sel.] = 0;*/
/*			norm2[{&avail_indices.}] = norm22[{&avail_indices.}] / sum( abs( norm22));*/
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


%macro proj_dpp_sampler_kernel_GS( K_set, dest, size=., random_state=23354);
	%local N rank tmp avail sampl;
	%let N = %size( &K_set. );
	%let rank = %trace( &K_set. );

	%if &size = . %then %let size = &rank.;

	%let avail = 1; 
	%do i=2 %to &N.;
		%let avail = &avail. 1;
	%end;

	%let tmp = %sysfunc( dosubl('
	data c( drop=k );	
		array COL[&size.] (&size.*0);
		do k=1 to &N.; 
			output;
		end;
	run;'));

	%let tmp = %sysfunc( dosubl('
		data norm2( keep=COL1 ); 
			set &K_set.; 
			array K_col[*] _numeric_; 
			COL1 = K_col[_N_] / &rank.; 
		run;
	'));
	%let tmp = %sysfunc( dosubl('
		data norm22( keep=COL1 ); 
			set &K_set.; 
			array K_col[*] _numeric_; 
			COL1 = K_col[_N_];
		run;
	'));

	%let tmp = %sysfunc( dosubl('
	data r(keep=x);
		call streaminit(&random_state.);
		do i=1 to &size.;
			x = rand("UNIFORM");
			output;
		end;
	run;
	'));
	

	%let sampl = ;
	%do it=1 %to &size.;
	
		%if &it <= &size. %then	
		%do;

			%let sel = %select_index_GS(&rank., &it.);

			%let sampl = &sampl. &sel. ;
			%let avail = %alter_avail(&avail., &sel.);
			%let avail_indices = %get_avail_indices( &avail.);
			
			%if &it ne &size. %then 
				%do;

			%adjust_probabilities( &K_set., &avail_indices., &sel., &it., &rank.);
			%end; 
		%end;
	%end;
	
	data &dest. (keep=Col1-Col&N.);
		array COL[&N.] (&N.*.);
		%do i=1 %to &size.;
			Col[&i.] = %scan( &sampl. , &i. ) ;
		%end;
	run;
	

%mend;


/*
%macro test();
	%do i=1 %to 2;
		%proj_dpp_sampler_kernel_GS(K, dest = res&i., size = ., random_state = 123&i.);
	%end;


	data full_res_GS;
		set %do i=1 %to 2; res&i. %end;;
	run;

	proc delete data = %do i=1 %to 2; res&i. %end;; quit;
%mend;

%test;
*/


%put _user_;
%proj_dpp_sampler_kernel_GS(K, dest = res1, size = 5, random_state = 1)
%proj_dpp_sampler_kernel_GS(K, dest = res2, size = 5, random_state = 11342)
%proj_dpp_sampler_kernel_GS(K, dest = res3, size = 5, random_state = 1512)
%proj_dpp_sampler_kernel_GS(K, dest = res4, size = 5, random_state = 4651)
%proj_dpp_sampler_kernel_GS(K, dest = res5, size = 5, random_state = 1543)
%proj_dpp_sampler_kernel_GS(K, dest = res6, size = 5, random_state = 134)
%proj_dpp_sampler_kernel_GS(K, dest = res7, size = 5, random_state = 1132)
%proj_dpp_sampler_kernel_GS(K, dest = res8, size = 5, random_state = 1543)
%proj_dpp_sampler_kernel_GS(K, dest = res9, size = 5, random_state = 1123)
%proj_dpp_sampler_kernel_GS(K, dest = res10, size = 5, random_state = 21123)
%proj_dpp_sampler_kernel_GS(K, dest = res11, size = 5, random_state = 2154)
%proj_dpp_sampler_kernel_GS(K, dest = res12, size = 5, random_state = 211342)
%proj_dpp_sampler_kernel_GS(K, dest = res13, size = 5, random_state = 21512)
%proj_dpp_sampler_kernel_GS(K, dest = res14, size = 5, random_state = 24651)
%proj_dpp_sampler_kernel_GS(K, dest = res15, size = 5, random_state = 21543)
%proj_dpp_sampler_kernel_GS(K, dest = res16, size = 5, random_state = 2134)
%proj_dpp_sampler_kernel_GS(K, dest = res17, size = 5, random_state = 21132)
%proj_dpp_sampler_kernel_GS(K, dest = res18, size = 5, random_state = 21543)
%proj_dpp_sampler_kernel_GS(K, dest = res19, size = 5, random_state = 21123)
%proj_dpp_sampler_kernel_GS(K, dest = res20, size = 5, random_state = 12123)

data full_res_GS_sized;
		set res1 res2 res3 res4 res5 res6 res7 res8 res9 res10 
			res11 res12 res13 res14 res15 res16 res17 res18 res19 res20;
	run;

%get_emp_kernel_from_sample(full_res_GS_sized);



%alter_kernel_for_comparison(K, K_com);

%plot_matrix(K_com);
%plot_matrix(emp_ker_from_full_res_GS_sized);

