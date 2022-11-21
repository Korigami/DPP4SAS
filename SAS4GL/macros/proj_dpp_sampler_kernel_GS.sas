
%macro select_index_GS(rank, it);
	
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
								call symputx("sel", k, "G");
								stop;
							run;
						')
				);
		&sel.
%mend;



%macro adjust_probabilities_GS_exact( K_set, avail_indices, sel, it, rank);
	%let tmp = %sysfunc( dosubl('
	proc iml;
			use &K_set.; 	read all  into K; 		close;
			use C;     		read all  into C; 		close;
			use norm2;		read all into norm2;	close;
			if &it. > 1 then do;
			C[{&avail_indices.}, &it.] 
					= (
						K[{&avail_indices.}, &sel.] 
						- C[{&avail_indices.}, 1:&it.-1]*C[&sel., 1:&it.-1]`
					   )/ sqrt( norm2[&sel.]);
			end;
			else do;
			C[{&avail_indices.}, &it.]  = K[{&avail_indices.}, &sel.] /sqrt( norm2[&sel.]);
			end;
			norm2[{&avail_indices.}] 
					= norm2[{&avail_indices.}] 
						- C[{&avail_indices.}, &it.]##2;
			norm2[&sel.] = 0;
			norm2[{&avail_indices.}] = norm2[{&avail_indices.}] / sum( abs( norm2));
			create C from C;
					append from C;
			close C; 
			create norm2 from norm2;
					append from norm2;
			close norm2; 
	quit;
	'));
%mend;


%macro proj_dpp_sampler_kernel_GS( K_set, dest, size=., random_state=23354);

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
			
			%if &it ne &size. %then %do;

			%adjust_probabilities_GS_exact( &K_set., &avail_indices., &sel., &it., &rank.);
			%end
		%end; 
		%end;
	%end;
	
	%put &sampl.;

		data &dest. (keep=sample);
			%do i=1 %to &size.;
				sample = %scan( &sampl. , &i. ) ;
				output;
			%end;
		run;
	

%mend;



/*
%proj_dpp_sampler_kernel_GS(K, sampl);
*/
