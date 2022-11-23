
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


%macro proj_dpp_sampler_eig_GS(eig_vecs, sample, size, random_state);
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

	data &sample.;
		array COL[&N.] (&N.*.);
		%do i=1 %to &size.;
			COL[&i.] = %scan(&sampl., &i.);
		%end;
		
	run;

	
%mend;

