%macro select_index_KuTa12(rank, it);
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

%macro adjust_probabilities_eig_KuTa12( eig_vecs,  sel, it, rank, sampl);
	%let tmp = %sysfunc( dosubl('
	proc iml;
		use &eig_vecs.; read all  into V; 		close;
		k = loc(V[&sel.,] ^= 0)[1];
		tmp = (V[, k] / (V[&sel, k])*V[&sel., ]);
		V = V- tmp;
		idx = remove(1:ncol(V),k);
		call qr(q, r, piv, lindep, V[,idx]);
		m = min(ncol(q), ncol(r));
		I = diag(j(ncol(q),1))[,1:ncol(r)];
		q = q*I;
		V = q*q`;
		norm2 = shape(0, nrow(V), 1);
		do i=1 to nrow(V);
				norm2[i] = V[i,]*V[i,]`;
		end;
		norm2[{&sampl.}] =0;
		norm2 = abs(norm2)/(&rank. - &it.);
		create &eig_vecs. from V;
						append from V;
		close &eig_vecs.;
		create norm2 from norm2;
						append from norm2;
		close norm2; 
	quit;
	'));
%mend;


%macro proj_dpp_sampler_eig_KuTa12(eig_vecs, sample, size=., random_state=.);
  %local N rank avail sampl avail_indices;

	%let N = %size(&eig_vecs.);
	%let rank = %numeric_var_count(&eig_vecs.);
	%put &=rank;

	%if &size. = . %then
		%let size = &rank.;

	%if &random_state = . %then %let random_state = %random_int();

	proc iml;
		use &eig_vecs.; read all into V; close;
		norm2 = shape(0, &N., 1);
		do i=1 to &N.;
			norm2[i] = V[i,]*V[i,]`;
		end;
		norm2 = norm2/&rank.;

		create norm2 from norm2; append from norm2; close norm2;
	quit;


	%generate_random(&size., random_state = &random_state.);	

	
	%let sampl = ;
	%do it=1 %to &size.; 
	
			%let sel = %select_index_KuTa12(&rank., &it.);
			%let sampl = &sampl. &sel. ;
			
			%if &it ne &size. %then 
				%do;
					%adjust_probabilities_eig_KuTa12( &eig_vecs.,  &sel., &it., &rank., &sampl.);
				%end;
	%end;

	data &sample.;
		array COL[&N.] (&N.*.);
		%do i=1 %to &size.;
			COL[&i.] = %scan(&sampl., &i.);
		%end;
		
	run;
	
	proc datasets nolist lib=work;
	 delete norm2 R;
	quit;

	
%mend;

