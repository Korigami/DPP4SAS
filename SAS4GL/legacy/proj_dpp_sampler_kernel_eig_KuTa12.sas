%macro adjust_probabilities_eig_KuTa12( eig_vecs,  sel, it, rank);
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

		norm2 = abs(norm2)/(&rank. - &it.);
		create norm2 from norm2;
						append from norm2;
		close norm2; 
	quit;
	'));
%mend;


%macro proj_dpp_sampler_eig_KuTa12(eig_vecs, sample, size, random_state);
  %local N rank avail sampl avail_indices;

	%let N = %size(&eig_vecs.);
	%let rank = %numeric_var_count(&eig_vecs.);

	%if &size. = . %then
		%let size = &rank.;


	proc iml;
		use &eig_vecs.; read all into V; close;
		norm2 = shape(0, &N., 1);
		do i=1 to &N.;
			norm2[i] = V[i,]*V[i,]`;
		end;

		create norm2 from norm2; append from norm2; close norm2;
	quit;


	%generate_random(&size., &random_state.);	

	
	%let sampl = ;
	%do it=1 %to &size.; 
	
			%let sel = %select_index_GS(&rank., &it.);
			%let sampl = &sampl. &sel. ;
			
			%if &it ne &size. %then 
				%do;
					%adjust_probabilities_eig_KuTa12( &eig_vecs.,  &sel., &it., &rank.);
				%end;
	%end;

	data &sample.;
		array COL[&N.] (&N.*.);
		%do i=1 %to &size.;
			COL[&i.] = %scan(&sampl., &i.);
		%end;
		
	run;

	
%mend;

%get_eigendecomposition(K, eig_vals, eig_vecs);
%proj_dpp_sampler_eig_KuTa12(eig_vecs, sample, size=., random_state=1);



