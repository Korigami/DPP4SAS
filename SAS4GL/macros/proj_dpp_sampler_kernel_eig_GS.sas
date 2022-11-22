
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
			norm22[&sel.] = 0;
			norm2[&sel.] = 0;
			norm2[{&avail_indices.}] = norm22[{&avail_indices.}] / sum( abs( norm22));

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
