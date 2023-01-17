
%macro select_index_Chol(rank, it);
	%local sel;
	%let sell = %sysfunc(dosubl('
							data _null_; 
								a=&it.;
								set R point=a;
								cum_sum = 0;
								k=&it-1;
								d_point = &it;
								do while(x > cum_sum);
									set d point=d_point;
									cum_sum + abs( col1  / (&rank. - &it. + 1) );
									k + 1;
									d_point + 1;
								end;
								call symputx("sel", k, "L");
								stop;
							run;
						')
				);
		&sel.
%mend;

%macro adjust_probabilities_chol(sel, it, size);
	%let tmp = %sysfunc(dosubl('
proc iml;
	use A; 	read all  into A; 	close;
	if &sel.+1 <= &size. then do;
		A[&sel.+1:&size.,{&it.,&sel.}] = A[&sel.+1:&size.,{&sel., &it.}];
	end;
	A[{&it., &sel.},{&it., &sel.}] = A[{&sel., &it.},{&sel., &it.}];
	if &it.> 1 then do;
		A[{&it., &sel.}, 1:&it.-1] = A[{&sel., &it.}, 1:&it.-1];

	end;
	use orig_indices; 	read all  into orig_indices; 	close;
	orig_indices[{&it., &sel.},] = orig_indices[{&sel.,&it.},];
	use d; 	read all  into d; 	close;
	d[{&it., &sel.}] = d[{&sel., &it.}];
	A[&it., &it.] = sqrt(d[&it.]);

	if &it > 1 then do;
		A[&it. + 1:&size., &it.] = A[&it. + 1:&size., &it.] - A[&it. + 1:&size., 1:&It.-1]*A[&it., 1:&it.-1]`;
	end;
	A[&it. + 1:&size., &it.] = A[&it. + 1:&size., &it.]/A[&it., &it.];
	d[&it + 1:&size.] = d[&it + 1:&size.] - A[&it. + 1:&size., &it.]##2;

	create A from A;
					append from A;
	close A;
	create orig_indices from orig_indices;
					append from orig_indices;
	close orig_indices;
	create d from d;
					append from d;
	close d;

quit;
'));
%mend;


%macro proj_dpp_sampler_kernel_Chol(K_set, dest, size=., random_state=.);
	%let N = %size( &K_set. );
	%let rank = %trace( &K_set. );
	%put &=rank;
	%if &size = . %then %let size = &rank.;
	%if &random_state = . %then %let random_state = %random_int();
	%let tmp = %sysfunc(dosubl('
		data A;
			set &K_set.;
		run;
	'));

	%let tmp = %sysfunc( dosubl('
		data d( keep=COL1 ); 
			set &K_set.; 
			array K_col[*] _numeric_; 
			COL1 = K_col[_N_];
		run;
	'));

/*	%let tmp = %sysfunc( dosubl('*/
/*	data r(keep=x);*/
/*		call streaminit(&random_state.);*/
/*		do i=1 to &size.;*/
/*			x = rand("UNIFORM");*/
/*			output;*/
/*		end;*/
/*	run;*/
/*	'));*/
	%generate_random(&size., random_state = &random_state.);

	%let tmp = %sysfunc(dosubl('
		data orig_indices;
			do i=1 to &N.;
				output;
			end;
		run;
	'));

	%let sampl = ;
	%do it=1 %to &size.;
	
		%if &it <= &size. %then	
			%do;
				%let sel = %select_index_Chol(&rank., &it.);
				
				%if &it ne &size. %then %do;
					%adjust_probabilities_chol( &sel., &it., &size.);
				%end;
		%end;
	%end;

/*	%let tmp = %sysfunc( dosubl('*/
/*	data &dest.;*/
/*		set orig_indices;*/
/*		output;*/
/*		if _N_ = &rank. then stop;*/
/*	run;*/

	data &dest. (keep=Col1-Col&N.);
		array COL[&N.] (&N.*.);
		%do i=1 %to &size.;
			j = &i;
			set orig_indices(rename =(col1 = x)) point=j;
			put x=;
			Col[&i.] = x;
			put Col&i.=;
		%end;
		output;
		stop;
	run;
/*	'));*/
	proc datasets nolist lib=work;
	 delete d r orig_indices A;
	quit;


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


%macro proj_dpp_sampler_kernel_GS( K_set, dest, size=., random_state=.);
	%local N rank tmp avail sampl;
	%let N = %size( &K_set. );
	%let rank = %trace( &K_set. );

	%if &size. = . %then %let size = &rank.;
	%if &random_state. = . %then %do;
		%let random_state = %random_int();
	%end;


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

/*	%let tmp = %sysfunc( dosubl('*/
/*	data r(keep=x);*/
/*		call streaminit(&random_state.);*/
/*		do i=1 to &size.;*/
/*			x = rand("UNIFORM");*/
/*			output;*/
/*		end;*/
/*	run;*/
/*	'));*/
	%generate_random(&size., random_state = &random_state.);

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

	proc datasets nolist lib=work;
	 delete norm2 norm22 C r;
	quit;
	

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



%macro proj_dpp_sampler_eig_GS(eig_vecs, dest, size=., random_state=.);
  	%local N rank avail sampl avail_indices;

	%let N = %size(&eig_vecs.);
	%let rank = %numeric_var_count(&eig_vecs.);

	%if &size. = . %then
		%let size = &rank.;

	%if &random_state = . %then %let random_state = %random_int();

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
		norm2 = norm2/&rank.;

		create norm2 from norm2; append from norm2; close norm2;
		create norm22 from norm22; append from norm22; close norm22;
		create C from C; append from C; close C;
	quit;


	%generate_random(&size., random_state = &random_state.);	

	
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
	proc datasets nolist lib=work;
	 delete norm2 norm22 C R;
	quit;

	
%mend;


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

%macro select_index_Schur(rank, it);
	%local sel;
	%let sell = %sysfunc(dosubl('
							data _null_; 
								a=&it.;
								set R point=a;
								cum_sum = 0;
								k=0;
								do while(x > cum_sum);
									set Schur_comp;
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




%macro adjust_probabilities_Schur(K_set, avail_indices, sel, it, rank, sampl,sampl_new);
	%let i = %scan(&sampl., 1); 
	%put it &it.;
	%put sel &sel.;


	proc iml;
		use &K_set.; 	read all  into K; 		close;
		use K_inv;     		read all  into K_inv; 		close;
		use Schur_comp;		read all into Schur_comp;	close;

		if &it.=1 then do;
			K_inv[1, 1] = 1.0 / K[&sel., &sel.];


		end;
		else if &it.=2 then do;
		
			a = K[&sel., &sel.];
	    	b = -K[ &sel. , &i.];
			c = -K[ &i. , &i.];
			temp1 = a||b;
			temp2 = b||c;
			tempK_inv = temp1//temp2;
			K_inv[{1,2},{1,2}] = tempK_inv/(K[&i., &i.] * K[&sel., &sel.] - K[&sel., &i.]##2) ;
	
		end;



		else do;
			
			temp = K_inv[1:&it.-1, 1:&it.-1] * K[{&sampl.},&sel.];
			schur_sel = K[&sel., &sel.] - K[&sel. , {&sampl.}]*temp;
/* 			print schur_sel;*/
			temp2 = temp/schur_sel;

			K_inv[1:&it.-1, 1:&it.-1] = K_inv[1:&it.-1, 1:&it.-1] + temp*temp2`;
			K_inv[1:&it.-1, &it.] = -temp2 ; 
            K_inv[&it., 1:&it.-1] = K_inv[1:&it.-1, &it.]`;
            K_inv[&it., &it.] = 1/schur_sel;
			


			end;

		K_iY = K[{&avail_indices.},{&sampl_new.}];
/*		print K_iY;*/
		tmp = K_iY * K_inv[1:&it.,1:&it.];
/*		print tmp;*/
		tmp2 = tmp * K_iY`;
/*		print tmp2;*/
		tmp3 = tmp2[,+];
/*		print tmp3;*/
		Schur_comp[{&avail_indices.}] = (vecdiag(K)[{&avail_indices.}] - tmp3);
		Schur_comp[&sel.]=0;
		Schur_comp[{&avail_indices.}] = abs(Schur_comp[{&avail_indices.}]) / (&rank. - &it.);

		/*Schur_comp[{&avail_indices.}]*/


/*
        K_iY = K[np.ix_(avail, sampl[:it + 1])]
        schur_comp[avail] = vecdiag(K)[{&avail.}]\
                        - inner1d(K_iY.dot(K_inv[:it+1, :it+1]), K_iY, axis=1)
*/


		create K_inv from K_inv;
			append from K_inv;
		close K_inv; 
		create Schur_comp from Schur_comp;
			append from Schur_comp;
		close Schur_comp;
		



	quit;
%mend;



%macro proj_dpp_sampler_kernel_Schur(K_set, dest, size=., random_state=.);
	%let N = %size( &K_set. );
	%let rank = %trace( &K_set. );

	%if &size = . %then %let size = &rank.;

	%if &random_state = . %then %let random_state = %random_int();

	%let avail = 1; 
	%do i=2 %to &N.;
		%let avail = &avail. 1;
	%end;

	%let tmp = %sysfunc( dosubl('
		data schur_comp( keep=COL1 ); 
			set &K_set.; 
			array K_col[*] _numeric_; 
			COL1 = K_col[_N_] / &rank.; 
		run;'));


	%let tmp = %sysfunc( dosubl('
	data K_inv( drop=k );array COL[&size.] (&size.*0);do k=1 to &size.; output;end;run;'));
	

/*		%let tmp = %sysfunc( dosubl('*/
/*	data r(keep=x);call streaminit(&random_state.);do i=1 to &size.;x = rand("UNIFORM");output;end;run;'));*/
	%generate_random(&size., random_state = &random_state.);


	%let sampl = ;
	%do it=1 %to &size.;
		/*pick a new item to add to sample*/
		%let sel = %select_index_Schur(&rank., &it.);
		%let i = %scan(&sampl.,1);
		%let avail = %alter_avail(&avail., &sel.);
		%put avail &avail.;
	 	%let avail_indices = %get_avail_indices( &avail.);
		%put avail_indices &avail_indices.;
		%let sampl_new = &sampl. &sel. ;
		%if &it ne &size. %then %do;
			%adjust_probabilities_Schur( &K_set., &avail_indices., &sel., &it., &rank., &sampl., &sampl_new.);
		%end;
		%let sampl = &sampl_new. ;
		

		
	%end;

	data &dest. (keep=Col1-Col&N.);
		array COL[&N.] (&N.*.);
		%do i=1 %to &size.;
			Col[&i.] = %scan( &sampl., &i. );
		%end;
	run;

	proc datasets nolist lib=work;
	 delete schur_comp K_inv r ;
	quit;

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
	proc datasets nolist lib=work;
	 delete R;
	quit;
%mend;


%macro dpp_sampler_generic_kernel(K_set, dest, random_state=1, size=.);
	%if &size. = . %then %let size = %size(&K_set.);
	%generate_random(&size.);	
	%sampler_generic(&K_set., &dest., &size.);
%mend;


options minoperator mindelimiter=',';


%macro proj_dpp_sampler_kernel(kernel,
								dest ,
								mode=GS ,
								size=.,
								random_state=.);

%if &random_state. = . %then %let random_state = %random_int();

%if &size. ne . %then %do;
	%let rank = %trace( &kernel. );
	%if &size. > &rank. %then %do;
		%put ERROR: size k=&size. > rank=&rank.; 
		%goto exit;

	%end;
%end;

%if &mode.=GS %then %do;
	%proj_dpp_sampler_kernel_GS(&kernel.,&dest., size = &size., random_state=&random_state.);
%end; 
%else %if &mode. = Chol %then %do;
	%proj_dpp_sampler_kernel_Chol(&kernel.,&dest., size = &size., random_state=&random_state.);
%end;

%else %if &mode. = Schur %then %do;
	%proj_dpp_sampler_kernel_Schur(&kernel.,&dest., size = &size., random_state=&random_state.);
%end;

%else %do;
	%let EMEXCEPTIONSTRING = ERROR;
	%put &em_codebar;
	%put ERROR: Invalid sampling mode, choose from: GS, Chol, Schur. Given: &mode.; 
	%put &em_codebar;
	%goto doendm;
%end;
%doendm:
%mend;

%macro proj_dpp_sampler_eig(eig_vecs,
								dest ,
								mode=GS ,
								size=.,
								random_state=.);

%if &random_state. = . %then %let random_state = %random_int();


%if &mode.=GS %then %do;
	%proj_dpp_sampler_eig_GS(&eig_vecs., &dest., size = &size., random_state = &random_state.);
%end; 
%else %if &mode. = GS_bis %then %do;
	%proj_dpp_sampler_eig_GS_bis(&eig_vecs., &dest., size = &size., random_state = &random_state.);
%end;

%else %if &mode. = KuTa12 %then %do;
	%proj_dpp_sampler_eig_KuTa12(&eig_vecs., &dest., size = &size., random_state = &random_state.);
%end;

%else %do;
	%let EMEXCEPTIONSTRING = ERROR;
	%put &em_codebar;
	%put ERROR: Invalid sampling mode, choose from: GS, GS_bis, KuTa12. Given: &mode.; 
	%put &em_codebar;
	%goto doendm;
%end;
%doendm:
%mend;

%macro K_eig_vecs(eig_vals,eig_vecs);
proc iml;
	use &eig_vals.; read all into e_vals; close;
	use &eig_vecs.; read all into e_vecs; close;
	e_vals = e_vals`;
	V = {};
	do i=1 to ncol(e_vals);
		if e_vals[i] > 0.5 then V = V || i;
	end;
	V = e_vecs[1:nrow(e_vecs),V];
	create proj_dpp_V from V; append from V;
quit;
%mend;

%macro sample_exact(kernel_type ,K_kernel = ., L_kernel = .,K_e_vals=., K_e_vecs=.,L_e_vals=., L_e_vecs=., dest = sample, mode=GS, random_state=., projection=N, size=.);
/*%init(&kernel_type. ,K_kernel = &K_kernel., L_kernel = &L_kernel.,K_e_vals=&K_e_vals., K_e_vecs=&K_e_vecs.,L_e_vals=&L_e_vals., L_e_vecs=&L_e_vecs., dest = &dest., mode=&mode., random_state=&random_state., projection=&projection., size=&size.);*/

%if &random_state. = . %then %let random_state = %random_int();
%if &mode. = Schur %then %do;
	%put NOTE: Schur;
	%if &kernel_type. = correlation & &projection. = Y %then %do;
		%computeK(K_kernel=K_kernel., L_kernel = &L_kernel. ,L_e_vals=&L_e_vals., L_e_vecs=&L_e_vecs.,K_e_vals=&K_e_vals., K_e_vecs=&K_e_vecs.);
		%if &K_kernel. = . %then %let K_kernel = K;
		%proj_dpp_sampler_kernel(&K_kernel., &dest.,mode=&mode., random_state = &random_state., size = &size.);
	%end;
	%else %do;
		%let EMEXCEPTIONSTRING = ERROR;
		%put &em_codebar;
		%put ERROR: Schur sampling mode is only available for kernel_type=correlation and projection=1. Given kernel_type=&kernel_type. and projection=&projection.;
		%put &em_codebar;
		%goto doendm;
	%end;

%end;

%else %if &mode. = Chol %then %do;
	
	%computeK(K_kernel=&K_kernel., L_kernel = &L_kernel. ,L_e_vals=&L_e_vals., L_e_vecs=&L_e_vecs.,K_e_vals=&K_e_vals., K_e_vecs=&K_e_vecs.);
	%if &K_kernel. = . %then %let K_kernel = K;
	%if &kernel_type. = correlation & &projection. = Y %then %do;
		%put NOTE: Chol proj;
		%proj_dpp_sampler_kernel(&K_kernel., &dest.,mode=&mode., random_state = &random_state.,  size = &size.);
	%end;
	%else %do;
		%put NOTE: Chol;
		%dpp_sampler_generic_kernel(&K_kernel., &dest., random_state = &random_state.,  size = .);
	%end;
%end;



%else %if &K_e_vals. ne . %then %do;
	%put NOTE: K_e_vals ne ., eig sampler;
	%if &kernel_type. = correlation & &projection. = Y %then %do;
		%K_eig_vecs(&K_e_vals.,&K_e_vecs.);
	%end;
	%else %do;
		%dpp_eig_vecs_selector(&K_e_vals.,&K_e_vecs., random_state = &random_state.);
	%end;
	%proj_dpp_sampler_eig(proj_dpp_V, &dest., mode = &mode.,random_state = &random_state.,  size = &size.);
/*	proc datasets nolist lib=work;*/
/*		delete proj_dpp_V;*/
/*	run;*/
	
%end;

%else %if &L_e_vals. ne . %then %do;
	%put NOTE: L_e_vals ne ., eig sampler;
	%compute_K_evals_from_L_evals(&L_e_vals.);
	%let K_e_vals = K_e_vals;
	%let K_e_vecs = &L_e_vecs.;
	%if &kernel_type. = correlation & &projection. = Y %then %do;
		%K_eig_vecs(&K_e_vals.,&K_e_vecs.);
	%end;
	%else %do;
		%dpp_eig_vecs_selector(&K_e_vals.,&K_e_vecs., random_state = &random_state.);
	%end;
	%proj_dpp_sampler_eig(proj_dpp_V, &dest., mode = &mode.,random_state = &random_state.,  size = &size.);
	proc datasets nolist lib=work;
		delete proj_dpp_V K_e_vals;
	run;
%end;

%else %if (&K_kernel. ne .) and (&projection. = Y) %then %do;
	%put NOTE: K_kernel ne . and  p=1;
	%proj_dpp_sampler_kernel(&K_kernel., &dest.,mode=&mode., random_state = &random_state.,  size = &size.);
%end;
%else %if &K_kernel. ne .  %then %do;
	%put NOTE: K_kernel ne . ;
	%get_eigendecomposition(&K_kernel., eig_vals, eig_vecs);
	%let  K_e_vals = eig_vals;
	%let  K_e_vecs = eig_vecs;
	%dpp_eig_vecs_selector(&K_e_vals.,&K_e_vecs., random_state = &random_state.);
	%proj_dpp_sampler_eig(proj_dpp_V, &dest., mode = &mode.,random_state = &random_state.,  size = &size.);
	proc datasets nolist lib=work;
		delete proj_dpp_V eig_vals eig_vecs;
	run;
%end;
%else %if &L_kernel. ne . %then %do;
	%put NOTE: L_kernel ne .;
	%get_eigendecomposition(&L_kernel., eig_vals, eig_vecs);
	%let  L_e_vals = eig_vals;
	%let  L_e_vecs = eig_vecs;

	%compute_K_evals_from_L_evals(&L_e_vals.);
	%let K_e_vals = K_e_vals;
	%let K_e_vecs = &L_e_vecs.;
	%if &kernel_type. = correlation & &projection.=Y %then %do;
		%K_eig_vecs(&K_e_vals.,&K_e_vecs.);
	%end;
	%else %do;
		%dpp_eig_vecs_selector(&K_e_vals.,&K_e_vecs., random_state = &random_state.);
	%end;
	%proj_dpp_sampler_eig(proj_dpp_V, &dest., mode = &mode.,random_state = &random_state.,  size = &size.);
	proc datasets nolist lib=work;
		delete proj_dpp_V K_e_vals eig_vals eig_vecs;
	run;
%end;

%else %do;
	%let EMEXCEPTIONSTRING = ERROR;
	%put &em_codebar;
	%put ERROR: None of the available samplers could be used based on the current DPP representation.;
	%put &em_codebar;
	%goto doendm;
%end;




proc datasets nolist lib=work;
	delete computed_K;
run;

%doendm:
%mend;

%macro compute_K_from_eig(K_e_vals, K_e_vecs);
	proc iml;
		use &K_e_vals.; 	read all  into K_e_vals; 		close;
		use &K_e_vecs.; 	read all  into K_e_vecs; 		close;
		if dimension(K_e_vals)[1] ^= 1 then K_e_vals = K_e_vals`;
		K = (K_e_vecs # K_e_vals)*K_e_vecs`;
		create K from K;
				append from K;
		close K;
	quit;
%mend;
%macro compute_K_evals_from_L_evals(L_e_vals);
	proc iml;
		use &L_e_vals.; 	read all  into L_e_vals; 		close;
		if dimension(L_e_vals)[1] ^= 1 then L_e_vals = L_e_vals`;
		K_e_vals = L_e_vals / (1 + L_e_vals);
		create K_e_vals from K_e_vals;
				append from K_e_vals;
		close K_e_vals;
	quit;
%mend;

%macro computeK(K_kernel = ., L_kernel = .,K_e_vals=., K_e_vecs=.,L_e_vals=., L_e_vecs=.);
	%if &K_kernel. ne . %then %do;
		%put &=K_kernel;
		%put NOTE: K already computed;
		%goto exit;
	%end;
	%else %if (&K_e_vals. ne .) & (&K_e_vecs. ne .) %then %do;
		%compute_K_from_eig(&K_e_vals., &K_e_vecs.);
		%let K_kernel = K;
	%end;
	%else %if &L_e_vals. ne . %then %do;
		%put liczy K z L;
		%compute_K_evals_from_L_evals(&L_e_vals.);
		%compute_K_from_eig(K_e_vals, &L_e_vecs.);
		proc datasets nolist lib=work;
		delete K_e_vals;
		run;
		%let K_kernel = K;
	%end;
	%else %if &L_kernel. ne . %then %do;
		%compute_K(&L_kernel., K);
		%let K_kernel = K;
	%end;
	%exit:
%mend;

%macro dpp_eig_vecs_selector(eig_vals, eig_vecs,random_state=.);
	proc iml;
		use &eig_vals.; 	read all  into eig_vals; 		close;
		CALL STREAMINIT(&random_state.);
		if dimension(eig_vals)[1] ^= 1 then eig_vals = eig_vals`;
		y = {};
		do i = 1 to dimension(eig_vals)[2];
			x = rand("uniform");
			if x < eig_vals[i] then do;
				y = y || i;
			end;
		end;
		
		use &eig_vecs.; 	read all  into eig_vecs; 	close;
		eig_vecs = eig_vecs[,y];
		create proj_dpp_V from eig_vecs; append from eig_vecs;
	quit;
%mend;










