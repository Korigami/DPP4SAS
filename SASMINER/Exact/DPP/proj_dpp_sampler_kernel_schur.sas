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


/*%matrix_K(10, 8);*/
/*%proj_dpp_sampler_kernel_Schur(K, test1);*/



