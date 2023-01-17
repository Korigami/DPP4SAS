%macro matrix_K(size, rank);

proc iml;
	call streaminit(0);
	x = j(&rank.*&size.,1);
	do i=1 to &rank.*&size.;
		x[i] = rand('uniform');
	end;
	x = shape(x,&size.,&rank.);
	call qr(q, r, piv, lindep, x);
	I = diag(j(&size.,1))[,1:&rank.];
	q = q*I;
	K = q*q`;
	create K from K;
			append from K;
	close K; 

quit;
%mend;




%macro size(set);
	%local k rc;

	%let rc = %sysfunc(dosubl('
			data _null_; 
				set &set. end=end;
				if _N_ = 1 then stop;
			run; 
			%let k =&sysnobs;
		'));
    &k
%mend;

%macro choice_bounded_wo_repl(bound, size, random_state);
	%local tmp;
	%let tmp = %sysfunc(dosubl('
		data r(keep = x);
			call streaminit(&random_state.);
			array R[&bound.] _temporary_ (&bound.*0);
			do while(no_of_sel < &size.);
				x = round(&bound.*rand("UNIFORM")+0.5);
				if R[x] ne 1 then
					do;
						R[x] = 1;
						no_of_sel + 1;
						output;
					end;
			end;
		run;
		'));
%mend;



%macro det_from_submatrix(matrix, indices);
	%local det tmp;
	
	%let tmp = %sysfunc(dosubl('
	proc iml;
		use &matrix.; read all into M; close;
	   	d = det(M[{&indices.},{&indices.}]);
		call symputx("det", char(d));
	quit;'));

	&det.
%mend;


%macro intersect1d(ar1, ar2); /* we assume that both ar1 and ar2 only contain variable x */
	%local intersection tmp;
	%let tmp = %sysfunc(dosubl('
		proc sql noprint;
			select A.x 
				into: intersection separated by ' '
			from &ar1. as A, &ar2. as B
			where A.x = B.x;
		quit;
	'));
	&intersection.
%mend;

%macro get_var_into_macrovar(dataset, variable);
	%local tmp result;

	%let tmp = %sysfunc(dosubl('
		proc sql noprint;
			select &variable.
				into: result separated by ' '
			from &dataset.;
		quit;
	'));

	&result.
%mend;




%macro initialize_AED_sampler(kernel, dest, random_state);

	%local N s0 det_s0 nb_trials tol T tmp;

	%let N = %size(&kernel.);
	
	%let tmp = %sysfunc(dosubl('
		data ground_set(keep = x);
			do i=1 to &N.;
				x =  i;
				output;
			end;	
		run;
	'));
	
	%let det_s0 = 0;
	%let nb_trials = 100;
	%let tol = 0.00000001;
	
	%do i=1 %to &nb_trials.;
		%if &det_s0. > &tol. %then 
			%do;
				%let tmp = %sysfunc(dosubl('
					proc delete data = ground_set R; quit;
					'));
					%let num_obs = %length(&s0.) - %length(%sysfunc(compress(&s0.))) + 1;
				data &dest. (keep=sample);
				%do i=1 %to &num_obs.;
				sample = %scan( &s0. , &i. ) ;
				output;
			%end;
		run;
				%return;
			%end;
		%else
			%do;
				%choice_bounded_wo_repl(%eval(2*&N.), &N., &random_state.);
				%let s0 = %intersect1d(ground_set, R);
			   	%let det_s0 = %det_from_submatrix(&kernel., &s0.);
			%end;

	%end; 

%mend;







%macro initialize_AD_and_E_sampler(kernel, size, dest, random_state);
   
	%local N s0 det_s0 nb_trials tol T tmp new_size;

	%let N = %size(&kernel.);
	
	%let tmp = %sysfunc(dosubl('
		data ground_set(keep = x);
			do i=1 to &N.;
				x = i;
				output;
			end;	
		run;
	'));
	
	%let det_s0 = 0;
	%let nb_trials = 100;
	%let tol = 0.00000001;
    
	%do i=1 %to &nb_trials.;
		%if &det_s0. > &tol. %then
			%do;
				%let tmp = %sysfunc(dosubl('
					proc delete data = ground_set R; quit;
					'));
										%let num_obs = %length(&s0.) - %length(%sysfunc(compress(&s0.))) + 1;
				data &dest. (keep=sample);
				%do i=1 %to &num_obs.;
					sample = %scan( &s0. , &i. ) ;
				output;
				%end;
				run;

				%return;
			%end; 
		%else
			%do;
				%if &size. ne . %then 
					%let new_size = &size.;
				%else 
					%let new_size = %sysfunc(round( 
									%sysevalf( %eval(&N.+1)*%sysfunc(rand(UNIFORM))+0.5  ) 
								));

				%choice_bounded_wo_repl(&N., &new_size., &random_state.);
				%let s0 = %get_var_into_macrovar(R, x);
				%let det_s0 = %det_from_submatrix(&kernel., &s0.);
				%put &s0;
			%end;
	%end;

%mend;





options mprint;


%macro RandBetween(min, max);
	%let x = %sysfunc(ranuni(0));
   %sysevalf((&min + %sysfunc(floor((1+&max-&min)*&x.))))
%mend;



%macro sampl_without_given(N, s, random_state);
%let ids = ;
%let s_size = %sysevalf(%length(&s.)-%length(%sysfunc(compress(&s.))) + 1); 
%let pom = 0;
%do i=1 %to &N.;
	%do j=1 %to &s_size.;
		%if %scan(&s., &j.) = &i. %then %do;
			%let pom=1;
		%end;

	%end;
	%if &pom. = 0 %then %do;
		%let ids = &ids. &i.;
	%end;
	%else %do;
	%let pom = 0;
	%end;
%end;

%let size_wo = %sysevalf(%length(&ids.)-%length(%sysfunc(compress(&ids.))) + 1); 
%let id  = %sysfunc(round(%sysevalf(&size_wo*%sysfunc(ranuni(&random_state.))+0.5)));
%let out =  %scan(&ids., &id.);
&out.
%mend;





%macro add_exchange_delete_sampler(kernel, s_init, dest, random_state=123, nb_iter=10);

	%let N = %size(&kernel.);
	%let tmp = %sysfunc(dosubl('
		data ground_set(keep = x);
			do i=1 to &N.;
				x =  i;
				output;
			end;	
		run;
	'));

	%let s0 = %get_var_into_macrovar(&s_init., sample);
	%let det_s0 = %det_from_submatrix(&kernel., &s0.);
	%let size_s0 = %sysevalf(%length(&s0.)-%length(%sysfunc(compress(&s0.))) + 1);
	%put &size_s0;
	%put inicjalne &s0.;

	%do idx=1 %to &nb_iter. ;
		%put &idx.;
		%let s1 = &s0.;
		/*wybór indeksu do ewentualnego usuniêcia*/
		%let s_ind = %sysfunc(round(%sysevalf(&size_s0.*%sysfunc(ranuni(&random_state))+0.5)));
		/*losowy indeks spoœród nienale¿¹cych do s0*/
		
		%let t = %sampl_without_given(&N., &s0., &random_state.);
		%let U = %sysfunc(ranuni(&random_state.));
		%let ratio = %sysevalf(&size_S0. / &N.); 
		%put &U.;
		
		%if &U. < %sysevalf(0.5 * (1 - &ratio.)**2) %then %do;

			%put ADD;
			%let s1 = &s1. &t.;
			%let det_s1 = %det_from_submatrix(&kernel., &s1.);
			%let rand_num = %sysfunc(ranuni(&random_state));
			%if &rand_num. < %sysevalf(&det_s1./&det_s0. * (&size_s0. + 1)/(&N. - &size_s0)) %then %do;
				%let s0 = &s1.;
				%let det_s0 = &det_s1;
				%let size_s0 = %sysevalf(&size_s0 + 1);
			%end;


		%end;
		%else %if (%sysevalf(0.5 * (1-&ratio.)**2) <= &U.) and ( &U. < %sysevalf(0.5 * (1-&ratio.))) %then %do;
			%put EXCHANGE;
			%put indeks &s_ind.;
			%let el_to_del = %sysfunc(scan(&s1.,&s_ind));
			%put element &el_to_del.;
			%let s_pom = ;
			%put &s1.;
			
			%do it = 1 %to &size_s0.;
				%if (%sysfunc(scan(&s1. , &it.)) ne  &el_to_del.) %then %do;
					%let s_pom = &s_pom. %sysfunc(scan(&s1.,&it.)); 
					%end;
			%end;
			%let s1 = &s_pom.;
			%put po odjeciu &s1.;
			%let s1 = &s1. &t.;
			%put po dodaniu &s1.;
			%let det_s1 = %det_from_submatrix(&kernel., &s1.);
			%let rand_num = %sysfunc(ranuni(&random_state));
			%if &rand_num. < %sysevalf(&det_s1./&det_s0.) %then %do;
				%let s0 = &s1.;
				%let det_s0 = &det_s1;
			%end;		
		%end;
		%else %if (%sysevalf(0.5*(1-&ratio.))<=&U.) and (&U. < %sysevalf(0.5 * (&ratio.**2 + (1 - &ratio.)))) %then %do;
			%put DELETE;
			%put indeks &s_ind.;
			%let el_to_del = %sysfunc(scan(&s1.,&s_ind));
			%put element &el_to_del.;
			%let s_pom = ;
			
				%do it = 1 %to &size_s0.;
					%if (%sysfunc(scan(&s1., &it.)) ne  &el_to_del.) %then %do;
						%let s_pom = &s_pom. %sysfunc(scan(&s1.,&it.)); 
					%end;
				%end;
			
			%let s1 = &s_pom.;
			%put po odjeciu &s1.;
			%let det_s1 = %det_from_submatrix(&kernel., &s1.);
			%let rand_num = %sysfunc(ranuni(&random_state));
			%if &rand_num. < %sysevalf(&det_s1./&det_s0. * &size_s0. / (&N. - (&size_s0-1))) %then %do;
				%let s0 = &s1.;
				%let det_s0 = &det_s1;
				%let size_s0 = %sysevalf(&size_s0 - 1);
				
			%end;	
		%end;

		%else %do;
			%put NOTHING;
			%let s0 = &s0.;

		%end;
		%put &s0.;
	%end;

			data &dest. (keep=sample);
			%do i=1 %to &size_s0.;
				sample = %scan( &s0. , &i. ) ;
				output;
			%end;
		run;




%mend;




%macro add_delete_sampler(kernel, s_init, dest, random_state=123, nb_iter=10);
	%let N = %size(&kernel.);
	%let s0 = %get_var_into_macrovar(&s_init., sample);
	%let det_s0 = %det_from_submatrix(&kernel., &s0.);
	%let size_s0 = %sysevalf(%length(&s0.)-%length(%sysfunc(compress(&s0.))) + 1);
	%do idx=1 %to &nb_iter. ;
		%let rand_num = %sysfunc(ranuni(&random_state));
		%if &rand_num. < 0.5 %then %do;
			%let s1= &s0.;
			%let size_s1 = %sysevalf(%length(&s1.)-%length(%sysfunc(compress(&s1.))) + 1);
			%let s = %sysfunc(round(%sysevalf(&N.*%sysfunc(ranuni(&random_state.))+0.5)));
			%let tmp = 0;
			%put s &s.;
			%do it=1 %to &size_s1.;
				%put %sysfunc(scan(&s1., &it.));
				%if %sysfunc(scan(&s1., &it.))=&s. %then %do;
					%let tmp = 1;
					%put hej;
					%end;
					
			%end;
			%put tmp &tmp.;
			%if &tmp. = 1 %then %do;
				%put remove;
			/*remove s*/
				
				%let s_pom = ;
			
				%do it = 1 %to &size_s1.;
				%put %sysfunc(scan(&s1., &it.)) ^=  &s.; 
					%if (%sysfunc(scan(&s1., &it.)) ^=  &s.) %then %do;
						%let s_pom = &s_pom. %sysfunc(scan(&s1.,&it.)); 
					%end;
				%end;
				%let s1 = &s_pom.;

			
			%end;
			%else %do;
			/*add s*/
				%let s1 = &s1. &s.;

			%end;
			%let det_s1 = %det_from_submatrix(&kernel., &s1.);
			%let rand_num2 = %sysfunc(ranuni(&random_state));

			%if &rand_num2. < %sysevalf(&det_s1./&det_s0.) %then %do;
				%let s0 = &s1.;
				%put zmiana;
				%let det_s0 = &det_s1;
			%end;

			%else %do;
				%let s0 = &s0.;
			%end;


		%end;
		%put &s0.;
	
	%end;
	%put &s0.;
	%let size_s0 = %sysevalf(%length(&s0.)-%length(%sysfunc(compress(&s0.))) + 1);

	data &dest. (keep=sample);
		%do i=1 %to &size_s0.;
			sample = %scan( &s0. , &i. ) ;
			output;
		%end;
	run;


%mend;
options mprint;



%macro basic_exchange_sampler(kernel, s_init, dest, random_state=123, nb_iter=10);
	%let N = %size(&kernel.);
	%let s0 = %get_var_into_macrovar(&s_init., sample);
	%let det_s0 = %det_from_submatrix(&kernel., &s0.);
	%let size_s0 = %sysevalf(%length(&s0.)-%length(%sysfunc(compress(&s0.))) + 1);
	%do idx=1 %to &nb_iter. ;
		%let rand_num = %sysfunc(ranuni(&random_state));
		%if &rand_num. < 0.5 %then %do;
			%let s1= &s0.;
			%let s_ind = %sysfunc(round(%sysevalf(&size_s0.*%sysfunc(ranuni(&random_state))+0.5)));
			%let t = %sampl_without_given(&N., &s0., &random_state.);
			%let el_to_del = %sysfunc(scan(&s1.,&s_ind));
			%let s_pom = ;
			%do it = 1 %to &size_s0.;
				%if (%sysfunc(scan(&s1. , &it.)) ne  &el_to_del.) %then %do;
					%let s_pom = &s_pom. %sysfunc(scan(&s1.,&it.)); 
					%end;
			%end;
			%let s1 = &s_pom.;
			%put po odjeciu &s1.;
			%let s1 = &s1. &t.;
			%put po dodaniu &s1.;
			%let det_s1 = %det_from_submatrix(&kernel., &s1.);
			%let rand_num2 = %sysfunc(ranuni(&random_state));

			%if &rand_num2. < %sysevalf(&det_s1./&det_s0.) %then %do;
				%let s0 = &s1.;
				%put zmiana;
				%let det_s0 = &det_s1;
			%end;

			%else %do;
				%let s0 = &s0.;
			%end;




		%end;
		%put &s0.;
	%end;
	%put &s0.;
	data &dest. (keep=sample);
		%do i=1 %to &size_s0.;
			sample = %scan( &s0. , &i. ) ;
			output;
		%end;
	run;
%mend;

/*
%matrix_K(10, 8);

%initialize_AED_sampler(K,s0, 12321564);
%add_exchange_delete_sampler(K,s0,dest1,random_state=123, nb_iter=10);

%initialize_AD_and_E_sampler(K, 5, s0, 1233);
%add_delete_sampler(K,s0,dest2);
%basic_exchange_sampler(K,s0,dest3);
*/
