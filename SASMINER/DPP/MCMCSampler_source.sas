options minoperator mindelimiter=',';

%macro choice_bounded_wo_repl(bound, size, random_state);
	%local tmp;
	%let tmp = %sysfunc(dosubl('
		data r(keep = x);
			call streaminit(&random_state.);
			array R[&bound.] _temporary_ (&bound.*0);
			do while(no_of_sel < &size.);
				x = floor(&bound.*rand("UNIFORM")+1);
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




%macro initialize_AED_sampler(kernel, dest, size = .,  random_state=.);
	

	%local N s0 det_s0 nb_trials tol T tmp uni_iter;
	%let uni_iter = 0;
	%if &random_state = . %then %let random_state = %random_int();

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
		%if %sysevalf(&det_s0. > &tol.) %then 
			%do;
				%let tmp = %sysfunc(dosubl('
					proc delete data = ground_set R; quit;
					'));
				%let num_obs = %sysevalf(%length(&s0.) - %length(%sysfunc(compress(&s0.))) + 1);
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
				%choice_bounded_wo_repl(%eval(2*&N.), &N., %eval(&random_state. * &i.));
				%let s0 = %intersect1d(ground_set, R);
				%put &=s0;
			   	%let det_s0 = %det_from_submatrix(&kernel., &s0.);
			%end;
/**/
	%end;
	%put ERROR: Initialization terminated unsuccessfully. After &nb_trials. random trials, no initial set S0 satisfies det L_S0 > &tol.;
%exit:

%mend;







%macro initialize_AD_and_E_sampler(kernel, dest, size = ., random_state=.);
   
	%local N s0 det_s0 nb_trials tol T tmp new_size uni_iter;
	%let uni_iter = 0;

	%if &random_state = . %then %let random_state = %random_int();

	%let N = %size(&kernel.);
	%put &=n;
	
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
		%put ---------------------------------;
		%put NOTE: Iteration number &i.;
		%put &=det_s0 &=tol;
		%if %sysevalf(&det_s0. > &tol.) %then
			%do;
				%let tmp = %sysfunc(dosubl('
					proc delete data = ground_set R; quit;
					'));
				%let num_obs = %sysevalf(%length(&s0.) - %length(%sysfunc(compress(&s0.))) + 1);
				%put check &=num_obs;
				%put &=dest;
				data &dest. (keep=sample);
					%do j=1 %to &num_obs.;
						sample = %scan( &s0. , &j. ) ;
						output;
					%end;
				run;

				%return;
			%end; 
		%else
			%do;
				%put ----3-----;
				%put &=size;
				%if &size. > 0  %then %do;
					%let new_size = &size.;
				%end;
				%else %do;
					%let uni_iter = %eval(&uni_iter + 1); 
					%let r_u = %random_uni(%eval(&uni_iter +&random_state.));
					%let new_size = %sysfunc(floor( 
									%sysevalf( %eval(&N.)*&r_u.+1  ) 
								));
				%end;

				%put &=new_size;
				%choice_bounded_wo_repl(&N., &new_size., %eval(&random_state.+&i.));
				%let s0 = %get_var_into_macrovar(R, x);
				%put ---------------------------------;
				%put &=s0;
				%let det_s0 = %det_from_submatrix(&kernel., &s0.);
		%end;
	%end;
	%exit:
%mend;








%macro RandBetween(min, max);
	%let x = %sysfunc(ranuni(0));
   %sysevalf((&min + %sysfunc(floor((1+&max-&min)*&x.))))
%mend;



%macro sampl_without_given(N, s, random_state);
%local uni_ter;
%let uni_iter = 0;
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
%let uni_iter = %eval(&uni_iter + 1); 
%let r_u = %random_uni(%eval(&uni_iter +&random_state.));
%let id  = %sysfunc(floor(%sysevalf(&size_wo*&r_u.+1)));
%let out =  %scan(&ids., &id.);
&out.
%mend;





%macro add_exchange_delete_sampler(kernel, s_init, dest, random_state=., nb_iter=10);
	%local uni_ter;
	%let uni_iter = 0;
	%let N = %size(&kernel.);
	%let tmp = %sysfunc(dosubl('
		data ground_set(keep = x);
			do i=1 to &N.;
				x =  i;
				output;
			end;	
		run;
	'));
	%if &random_state = . %then %let random_state = %random_int();
	%let s0 = %get_var_into_macrovar(&s_init., sample);
	%let det_s0 = %det_from_submatrix(&kernel., &s0.);
	%let size_s0 = %sysevalf(%length(&s0.)-%length(%sysfunc(compress(&s0.))) + 1);
	%put &=size_s0;
	%put inicjalne &s0.;

	%do idx=1 %to &nb_iter. ;
/*		%put &=idx.;*/
		%let s1 = &s0.;
		/*selection of index tfor possible removal */
		%let uni_iter = %eval(&uni_iter + 1); 
		%let r_u = %random_uni(%eval(&uni_iter +&random_state.));
		%let s_ind = %sysfunc(floor(%sysevalf(&size_s0.*&r_u.+1)));
		/*random index outside of s0*/
		
		%let t = %sampl_without_given(&N., &s0., &random_state.);
/*		%put &=random_state.;*/
		%let uni_iter = %eval(&uni_iter + 1); 
		%let U = %random_uni(%eval(&uni_iter +&random_state.));
/*		%let U = %sysfunc(ranuni(&random_state.));*/
		%let ratio = %sysevalf(&size_S0. / &N.); 
/*		%put &=U.;*/
		
		%if %sysevalf(&U. < %sysevalf(0.5 * (1 - &ratio.)**2)) %then %do;

			%let s1 = &s1. &t.;
			%let det_s1 = %det_from_submatrix(&kernel., &s1.);
			%let uni_iter = %eval(&uni_iter + 1); 
			%let rand_num =  %random_uni(%eval(&uni_iter +&random_state.));
/*			%put &=rand_num;*/

			%let tmp = %sysevalf( &det_s1./&det_s0. * (&size_s0. + 1)/(&N. - &size_s0));
/*			%put eq = &tmp.;*/
			%if %sysevalf(&rand_num. < &tmp.) %then %do;
				%put ADD;
				%put OLD &=s0;
/*				%put element &t.;*/
				%let s0 = &s1.;
				%let det_s0 = &det_s1;
				%let size_s0 = %sysevalf(&size_s0 + 1);
			%end;
			%else %do;
				%PUT NOTHING;
			%end;


		%end;
		%else %if (%sysevalf(0.5 * (1-&ratio.)**2) <= &U.) and ( &U. < %sysevalf(0.5 * (1-&ratio.))) %then %do;
			%let el_to_del = %sysfunc(scan(&s1.,&s_ind));
			%let s_pom = ;
			
			%do it = 1 %to &size_s0.;
				%if (%sysfunc(scan(&s1. , &it.)) ne  &el_to_del.) %then %do;
					%let s_pom = &s_pom. %sysfunc(scan(&s1.,&it.)); 
					%end;
			%end;
			%let s1 = &s_pom.;
/*			%put po odjeciu &s1.;*/
			%let s1 = &s1. &t.;
/*			%put po dodaniu &s1.;*/
			%let det_s1 = %det_from_submatrix(&kernel., &s1.);
			%let uni_iter = %eval(&uni_iter + 1); 
			%let rand_num =  %random_uni(%eval(&uni_iter +&random_state.));
/*			%put &=rand_num;*/
/*			%put eq = %sysevalf(&det_s1./&det_s0.);*/
			%if %sysevalf(&rand_num. < %sysevalf(&det_s1./&det_s0.)) %then %do;
				%put EXCHANGE;
				%put OLD &=s0;
/*				%put indeks &s_ind.;*/
/*				%put element &el_to_del.;*/
/*				%put &s1.;*/
				%let s0 = &s1.;
				%let det_s0 = &det_s1;
			%end;
			%else %do;
				%PUT NOTHING;
			%end;	
		%end;
		%else %if (%sysevalf(0.5*(1-&ratio.))<=&U.) and (&U. < %sysevalf(0.5 * (&ratio.**2 + (1 - &ratio.)))) %then %do;
/*			%put indeks &s_ind.;*/
			%let el_to_del = %sysfunc(scan(&s1.,&s_ind));
/*			%put element &el_to_del.;*/
			%let s_pom = ;
			
				%do it = 1 %to &size_s0.;
					%if (%sysfunc(scan(&s1., &it.)) ne  &el_to_del.) %then %do;
						%let s_pom = &s_pom. %sysfunc(scan(&s1.,&it.)); 
					%end;
				%end;
			
			%let s1 = &s_pom.;
/*			%put po odjeciu &s1.;*/
			%let det_s1 = %det_from_submatrix(&kernel., &s1.);
			%let uni_iter = %eval(&uni_iter + 1); 
			%let rand_num =  %random_uni(%eval(&uni_iter +&random_state.));
/*			%put &=rand_num;*/

/*			%put eq = %sysevalf(&det_s1./&det_s0. * &size_s0. / (&N. - (&size_s0-1)));*/
			%if %sysevalf(&rand_num. < %sysevalf(&det_s1./&det_s0. * &size_s0. / (&N. - (&size_s0-1)))) %then %do;
				%put DELETE;
				%put OLD &=s0;
				%let s0 = &s1.;
				%let det_s0 = &det_s1;
				%let size_s0 = %sysevalf(&size_s0 - 1);
				
			%end;
			%else %do;
				%PUT NOTHING;
			%end;	
		%end;

		%else %do;
			%put NOTHING;
			%let s0 = &s0.;

		%end;
		%put &=s0.;
	%end;

			data &dest. (keep=sample);
			%do i=1 %to &size_s0.;
				sample = %scan( &s0. , &i. ) ;
				output;
			%end;
		run;




%mend;




%macro add_delete_sampler(kernel, s_init, dest, random_state=., nb_iter=10);
	%local uni_ter;
	%let uni_iter = 0;
	%let N = %size(&kernel.);
	%if &random_state = . %then %let random_state = %random_int();
	%let s0 = %get_var_into_macrovar(&s_init., sample);
	%let det_s0 = %det_from_submatrix(&kernel., &s0.);
	%let size_s0 = %sysevalf(%length(&s0.)-%length(%sysfunc(compress(&s0.))) + 1);
	%do idx=1 %to &nb_iter. ;
		%let uni_iter = %eval(&uni_iter + 1); 
		%let rand_num =  %random_uni(%eval(&uni_iter +&random_state.));
		%if %sysevalf(&rand_num. < 0.5) %then %do;
			%let s1= &s0.;
			%let size_s1 = %sysevalf(%length(&s1.)-%length(%sysfunc(compress(&s1.))) + 1);
			%let uni_iter = %eval(&uni_iter + 1); 
			%let r_u =  %random_uni(%eval(&uni_iter +&random_state.));
			%let s = %sysfunc(floor(%sysevalf(&N.*&r_u.+1)));
			%let tmp = 0;
			%put s &s.;
			%do it=1 %to &size_s1.;
/*				%put %sysfunc(scan(&s1., &it.));*/
				%if %sysfunc(scan(&s1., &it.))=&s. %then %do;
					%let tmp = 1;
/*					%put hej;*/
					%end;
					
			%end;
/*			%put tmp &tmp.;*/
			%if &tmp. = 1 %then %do;
/*				%put remove;*/
			/*remove s*/
				
				%let s_pom = ;
			
				%do it = 1 %to &size_s1.;
/*				%put %sysfunc(scan(&s1., &it.)) ^=  &s.; */
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
			%let uni_iter = %eval(&uni_iter + 1); 
			%let rand_num2 = %random_uni(%eval(&uni_iter +&random_state.));



			%if %sysevalf(&rand_num2. < %sysevalf(&det_s1./&det_s0.)) %then %do;
				%if &tmp. = 1 %then %put DELETE;
				%else %do;
				%put ADD;
				%end;
				%put OLD &=s0;
				%let s0 = &s1.;
				%put &=s0;
				%put zmiana;
				%let det_s0 = &det_s1;
			%end;

			%else %do;
				%put NOTHING;
				%put &s0.;
				%let s0 = &s0.;
			%end;


		%end;
		%else %do;
		%put NOTHING;
		%put &s0.;
		%end;
	
	%end;
	%put &s0.;
	%let size_s0 = %sysevalf(%length(&s0.)-%length(%sysfunc(compress(&s0.))) + 1);

	%put output creation &=dest;
	data &dest. (keep=sample);
		%do i=1 %to &size_s0.;
			sample = %scan( &s0. , &i. ) ;
			output;
		%end;
	run;


%mend;




%macro basic_exchange_sampler(kernel, s_init, dest, random_state=., nb_iter=10);
	%local uni_ter;
	%let uni_iter = 0;
	%let N = %size(&kernel.);
	%if &random_state = . %then %let random_state = %random_int();
	%let s0 = %get_var_into_macrovar(&s_init., sample);
	%let det_s0 = %det_from_submatrix(&kernel., &s0.);
	%let size_s0 = %sysevalf(%length(&s0.)-%length(%sysfunc(compress(&s0.))) + 1);
	%do idx=1 %to &nb_iter. ;
		%let uni_iter = %eval(&uni_iter + 1); 
		%let rand_num = %random_uni(%eval(&uni_iter +&random_state.));
		%if %sysevalf(&rand_num. < 0.5) %then %do;
			%let s1= &s0.;
			%let uni_iter = %eval(&uni_iter + 1); 
			%let r_u = %random_uni(%eval(&uni_iter +&random_state.));
			%let s_ind = %sysfunc(floor(%sysevalf(&size_s0.*&r_u.+1)));
			%let t = %sampl_without_given(&N., &s0., &random_state.);
			%let el_to_del = %sysfunc(scan(&s1.,&s_ind));
			%let s_pom = ;
			%do it = 1 %to &size_s0.;
				%if (%sysfunc(scan(&s1. , &it.)) ne  &el_to_del.) %then %do;
					%let s_pom = &s_pom. %sysfunc(scan(&s1.,&it.)); 
					%end;
			%end;
			%let s1 = &s_pom.;
/*			%put po odjeciu &s1.;*/
			%let s1 = &s1. &t.;
/*			%put po dodaniu &s1.;*/
			%let det_s1 = %det_from_submatrix(&kernel., &s1.);
			%let uni_iter = %eval(&uni_iter + 1); 
			%let rand_num2 = %random_uni(%eval(&uni_iter +&random_state.));

			%if %sysevalf(&rand_num2. < %sysevalf(&det_s1./&det_s0.)) %then %do;
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



%macro sample_mcmc(K=., s_init=.,dest = sample, sampling_mode = .,kernel_type = .,projection = ., random_state=., nb_iter=10, size=., input_initial = N);

	%if &sampling_mode. in (AED, AD, E) %then %do;
		%if &kernel_type. = likelihood %then %do;
			
			%dpp_sampler_mcmc(K=&K., mode=&sampling_mode., random_state=&random_state., s_init=&s_init., nb_iter=&nb_iter., T_max=., size2 = &size., dest=&DEST.,input_initial = &input_initial.);
		%end;
		%else %do;
			%compute_L(&K., comp_L);
			%if (%sysfunc(exist(comp_L))) %then %do;
/*				%let comp_L = comp_L;*/
/*				%put dpp_sampler_mcmc(K=&K., mode=&sampling_mode., random_state=&random_state., s_init=&s_init., nb_iter=&nb_iter., T_max=., size2 = &size., dest=&DEST.,input_initial = &input_initial.);*/
				
				%dpp_sampler_mcmc(dest=&DEST.,K=comp_L , mode=&sampling_mode., random_state=&random_state., s_init=&s_init., nb_iter=&nb_iter., T_max=., size2 = &size.,input_initial = &input_initial.);
					proc datasets nolist lib=work;
						delete comp_L;
					run;
			%end;
			%else %do;
				%put ERROR: K data set must have size = rank, to be computed to L;
				%goto exit;
			%end;
		%end;
	%end;
	%else %do;
		%put ERROR: Invalide mode;
	%end;
%exit:
%mend;


%macro dpp_sampler_mcmc(K=., mode=AED, random_state=., s_init=., nb_iter=10, T_max=.,size2 = ., dest=sample, input_initial = N);
%put ----2-----;
%put &=size2;
%let s = 0;
	%if &s_init. ne . %then %do;
		data s0;
			set &s_init.;
			array n[*] _numeric_;
			sample = n[1];
			keep sample;
		run;
		%let s_init = s0;
	%end; 

	%if &mode. = AED %then %do;
		%if &s_init. = . %then %do;
			%let s = 1;
			%initialize_AED_sampler(&K.,s0,size = &size2., random_state=&random_state. );
			%let s_init = s0;
			%put -----------+++++++++++++++++++++_________________;
		%end;
		%add_exchange_delete_sampler(&K.,&s_init.,&dest.,random_state=&random_state., nb_iter=&nb_iter.);
	%end;

	%if &mode. = AD %then %do;
		%if &s_init. = . %then %do;
			%let s = 1;
			
			%initialize_AD_and_E_sampler(&K., s0,size = &size2., random_state = &random_state.);
			%let s_init = s0;

		%end;
		%add_delete_sampler(&K.,&s_init.,&dest.,random_state=&random_state., nb_iter=&nb_iter.);
	%end;

	%if &mode. = E %then %do;
		%if &s_init. = . %then %do;
			%let s = 1;
			%initialize_AD_and_E_sampler(&K., s0,size = &size2.,  random_state = &random_state.);
			%let s_init = s0;
		%end;
		%basic_exchange_sampler(&K.,&s_init.,&dest.,random_state=&random_state., nb_iter=&nb_iter.);
	%end;
	proc datasets nolist lib=work;
		delete s0;
	run;
%exit:
%mend;





