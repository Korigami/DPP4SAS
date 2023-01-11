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



%macro sample_poissonized_plancherel(theta, dest=work.poissonized_sample,random_state=.);


%if &theta <=0  %then %do;
	%put Theta must an integer greater then 0; 
	%goto exit;

%end;

%if &random_state=. %then %do;
	%let random_state=0;

%end;



data _null_;
	call streaminit(&random_state);
	N=rand('poisson',&theta.);
	call symputx("N", N);

run;



data sigma;
call streaminit(&random_state);
do i=1 to &N.;
num=i;
rand = rand("uniform",0,1);
output;
end;
drop i ;

run;



proc sort data=sigma;
by rand;

run;

%let sigma = %get_var_into_macrovar(sigma, num);

/*RSK*/
%let num_row=1;
%let P1 = ;


%do i=1 %to &N.;
	%let x = %sysfunc(scan(&sigma., &i.));
	%let if_leave=0;
	
	%do j=1 %to &num_row.;

		%if %sysevalf(%length(&&P&j.)) =0 %then %let size&j =1;
 
		%else %let size&j = %sysevalf(%length(&&P&j.)-%length(%sysfunc(compress(&&P&j.))) + 1); 

		%if &x. >= %sysfunc(scan(&&P&j., &&size&j)) %then %do;

			%let P&j. = &&P&j. &x.;

			%let if_leave=1;
			%goto leave;
		%end;


		%else %do;

			%let pom = 0;
			%let id = 1;
			%let prev=1;
			%let P_tmp = ;
			%do %while (&pom.=0);
				%let element = %sysfunc(scan(&&P&j., &id.));

				%if &x. > &element. %then %do;
					%let id = %sysevalf(&id. + 1);
				%end;
				%else %do;
					%do k=1 %to &&size&j;
						%if &k.=&id. %then %do;
							%let P_tmp = &P_tmp &x;
						%end;
						%else %do;
							%let P_tmp = &P_tmp %sysfunc(scan(&&P&j., &k));

						%end;

					%end;
					%let x = &element;
					%let P&j = &P_tmp;
					%let pom = 1;

				%end;

				
				
			%end;
				

		%end;

		

	%end;

	%if &if_leave=0 %then %do;
		%let num_row = %sysevalf(&num_row + 1);
		%let P&num_row = &x;

	%end;


	%leave:



%end;

%put &num_row;
%let young_diag=;
%let sampl = ;
%do number=1 %to &num_row; 

	%put p_num &number;
	%put &&P&number;
	%let young_diag_num = /*&young_diag. */%sysevalf(%length(&&P&number.)-%length(%sysfunc(compress(&&P&number.))) + 1);
	%let tmp_num= /*&tmp_arr */%sysevalf(&number. - 0.5);
	%let sampl = &sampl %sysevalf(&young_diag_num. - &tmp_num);
	%let young_diag = &young_diag &young_diag_num;

%end;
%put &=young_diag.;
%put &=sampl;

%let sampl_size = %sysevalf(%length(&sampl.)-%length(%sysfunc(compress(&sampl.))) + 1);


data &dest. (keep=sample);
	%do i=1 %to &sampl_size;
		young = %scan( &young_diag., &i. ) ;
		tmp_num = &i. - 0.5;
		sample = young -tmp_num;
		output;
	%end;

run;

%mend;

