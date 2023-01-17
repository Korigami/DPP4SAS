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


%macro carries_sample(base,size=100, dest=work.carries_sample, random_state=.);

%if &random_state=. %then %let random_state=0;


data rand_ints;
call streaminit(&random_state);
retain cum_sum;
do i=1 to &size.;
x= floor(&base. * rand("UNIFORM")) + 1 /* rand('integer',1,&base.); */ ;
cum_sum +x;
rest = mod(cum_sum, &base.);
output;
end;
drop i;

run;

%let rests = %get_var_into_macrovar(rand_ints,rest);


%let carries = ;

%do i = 2 %to &size.;
%if %sysfunc(scan(&rests., %sysevalf(&i.-1)))>
			%sysfunc(scan(&rests.,&i.))%then %do;

	%let carries = &carries. %sysevalf(&i.-1);
	
%end;
%end;
%put &carries.;
%let carries_size=%sysevalf(%length(&carries)-%length(%sysfunc(compress(&carries.))) + 1);


data &dest. (keep=sample);
	%do i=1 %to &carries_size.;
		sample = %scan( &carries., &i. ) ;
		output;
	%end;

run;
%mend;



%macro descent_sampler(size=100, dest=work.descent_sample,  random_state=.);

%if &random_state=. %then %let random_state=0;


data permut;
call streaminit(&random_state);
do i=1 to &size.;
num=i;
rand = rand("uniform",0,1);
output;
end;
drop i ;

run;



proc sort data=permut;
by rand;

run;

%let permut = %get_var_into_macrovar(permut, num);
%put &permut.;

%let descent = ;

%do i = 2 %to &size.;
%if %sysfunc(scan(&permut., %sysevalf(&i.-1)))>
			%sysfunc(scan(&permut.,&i.))%then %do;

	%let descent = &descent. %sysevalf(&i.-1);
	
%end;
%end;
%put &descent.;


%let descent_size=%sysevalf(%length(&descent.)-%length(%sysfunc(compress(&descent.))) + 1);


data &dest. (keep=sample);
	%do i=1 %to &descent_size;
		sample = %scan( &descent., &i. ) ;
		output;
	%end;

run;


%mend;




%macro virtual_descent_sampler(x0=0.5, size=100, dest= work.virt_descent_sample, random_state=.);

%if &random_state=. %then %let random_state=0;


data permut;
call streaminit(&random_state);
do i=1 to &size.+1;
num=i;
rand = rand("uniform",0,1);
output;
end;
drop i ;

run;



proc sort data=permut;
by rand;

run;

%let permut = %get_var_into_macrovar(permut, num);

%let X = ;

%do i = 2 %to &size.+1;
%if %sysfunc(scan(&permut.,&i.-1))>
			%sysfunc(scan(&permut.,&i.))%then %do;

	%let X = &X. 1;
	
%end;
%else %do;
	%let X = &X. 0;

%end;
%end;
%put &X.;


data Y;
call streaminit(&random_state);
do i=1 to &size.+1;
binom=rand('BINOMIAL', &x0., 2);
if binom ^= 1 then y=1;
else y=0;
output;
end;

run;

%let Y = %get_var_into_macrovar(Y, y);

%let virtual_descent = ;
%do i=1 %to &size.;
/*if (~Y[i] and Y[i + 1]) or (~Y[i] and ~Y[i + 1] and X[i])*/
%if (%sysfunc(scan(&Y,&i))=0 & %sysfunc(scan(&Y,&i.+1))=1) | (%sysfunc(scan(&Y,&i))=0 & %sysfunc(scan(&Y,&i.+1))=0 & %sysfunc(scan(&X,&i.))=1) %then %do;
%let virtual_descent = &virtual_descent. %sysevalf(&i-1);
%end;


%end;
%put &virtual_descent.;

%let virt_descent_size=%sysevalf(%length(&virtual_descent.)-%length(%sysfunc(compress(&virtual_descent.))) + 1);


data &dest. (keep=sample);
	%do i=1 %to &virt_descent_size;
		sample = %scan( &virtual_descent., &i. ) ;
		output;
	%end;

run;


%mend;


%macro stat_1_dep_sampler(size=100, dest=work.stat_sample, mode=descent, base=., x0=0.5, random_state=.);

%if ((&mode ^= descent) & (&mode ^= carries) & (&mode ^= virtual)) %then %do;
	%put Mode must be one of the following: 'descent', 'carries' or 'virtual';
%end;

%if &mode=carries %then %do;
	%if &base=. | &base < 2 %then %do;
		%put Base must be provided (integer > 1);
	%end;
	%else %do;
		%carries_sample(base=&base., size=&size., dest=&dest., random_state=&random_state.);
	%end;

%end;

%else %if &mode=descent %then %do;
	%descent_sampler(size=&size., dest=&dest., random_state=&random_state.);	
%end;

%else %if &mode=virtual %then %do;
	%if ((&x0=.) | (&x0<0) | (&x0>1)) %then %do;
		%put x_0 must be provided : 0<x0<1;

	%end;
	%else %do;
		%virtual_descent_sampler(x0=&x0, size=&size, dest=&dest., random_state=&random_state);
	%end;
%end;



%mend;
