



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
