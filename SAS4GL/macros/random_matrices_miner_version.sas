%macro hermite_sampler_full(result_eigvals, size=10, beta=2, normalize=1, random_state=1618);

	proc iml;
		call streaminit(&random_state.);
		
		result_eigvals = shape(0, &size., 1);

		if &beta. = 1 then
			do;
				A = RandNormal( &size., shape(0, &size., 1), I(&size.));
				result_eigvals = eigval( A + A` )/sqrt(2) ;
			end;

		else if &beta. = 2 then
			do;
				A = RandNormal( &size., shape(0, &size., 1), I(&size.));
				B = RandNormal( &size., shape(0, &size., 1), I(&size.));

				Re_mat = A + A`;
				Im_mat = B - B`;

				collect_matr = shape(0, 2 * &size., 2 * &size.);
				/* first column */
				collect_matr[1:&size.,      1:&size.] = Re_mat;
				collect_matr[&size.+1:2*&size.,1:&size.] = Im_mat;
				/* second column */
				collect_matr[1:&size.,         &size.+1:2*&size.] = - Im_mat;
				collect_matr[&size.+1:2*&size.,&size.+1:2*&size.] = Re_mat;

				result_eigvals = eigval(collect_matr)[{%do i=1 %to &size.; %eval(2*&i.) %end;}] / sqrt(2); /* eigevalues get doubled --- we select only one from each pair */ 
			end;

		else if &beta. = 4 then
			do;
				A = RandNormal( &size., shape(0, &size., 1), I(&size.));
				B = RandNormal( &size., shape(0, &size., 1), I(&size.));
				C = RandNormal( &size., shape(0, &size., 1), I(&size.));
				D = RandNormal( &size., shape(0, &size., 1), I(&size.));

				Re_diag  =  A + A`;
				Re_ndiag = -C + C;
				Im_diag  =  B - B`;
				Im_ndiag =  D - D`;

				collect_matr = shape(0, 4 * &size., 4 * &size.);

				/* first column */
				collect_matr[1:&size.,            1:&size.] = Re_diag;
				collect_matr[&size.+1:2*&size.,   1:&size.] = Re_ndiag;
				collect_matr[2*&size.+1:3*&size., 1:&size.] = Im_diag;
				collect_matr[3*&size.+1:4*&size., 1:&size.] = Im_ndiag;
				/* second column */
				collect_matr[1:&size.,            &size.+1:2*&size.] = - Re_ndiag;
				collect_matr[&size.+1:2*&size.,   &size.+1:2*&size.] = Re_diag;
				collect_matr[2*&size.+1:3*&size., &size.+1:2*&size.] = Im_ndiag;
				collect_matr[3*&size.+1:4*&size., &size.+1:2*&size.] = - Im_diag;
				/* third column */
				collect_matr[1:&size.,            2*&size.+1:3*&size.] = -Im_diag;
				collect_matr[&size.+1:2*&size.,   2*&size.+1:3*&size.] = -Im_ndiag;
				collect_matr[2*&size.+1:3*&size., 2*&size.+1:3*&size.] = Re_diag;
				collect_matr[3*&size.+1:4*&size., 2*&size.+1:3*&size.] = Re_ndiag;
				/* third column */
				collect_matr[1:&size.,            3*&size.+1:4*&size.] = -Im_ndiag;
				collect_matr[&size.+1:2*&size.,   3*&size.+1:4*&size.] = Im_diag;
				collect_matr[2*&size.+1:3*&size., 3*&size.+1:4*&size.] = -Re_ndiag;
				collect_matr[3*&size.+1:4*&size., 3*&size.+1:4*&size.] = Re_diag;
		
				result_eigvals = eigval(collect_matr)[{%do i=1 %to &size.; %eval(4*&i.) %end;}] / sqrt(2); /* eigevalues get quadrupled --- we select only one from each four */
			end;

		else 
			do;
				print "Beta has to be 1, 2 or 4.";
				result_eigvals = .;
			end;

		if &normalize. then
			do;
				
				result_eigvals = result_eigvals / sqrt( &beta. * &size. );
			end;
		
		create &result_eigvals. from result_eigvals; 
					append from result_eigvals; 

	quit;
%mend;


%macro semi_circle_law(x, R=2.0);
	%local ret;
	%let ret = %sysevalf(
					%sysevalf( 
						2/( %sysevalf( %sysfunc(constant(PI)) * %sysevalf(&R.*&R.) ) ) )
				* %sysfunc(sqrt( %sysfunc(max( %sysevalf( %sysevalf(&R.*&R.) - %sysevalf(&x.*&x.)), 0 )) ))
				);
	&ret.
%mend;

%macro mu_ref_normal_sampler_tridiag(result_eigvals, size=10, loc=0.0, scale=1.0, beta=2, normalize=1, random_state=1618);

	proc iml;
		call streaminit(&random_state.);
		if &beta. <= 0 then
			do;
        		print "`beta` must be positive. Given: &beta.";
			end;
		else 
			do;
				M = diag( RandNormal(1, shape(&loc., &size., 1), I(&size.))); /* diagonal */
				do i = 1 to &size.-1;
					x = sqrt( &scale. * &scale. * rand("GAMMA", &beta. * (&size. - i)/2));
					M[i, i+1] = x;
					M[i+1, i] = x;
				end;

				result_eigvals = eigval(M);

				if &normalize. then
					do;
						result_eigvals = (result_eigvals - &loc.) / (sqrt(0.5 * &beta. * &size.) * &scale.) ;
					end; 

				create &result_eigvals. from result_eigvals; 
					append from result_eigvals; 				
			end;
	quit;
%mend;


/* 
	LAGUERRE
*/

options mprint; 
%macro marcenko_pasteur_law(x, M, N, sigma=1.0); /* M >= N */
	%local c Lm Lp res;
    %let c = %sysevalf( &N. / &M.);
	%let Lm = %sysevalf(%sysevalf(&sigma. * (%sysevalf(1 - %sysfunc(sqrt(&c.)) )) )**2);
	%let Lp = %sysevalf(%sysevalf(&sigma. * (%sysevalf(1 + %sysfunc(sqrt(&c.)) )) )**2);


	%let res = %sysevalf(
					%sysevalf(
						%sysfunc(sqrt( %sysfunc(max( (&Lp.-&x.)*(&x.-&Lm.),0 )) )) 
		 				/ %sysevalf(&c.*&x.)
					) 
					/ %sysevalf(
						%sysevalf( 2*%sysfunc(constant(PI)) )
						*%sysevalf(&sigma.**2)
						)
				);
	&res.
%mend;


%macro laguerre_sampler_full(result_eigvals, M, N, beta=2, normalize=1, random_state=1618); /* maybe change N to size */

	proc iml;
		call streaminit(&random_state.);

		result_eigvals = shape(0, &N., 1);
		if &beta. = 1 then
			do;
				A = RandNormal( &N., shape(0, &M., 1), I(&M.));
				result_eigvals = eigval(A * A`);
			end;

		else if &beta. = 2 then
			do;
				A = RandNormal( &N., shape(0, &M., 1), I(&M.));
				B = RandNormal( &N., shape(0, &M., 1), I(&M.));
				Re_mat = A * A` + B * B`;
				Im_mat = B * A` - A * B`;

				collect_matr = shape(0, 2 * &N., 2 * &N.);
				/* first column */
				collect_matr[1:&N.,      1:&N.] = Re_mat;
				collect_matr[&N.+1:2*&N.,1:&N.] = Im_mat;
				/* second column */
				collect_matr[1:&N.,      &N.+1:2*&N.] = - Im_mat;
				collect_matr[&N.+1:2*&N.,&N.+1:2*&N.] = Re_mat;

				result_eigvals = eigval(collect_matr)[{%do i=1 %to &N.; %eval(2*&i.) %end;}]; /* eigevalues get doubled --- we select only one from each pair */
			end;

		else if &beta. = 4 then
			do;
				A = RandNormal( &N., shape(0, &M., 1), I(&M.));
				B = RandNormal( &N., shape(0, &M., 1), I(&M.));
				C = RandNormal( &N., shape(0, &M., 1), I(&M.));
				D = RandNormal( &N., shape(0, &M., 1), I(&M.));

				Re_diag  = A*A` + B*B` + C*C` + D*D`;
				Re_ndiag = A*C` + D*B` - C*A` - B*D`;
				Im_diag  = B*A` - A*B` + D*C` - C*D`;
				Im_ndiag = D*A` + C*B` - B*C` - A*D`;

				collect_matr = shape(0, 4 * &N., 4 * &N.);

				/* first column */
				collect_matr[1:&N.,         1:&N.] = Re_diag;
				collect_matr[&N.+1:2*&N.,   1:&N.] = Re_ndiag;
				collect_matr[2*&N.+1:3*&N., 1:&N.] = Im_diag;
				collect_matr[3*&N.+1:4*&N., 1:&N.] = Im_ndiag;
				/* second column */
				collect_matr[1:&N.,         &N.+1:2*&N.] = - Re_ndiag;
				collect_matr[&N.+1:2*&N.,   &N.+1:2*&N.] = Re_diag;
				collect_matr[2*&N.+1:3*&N., &N.+1:2*&N.] = Im_ndiag;
				collect_matr[3*&N.+1:4*&N., &N.+1:2*&N.] = - Im_diag;
				/* third column */
				collect_matr[1:&N.,         2*&N.+1:3*&N.] = -Im_diag;
				collect_matr[&N.+1:2*&N.,   2*&N.+1:3*&N.] = -Im_ndiag;
				collect_matr[2*&N.+1:3*&N., 2*&N.+1:3*&N.] = Re_diag;
				collect_matr[3*&N.+1:4*&N., 2*&N.+1:3*&N.] = Re_ndiag;
				/* fourth column */
				collect_matr[1:&N.,         3*&N.+1:4*&N.] = -Im_ndiag;
				collect_matr[&N.+1:2*&N.,   3*&N.+1:4*&N.] = Im_diag;
				collect_matr[2*&N.+1:3*&N., 3*&N.+1:4*&N.] = -Re_ndiag;
				collect_matr[3*&N.+1:4*&N., 3*&N.+1:4*&N.] = Re_diag;
		
				result_eigvals = eigval(collect_matr)[{%do i=1 %to &N.; %eval(4*&i.) %end;}]; /* eigevalues get quadrupled --- we select only one from each four */
				
			end;

		else 
			do;
				print "Beta has to be 1, 2 or 4.";
				result_eigvals = .;
			end;

		if &normalize. then
			do;
				result_eigvals = result_eigvals / (&beta. * &M.); 
			end;

		create &result_eigvals. from result_eigvals; 
					append from result_eigvals; 
	quit;
%mend;


%macro mu_ref_gamma_sampler_tridiag(result_eigvals, shape=1.0, scale=1.0, beta=2, size=10, normalize=1, random_state=1618);

	proc iml;
		call streaminit(&random_state.);

		result_eigvals = shape(0, &size., 1);
		if &beta. <= 0 then
			do;
        		print "`beta` must be positive. Given: &beta.";
				result_eigvals = .;
			end;
		else 
			do;
				xi_odd = shape(0, &size., 1);
				xi_even = shape(0, &size., 1);

				do i = 1 to &size.-1;
					xi_odd[i]    = &scale. * rand("GAMMA", &beta. * (&size. - i)/2 + &shape.);
					xi_even[i+1] = &scale. * rand("GAMMA", &beta. * (&size. - i)/2);
				end;
				xi_odd[&size.]   = &scale. * rand("GAMMA", &shape.);

				M = diag( xi_odd + xi_even);

				do i = 1 to &size.-1;
					M[i+1,i] = sqrt(xi_odd[i] * xi_even[i+1]);
					M[i,i+1] = sqrt(xi_odd[i] * xi_even[i+1]);
				end;

				result_eigvals = eigval(M)[,1];
			
				if &normalize. then
					do; /* potencjalnie cos nie tak */
						result_eigvals = result_eigvals / (0.5 * &scale. * &beta. * (2 / &beta. * &shape. + &size. - 1));
					end;

				create &result_eigvals. from result_eigvals; 
					append from result_eigvals; 		
			end;
	quit;
%mend;


/*
	JACOBI
*/


%macro jacobi_sampler_full(result_eigvals, M_1, M_2, N, beta=2, normalize=1, random_state=1618); /* normalization does nothing here */
	
	proc iml;
		call streaminit(&random_state.);
		
		if &beta. = 1 then
			do;
				X = RandNormal( &N., shape(0, &M_1., 1), I(&M_1.));
				Y = RandNormal( &N., shape(0, &M_2., 1), I(&M_2.));

				X_tmp = X*X`;
				Y_tmp = Y*Y`;

				res_tmp = X_tmp * inv( X_tmp + Y_tmp );
				result_eigvals = eigval(res_tmp)[,1];
				create &result_eigvals. from result_eigvals; 
					append from result_eigvals; 
			end;
		else if &beta. = 2 then
			do;
				A_X = RandNormal( &N., shape(0, &M_1., 1), I(&M_1.));
				B_X = RandNormal( &N., shape(0, &M_1., 1), I(&M_1.));
				A_Y = RandNormal( &N., shape(0, &M_2., 1), I(&M_2.));
				B_Y = RandNormal( &N., shape(0, &M_2., 1), I(&M_2.));

				A_X_tmp = A_X*A_X` + B_X*B_X`;
				B_X_tmp = B_X*A_X` - A_X*B_X`;
				A_Y_tmp = A_Y*A_Y` + B_Y*B_Y`;
				B_Y_tmp = B_Y*A_Y` - A_Y*B_Y`;

				A_to_solve = shape(0, 2*&N., 2*&N.);
				B_to_solve = shape(0, 2*&N., &N.);

				/* first column */
				A_to_solve[1:&N.,             1:&N.] = A_X_tmp` + A_Y_tmp`;
				A_to_solve[&N.+1:2*&N.,       1:&N.] = B_X_tmp` + B_Y_tmp`;
				/* second column */
				A_to_solve[1:&N.,       &N.+1:2*&N.] = -(B_X_tmp` + B_Y_tmp`);
				A_to_solve[&N.+1:2*&N., &N.+1:2*&N.] = A_X_tmp` + A_Y_tmp`;

				B_to_solve[1:&N.,       1:&N.] = A_X_tmp`;
				B_to_solve[&N.+1:2*&N., 1:&N.] = B_X_tmp`;

				quot = solve(A_to_solve, B_to_solve)`;

				collect_mat = shape(0, 2*&N., 2*&N.);
				/* first column */
				collect_mat[1:&N.,             1:&N.] = quot[1:&N.,       1:&N.];
				collect_mat[&N.+1:2*&N.,       1:&N.] = quot[1:&N., &N.+1:2*&N.];
				/* second column */
				collect_mat[1:&N.,       &N.+1:2*&N.] = -quot[1:&N., &N.+1:2*&N.];
				collect_mat[&N.+1:2*&N., &N.+1:2*&N.] = quot[1:&N.,       1:&N.];

				result_eigvals = eigval(collect_mat)[{%do i=1 %to &N.; %eval(2*&i.) %end;},1];
				
				create &result_eigvals. from result_eigvals; 
					append from result_eigvals; 
			end;

		else if &beta. = 4 then
			do;
				A_X = RandNormal( &N., shape(0, &M_1., 1), I(&M_1.));
				B_X = RandNormal( &N., shape(0, &M_1., 1), I(&M_1.));
				C_X = RandNormal( &N., shape(0, &M_1., 1), I(&M_1.));
				D_X = RandNormal( &N., shape(0, &M_1., 1), I(&M_1.));
				A_Y = RandNormal( &N., shape(0, &M_2., 1), I(&M_2.));
				B_Y = RandNormal( &N., shape(0, &M_2., 1), I(&M_2.));
				C_Y = RandNormal( &N., shape(0, &M_2., 1), I(&M_2.));
				D_Y = RandNormal( &N., shape(0, &M_2., 1), I(&M_2.));

				/* Changes in sign are because of transposition */
				Re_diag_X_t  =  A_X*A_X` + B_X*B_X` + C_X*C_X` + D_X*D_X`;
				Re_ndiag_X_t = -A_X*C_X` - D_X*B_X` + C_X*A_X` + B_X*D_X`;
				Im_diag_X_t  = -B_X*A_X` + A_X*B_X` - D_X*C_X` + C_X*D_X`;
				Im_ndiag_X_t = -D_X*A_X` - C_X*B_X` + B_X*C_X` + A_X*D_X`;

				Re_diag_tot_t  = Re_diag_X_t  + A_Y*A_Y` + B_Y*B_Y` + C_Y*C_Y` + D_Y*D_Y`;
				Re_ndiag_tot_t = Re_ndiag_X_t - A_Y*C_Y` - D_Y*B_Y` + C_Y*A_Y` + B_Y*D_Y`;
				Im_diag_tot_t  = Im_diag_X_t  - B_Y*A_Y` + A_Y*B_Y` - D_Y*C_Y` + C_Y*D_Y`;
				Im_ndiag_tot_t = Im_ndiag_X_t - D_Y*A_Y` - C_Y*B_Y` + B_Y*C_Y` + A_Y*D_Y`;

				A_to_solve = shape(0, 4*&N., 4*&N.);
				B_to_solve = shape(0, 4*&N., 2*&N.);

				/* first column */
				A_to_solve[1:&N.,         1:&N.] = Re_diag_tot_t;
				A_to_solve[&N.+1:2*&N.,   1:&N.] = - Re_ndiag_tot_t;
				A_to_solve[2*&N.+1:3*&N., 1:&N.] = Im_diag_tot_t;
				A_to_solve[3*&N.+1:4*&N., 1:&N.] = Im_ndiag_tot_t;
				/* second column */
				A_to_solve[1:&N.,         &N.+1:2*&N.] =   Re_ndiag_tot_t;
				A_to_solve[&N.+1:2*&N.,   &N.+1:2*&N.] =   Re_diag_tot_t;
				A_to_solve[2*&N.+1:3*&N., &N.+1:2*&N.] =   Im_ndiag_tot_t;
				A_to_solve[3*&N.+1:4*&N., &N.+1:2*&N.] = - Im_diag_tot_t;
				/* third column */
				A_to_solve[1:&N.,         2*&N.+1:3*&N.] = - Im_diag_tot_t;
				A_to_solve[&N.+1:2*&N.,   2*&N.+1:3*&N.] = - Im_ndiag_tot_t;
				A_to_solve[2*&N.+1:3*&N., 2*&N.+1:3*&N.] =   Re_diag_tot_t;
				A_to_solve[3*&N.+1:4*&N., 2*&N.+1:3*&N.] = - Re_ndiag_tot_t;
				/* fourth column */
				A_to_solve[1:&N.,         3*&N.+1:4*&N.] = - Im_ndiag_tot_t;
				A_to_solve[&N.+1:2*&N.,   3*&N.+1:4*&N.] =   Im_diag_tot_t;
				A_to_solve[2*&N.+1:3*&N., 3*&N.+1:4*&N.] = - Re_ndiag_tot_t;
				A_to_solve[3*&N.+1:4*&N., 3*&N.+1:4*&N.] =   Re_diag_tot_t;
	
				/* first column */
				B_to_solve[1:&N.,         1:&N.] =   Re_diag_X_t;
				B_to_solve[&N.+1:2*&N.,   1:&N.] = - Re_ndiag_X_t;
				B_to_solve[2*&N.+1:3*&N., 1:&N.] =   Im_diag_X_t;
				B_to_solve[3*&N.+1:4*&N., 1:&N.] =   Im_ndiag_X_t;
				/* second column */
				B_to_solve[1:&N.,         &N.+1:2*&N.] =   Re_ndiag_X_t;
				B_to_solve[&N.+1:2*&N.,   &N.+1:2*&N.] =   Re_diag_X_t;
				B_to_solve[2*&N.+1:3*&N., &N.+1:2*&N.] =   Im_ndiag_X_t;
				B_to_solve[3*&N.+1:4*&N., &N.+1:2*&N.] = - Im_diag_X_t;

				quot = solve(A_to_solve, B_to_solve);

				collect_mat = shape(0, 4*&N., 4*&N.);
				/* first column */
				collect_mat[1:2*&N.,             1:2*&N.] = quot[1:2*&N.,       1:2*&N.];
				collect_mat[2*&N.+1:4*&N.,       1:2*&N.] = quot[2*&N.+1:4*&N., 1:2*&N.];
				/* second column */
				collect_mat[1:2*&N.,       2*&N.+1:4*&N.] = -quot[2*&N.+1:4*&N., 1:2*&N.];
				collect_mat[2*&N.+1:4*&N., 2*&N.+1:4*&N.] = quot[1:2*&N.,        1:2*&N.];


				result_eigvals = eigval(collect_mat)[{%do i=1 %to &N.; %eval(4*&i.) %end;},1];
				
				create &result_eigvals. from result_eigvals; 
					append from result_eigvals; 
			end;

		else 
			do;
				print "Beta has to be 1, 2 or 4.";
			end;

	quit;
%mend;


%macro wachter_law(x, M_1, M_2, N); /* M_1, M_2>=N */
    %local a b Lm Mp ret;
	%let a = %sysevalf( &M_1. / &N.);
	%let b = %sysevalf( &M_2. / &N.);

    %let Lm = %sysevalf(
		%sysevalf( 
			%sysfunc(sqrt(
				%sysevalf(&a. 
					* %sysevalf(&a. 
						+ %sysevalf(&b. - 1)
					)
				)
			)) 
		- %sysfunc(sqrt(&b.))
		) / %sysevalf(&a. + &b.))**2;

    %let Lp = %sysevalf(
		%sysevalf( 
			%sysfunc(sqrt(
				%sysevalf(&a. 
					* %sysevalf(&a. 
						+ %sysevalf(&b. - 1)
					)
				)
			)) 
		+ %sysfunc(sqrt(&b.))
		) / %sysevalf(&a. + &b.))**2;

    %let ret = %sysevalf(
					%sysevalf(
						%sysevalf(
							%sysevalf(&a.+&b.)
							/ %sysevalf(2 * %sysfunc(constant(PI)) ) 
						)
						* %sysevalf(1/ %sysevalf(&x.*%sysevalf(1-&x.) ) ) 
					)
			* %sysfunc(sqrt( %sysfunc(max( %sysevalf( %sysevalf(&Lp.-&x.)*%sysevalf(&x.-&Lm.) ) ,0 )) )) 
			);
	&ret.
%mend;


%macro mu_ref_beta_sampler_tridiag(result_eigvals, a, b, beta=2, size=10, normalize=1, random_state=1618); /* normalization does nothing here */
	
	proc iml;
		call streaminit(&random_state.);

			if &beta. <= 0 then
			do;
        		print "`beta` must be positive. Given: &beta.";
			end;
		else 
			do;
				c_odd  = shape(0, &size., 1);
				c_even = shape(0, &size., 1);
				do i=1 to &size.-1;
					c_odd[i]    = rand("BETA", &beta.*(&size. - i)/2 + &a., &beta.*(&size. - i)/2 + &b.);
					c_even[i+1] = rand("BETA", &beta.*(&size. - i)/2, &beta.*(&size. - i - 1)/2 + &b. + &a.);
				end;
				c_odd[&size.] = rand("BETA", &a., &b.);

				xi_odd = (1 - c_even) # c_odd;
				xi_even = shape(0, &size., 1);
				xi_even[2:&size.] = (1 - c_odd[1:&size.-1]) # c_even[2:&size.];

				alpha_coef = xi_even + xi_odd;
				beta_coef = sqrt(xi_odd[1:&size.-1] # xi_even[2:&size.]);
	
				M = diag( alpha_coef);
				do i=1 to &size.-1;
					M[i, i+1] = beta_coef[i];
					M[i+1, i] = beta_coef[i];
				end;

				result_eigvals = eigval(M)[,1];
				create &result_eigvals. from result_eigvals; 
					append from result_eigvals; 
			end;
	quit;
%mend;

/* 
	CIRCULAR
*/


%macro circular_sampler_full(result_eigvals, N, beta=2, haar_mode='QR', normalize=1, heurestic_fix=1, random_state=1618); /*normalization does nothing here */
    	
	proc iml;
		call streaminit(&random_state.);

		if "%upcase(&haar_mode.)" = "HERMITE" then
			do;

				if &beta. = 1 then
					do;
						A = RandNormal( &N., shape(0, &N., 1), I(&N.));
						eigvecs = eigvec( A + A` );
						result_eigvals = eigval( eigvecs);

						create &result_eigvals. from result_eigvals; 
							append from result_eigvals; 
					end;

				else if &beta. = 2 then
					do;
						A = RandNormal( &N., shape(0, &N., 1), I(&N.));
						B = RandNormal( &N., shape(0, &N., 1), I(&N.));

						Re_mat = A + A`;
						Im_mat = B - B`;

						collect_matr = shape(0, 2 * &N., 2 * &N.);
						/* first column */
						collect_matr[1:&N.,      1:&N.] = Re_mat;
						collect_matr[&N.+1:2*&N.,1:&N.] = Im_mat;
						/* second column */
						collect_matr[1:&N.,      &N.+1:2*&N.] = - Im_mat;
						collect_matr[&N.+1:2*&N.,&N.+1:2*&N.] = Re_mat;

						eigvecs = eigvec(collect_matr);
						result_eigvals = eigval(eigvecs); 

						if &heurestic_fix. then
							do;
								eigvals_fixed_heur = shape(0, &N., 2);
								probs_for_heurestic = shape(0, &N., 1);
								call randgen(probs_for_heurestic, "Uniform");	
								
								

								do i=1 to &N.;
									if probs_for_heurestic[i] < 0.5 then
										eigvals_fixed_heur[i,] = result_eigvals[2*i-1,];
									else
										eigvals_fixed_heur[i,] = result_eigvals[2*i,];
								end;

								create &result_eigvals. from eigvals_fixed_heur; 
									append from eigvals_fixed_heur;  
							end;
						else
							do;
								create &result_eigvals. from result_eigvals; 
									append from result_eigvals;  
							end;

						
					end;

				else if &beta. = 4 then
					do;
						A = RandNormal( &N., shape(0, &N., 1), I(&N.));
						B = RandNormal( &N., shape(0, &N., 1), I(&N.));
						C = RandNormal( &N., shape(0, &N., 1), I(&N.));
						D = RandNormal( &N., shape(0, &N., 1), I(&N.));

						Re_diag  =  A + A`;
						Re_ndiag = -C + C;
						Im_diag  =  B - B`;
						Im_ndiag =  D - D`;

						collect_matr = shape(0, 4 * &N., 4 * &N.);

						/* first column */
						collect_matr[1:&N.,         1:&N.] = Re_diag;
						collect_matr[&N.+1:2*&N.,   1:&N.] = Re_ndiag;
						collect_matr[2*&N.+1:3*&N., 1:&N.] = Im_diag;
						collect_matr[3*&N.+1:4*&N., 1:&N.] = Im_ndiag;
						/* second column */
						collect_matr[1:&N.,         &N.+1:2*&N.] = - Re_ndiag;
						collect_matr[&N.+1:2*&N.,   &N.+1:2*&N.] = Re_diag;
						collect_matr[2*&N.+1:3*&N., &N.+1:2*&N.] = Im_ndiag;
						collect_matr[3*&N.+1:4*&N., &N.+1:2*&N.] = - Im_diag;
						/* third column */
						collect_matr[1:&N.,         2*&N.+1:3*&N.] = -Im_diag;
						collect_matr[&N.+1:2*&N.,   2*&N.+1:3*&N.] = -Im_ndiag;
						collect_matr[2*&N.+1:3*&N., 2*&N.+1:3*&N.] = Re_diag;
						collect_matr[3*&N.+1:4*&N., 2*&N.+1:3*&N.] = Re_ndiag;
						/* third column */
						collect_matr[1:&N.,         3*&N.+1:4*&N.] = -Im_ndiag;
						collect_matr[&N.+1:2*&N.,   3*&N.+1:4*&N.] = Im_diag;
						collect_matr[2*&N.+1:3*&N., 3*&N.+1:4*&N.] = -Re_ndiag;
						collect_matr[3*&N.+1:4*&N., 3*&N.+1:4*&N.] = Re_diag;


						eigvecs = eigvec(collect_matr);	
						result_eigvals = eigval(eigvecs); 

						if &heurestic_fix. then
							do;
								eigvals_fixed_heur = shape(0, &N., 2);
								probs_for_heurestic = shape(0, &N., 1);
								call randgen(probs_for_heurestic, "Uniform");	
								
								

								do i=1 to &N.;
									if probs_for_heurestic[i] < 0.25 then
										eigvals_fixed_heur[i,] = result_eigvals[2*i-1,];

									else if probs_for_heurestic[i] < 0.5 then
										eigvals_fixed_heur[i,] = result_eigvals[2*i,];

									else if probs_for_heurestic[i] < 0.75 then
										eigvals_fixed_heur[i,] = result_eigvals[4*&N. - 2*i+1,];

									else
										eigvals_fixed_heur[i,] = result_eigvals[4*&N. - 2*i+2,];
								end;

								create &result_eigvals. from eigvals_fixed_heur; 
									append from eigvals_fixed_heur;  
							end;
						else
							do;
								create &result_eigvals. from result_eigvals; 
									append from result_eigvals;  
							end;
					end;

				else 
					do;
						print "For Hermite version Beta has to be 1, 2 or 4.";
					end;
			end;

		else if "%upcase(&haar_mode.)" = "QR" then
			do;
				if &beta. = 1 then
					do;
						A = RandNormal( &N., shape(0, &N., 1), I(&N.));
						call qr(Q, R, piv, indep, A);
						D = vecdiag(R); /* should there be a transposition here? */
      					U = Q # sign(D); 
						result_eigvals = eigval(U); 
						
						create &result_eigvals. from result_eigvals; 
							append from result_eigvals; 
					end;

				else 
					do;
						print "`QR` method not provided for `beta` ither than 1";
					end;
			end;

		else 
			do;
				print "Wrong haar_mode supplied --- only `Hermite` and `QR` are available";
			end;
	quit;
%mend;


%macro mu_ref_unif_unit_circle_sampler(result_eigvals, beta=2, size=10, normalize=1, heurestic_fix=1, random_state=1618); /* normalization does nothing here */

	proc iml;
		call streaminit(&random_state.);
		if &beta. = int(&beta.) && sign(&beta.) = 1 then
			do;
				alpha_re = shape(0, &size., 1);
				alpha_im = shape(0, &size., 1);
				do i=1 to &size.;
					nu = 1 + &beta. * (&size. - i);
					gauss_vec = RandNormal(1, shape(0, nu+1, 1), I(nu+1));
					alpha_re[i] = gauss_vec[1] / norm(gauss_vec);
					alpha_im[i] = gauss_vec[2] / norm(gauss_vec);
				end;

				rho = sqrt( 1 - alpha_re[1:&size.-1]##2 - alpha_im[1:&size.-1]##2 );

				L = shape(0, 2*(&size.), 2*(&size.));
				M = shape(0, 2*(&size.), 2*(&size.));
				
				M[1,1] = 1;
				M[&size.+1,&size.+1] = 1;
				if mod(&size., 2) = 0 then 
					do;
						do i=1 to int(&size./2);
							/* L */

							/* over diagonal */
							L[2*i-1, 2*i] = rho[2*i-1];
							L[2*i-1 + &size., 2*i + &size.] = rho[2*i-1];
							/* diagonal */
							L[2*i-1, 2*i-1] = alpha_re[2*i-1];
							L[2*i, 2*i] = -alpha_re[2*i-1];
							L[2*i-1 + &size., 2*i-1 + &size.] = alpha_re[2*i-1];
							L[2*i + &size., 2*i + &size.] = -alpha_re[2*i-1];

							L[2*i-1 + &size., 2*i-1] = -alpha_im[2*i-1];
							L[2*i + &size., 2*i] = -alpha_im[2*i-1];
							L[2*i-1 , 2*i-1 + &size.] = alpha_im[2*i-1];
							L[2*i, 2*i + &size.] = alpha_im[2*i-1];
							/* under diagonal */
							L[2*i, 2*i-1] = rho[2*i-1];
							L[2*i + &size., 2*i-1 + &size.] = rho[2*i-1];
						end;
						
						
						do i=1 to int(&size./2)-1;
							/* M */

							/* over diagonal */
							M[2*i, 2*i+1] = rho[2*i];
							M[2*i + &size., 2*i + &size.+1] = rho[2*i];
							/* diagonal */
							M[2*i, 2*i] = alpha_re[2*i];
							M[2*i+1, 2*i+1] = -alpha_re[2*i];
							M[2*i + &size., 2*i + &size.] = alpha_re[2*i];
							M[2*i+1 + &size., 2*i+1 + &size.] = -alpha_re[2*i];

							M[2*i+ &size., 2*i] = -alpha_im[2*i];
							M[2*i+1+ &size., 2*i+1] = -alpha_im[2*i];
							M[2*i, 2*i + &size.] = alpha_im[2*i];
							M[2*i+1, 2*i+1 + &size.] = alpha_im[2*i];
							/* under diagonal */
							M[2*i+1, 2*i] = rho[2*i];
							M[2*i+1 + &size., 2*i + &size.] = rho[2*i];
						end;
						M[&size.,&size.] = alpha_re[&size.];
						M[2*(&size.), 2*(&size.)] = alpha_re[&size.];
						M[&size., 2*(&size.)] = alpha_im[&size.];
						M[2*(&size.), &size.] = -alpha_im[&size.];
					end;
				else
					do;
						/* L */
						L[&size., &size.] = alpha_re[&size.];
						L[2*(&size.), 2*(&size.)] = alpha_re[&size.];
						L[&size., 2*(&size.)] = alpha_im[&size.];
						L[2*(&size.), &size.] = -alpha_im[&size.];

						do i=1 to int(&size./2);
							/* over diagonal */
							L[2*i-1, 2*i] = rho[2*i-1];
							L[2*i-1 + &size., 2*i + &size.] = rho[2*i-1];
							/* diagonal */
							L[2*i-1, 2*i-1] = alpha_re[2*i-1];
							L[2*i, 2*i] = -alpha_re[2*i-1];
							L[2*i-1 + &size., 2*i-1 + &size.] = alpha_re[2*i-1];
							L[2*i + &size., 2*i + &size.] = -alpha_re[2*i-1];

							L[2*i-1 + &size., 2*i-1] = -alpha_im[2*i-1];
							L[2*i + &size., 2*i] = -alpha_im[2*i-1];
							L[2*i-1 , 2*i-1 + &size.] = alpha_im[2*i-1];
							L[2*i, 2*i + &size.] = alpha_im[2*i-1];
							/* under diagonal */
							L[2*i, 2*i-1] = rho[2*i-1];
							L[2*i + &size., 2*i-1 + &size.] = rho[2*i-1];

							/* M */

							/* over diagonal */
							M[2*i, 2*i+1] = rho[2*i];
							M[2*i + &size., 2*i + &size.+1] = rho[2*i];
							/* diagonal */
							M[2*i, 2*i] = alpha_re[2*i];
							M[2*i+1, 2*i+1] = -alpha_re[2*i];
							M[2*i + &size., 2*i + &size.] = alpha_re[2*i];
							M[2*i+1 + &size., 2*i + &size.+1] = -alpha_re[2*i];

							M[2*i + &size., 2*i] = -alpha_im[2*i];
							M[2*i+1 + &size., 2*i+1] = -alpha_im[2*i];
							M[2*i , 2*i + &size.] = alpha_im[2*i];
							M[2*i+1, 2*i + &size.+1] = alpha_im[2*i];
							/* under diagonal */
							M[2*i+1, 2*i] = rho[2*i];
							M[2*i+1 + &size., 2*i + &size.] = rho[2*i];
						end;
					end;

				result_eigvals = eigval(L*M); 

				if &heurestic_fix. then
							do;
								eigvals_fixed_heur = shape(0, &size., 2);
								probs_for_heurestic = shape(0, &size., 1);
								call randgen(probs_for_heurestic, "Uniform");	
								
								

								do i=1 to &size.;
									if probs_for_heurestic[i] < 0.5 then
										eigvals_fixed_heur[i,] = result_eigvals[2*i-1,];
									else
										eigvals_fixed_heur[i,] = result_eigvals[2*i,];
								end;

								create &result_eigvals. from eigvals_fixed_heur; 
									append from eigvals_fixed_heur;  
							end;
				else
							do;
								create &result_eigvals. from result_eigvals; 
									append from result_eigvals;  
							end;  
			end;
		else
			do;
				print "`beta` must be positive integer. Given: &beta.";
			end;	
	quit;
%mend;



/*
	GINIBRE
*/


%macro ginibre_sampler_full(result_eigvals, N, normalize=1, heurestic_fix=1, random_state=1618);
  
	proc iml;
		call streaminit(&random_state.);
		A = RandNormal( &N., shape(0, &N., 1), I(&N.));
		B = RandNormal( &N., shape(0, &N., 1), I(&N.));

    	collect_matrix = shape(0, 2*&N., 2*&N.);
		collect_matrix[1:&N., 1:&N.] = A;
		collect_matrix[&N.+1:2*&N., 1:&N.] = B;
		collect_matrix[1:&N.,&N.+1:2*&N.] = -B;
		collect_matrix[&N.+1:2*&N.,&N.+1:2*&N.] = A;

		result_eigvals = eigval(collect_matrix)/sqrt(2);
		if &normalize. then
			do;
				result_eigvals = result_eigvals / sqrt(&N.);
			end;

		if &heurestic_fix. then
			do;
				eigvals_fixed_heur = shape(0, &N., 2);
				probs_for_heurestic = shape(0, &N., 1);
				call randgen(probs_for_heurestic, "Uniform");	
												
				do i=1 to &N.;
					if probs_for_heurestic[i] < 0.5 then
						eigvals_fixed_heur[i,] = result_eigvals[2*i-1,];
					else
						eigvals_fixed_heur[i,] = result_eigvals[2*i,];
				end;

				create &result_eigvals. from eigvals_fixed_heur; 
					append from eigvals_fixed_heur;  
			end;
		else
			do;
				create &result_eigvals. from result_eigvals; 
					append from result_eigvals; 
			end;
		
	quit;
%mend;


%macro sample_from_beta_ensemble_full(
	result_eigvals, 
	ensemble_version, 
	M_1, M_2, /* M variables only for Laguerre (first) and Jacobi (both) */
	size=10, 
	beta=2, 
	normalize=1, 
	haar_mode="Hermite", /* haar_mode only available for circular ensemble */
	heurestic_fix=1, /* heurestic_fix only available for circular ensemble */
	random_state=1618); 
	
	%if %upcase(&ensemble_version.) = HERMITE %then
		%hermite_sampler_full(
			result_eigvals=&result_eigvals.,
			size=&size., 
			beta=&beta., 
			normalize=&normalize.,
			random_state=&random_state.);

	%else %if %upcase(&ensemble_version.) = LAGUERRE %then
		%laguerre_sampler_full(
			result_eigvals=&result_eigvals.,
			M=&M_1., 
			N=&size., 
			beta=&beta., 
			normalize=&normalize., 
			random_state=&random_state.); 

	%else %if %upcase(&ensemble_version.) = JACOBI %then
		%jacobi_sampler_full(
			result_eigvals=&result_eigvals., 
			M_1=&M_1., 
			M_2=&M_2., 
			N=&size., 
			beta=&beta., 
			normalize=&normalize., 
			random_state=&random_state.); 

	%else %if %upcase(&ensemble_version.) = CIRCULAR %then
		%circular_sampler_full(
			result_eigvals=&result_eigvals., 
			N=&size., 
			beta=&beta., 
			haar_mode=&haar_mode., 
			normalize=&normalize., 
			heurestic_fix=&heurestic_fix., 
			random_state=&random_state.); /*normalization does nothing here */
    
	%else %if %upcase(&ensemble_version.) = GINIBRE %then
		%ginibre_sampler_full(
			result_eigvals=&result_eigvals., 
			N=&size., 
			normalize=&normalize., 
			heurestic_fix=&heurestic_fix., 
			random_state=&random_state.);

	%else
		%put Wrong type supplied - only Hermite, Laguerre, Jacobi, and Circular are available.;

%mend;


%macro sample_from_beta_ensemble_banded(
	result_eigvals, 
	ensemble_version, 
	size=10, 
	beta=2, 
	/* distribution specific parameters */
	loc=0.0, 
	scale=1.0, 
	shape = 1.0,
	a = 1.0,
	b = 1.0,
	/* parameters for the algorithm */
	normalize=1, 
	heurestic_fix=1,
	random_state=1618);
	
	%if %upcase(&ensemble_version.) = HERMITE %then
		%do;
			%mu_ref_normal_sampler_tridiag(
				result_eigvals = &result_eigvals., 
				size=&size., 
				loc=&loc., 
				scale=&scale., 
				beta=&beta., 
				normalize=&normalize.,
				random_state=&random_state.);
		%end;
		
	%else %if %upcase(&ensemble_version.) = LAGUERRE %then
		%do;
			%mu_ref_gamma_sampler_tridiag(
				result_eigvals=&result_eigvals., 
				shape=&shape.,  
				scale=&scale., 
				size=&size., 
				beta=&beta., 
				normalize=&normalize.,
				random_state=&random_state.);
		%end;

	%else %if %upcase(&ensemble_version.) = JACOBI %then
		%do;
			%mu_ref_beta_sampler_tridiag(
				result_eigvals=&result_eigvals., 
				a=&a., 
				b=&b., 
				beta=&beta.,
				size=&size., 
				normalize=&normalize., 
				random_state=&random_state.); /* normalization does nothing here */
		%end;
	
	%else %if %upcase(&ensemble_version.) = CIRCULAR %then
		%do;
			%mu_ref_unif_unit_circle_sampler(
				result_eigvals=&result_eigvals., 
				beta=&beta., 
				size=&size., 
				normalize=&normalize., 
				heurestic_fix=&heurestic_fix.,
				random_state=&random_state.); /* normalization does nothing here */
		%end;

	%else
		%do;
			%put Wrong type supplied - only Hermite, Laguerre, Jacobi, and Circular are available.;
		%end;

%mend;


