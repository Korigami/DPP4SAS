

proc iml;

start _mu_ref_normal_sampler_tridiag(size=10, loc=0.0, scale=1.0, beta=2, normalize=1, random_state=1618);

	call streaminit(random_state);
		if beta <= 0 then
			do;
				result_eigvals = .;
        		print "`beta` must be positive.";
			end;
		else 
			do;
				M = diag( RandNormal(1, shape(loc, size, 1), I(size))); /* diagonal */
				do i = 1 to size-1;
					x = sqrt( scale * scale * rand("GAMMA", beta * (size - i)/2));
					M[i, i+1] = x;
					M[i+1, i] = x;
				end;

				result_eigvals = eigval(M);

				if normalize then
					do;
						result_eigvals = (result_eigvals - loc) / (sqrt(0.5 * beta * size) * scale) ;
					end;  		
			end;
	return result_eigvals;

finish;

	start _hermite_sampler_full(size=10, beta=2, normalize=1, random_state=1618);

		call streaminit(random_state);
		result_eigvals = shape(0, size, 1);

		if beta = 1 then
			do;
				A = RandNormal( size, shape(0, size, 1), I(size));
				result_eigvals = eigval( A + A` )/sqrt(2) ;
			end;

		else if beta = 2 then
			do;
				A = RandNormal(size, shape(0, size, 1), I(size));
				B = RandNormal(size, shape(0, size, 1), I(size));

				Re_mat = A + A`;
				Im_mat = B - B`;

				collect_matr = shape(0, 2 * size, 2 * size);

				ind1 = 1:size;
				ind2 = (size+1):(2*size);
				/* first column */
				collect_matr[ind1, ind1] = Re_mat;
				collect_matr[ind2, ind1] = Im_mat;
				/* second column */
				collect_matr[ind1, ind2] = - Im_mat;
				collect_matr[ind2, ind2] = Re_mat;

				chosen_indices = do(2, 2*size, 2);
				result_eigvals = eigval(collect_matr)[chosen_indices] / sqrt(2); /* eigevalues get doubled --- we select only one from each pair */ 
			end;

		else if beta = 4 then
			do;
				A = RandNormal(size, shape(0, size, 1), I(size));
				B = RandNormal(size, shape(0, size, 1), I(size));
				C = RandNormal(size, shape(0, size, 1), I(size));
				D = RandNormal(size, shape(0, size, 1), I(size));

				Re_diag  =  A + A`;
				Re_ndiag = -C + C;
				Im_diag  =  B - B`;
				Im_ndiag =  D - D`;

				collect_matr = shape(0, 4 * size, 4 * size);

				ind1 = 1:size;
				ind2 = (size+1):(2*size);
				ind3 = (2*size+1):(3*size);
				ind4 = (3*size+1):(4*size);
				/* first column */
				collect_matr[ind1, ind1] = Re_diag;
				collect_matr[ind2, ind1] = Re_ndiag;
				collect_matr[ind3, ind1] = Im_diag;
				collect_matr[ind4, ind1] = Im_ndiag;
				/* second column */
				collect_matr[ind1, ind2] = - Re_ndiag;
				collect_matr[ind2, ind2] = Re_diag;
				collect_matr[ind3, ind2] = Im_ndiag;
				collect_matr[ind4, ind2] = - Im_diag;
				/* third column */
				collect_matr[ind1, ind3] = -Im_diag;
				collect_matr[ind2, ind3] = -Im_ndiag;
				collect_matr[ind3, ind3] = Re_diag;
				collect_matr[ind4, ind3] = Re_ndiag;
				/* third column */
				collect_matr[ind1, ind4] = -Im_ndiag;
				collect_matr[ind2, ind4] = Im_diag;
				collect_matr[ind3, ind4] = -Re_ndiag;
				collect_matr[ind4, ind4] = Re_diag;

				chosen_indices = do(4, 4*size, 4);
				result_eigvals = eigval(collect_matr)[chosen_indices] / sqrt(2); /* eigevalues get quadrupled --- we select only one from each four */
			end;

		else 
			do;
				print "Beta has to be 1, 2 or 4.";
				result_eigvals = .;
			end;

		if normalize then
			do;
				
				result_eigvals = result_eigvals / sqrt( beta * size );
			end;
		
		return result_eigvals;
	finish;

	start _laguerre_sampler_full(M, N=10, beta=2, normalize=1, random_state=1618); /* maybe change N to size */

		call streaminit(random_state);

		result_eigvals = shape(0, N, 1);
		if beta = 1 then
			do;
				A = RandNormal( N, shape(0, M, 1), I(M));
				result_eigvals = eigval(A * A`);
			end;

		else if beta = 2 then
			do;
				A = RandNormal( N, shape(0, M, 1), I(M));
				B = RandNormal( N, shape(0, M, 1), I(M));
				Re_mat = A * A` + B * B`;
				Im_mat = B * A` - A * B`;

				collect_matr = shape(0, 2 * N, 2 * N);

				ind1 = 1:N;
				ind2 = (N+1):(2*N);
				/* first column */
				collect_matr[ind1, ind1] = Re_mat;
				collect_matr[ind2, ind1] = Im_mat;
				/* second column */
				collect_matr[ind1, ind2] = - Im_mat;
				collect_matr[ind2, ind2] = Re_mat;

				chosen_indices = do(2, 2*N, 2);
				result_eigvals = eigval(collect_matr)[chosen_indices]; /* eigevalues get doubled --- we select only one from each pair */
			end;

		else if beta = 4 then
			do;
				A = RandNormal( N, shape(0, M, 1), I(M));
				B = RandNormal( N, shape(0, M, 1), I(M));
				C = RandNormal( N, shape(0, M, 1), I(M));
				D = RandNormal( N, shape(0, M, 1), I(M));

				Re_diag  = A*A` + B*B` + C*C` + D*D`;
				Re_ndiag = A*C` + D*B` - C*A` - B*D`;
				Im_diag  = B*A` - A*B` + D*C` - C*D`;
				Im_ndiag = D*A` + C*B` - B*C` - A*D`;

				collect_matr = shape(0, 4 * N, 4 * N);

				ind1 = 1:N;
				ind2 = (N+1):(2*N);
				ind3 = (2*N+1):(3*N);
				ind4 = (3*N+1):(4*N);
				/* first column */
				collect_matr[ind1, ind1] = Re_diag;
				collect_matr[ind2, ind1] = Re_ndiag;
				collect_matr[ind3, ind1] = Im_diag;
				collect_matr[ind4, ind1] = Im_ndiag;
				/* second column */
				collect_matr[ind1, ind2] = - Re_ndiag;
				collect_matr[ind2, ind2] = Re_diag;
				collect_matr[ind3, ind2] = Im_ndiag;
				collect_matr[ind4, ind2] = - Im_diag;
				/* third column */
				collect_matr[ind1, ind3] = -Im_diag;
				collect_matr[ind2, ind3] = -Im_ndiag;
				collect_matr[ind3, ind3] = Re_diag;
				collect_matr[ind4, ind3] = Re_ndiag;
				/* fourth column */
				collect_matr[ind1, ind4] = -Im_ndiag;
				collect_matr[ind2, ind4] = Im_diag;
				collect_matr[ind3, ind4] = -Re_ndiag;
				collect_matr[ind4, ind4] = Re_diag;
				
				chosen_indices = do(4, 4*N, 4);		
				result_eigvals = eigval(collect_matr)[chosen_indices]; /* eigevalues get quadrupled --- we select only one from each four */
				
			end;

		else 
			do;
				print "Beta has to be 1, 2 or 4.";
				result_eigvals = .;
			end;

		if normalize then
			do;
				result_eigvals = result_eigvals / (beta * M); 
			end;

		return result_eigvals; 
	finish;

	start _mu_ref_gamma_sampler_tridiag(shape=1.0, scale=1.0, beta=2, size=10, normalize=1, random_state=1618);

		call streaminit(random_state);

		result_eigvals = shape(0, size, 1);
		if beta <= 0 then
			do;
        		print "`beta` must be positive.";
				result_eigvals = .;
			end;
		else 
			do;
				xi_odd  = shape(0, size, 1);
				xi_even = shape(0, size, 1);

				do i = 1 to size-1;
					xi_odd[i]    = scale * rand("GAMMA", beta * (size - i)/2 + shape);
					xi_even[i+1] = scale * rand("GAMMA", beta * (size - i)/2);
				end;
				xi_odd[size]   = scale * rand("GAMMA", shape);

				M = diag( xi_odd + xi_even);

				do i = 1 to size-1;
					M[i+1,i] = sqrt(xi_odd[i] * xi_even[i+1]);
					M[i,i+1] = sqrt(xi_odd[i] * xi_even[i+1]);
				end;

				result_eigvals = eigval(M)[,1];
			
				if normalize then
					do; /* potencjalnie cos nie tak */
						result_eigvals = result_eigvals / (0.5 *scale * beta * (2 / beta * shape + size - 1));
					end;	
			end;

		return result_eigvals; 	
	finish;

	start _jacobi_sampler_full(M_1, M_2, N=10, beta=2, normalize=1, random_state=1618); /* normalization does nothing here */
	
		call streaminit(random_state);

		
		if beta = 1 then
			do;
				X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				Y = RandNormal( N, shape(0, M_2, 1), I(M_2));

				X_tmp = X*X`;
				Y_tmp = Y*Y`;

				res_tmp = X_tmp * inv( X_tmp + Y_tmp );
				result_eigvals = eigval(res_tmp)[,1];

				return result_eigvals; 
			end;

		else if beta = 2 then
			do;
				A_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				B_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				A_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));
				B_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));

				A_X_tmp = A_X*A_X` + B_X*B_X`;
				B_X_tmp = B_X*A_X` - A_X*B_X`;
				A_Y_tmp = A_Y*A_Y` + B_Y*B_Y`;
				B_Y_tmp = B_Y*A_Y` - A_Y*B_Y`;

				A_to_solve = shape(0, 2*N, 2*N);
				B_to_solve = shape(0, 2*N, N);


				ind1 = 1:N;
				ind2 = (N+1):(2*N);
				/* first column */
				A_to_solve[ind1, ind1] = A_X_tmp` + A_Y_tmp`;
				A_to_solve[ind2, ind1] = B_X_tmp` + B_Y_tmp`;
				/* second column */
				A_to_solve[ind1, ind2] = -(B_X_tmp` + B_Y_tmp`);
				A_to_solve[ind2, ind2] = A_X_tmp` + A_Y_tmp`;

				B_to_solve[ind1, ind1] = A_X_tmp`;
				B_to_solve[ind2, ind1] = B_X_tmp`;

				quot = solve(A_to_solve, B_to_solve)`;

				collect_mat = shape(0, 2*N, 2*N);
				/* first column */
				collect_mat[ind1, ind1] = quot[ind1, ind1];
				collect_mat[ind2, ind1] = quot[ind1, ind2];
				/* second column */
				collect_mat[ind1, ind2] = -quot[ind1, ind2];
				collect_mat[ind2, ind2] =  quot[ind1, ind1];

				chosen_indices = do(2, 2*N, 2);
				result_eigvals = eigval(collect_mat)[chosen_indices, 1];
				
				return result_eigvals; 
			end;

		else if beta = 4 then
			do;
				A_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				B_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				C_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				D_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				A_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));
				B_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));
				C_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));
				D_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));

				/* Changes in sign are because of transposition */
				Re_diag_X_t  =  A_X*A_X` + B_X*B_X` + C_X*C_X` + D_X*D_X`;
				Re_ndiag_X_t = -A_X*C_X` - D_X*B_X` + C_X*A_X` + B_X*D_X`;
				Im_diag_X_t  = -B_X*A_X` + A_X*B_X` - D_X*C_X` + C_X*D_X`;
				Im_ndiag_X_t = -D_X*A_X` - C_X*B_X` + B_X*C_X` + A_X*D_X`;

				Re_diag_tot_t  = Re_diag_X_t  + A_Y*A_Y` + B_Y*B_Y` + C_Y*C_Y` + D_Y*D_Y`;
				Re_ndiag_tot_t = Re_ndiag_X_t - A_Y*C_Y` - D_Y*B_Y` + C_Y*A_Y` + B_Y*D_Y`;
				Im_diag_tot_t  = Im_diag_X_t  - B_Y*A_Y` + A_Y*B_Y` - D_Y*C_Y` + C_Y*D_Y`;
				Im_ndiag_tot_t = Im_ndiag_X_t - D_Y*A_Y` - C_Y*B_Y` + B_Y*C_Y` + A_Y*D_Y`;

				A_to_solve = shape(0, 4*N, 4*N);
				B_to_solve = shape(0, 4*N, 2*N);

				ind1 = 1:N; 
				ind2 = (N+1):(2*N);
 				ind3 = (2*N+1):(3*N);
				ind4 = (3*N+1):(4*N);
				/* first column */
				A_to_solve[ind1, ind1] = Re_diag_tot_t;
				A_to_solve[ind2, ind1] = - Re_ndiag_tot_t;
				A_to_solve[ind3, ind1] = Im_diag_tot_t;
				A_to_solve[ind4, ind1] = Im_ndiag_tot_t;
				/* second column */
				A_to_solve[ind1, ind2] =   Re_ndiag_tot_t;
				A_to_solve[ind2, ind2] =   Re_diag_tot_t;
				A_to_solve[ind3, ind2] =   Im_ndiag_tot_t;
				A_to_solve[ind4, ind2] = - Im_diag_tot_t;
				/* third column */
				A_to_solve[ind1, ind3] = - Im_diag_tot_t;
				A_to_solve[ind2, ind3] = - Im_ndiag_tot_t;
				A_to_solve[ind3, ind3] =   Re_diag_tot_t;
				A_to_solve[ind4, ind3] = - Re_ndiag_tot_t;
				/* fourth column */
				A_to_solve[ind1, ind4] = - Im_ndiag_tot_t;
				A_to_solve[ind2, ind4] =   Im_diag_tot_t;
				A_to_solve[ind3, ind4] = - Re_ndiag_tot_t;
				A_to_solve[ind4, ind4] =   Re_diag_tot_t;
	
				/* first column */
				B_to_solve[ind1, ind1] =   Re_diag_X_t;
				B_to_solve[ind2, ind1] = - Re_ndiag_X_t;
				B_to_solve[ind3, ind1] =   Im_diag_X_t;
				B_to_solve[ind4, ind1] =   Im_ndiag_X_t;
				/* second column */
				B_to_solve[ind1, ind2] =   Re_ndiag_X_t;
				B_to_solve[ind2, ind2] =   Re_diag_X_t;
				B_to_solve[ind3, ind2] =   Im_ndiag_X_t;
				B_to_solve[ind4, ind2] = - Im_diag_X_t;

				quot = solve(A_to_solve, B_to_solve);

				collect_mat = shape(0, 4*N, 4*N);
				Lind1 = 1:(2*N);
				Lind2 = (2*N+1):(4*N);


				/* first column */
				collect_mat[Lind1, Lind1] = quot[Lind1, Lind1];
				collect_mat[Lind2, Lind1] = quot[Lind2, Lind1];
				/* second column */
				collect_mat[Lind1, Lind2] = -quot[Lind2, Lind1];
				collect_mat[Lind2, Lind2] = quot[Lind1, Lind1];

				
				chosen_indices = do(4, 4*N, 4);
				result_eigvals = eigval(collect_mat)[chosen_indices, 1];
				
				return result_eigvals; 
			end;

		else 
			do;
				print "Beta has to be 1, 2 or 4.";
				result_eigvals = .;
				return result_eigvals;
			end;

	finish;

	start _mu_ref_beta_sampler_tridiag(a, b, beta=2, size=10, normalize=1, random_state=1618); /* normalization does nothing here */
	
	
		call streaminit(random_state);

		if beta <= 0 then
			do;
        		print "`beta` must be positive.";
				result_eigvals = .;
				return result_eigvals;
			end;
		else 
			do;
				c_odd  = shape(0, size, 1);
				c_even = shape(0, size, 1);
				do i=1 to size-1;
					c_odd[i]    = rand("BETA", beta*(size - i)/2 + a, beta*(size - i)/2 + b);
					c_even[i+1] = rand("BETA", beta*(size - i)/2, beta*(size - i - 1)/2 + b + a);
				end;
				c_odd[size] = rand("BETA", a, b);

				xi_odd = (1 - c_even) # c_odd;
				xi_even = shape(0, size, 1);
				xi_even[2:size] = (1 - c_odd[1:size-1]) # c_even[2:size];

				alpha_coef = xi_even + xi_odd;
				beta_coef = sqrt(xi_odd[1:size-1] # xi_even[2:size]);
	
				M = diag( alpha_coef);
				do i=1 to size-1;
					M[i, i+1] = beta_coef[i];
					M[i+1, i] = beta_coef[i];
				end;

				result_eigvals = eigval(M)[,1];
				return result_eigvals;
			end;
	finish;


	start _circular_sampler_full(N, beta=2, haar_mode='QR', normalize=1, heurestic_fix=1, random_state=1618); /*normalization does nothing here */
    	
		call streaminit(random_state);

		if haar_mode = "Hermite" then
			do;

				if beta = 1 then
					do;
						A = RandNormal( N, shape(0, N, 1), I(N));
						eigvecs = eigvec( A + A` );
						result_eigvals = eigval( eigvecs);
						return result_eigvals;
					end;

				else if beta = 2 then
					do;
						A = RandNormal( N, shape(0, N, 1), I(N));
						B = RandNormal( N, shape(0, N, 1), I(N));

						Re_mat = A + A`;
						Im_mat = B - B`;

						collect_matr = shape(0, 2 * N, 2 * N);
						ind1 = 1:N;
						ind2 = (N+1):(2*N);
						/* first column */
						collect_matr[ind1, ind1] = Re_mat;
						collect_matr[ind2, ind1] = Im_mat;
						/* second column */
						collect_matr[ind1, ind2] = - Im_mat;
						collect_matr[ind2, ind2] = Re_mat;

						eigvecs = eigvec(collect_matr);
						result_eigvals = eigval(eigvecs); 


						if heurestic_fix then
							do;
								eigvals_fixed_heur = shape(0, N, 2);
								probs_for_heurestic = shape(0, N, 1);
								call randgen(probs_for_heurestic, "Uniform");	

								do i=1 to N;
									if probs_for_heurestic[i] < 0.5 then
										eigvals_fixed_heur[i,] = result_eigvals[2*i-1,];
									else
										eigvals_fixed_heur[i,] = result_eigvals[2*i,];
								end;

								return eigvals_fixed_heur;  
							end;
						else
							return result_eigvals;		
					end;

				else if beta = 4 then
					do;
						A = RandNormal( N, shape(0, N, 1), I(N));
						B = RandNormal( N, shape(0, N, 1), I(N));
						C = RandNormal( N, shape(0, N, 1), I(N));
						D = RandNormal( N, shape(0, N, 1), I(N));

						Re_diag  =  A + A`;
						Re_ndiag = -C + C;
						Im_diag  =  B - B`;
						Im_ndiag =  D - D`;

						collect_matr = shape(0, 4 * N, 4 * N);

						ind1 = 1:N;
						ind2 = (N+1):(2*N);
						ind3 = (2*N+1):(3*N);
						ind4 = (3*N+1):(4*N);
						/* first column */
						collect_matr[ind1, ind1] = Re_diag;
						collect_matr[ind2, ind1] = Re_ndiag;
						collect_matr[ind3, ind1] = Im_diag;
						collect_matr[ind4, ind1] = Im_ndiag;
						/* second column */
						collect_matr[ind1, ind2] = - Re_ndiag;
						collect_matr[ind2, ind2] = Re_diag;
						collect_matr[ind3, ind2] = Im_ndiag;
						collect_matr[ind4, ind2] = - Im_diag;
						/* third column */
						collect_matr[ind1, ind3] = -Im_diag;
						collect_matr[ind2, ind3] = -Im_ndiag;
						collect_matr[ind3, ind3] = Re_diag;
						collect_matr[ind4, ind3] = Re_ndiag;
						/* fourth column */
						collect_matr[ind1, ind4] = -Im_ndiag;
						collect_matr[ind2, ind4] = Im_diag;
						collect_matr[ind3, ind4] = -Re_ndiag;
						collect_matr[ind4, ind4] = Re_diag;


						eigvecs = eigvec(collect_matr);	
						result_eigvals = eigval(eigvecs); 

						if heurestic_fix then
							do;
								eigvals_fixed_heur = shape(0, N, 2);
								probs_for_heurestic = shape(0, N, 1);
								call randgen(probs_for_heurestic, "Uniform");	
								
								do i=1 to N;
									if probs_for_heurestic[i] < 0.25 then
										eigvals_fixed_heur[i,] = result_eigvals[2*i-1,];

									else if probs_for_heurestic[i] < 0.5 then
										eigvals_fixed_heur[i,] = result_eigvals[2*i,];

									else if probs_for_heurestic[i] < 0.75 then
										eigvals_fixed_heur[i,] = result_eigvals[4*N - 2*i+1,];

									else
										eigvals_fixed_heur[i,] = result_eigvals[4*N - 2*i+2,];
								end;

								return eigvals_fixed_heur; 
							end;
						else
							do;
								return result_eigvals; 
							end;
					end;

				else 
					do;
						print "For Hermite version Beta has to be 1, 2 or 4.";
						return .;
					end;
			end;

		else if haar_mode = "QR" then
			do;
				if beta = 1 then
					do;
						A = RandNormal( N, shape(0, N, 1), I(N));
						call qr(Q, R, piv, indep, A);
						D = vecdiag(R); /* should there be a transposition here? */
      					U = Q # sign(D); 
						result_eigvals = eigval(U); 
						
						return result_eigvals; 
					end;

				else 
					do;
						print "`QR` method not provided for `beta` ither than 1";
						return .;
					end;
			end;

		else 
			do;
				print "Wrong haar_mode supplied --- only `Hermite` and `QR` are available";
				return .;
			end;
	finish;


	start _mu_ref_unif_unit_circle_sampler(beta=2, size=10, normalize=1, heurestic_fix=1, random_state=1618); /* normalization does nothing here */

	
		call streaminit(random_state);
		if beta = int(beta) && sign(beta) = 1 then
			do;
				alpha_re = shape(0, size, 1);
				alpha_im = shape(0, size, 1);
				do i=1 to size;
					nu = 1 + beta * (size - i);
					gauss_vec = RandNormal(1, shape(0, nu+1, 1), I(nu+1));
					alpha_re[i] = gauss_vec[1] / norm(gauss_vec);
					alpha_im[i] = gauss_vec[2] / norm(gauss_vec);
				end;

				rho = sqrt( 1 - alpha_re[1:size-1]##2 - alpha_im[1:size-1]##2 );

				L = shape(0, 2*size, 2*size);
				M = shape(0, 2*size, 2*size);
				
				M[1,1] = 1;
				M[size+1,size+1] = 1;
				if mod(size, 2) = 0 then 
					do;
						do i=1 to int(size/2);
							/* L */

							/* over diagonal */
							L[2*i-1, 2*i] = rho[2*i-1];
							L[2*i-1 + size, 2*i + size] = rho[2*i-1];
							/* diagonal */
							L[2*i-1, 2*i-1] = alpha_re[2*i-1];
							L[2*i, 2*i] = -alpha_re[2*i-1];
							L[2*i-1 + size, 2*i-1 + size] = alpha_re[2*i-1];
							L[2*i + size, 2*i + size] = -alpha_re[2*i-1];

							L[2*i-1 + size, 2*i-1] = -alpha_im[2*i-1];
							L[2*i + size, 2*i] = -alpha_im[2*i-1];
							L[2*i-1 , 2*i-1 + size] = alpha_im[2*i-1];
							L[2*i, 2*i + size] = alpha_im[2*i-1];
							/* under diagonal */
							L[2*i, 2*i-1] = rho[2*i-1];
							L[2*i + size, 2*i-1 + size] = rho[2*i-1];
						end;
						
						
						do i=1 to int(size/2)-1;
							/* M */

							/* over diagonal */
							M[2*i, 2*i+1] = rho[2*i];
							M[2*i + size, 2*i + size+1] = rho[2*i];
							/* diagonal */
							M[2*i, 2*i] = alpha_re[2*i];
							M[2*i+1, 2*i+1] = -alpha_re[2*i];
							M[2*i + size, 2*i + size] = alpha_re[2*i];
							M[2*i+1 + size, 2*i+1 + size] = -alpha_re[2*i];

							M[2*i+ size, 2*i] = -alpha_im[2*i];
							M[2*i+1+ size, 2*i+1] = -alpha_im[2*i];
							M[2*i, 2*i + size] = alpha_im[2*i];
							M[2*i+1, 2*i+1 + size] = alpha_im[2*i];
							/* under diagonal */
							M[2*i+1, 2*i] = rho[2*i];
							M[2*i+1 + size, 2*i + size] = rho[2*i];
						end;
						M[size,size] = alpha_re[size];
						M[2*size, 2*size] = alpha_re[size];
						M[size, 2*size] = alpha_im[size];
						M[2*size, size] = -alpha_im[size];
					end;
				else
					do;
						/* L */
						L[size, size] = alpha_re[size];
						L[2*size, 2*size] = alpha_re[size];
						L[size, 2*size] = alpha_im[size];
						L[2*size, size] = -alpha_im[size];

						do i=1 to int(size/2);
							/* over diagonal */
							L[2*i-1, 2*i] = rho[2*i-1];
							L[2*i-1 + size, 2*i + size] = rho[2*i-1];
							/* diagonal */
							L[2*i-1, 2*i-1] = alpha_re[2*i-1];
							L[2*i, 2*i] = -alpha_re[2*i-1];
							L[2*i-1 + size, 2*i-1 + size] = alpha_re[2*i-1];
							L[2*i + size, 2*i + size] = -alpha_re[2*i-1];

							L[2*i-1 + size, 2*i-1] = -alpha_im[2*i-1];
							L[2*i + size, 2*i] = -alpha_im[2*i-1];
							L[2*i-1 , 2*i-1 + size] = alpha_im[2*i-1];
							L[2*i, 2*i + size] = alpha_im[2*i-1];
							/* under diagonal */
							L[2*i, 2*i-1] = rho[2*i-1];
							L[2*i + size, 2*i-1 + size] = rho[2*i-1];

							/* M */

							/* over diagonal */
							M[2*i, 2*i+1] = rho[2*i];
							M[2*i + size, 2*i + size+1] = rho[2*i];
							/* diagonal */
							M[2*i, 2*i] = alpha_re[2*i];
							M[2*i+1, 2*i+1] = -alpha_re[2*i];
							M[2*i + size, 2*i + size] = alpha_re[2*i];
							M[2*i+1 + size, 2*i + size+1] = -alpha_re[2*i];

							M[2*i + size, 2*i] = -alpha_im[2*i];
							M[2*i+1 + size, 2*i+1] = -alpha_im[2*i];
							M[2*i , 2*i + size] = alpha_im[2*i];
							M[2*i+1, 2*i + size+1] = alpha_im[2*i];
							/* under diagonal */
							M[2*i+1, 2*i] = rho[2*i];
							M[2*i+1 + size, 2*i + size] = rho[2*i];
						end;
					end;

				result_eigvals = eigval(L*M); 

				if heurestic_fix then
							do;
								eigvals_fixed_heur = shape(0, size, 2);
								probs_for_heurestic = shape(0, size, 1);
								call randgen(probs_for_heurestic, "Uniform");	
								
								

								do i=1 to size;
									if probs_for_heurestic[i] < 0.5 then
										eigvals_fixed_heur[i,] = result_eigvals[2*i-1,];
									else
										eigvals_fixed_heur[i,] = result_eigvals[2*i,];
								end;

								return eigvals_fixed_heur; 
							end;
				else
							do;
								return result_eigvals; 
							end;  
			end;
		else
			do;
				print "`beta` must be positive integer.";
				return .;
			end;	
	
	finish;

	start _ginibre_sampler_full(N, normalize=1, heurestic_fix=1, random_state=1618);
  
	
		call streaminit(random_state);
		A = RandNormal( N, shape(0, N, 1), I(N));
		B = RandNormal( N, shape(0, N, 1), I(N));

    	collect_matrix = shape(0, 2*N, 2*N);
		ind1 = 1:N; ind2 = (N+1):(2*N);
		collect_matrix[ind1, ind1] = A;
		collect_matrix[ind2, ind1] = B;
		collect_matrix[ind1, ind2] = -B;
		collect_matrix[ind2, ind2] = A;

		result_eigvals = eigval(collect_matrix)/sqrt(2);

		if normalize then
			do;
				result_eigvals = result_eigvals / sqrt(N);
			end;

		if heurestic_fix then
			do;
				eigvals_fixed_heur = shape(0, N, 2);
				probs_for_heurestic = shape(0, N, 1);
				call randgen(probs_for_heurestic, "Uniform");	
												
				do i=1 to N;
					if probs_for_heurestic[i] < 0.5 then
						eigvals_fixed_heur[i,] = result_eigvals[2*i-1,];
					else
						eigvals_fixed_heur[i,] = result_eigvals[2*i,];
				end;

				return eigvals_fixed_heur; 
			end;
		else
			do;
				return result_eigvals;
			end;

		
	finish;


start sample_from_beta_ensemble_full( 
	ensemble_version, 
	M_1, M_2, /* M variables only for Laguerre (first) and Jacobi (both) */
	size=10, 
	beta=2, 
	normalize=1, 
	haar_mode="Hermite", /* haar_mode only available for circular ensemble */
	heurestic_fix=1, /* heurestic_fix only available for circular ensemble */
	random_state=1618); 
	
	if ensemble_version = "Hermite" then
		return _hermite_sampler_full(size, beta, normalize, random_state);

	else if ensemble_version = "Laguerre" then
		return _laguerre_sampler_full(M_1, size, beta, normalize, random_state); 

	else if ensemble_version = "Jacobi" then
		return _jacobi_sampler_full(M_1, M_2, size, beta, normalize, random_state); 

	else if ensemble_version = "Circular" then
		return _circular_sampler_full(size, beta, haar_mode, normalize, heurestic_fix, random_state); /*normalization does nothing here */
    
	else if ensemble_version = "Ginibre" then
		return _ginibre_sampler_full(size, normalize, heurestic_fix, random_state);

	else
		do;
			print "Wrong type supplied - only Hermite, Laguerre, Jacobi, and Circular are available.";
			return .;
		end;

finish;

start sample_from_beta_ensemble_banded(
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
	
	if ensemble_version = "Hermite" then 
			return _mu_ref_normal_sampler_tridiag(size, loc, scale, beta, normalize, random_state);
		
	else if ensemble_version = "Laguerre" then
			return _mu_ref_gamma_sampler_tridiag(shape, scale, beta, size, normalize, random_state);

	else if ensemble_version = "Jacobi" then
			return _mu_ref_beta_sampler_tridiag(a, b, beta, size, normalize, random_state); /* normalization does nothing here */
	
	else if ensemble_version = "Circular" then
			return _mu_ref_unif_unit_circle_sampler(beta, size, normalize, heurestic_fix, random_state); /* normalization does nothing here */

	else
 		do;
			print "Wrong type supplied - only Hermite, Laguerre, Jacobi, and Circular are available.";
			return .;
		end;

finish;

	ensemble_version = "Circular";
	size=10;
	beta=2;
	loc=0.0;
	scale=1.0;
	shape = 1.0;
	a = 1.0;
	b = 1.0;
	normalize=0;
	heurestic_fix=0;
	random_state=1618;

	sample = sample_from_beta_ensemble_banded(
	ensemble_version, 
	size, 
	beta, 
	/* distribution specific parameters */
	loc, 
	scale, 
	shape,
	a,
	b,
	/* parameters for the algorithm */
	normalize, 
	heurestic_fix,
	random_state);

	run scatter(sample[,1], sample[,2]);
quit;




proc iml;

start _mu_ref_normal_sampler_tridiag(size=10, loc=0.0, scale=1.0, beta=2, normalize=1, random_state=1618);

	call streaminit(random_state);
		if beta <= 0 then
			do;
				result_eigvals = .;
        		print "`beta` must be positive.";
			end;
		else 
			do;
				M = diag( RandNormal(1, shape(loc, size, 1), I(size))); /* diagonal */
				do i = 1 to size-1;
					x = sqrt( scale * scale * rand("GAMMA", beta * (size - i)/2));
					M[i, i+1] = x;
					M[i+1, i] = x;
				end;

				result_eigvals = eigval(M);

				if normalize then
					do;
						result_eigvals = (result_eigvals - loc) / (sqrt(0.5 * beta * size) * scale) ;
					end;  		
			end;
	return result_eigvals;

finish;

	start _hermite_sampler_full(size=10, beta=2, normalize=1, random_state=1618);

		call streaminit(random_state);
		result_eigvals = shape(0, size, 1);

		if beta = 1 then
			do;
				A = RandNormal( size, shape(0, size, 1), I(size));
				result_eigvals = eigval( A + A` )/sqrt(2) ;
			end;

		else if beta = 2 then
			do;
				A = RandNormal(size, shape(0, size, 1), I(size));
				B = RandNormal(size, shape(0, size, 1), I(size));

				Re_mat = A + A`;
				Im_mat = B - B`;

				collect_matr = shape(0, 2 * size, 2 * size);

				ind1 = 1:size;
				ind2 = (size+1):(2*size);
				/* first column */
				collect_matr[ind1, ind1] = Re_mat;
				collect_matr[ind2, ind1] = Im_mat;
				/* second column */
				collect_matr[ind1, ind2] = - Im_mat;
				collect_matr[ind2, ind2] = Re_mat;

				chosen_indices = do(2, 2*size, 2);
				result_eigvals = eigval(collect_matr)[chosen_indices] / sqrt(2); /* eigevalues get doubled --- we select only one from each pair */ 
			end;

		else if beta = 4 then
			do;
				A = RandNormal(size, shape(0, size, 1), I(size));
				B = RandNormal(size, shape(0, size, 1), I(size));
				C = RandNormal(size, shape(0, size, 1), I(size));
				D = RandNormal(size, shape(0, size, 1), I(size));

				Re_diag  =  A + A`;
				Re_ndiag = -C + C;
				Im_diag  =  B - B`;
				Im_ndiag =  D - D`;

				collect_matr = shape(0, 4 * size, 4 * size);

				ind1 = 1:size;
				ind2 = (size+1):(2*size);
				ind3 = (2*size+1):(3*size);
				ind4 = (3*size+1):(4*size);
				/* first column */
				collect_matr[ind1, ind1] = Re_diag;
				collect_matr[ind2, ind1] = Re_ndiag;
				collect_matr[ind3, ind1] = Im_diag;
				collect_matr[ind4, ind1] = Im_ndiag;
				/* second column */
				collect_matr[ind1, ind2] = - Re_ndiag;
				collect_matr[ind2, ind2] = Re_diag;
				collect_matr[ind3, ind2] = Im_ndiag;
				collect_matr[ind4, ind2] = - Im_diag;
				/* third column */
				collect_matr[ind1, ind3] = -Im_diag;
				collect_matr[ind2, ind3] = -Im_ndiag;
				collect_matr[ind3, ind3] = Re_diag;
				collect_matr[ind4, ind3] = Re_ndiag;
				/* third column */
				collect_matr[ind1, ind4] = -Im_ndiag;
				collect_matr[ind2, ind4] = Im_diag;
				collect_matr[ind3, ind4] = -Re_ndiag;
				collect_matr[ind4, ind4] = Re_diag;

				chosen_indices = do(4, 4*size, 4);
				result_eigvals = eigval(collect_matr)[chosen_indices] / sqrt(2); /* eigevalues get quadrupled --- we select only one from each four */
			end;

		else 
			do;
				print "Beta has to be 1, 2 or 4.";
				result_eigvals = .;
			end;

		if normalize then
			do;
				
				result_eigvals = result_eigvals / sqrt( beta * size );
			end;
		
		return result_eigvals;
	finish;

	start _laguerre_sampler_full(M, N=10, beta=2, normalize=1, random_state=1618); /* maybe change N to size */

		call streaminit(random_state);

		result_eigvals = shape(0, N, 1);
		if beta = 1 then
			do;
				A = RandNormal( N, shape(0, M, 1), I(M));
				result_eigvals = eigval(A * A`);
			end;

		else if beta = 2 then
			do;
				A = RandNormal( N, shape(0, M, 1), I(M));
				B = RandNormal( N, shape(0, M, 1), I(M));
				Re_mat = A * A` + B * B`;
				Im_mat = B * A` - A * B`;

				collect_matr = shape(0, 2 * N, 2 * N);

				ind1 = 1:N;
				ind2 = (N+1):(2*N);
				/* first column */
				collect_matr[ind1, ind1] = Re_mat;
				collect_matr[ind2, ind1] = Im_mat;
				/* second column */
				collect_matr[ind1, ind2] = - Im_mat;
				collect_matr[ind2, ind2] = Re_mat;

				chosen_indices = do(2, 2*N, 2);
				result_eigvals = eigval(collect_matr)[chosen_indices]; /* eigevalues get doubled --- we select only one from each pair */
			end;

		else if beta = 4 then
			do;
				A = RandNormal( N, shape(0, M, 1), I(M));
				B = RandNormal( N, shape(0, M, 1), I(M));
				C = RandNormal( N, shape(0, M, 1), I(M));
				D = RandNormal( N, shape(0, M, 1), I(M));

				Re_diag  = A*A` + B*B` + C*C` + D*D`;
				Re_ndiag = A*C` + D*B` - C*A` - B*D`;
				Im_diag  = B*A` - A*B` + D*C` - C*D`;
				Im_ndiag = D*A` + C*B` - B*C` - A*D`;

				collect_matr = shape(0, 4 * N, 4 * N);

				ind1 = 1:N;
				ind2 = (N+1):(2*N);
				ind3 = (2*N+1):(3*N);
				ind4 = (3*N+1):(4*N);
				/* first column */
				collect_matr[ind1, ind1] = Re_diag;
				collect_matr[ind2, ind1] = Re_ndiag;
				collect_matr[ind3, ind1] = Im_diag;
				collect_matr[ind4, ind1] = Im_ndiag;
				/* second column */
				collect_matr[ind1, ind2] = - Re_ndiag;
				collect_matr[ind2, ind2] = Re_diag;
				collect_matr[ind3, ind2] = Im_ndiag;
				collect_matr[ind4, ind2] = - Im_diag;
				/* third column */
				collect_matr[ind1, ind3] = -Im_diag;
				collect_matr[ind2, ind3] = -Im_ndiag;
				collect_matr[ind3, ind3] = Re_diag;
				collect_matr[ind4, ind3] = Re_ndiag;
				/* fourth column */
				collect_matr[ind1, ind4] = -Im_ndiag;
				collect_matr[ind2, ind4] = Im_diag;
				collect_matr[ind3, ind4] = -Re_ndiag;
				collect_matr[ind4, ind4] = Re_diag;
				
				chosen_indices = do(4, 4*N, 4);		
				result_eigvals = eigval(collect_matr)[chosen_indices]; /* eigevalues get quadrupled --- we select only one from each four */
				
			end;

		else 
			do;
				print "Beta has to be 1, 2 or 4.";
				result_eigvals = .;
			end;

		if normalize then
			do;
				result_eigvals = result_eigvals / (beta * M); 
			end;

		return result_eigvals; 
	finish;

	start _mu_ref_gamma_sampler_tridiag(shape=1.0, scale=1.0, beta=2, size=10, normalize=1, random_state=1618);

		call streaminit(random_state);

		result_eigvals = shape(0, size, 1);
		if beta <= 0 then
			do;
        		print "`beta` must be positive.";
				result_eigvals = .;
			end;
		else 
			do;
				xi_odd  = shape(0, size, 1);
				xi_even = shape(0, size, 1);

				do i = 1 to size-1;
					xi_odd[i]    = scale * rand("GAMMA", beta * (size - i)/2 + shape);
					xi_even[i+1] = scale * rand("GAMMA", beta * (size - i)/2);
				end;
				xi_odd[size]   = scale * rand("GAMMA", shape);

				M = diag( xi_odd + xi_even);

				do i = 1 to size-1;
					M[i+1,i] = sqrt(xi_odd[i] * xi_even[i+1]);
					M[i,i+1] = sqrt(xi_odd[i] * xi_even[i+1]);
				end;

				result_eigvals = eigval(M)[,1];
			
				if normalize then
					do; /* potencjalnie cos nie tak */
						result_eigvals = result_eigvals / (0.5 *scale * beta * (2 / beta * shape + size - 1));
					end;	
			end;

		return result_eigvals; 	
	finish;

	start _jacobi_sampler_full(M_1, M_2, N=10, beta=2, normalize=1, random_state=1618); /* normalization does nothing here */
	
		call streaminit(random_state);

		
		if beta = 1 then
			do;
				X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				Y = RandNormal( N, shape(0, M_2, 1), I(M_2));

				X_tmp = X*X`;
				Y_tmp = Y*Y`;

				res_tmp = X_tmp * inv( X_tmp + Y_tmp );
				result_eigvals = eigval(res_tmp)[,1];

				return result_eigvals; 
			end;

		else if beta = 2 then
			do;
				A_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				B_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				A_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));
				B_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));

				A_X_tmp = A_X*A_X` + B_X*B_X`;
				B_X_tmp = B_X*A_X` - A_X*B_X`;
				A_Y_tmp = A_Y*A_Y` + B_Y*B_Y`;
				B_Y_tmp = B_Y*A_Y` - A_Y*B_Y`;

				A_to_solve = shape(0, 2*N, 2*N);
				B_to_solve = shape(0, 2*N, N);


				ind1 = 1:N;
				ind2 = (N+1):(2*N);
				/* first column */
				A_to_solve[ind1, ind1] = A_X_tmp` + A_Y_tmp`;
				A_to_solve[ind2, ind1] = B_X_tmp` + B_Y_tmp`;
				/* second column */
				A_to_solve[ind1, ind2] = -(B_X_tmp` + B_Y_tmp`);
				A_to_solve[ind2, ind2] = A_X_tmp` + A_Y_tmp`;

				B_to_solve[ind1, ind1] = A_X_tmp`;
				B_to_solve[ind2, ind1] = B_X_tmp`;

				quot = solve(A_to_solve, B_to_solve)`;

				collect_mat = shape(0, 2*N, 2*N);
				/* first column */
				collect_mat[ind1, ind1] = quot[ind1, ind1];
				collect_mat[ind2, ind1] = quot[ind1, ind2];
				/* second column */
				collect_mat[ind1, ind2] = -quot[ind1, ind2];
				collect_mat[ind2, ind2] =  quot[ind1, ind1];

				chosen_indices = do(2, 2*N, 2);
				result_eigvals = eigval(collect_mat)[chosen_indices, 1];
				
				return result_eigvals; 
			end;

		else if beta = 4 then
			do;
				A_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				B_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				C_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				D_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				A_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));
				B_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));
				C_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));
				D_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));

				/* Changes in sign are because of transposition */
				Re_diag_X_t  =  A_X*A_X` + B_X*B_X` + C_X*C_X` + D_X*D_X`;
				Re_ndiag_X_t = -A_X*C_X` - D_X*B_X` + C_X*A_X` + B_X*D_X`;
				Im_diag_X_t  = -B_X*A_X` + A_X*B_X` - D_X*C_X` + C_X*D_X`;
				Im_ndiag_X_t = -D_X*A_X` - C_X*B_X` + B_X*C_X` + A_X*D_X`;

				Re_diag_tot_t  = Re_diag_X_t  + A_Y*A_Y` + B_Y*B_Y` + C_Y*C_Y` + D_Y*D_Y`;
				Re_ndiag_tot_t = Re_ndiag_X_t - A_Y*C_Y` - D_Y*B_Y` + C_Y*A_Y` + B_Y*D_Y`;
				Im_diag_tot_t  = Im_diag_X_t  - B_Y*A_Y` + A_Y*B_Y` - D_Y*C_Y` + C_Y*D_Y`;
				Im_ndiag_tot_t = Im_ndiag_X_t - D_Y*A_Y` - C_Y*B_Y` + B_Y*C_Y` + A_Y*D_Y`;

				A_to_solve = shape(0, 4*N, 4*N);
				B_to_solve = shape(0, 4*N, 2*N);

				ind1 = 1:N; 
				ind2 = (N+1):(2*N);
 				ind3 = (2*N+1):(3*N);
				ind4 = (3*N+1):(4*N);
				/* first column */
				A_to_solve[ind1, ind1] = Re_diag_tot_t;
				A_to_solve[ind2, ind1] = - Re_ndiag_tot_t;
				A_to_solve[ind3, ind1] = Im_diag_tot_t;
				A_to_solve[ind4, ind1] = Im_ndiag_tot_t;
				/* second column */
				A_to_solve[ind1, ind2] =   Re_ndiag_tot_t;
				A_to_solve[ind2, ind2] =   Re_diag_tot_t;
				A_to_solve[ind3, ind2] =   Im_ndiag_tot_t;
				A_to_solve[ind4, ind2] = - Im_diag_tot_t;
				/* third column */
				A_to_solve[ind1, ind3] = - Im_diag_tot_t;
				A_to_solve[ind2, ind3] = - Im_ndiag_tot_t;
				A_to_solve[ind3, ind3] =   Re_diag_tot_t;
				A_to_solve[ind4, ind3] = - Re_ndiag_tot_t;
				/* fourth column */
				A_to_solve[ind1, ind4] = - Im_ndiag_tot_t;
				A_to_solve[ind2, ind4] =   Im_diag_tot_t;
				A_to_solve[ind3, ind4] = - Re_ndiag_tot_t;
				A_to_solve[ind4, ind4] =   Re_diag_tot_t;
	
				/* first column */
				B_to_solve[ind1, ind1] =   Re_diag_X_t;
				B_to_solve[ind2, ind1] = - Re_ndiag_X_t;
				B_to_solve[ind3, ind1] =   Im_diag_X_t;
				B_to_solve[ind4, ind1] =   Im_ndiag_X_t;
				/* second column */
				B_to_solve[ind1, ind2] =   Re_ndiag_X_t;
				B_to_solve[ind2, ind2] =   Re_diag_X_t;
				B_to_solve[ind3, ind2] =   Im_ndiag_X_t;
				B_to_solve[ind4, ind2] = - Im_diag_X_t;

				quot = solve(A_to_solve, B_to_solve);

				collect_mat = shape(0, 4*N, 4*N);
				Lind1 = 1:(2*N);
				Lind2 = (2*N+1):(4*N);


				/* first column */
				collect_mat[Lind1, Lind1] = quot[Lind1, Lind1];
				collect_mat[Lind2, Lind1] = quot[Lind2, Lind1];
				/* second column */
				collect_mat[Lind1, Lind2] = -quot[Lind2, Lind1];
				collect_mat[Lind2, Lind2] = quot[Lind1, Lind1];

				
				chosen_indices = do(4, 4*N, 4);
				result_eigvals = eigval(collect_mat)[chosen_indices, 1];
				
				return result_eigvals; 
			end;

		else 
			do;
				print "Beta has to be 1, 2 or 4.";
				result_eigvals = .;
				return result_eigvals;
			end;

	finish;

	start _mu_ref_beta_sampler_tridiag(a, b, beta=2, size=10, normalize=1, random_state=1618); /* normalization does nothing here */
	
	
		call streaminit(random_state);

		if beta <= 0 then
			do;
        		print "`beta` must be positive.";
				result_eigvals = .;
				return result_eigvals;
			end;
		else 
			do;
				c_odd  = shape(0, size, 1);
				c_even = shape(0, size, 1);
				do i=1 to size-1;
					c_odd[i]    = rand("BETA", beta*(size - i)/2 + a, beta*(size - i)/2 + b);
					c_even[i+1] = rand("BETA", beta*(size - i)/2, beta*(size - i - 1)/2 + b + a);
				end;
				c_odd[size] = rand("BETA", a, b);

				xi_odd = (1 - c_even) # c_odd;
				xi_even = shape(0, size, 1);
				xi_even[2:size] = (1 - c_odd[1:size-1]) # c_even[2:size];

				alpha_coef = xi_even + xi_odd;
				beta_coef = sqrt(xi_odd[1:size-1] # xi_even[2:size]);
	
				M = diag( alpha_coef);
				do i=1 to size-1;
					M[i, i+1] = beta_coef[i];
					M[i+1, i] = beta_coef[i];
				end;

				result_eigvals = eigval(M)[,1];
				return result_eigvals;
			end;
	finish;


	start _circular_sampler_full(N, beta=2, haar_mode='QR', normalize=1, heurestic_fix=1, random_state=1618); /*normalization does nothing here */
    	
		call streaminit(random_state);

		if haar_mode = "Hermite" then
			do;

				if beta = 1 then
					do;
						A = RandNormal( N, shape(0, N, 1), I(N));
						eigvecs = eigvec( A + A` );
						result_eigvals = eigval( eigvecs);
						return result_eigvals;
					end;

				else if beta = 2 then
					do;
						A = RandNormal( N, shape(0, N, 1), I(N));
						B = RandNormal( N, shape(0, N, 1), I(N));

						Re_mat = A + A`;
						Im_mat = B - B`;

						collect_matr = shape(0, 2 * N, 2 * N);
						ind1 = 1:N;
						ind2 = (N+1):(2*N);
						/* first column */
						collect_matr[ind1, ind1] = Re_mat;
						collect_matr[ind2, ind1] = Im_mat;
						/* second column */
						collect_matr[ind1, ind2] = - Im_mat;
						collect_matr[ind2, ind2] = Re_mat;

						eigvecs = eigvec(collect_matr);
						result_eigvals = eigval(eigvecs); 


						if heurestic_fix then
							do;
								eigvals_fixed_heur = shape(0, N, 2);
								probs_for_heurestic = shape(0, N, 1);
								call randgen(probs_for_heurestic, "Uniform");	

								do i=1 to N;
									if probs_for_heurestic[i] < 0.5 then
										eigvals_fixed_heur[i,] = result_eigvals[2*i-1,];
									else
										eigvals_fixed_heur[i,] = result_eigvals[2*i,];
								end;

								return eigvals_fixed_heur;  
							end;
						else
							return result_eigvals;		
					end;

				else if beta = 4 then
					do;
						A = RandNormal( N, shape(0, N, 1), I(N));
						B = RandNormal( N, shape(0, N, 1), I(N));
						C = RandNormal( N, shape(0, N, 1), I(N));
						D = RandNormal( N, shape(0, N, 1), I(N));

						Re_diag  =  A + A`;
						Re_ndiag = -C + C;
						Im_diag  =  B - B`;
						Im_ndiag =  D - D`;

						collect_matr = shape(0, 4 * N, 4 * N);

						ind1 = 1:N;
						ind2 = (N+1):(2*N);
						ind3 = (2*N+1):(3*N);
						ind4 = (3*N+1):(4*N);
						/* first column */
						collect_matr[ind1, ind1] = Re_diag;
						collect_matr[ind2, ind1] = Re_ndiag;
						collect_matr[ind3, ind1] = Im_diag;
						collect_matr[ind4, ind1] = Im_ndiag;
						/* second column */
						collect_matr[ind1, ind2] = - Re_ndiag;
						collect_matr[ind2, ind2] = Re_diag;
						collect_matr[ind3, ind2] = Im_ndiag;
						collect_matr[ind4, ind2] = - Im_diag;
						/* third column */
						collect_matr[ind1, ind3] = -Im_diag;
						collect_matr[ind2, ind3] = -Im_ndiag;
						collect_matr[ind3, ind3] = Re_diag;
						collect_matr[ind4, ind3] = Re_ndiag;
						/* fourth column */
						collect_matr[ind1, ind4] = -Im_ndiag;
						collect_matr[ind2, ind4] = Im_diag;
						collect_matr[ind3, ind4] = -Re_ndiag;
						collect_matr[ind4, ind4] = Re_diag;


						eigvecs = eigvec(collect_matr);	
						result_eigvals = eigval(eigvecs); 

						if heurestic_fix then
							do;
								eigvals_fixed_heur = shape(0, N, 2);
								probs_for_heurestic = shape(0, N, 1);
								call randgen(probs_for_heurestic, "Uniform");	
								
								do i=1 to N;
									if probs_for_heurestic[i] < 0.25 then
										eigvals_fixed_heur[i,] = result_eigvals[2*i-1,];

									else if probs_for_heurestic[i] < 0.5 then
										eigvals_fixed_heur[i,] = result_eigvals[2*i,];

									else if probs_for_heurestic[i] < 0.75 then
										eigvals_fixed_heur[i,] = result_eigvals[4*N - 2*i+1,];

									else
										eigvals_fixed_heur[i,] = result_eigvals[4*N - 2*i+2,];
								end;

								return eigvals_fixed_heur; 
							end;
						else
							do;
								return result_eigvals; 
							end;
					end;

				else 
					do;
						print "For Hermite version Beta has to be 1, 2 or 4.";
						return .;
					end;
			end;

		else if haar_mode = "QR" then
			do;
				if beta = 1 then
					do;
						A = RandNormal( N, shape(0, N, 1), I(N));
						call qr(Q, R, piv, indep, A);
						D = vecdiag(R); /* should there be a transposition here? */
      					U = Q # sign(D); 
						result_eigvals = eigval(U); 
						
						return result_eigvals; 
					end;

				else 
					do;
						print "`QR` method not provided for `beta` ither than 1";
						return .;
					end;
			end;

		else 
			do;
				print "Wrong haar_mode supplied --- only `Hermite` and `QR` are available";
				return .;
			end;
	finish;


	start _mu_ref_unif_unit_circle_sampler(beta=2, size=10, normalize=1, heurestic_fix=1, random_state=1618); /* normalization does nothing here */

	
		call streaminit(random_state);
		if beta = int(beta) && sign(beta) = 1 then
			do;
				alpha_re = shape(0, size, 1);
				alpha_im = shape(0, size, 1);
				do i=1 to size;
					nu = 1 + beta * (size - i);
					gauss_vec = RandNormal(1, shape(0, nu+1, 1), I(nu+1));
					alpha_re[i] = gauss_vec[1] / norm(gauss_vec);
					alpha_im[i] = gauss_vec[2] / norm(gauss_vec);
				end;

				rho = sqrt( 1 - alpha_re[1:size-1]##2 - alpha_im[1:size-1]##2 );

				L = shape(0, 2*size, 2*size);
				M = shape(0, 2*size, 2*size);
				
				M[1,1] = 1;
				M[size+1,size+1] = 1;
				if mod(size, 2) = 0 then 
					do;
						do i=1 to int(size/2);
							/* L */

							/* over diagonal */
							L[2*i-1, 2*i] = rho[2*i-1];
							L[2*i-1 + size, 2*i + size] = rho[2*i-1];
							/* diagonal */
							L[2*i-1, 2*i-1] = alpha_re[2*i-1];
							L[2*i, 2*i] = -alpha_re[2*i-1];
							L[2*i-1 + size, 2*i-1 + size] = alpha_re[2*i-1];
							L[2*i + size, 2*i + size] = -alpha_re[2*i-1];

							L[2*i-1 + size, 2*i-1] = -alpha_im[2*i-1];
							L[2*i + size, 2*i] = -alpha_im[2*i-1];
							L[2*i-1 , 2*i-1 + size] = alpha_im[2*i-1];
							L[2*i, 2*i + size] = alpha_im[2*i-1];
							/* under diagonal */
							L[2*i, 2*i-1] = rho[2*i-1];
							L[2*i + size, 2*i-1 + size] = rho[2*i-1];
						end;
						
						
						do i=1 to int(size/2)-1;
							/* M */

							/* over diagonal */
							M[2*i, 2*i+1] = rho[2*i];
							M[2*i + size, 2*i + size+1] = rho[2*i];
							/* diagonal */
							M[2*i, 2*i] = alpha_re[2*i];
							M[2*i+1, 2*i+1] = -alpha_re[2*i];
							M[2*i + size, 2*i + size] = alpha_re[2*i];
							M[2*i+1 + size, 2*i+1 + size] = -alpha_re[2*i];

							M[2*i+ size, 2*i] = -alpha_im[2*i];
							M[2*i+1+ size, 2*i+1] = -alpha_im[2*i];
							M[2*i, 2*i + size] = alpha_im[2*i];
							M[2*i+1, 2*i+1 + size] = alpha_im[2*i];
							/* under diagonal */
							M[2*i+1, 2*i] = rho[2*i];
							M[2*i+1 + size, 2*i + size] = rho[2*i];
						end;
						M[size,size] = alpha_re[size];
						M[2*size, 2*size] = alpha_re[size];
						M[size, 2*size] = alpha_im[size];
						M[2*size, size] = -alpha_im[size];
					end;
				else
					do;
						/* L */
						L[size, size] = alpha_re[size];
						L[2*size, 2*size] = alpha_re[size];
						L[size, 2*size] = alpha_im[size];
						L[2*size, size] = -alpha_im[size];

						do i=1 to int(size/2);
							/* over diagonal */
							L[2*i-1, 2*i] = rho[2*i-1];
							L[2*i-1 + size, 2*i + size] = rho[2*i-1];
							/* diagonal */
							L[2*i-1, 2*i-1] = alpha_re[2*i-1];
							L[2*i, 2*i] = -alpha_re[2*i-1];
							L[2*i-1 + size, 2*i-1 + size] = alpha_re[2*i-1];
							L[2*i + size, 2*i + size] = -alpha_re[2*i-1];

							L[2*i-1 + size, 2*i-1] = -alpha_im[2*i-1];
							L[2*i + size, 2*i] = -alpha_im[2*i-1];
							L[2*i-1 , 2*i-1 + size] = alpha_im[2*i-1];
							L[2*i, 2*i + size] = alpha_im[2*i-1];
							/* under diagonal */
							L[2*i, 2*i-1] = rho[2*i-1];
							L[2*i + size, 2*i-1 + size] = rho[2*i-1];

							/* M */

							/* over diagonal */
							M[2*i, 2*i+1] = rho[2*i];
							M[2*i + size, 2*i + size+1] = rho[2*i];
							/* diagonal */
							M[2*i, 2*i] = alpha_re[2*i];
							M[2*i+1, 2*i+1] = -alpha_re[2*i];
							M[2*i + size, 2*i + size] = alpha_re[2*i];
							M[2*i+1 + size, 2*i + size+1] = -alpha_re[2*i];

							M[2*i + size, 2*i] = -alpha_im[2*i];
							M[2*i+1 + size, 2*i+1] = -alpha_im[2*i];
							M[2*i , 2*i + size] = alpha_im[2*i];
							M[2*i+1, 2*i + size+1] = alpha_im[2*i];
							/* under diagonal */
							M[2*i+1, 2*i] = rho[2*i];
							M[2*i+1 + size, 2*i + size] = rho[2*i];
						end;
					end;

				result_eigvals = eigval(L*M); 

				if heurestic_fix then
							do;
								eigvals_fixed_heur = shape(0, size, 2);
								probs_for_heurestic = shape(0, size, 1);
								call randgen(probs_for_heurestic, "Uniform");	
								
								

								do i=1 to size;
									if probs_for_heurestic[i] < 0.5 then
										eigvals_fixed_heur[i,] = result_eigvals[2*i-1,];
									else
										eigvals_fixed_heur[i,] = result_eigvals[2*i,];
								end;

								return eigvals_fixed_heur; 
							end;
				else
							do;
								return result_eigvals; 
							end;  
			end;
		else
			do;
				print "`beta` must be positive integer.";
				return .;
			end;	
	
	finish;

	start _ginibre_sampler_full(N, normalize=1, heurestic_fix=1, random_state=1618);
  
	
		call streaminit(random_state);
		A = RandNormal( N, shape(0, N, 1), I(N));
		B = RandNormal( N, shape(0, N, 1), I(N));

    	collect_matrix = shape(0, 2*N, 2*N);
		ind1 = 1:N; ind2 = (N+1):(2*N);
		collect_matrix[ind1, ind1] = A;
		collect_matrix[ind2, ind1] = B;
		collect_matrix[ind1, ind2] = -B;
		collect_matrix[ind2, ind2] = A;

		result_eigvals = eigval(collect_matrix)/sqrt(2);

		if normalize then
			do;
				result_eigvals = result_eigvals / sqrt(N);
			end;

		if heurestic_fix then
			do;
				eigvals_fixed_heur = shape(0, N, 2);
				probs_for_heurestic = shape(0, N, 1);
				call randgen(probs_for_heurestic, "Uniform");	
												
				do i=1 to N;
					if probs_for_heurestic[i] < 0.5 then
						eigvals_fixed_heur[i,] = result_eigvals[2*i-1,];
					else
						eigvals_fixed_heur[i,] = result_eigvals[2*i,];
				end;

				return eigvals_fixed_heur; 
			end;
		else
			do;
				return result_eigvals;
			end;

		
	finish;


start sample_from_beta_ensemble_full( 
	ensemble_version, 
	M_1, M_2, /* M variables only for Laguerre (first) and Jacobi (both) */
	size=10, 
	beta=2, 
	normalize=1, 
	haar_mode="Hermite", /* haar_mode only available for circular ensemble */
	heurestic_fix=1, /* heurestic_fix only available for circular ensemble */
	random_state=1618); 
	
	if ensemble_version = "Hermite" then
		return _hermite_sampler_full(size, beta, normalize, random_state);

	else if ensemble_version = "Laguerre" then
		return _laguerre_sampler_full(M_1, size, beta, normalize, random_state); 

	else if ensemble_version = "Jacobi" then
		return _jacobi_sampler_full(M_1, M_2, size, beta, normalize, random_state); 

	else if ensemble_version = "Circular" then
		return _circular_sampler_full(size, beta, haar_mode, normalize, heurestic_fix, random_state); /*normalization does nothing here */
    
	else if ensemble_version = "Ginibre" then
		return _ginibre_sampler_full(size, normalize, heurestic_fix, random_state);

	else
		do;
			print "Wrong type supplied - only Hermite, Laguerre, Jacobi, and Circular are available.";
			return .;
		end;

finish;

start sample_from_beta_ensemble_banded(
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
	
	if ensemble_version = "Hermite" then 
			return _mu_ref_normal_sampler_tridiag(size, loc, scale, beta, normalize, random_state);
		
	else if ensemble_version = "Laguerre" then
			return _mu_ref_gamma_sampler_tridiag(shape, scale, beta, size, normalize, random_state);

	else if ensemble_version = "Jacobi" then
			return _mu_ref_beta_sampler_tridiag(a, b, beta, size, normalize, random_state); /* normalization does nothing here */
	
	else if ensemble_version = "Circular" then
			return _mu_ref_unif_unit_circle_sampler(beta, size, normalize, heurestic_fix, random_state); /* normalization does nothing here */

	else
 		do;
			print "Wrong type supplied - only Hermite, Laguerre, Jacobi, and Circular are available.";
			return .;
		end;

finish;

	ensemble_version = "Circular";
	size=10;
	beta=20;
	loc=0.0;
	scale=1.0;
	shape = 1.0;
	a = 1.0;
	b = 1.0;
	normalize=0;
	heurestic_fix=0;
	random_state=1618;

	sample = sample_from_beta_ensemble_banded(
	ensemble_version, 
	size, 
	beta, 
	/* distribution specific parameters */
	loc, 
	scale, 
	shape,
	a,
	b,
	/* parameters for the algorithm */
	normalize, 
	heurestic_fix,
	random_state);

	run scatter(sample[,1], sample[,2]);
quit;





proc iml;

start _mu_ref_normal_sampler_tridiag(size=10, loc=0.0, scale=1.0, beta=2, normalize=1, random_state=1618);

	call streaminit(random_state);
		if beta <= 0 then
			do;
				result_eigvals = .;
        		print "`beta` must be positive.";
			end;
		else 
			do;
				M = diag( RandNormal(1, shape(loc, size, 1), I(size))); /* diagonal */
				do i = 1 to size-1;
					x = sqrt( scale * scale * rand("GAMMA", beta * (size - i)/2));
					M[i, i+1] = x;
					M[i+1, i] = x;
				end;

				result_eigvals = eigval(M);

				if normalize then
					do;
						result_eigvals = (result_eigvals - loc) / (sqrt(0.5 * beta * size) * scale) ;
					end;  		
			end;
	return result_eigvals;

finish;

	start _hermite_sampler_full(size=10, beta=2, normalize=1, random_state=1618);

		call streaminit(random_state);
		result_eigvals = shape(0, size, 1);

		if beta = 1 then
			do;
				A = RandNormal( size, shape(0, size, 1), I(size));
				result_eigvals = eigval( A + A` )/sqrt(2) ;
			end;

		else if beta = 2 then
			do;
				A = RandNormal(size, shape(0, size, 1), I(size));
				B = RandNormal(size, shape(0, size, 1), I(size));

				Re_mat = A + A`;
				Im_mat = B - B`;

				collect_matr = shape(0, 2 * size, 2 * size);

				ind1 = 1:size;
				ind2 = (size+1):(2*size);
				/* first column */
				collect_matr[ind1, ind1] = Re_mat;
				collect_matr[ind2, ind1] = Im_mat;
				/* second column */
				collect_matr[ind1, ind2] = - Im_mat;
				collect_matr[ind2, ind2] = Re_mat;

				chosen_indices = do(2, 2*size, 2);
				result_eigvals = eigval(collect_matr)[chosen_indices] / sqrt(2); /* eigevalues get doubled --- we select only one from each pair */ 
			end;

		else if beta = 4 then
			do;
				A = RandNormal(size, shape(0, size, 1), I(size));
				B = RandNormal(size, shape(0, size, 1), I(size));
				C = RandNormal(size, shape(0, size, 1), I(size));
				D = RandNormal(size, shape(0, size, 1), I(size));

				Re_diag  =  A + A`;
				Re_ndiag = -C + C;
				Im_diag  =  B - B`;
				Im_ndiag =  D - D`;

				collect_matr = shape(0, 4 * size, 4 * size);

				ind1 = 1:size;
				ind2 = (size+1):(2*size);
				ind3 = (2*size+1):(3*size);
				ind4 = (3*size+1):(4*size);
				/* first column */
				collect_matr[ind1, ind1] = Re_diag;
				collect_matr[ind2, ind1] = Re_ndiag;
				collect_matr[ind3, ind1] = Im_diag;
				collect_matr[ind4, ind1] = Im_ndiag;
				/* second column */
				collect_matr[ind1, ind2] = - Re_ndiag;
				collect_matr[ind2, ind2] = Re_diag;
				collect_matr[ind3, ind2] = Im_ndiag;
				collect_matr[ind4, ind2] = - Im_diag;
				/* third column */
				collect_matr[ind1, ind3] = -Im_diag;
				collect_matr[ind2, ind3] = -Im_ndiag;
				collect_matr[ind3, ind3] = Re_diag;
				collect_matr[ind4, ind3] = Re_ndiag;
				/* third column */
				collect_matr[ind1, ind4] = -Im_ndiag;
				collect_matr[ind2, ind4] = Im_diag;
				collect_matr[ind3, ind4] = -Re_ndiag;
				collect_matr[ind4, ind4] = Re_diag;

				chosen_indices = do(4, 4*size, 4);
				result_eigvals = eigval(collect_matr)[chosen_indices] / sqrt(2); /* eigevalues get quadrupled --- we select only one from each four */
			end;

		else 
			do;
				print "Beta has to be 1, 2 or 4.";
				result_eigvals = .;
			end;

		if normalize then
			do;
				
				result_eigvals = result_eigvals / sqrt( beta * size );
			end;
		
		return result_eigvals;
	finish;

	start _laguerre_sampler_full(M, N=10, beta=2, normalize=1, random_state=1618); /* maybe change N to size */

		call streaminit(random_state);

		result_eigvals = shape(0, N, 1);
		if beta = 1 then
			do;
				A = RandNormal( N, shape(0, M, 1), I(M));
				result_eigvals = eigval(A * A`);
			end;

		else if beta = 2 then
			do;
				A = RandNormal( N, shape(0, M, 1), I(M));
				B = RandNormal( N, shape(0, M, 1), I(M));
				Re_mat = A * A` + B * B`;
				Im_mat = B * A` - A * B`;

				collect_matr = shape(0, 2 * N, 2 * N);

				ind1 = 1:N;
				ind2 = (N+1):(2*N);
				/* first column */
				collect_matr[ind1, ind1] = Re_mat;
				collect_matr[ind2, ind1] = Im_mat;
				/* second column */
				collect_matr[ind1, ind2] = - Im_mat;
				collect_matr[ind2, ind2] = Re_mat;

				chosen_indices = do(2, 2*N, 2);
				result_eigvals = eigval(collect_matr)[chosen_indices]; /* eigevalues get doubled --- we select only one from each pair */
			end;

		else if beta = 4 then
			do;
				A = RandNormal( N, shape(0, M, 1), I(M));
				B = RandNormal( N, shape(0, M, 1), I(M));
				C = RandNormal( N, shape(0, M, 1), I(M));
				D = RandNormal( N, shape(0, M, 1), I(M));

				Re_diag  = A*A` + B*B` + C*C` + D*D`;
				Re_ndiag = A*C` + D*B` - C*A` - B*D`;
				Im_diag  = B*A` - A*B` + D*C` - C*D`;
				Im_ndiag = D*A` + C*B` - B*C` - A*D`;

				collect_matr = shape(0, 4 * N, 4 * N);

				ind1 = 1:N;
				ind2 = (N+1):(2*N);
				ind3 = (2*N+1):(3*N);
				ind4 = (3*N+1):(4*N);
				/* first column */
				collect_matr[ind1, ind1] = Re_diag;
				collect_matr[ind2, ind1] = Re_ndiag;
				collect_matr[ind3, ind1] = Im_diag;
				collect_matr[ind4, ind1] = Im_ndiag;
				/* second column */
				collect_matr[ind1, ind2] = - Re_ndiag;
				collect_matr[ind2, ind2] = Re_diag;
				collect_matr[ind3, ind2] = Im_ndiag;
				collect_matr[ind4, ind2] = - Im_diag;
				/* third column */
				collect_matr[ind1, ind3] = -Im_diag;
				collect_matr[ind2, ind3] = -Im_ndiag;
				collect_matr[ind3, ind3] = Re_diag;
				collect_matr[ind4, ind3] = Re_ndiag;
				/* fourth column */
				collect_matr[ind1, ind4] = -Im_ndiag;
				collect_matr[ind2, ind4] = Im_diag;
				collect_matr[ind3, ind4] = -Re_ndiag;
				collect_matr[ind4, ind4] = Re_diag;
				
				chosen_indices = do(4, 4*N, 4);		
				result_eigvals = eigval(collect_matr)[chosen_indices]; /* eigevalues get quadrupled --- we select only one from each four */
				
			end;

		else 
			do;
				print "Beta has to be 1, 2 or 4.";
				result_eigvals = .;
			end;

		if normalize then
			do;
				result_eigvals = result_eigvals / (beta * M); 
			end;

		return result_eigvals; 
	finish;

	start _mu_ref_gamma_sampler_tridiag(shape=1.0, scale=1.0, beta=2, size=10, normalize=1, random_state=1618);

		call streaminit(random_state);

		result_eigvals = shape(0, size, 1);
		if beta <= 0 then
			do;
        		print "`beta` must be positive.";
				result_eigvals = .;
			end;
		else 
			do;
				xi_odd  = shape(0, size, 1);
				xi_even = shape(0, size, 1);

				do i = 1 to size-1;
					xi_odd[i]    = scale * rand("GAMMA", beta * (size - i)/2 + shape);
					xi_even[i+1] = scale * rand("GAMMA", beta * (size - i)/2);
				end;
				xi_odd[size]   = scale * rand("GAMMA", shape);

				M = diag( xi_odd + xi_even);

				do i = 1 to size-1;
					M[i+1,i] = sqrt(xi_odd[i] * xi_even[i+1]);
					M[i,i+1] = sqrt(xi_odd[i] * xi_even[i+1]);
				end;

				result_eigvals = eigval(M)[,1];
			
				if normalize then
					do; /* potencjalnie cos nie tak */
						result_eigvals = result_eigvals / (0.5 *scale * beta * (2 / beta * shape + size - 1));
					end;	
			end;

		return result_eigvals; 	
	finish;

	start _jacobi_sampler_full(M_1, M_2, N=10, beta=2, normalize=1, random_state=1618); /* normalization does nothing here */
	
		call streaminit(random_state);

		
		if beta = 1 then
			do;
				X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				Y = RandNormal( N, shape(0, M_2, 1), I(M_2));

				X_tmp = X*X`;
				Y_tmp = Y*Y`;

				res_tmp = X_tmp * inv( X_tmp + Y_tmp );
				result_eigvals = eigval(res_tmp)[,1];

				return result_eigvals; 
			end;

		else if beta = 2 then
			do;
				A_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				B_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				A_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));
				B_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));

				A_X_tmp = A_X*A_X` + B_X*B_X`;
				B_X_tmp = B_X*A_X` - A_X*B_X`;
				A_Y_tmp = A_Y*A_Y` + B_Y*B_Y`;
				B_Y_tmp = B_Y*A_Y` - A_Y*B_Y`;

				A_to_solve = shape(0, 2*N, 2*N);
				B_to_solve = shape(0, 2*N, N);


				ind1 = 1:N;
				ind2 = (N+1):(2*N);
				/* first column */
				A_to_solve[ind1, ind1] = A_X_tmp` + A_Y_tmp`;
				A_to_solve[ind2, ind1] = B_X_tmp` + B_Y_tmp`;
				/* second column */
				A_to_solve[ind1, ind2] = -(B_X_tmp` + B_Y_tmp`);
				A_to_solve[ind2, ind2] = A_X_tmp` + A_Y_tmp`;

				B_to_solve[ind1, ind1] = A_X_tmp`;
				B_to_solve[ind2, ind1] = B_X_tmp`;

				quot = solve(A_to_solve, B_to_solve)`;

				collect_mat = shape(0, 2*N, 2*N);
				/* first column */
				collect_mat[ind1, ind1] = quot[ind1, ind1];
				collect_mat[ind2, ind1] = quot[ind1, ind2];
				/* second column */
				collect_mat[ind1, ind2] = -quot[ind1, ind2];
				collect_mat[ind2, ind2] =  quot[ind1, ind1];

				chosen_indices = do(2, 2*N, 2);
				result_eigvals = eigval(collect_mat)[chosen_indices, 1];
				
				return result_eigvals; 
			end;

		else if beta = 4 then
			do;
				A_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				B_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				C_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				D_X = RandNormal( N, shape(0, M_1, 1), I(M_1));
				A_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));
				B_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));
				C_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));
				D_Y = RandNormal( N, shape(0, M_2, 1), I(M_2));

				/* Changes in sign are because of transposition */
				Re_diag_X_t  =  A_X*A_X` + B_X*B_X` + C_X*C_X` + D_X*D_X`;
				Re_ndiag_X_t = -A_X*C_X` - D_X*B_X` + C_X*A_X` + B_X*D_X`;
				Im_diag_X_t  = -B_X*A_X` + A_X*B_X` - D_X*C_X` + C_X*D_X`;
				Im_ndiag_X_t = -D_X*A_X` - C_X*B_X` + B_X*C_X` + A_X*D_X`;

				Re_diag_tot_t  = Re_diag_X_t  + A_Y*A_Y` + B_Y*B_Y` + C_Y*C_Y` + D_Y*D_Y`;
				Re_ndiag_tot_t = Re_ndiag_X_t - A_Y*C_Y` - D_Y*B_Y` + C_Y*A_Y` + B_Y*D_Y`;
				Im_diag_tot_t  = Im_diag_X_t  - B_Y*A_Y` + A_Y*B_Y` - D_Y*C_Y` + C_Y*D_Y`;
				Im_ndiag_tot_t = Im_ndiag_X_t - D_Y*A_Y` - C_Y*B_Y` + B_Y*C_Y` + A_Y*D_Y`;

				A_to_solve = shape(0, 4*N, 4*N);
				B_to_solve = shape(0, 4*N, 2*N);

				ind1 = 1:N; 
				ind2 = (N+1):(2*N);
 				ind3 = (2*N+1):(3*N);
				ind4 = (3*N+1):(4*N);
				/* first column */
				A_to_solve[ind1, ind1] = Re_diag_tot_t;
				A_to_solve[ind2, ind1] = - Re_ndiag_tot_t;
				A_to_solve[ind3, ind1] = Im_diag_tot_t;
				A_to_solve[ind4, ind1] = Im_ndiag_tot_t;
				/* second column */
				A_to_solve[ind1, ind2] =   Re_ndiag_tot_t;
				A_to_solve[ind2, ind2] =   Re_diag_tot_t;
				A_to_solve[ind3, ind2] =   Im_ndiag_tot_t;
				A_to_solve[ind4, ind2] = - Im_diag_tot_t;
				/* third column */
				A_to_solve[ind1, ind3] = - Im_diag_tot_t;
				A_to_solve[ind2, ind3] = - Im_ndiag_tot_t;
				A_to_solve[ind3, ind3] =   Re_diag_tot_t;
				A_to_solve[ind4, ind3] = - Re_ndiag_tot_t;
				/* fourth column */
				A_to_solve[ind1, ind4] = - Im_ndiag_tot_t;
				A_to_solve[ind2, ind4] =   Im_diag_tot_t;
				A_to_solve[ind3, ind4] = - Re_ndiag_tot_t;
				A_to_solve[ind4, ind4] =   Re_diag_tot_t;
	
				/* first column */
				B_to_solve[ind1, ind1] =   Re_diag_X_t;
				B_to_solve[ind2, ind1] = - Re_ndiag_X_t;
				B_to_solve[ind3, ind1] =   Im_diag_X_t;
				B_to_solve[ind4, ind1] =   Im_ndiag_X_t;
				/* second column */
				B_to_solve[ind1, ind2] =   Re_ndiag_X_t;
				B_to_solve[ind2, ind2] =   Re_diag_X_t;
				B_to_solve[ind3, ind2] =   Im_ndiag_X_t;
				B_to_solve[ind4, ind2] = - Im_diag_X_t;

				quot = solve(A_to_solve, B_to_solve);

				collect_mat = shape(0, 4*N, 4*N);
				Lind1 = 1:(2*N);
				Lind2 = (2*N+1):(4*N);


				/* first column */
				collect_mat[Lind1, Lind1] = quot[Lind1, Lind1];
				collect_mat[Lind2, Lind1] = quot[Lind2, Lind1];
				/* second column */
				collect_mat[Lind1, Lind2] = -quot[Lind2, Lind1];
				collect_mat[Lind2, Lind2] = quot[Lind1, Lind1];

				
				chosen_indices = do(4, 4*N, 4);
				result_eigvals = eigval(collect_mat)[chosen_indices, 1];
				
				return result_eigvals; 
			end;

		else 
			do;
				print "Beta has to be 1, 2 or 4.";
				result_eigvals = .;
				return result_eigvals;
			end;

	finish;

	start _mu_ref_beta_sampler_tridiag(a, b, beta=2, size=10, normalize=1, random_state=1618); /* normalization does nothing here */
	
	
		call streaminit(random_state);

		if beta <= 0 then
			do;
        		print "`beta` must be positive.";
				result_eigvals = .;
				return result_eigvals;
			end;
		else 
			do;
				c_odd  = shape(0, size, 1);
				c_even = shape(0, size, 1);
				do i=1 to size-1;
					c_odd[i]    = rand("BETA", beta*(size - i)/2 + a, beta*(size - i)/2 + b);
					c_even[i+1] = rand("BETA", beta*(size - i)/2, beta*(size - i - 1)/2 + b + a);
				end;
				c_odd[size] = rand("BETA", a, b);

				xi_odd = (1 - c_even) # c_odd;
				xi_even = shape(0, size, 1);
				xi_even[2:size] = (1 - c_odd[1:size-1]) # c_even[2:size];

				alpha_coef = xi_even + xi_odd;
				beta_coef = sqrt(xi_odd[1:size-1] # xi_even[2:size]);
	
				M = diag( alpha_coef);
				do i=1 to size-1;
					M[i, i+1] = beta_coef[i];
					M[i+1, i] = beta_coef[i];
				end;

				result_eigvals = eigval(M)[,1];
				return result_eigvals;
			end;
	finish;


	start _circular_sampler_full(N, beta=2, haar_mode='QR', normalize=1, heurestic_fix=1, random_state=1618); /*normalization does nothing here */
    	
		call streaminit(random_state);

		if haar_mode = "Hermite" then
			do;

				if beta = 1 then
					do;
						A = RandNormal( N, shape(0, N, 1), I(N));
						eigvecs = eigvec( A + A` );
						result_eigvals = eigval( eigvecs);
						return result_eigvals;
					end;

				else if beta = 2 then
					do;
						A = RandNormal( N, shape(0, N, 1), I(N));
						B = RandNormal( N, shape(0, N, 1), I(N));

						Re_mat = A + A`;
						Im_mat = B - B`;

						collect_matr = shape(0, 2 * N, 2 * N);
						ind1 = 1:N;
						ind2 = (N+1):(2*N);
						/* first column */
						collect_matr[ind1, ind1] = Re_mat;
						collect_matr[ind2, ind1] = Im_mat;
						/* second column */
						collect_matr[ind1, ind2] = - Im_mat;
						collect_matr[ind2, ind2] = Re_mat;

						eigvecs = eigvec(collect_matr);
						result_eigvals = eigval(eigvecs); 


						if heurestic_fix then
							do;
								eigvals_fixed_heur = shape(0, N, 2);
								probs_for_heurestic = shape(0, N, 1);
								call randgen(probs_for_heurestic, "Uniform");	

								do i=1 to N;
									if probs_for_heurestic[i] < 0.5 then
										eigvals_fixed_heur[i,] = result_eigvals[2*i-1,];
									else
										eigvals_fixed_heur[i,] = result_eigvals[2*i,];
								end;

								return eigvals_fixed_heur;  
							end;
						else
							return result_eigvals;		
					end;

				else if beta = 4 then
					do;
						A = RandNormal( N, shape(0, N, 1), I(N));
						B = RandNormal( N, shape(0, N, 1), I(N));
						C = RandNormal( N, shape(0, N, 1), I(N));
						D = RandNormal( N, shape(0, N, 1), I(N));

						Re_diag  =  A + A`;
						Re_ndiag = -C + C;
						Im_diag  =  B - B`;
						Im_ndiag =  D - D`;

						collect_matr = shape(0, 4 * N, 4 * N);

						ind1 = 1:N;
						ind2 = (N+1):(2*N);
						ind3 = (2*N+1):(3*N);
						ind4 = (3*N+1):(4*N);
						/* first column */
						collect_matr[ind1, ind1] = Re_diag;
						collect_matr[ind2, ind1] = Re_ndiag;
						collect_matr[ind3, ind1] = Im_diag;
						collect_matr[ind4, ind1] = Im_ndiag;
						/* second column */
						collect_matr[ind1, ind2] = - Re_ndiag;
						collect_matr[ind2, ind2] = Re_diag;
						collect_matr[ind3, ind2] = Im_ndiag;
						collect_matr[ind4, ind2] = - Im_diag;
						/* third column */
						collect_matr[ind1, ind3] = -Im_diag;
						collect_matr[ind2, ind3] = -Im_ndiag;
						collect_matr[ind3, ind3] = Re_diag;
						collect_matr[ind4, ind3] = Re_ndiag;
						/* fourth column */
						collect_matr[ind1, ind4] = -Im_ndiag;
						collect_matr[ind2, ind4] = Im_diag;
						collect_matr[ind3, ind4] = -Re_ndiag;
						collect_matr[ind4, ind4] = Re_diag;


						eigvecs = eigvec(collect_matr);	
						result_eigvals = eigval(eigvecs); 

						if heurestic_fix then
							do;
								eigvals_fixed_heur = shape(0, N, 2);
								probs_for_heurestic = shape(0, N, 1);
								call randgen(probs_for_heurestic, "Uniform");	
								
								do i=1 to N;
									if probs_for_heurestic[i] < 0.25 then
										eigvals_fixed_heur[i,] = result_eigvals[2*i-1,];

									else if probs_for_heurestic[i] < 0.5 then
										eigvals_fixed_heur[i,] = result_eigvals[2*i,];

									else if probs_for_heurestic[i] < 0.75 then
										eigvals_fixed_heur[i,] = result_eigvals[4*N - 2*i+1,];

									else
										eigvals_fixed_heur[i,] = result_eigvals[4*N - 2*i+2,];
								end;

								return eigvals_fixed_heur; 
							end;
						else
							do;
								return result_eigvals; 
							end;
					end;

				else 
					do;
						print "For Hermite version Beta has to be 1, 2 or 4.";
						return .;
					end;
			end;

		else if haar_mode = "QR" then
			do;
				if beta = 1 then
					do;
						A = RandNormal( N, shape(0, N, 1), I(N));
						call qr(Q, R, piv, indep, A);
						D = vecdiag(R); /* should there be a transposition here? */
      					U = Q # sign(D); 
						result_eigvals = eigval(U); 
						
						return result_eigvals; 
					end;

				else 
					do;
						print "`QR` method not provided for `beta` ither than 1";
						return .;
					end;
			end;

		else 
			do;
				print "Wrong haar_mode supplied --- only `Hermite` and `QR` are available";
				return .;
			end;
	finish;


	start _mu_ref_unif_unit_circle_sampler(beta=2, size=10, normalize=1, heurestic_fix=1, random_state=1618); /* normalization does nothing here */

	
		call streaminit(random_state);
		if beta = int(beta) && sign(beta) = 1 then
			do;
				alpha_re = shape(0, size, 1);
				alpha_im = shape(0, size, 1);
				do i=1 to size;
					nu = 1 + beta * (size - i);
					gauss_vec = RandNormal(1, shape(0, nu+1, 1), I(nu+1));
					alpha_re[i] = gauss_vec[1] / norm(gauss_vec);
					alpha_im[i] = gauss_vec[2] / norm(gauss_vec);
				end;

				rho = sqrt( 1 - alpha_re[1:size-1]##2 - alpha_im[1:size-1]##2 );

				L = shape(0, 2*size, 2*size);
				M = shape(0, 2*size, 2*size);
				
				M[1,1] = 1;
				M[size+1,size+1] = 1;
				if mod(size, 2) = 0 then 
					do;
						do i=1 to int(size/2);
							/* L */

							/* over diagonal */
							L[2*i-1, 2*i] = rho[2*i-1];
							L[2*i-1 + size, 2*i + size] = rho[2*i-1];
							/* diagonal */
							L[2*i-1, 2*i-1] = alpha_re[2*i-1];
							L[2*i, 2*i] = -alpha_re[2*i-1];
							L[2*i-1 + size, 2*i-1 + size] = alpha_re[2*i-1];
							L[2*i + size, 2*i + size] = -alpha_re[2*i-1];

							L[2*i-1 + size, 2*i-1] = -alpha_im[2*i-1];
							L[2*i + size, 2*i] = -alpha_im[2*i-1];
							L[2*i-1 , 2*i-1 + size] = alpha_im[2*i-1];
							L[2*i, 2*i + size] = alpha_im[2*i-1];
							/* under diagonal */
							L[2*i, 2*i-1] = rho[2*i-1];
							L[2*i + size, 2*i-1 + size] = rho[2*i-1];
						end;
						
						
						do i=1 to int(size/2)-1;
							/* M */

							/* over diagonal */
							M[2*i, 2*i+1] = rho[2*i];
							M[2*i + size, 2*i + size+1] = rho[2*i];
							/* diagonal */
							M[2*i, 2*i] = alpha_re[2*i];
							M[2*i+1, 2*i+1] = -alpha_re[2*i];
							M[2*i + size, 2*i + size] = alpha_re[2*i];
							M[2*i+1 + size, 2*i+1 + size] = -alpha_re[2*i];

							M[2*i+ size, 2*i] = -alpha_im[2*i];
							M[2*i+1+ size, 2*i+1] = -alpha_im[2*i];
							M[2*i, 2*i + size] = alpha_im[2*i];
							M[2*i+1, 2*i+1 + size] = alpha_im[2*i];
							/* under diagonal */
							M[2*i+1, 2*i] = rho[2*i];
							M[2*i+1 + size, 2*i + size] = rho[2*i];
						end;
						M[size,size] = alpha_re[size];
						M[2*size, 2*size] = alpha_re[size];
						M[size, 2*size] = alpha_im[size];
						M[2*size, size] = -alpha_im[size];
					end;
				else
					do;
						/* L */
						L[size, size] = alpha_re[size];
						L[2*size, 2*size] = alpha_re[size];
						L[size, 2*size] = alpha_im[size];
						L[2*size, size] = -alpha_im[size];

						do i=1 to int(size/2);
							/* over diagonal */
							L[2*i-1, 2*i] = rho[2*i-1];
							L[2*i-1 + size, 2*i + size] = rho[2*i-1];
							/* diagonal */
							L[2*i-1, 2*i-1] = alpha_re[2*i-1];
							L[2*i, 2*i] = -alpha_re[2*i-1];
							L[2*i-1 + size, 2*i-1 + size] = alpha_re[2*i-1];
							L[2*i + size, 2*i + size] = -alpha_re[2*i-1];

							L[2*i-1 + size, 2*i-1] = -alpha_im[2*i-1];
							L[2*i + size, 2*i] = -alpha_im[2*i-1];
							L[2*i-1 , 2*i-1 + size] = alpha_im[2*i-1];
							L[2*i, 2*i + size] = alpha_im[2*i-1];
							/* under diagonal */
							L[2*i, 2*i-1] = rho[2*i-1];
							L[2*i + size, 2*i-1 + size] = rho[2*i-1];

							/* M */

							/* over diagonal */
							M[2*i, 2*i+1] = rho[2*i];
							M[2*i + size, 2*i + size+1] = rho[2*i];
							/* diagonal */
							M[2*i, 2*i] = alpha_re[2*i];
							M[2*i+1, 2*i+1] = -alpha_re[2*i];
							M[2*i + size, 2*i + size] = alpha_re[2*i];
							M[2*i+1 + size, 2*i + size+1] = -alpha_re[2*i];

							M[2*i + size, 2*i] = -alpha_im[2*i];
							M[2*i+1 + size, 2*i+1] = -alpha_im[2*i];
							M[2*i , 2*i + size] = alpha_im[2*i];
							M[2*i+1, 2*i + size+1] = alpha_im[2*i];
							/* under diagonal */
							M[2*i+1, 2*i] = rho[2*i];
							M[2*i+1 + size, 2*i + size] = rho[2*i];
						end;
					end;

				result_eigvals = eigval(L*M); 

				if heurestic_fix then
							do;
								eigvals_fixed_heur = shape(0, size, 2);
								probs_for_heurestic = shape(0, size, 1);
								call randgen(probs_for_heurestic, "Uniform");	
								
								

								do i=1 to size;
									if probs_for_heurestic[i] < 0.5 then
										eigvals_fixed_heur[i,] = result_eigvals[2*i-1,];
									else
										eigvals_fixed_heur[i,] = result_eigvals[2*i,];
								end;

								return eigvals_fixed_heur; 
							end;
				else
							do;
								return result_eigvals; 
							end;  
			end;
		else
			do;
				print "`beta` must be positive integer.";
				return .;
			end;	
	
	finish;

	start _ginibre_sampler_full(N, normalize=1, heurestic_fix=1, random_state=1618);
  
	
		call streaminit(random_state);
		A = RandNormal( N, shape(0, N, 1), I(N));
		B = RandNormal( N, shape(0, N, 1), I(N));

    	collect_matrix = shape(0, 2*N, 2*N);
		ind1 = 1:N; ind2 = (N+1):(2*N);
		collect_matrix[ind1, ind1] = A;
		collect_matrix[ind2, ind1] = B;
		collect_matrix[ind1, ind2] = -B;
		collect_matrix[ind2, ind2] = A;

		result_eigvals = eigval(collect_matrix)/sqrt(2);

		if normalize then
			do;
				result_eigvals = result_eigvals / sqrt(N);
			end;

		if heurestic_fix then
			do;
				eigvals_fixed_heur = shape(0, N, 2);
				probs_for_heurestic = shape(0, N, 1);
				call randgen(probs_for_heurestic, "Uniform");	
												
				do i=1 to N;
					if probs_for_heurestic[i] < 0.5 then
						eigvals_fixed_heur[i,] = result_eigvals[2*i-1,];
					else
						eigvals_fixed_heur[i,] = result_eigvals[2*i,];
				end;

				return eigvals_fixed_heur; 
			end;
		else
			do;
				return result_eigvals;
			end;

		
	finish;


start sample_from_beta_ensemble_full( 
	ensemble_version, 
	M_1, M_2, /* M variables only for Laguerre (first) and Jacobi (both) */
	size=10, 
	beta=2, 
	normalize=1, 
	haar_mode="Hermite", /* haar_mode only available for circular ensemble */
	heurestic_fix=1, /* heurestic_fix only available for circular ensemble */
	random_state=1618); 
	
	if ensemble_version = "Hermite" then
		return _hermite_sampler_full(size, beta, normalize, random_state);

	else if ensemble_version = "Laguerre" then
		return _laguerre_sampler_full(M_1, size, beta, normalize, random_state); 

	else if ensemble_version = "Jacobi" then
		return _jacobi_sampler_full(M_1, M_2, size, beta, normalize, random_state); 

	else if ensemble_version = "Circular" then
		return _circular_sampler_full(size, beta, haar_mode, normalize, heurestic_fix, random_state); /*normalization does nothing here */
    
	else if ensemble_version = "Ginibre" then
		return _ginibre_sampler_full(size, normalize, heurestic_fix, random_state);

	else
		do;
			print "Wrong type supplied - only Hermite, Laguerre, Jacobi, and Circular are available.";
			return .;
		end;

finish;

start sample_from_beta_ensemble_banded(
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
	
	if ensemble_version = "Hermite" then 
			return _mu_ref_normal_sampler_tridiag(size, loc, scale, beta, normalize, random_state);
		
	else if ensemble_version = "Laguerre" then
			return _mu_ref_gamma_sampler_tridiag(shape, scale, beta, size, normalize, random_state);

	else if ensemble_version = "Jacobi" then
			return _mu_ref_beta_sampler_tridiag(a, b, beta, size, normalize, random_state); /* normalization does nothing here */
	
	else if ensemble_version = "Circular" then
			return _mu_ref_unif_unit_circle_sampler(beta, size, normalize, heurestic_fix, random_state); /* normalization does nothing here */

	else
 		do;
			print "Wrong type supplied - only Hermite, Laguerre, Jacobi, and Circular are available.";
			return .;
		end;

finish;

	ensemble_version = "Circular";
	size=10;
	beta=200;
	loc=0.0;
	scale=1.0;
	shape = 1.0;
	a = 1.0;
	b = 1.0;
	normalize=0;
	heurestic_fix=0;
	random_state=1618;

	sample = sample_from_beta_ensemble_banded(
	ensemble_version, 
	size, 
	beta, 
	/* distribution specific parameters */
	loc, 
	scale, 
	shape,
	a,
	b,
	/* parameters for the algorithm */
	normalize, 
	heurestic_fix,
	random_state);

	run scatter(sample[,1], sample[,2]);
quit;
