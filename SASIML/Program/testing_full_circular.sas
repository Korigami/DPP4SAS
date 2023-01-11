proc iml;
    
	package load dppsampl;
	
	ensemble_version = "Circular";
	size=10;
	beta=4;

	M_1=10; M_2 = 10;
	haar_mode="Hermite";


	normalize=0;
	heurestic_fix=0;
	random_state=1618;

	sample = sample_from_beta_ensemble_full(
	ensemble_version, 
	M_1, M_2, /* M variables only for Laguerre (first) and Jacobi (both) */
	size, 
	beta, 
	normalize, 
	haar_mode, /* haar_mode only available for circular ensemble */
	heurestic_fix, /* heurestic_fix only available for circular ensemble */
	random_state); 

	run scatter(sample[,1], sample[,2]);
quit;
