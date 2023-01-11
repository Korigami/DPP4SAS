
proc iml;
    
	package load dppsampl;
	
	ensemble_version = "Jacobi";
	size=1000;
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

	run histogram(sample);
quit;
