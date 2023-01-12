%include 'path to random_matrices_fixes.sas';
%sample_from_beta_ensemble_banded(
	result_eigvals= as, 
	ensemble_version =Hermite, 
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

	%mend;
