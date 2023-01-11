proc iml;
	package load dppsampl;

	sample = poiss_planch_sample(1000, 1618);
	run histogram(sample);
quit;
