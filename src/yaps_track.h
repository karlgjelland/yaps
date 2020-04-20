	DATA_ARRAY(H);			// Position of hydros 
	DATA_ARRAY(toa);   		// Time of arrival at hydro. One row per buoy, one column per ping
	DATA_INTEGER(nh);
	DATA_INTEGER(np);
	DATA_STRING(pingType);
	DATA_SCALAR(bi_epsilon);    
	DATA_SCALAR(bi_penalty);
	DATA_SCALAR(rbi_min);
	DATA_SCALAR(rbi_max);
	DATA_STRING(ss_data_what);	// 'est' = estimate SS-data; 'data' = Use SS-data
	DATA_VECTOR(ss_data);		// Vector of SS-data if used - length(ss_data) = np
	DATA_IVECTOR(ss_idx);
	DATA_INTEGER(n_ss);
	DATA_SCALAR(approxBI);
	DATA_VECTOR(Edist);
	DATA_VECTOR(biTable);

	DATA_STRING(how_3d);
	DATA_VECTOR(z_vec);
	
	
	PARAMETER_VECTOR(X);	//Position at time of ping
	PARAMETER_VECTOR(Y);	//Position at time of ping
	PARAMETER_VECTOR(top);		// Estimated time of pings
	PARAMETER_VECTOR(ss);		// Estimated speed of sound

	PARAMETER(logD_xy);    		// Diffusivity of fish
	Type D_xy = exp(logD_xy);
	
	PARAMETER(logSigma_bi);		// Sigma for burst interval
	Type sigma_bi = exp(logSigma_bi);

	PARAMETER(logD_v);    		// Diffusivity of sound speed
	Type D_v = exp(logD_v);
	
	PARAMETER(logSigma_toa);	// Sigma TimeOfArrival
	Type sigma_toa = exp(logSigma_toa);

	PARAMETER(logScale);		// scale-parameter for t-dist
	Type scale = exp(logScale);

	PARAMETER(log_t_part);		// t-part of mixture model 
	Type t_part = exp(log_t_part);
	Type G_part = Type(1.0) - t_part; //Gaussian part of mixture model

	PARAMETER_VECTOR(tag_drift);

	array<Type> mu_toa(nh,np);  // mu-matrix
	array<Type> dist(nh,np);	// dist-matrix
	vector<Type> top_pbi(np);

	Type nll = 0.0;

	if(pingType == "pbi"){
		top_pbi(0) = Type(0.0);
		for(int i = 1; i < np; ++i)	{
			top_pbi(i) = top_pbi(i-1) + biTable(i-1);
		}
		for(int i = 0; i < np; ++i)	{
			nll -= dnorm(top(i), top_pbi(i) + tag_drift(i), Type(1E-6), true);
		}
	} else {
		for(int i = 0; i < np; ++i)	{
			nll -= dnorm(tag_drift(i), Type(0), Type(10), true);
		}
	}


	for(int i=0; i<np; ++i) //iterate pings
	{
		for(int h=0; h<nh; ++h){ //iterate hydros
			if(!isNA(toa(h,i))){ //ignore NA's...
				if(how_3d == "none"){
					dist(h,i) = sqrt((H(h,0)-X(i))*(H(h,0)-X(i)) + (H(h,1)-Y(i))*(H(h,1)-Y(i)));
				} else if(how_3d == "data"){
					dist(h,i) = sqrt((H(h,0)-X(i))*(H(h,0)-X(i)) + (H(h,1)-Y(i))*(H(h,1)-Y(i)) + (H(h,2)-z_vec(i))*(H(h,2)-z_vec(i)));
				}
				if(ss_data_what == "est"){
					mu_toa(h,i) = top(i) +  dist(h,i)/ss(ss_idx(i));
				} else {
					mu_toa(h,i) = top(i) +  dist(h,i)/ss_data(i);
				}
				Type eps = toa(h,i) - mu_toa(h,i);
				
				nll -= Edist(0) * dnorm(eps, Type(0), sigma_toa, true); 					//Gaussian part					
				
				nll -= Edist(1) * log( G_part * dnorm(eps, Type(0),sigma_toa,false) + 		//Gaussian part
						t_part * dt(eps/scale, Type(3.0), false) );					//t part
				
				nll -= Edist(2) * log(dt(eps/scale, Type(3.0), false)/scale);

			}
		}
	}
	
	// Needed to ensure positive definite Hessian...
	nll -= dnorm(log_t_part, Type(0), Type(25), true);
	nll -= dnorm(logSigma_toa, Type(0), Type(25), true);
	nll -= dnorm(logSigma_bi, Type(0), Type(25), true);
	nll -= dnorm(logScale, Type(0), Type(25), true);
	nll -= dnorm(logSigma_bi, Type(0), Type(25), true);
	nll -= dnorm(logD_v, Type(0), Type(25), true);

	//position component
	nll -= dnorm(X(0),Type(0),Type(10000),true);
	nll -= dnorm(Y(0),Type(0),Type(10000),true);
	for(int i=1; i<np; ++i)	{
		nll -= dnorm(X(i), X(i-1),sqrt(2*D_xy*(top(i) - top(i-1))),true);	
		nll -= dnorm(Y(i), Y(i-1),sqrt(2*D_xy*(top(i) - top(i-1))),true);
	}
	
	// //speed of sound component
	nll -= dnorm(ss(0),Type(1475.0),Type(40.0),true);		
	for(int i = 1; i < n_ss; ++i){
		nll -= dnorm(ss(i), ss(i-1),sqrt(2*D_v), true);
		nll -= dnorm(ss(i),Type(1475.0),Type(100.0),true);		
	}
	
	//burst interval component
	if(pingType == "sbi"){
		nll -= dnorm(top(0),Type(0.0),Type(4.0),true);
		nll -= dnorm(top(1),Type(approxBI),Type(4.0),true);
		for(int i = 2; i < np; ++i)	{
			nll -= dnorm(top(i)-2*top(i-1)+top(i-2), Type(0),sigma_bi, true);
		}
	} else if(pingType == "rbi"){
		nll -= dnorm(top(0),Type(0.0),Type(4.0),true);
		for(int i = 1; i < np; ++i)	{
			nll -= dnorm(top(i), top(i-1) + (rbi_max - rbi_min)/2, (rbi_max-rbi_min)/2, true);
			nll += bi_penalty * (softplus((top(i) - top(i-1)) - rbi_max, bi_epsilon) + softplus(rbi_min - (top(i) - top(i-1)), bi_epsilon));
		}
	} else if(pingType == "pbi"){
		nll -= dnorm(tag_drift(0), Type(0.0), Type(4), true);
		nll -= dnorm(tag_drift(1), Type(0.0), Type(4), true);
		for(int i = 2; i < np; ++i)	{
			nll -= dnorm(tag_drift(i)-2*tag_drift(i-1)+tag_drift(i-2), Type(0), sigma_bi, true);
		}
	}


	REPORT(mu_toa);
	
	return nll;