

using namespace std;
using namespace std::chrono;




// functions to run the squaring procedure: implementations at the end
void read_parameters();  // reads input parameters from file pdm.inp
bool check_global_variables(); // checks if there's nothing crazy around
double asympbessel_pf(double n, double x); // computes the pre-factor of asymptotic Bessel functions
double pot(double r); //this is where the physics of your problem comes in
double interaction(double r1, double r2); // numerical integration of the potential
double free_pw(double r1, double r2, int l, double tau, int dim); // calculates expression for free partial waves
double correction_integral(double r1, double r2, double tau); // calculates correction for truncating point, see documentation

// HERE WE GO
int main()
{
	auto start = high_resolution_clock::now();		

	read_parameters();
	bool go;
	go = check_global_variables();
	if(!go)
	{
		return -1;	
	}

	int h, i, j, k, l;	
	// initialize vectors!
	vector<double> tau;	//time slice at each squaring
	vector<double> fc;	//pre-factor on gaussian propagator at each squaring
	vector<double> gc;	//factor inside gaussian at each squaring
	vector<double> twl;	//thermal wavelength
	for(i=0; i <= NSQ_; i++)
	{
		tau.push_back(0.);
		fc.push_back(0.);
		gc.push_back(0.);
		twl.push_back(0.);
	}
	for(i=0; i <= NSQ_; i++)
	{
		tau[i] = TAU_MAX_ * pow(2.0, double(i));
		fc[i] = pow(4. * PI * LAMBDA_ * tau[i], -1.5);
		gc[i] = 1. / (4. * PI * LAMBDA_ * tau[i]);
		twl[i] = 2. * sqrt(2. * LAMBDA_ * tau[i]);
	}

	// initialize grids
	// spatial first
	vector<double> grid;
	for(i=1; i<=NGRID_; i++)		// notice that here we could start with 0, but that runs into trouble for 3D
	{
		grid.push_back(i * DR_);
	}
	// angular grid
	vector<double> agrid;
	for(i=0; i<NA_; i++)
		agrid.push_back(i * DA_);

	
	// now initialize a matrix that will receive the partial wave summation
	int nad;
	int nau;
	nad = int(1. / DA_) - NA_;
	nau = int(1. / DA_) - 1;
	vector< vector < vector<double> > > pws;
	vector< double > zero1d;
	vector< vector < double > > zero2d;
	for(i=0; i<NA_; i++)
		zero1d.push_back(0.);
	for(i=0; i<NGRID_; i++)
		zero2d.push_back(zero1d);
	for(i=0; i<NGRID_; i++)
		pws.push_back(zero2d);

	//initialize vector for partial waves
	vector< vector< vector< double > > > pw;
	vector< vector< vector< double > > > pw2;
	zero1d.clear();
	for(i=0;i<NGRID_;i++)
		zero1d.push_back(0.);
	zero2d.clear();
	for(i=0;i<NGRID_;i++)
		zero2d.push_back(zero1d);
	for(i=0; i<NW_; i++) {
		pw.push_back(zero2d);
		pw2.push_back(zero2d);
	}

	// define some variables that will be useful
	double r1, r2;		// radial distances between the pair
	double r12;		// r1 * r2;
	double dm;		// will receive the density matrix value
	double dmf;		// will receive the free density matrix
	double arg_exp;		//will receive the function
	double it;

	// start high temperature expression, free-particle Ceperley1995 Eq. 4.42
	for(k=0; k<NW_; k++)
	{
		for(i=0; i<NGRID_; i++)
		{
			for(j=0; j<NGRID_; j++)
			{
				r1 = grid[i];
				r2 = grid[j];
				pw[k][i][j] = free_pw(r1,r2,k,tau[0],NDIM_);
			}
		}
	}
	// calculate interactions separately since it does not depend on the order of the partial wave
	// this saves quite a bit of time
	for(i=0; i<NGRID_; i++)
	{
		for(j=0; j<NGRID_; j++)
		{
			r1 = grid[i];
			r2 = grid[j];
			arg_exp = -tau[0] * interaction(r1,r2);
			it = exp(arg_exp);
			for(k=0; k<NW_; k++)
			{
				pw[k][i][j] = pw[k][i][j] * it;
			}
		}
	}


	// SQUARE IT NOW
	double qp;
	double quad;
	quad = 0.;
	int ncut;
	int nsml;
	double reff;
	double c1,c2;

	for(h=1; h<=NSQ_; h++)							// loop over squaring
	{
		ncut = NGRID_ - int(twl[h] / DR_);
		nsml = int(twl[h] / DR_) + 1;

		for(l=0; l<NW_; l++)						// loop over waves
		{
			for(i=0; i<NGRID_; i++)					// loop over spatial grid 1
			{
				for(j=0; j<NGRID_; j++)				// loop over spatial grid 2
				{
					r1 = grid[i];
					r2 = grid[j];
					for(k=0; k<NGRID_; k++)			// loop over spatial grid 3 (here's the actual matrix multiplication)
					{
						qp = pw[l][i][k] * pw[l][j][k] * DR_;
						quad = quad + qp;
					}


/*
					// include correction RGStorer1968 eq.3.9
					c1 = pw[i][NGRID_-1][l]*pw[NGRID_-1][j][l];
					c2 = free_pw(r1,NGRID_*DR_,l,tau[h],NDIM_) * free_pw(NGRID_*DR_,r2,l,tau[h],NDIM_);
					if (c1*c2==0)
						dm = 0;
					else
					{
						dm = c1/c2;
						dm = dm * exp(-tau[h-1]*(pot(r1)+pot(r2)));
						dm = dm * correction_integral(r1,r2,tau[h]);
					}
					quad = quad + dm;
*/
			
					pw2[l][i][j] = quad;
					quad = 0.;
				}
			}


			// CORRECTION NEAR TRUNCATION POINT
			// presumably potential->0 
			for(i=ncut; i<NGRID_; i++)
			{
				for(j=ncut; j<NGRID_; j++)
				{
					r1 = grid[i];
					r2 = grid[j];
					dm = free_pw(r1,r2,l,tau[h],NDIM_);
					//interaction
					dm = dm * exp(-tau[h-1]*(pot(r1) + pot(r2)));  // primitive approximation
					pw2[l][i][j] = dm;
				}
			}



			//reset values for next convolution
			for(i=0; i<NGRID_; i++)
			{	
				for(j=0; j<NGRID_; j++)
				{
					pw[l][i][j] = pw2[l][i][j];
					pw2[l][i][j] = 0.;

				}
			}


		}
	}


	// print waves
	ofstream partwave;
	string fname;
	if(PRINT_PW_)
	{
		for(k=0; k<NW_; k++)
		{
			fname = "pwl" + to_string(k) + ".csv";
			partwave.open(fname.c_str(), ios::out);
			for(i=0; i<NGRID_/RN_; i++)
			{
				for(j=0; j<NGRID_/RN_; j++)
				{
					partwave << grid[i*RN_] << "," << grid[j*RN_] << "," << pw[k][i*RN_][j*RN_] << endl;
				}
			}
			partwave.close();
		}
	}




	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	// NOW PERFORM PARTIAL WAVE SUM
	double darg;
	for(k=0; k<NA_; k++)
	{
		darg = agrid[k];
		for(i=1; i<NGRID_; i++)
		{	
			for(j=1; j<NGRID_; j++)
			{
				r12 = grid[i]*grid[j];

				if(NDIM_==2)
				{
					dm = 1. / (2.*PI*sqrt(r12));
					pws[i][j][k] = pws[i][j][k] + pw[0][i][j] * dm;
					for(l=1; l<NW_; l++)
					{
						pws[i][j][k] = pws[i][j][k] + dm * 2. * pw[l][i][j] * cos(l*darg);
					}
				}
				else if(NDIM_==3)
				{
					dm = 1. / (4. * PI * r12);
					pws[i][j][k] = pws[i][j][k] + pw[0][i][j] * dm;
					for(l=1; l<NW_; l++)
					{
						pws[i][j][k] = pws[i][j][k] + dm * (2.*l + 1)*pw[l][i][j] * gsl_sf_legendre_Pl(l, cos(darg) );
					}

				} 
			}
		}
	}

	// PRINT IT
	ofstream pdm;
	fname = "pdm.csv";
	pdm.open(fname.c_str(), ios::out);
	for(i=0;i<NGRID_/RN_;i++)
	{
		for(j=0;j<NGRID_/RN_;j++)
		{
			for(k=0;k<NA_;k++)
			{
				pdm << grid[i*RN_] << "," << grid[j*RN_] << "," << agrid[k] << "," << pws[i*RN_][j*RN_][k] << endl;
			}
		}
	}
	pdm.close();

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);

	cout << "Calculation completed successfully. \n Execution time was " << duration.count() << " s. \nSee you later!" << endl;
	return 0;
}








double correction_integral(double r1, double r2, double tau)   // calculates correction for truncating point, see documentation
{
	double c;
	double arg_erf;
	c = 4. * sqrt(2.* PI * LAMBDA_ * tau);
	c = 1./c;
	c = c * exp( (-3.*r1*r1 -3.*r2*r2 + 2.*r1*r2) / (16.*LAMBDA_*tau) );
	arg_erf = 0.5*(r1+r2) - NGRID_*DR_;
	arg_erf = arg_erf / sqrt(2.*LAMBDA_*tau);
	c = c * (1 - gsl_sf_erf(arg_erf));
	return c;
}
