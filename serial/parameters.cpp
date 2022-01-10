#include "parameters.h"

using namespace std;

Parameters::Parameters()
{
        this->lambda = 0;
        this->potId = 0;
        this->nGrid = 0;
        this->dr = 0;
        this->na = 0;
        this->da = 0;
        this->nw = 0;
        this->nSq = 0;
        this->temp = 0;
        this->tau = 0;
        this->nDim = 0;
        this->printPw = false;
        this->printRn = 1;
}

Parameters::~Parameters()
{}


void Parameters::checkValues() {
        if(this->lambda == 0) {
            this->lambda = 1.0;
            cout << "lambda value has not been provided, using default value: " << this->lambda << endl;
        }

        if(this->potId == 0) {
            cout << "Free particle potential potId = 0" << endl;
        }
        else if (this->potId < 0 || this->potId > 3) {
            cout << "potId not implemented (0 through 3 only)! Using free praticle expression" << endl;
            this->potId = 0;
        }

        if(this->nGrid == 0) {
            this->nGrid = 100;
            cout << "nGrid value has not been provided, using default value: " << this->nGrid << endl;
        }

        if(this->dr == 0) {
            this->dr = 0.01;
            cout << "dr value has not been provided, using default value: " << this->dr << endl;
        }

        if(this->na == 0) {
            this->na = 10;
            this->da = PI / double(this->na);
            cout << "na value has not been provided, using default value: " << this->na << endl;
        }

        if(this->nw == 0) {
            this->nw = 1;
            cout << "nw value has not been provided, using default value: " << this->nw << endl;
        }

        if(this->nSq == 0) {
            cout << "nSq value has not been provided, partial waves will not be squared" << endl;
        }

        if(this->temp == 0) {
            this->temp = 10.0;
            this->tau = 1.0 / this->temp;
            cout << "tHigh value has not been provided, using default value: " << this->temp << endl;
        }

        if(this->nDim == 0) {
            this->nDim = 2;
            cout << "nDim value has not been provided, using default value: " << this->nDim << endl;
        }
}


int Parameters::readInput() {
    ifstream input;
	string fname;
	fname="pdm.inp";
	string command;
	
	input.open(fname.c_str(), ios::in);
	if (input.fail())
	{
		cout << "Input file pdm.inp could not be found!" << endl;
		return -1;
	}

	string line;
	string temp;
	while(getline(input,line))
	{
		istringstream row(line);
		row >> command;
		if( command=="PHYSI" )
		{
			if(not (row >> this->lambda >> this->nDim >> this->potId)) {
				cout << command << " has missing values in the input file pdm.inp" << endl;
                return -1;
            }
			else if (row >> temp) {
				cout << command << " has too many values in the input file pdm.inp" << endl;
                return -1;
            }
		}
		else if ( command=="RGRID" )
		{
			if(not (row >> this->nGrid >> this->dr)) {
				cout << command << " has missing values in the input file pdm.inp" << endl;
                return -1;
            }
			else if (row >> temp) {
				cout << command << " has too many values in the input file pdm.inp" << endl;
                return -1;
            }
		}
		else if ( command=="AGRID" )
		{
			if(not (row >> this->na)) {
				cout << command << " has missing values in the input file pdm.inp" << endl;
                return -1;
            }
			else if (row >> temp) {
				cout << command << " has too many values in the input file pdm.inp" << endl;
                return -1;
            }
			this->da = PI / double(this->na);
		}
		else if ( command=="NWAVE" )
		{
			if(not (row >> this->nw)) {
				cout << command << " has missing values in the input file pdm.inp" << endl;
                return -1;
            }
			else if (row >> temp) {
				cout << command << " has too many values in the input file pdm.inp" << endl;
                return -1;
            }
		}
		else if ( command=="SQUAR")
		{
			if(not (row >> this->nSq >> this->temp)) {
				cout << command << " has missing values in the input file pdm.inp" << endl;
                return -1;
            }
			else if (row >> temp) {
				cout << command << " has too many values in the input file pdm.inp" << endl;
                return -1;
            }
			this->tau = 1. / this->temp;
		}
		else if ( command=="PRWAV" )
		{	
			this->printPw = true;
		}
		else if ( command=="PRNOR")
		{
			if(not (row >> this->printRn)) {
				cout << command << " has a mssing argument in pdm.inp" << endl;
                return -1;
            }
			else if(row >> temp) {
				cout << command << " has too many arguments in pdm.inp" << endl;
                return -1;
            }
		}
		else
		{
			cout << command << " not an existing input command! Check pdm.inp" << endl;
			return -1;
		}
	}
	input.close();

    // now check to see if everything that is needed has been provided
    // if not, use defaults
    this->checkValues();
    cout << flush;

	return 0;
}