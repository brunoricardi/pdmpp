#include "parameters.h"

using namespace std;

parameters::parameters()
{}

parameters::~parameters()
{}

int parameters::readInput() {
    ifstream input;
	string fname;
	fname="pdm.inp";
	string command;
	double read_d;
	int read_i;
	
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
                return -1
            }
			else if (row >> temp) {
				cout << command << " has too many values in the input file pdm.inp" << endl;
                return -1;
            }
		}
		else if ( command=="SQUAR")
		{
			if(not (row >> this->nSq >> this->tHigh)) {
				cout << command << " has missing values in the input file pdm.inp" << endl;
                return -1;
            }
			else if (row >> temp) {
				cout << command << " has too many values in the input file pdm.inp" << endl;
                return -1;
            }
			this->tauMax = 1. / tHigh;
		}
		else if ( command=="PRWAV" )
		{	
			this->printPw = true;
		}
		else if ( command=="PRNOR")
		{
			if(not (row >> this->rn)) {
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
	return 0;
}