#include <vector>
#include <stdio.h>
#include <iostream>

using namespace std;

int main() {
	vector <double> zero1d;
	vector <vector <double>> zero2d;
//	int NW=1;
//	int NGRID=3;
	int i,j;
/*
	for(i=0;i<NW; i++)
		zero1d.push_back(0.);
	cout << "zero1d size: " << zero1d.size() << endl; 
	for(i=0;i<NGRID;i++)
		zero2d.push_back(zero1d);
	cout << "zero2d size: " << zero2d.size() << endl;

	for(i=0; i<NGRID; i++) {
		for(j=0; j<NGRID; j++) {
			cout << i << "\t" << j << "\t" << zero2d[i][j] << endl;
		}
	}

*/
/*
	int R = 3;
	int C = 2;
	vector<int> v(C);
	vector<vector<int>> mat(R,v);

	for(i=0; i<5; i++) {
		for(j=0; j<C; j++) {
			cout << i << "\t" << j << "\t" << mat[i][j] << endl;
		}
	}
*/
/*
	int l;
	int NW=2;
	int NGRID=3;
	vector<int> v(NW,0);
	vector<vector<int>> m(NGRID,v);
	vector<vector<vector<int>>> t(NGRID,m);

	for(i=0; i<NGRID; i++) {
		for(j=0; j<NGRID; j++) {
			for(l=0; l<5; l++) {
				cout << i << "\t" << j << "\t" << l << "\t" << t[i][j][l] << endl;
			}
		}
	}
*/

	typedef vector<int> v1d;
	typedef vector<v1d> v2d;
	typedef vector<v2d> v3d;
	int NGRID=3,NW=2;
	int l;
	v3d t(NGRID, v2d(NGRID, v1d(NW,0)));
	for(i=0; i<NGRID; i++) {
		for(j=0; j<NGRID; j++) {
			for(l=0; l<NW; l++) {
				cout << i << "\t" << j << "\t" << l << "\t" << t[i][j][l] << endl;
			}
		}
	}
	v1d vec(2,0);
	//vec.assign(2,0);
	vec[7] = 0;
	
/*
	vector<int> ve(5);
	for(i=0; i<8; i++)
		cout << ve[i] << endl;
*/
/*
	int jim[NGRID][NGRID][NW] = {0};
	for(i=0;i<NGRID;i++) {
		for(j=0;j<NGRID;j++) {
			for(l=0;l<NW;l++)
				cout << jim[i][j][l] << endl;
		}
	}

*/

	
	return 0;
}
