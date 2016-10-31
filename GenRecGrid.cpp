// This script generate a mxn rectangular lattice filled with triangular mesh 
#include <iostream>
#include <fstream>
#include <cmath>
#include <igl/readOBJ.h>
#include <igl/viewer/Viewer.h>
using namespace std;
using namespace Eigen;

int main() {
	int m, n, NSq;
	double space;
	
	cout << "This script will generate a mxn rectangular lattice" << endl;
	cout << "Please enter m:" << endl;
        cin >> m;
	cout << "Please enter n:" << endl;
        cin >> n;
	cout << "Please enter spacing:" << endl;
        cin >> space;
	
	NSq=(m-1)*(n-1);
	MatrixXi F(NSq*2,3);
	MatrixXd V(m*n,3);
	for (int i=0; i<m ; i++){
		for (int j=0; j<n; j++){
			V.row(i*n+j) << (j%n)*space, i*space, 0;
		}
	}
	for (int i=0;i<m-1;i++){
		for (int j=1; j<n; j++){
			int flag;
			flag=(i+j)%2;
			F.row(i*(n-1)+j-1) << i*n+j-1, i*n+j, i*n+j+n+flag-1;	
			if (flag==1){
				F.row(i*(n-1)+j-1+NSq) << i*n+j-1, i*n+j+n, i*n+j+n-1;
			}
			else {
				F.row(i*(n-1)+j-1+NSq) << i*n+j, i*n+j+n, i*n+j+n-1; 
			}
		}
	}
//	cout << V << endl <<"F"<< endl << F << endl;
	igl::writeOBJ("RecLattice.obj",V,F);
	return 0;
}
