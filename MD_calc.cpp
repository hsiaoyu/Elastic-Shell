// This file generate the machine direction in barycentric coordinate of each face and store in MDinFace.dmat
// Read out MDinFace.dmat using readDMAT and each row of the matrix is the barycentric coordinate representation of MD 
#include <iostream>
#include <fstream>
#include <cmath>
#include <igl/readOBJ.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include "MD_calc.h"
using namespace std;
using namespace Eigen;
MatrixXd V;// ENbar is the edge normal of Vbar
MatrixXi F; // Eash row of EdgeNV(3*NTri,4) stores the index off the four vertices that are related to one edge
int main (){
    MatrixXd GlobalMD;
    // Take special care with the GlobalMD.dmat format!
    igl::readDMAT("GlobalMD.dmat", GlobalMD);
    // V is the rest state of the paper and each faces should be inplane with the corresponding MD
    igl::readOBJ("V.obj", V, F);
    MatrixXd MachineDirection(F.rows(),2);
    MachineDirection = CalcMD(GlobalMD);
    igl::writeDMAT("MDinFace.dmat",MachineDirection,1);
   // ofstream OutFile("MDinFace.dmat");
   // OutFile << 2 << " " << F.rows() << endl << MachineDirection.transpose() << endl;
   /* OutFile.close();
    cout << MachineDirection << endl;
    MatrixXd tmp;
    igl::readDMAT("MDinFace.dmat",tmp);
    cout << tmp << endl;
   */return 0;
}

//CalcMD calculates GlobalMD in barycentric coordinate in each face
MatrixXd CalcMD(MatrixXd GlobalMD){
	Vector3d e1, e2;
	Matrix2d A;
	MatrixXd MD(F.rows(),2);
        Vector2d GMD;
//	GMD is defined as follow if MD is uniform throughout the paper
	GMD << GlobalMD(0,0), GlobalMD(0,1);
	for (int i=0 ; i<F.rows(); i++){
        /*	// If MD is not uniform throughout the paper use the below definition for GMD
		GMD << GlobalMD(i,0), GlobalMD(i,1);
	*/
		e1 = V.row(F(i,1))-V.row(F(i,0));
        	e2 = V.row(F(i,2))-V.row(F(i,1));
		A << e2(1), -e2(0), -e1(1), e1(0);
		A *= 1/(e1(0)*e2(1)-e1(1)*e2(0));
		Vector2d tmp;
		tmp = A*GMD;
		// MD is normalized so that IAbar would be easier to calculate by unitary transform
		MD.row(i) = tmp.transpose()/tmp.norm();
//		cout << MD(i,0)*e1+MD(i,1)*e2 << endl << "GMD" << endl <<  GMD << " " << GlobalMD(0,2) << endl;
	}
	return MD;
}

