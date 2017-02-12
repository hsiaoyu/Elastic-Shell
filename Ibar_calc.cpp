// This file generate IAbar from a known shrinkage behavior and precalculated machine direction
#include <iostream>
#include <fstream>
#include <cmath>
#include <igl/per_face_normals.h>
#include <igl/readOBJ.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include "MD_calc.h"
using namespace std;
using namespace Eigen;
MatrixXd V;// ENbar is the edge normal of Vbar
MatrixXi F; // Eash row of EdgeNV(3*NTri,4) stores the index off the four vertices that are related to one edge
int main (){
    igl::readOBJ("V.obj", V, F);
    MatrixXd MachineDirection(F.rows(),2);
    igl::readDMAT("MDinFace.dmat", MachineDirection);
    VectorXd MoistureLevel(F.rows());
    // MoistureLevel stores the moist of face i in element i
    igl::readDMAT("MoistureLevel.dmat", MoistureLevel);
    // V is the rest state of the paper and each faces should be inplane with the corresponding MD
    ifstream InFile("calcIAbar.txt");
    double t, ShrinkCoeffMD, ShrinkCoeffCD;
    InFile >> t;
    InFile >> ShrinkCoeffMD;
    InFile >> ShrinkCoeffCD;
    MatrixXd Itot;
    Itot=CalcIbar(MachineDirection, t , ShrinkCoeffMD, ShrinkCoeffCD, MoistureLevel);
    igl::writeDMAT("Ibar.dmat",Itot,1)
    return 0;
}

//CalcIAbar calculates IAbar and IBbar based on the moisture distribution assuming that only upper surface changes and the lower surface remains intact
MatrixXd CalcIbar(MatrixXd MachineDirection, double t, double ShrinkMD, double ShrinkCD, VectorXd MoistureLevel){
	int NTri = F.rows();
	MatrixXd FN, EN, IAtot(2*NTri,2), IBtot(2*NTri,2), Itot(4*NTri,2);
    	igl::per_face_normals(V, F, FN);
    	EN=Enormal(FN);
	Matrix2d Rot;
	Rot << 0, -1, 1, 0;
	for (i=0; i<NTri; i++){
		Matrix2d IA, IB, tmp1, tmp2, Utransform, MShrink;
		double c1, c2;
		Vector2d MD, CD;
		IB=I2(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)),EN.row(3*i),EN.row(3*i+1),EN.row(3*i+2));
		IA=I1(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)));
		tmp1=IA+t*IB/2;
		tmp2=IA-t*IB/2;
		MD=MachineDirection.row(i).transpose();
		CD=Rot*MD;
		cout << "MD*CD=" << MD*CD << endl;
		c1=ShrinkMD*MD.transpose()*tmp1*MD*MoistureLevel(i);
		c2=ShrinkCD*CD.transpose()*tmp1*CD*MoistureLevel(i);
		Utransform << MD(0), MD(1), CD(0), CD(1); 
		//Utransform *= 1/(MD(0)*CD(1)-MD(1)*CD(0));
		MShrink << c1, 0, 0, c2;
		tmp1= Utransform.transpose()*MShrink*Utranform;
		IAtot.block(2*i,0,2,2)=0.5*(tmp1+tmp2);
		IBtot.block(2*i,0,2,2)=0.5*(tmp1-tmp2)/t;
	}
	Itot << IAtot, IBtot;
	return Itot;
}

Matrix2d I1(Vector3d v1d, Vector3d v2d, Vector3d v3d){
	double q1, q2, q3;
	Matrix2d I;
	q1=(v2d-v1d).squaredNorm();
	q2=(v3d-v2d).squaredNorm();
	q3=(v1d-v3d).squaredNorm();
	I=QMatrix(q1,q2,q3);
	return I;	
}

Matrix2d I2(Vector3d v1, Vector3d v2, Vector3d v3, Vector3d n1, Vector3d n2, Vector3d n3){
	Vector3d e1, e2, e3;
	double q1, q2, q3;
	Matrix2d I;
	e1=v2-v1;
	e2=v3-v2;
	e3=v1-v3;
	q1=2*e1.dot(n3-n2);
	q2=2*e2.dot(n1-n3);
	q3=2*e3.dot(n2-n1);
	I=QMatrix(q1,q2,q3);
	return I;
}

MatrixXd Enormal(MatrixXd FNf){	// FN contains the face normals and EdgeF contains the neighboring faces of each edge
	MatrixXd EdgeN(EdgeF.size(),3);
	int i;
	for(i=0; i<EdgeF.size(); i++){
		if (EdgeF(i)==-1){
			EdgeN.row(i) << FNf.row(i/3);
		}
		else {
			EdgeN.row(i)=FNf.row(i/3)+FNf.row(EdgeF(i));
			EdgeN.row(i)=EdgeN.row(i).normalized();
		}	
	}
	return EdgeN;
} 
