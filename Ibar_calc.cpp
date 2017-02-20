// This file generate IAbar from a known shrinkage behavior and precalculated machine direction
#include <iostream>
#include <fstream>
#include <cmath>
#include <igl/per_face_normals.h>
#include <igl/readOBJ.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include "Ibar_calc.h"
using namespace std;
using namespace Eigen;
MatrixXd V;// ENbar is the edge normal of Vbar
MatrixXi F; // Eash row of EdgeNV(3*NTri,4) stores the index off the four vertices that are related to one edge
VectorXi EdgeF;
Matrix2d CanonM1, CanonM2, CanonM3;
double area=0.5;

int main (){
    MatrixXd GlobalMD;
    // Take special care with the GlobalMD.dmat format! GlobalMD should be inplane with each faces.
    igl::readOBJ("Vbar.obj", V, F);
    ifstream Infile("Input_Filename.txt");
    string str;
    getline(Infile,str); 
    const char *cstr = str.c_str();
    igl::readDMAT(cstr, GlobalMD);
    MatrixXd MachineDirection(2*F.rows(),2);
    MachineDirection = CalcMD(GlobalMD);
    //---------------------------------------------
    VectorXd MoistureLevel(F.rows());
    // MoistureLevel stores the moist of face i in element i
    getline(Infile,str); 
    const char *cstr1 = str.c_str();
    igl::readDMAT(cstr1, MoistureLevel);
    // V is the rest state of the paper and each faces should be inplane with the corresponding MD
    getline(Infile,str); 
    const char *cstr2 = str.c_str();
    ifstream InFile(cstr2);
    double t, ShrinkCoeffMD, ShrinkCoeffCD;
    InFile >> t;
    InFile >> ShrinkCoeffMD;
    InFile >> ShrinkCoeffCD;
    CanonM3 << 1, 0, 0, 0;
    CanonM2 << 1, 1, 1 ,1;
    CanonM1 << 0, 0, 0, 1;
    EdgeF=NeighborF();
    MatrixXd Itot;
    Itot=CalcIbar(MachineDirection, t , ShrinkCoeffMD, ShrinkCoeffCD, MoistureLevel);
    getline(Infile,str); 
    const char *cstr3 = str.c_str();
    igl::writeDMAT(cstr3,Itot,1);
   
/* MatrixXd Icheck;
    igl::readDMAT("Ibar.dmat", Icheck);
    cout << Itot << endl << "Check" << endl << Icheck << endl;
*/
    return 0;
}

//CalcIAbar calculates IAbar and IBbar based on the moisture distribution assuming that only upper surface changes and the lower surface remains intact
MatrixXd CalcIbar(MatrixXd MachineDirection, double t, double ShrinkMD, double ShrinkCD, VectorXd MoistureLevel){
	int NTri = F.rows();
	MatrixXd FN, EN, IAtot(2*NTri,2), IBtot(2*NTri,2), Itot(4*NTri,2);
    	igl::per_face_normals(V, F, FN);
    	EN=Enormal(FN);
	for (int i=0; i<NTri; i++){
		Matrix2d IA, IB, tmp1, tmp2, U, Utransform, MShrink;
		double c1, c2, c3, ctemp;
		Vector2d MD, CD;
		IB=I2(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)),EN.row(3*i),EN.row(3*i+1),EN.row(3*i+2));
		IA=I1(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)));
		tmp1=IA+t*IB;
		tmp2=IA-t*IB;
		//Canonical e1=(1,0), e2=(-1,1)
		MD<<MachineDirection(i,0)-MachineDirection(i,1),MachineDirection(i,1);
		CD<<MachineDirection(i+NTri,0)-MachineDirection(i+NTri,1),MachineDirection(i+NTri,1);
		ctemp=MD.transpose()*tmp1*MD;
		c1=(1.0-ShrinkMD*MoistureLevel(i))*(1.0-ShrinkMD*MoistureLevel(i))*ctemp;
		ctemp=CD.transpose()*tmp1*MD;
		c2=(1.0-ShrinkCD*MoistureLevel(i))*(1.0-ShrinkMD*MoistureLevel(i))*ctemp;
		ctemp=CD.transpose()*tmp1*CD;
		c3=(1.0-ShrinkCD*MoistureLevel(i))*(1.0-ShrinkCD*MoistureLevel(i))*ctemp;
		//c2=ShrinkCD*CD.transpose()*tmp1*CD*MoistureLevel(i);
		U << MD(0), CD(0), MD(1), CD(1);
		Utransform=U.inverse(); 
		//Utransform *= 1/(MD(0)*CD(1)-MD(1)*CD(0));
		MShrink << c1, c2, c2, c3;
		tmp1= Utransform.transpose()*MShrink*Utransform;
		IAtot.block(2*i,0,2,2)=0.5*(tmp1+tmp2);
		IBtot.block(2*i,0,2,2)=0.5*(tmp1-tmp2)/t;
	}
	Itot << IAtot, IBtot;
	return Itot;
}

//CalcMD calculates GlobalMD in barycentric coordinate in each face
MatrixXd CalcMD(MatrixXd GlobalMD){
	Vector3d e1, e2;
	Matrix2d A;
	MatrixXd MD(2*F.rows(),2);
        Vector2d GMD,GCD;
//	GMD is defined as follow if MD is uniform throughout the paper
	GMD << GlobalMD(0,0), GlobalMD(0,1);
        Matrix2d Rot;
	Rot << 0,-1,1,0;
        GCD = Rot*GMD;
	for (int i=0 ; i<F.rows(); i++){
        /*	// If MD is not uniform throughout the paper use the below definition for GMD
		GMD << GlobalMD(i,0), GlobalMD(i,1);
	*/
		e1 = V.row(F(i,1))-V.row(F(i,0));
        	e2 = V.row(F(i,2))-V.row(F(i,1));
		A << e2(1), -e2(0), -e1(1), e1(0);
		A *= 1/(e1(0)*e2(1)-e1(1)*e2(0));
		Vector2d tmp1, tmp2;
		tmp1 = A*GMD;
                tmp2 = A*GCD;
		// MD is normalized so that IAbar would be easier to calculate by unitary transform
		MD.row(i) = tmp1.transpose()/tmp1.norm();
		MD.row(i+F.rows()) = tmp2.transpose()/tmp2.norm();
//		cout << MD(i,0)*e1+MD(i,1)*e2 << endl << "GMD" << endl <<  GMD << " " << GlobalMD(0,2) << endl;
	}
	return MD;
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

VectorXi NeighborF( ){
	int NTri=F.rows(),i,j,m,n,flag;
	VectorXi EdgeFf(3*NTri);	// EdgeF records the neighboring face of an edge in the tri. mesh in ccw order, #n edge in  face #i is Edge(3*i+n) 
	EdgeFf=-1*EdgeFf.setOnes(3*NTri);
	for (i=0; i< NTri; i++){
		for (j=0; j<3; j++){
			m=0;
			flag=0;
			while (m<NTri && flag<2){
				flag=0;
				if (m!=i){	//avoid self counting
					for(n=0; n<3; n++){
						if (F(i,j)==F(m,n) || F(i,(j+1)%3)==F(m,n)){
							flag++;
						}
					}
				}
				m++;
			}
		
			if (flag==2){
				EdgeFf(3*i+j)=m-1;
			}
		}
	}
	return EdgeFf;
}

Matrix2d QMatrix(double q1, double q2, double q3){
	Matrix2d Q;
	double c;
	c=-1/(8*area*area);
	Q=c*((q1-q2-q3)*CanonM1+(q2-q3-q1)*CanonM2+(q3-q1-q2)*CanonM3);
	return Q;
}
