// This file could calculate energy correctly!   9/12/2016
#include <iostream>
#include <fstream>
#include <cmath>
#include <igl/per_face_normals.h>
#include <igl/readOBJ.h>
#include <igl/viewer/Viewer.h>
#include <igl/per_edge_normals.h>
#include "func.h"
using namespace std;
using namespace Eigen;

int main() {
	int i, j, NTri,NNode,NLen,NWid;
	double L, midL, midW, r, Energy1, constE, nu, t;
	Vector3d E(0,0,0), tmpE;
	constE=1.7242e7;
	nu=0.3;
	t=0.01;
	Vector2d newV;
	Vector3d Min, Max;
	MatrixXd Vbar;
	MatrixXi F;
	igl::readOBJ("rectangle_slit.obj", Vbar, F);
	
	//igl::viewer::Viewer viewer;
	//viewer.data.set_mesh(V, F);
	NTri=F.rows();
	NNode=Vbar.rows();	
	MatrixXd V0(NNode,2), V(NNode,3),Epnt(3*NTri,3), Epnt2(3*NTri,3), EN(3*NTri,3), FN(NTri,3), ENbar(3*NTri,3), FNbar(NTri,3);
	VectorXi EdgeF(3*NTri);
	EdgeF=NeighborF(F);
	Min=Vbar.colwise().minCoeff();
	Max=Vbar.colwise().maxCoeff();
	L=(Max(0)-Min(0));
	midL=(Max(0)+Min(0))/2;
	midW=(Max(1)+Min(1))/2;
	r=L/M_PI;
/*	cout << "L is " << L << " r is " << r << " mid is " << midL << endl;
	newV=FoldCircle(midL, midL, r);
	cout << newV << endl;
	newV=FoldCircle(midL, midL-(L/2), r);
	cout << newV << endl;
*/
	for (i=0; i<NNode; i++){
		V0.row(i) << Vbar(i,0), Vbar(i,1);
		V.row(i) << Vbar(i,0), (Vbar(i,1)-midW)*1.2 , Vbar(i,2);
	/*	// Bend into a cylinder
		newV = FoldCircle(midL, Vbar(i,0), r);
		V.row(i)  << newV(0), Vbar(i,1), newV(1);
	*/
	}
	
//	igl::viewer::Viewer viewer;
//	viewer.data.set_mesh(V, F);	
	igl::per_face_normals(Vbar, F, FNbar);		
	igl::per_face_normals(V, F, FN);		
	ENbar=Enormal(FNbar, EdgeF);
	EN=Enormal(FN, EdgeF);

/*	for (i=0; i<NTri; i++){
		for(j=0; j<3; j++){
			Epnt.row(3*i+j)=0.5*(V.row(F(i,j))+V.row(F(i,(j+1)%3)));
			Epnt2.row(3*i+j)=Epnt.row(3*i+j)+EN.row(3*i+j);
	//		viewer.data.add_edges(Epnt.row(3*i+j),Epnt2.row(3*i+j),RowVector3d(1,0,0));
		}
	}
*/
//	viewer.launch();
	for(i=0; i<NTri; i++){
	tmpE=Energy(constE, nu, t, V0.row(F(i,0)), V0.row(F(i,1)), V0.row(F(i,2)), Vbar.row(F(i,0)), Vbar.row(F(i,1)), Vbar.row(F(i,2)), V.row(F(i,0)), V.row(F(i,1)), V.row(F(i,2)), ENbar.row(3*i), ENbar.row(3*i+1), ENbar.row(3*i+2), EN.row(3*i), EN.row(3*i+1), EN.row(3*i+2));
        E=E+tmpE;
	}
	
	ofstream myfile;
	myfile.open("output_shell3.txt");
	myfile << "Rectangular plate with slit in the middle." << endl <<"Shift to middle and stretch in Y direction by 1.2*y_position" << endl; 
	myfile << "E1=" << E(0) << " E2=" << E(1) << " Etot=" << E(2) << endl;
	myfile.close();

	return 0;
}


//E1 computes the first energy term, v0 for rest 2D mesh, v for deformed 3D mesh, vbar for rest 3D vertices
Vector3d Energy(double constE, double nu, double t, Vector2d v01, Vector2d v02, Vector2d v03, Vector3d vbar1, Vector3d vbar2, Vector3d vbar3, Vector3d v1, Vector3d v2, Vector3d v3, Vector3d nbar1, Vector3d nbar2, Vector3d nbar3, Vector3d n1, Vector3d n2, Vector3d n3){
	Vector2d e1, e2, e3;
	Vector3d e1_3D, e2_3D, E;
	Matrix2d IA, IAbar, A, IB, IBbar, B;
	double C, E1, E2, Etot, area;
	e1=v02-v01;
	e2=v03-v02;
	e3=v01-v03;
	e1_3D << e1,0;
	e2_3D << e2,0;
	area=0.5*e1_3D.cross(e2_3D).norm();
	IA=I1(v01,v02,v03,v1,v2,v3);
	IAbar=I1(v01,v02,v03,vbar1,vbar2,vbar3);
	IB=I2(e1,e2,e3,v1,v2,v3,n1,n2,n3);
	IBbar=I2(e1,e2,e3,vbar1,vbar2,vbar3,nbar1,nbar2,nbar3);
	A=IAbar.inverse()*(IA-IAbar);
	B=IAbar.inverse()*(IB-IBbar);
	C=constE/(1-nu*nu);
	E1=area*C*t*(nu*pow(A.trace(),2.0)+(1-nu)*((A*A).trace()))/8;	//Add area!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	E2=area*C*pow(t,3)*(nu*pow(B.trace(),2.0)+(1-nu)*((B*B).trace()))/24;	
	Etot=E1+E2;
	E << E1, E2, Etot;
	return E; 
}

Vector2d FoldCircle(double mid, double xpos, double r){
	Vector2d v;
	double theta;
	theta=(xpos-mid)/r;
	v << mid+r*sin(theta), r*cos(theta);
	return v;
}

VectorXi NeighborF(MatrixXi F){
	int NTri=F.rows(),i,j,m,n,flag;
	VectorXi EdgeF(3*NTri);	// EdgeF records the neighboring face of an edge in the tri. mesh in ccw order, #n edge in  face #i is Edge(3*i+n) 
	EdgeF=-1*EdgeF.setOnes(3*NTri);
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
				EdgeF(3*i+j)=m-1;
			}
		}
	}
	return EdgeF;
}
//Enormal calculates the normal of the edges by averaging over neighboring faces and the #n edge in #i face has normal EdgeN(3*i+n)
MatrixXd Enormal(MatrixXd FN, VectorXi EdgeF){	// FN contains the face normals and EdgeF contains the neighboring faces of each edge
	MatrixXd EdgeN(EdgeF.size(),3);
	int i;
	for(i=0; i<EdgeF.size(); i++){
		if (EdgeF(i)==-1){
			EdgeN.row(i) << FN.row(i/3);
		}
		else {
			EdgeN.row(i)=FN.row(i/3)+FN.row(EdgeF(i));
			EdgeN.row(i)=EdgeN.row(i).normalized();
		}	
	}
	return EdgeN;
} 

//I1 computes the matrix of the first fundamental form based on the vertices of the triangle mesh
Matrix2d I1(Vector2d v1, Vector2d v2, Vector2d v3, Vector3d v1d, Vector3d v2d, Vector3d v3d){
	Vector2d e1, e2, e3;
	double q1, q2, q3;
	Matrix2d I;
	e1=v2-v1;
	e2=v3-v2;
	e3=v1-v3;
	q1=(v2d-v1d).squaredNorm();
	q2=(v3d-v2d).squaredNorm();
	q3=(v1d-v3d).squaredNorm();
	I=QMatrix(q1,q2,q3,e1,e2,e3);
	return I;	
}

Matrix2d I2(Vector2d e01, Vector2d e02, Vector2d e03, Vector3d v1, Vector3d v2, Vector3d v3, Vector3d n1, Vector3d n2, Vector3d n3){
	Vector3d e1, e2, e3;
	double q1, q2, q3;
	Matrix2d I;
	e1=v2-v1;
	e2=v3-v2;
	e3=v1-v3;
	q1=2*e1.dot(n3-n2);
	q2=2*e2.dot(n1-n3);
	q3=2*e3.dot(n2-n1);
	I=QMatrix(q1,q2,q3,e01,e02,e03);
	return I;
}

//QMatrix returns the quadratic function Q in discrete triangular mesh as a matrix (operator)
Matrix2d QMatrix(double q1, double q2, double q3, Vector2d e1, Vector2d e2, Vector2d e3){
	Vector2d e1n, e2n, e3n;	
	Vector3d e1_3D,e2_3D;
	Matrix2d Rot, Q;
	double c, area;
	Rot << 0,1,-1,0;
	e1n=Rot*e1;
	e2n=Rot*e2;
	e3n=Rot*e3;
	e1_3D << e1,0;
	e2_3D << e2,0;
	area=0.5*e1_3D.cross(e2_3D).norm();
	c=-1/(8*area*area);
	Q=c*((q1-q2-q3)*e1n*e1n.transpose()+(q2-q3-q1)*e2n*e2n.transpose()+(q3-q1-q2)*e3n*e3n.transpose());
	return Q;
}


//edge computes the edge vector given the vertices of the triangle   
Vector3d edge(Vector3d a, Vector3d b){
        Vector3d c=b-a;
        return c;
}

//rMatrix computes the rotation matrix to generate the dual edges, e1 e2 in ccw ordering
Matrix3d rMatrix(Vector3d e1, Vector3d e2){
	Vector3d axis=e1.cross(e2);
	//axis=axis.normalize();
	Matrix3d rot;	//rot is the planar rotation matrix wrt to the axis
	rot=AngleAxisd(-0.5*M_PI,axis.normalized());
	return rot;
}

