// Last change 9/12/2016 
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
	MatrixXd Vbar, dI(18,2);
	MatrixXi F;
	igl::readOBJ("rectangle_slit.obj", Vbar, F);
	
	
	NTri=F.rows();
	NNode=Vbar.rows();	
	MatrixXd FFtot(NNode,3),V0(NNode,2), V(NNode,3),Epnt(3*NTri,3), Epnt2(3*NTri,3), EN(3*NTri,3), FN(NTri,3), ENbar(3*NTri,3), FNbar(NTri,3), Vnew(NNode,3);
	VectorXi EdgeF(3*NTri);
	EdgeF=NeighborF(F);
	cout << "test" << endl;
	Min=Vbar.colwise().minCoeff();
	Max=Vbar.colwise().maxCoeff();
	L=(Max(0)-Min(0));
	midL=(Max(0)+Min(0))/2;
	midW=(Max(1)+Min(1))/2;
	r=L/M_PI;
	cout << "L is " << L << " r is " << r << " mid is " << midL << endl;
/*	newV=FoldCircle(midL, midL, r);
	cout << newV << endl;
	newV=FoldCircle(midL, midL-(L/2), r);
	cout << newV << endl;
*/
	//--------------Streching the rectangle--------------------------
	for (i=0; i<NNode; i++){
		V0.row(i) << Vbar(i,0), Vbar(i,1);
		V.row(i) << Vbar(i,0), (Vbar(i,1)-midW)*1.2 , Vbar(i,2);
	}
	//----------------------------------------------------------------
	igl::viewer::Viewer viewer;
	viewer.data.set_mesh(V, F);
	Matrix3d FF;
	FFtot= MatrixXd::Zero(NNode,3);
	for(i=0; i<NTri; i++){
	FF=Force(constE, nu, t, V0.row(F(i,0)), V0.row(F(i,1)), V0.row(F(i,2)), Vbar.row(F(i,0)), Vbar.row(F(i,1)), Vbar.row(F(i,2)), V.row(F(i,0)), V.row(F(i,1)), V.row(F(i,2)));
//	cout << FF << endl;
		for (j=0; j<3; j++){
			FFtot.row(F(i,j))=FFtot.row(F(i,j))+FF.row(j);
			//cout << FF << endl;
		}
	}
	
	Vnew= V+FFtot/50000;
	cout << FFtot << endl;
	
//	cout << Vnew << endl;
	viewer.data.add_edges(V,Vnew,RowVector3d(1,0,0));
	viewer.launch();
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

Matrix3d Force(double constE, double nu, double t, Vector2d v01, Vector2d v02, Vector2d v03, Vector3d vbar1, Vector3d vbar2, Vector3d vbar3, Vector3d v1, Vector3d v2, Vector3d v3){
	// F(i,j) returns the force of the i-th vertex and the j-th coordinate
	MatrixXd dI(18,2);
	Vector2d e1, e2, e3;
	Vector3d e1_3D, e2_3D;
	Matrix2d IA, IAbar, A, tmp;
	Matrix3d F1; 
	int i, j;
	double c1, c2, area, C;
	e1=v02-v01;
	e2=v03-v02;
	e3=v01-v03;
	e1_3D << e1,0;
	e2_3D << e2,0;
	area=0.5*e1_3D.cross(e2_3D).norm();
	dI=DelI(v01,v02,v03,v1,v2,v3);
	IA=I1(v01,v02,v03,v1,v2,v3);
	IAbar=I1(v01,v02,v03,vbar1,vbar2,vbar3);
	A=IAbar.inverse()*(IA-IAbar);	
	C=constE/(1-nu*nu);
	c1=C*nu*area/8; //////////////////////////////////////Check~~~~~~~~~~~
	c2=C*(1-nu)*area/8;
	for (i=0; i<3; i++){
		for(j=0; j<3; j++){
			tmp=dI.block(6*i+2*j,0,2,2);
			F1(i,j)=-2*t*(c1*A.trace()*tmp.trace()+c2*(A*tmp).trace()); // Force due to the fisrt fundamental term
		}
	}
	return F1;
}

MatrixXd DelI(Vector2d v01, Vector2d v02, Vector2d v03, Vector3d v1, Vector3d v2, Vector3d v3){
	Matrix2d E1, E2, E3, Rot, tmp;
	MatrixXd T1(6,2), T2(6,2), T3(6,2), dI(18,2);
	Vector2d e1,e2,e3,e1n, e2n, e3n;
	Vector3d e1_3D, e2_3D;
	double area, c;
	int i;
	e1=v02-v01;
	e2=v03-v02;
	e3=v01-v03;
	Rot << 0,1,-1,0;
	e1_3D << e1,0;
	e2_3D << e2,0;
	area=0.5*e1_3D.cross(e2_3D).norm();
	e1n=Rot*e1;
	e2n=Rot*e2;
	e3n=Rot*e3;
	E1=e1n*e1n.transpose();
	E2=e2n*e2n.transpose();
	E3=e3n*e3n.transpose();
	for (i=0; i<3; i++){
		T1.block(2*i,0,2,2)=2*(v3(i)-v2(i))*E1+2*(v3(i)+v2(i)-2*v1(i))*E2+2*(v2(i)-v3(i))*E3;
		T2.block(2*i,0,2,2)=2*(v3(i)-v1(i))*E1+2*(v1(i)-v3(i))*E2+2*(v1(i)+v3(i)-2*v2(i))*E3;
		T3.block(2*i,0,2,2)=2*(v1(i)+v2(i)-2*v3(i))*E1+2*(v1(i)-v2(i))*E2+2*(v2(i)-v1(i))*E3;
	}
	dI << T1, T2, T3;
	c=-1/(8*area*area);
	dI=c*dI;
	return dI;
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

