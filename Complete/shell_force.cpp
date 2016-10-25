// Last change 10/17/2016 
#include <iostream>
#include <fstream>
#include <cmath>
#include <igl/per_face_normals.h>
#include <igl/readOBJ.h>
#include <igl/viewer/Viewer.h>
#include <igl/per_edge_normals.h>
#include "func_Canon.h"
using namespace std;
using namespace Eigen;
Matrix2d CanonM1, CanonM2, CanonM3;
double area; // area of the canonical triangle

int main() {
	CanonM3 << 1, 0, 0, 0;
	CanonM2 << 1, 1, 1 ,1;
	CanonM1 << 0, 0, 0, 1;
        area = 0.5;
	int i, j, tmp, NTri,NNode,NLen,NWid;
	double L, midL, midW, r, Energy, constE, nu, t, dE1, dE2;
	Vector3d Ei, Ef;
	constE=1.7242e7;
	nu=0.2;
	t=0.8;
	Vector2d newV;
	Vector3d Min, Max;
	MatrixXd Vbar, dI(18,2), V, dV, V2;
	MatrixXi F;
	igl::readOBJ("Vbar.obj", Vbar, F);
	igl::readOBJ("V.obj",V,F); 
	NTri=F.rows();
	NNode=Vbar.rows();	
	MatrixXd FF(2*NNode,3), Epnt(3*NTri,3), Epnt2(3*NTri,3), EN(3*NTri,3), FN(NTri,3), ENbar(3*NTri,3), FNbar(NTri,3), Vnew(NNode,3);
	VectorXi EdgeF(3*NTri);
	EdgeF=NeighborF(F);

//------------------- Finite difference test for the derivative of energy w.r.t the vertices(Force)---------------
        igl::readOBJ("V2.obj",V2,F); //  V is the perturb state and V2 is the initial state
	MatrixXi EdgeNV(3*NTri,4);
	MatrixXd EN2(3*NTri,3), FN2(NTri,3);
	MatrixXd FFtot(2*NNode,3);
	EdgeNV=VofEdgeN(F,EdgeF);
	igl::per_face_normals(Vbar, F, FNbar);
	igl::per_face_normals(V, F, FN);
	igl::per_face_normals(V2, F, FN2);
	ENbar=Enormal(FNbar,EdgeF);		
	EN=Enormal(FN,EdgeF);		
	EN2=Enormal(FN2,EdgeF);		
	Ei=TotEnergy(constE, nu, t, F, Vbar, V2, ENbar, EN2);
	Ef=TotEnergy(constE, nu, t, F, Vbar, V, ENbar, EN);
	cout << "Theoretical dE" << endl << "dE1       dE2" << endl <<  Ef(0)-Ei(0) << "       "<< Ef(1)-Ei(1) << endl;
	FFtot=Force(constE, nu, t, F, Vbar, V2, ENbar, EN2, EdgeNV);   //Change the output for Force to FF << FF1, FF2;
	dV=V-V2;
	dE1=0;
	dE2=0;
	for (i=0; i<NNode; i++){
		dE1+=-FFtot.row(i).dot(dV.row(i));
		dE2+=-FFtot.row(i+NNode).dot(dV.row(i));
	}
	cout << "Finite Difference dE" << endl << "dE1         dE2" << endl << dE1 << "       " << dE2 << endl;
//-------------------------------------------------------------------------------------------------------------

/*//---------------Fintie difference testing for dEN-----------------------------------------------------
	MatrixXd dENM(9*NTri,12), dEN(3*NTri,3);
	dENM=dEnormM(EdgeNV, Vbar);
	dEN.setZero();
	for (i=0;i<3*NTri;i++){
		for (j=0; j<4; j++){
			tmp=EdgeNV(i,j);
			if (tmp==-1){
				tmp=0;
			}
			dEN.row(i)=dEN.row(i)+(dENM.block(3*i,3*j,3,3)*dV.row(tmp).transpose()).transpose();
		}
	}
	
	MatrixXd EN1(3*NTri,3), FN1(NTri,3);
	igl::per_face_normals(Vbar, F, FN);
	igl::per_face_normals(V1, F, FN1);
	EN=Enormal(FN,EdgeF);		
	EN1=Enormal(FN1,EdgeF);		
	cout << "EN1" << endl << EN1 << endl << "EN" << endl << EN << endl << "EN1-EN" << endl << EN1-EN << endl;
	cout << "dEN" << endl << dEN << endl;
*///----------------------------------------------------------------------------------------------------------
/*//----------------------Calculate Total Force------------------------------------------------------------
	Matrix3d FF;
	MatrixXd FFtot1(NNode,3);
	FFtot1= MatrixXd::Zero(NNode,3);
	for(i=0; i<NTri; i++){
	FF=Force(constE, nu, t, V0.row(F(i,0)), V0.row(F(i,1)), V0.row(F(i,2)), Vbar.row(F(i,0)), Vbar.row(F(i,1)), Vbar.row(F(i,2)), V.row(F(i,0)), V.row(F(i,1)), V.row(F(i,2)));
//	cout << FF << endl;
		for (j=0; j<3; j++){
			FFtot1.row(F(i,j))=FFtot1.row(F(i,j))+FF.row(j);
			//cout << FF << endl;
		}
	}
	cout << "Correct F1" << endl << FFtot1 << endl;
*///--------------------------------------------------------------------------------------------------------

/*
//----------------------Output Force-----------------------------------------------------------	
	Vnew= V+FFtot/10;
//	cout << FFtot << endl;
//----------------------------------------------------------------------------------------------	
	
//-----------------Show Force in viewer-----------------------------------------------
	igl::viewer::Viewer viewer;
	viewer.data.set_mesh(V, F);
	viewer.data.add_edges(V,Vnew,RowVector3d(1,0,0));
	viewer.launch();
//--------------------------------------------------------------------------------------
*/
	return 0;
}

Vector3d TotEnergy(double constE, double nu, double t, MatrixXi F, MatrixXd Vbar, MatrixXd V, MatrixXd ENbar, MatrixXd EN){
	int i, j, k, NTri, NNode;
	Vector3d E;
	NTri=F.rows();
	E << 0, 0, 0;
	for (i=0; i<NTri; i++) {
		E+=Energy(constE, nu , t, Vbar.row(F(i,0)), Vbar.row(F(i,1)), Vbar.row(F(i,2)), V.row(F(i,0)), V.row(F(i,1)), V.row(F(i,2)), ENbar.row(3*i), ENbar.row(3*i+1), ENbar.row(3*i+2),  EN.row(3*i), EN.row(3*i+1), EN.row(3*i+2));
	}
	return E;
}


//E1 computes the first energy term, v0 for rest 2D mesh, v for deformed 3D mesh, vbar for rest 3D vertices
Vector3d Energy(double constE, double nu, double t, Vector3d vbar1, Vector3d vbar2, Vector3d vbar3, Vector3d v1, Vector3d v2, Vector3d v3, Vector3d nbar1, Vector3d nbar2, Vector3d nbar3, Vector3d n1, Vector3d n2, Vector3d n3){
	Vector3d E;
	Matrix2d IA, IAbar, A, IB, IBbar, B;
	double C, E1, E2, Etot;
	IA=I1(v1,v2,v3);
	IAbar=I1(vbar1,vbar2,vbar3);
	IB=I2(v1,v2,v3,n1,n2,n3);
	IBbar=I2(vbar1,vbar2,vbar3,nbar1,nbar2,nbar3);
	A=IAbar.inverse()*(IA-IAbar);
	B=IAbar.inverse()*(IB-IBbar);
	C=constE/(1-nu*nu);
	E1=area*C*t*(nu*pow(A.trace(),2.0)+(1-nu)*((A*A).trace()))/8;	//Add area!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	E2=area*C*pow(t,3)*(nu*pow(B.trace(),2.0)+(1-nu)*((B*B).trace()))/24;	
	Etot=E1+E2;
	E << E1, E2, Etot;
	return E; 
}

// Return the force on every vertices
MatrixXd Force(double constE, double nu, double t, MatrixXi F, MatrixXd Vbar, MatrixXd V, MatrixXd ENbar, MatrixXd EN, MatrixXi EdgeNV){
	RowVector3d dval1, dval2, dval3;
	Matrix2d IA, IAbar, A, IB, IBbar, B, Rot, tmp, Inv;
	Matrix3d Ed;
	int i, j, k, NTri, NNode;
	double c1, c2, C, dval;
	NTri=F.rows();
	NNode=V.rows();
	MatrixXd dN, dI(18,2), Tr(2,3), FF(2*NNode,3), FF2(NNode,3), FF1(NNode,3);
	dN=dEnormM(EdgeNV, V); //dN is the derivative of egde normal w.r.t. the four vectices cooresponding to the edge
	FF2.setZero();
	FF1.setZero();
	for (i=0; i<NTri; i++) {
		Ed << V.row(F(i,1))-V.row(F(i,0)), V.row(F(i,2))-V.row(F(i,1)), V.row(F(i,0))-V.row(F(i,2)); //Ed.row(i) is the ith  edge vector for deformed triangle
		C=constE/(1-nu*nu);
		c1=C*nu*area/24;
		c2=C*(1-nu)*area/24;
		IB=I2(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)),EN.row(3*i),EN.row(3*i+1),EN.row(3*i+2));
		IBbar=I2(Vbar.row(F(i,0)),Vbar.row(F(i,1)),Vbar.row(F(i,2)),ENbar.row(3*i),ENbar.row(3*i+1),ENbar.row(3*i+2));
		IA=I1(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)));
		IAbar=I1(Vbar.row(F(i,0)),Vbar.row(F(i,1)),Vbar.row(F(i,2)));
		Inv=IAbar.inverse();
                B=Inv*(IB-IBbar);
		A=Inv*(IA-IAbar);	
		dI=DelI(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)));
		Tr << B.trace()*(Inv*CanonM1).trace(), B.trace()*(Inv*CanonM2).trace(), B.trace()*(Inv*CanonM3).trace(), (B*Inv*CanonM1).trace(), (B*Inv*CanonM2).trace(), (B*Inv*CanonM3).trace();
		for (j=0; j<3; j++){
			for (k=0; k<4; k++){
				if (EdgeNV(3*i+j,k)!=-1){
					dval1=2*Ed.row((j+1)%3)*dN.block(9*i+3*j,3*k,3,3);	// orginate from the term 2<dn(j),e(j+1)>, dn(j)(dV(EdgeNV(j,k)))
					dval2=-2*Ed.row((j+2)%3)*dN.block(9*i+3*j,3*k,3,3); // originate form the term -2<dn(j),e(j-1)>
					FF2.row(EdgeNV(3*i+j,k))+= -2*pow(t,3)*(c1*(Tr(0,(j+1)%3)-Tr(0,j)-Tr(0,(j+2)%3))+c2*(Tr(1,(j+1)%3)-Tr(1,j)-Tr(1,(j+2)%3)))*dval1;
					FF2.row(EdgeNV(3*i+j,k))+= -2*pow(t,3)*(c1*(Tr(0,(j+2)%3)-Tr(0,j)-Tr(0,(j+1)%3))+c2*(Tr(1,(j+2)%3)-Tr(1,j)-Tr(1,(j+1)%3)))*dval2;
				}
			}
			dval3=2*(EN.row(3*i+((j+2)%3))-EN.row(3*i+((j+1)%3))); // <n3-n2,dv1>, <n3-n3,dv2>
			FF2.row(F(i,j))+= 2*pow(t,3)*(c1*(Tr(0,j)-Tr(0,(j+1)%3)-Tr(0,(j+2)%3))+c2*(Tr(1,j)-Tr(1,(j+1)%3)-Tr(1,(j+2)%3)))*dval3;// dQ(j)(dV(j))
			FF2.row(F(i,(j+1)%3))+= -2*pow(t,3)*(c1*(Tr(0,j)-Tr(0,(j+1)%3)-Tr(0,(j+2)%3))+c2*(Tr(1,j)-Tr(1,(j+1)%3)-Tr(1,(j+2)%3)))*dval3;//dQ(j)(dV(j+1))
		}
		for (j=0; j<3; j++){
			for(k=0; k<3; k++){
				tmp=Inv*dI.block(6*j+2*k,0,2,2);
				FF1(F(i,j),k)+=-6*t*(c1*A.trace()*tmp.trace()+c2*(A*tmp).trace()); // Force due to the fisrt fundamental term
			}	
		}
	}
	C=-1/(8*area*area);
	FF2=C*FF2;
	FF<<FF1,FF2;
	//FF = FF1+FF2;
	return FF;
}

MatrixXd DelI(Vector3d v1, Vector3d v2, Vector3d v3){ // return delta I w.r.t to the change in Vi, V2, V3
	MatrixXd T1(6,2), T2(6,2), T3(6,2), dI(18,2);
	double c;
	int i;
	for (i=0; i<3; i++){
		T1.block(2*i,0,2,2)=2*(v3(i)-v2(i))*CanonM1+2*(v3(i)+v2(i)-2*v1(i))*CanonM2+2*(v2(i)-v3(i))*CanonM3;
		T2.block(2*i,0,2,2)=2*(v3(i)-v1(i))*CanonM1+2*(v1(i)-v3(i))*CanonM2+2*(v1(i)+v3(i)-2*v2(i))*CanonM3;
		T3.block(2*i,0,2,2)=2*(v1(i)+v2(i)-2*v3(i))*CanonM1+2*(v1(i)-v2(i))*CanonM2+2*(v2(i)-v1(i))*CanonM3;
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

MatrixXi VofEdgeN(MatrixXi F, VectorXi EdgeF){  // the related vertices to an edge normal
	int NEdge=EdgeF.size(), i ,j, NTri=F.rows(), m;
	MatrixXi EdgeNV(NEdge,4);
	EdgeNV.col(3).setConstant(-1);
	for (i=0; i<NTri; i++){
		for (j=0; j<3; j++){
			EdgeNV(3*i+j,0)= F(i,j);
			EdgeNV(3*i+j,1)= F(i,(j+1)%3);
			EdgeNV(3*i+j,2)= F(i,(j+2)%3);
			if (EdgeF(3*i+j)!= -1){
				for (m=0; m<3; m++){
					if(F(i,j)!=F(EdgeF(3*i+j),m)&&F(i,(j+1)%3)!=F(EdgeF(3*i+j),m)){
						EdgeNV(3*i+j,3)=F(EdgeF(3*i+j),m);
					}
				}
			}	
		}
	} 
	return EdgeNV;
}

Matrix3d CrossM(Vector3d v){ // return the matrix A such that axb=Ab, A is the cross product matrix
	Matrix3d M;
	M << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
	return M;
}

Matrix3d dNormVecM(Vector3d v){ // return the matrix M such that for F(x)=x/|x|, dF=M*dx
	Matrix3d M, I;
	double a=1/v.norm();
	I.setIdentity();
	M=a*I-pow(a,3)*v*v.transpose();
	return M;
}

MatrixXd dEnormM(MatrixXi EdgeNV, MatrixXd V){ // return the derivative of each edge normal w. the vertices related to the normal of edge
	// EdgeNV is the matrix that contains the related vertices to an edge, matrix V is just the vertice coordinate matrix
	int i, NEdge=EdgeNV.rows();
	Matrix3d Mn, M1, M2, c1_1, c2_1, c1_2, c2_2;
	MatrixXd M(3*NEdge,12);   // M.block(3*i,3*j,3,3) corresponse to the derivative of the i-th edge with EdgeV(i,j) component
	Vector3d e1_1, e2_1, e1_2, e2_2;
	for (i=0; i<NEdge; i++){
		e1_1=V.row(EdgeNV(i,1))-V.row(EdgeNV(i,0));   // edg1 1 in face 1 (current triangle)
		e2_1=V.row(EdgeNV(i,2))-V.row(EdgeNV(i,1));
		e1_2=V.row(EdgeNV(i,0))-V.row(EdgeNV(i,1));   // edge 1 in face 2 (neighboring triangle)
		if (EdgeNV(i,3)==-1){	// no neighboring face for Edge(i)
			Mn.setIdentity();
			M2.setZero();
			e2_2.setZero();
		}
		else {
			e2_2=V.row(EdgeNV(i,3))-V.row(EdgeNV(i,0));
			Mn=dNormVecM(e1_1.cross(e2_1).normalized()+e1_2.cross(e2_2).normalized());
			M2=dNormVecM(e1_2.cross(e2_2));
		}
		c1_1=CrossM(e1_1);
		c2_1=CrossM(e2_1);
		c1_2=CrossM(e1_2);
		c2_2=CrossM(e2_2);
		M1=dNormVecM(e1_1.cross(e2_1));
		M.block(3*i,0,3,3)=Mn*(M1*c2_1-M2*(c1_2+c2_2));
		M.block(3*i,3,3,3)=Mn*(-M1*(c1_1+c2_1)+M2*c2_2);
		M.block(3*i,6,3,3)=Mn*M1*c1_1;
		M.block(3*i,9,3,3)=Mn*M2*c1_2;
	}
	return M;
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

//QMatrix returns the quadratic function Q in discrete triangular mesh as a matrix (operator)
Matrix2d QMatrix(double q1, double q2, double q3){
	Matrix2d Q;
	double c;
	c=-1/(8*area*area);
	Q=c*((q1-q2-q3)*CanonM1+(q2-q3-q1)*CanonM2+(q3-q1-q2)*CanonM3);
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

