// This file simulate shinkrage or inflation of paper based on the defined IAbar and IBbar
// The final result of V is the stable shape of the thin shell based on the defined shinkrage and inflation rule
#include <iostream>
#include <fstream>
#include <cmath>
#include <igl/per_face_normals.h>
#include <igl/readOBJ.h>
#include <igl/readDMAT.h>
#include <igl/per_edge_normals.h>
#include <igl/viewer/Viewer.h>
#include "burning_static.h"
#include <lbfgs.h>
using namespace std;
using namespace Eigen;
Matrix2d CanonM1, CanonM2, CanonM3;
double area=0.5, constE, nu, t, delt; // area of the canonical triangle, constS is the shrinkage constant
MatrixXd V, Vbar, ENbar, Vfix, Ibartot;// ENbar is the edge normal of Vbar
MatrixXi F, EdgeNV; // Eash row of EdgeNV(3*NTri,4) stores the index off the four vertices that are related to one edge
VectorXi EdgeF;// EdgeF(3*NTri) stores the index of the  adjecent face of the edge other than the face that is indicated by the vector index 
int NNode; // NNode is the total # of nodes
int* FixIndex;
Vector3d V_ini1, V_ini2; // Vini1, Vini2 are the fixed points of V when enabled gravity

class objective_function
{
protected:
    lbfgsfloatval_t *m_x;

public:
    objective_function() : m_x(NULL)
    {
    }

    virtual ~objective_function()
    {
        if (m_x != NULL) {
            lbfgs_free(m_x);
            m_x = NULL;
        }
    }

    int run(int Nvar)
    {
        lbfgsfloatval_t fx;
        lbfgsfloatval_t *m_x = lbfgs_malloc(Nvar);

        if (m_x == NULL) {
            printf("ERROR: Failed to allocate a memory block for variables.\n");
            return 1;
        }
	
	int nskip=0;
        for (int i = 0; i < NNode; i++) {
		if(i==FixIndex[nskip]){
			nskip++;
		}	
		else {
			m_x[3*(i-nskip)] = V(i,0);
			m_x[3*(i-nskip)+1] = V(i,1);
			m_x[3*(i-nskip)+2] = V(i,2);
		}
        }


	// The fifth component is the maximum iternation allowed
        lbfgs_parameter_t _defparam = {
        6, 1e-5, 0, 1e-5,
        50000, LBFGS_LINESEARCH_DEFAULT, 40,
        1e-20, 1e20, 1e-4, 0.9, 0.9, 1.0e-16,
        0.0, 0, -1,
        };


	// last argument=NULL uses defualt value
        int ret = lbfgs(Nvar, m_x, &fx, _evaluate, _progress, this, &_defparam);

        printf("L-BFGS optimization terminated with status code = %d\n", ret);
        printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, m_x[0], m_x[1]);
        
        return ret;
    }

protected:
    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        return reinterpret_cast<objective_function*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        lbfgsfloatval_t fx = 0.0;

	int nskip=0;
        for (int i = 0; i < NNode; i++) {
		if(i==FixIndex[nskip]){
			V.row(i)=Vfix.row(nskip);
			nskip++;
		}	
		else {
			V(i,0)=x[3*(i-nskip)];
			V(i,1)=x[3*(i-nskip)+1];
			V(i,2)=x[3*(i-nskip)+2];
		}
        }
        MatrixXd FF;
	FF = -1*Force(); 
	nskip = 0;
        for (int i = 0; i < NNode; i++) {
		if(i==FixIndex[nskip]){
			nskip++;
		}	
		else {
			g[3*(i-nskip)] = FF(i,0);
			g[3*(i-nskip)+1] = FF(i,1);
			g[3*(i-nskip)+2] = FF(i,2);
		}
        }
        fx=TotEnergy();
        return fx;
    }

    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        return reinterpret_cast<objective_function*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        printf("Iteration %d:\n", k);
        cout << fx << endl;
	return 0;
    }
};

int main (){
//int main(int argc, char *argv){
    CanonM3 << 1, 0, 0, 0;
    CanonM2 << 1, 1, 1 ,1;
    CanonM1 << 0, 0, 0, 1;
    int i, j, tmp, NTri, Nfix;
    Vector3d Ei, Ef;
    ifstream Infile("Input_Filename.txt");
    string str;
    getline(Infile,str); 
    getline(Infile,str); 
    getline(Infile,str); 
    getline(Infile,str); 
    const char *cstr = str.c_str();
    igl::readDMAT(cstr, Ibartot);
    ifstream InFile;
    getline(Infile,str); 
    const char *cstr2 = str.c_str();
    InFile.open(cstr2);
    InFile >> constE;
    InFile >> nu;
    InFile >> t;
 //   igl::readOBJ("Vbar.obj", Vbar, F);
    igl::readOBJ("V.obj",V,F);
    InFile >> Nfix;
    FixIndex = new int [Nfix];  //Create an array of size Nfix
    NTri=F.rows();
    NNode=V.rows();	
    MatrixXd Vtemp(Nfix,3);
    for (i=0; i<Nfix; i++){
   	 InFile >> FixIndex[i];
	 Vtemp.row(i)=V.row(FixIndex[i]);
    }
    if (i==0){
	FixIndex[0]=-1;
    }
    cout << FixIndex[0] << endl;
    Vfix=Vtemp;
    MatrixXd FFtot(NNode,3), EN(3*NTri,3), FN(NTri,3), FNbar(NTri,3), Vnew(NNode,3);
    EdgeF=NeighborF();
    EdgeNV=VofEdgeN();
//   igl::per_face_normals(Vbar, F, FNbar);
//    ENbar=Enormal(FNbar);
    
    objective_function obj;
    obj.run(3*(NNode-Nfix));
  
    getline(Infile,str); 
    const char *cstr1 = str.c_str();
    cout << cstr1 << endl;
    igl::writeOBJ(cstr1,V,F);
   
  // igl::viewer::Viewer viewer;
   //viewer.data.set_mesh(V, F);
   //viewer.launch();
return 0;
}


VectorXd VMass(double rho){
	int i;
        double areaT;
        Vector3d e1 , e2;
        VectorXd M(V.rows());
        M.setZero();
	for (i=0; i<F.rows(); i++){
        	e1 = V.row(F(i,1))-V.row(F(i,0));
        	e2 = V.row(F(i,2))-V.row(F(i,1));
		areaT=0.5*e1.cross(e2).norm();
                M(F(i,0))+=rho*areaT/3;
                M(F(i,1))+=rho*areaT/3;
                M(F(i,2))+=rho*areaT/3;
	}
	return M;
}

double TotEnergy( ){
	int i, j, k, NTri;
        double E;
	NTri=F.rows();
	MatrixXd FN(NTri,3), EN(3*NTri,3);
	E = 0;
	igl::per_face_normals(V, F, FN);
	EN=Enormal(FN);
	for (i=0; i<NTri; i++) {
		E+=Energy(constE, nu , t, i, V.row(F(i,0)), V.row(F(i,1)), V.row(F(i,2)), EN.row(3*i), EN.row(3*i+1), EN.row(3*i+2));
	}
	return E;
}


//E1 computes the first energy term, v0 for rest 2D mesh, v for deformed 3D mesh, vbar for rest 3D vertices
double Energy(double constE, double nu, double t, int FaceIndex, Vector3d v1, Vector3d v2, Vector3d v3,Vector3d n1, Vector3d n2, Vector3d n3){
	Vector3d E;
	Matrix2d IA, IAbar, A, IB, IBbar, B;
	double C, E1, E2, Etot;
	IA=I1(v1,v2,v3);
	IB=I2(v1,v2,v3,n1,n2,n3);
	IAbar=Ibartot.block(2*FaceIndex,0,2,2);
	IBbar=Ibartot.block(2*(FaceIndex+F.rows()),0,2,2);
	A=IAbar.inverse()*(IA-IAbar);
	B=IAbar.inverse()*(IB-IBbar);
	C=constE/(1-nu*nu);
	E1=area*C*t*(nu*pow(A.trace(),2.0)+(1-nu)*((A*A).trace()))/8;	//Add area!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	E2=area*C*pow(t,3)*(nu*pow(B.trace(),2.0)+(1-nu)*((B*B).trace()))/24;	
	Etot=E1+E2;
        E << E1, E2, Etot;
	return Etot; 
}

// Return the force on every vertices
MatrixXd Force(){
	RowVector3d dval1, dval2, dval3;
	Matrix2d IA, IAbar, A, IB, IBbar, B, Rot, tmp, Inv;
	Matrix3d Ed;
	int i, j, k, NTri;
	double c1, c2, C, dval;
	NTri=F.rows();
	MatrixXd dN, dI(18,2), Tr(2,3), FF(NNode,3), FF2(NNode,3), FF1(NNode,3), FN(NTri,3), EN(3*NTri,3);
	igl::per_face_normals(V, F, FN);
	EN=Enormal(FN);	
	dN=dEnormM(); //dN is the derivative of egde normal w.r.t. the four vectices cooresponding to the edge
        FF2.setZero();
	FF1.setZero();
	/*IAbar.setIdentity();
	IAbar=(1.+constS)*IAbar/2;
	IBbar.setIdentity();
	IBbar=(constS-1.)*IBbar/(2*t);*/
	for (i=0; i<NTri; i++) {
		Ed << V.row(F(i,1))-V.row(F(i,0)), V.row(F(i,2))-V.row(F(i,1)), V.row(F(i,0))-V.row(F(i,2)); //Ed.row(i) is the ith  edge vector for deformed triangle
		C=constE/(1-nu*nu);
		c1=C*nu*area/24;
		c2=C*(1-nu)*area/24;
		IB=I2(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)),EN.row(3*i),EN.row(3*i+1),EN.row(3*i+2));
		IA=I1(V.row(F(i,0)),V.row(F(i,1)),V.row(F(i,2)));
		IAbar=Ibartot.block(2*i,0,2,2);
		IBbar=Ibartot.block(2*(i+F.rows()),0,2,2);
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
               //cout << FF1 << endl << FF2 << endl;
	}
	C=-1/(8*area*area);
	FF2=C*FF2;
	//FF<<FF1,FF2;
	FF = FF1+FF2;
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
	Vector2d vec;
	double theta;
	theta=(xpos-mid)/r;
	vec << mid+r*sin(theta), r*cos(theta);
	return vec;
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

MatrixXi VofEdgeN( ){  // the related vertices to an edge normal
	int NEdge=EdgeF.size(), i ,j, NTri=F.rows(), m;
	MatrixXi EdgeNVf(NEdge,4);
	EdgeNVf.col(3).setConstant(-1);
	for (i=0; i<NTri; i++){
		for (j=0; j<3; j++){
			EdgeNVf(3*i+j,0)= F(i,j);
			EdgeNVf(3*i+j,1)= F(i,(j+1)%3);
			EdgeNVf(3*i+j,2)= F(i,(j+2)%3);
			if (EdgeF(3*i+j)!= -1){
				for (m=0; m<3; m++){
					if(F(i,j)!=F(EdgeF(3*i+j),m)&&F(i,(j+1)%3)!=F(EdgeF(3*i+j),m)){
						EdgeNVf(3*i+j,3)=F(EdgeF(3*i+j),m);
					}
				}
			}	
		}
	} 
	return EdgeNVf;
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

MatrixXd dEnormM( ){ // return the derivative of each edge normal w. the vertices related to the normal of edge
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

