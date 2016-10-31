#include <iostream>
#include <cmath>
#include <Eigen/Geometry>
#include <igl/per_face_normals.h>
#include <igl/per_edge_normals.h>
using namespace Eigen;

bool pre_draw(igl::viewer::Viewer& viewer);
void Sim();
VectorXd VMass(double);
Vector3d edgeV(Vector3d, Vector3d);
Matrix3d rMatrix(Vector3d, Vector3d);
Vector3d TotEnergy();
Vector3d Energy(double, double, double, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d);
Matrix2d QMatrix(double, double, double);
MatrixXd Force();
MatrixXd DelI(Vector3d, Vector3d, Vector3d);
VectorXi NeighborF();
MatrixXi VofEdgeN(); 
Matrix3d CrossM(Vector3d);
Matrix3d dNormVecM(Vector3d);
MatrixXd dEnormM();
Vector2d FoldCircle(double,  double, double);
MatrixXd Enormal(MatrixXd);
Matrix2d I2(Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d);
Matrix2d I1(Vector3d, Vector3d, Vector3d);
/*double elength(Eigen::Vector3d a, Eigen::Vector3d b);
int testtest();*/
