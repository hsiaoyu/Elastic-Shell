#include <iostream>
#include <cmath>
#include <Eigen/Geometry>
#include <igl/per_face_normals.h>
#include <igl/per_edge_normals.h>
using namespace Eigen;

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier);
VectorXd VMass(double, MatrixXi, MatrixXd);
Vector3d edgeV(Vector3d, Vector3d);
Matrix3d rMatrix(Vector3d, Vector3d);
Vector3d TotEnergy(double, double, double, MatrixXi, MatrixXd, MatrixXd, MatrixXd, MatrixXd);
Vector3d Energy(double, double, double, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d);
Matrix2d QMatrix(double, double, double);
MatrixXd Force(double, double, double, MatrixXi, MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXi);
MatrixXd DelI(Vector3d, Vector3d, Vector3d);
VectorXi NeighborF(MatrixXi);
MatrixXi VofEdgeN(MatrixXi, VectorXi); 
Matrix3d CrossM(Vector3d);
Matrix3d dNormVecM(Vector3d);
MatrixXd dEnormM(MatrixXi, MatrixXd);
Vector2d FoldCircle(double,  double, double);
MatrixXd Enormal(MatrixXd, VectorXi);
Matrix2d I2(Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d);
Matrix2d I1(Vector3d, Vector3d, Vector3d);
/*double elength(Eigen::Vector3d a, Eigen::Vector3d b);
int testtest();*/
