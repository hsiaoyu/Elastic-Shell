#include <iostream>
#include <cmath>
#include <Eigen/Geometry>
#include <igl/per_face_normals.h>
#include <igl/per_edge_normals.h>
using namespace Eigen;

Matrix2d test(int);
Vector3d edgeV(Vector3d, Vector3d);
Matrix3d rMatrix(Vector3d, Vector3d);
Vector3d TotEnergy(double, double, double, MatrixXi, MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd);
Vector3d Energy(double, double, double, Vector2d, Vector2d, Vector2d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d);
Matrix2d QMatrix(double, double, double, Vector2d, Vector2d, Vector2d);
MatrixXd Force(double, double, double, MatrixXi, MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXd, MatrixXi);
MatrixXd DelI(Vector2d, Vector2d, Vector2d, Vector3d, Vector3d, Vector3d);
VectorXi NeighborF(MatrixXi);
MatrixXi VofEdgeN(MatrixXi, VectorXi); 
Matrix3d CrossM(Vector3d);
Matrix3d dNormVecM(Vector3d);
MatrixXd dEnormM(MatrixXi, MatrixXd);
Vector2d FoldCircle(double,  double, double);
MatrixXd Enormal(MatrixXd, VectorXi);
Matrix2d I2(Vector2d, Vector2d, Vector2d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d);
Matrix2d I1(Vector2d, Vector2d, Vector2d, Vector3d, Vector3d, Vector3d );
/*double elength(Eigen::Vector3d a, Eigen::Vector3d b);
int testtest();*/
