#include <iostream>
#include <cmath>
#include <Eigen/Geometry>
using namespace Eigen;

MatrixXd CalcIbar(MatrixXd, double, double, double, VectorXd);
Matrix2d I1(Vector3d, Vector3d, Vector3d);
Matrix2d I2(Vector3d, Vector3d, Vector3d, Vector3d, Vector3d, Vector3d);
MatrixXd Enormal(MatrixXd);
VectorXi NeighborF( );
Matrix2d QMatrix(double, double, double);
