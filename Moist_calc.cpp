// This file generate the machine direction in barycentric coordinate of each face and store in MDinFace.dmat
// Read out MDinFace.dmat using readDMAT and each row of the matrix is the barycentric coordinate representation of MD 
#include <iostream>
#include <fstream>
#include <cmath>
#include <igl/readOBJ.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
using namespace std;
using namespace Eigen;
MatrixXd V;// ENbar is the edge normal of Vbar
MatrixXi F; // Eash row of EdgeNV(3*NTri,4) stores the index off the four vertices that are related to one edge
int main (){
    igl::readOBJ("Vbar.obj", V, F);
    VectorXd Moisture(F.rows());
    double moist;
    cout << "What is the uniform moisture difference(%):" << endl;
    cin >> moist;
    for(int i=0; i< F.rows(); i++){
	Moisture(i)=moist;
    }    
    igl::writeDMAT("Input_MoistureLevel.dmat",Moisture,1);
    return 0;
}
