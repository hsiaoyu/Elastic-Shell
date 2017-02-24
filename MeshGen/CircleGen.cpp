// This file generate the points on the circle
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;
#define PI 3.14159265
int main() {
    double rad, theta;
    int N,i;
    cout << "Enter radius" << endl;
    cin >> rad;
    cout << "Enter number of points on the circle" << endl;
    cin >> N;
    ofstream Outfile("Circle.poly");
    Outfile << N << " 2 0 0" << endl;
    for (i=0; i<N; i++){
	theta=2*PI*i/N;
        Outfile << i+1 << " " << rad*cos(theta) << " " << rad*sin(theta) << endl;
    }
    Outfile << N << " 0" << endl;
    for (i=1; i<N; i++){
         Outfile << i << " " << i << " " << i+1 << endl;
    }
    Outfile << N << " " << N << " " << 1 << endl;  
    Outfile << "0" << endl;
    return 0;
}
