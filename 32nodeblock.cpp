#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

int main(){
    const double T_west = 100.0; // west wall temperature
    const double L = 1.0; //length of rod, meter
    const double k = 1.0; //watts per meter kelvin
    const int N = 32; //discretization qty - dimension of array
    const double h = 100.0, Tinf = 25.0;
    const double dx = L/N;       
    const double C = h*Tinf / (k/dx + h/2.0);
    double b = (k/dx - (h/2))/(k/dx + (h/2));
    
    double x[N];
    for (int i = 0; i < N; ++i){
        x[i] = (i+0.5)*dx;
    }

    double q[N]; // watts per cubic meter, volumetric heat gen, 10x+5
    for (int i = 0; i < N; ++i){
        q[i] = (10.0*x[i])+5;
    }

    double rhs[N];
    double delx2k = dx*dx / k;
    double A[N][N]={}; //initialize Matrix A

    //interior cells
    for (int i = 0; i <= N-3; ++i){
        //left hand side
        A[i+1][i] = -1;
        A[i+1][i+1] = 2;
        A[i+1][i+2] = -1;
    }

    //right hand side
    for (int i = 0; i < N; ++i){
        rhs[i] = q[i]*delx2k;
        if (i == 0) rhs[i] += 2.0*T_west;
        else if (i == N-1) rhs[i] += C;
    }

    //first row - left wall dirichlet method
    A[0][0] = 3.0;
    A[0][1] = -1.0;
    A[N-1][N-1] = 2-b; //T32
    A[N-1][N-2] = -1.0; //T31
    
    std::vector<double> T(N, Tinf);     // initial guess
    int iter = 0;
    double change;
    do {
        change = 0.0;
        for (int i = 0; i < N; ++i) {
            double Tnew;
            if (i == 0)                 
                Tnew = (2*T_west + T[1] + q[0]*dx*dx/k) / 3.0;
            else if (i == N-1)          
                Tnew = (T[N-2] + C + q[N-1]*dx*dx/k) / (2.0 - b);
            else        
                Tnew = (T[i-1] + T[i+1] + q[i]*dx*dx/k) / 2.0;
            change = std::fmax(change, std::fabs(Tnew - T[i]));
            T[i] = Tnew;
        }
        ++iter;
    } while (change > 1e-10 && iter < 200000);

    std::cout << "GS sweeps: " << iter << std::endl;
       std::cout << "node      x         T\n";
       for (int i = 0; i < N; ++i)
           std::cout << std::setw(4) << i+1
                     << "  " << std::setw(7) << x[i]
                     << "  " << std::setw(10) << T[i] << std::endl;

       
       double maxres = 0.0;
       for (int i = 0; i < N; ++i) {
           double s = rhs[i];
           for (int j = 0; j < N; ++j) s -= A[i][j] * T[j];
           if (std::fabs(s) > maxres) maxres = std::fabs(s);
       }
       std::cout << "max |A*T - rhs| = " << maxres << std::endl;
       
    return 0;
}
