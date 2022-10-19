#include <cstdlib>
#include<iostream>
#include<math.h>

using namespace std;

const int n = 50;
const int lymda = 2;
const double a = 0;
const double b = 1;

double K(double x, double y) {
    return(x - y);
}

//
//        double middlepryam(double a, double b, double xi) {
//      double nn=1000,h,x,in=0,i;
//     h=(b-a)/nn;
//      x=a+(h/2);
//    // cout<<" x= "<<x<<endl;
//      while (x<=b-(h/2)) {
//                  in=in+(K(xi,x)*h);
//                  x=x+h;
//       
//      
//      }
//     
//    return in;   
//   }
//    

double middlepryam2(double a, double b, double a1, double b1) {
    double nn = 1000, h, h1, x, x1, in = 0, i;
    h = (b - a) / nn;
    h1 = (b1 - a1) / nn;
    x = a + (h / 2);
    x1 = a1 + (h1 / 2);
    //cout<<" x= "<<x<<endl;
    while (x1 <= b1 - (h1 / 2)) {
        while (x <= b - (h / 2)) {
            in = in + (K(x1, x) * h);
            x = x + h;


        }
        x1 = x1 + h1;
    }
    return in;
}

double del(int i, int j) {
        return(i == j);
}

double phi(double xi, int i) {
    double x[n + 1], h, s;
    int j;

    h = (b - a) / n;

    for (j = 0; j < n + 1; j++) {

        x[j] = a + j * h;

    }
    if ((xi >= x[i]) && (xi <= x[i + 1])) {
        s = 1;
    }
    else {
        s = 0;
    }
    return(s);
}

double u(double xi, double c[n]) {
    int i;
    double s = 0;
    for (i = 0; i < n; i++) {
        s = s + c[i] * phi(xi, i);

    }
    return(s);
}

void Gauss(int k, double Matrix[n][n + 1]) {
    if (Matrix[k][k] != 1) {
        double T = Matrix[k][k];
        for (int j = k; j < n + 1; j++) {
            Matrix[k][j] = Matrix[k][j] / T;
        }
    }
    for (int i = 0; i < n; i++) {
        if ((Matrix[i][k] != 0) && (i != k)) {
            double T = Matrix[i][k];
            Matrix[i][k] = 0;
            for (int j = k + 1; j < n + 1; j++) {
                Matrix[i][j] -= Matrix[k][j] * T;
            }
        }
    }
    if (k < n - 1) {
        Gauss(k + 1, Matrix);
    }
}

int main(int argc, char** argv) {

    double h, x[n + 1], xi[n], A[n][n + 1], c[n];
    int i, j, k;

    h = (b - a) / n;

    for (i = 0; i < n + 1; i++) {

        x[i] = a + i * h;
        // cout<< x[i]<<endl; 
    }
    for (i = 0; i < n; i++) {

        xi[i] = x[i] + (h / 2);
        //  cout<< xi[i]<<endl; 
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            //cout<<" j= "<<j<<" x (j)= "<<x[j]<<" x (j+1)= "<<x[j+1]<<endl;
            A[i][j] = del(i, j) - lymda * middlepryam2(x[j], x[j + 1], x[i], x[i + 1]);
        }
        A[i][n] = (xi[i] * xi[i]) - lymda * ((xi[i] / 3) - 0.25);
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n + 1; j++) {

            cout << A[i][j] << " ";
        }
        cout << endl;
    }

    Gauss(0, A);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n + 1; j++) {

            cout << A[i][j] << " ";
        }
        c[i] = A[i][n];
        cout << endl;
    }

    for (i = 0; i < n; i++) {
        cout << u(xi[i], c) << "  " << pow(xi[i], 2) << endl;

    }
    return 0;
}

