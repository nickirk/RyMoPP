#include <stdio.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <limits>
#include <math.h>
#include <string>
#include <sstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <complex>
#include <iomanip>  //for setprecision
#include <sys/stat.h>
#include <sys/types.h>
#include "integral.h"
#include "line.h"
using namespace std;
using namespace Eigen;
#define pi 3.141592653589793
//

MatrixXd InitHam0(vector<vector<int>>& WaveFuncsBasis, MatrixXd mat, int& N, int& DN, int& DN1,int& DL, bool FMagneticNumber){
/*Initialise H0. Read the eigen energies of the basis wave functions and put them into the diagonal
entries of H0.*/
  stringstream stream;
  double j=0;
  double energy0=0.;
  string name;
  int size=WaveFuncsBasis.size();
//Start the loop for reading wave function files
  for (int i=0; i<size; i++){
    //reading out paramters of wave function
    int n=WaveFuncsBasis[i][0];
    int l=WaveFuncsBasis[i][1];
    int m=WaveFuncsBasis[i][2];
    int dj=WaveFuncsBasis[i][3]; 
    if (dj<0) j=l-0.5; 
    else j=l+0.5;
    stream.str(string());
    stream << fixed << setprecision(1) << j;
    string J=stream.str();
    name="../Wave/wfkt_X_Rb87_n=" + to_string(n) +"_L=" + to_string(l) +"_J=" + J +"_step=0.001.txt";
    ifstream file;
    file.open(name);
    if (!file.good()){
      name="../Wave/wfkt_X_Rb87_n=" + to_string(n) +"_L=" + to_string(l) +"_J=0_step=0.001.txt";
      file.open(name);
      if (!file.good()){
        cout << "Reading " << name << " failed! Check if file exists or not!" << endl; 
        cout<< "Unable to open file! Reading Diagonal entries failed! Programme Stops!";
        abort();
      }
    }
    //eigen energy is put at the end of each wave function file, get the last line of the file
    string line=getLastLine(file);
    double energy=stod(line);
    if (n==N && l==0) energy0=energy;// Shift the whole spectra w.r.t E0, which is the eigen energy of the central eigen state.
    cout << "i= " << i<< endl;
    cout << "n,l,m,dj= " << n << ", "<<l<< ", " <<m <<", "<< dj << endl;
    cout << "energy= " << energy << endl;
    cout << "energy0= " << energy0 << endl;
    mat(i,i)=energy;
  }
  for (int t=0; t<size; t++){
    cout << "Before subtraction, H0(t,t)= " << mat(t,t) << endl;	
    mat(t,t)=mat(t,t)-energy0;	
    cout << "After subtraction, H0(t,t)= " << mat(t,t) << endl;	
  }
  return mat;
}

int findPosition(MatrixXd V, double x){
//This function is used to find the position of a value in a vector.
  int i;
  for (  i=0; i<V.rows(); i++){
    if (abs(x-V.col(0)(i))<=0.0001){
      return i;
    }
  }
  return i-1;
}

Vector2cd WaveFuncAngS(double theta, double phi, int n, int l, int m, int dj){
//Calculating the angular part of S wave scattering
  double ReY0, ImY0, ReY1, ImY1,mj;
  //Define the real and imaginary part of the two spherical
  //harmonics in a spinor
  double coef0=0.,coef1=0.;
  int ml,ml_1;
  Vector2cd Ang;
 //converting a iteration index m into m_l
  if((-l+m-1)<=-l) {
    ml=-l;
  }
  else ml=-l+m-1;
  if((-l+m)>=l) {
    ml_1=l;
  }
  else ml_1=-l+m;
  //calculate m_j from iteration index m
  mj= -(l+1)+m+1./2;
  SpherY(l, ml, theta, phi, ReY0, ImY0);
  SpherY(l, ml_1, theta, phi, ReY1, ImY1);
  complex<double> Y0(ReY0,ImY0);
  complex<double> Y1(ReY1,ImY1);

  if (l==0){
    coef0=sqrt((l+1./2+mj)/(2*(l+1./2))) ;
    coef1=sqrt((l+1./2-mj)/(2*(l+1./2))) ;
  }
  else if(dj>0) {
    coef0=sqrt((l+1./2+mj)/(2*(l+1./2))) ;
    coef1=sqrt((l+1./2-mj)/(2*(l+1./2))) ;
  }
  else {
    coef0=-sqrt((l-1./2-mj+1)/(2*(l-1./2)+2)) ;
    coef1=sqrt((l-1./2+mj+1)/(2*(l-1./2)+2)) ;
  }

  Ang(0)=coef0*Y0;
  Ang(1)=coef1*Y1;
  return Ang;
}

Vector2cd WaveFuncAngSNoSin(double theta, double phi, int n, int l, int m, int dj){
//Wave function without sin\theta. For calculating PWS3.
  double ReY0=0., ImY0=0., ReY1=0., ImY1=0.,mj;
  double coef0=0.,coef1=0.;
  int ml,ml_1;

  Vector2cd Ang;
  if ((-l+m-1)<=-l) {
    ml=-l;
  }
  else ml=-l+m-1;
  if ((-l+m)>=l) {
    ml_1=l;
  }
  else ml_1=-l+m;
  mj= -(l+1)+m+1./2;
  if(abs(ml_1)==1){
    SpherYNoSin(l, ml_1, theta, phi, ReY1, ImY1);
  }
  if (abs(ml)==1){
    SpherYNoSin(l, ml, theta, phi, ReY0, ImY0);
  }

  complex<double> Y0(ReY0,ImY0);
  complex<double> Y1(ReY1,ImY1);
  if (l==0){
    coef0=sqrt((l+1./2+mj)/(2*(l+1./2))) ;
    coef1=sqrt((l+1./2-mj)/(2*(l+1./2))) ;
  }// dj>0 means j=l+1./2
  else if(dj>0) {
    coef0=sqrt((l+1./2+mj)/(2*(l+1./2))) ;
    coef1=sqrt((l+1./2-mj)/(2*(l+1./2))) ;
  }
  else {// j=l-1./2
    coef0=-sqrt((l-1./2-mj+1)/(2*(l-1./2)+2)) ;
    coef1=sqrt((l-1./2+mj+1)/(2*(l-1./2)+2)) ;
  }
  Ang(0)=coef0*Y0;
  Ang(1)=coef1*Y1;
  return Ang;
}

MatrixXd ReadR(string NameofFile){
//Read the radial part of the wave function into a matrix, including the
//x and the X(x).
  string line;
  int row,col;
  ifstream tFile;
  tFile.open(NameofFile);
  int i;
  for ( i = 0; getline(tFile, line); i++)
  ;
  tFile.close();
  ifstream pFile;
  pFile.open(NameofFile);
  int maxRow=i-2;//becuase the lastline of the file is not wavefunction, but the energy. Ignore it.
  MatrixXd   mat(maxRow,2);
  if (pFile.is_open()){
    row=0;      // getline() returns the stream. testing the stream with while returns error such as EOF.
    while(getline(pFile, line) && row<maxRow){
      // so here we know that the read was a success and that line has valid data
      stringstream ss(line);
      col=0;
      double value;
      while(ss >> value){
       mat(row,col)=value;
       col++;
      }
      row++;
    }
    pFile.close();
  }
  else{
    cout << "Unable to open file. Read radial part failed, programme Stops.";
    abort();
  }
  return mat;
}

MatrixXd ReadH(string NameofFile){
//This is designed for quick testing. There is a option in main() to
//write out the full Hamiltonian into a file during a previous successful run,
//so if it is enabled, the programme will read the H.txt file instead of 
//calculating the full Hamiltonian again, which is very costy.
  string line;
  int row,col;
  ifstream tFile;
  tFile.open(NameofFile);
  int i;
  for ( i = 0; getline(tFile, line); i++)
  ;
  tFile.close();
  ifstream pFile;
  pFile.open(NameofFile);
  int maxRow=i;//becuase the lastline of the file is not wavefunction, but the energy. Ignore it.
  MatrixXd   mat(maxRow,maxRow);
  if (pFile.is_open()){
    row=0;// getline() returns the stream. testing the stream with while returns error such as EOF.
    while(getline(pFile, line) && row<maxRow){
     // so here we know that the read was a success and that line has valid data
      stringstream ss(line);
      col=0;
      double value;
      while(ss >> value){
        mat(row,col)=value;
        col++;
      }
      row++;
    }
    pFile.close();
  }
  else{
    cout << "Unable to open file. Read H failed, programme Stops.";
    abort();
  }
  return mat;
}

double AngularS(double theta, double phi, double x, int N, int n,int l,int m, int dj,int n1,int l1,int m1,int dj1){
//
  double neff=N-(3.1311804+0.1784/pow((N-3.1311804),2));
  double R=pow(x,2), k=-1./pow(neff,2)+2./R;
  //Fitting paramters for scattering length
  double a= 428.36,  b= -12892.4,c=231179, d = -1.97045e+06, e= 6.57298e+06, f  = 1.82505;
  double a2= 565.637,  b2=-6450.73,c2= 34658.6, d2 = -53720.5, e2= -93072.7, f2  = -17.1526;
  double ASS=0.0;
  double ATS=0.0;
  double xtmp=x;
  while (k<0.){
    xtmp=xtmp-0.001;
    R=pow(xtmp,2);
    k=-1./pow(neff,2)+2./R;
  }
  k=pow(k,1./2);
  ASS=a*k+b*pow(k,2)+c*pow(k,3) +d*pow(k,4) +e*pow(k,5)+f;
  ATS=a2*k+b2*pow(k,2)+c2*pow(k,3) +d2*pow(k,4) +e2*pow(k,5)+f2;
  bool ST = true;
  Vector2cd Ang0=WaveFuncAngS(theta, phi, n, l, m, dj);
  Vector2cd Ang1=WaveFuncAngS(theta, phi, n1, l1, m1, dj1);
  //If Singlet and Triplet option is on, calculate it
  if (ST){
    Ang1(0)=ATS*Ang1(0);
    Ang1(1)=0.5*(ASS+ATS)*Ang1(1);
  }
  else{
    Ang1=Ang1*ATS;
  }

  complex <double> tmp=Ang0.adjoint()*Ang1;//angular part of integral
  return tmp.real();
}

double DRadialDr(MatrixXd V, int i){
//Calculate \frac{dR}{dr}, remember the change of variables:
//x=\sqrt{r}, X(x)=x^{3/2}R(r)
  double step=V(1,0)-V(0,0);
  double dR=0.;
  double dRdr=0.;
  if (i==0) i=1;
  dRdr=-3./4*pow(V(i,0),-7./2)*V(i,1)+1./2*pow(V(i,0),-5./2)*(V(i,1)-V(i-1,1))/step;
  return dRdr;
}

Vector2cd WaveFuncAngPDerivativeTheta(double theta, double phi, int n, int l, int m, int dj){
//Calculate \frac{\partial \Psi}{\partial \theta}
  Vector2cd Ang;
  int ml,ml_1;
  double ReY0=0., ImY0=0., ReY1=0., ImY1=0.,mj;
  double GammaProduct=0., GammaProduct1=0.;
  complex<double> i(0.,1.);
  if((-l+m-1)<=-l) {
    ml=-l;
  }
  else ml=-l+m-1;
  if((-l+m)>=l) {
    ml_1=l;
  }
  else ml_1=-l+m;
  
  if (ml+1<=l){     //calculate Y^{m+1}_l
    SpherY(l, ml+1, theta, phi, ReY0, ImY0);
    GammaProduct=sqrt(double(tgamma(double(l-ml+1)))
      *double(tgamma(double(l+ml+2)))/(double(tgamma(double(l+ml+1)))*double(tgamma(double(l-ml)))));
  }
  complex<double>  Y0_1(ReY0,ImY0);
  
  if (ml_1+1<=l) {
    SpherY(l, ml_1+1, theta, phi, ReY1, ImY1);
    GammaProduct1=sqrt(tgamma(double(l-ml_1+1))
      *tgamma(double(l+ml_1+2))/(tgamma(double(l+ml_1+1))*tgamma(double(l-ml_1))));
  }
  complex<double> Y1_1(ReY1,ImY1);
  
  double coef0=0., coef1=0.;
  mj= -(l+1)+m+1./2;
  if (abs(ml)==1){
    SpherYNoSin(l, ml, theta, phi, ReY0, ImY0);
  }	
  complex<double> Y0(ReY0,ImY0);
  if (abs(ml_1)==1){
    SpherYNoSin(l, ml_1, theta, phi, ReY1, ImY1);
  }		
  complex<double> Y1(ReY1,ImY1);
  if (l==0){
    coef0=sqrt((l+1./2+mj)/(2*(l+1./2))) ;
    coef1=sqrt((l+1./2-mj)/(2*(l+1./2)));
  }// dj>0 means j=l+1./2
  else if(dj>0) {
    coef0=sqrt((l+1./2+mj)/(2*(l+1./2)));
    coef1=sqrt((l+1./2-mj)/(2*(l+1./2)));
  }
  else {// j=l-1./2
    coef0=-sqrt((l-1./2-mj+1)/(2*(l-1./2)+2));
    coef1=sqrt((l-1./2+mj+1)/(2*(l-1./2)+2)) ;
  }

  Ang(0)=(double(ml)*Y0+GammaProduct*Y0_1)*coef0;
  Ang(1)=(double(ml_1)*Y1+GammaProduct1*Y1_1)*coef1;
//	cout << "l, m, m1, ml, ml_1, mj= " << l << "," << m <<  "," << ml << "," << ml_1 << "," << mj << endl;
//	cout << "Y0 ,Y0_1  =" << Y0 << ",  "<< Y0_1 <<endl;
//	cout << "coef0, coef1= " << coef0<< "  " << coef1 << endl;
//	cout << "ml*1./tan(theta)*Y0=  " << double(ml)*Y0 << endl;
//	cout << "exp(-i*phi)=  "  << exp(-i*phi) << endl;
//	cout << "GammaProduct, GammaProduct1= " << GammaProduct <<"  " << GammaProduct1 << endl;
//	cout << "thetaD " << Ang << "\n" << endl;
  return Ang;
}

double AngularP(double theta, double phi, double x,double RadialP_R, double RadialP_dR, int N, int n,int l,int m,
 int dj,int n1,int l1,int m1,int dj1){
  Vector2cd Ang0=WaveFuncAngS(theta, phi, n, l, m, dj);
  Vector2cd Ang1=WaveFuncAngS(theta, phi, n1, l1, m1, dj1);
  Vector2cd PhiDer0=WaveFuncAngSNoSin(theta, phi, n, l, m, dj);
  Vector2cd PhiDer1=WaveFuncAngSNoSin(theta, phi, n1, l1, m1, dj1);
  Vector2cd ThetaDer0=WaveFuncAngPDerivativeTheta(theta, phi,  n, l, m, dj);
  Vector2cd ThetaDer1=WaveFuncAngPDerivativeTheta(theta, phi,  n1, l1, m1, dj1);
  double tmp1=0.;
  double tmp2=0.;
  double tmp3=0.;
  double r=pow(x,2);
  double neff=N-(3.1311804+0.1784/pow((N-3.1311804),2));
  double R=pow(x,2), k=-1./pow(neff,2)+2./R;
  double ASP=0.0;
  double ATP=0.0;
  double xtmp=x;
  //fitting parameters for the scattering length
  double a=-85.141,b=0.296747,c=-0.00291864,d=1.1992e-05,e=-1.66351e-08,f=470.161;
  double a2=-0.00954914,b2=0.266178,c2=-1.69726,d2=9.05673,e2=-27.0818,f2=-1.20294e-05;
  while (k<0.){
    xtmp=xtmp-0.001;
    R=pow(xtmp,2);
    k=-1./pow(neff,2)+2./R;
  }
  	
  k=pow(k,1./2);
  ASP=a/k+b/pow(k,2)+c/pow(k,3) +d/pow(k,4) +e/pow(k,5)+f;
  ATP=1./(a2*k+b2*pow(k,2)+c2*pow(k,3) +d2*pow(k,4) +e2*pow(k,5)+f2);
  int ml, ml_1,ml1,ml1_1;
  bool FST=true;
  if (FST){
    Ang0(0)=ATP*Ang0(0);
    Ang0(1)=1./2*Ang0(1)*(ASP+ATP);
    PhiDer0(0)=ATP*PhiDer0(0);
    PhiDer0(1)=1./2*PhiDer0(1)*(ASP+ATP);
    ThetaDer0(0)=ATP*ThetaDer0(0);
    ThetaDer0(1)=1./2*ThetaDer0(1)*(ASP+ATP);
  }
  if((-l+m-1)<=-l) {
    ml=-l;
  }
  else ml=-l+m-1;
  if((-l+m)>=l) {
    ml_1=l;
  }
  else ml_1=-l+m;
  if((-l1+m1-1)<=-l1) {
    ml1=-l1;
  }
  else ml1=-l1+m1-1;
  if((-l1+m1)>=l1) {
    ml1_1=l1;
  }
  else ml1_1=-l1+m1;
  complex <double> Ctmp1=Ang0.adjoint()*Ang1;
  tmp1=Ctmp1.real();
  complex<double> Ctmp2=ThetaDer0.adjoint()*ThetaDer1;
  tmp2=Ctmp2.real();
  PhiDer0(0)=PhiDer0(0)*double(ml);
  PhiDer0(1)=PhiDer0(1)*double(ml_1);
  PhiDer1(0)=PhiDer1(0)*double(ml1);
  PhiDer1(1)=PhiDer1(1)*double(ml_1);
  complex<double> Ctmp3=PhiDer0.adjoint()*PhiDer1;
  tmp3=Ctmp3.real();
  double total=tmp1*RadialP_dR+(tmp2+tmp3)*RadialP_R/pow(r,2);
  cout << "x= " << x << endl;
  cout << "RadialP_dR= " << RadialP_dR<< endl;
  cout << "RadialP_R/pow(r,2)= " << RadialP_R/pow(r,2)<< endl; 
  cout <<  "tmp1,tmp2,tmp3==" << tmp1<< "," << tmp2<< "," << tmp3 << endl;
  cout << "Total P=  " << total<< endl;
  return total;
}

double Integrate(MatrixXd* WaveFuncs, double x,int i, int j, int N, int n, int l, int m, int dj, int n1, int l1, int m1, int dj1, bool FPWave){
  double I_S=0., I_T=0., I=0.;
  int nt, np;
  int ml, ml_1,ml1,ml1_1;
  double mj,mj1,coef01,coef11,coef0,coef1;
  double ReY0, ImY0, ReY1, ImY1;
  string name0, name1;
  nt = 181; np = 6;
  
  I_S=AngularS(0., 0.,  x,  N,  n, l, m,  dj, n1, l1, m1, dj1);
  if((-l+m-1)<=-l) {
    ml=-l;
  }
  else ml=-l+m-1;
  if((-l+m)>=l) {
    ml_1=l;
  }
  else ml_1=-l+m;
  if((-l1+m1-1)<=-l1) {
    ml1=-l1;
  }
  else ml1=-l1+m1-1;
  if((-l1+m1)>=l1) {
    ml1_1=l1;
  }
  else ml1_1=-l1+m1;
  mj= -(l+1)+m+1./2;
  mj1=-(l1+1)+m1+1./2;
  if (l==0){
    coef0=sqrt((l+1./2+mj)/(2*(l+1./2))) ;
    coef1=sqrt((l+1./2-mj)/(2*(l+1./2))) ;
  }// dj>0 means j=l+1./2
  else if(dj>0) {
    coef0=sqrt((l+1./2+mj)/(2*(l+1./2))) ;
    coef1=sqrt((l+1./2-mj)/(2*(l+1./2))) ;
  }
  else {// j=l-1./2
    coef0=-sqrt((l-1./2-mj+1)/(2*(l-1./2)+2)) ;
    coef1=sqrt((l-1./2+mj+1)/(2*(l-1./2)+2)) ;
  }
  if (l1==0){
    coef01=sqrt((l1+1./2+mj1)/(2*(l1+1./2))) ;
    coef11=sqrt((l1+1./2-mj1)/(2*(l1+1./2))) ;
  }// dj>0 means j=l+1./2
  else if(dj1>0) {
    coef01=sqrt((l1+1./2+mj1)/(2*(l1+1./2))) ;
    coef11=sqrt((l1+1./2-mj1)/(2*(l1+1./2))) ;
  }
  else {// j=l-1./2
    coef01=-sqrt((l1-1./2-mj1+1)/(2*(l1-1./2)+2)) ;
    coef11=sqrt((l1-1./2+mj1+1)/(2*(l1-1./2)+2)) ;
  }
  cout << "In Integrate"<< endl;	
  int id = omp_get_thread_num(); 
  cout << "ID= " << id << endl;
  cout << "x= " << x << endl;
  cout << "l, l1, m, m1=  " << l << ",  " << l1 << ",  " << m << ",  " << m1 << endl;
  cout << "ml, ml_1, ml1, ml1_1  " << ml <<", " << ml_1 <<", " << ml1 << ", " << ml1_1 << endl;
  cout << " mj, mj_1= " << mj << ", " << mj1 << endl;
  
  
  MatrixXd& Radial0=WaveFuncs[i];
  MatrixXd& Radial1=WaveFuncs[j];
  double R0;
  double R1;
  int i0,k0;
  
  i0 = findPosition (Radial0, x);
  R0= Radial0(i0,1);
  k0 = findPosition (Radial1, x);
  R1 = Radial1(k0,1);
  if (FPWave){
    //If P wave scattering is on.
    double RadialP_dR=0.;
    double RadialP_R=0.;
    RadialP_dR=DRadialDr(Radial0,i0)*DRadialDr(Radial1,k0);
    RadialP_R=R0*R1/pow(x,3);//change from X(x) to R(r).
    I_T = AngularP(0.,0.,x,RadialP_R, RadialP_dR,N,n,l,m,dj,n1,l1,m1,dj1);
  }
//	 if (abs(I_S)< 10.e-16) I_S=0.;
  I=2*pi*I_S*R0*R1*pow(x,-3)+I_T*6*pi;
/*In the scattering theory, they adopted the normalization of spherical harmonics
 as <Y_m^l|Y_m'^l'>=4*pi/(2l+1)\delta_{mm'}\delta_{ll'}, but in the code we 
use the normalization <Y_m^l|Y_m'^l'>=\delta_{mm'}\delta_{ll'}, in his paper 
Creation of Polar and Nonpolar Ultra-Long-Range Rydberg Molecules, he added 
an (2l+1)/4pi to normalize the spherical harmonics, but here we are 
already normalized, so we don't have to do that.*/
  return I;
}

MatrixXd OffDiag(vector<vector<int>>& WaveFuncsBasis, MatrixXd* WaveFuncs,
 MatrixXd mat, double x, int N, int DN, int DN1, int DL, bool FMagneticNumber, bool FPWave){
  int maxm=1, maxm1=1;
  double I=0.;
  int size = WaveFuncsBasis.size();
  int i;
  #pragma omp parallel for schedule(dynamic,1) //collapse(2)
  //parallel computing
  for (i=0; i<size; i++){
    int id = omp_get_thread_num(); 
    int n=WaveFuncsBasis[i][0];
    int l=WaveFuncsBasis[i][1];
    int m=WaveFuncsBasis[i][2];
    int dj=WaveFuncsBasis[i][3];
    for (int j=i; j< size;j++){
      int n1=WaveFuncsBasis[j][0];
      int l1=WaveFuncsBasis[j][1];
      int m1=WaveFuncsBasis[j][2];
      int dj1=WaveFuncsBasis[j][3];
      I=Integrate(WaveFuncs, x, i, j, N, n, l, m, dj, n1, l1, m1, dj1, FPWave);
      cout << "In Offdiag" << endl;
      cout << "ID= " << id << endl;
      cout << "x= " << x << endl;
      cout << "(i,j)=" << "("<< i << "," << j <<")" << endl;
      cout << "(n,l,m,dj)=" << n << "," << l <<"," << m << ","<<dj << endl;
      cout << "(n1,l1,m1,dj1)=" <<n1<< ","<< l1 <<","<< m1 <<","<< dj1 << endl;
      cout << "I =             " << I << "\n" << endl;
      mat(i,j)=I;
      mat(j,i)=I;
    }
  }
  return mat;
}

int main(int argc, char const *argv[]) {
  int N=40, DN=4, DN1=1, DL=4,DM=1;//DM should smaller than l// maxm=1;
  int size=0,n,l,j;
  double xmin=22, xmax=52,I;
  double x=xmin, step=0.2;
  int loopnumber= int((xmax-xmin)/step);
  bool FReadH=false, FMagneticNumber=true, FPWave=true;
   /* control the magnetic
   Number choices. In principle this option should always keep as true.
   Otherwise it means there is no spin orbital coupling,
   which doesn't make much sense.*/

  bool FManualSetBasis=true;
  bool FManyFolds=true;
  string name0;
  omp_set_nested(1);
  omp_set_num_threads(3);
  vector<vector<int>> WaveFuncsBasis;
  vector<vector<int>> WaveFuncsBasistmp={};
  WaveFuncsBasistmp={{37,5},{36,5}};
  //Option for manual input some wave functions
  //In the future this can be writting into an interactive part for user to
  //input on the run.
  cout << "Generating Basis Functions And Writing Basis Info Into File... " << endl;
  ofstream myfile1;
  myfile1.open("WaveFuncsBasis");
  if (FManualSetBasis){
    int ltmp, lstart;
    for (int it=0; it < WaveFuncsBasistmp.size(); it++){
      if (FManyFolds){ 
        ltmp=WaveFuncsBasistmp[it][0]-1;
        lstart=WaveFuncsBasistmp[it][1];
      }
      else{
        ltmp=WaveFuncsBasistmp[it][1];
        lstart=0;
      }
      for (int l=lstart; l<=ltmp; l++){
        if (l==0){
          int mtmp=1,djtmp=1;
          vector<int> Basistmp;
          Basistmp.push_back(WaveFuncsBasistmp[it][0]);
          Basistmp.push_back(l);
          Basistmp.push_back(mtmp);
          Basistmp.push_back(djtmp);
          WaveFuncsBasis.push_back(Basistmp);
          myfile1 << WaveFuncsBasistmp[it][0] << " " << l << " " << mtmp << " "<< djtmp << endl;
          size++;
          continue;
        }
        //if (FMagneticNumber) maxm=2*(l+1);
        int DMtmp=DM;
        if (l<DM) DMtmp=l;
        for (int m=l-DMtmp; m<=l+DMtmp+1; m++){
          for (int dj=-1;dj<=1;dj+=2){ 
            vector<int> Basistmp;
            Basistmp.push_back(WaveFuncsBasistmp[it][0]);
            Basistmp.push_back(l);
            Basistmp.push_back(m);
            Basistmp.push_back(dj);
            WaveFuncsBasis.push_back(Basistmp);
            myfile1 << WaveFuncsBasistmp[it][0] << " " << l  << " " << m << " "<< dj << endl;
            size++;
          }
        }
      }
    }
  }
  for (int n=N-DN; n<=N+DN1; n++){
    for (int l=0; l<=DL;l++){
      if (l==0){
        int m=1,dj=1;
	vector<int> pars;
	pars.push_back(n);
	pars.push_back(l);
	pars.push_back(m);
	pars.push_back(dj);
	WaveFuncsBasis.push_back(pars);
	continue;
      }
      //if (FMagneticNumber) maxm=2*(l+1);
      int DMt=DM;
      if (l<DM) DMt=l;
      for (int m=l-DMt; m<=l+DMt+1; m++){
        for (int dj=-1;dj<=1;dj+=2){ 
	  vector<int> pars;
	  pars.push_back(n);
	  pars.push_back(l);
	  pars.push_back(m);
	  pars.push_back(dj);
	  WaveFuncsBasis.push_back(pars);
	  int nt=WaveFuncsBasis[size][0];
	  int lt=WaveFuncsBasis[size][1];
	  int mt=WaveFuncsBasis[size][2];
	  int djt=WaveFuncsBasis[size][3];
	  myfile1  << nt << " " << lt <<" " << mt << " "<<djt << endl;
	  size++;
	}
      }
    }
  }
  cout << "Done" << endl;
  myfile1.close();	
  cout << "Size of Matrix=" << WaveFuncsBasis.size() << endl;
  size=WaveFuncsBasis.size();
  MatrixXd *WaveFuncs = new MatrixXd[size];
  cout << "Starting to load wavefunctions into memory..."<< endl;
  for (int i=0; i<size; i++){
    int n=WaveFuncsBasis[i][0];
    int l=WaveFuncsBasis[i][1];
    int m=WaveFuncsBasis[i][2];
    int dj=WaveFuncsBasis[i][3];
    stringstream stream;
    double j=l+dj*0.5;
    if (l==0) j=0.5;
    stream << fixed << setprecision(1) << j;
    string J=stream.str();
    name0="../Wave/wfkt_X_Rb87_n="+ to_string(n)+"_L="+to_string(l) +"_J="+ J +"_step=0.001.txt";
    ifstream file0;
    file0.open(name0);
    if (!file0.good()){
      name0="../Wave/wfkt_X_Rb87_n=" + to_string(n) +"_L=" + to_string(l) +"_J=0_step=0.001.txt";
      file0.open(name0);
        if (!file0.good()){
      	  cout<< "File can't be found! Fatal error! Programme stops!"<< endl;
      	  return 0;
        }
        else file0.close();
    }
    file0.close();
    WaveFuncs[i]=ReadR(name0);
    cout << i << " wavefunctions loaded" << endl;
  }
  cout << "Done" << endl;
  string directory=to_string(size);
  char const *pchar = directory.c_str();
  struct stat sb;

  if (stat(pchar, &sb) == 0 && S_ISDIR(sb.st_mode)){
    cout << "Directory " << pchar <<" exists! EigenValue files will be written" 
      "into it.\n Be careful not to overwrite previous data in this directory"
      "if you still need them." << endl;
  }
  else{
    if(mkdir(pchar, 0777)==-1){//creating a directory
      cerr<<"Error :  "<<strerror(errno)<<endl;
      exit(1);
    }
  }

  for (int t=0; t<=loopnumber; t+=1){
  //Start the loop of calculating at different x position
    cout << "Calculating at x= " << x << endl;
    MatrixXd H0=MatrixXd::Zero(size,size);
    MatrixXd H1=MatrixXd::Zero(size,size);
    MatrixXd H=MatrixXd::Zero(size,size);
    
    if(FReadH){
    	H=ReadH("H.txt");
    }
    else {
      H0=InitHam0(WaveFuncsBasis, H0,N,DN,DN1,DL,FMagneticNumber);
      ofstream myfileH0;
      myfileH0.open("H0.txt");
      myfileH0 << H0 << endl;
      myfileH0.close();
      H1=OffDiag(WaveFuncsBasis, WaveFuncs, H1, x,  N, DN, DN1,DL, FMagneticNumber, FPWave);
      H=H0+H1;
    //	cout <<  H0 << "\n" <<endl;
    //	cout << H1 << "\n" << endl;
    //	cout << H << endl;
      ofstream myfileH;
      myfileH.open("H.txt");
      myfileH << H << endl;
      myfileH.close();
    }
    //use Solvers from lib Eigen to diagonalise matrix
    SelfAdjointEigenSolver<MatrixXd> eigensolver(H);
    if (eigensolver.info() != Success){
      cout << "Unable to diagonalize Hamiltonian!\n";
      abort();
    }
    else{
      cout << "Success in diagonalizing Hamiltonian,\n Writing EigenValues to file... \n";
      ofstream myfile;
      myfile.open ("./"+ to_string(size) + "/EigenValue"+to_string(x)+".txt");
      myfile <<setprecision(15)<< eigensolver.eigenvalues()<< endl;
      myfile.close();
      cout << setprecision(15)<< eigensolver.eigenvalues()<< endl;
      cout << "Writing into file succeeded!" << endl;
    }
    x+=step;
  }
  return 0;
}
