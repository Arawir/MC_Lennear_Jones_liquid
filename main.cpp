#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;
////////////////////////////////////////////////////////////////////////////////
#define N_max 10000  

double ro = 0.844;
double sigma = 1.0;
double sigma2 = sigma*sigma;
double eps = 1.0;
double L = 20;
double N = (int)(L*L*ro);
double Z = 2.5;
double Z2 = 6.25;
double A_r[N_max][2];
double step = 0.2;
double T_red = 0.71;

double h_step = 0.001;
int h_tab[100000];
bool hist_active=0;
long long int h_num=0;

fstream plik_1;
fstream plik_2; //linie zewnetrzne pudelka
fstream plik_3;
////////////////////////////////////////////////////////////////////////////////
double d_rand(){
  return (double)rand()/(double)RAND_MAX;
}
void Draw_box(){
  plik_2.open("2.txt",ios::out);
  plik_2 << "-0.5 -0.5" << endl
	 << "-0.5 " << L+0.5 << endl
	 << L+0.5 << " " << L+0.5 << endl
	 << L+0.5 << " -0.5" << endl
	 << "-0.5 -0.5";
  plik_2.close();
}

void Init(){
  srand( time(NULL) );
  Draw_box();
  for(int i=0; i<1000; i++) h_tab[i]=0; // zerujemy histogram
}

void Close(){

}

void Generate_0(){
  for(int i=0; i<N; i++){
    A_r[i][0] = (double)(i/(int)L);
    A_r[i][1] = (double)(i%(int)L);
  }
}

void Save_picture(){
  plik_1.open("1.txt", ios::out);
  for(int i=0; i<N; i++){
    plik_1 << A_r[i][0] << " " << A_r[i][1] << endl;
  }
  plik_1.close();
}

double R2(double x1, double y1, double x2, double y2){
  double dx = min( abs(x1-x2) , L-abs(x1-x2) );
  double dy = min( abs(y1-y2) , L-abs(y1-y2) );
  return  dx*dx+dy*dy;
}

double dU(double x0, double y0, double xn, double yn, int num){
  double U0=0, Un=0, rn, r0, tmp;

  for(int i=0; i<N; i++){
    if(i==num) continue;
    r0 = R2(x0, y0, A_r[i][0], A_r[i][1]);
    rn = R2(xn, yn, A_r[i][0], A_r[i][1]);

    if(hist_active){ h_tab[(int)(sqrt(r0)/h_step)]++; }
    if((rn>Z2)&&(r0>Z2)) continue; // obcięcie
    
    tmp = sigma2/r0;
    tmp = tmp*tmp*tmp;
    U0 = U0 + tmp*tmp - tmp;


    tmp = sigma2/rn;
    tmp = tmp*tmp*tmp;
    Un = Un + tmp*tmp - tmp ;
  }
  return 4.0*eps*(Un-U0)/T_red;
}

void MCS(){
  double xn, yn, dE;
  for(int i=0; i<N; i++){
    if(hist_active){ h_num++; }
    xn = A_r[i][0]+(d_rand()-0.5)*step;
    yn = A_r[i][1]+(d_rand()-0.5)*step;
    if(xn>L) xn = xn-L;
    if(xn<0) xn = xn+L;
    if(yn>L) yn = yn-L;
    if(yn<0) yn = yn+L;
    dE = dU(A_r[i][0], A_r[i][1], xn, yn, i);
    if(dE<0){ A_r[i][0] = xn; A_r[i][1] = yn; continue; }
    if(d_rand()<=exp(-dE)){ A_r[i][0] = xn; A_r[i][1] = yn; continue; }
    
  }
}

void Write_data(){
  cout << "////////////////////////////////////////" << endl
       << "  T_red = " << T_red << endl
       << "  Ro = " << ro << endl
       << "  Sigma = " << sigma << endl
       << "  Epsilon = " << eps << endl
       << "  R_cut = " << Z << endl
       << "  Step = " << step << endl
       << "  L = " << L << endl
       << "  N = " << N << endl
       << "  HIST_ACTIVE = " << hist_active << endl
       << "////////////////////////////////////////" << endl;
}

void Save_hist(){
  int I_max = (int)(L/2.0/h_step);
  double P=0;
  cout << h_num << endl;
  plik_3.open("3.txt",ios::out);
  for(int i=0; i<I_max; i++){
    P = M_PI*(h_step*(double)i*h_step*(double)i - h_step*(double)(i-1)*h_step*(double)(i-1));
    plik_3 << (double)i*h_step << "   " <<(double)h_tab[i]/P/(double)h_num/ro << endl;
  }
  plik_3.close();
}


////////////////////////////////////////////////////////////////////////////////
int main(){
  char a='0';
  Init();
  Generate_0();

  while(1){
    cout << "///////////////////////////////////////////////a0|s2|d3|f4|z|q|t 123|l 123"<< endl;
    cin >> a;
    if(a=='a'){
      for(int i=0; i<1; i++){ MCS(); }
      a='0';
    }
    if(a=='s'){
      for(int i=0; i<100; i++){ MCS(); }
      a='0';
    }
    if(a=='d'){
      for(int j=0; j<100; j++){
	for(int i=0; i<10; i++){ MCS(); }
	cout << j << endl;
      }
      a='0';
    }
    if(a=='f'){
      for(int j=0; j<100; j++){
	for(int i=0; i<100; i++){ MCS(); }
	cout << j << endl;
      }
      a='0';
    }
    if(a=='t'){
      cin >> T_red;
    }
    if(a=='l'){
      cin >> L;
      ro = N/L/L;
      Draw_box();
    }
    if(a=='r'){
      cin >> ro;
      L = sqrt(N/ro);
      Draw_box();
    }
    if(a=='w'){
      Write_data();
    }
    if(a=='h'){
      if(hist_active==1) hist_active =0;
      else  hist_active =1;
    }
    if(a=='y'){
      cin >> step;
    }
    if(a=='j'){
      for(int i=0; i<100000; i++) h_tab[i]=0;
      h_num=0;
    }
    if(a=='z'){ Save_picture(); a='0';}
    if(a=='x'){ Save_hist(); a='0';}
    if(a=='q'){ break; }

  }

  Close();
  return 0;
}
