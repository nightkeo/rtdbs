#include <windows.h>

#include "RayTracing.h"
#include <iostream>
#include <iomanip>
#include <complex>
#include "nr.h"
#include <limits>
#include <ctime>
#include "spline.h"

#define empty -100
//const double pi=3.1415926535897932384626433832795;
#define pi 3.1415926535897932384626433832795
#define m_e 0.910938215454545E-27
#define zar 4.8E-10
#define c 2.99792458E10

//**********************************************************************
//**********************************************************************
//**********************************************************************
class coord2{     //создадим класс координат описывающих параметры магнитных поверхностей

	///**************************************************************************************
	//    (R-R0)/rm=alpha*rho*cos(teta)-delta-gamma*sin^2(teta)
	//               -br*sin^2(2*teta)*cos(teta)
	//    (Z-Z0)/rm=alpha*ell*sin(teta)-bz*sin(2*teta)
	//    definitions:
	//    (R0,Z0) - magnetic axis position in meters
	//    rm      - minor radius in meters
	//    rho(Psi)={Rmax(Psi)-Rmin(Psi)}/2
	//    dimensionless functions:
	//    delta(rho),gamma(rho),ell(rho),br(rho),bz(rho)
	//**************************************************************************************

public :
	double shf[5],ell[5],trin[5],bx[5],bz[5],dpsi[5];
    double rm,x0,z0,x0_real,z0_real;

//______________________________________________________________________
	void set_pol(double *tmp1,double *tmp2,double *tmp3,double *tmp4,double *tmp5,double *tmp6)
	{
		int i;
		for(i=0;i<5;i++)shf[i]=tmp1[i];     //задаем шафранов сдвиг (как функция-полином от rho)
		for(i=0;i<5;i++) ell[i]=tmp2[i];	//Задаем элиптичность	(как функция-полином от rho)
		for(i=0;i<5;i++) trin[i]=tmp3[i];	//Задаем треугольность	(как функция-полином от rho)
		for(i=0;i<5;i++) bx[i]=tmp4[i];		//Задаем bx	(как функция-полином от rho)
		for(i=0;i<5;i++)bz[i]=tmp5[i];		//Задаем bz (как функция-полином от rho)
		for(i=0;i<5;i++)dpsi[i]=tmp6[i];	//Задаем произв. полоид потока	(как функция-полином от rho)
    }
	void set_R0(double &tmp1,double &tmp2,double &tmp3)
	{
	  rm=tmp1;				//Задаём коэффицент масштабирования
	  x0_real=tmp2;			//Задаем центр x0 магн оси в координатах токамака
	  z0_real=tmp3;			//Задаем центр z0 магн оси в координатах токамака
	  x0=tmp2/rm;			//нормируем x0 центр на коэффицент масштабирования
	  z0=tmp3/rm;			//нормируем z0 центр на коэффицент масштабирования
	}

//______________________________________________________________________
//возращаем значения(shf ell trin bx bz) как функциии rho

	inline double ret_shf(double rho)
	{
       return  (((shf[4]*rho+shf[3])*rho+shf[2])*rho+shf[1])*rho+shf[0];
	}
    inline double ret_ell(double rho)
	{
       return  (((ell[4]*rho+ell[3])*rho+ell[2])*rho+ell[1])*rho+ell[0];
	}

	inline double ret_trin(double rho)
	{
       return  (((trin[4]*rho+trin[3])*rho+trin[2])*rho+trin[1])*rho+trin[0];
	}

	inline double ret_bx(double rho)
	{
       return  (((bx[4]*rho+bx[3])*rho+bx[2])*rho+bx[1])*rho+bx[0];
	}
    inline double ret_bz(double rho)
	{
       return  (((bz[4]*rho+bz[3])*rho+bz[2])*rho+bz[1])*rho+bz[0];
	}

//______________________________________________
//возращаем производные значений(shf ell trin bx bz) как функциии rho

	inline double ret_dshf(double rho)
	{
       return  ((4.0*shf[4]*rho+3.0*shf[3])*rho+2.0*shf[2])*rho+shf[1];
	}
    inline double ret_dell(double rho)
	{
       return  ((4.0*ell[4]*rho+3.0*ell[3])*rho+2.0*ell[2])*rho+ell[1];
	}

	inline double ret_dtrin(double rho)
	{
       return  ((4.0*trin[4]*rho+3.0*trin[3])*rho+2.0*trin[2])*rho+trin[1];
	}

	inline double ret_dbx(double rho)
	{
       return  ((4.0*bx[4]*rho+3.0*bx[3])*rho+2.0*bx[2])*rho+bx[1];
	}
    inline double ret_dbz(double rho)
	{
       return  ((4.0*bz[4]*rho+3.0*bz[3])*rho+2.0*bz[2])*rho+bz[1];
	}

//______________________________________________
//возращаем вторые производные значений(shf ell trin bx bz) как функциии rho

	inline double ret_ddshf(double rho)
	{
       return  (12.0*shf[4]*rho+6.0*shf[3])*rho+2.0*shf[2];


	}
    inline double ret_ddell(double rho)
	{
       return  (12.0*ell[4]*rho+6.0*ell[3])*rho+2.0*ell[2];


	}

	inline double ret_ddtrin(double rho)
	{
       return  (12.0*trin[4]*rho+6.0*trin[3])*rho+2.0*trin[2];
	}

	inline double ret_ddbx(double rho)
	{
       return  (12.0*bx[4]*rho+6.0*bx[3])*rho+2.0*bx[2];


	}
    inline double ret_ddbz(double rho)
	{
       return  (12.0*bz[4]*rho+6.0*bz[3])*rho+2.0*bz[2];


	}

//______________________________________________
//возращаем значения первой и второй производной psi как функциии rho

	inline double ret_dpsi(double rho)
	{
		return (((dpsi[4]*rho+dpsi[3])*rho+dpsi[2])*rho+dpsi[1])*rho+dpsi[0];
	}

	inline double ret_ddpsi(double rho)
	{
       return  ((4.0*dpsi[4]*rho+3.0*dpsi[3])*rho+2.0*dpsi[2])*rho+dpsi[1];
	}
 }; //на этом форм. класса кооридант завершено.

 //**********************************************************************
 //**********************************************************************
 //**********************************************************************
 //используем namespace чтобы не использовать теже переменные как и в других файлах

namespace stdray{

bool flag3;
bool flag33;
bool flag;
bool flag2;

bool ret; //если 0 Ray tracing не выполнился из-за ошибки

ofstream out;
ofstream cutoff;
ofstream out3;
ofstream ray;
ifstream initial; //начальные условия для траектории xs,ys,zs,nxs,nys,nzs,F (читаем только F)
ofstream raylog;
ifstream btor;
ofstream ne_log;

int inversT; //Инверсия времени
int moda; //=0 or 1  мода волны(нормальная или ненормальная)
double alpha;
double F;  // F=2*pi*f  f(Гц)
double Btor_0; // в системе си! при R0=Rcenter смотри g файл!
double R_center;
//double plotn_a;              //чтобы н границе не была большая или маленькая плотность
//double plotn0,W_s,W_delt;   //описание плотности !!! в сгс! (чч/cm^2)
double rho_temp, n_temp;
vector<double> rho_array,ne_array;
double n1y0;
double n1z0;
double w_start,theta_start;
double signBpol=1.0;

coord2 cor;
}

using namespace stdray;

 //**********************************************************************
 //**********************************************************************
 //**********************************************************************

 tk::spline ne_spline;


bool DlgFormRayTracing ()  //при нажатие кнопки вызывается
{
	std::cout << "alive" << std::endl;
	int ncoeff;
	double rm,x0,z0;  //double E-38..E38
	double shafr[5],ellip[5],triang[5],bx_m[5],bz_m[5],dpsif[5],*p;
	double zero;



	time_t starttime,endtime;	// Посмотрим сколько мы времени тратим на ray tracing
	starttime=time(NULL);		//

	ret=1;  // Если этот указатель обнулится значала ray tracing не выполнился из-за ошибки в derivs!!!

	flag3=1;
	flag33=1;
	flag=0;
	flag2=0;

	inversT=empty; //Инверсия времени
	moda=empty; //=0 or 1  мода волны(нормальная или ненормальная)
	//**************************************************************************************
	alpha=empty;
	F=empty;  // F=2*pi*f  f(Гц)

	//plotn_a=empty; //чтобы на границе не была большая или маленькая плотность

	//plotn0=empty,W_s=empty,W_delt=empty;   //описание плотности !!! в сгс! (чч/cm^2)

	n1y0=empty;
	n1z0=empty;

	w_start=empty,theta_start=empty;

	out.open("out1.txt");
	cutoff.open("cutoff.txt");
	out3.open("out3.txt");
	ray.open("ray_kuk.txt");
	initial.open("initial_data.txt");
	raylog.open("raytracing.log");
	btor.open("Btor.txt");
	ne_log.open("ne_log_kuk.txt");



	if(!initial){
		//MessageBox(hwnd,"don't exist initial_data.txt! Please press OK","Bad news",0);
		return 0;
	}

	///****
	//1) создадим экземпляр координат опис магн поверхн. коэффиценты полиномов возьмём из файла input_test.dat
	//Аппроксимация psi берется из файла  input_test2.dat

	ifstream in("input_test.dat");  //коэффиценты аппроксимации магнитной конфигурации
	ifstream in2("input_test2.dat"); //коэффиценты аппроксимации dpsi_dro + w_start и rho start
	ifstream in3("input_test3.dat"); //коээфиценты аппроксимации плотности

	if(!in){
		//MessageBox(hwnd,"don't exist in.dat! Please press OK","Bad news",0);
		return 0;}
	if(!in2){
		//MessageBox(hwnd,"don't exist in2.dat! Please press OK","Bad news",0);
		return 0;}
	if(!in3){
		//MessageBox(hwnd,"don't exist in3.dat! Please press OK","Bad news",0);
		return 0;}

	//________________________________________________________
	in >> ncoeff;                                                     //read ncoeff,rm,x0,z0
	in>> rm >>x0>>z0;
	//________________________________________________________
	p=shafr;
	in>> *p ; p++; in>> *p; p++; in>> *p; p++; in>> *p; p++; in>> *p; //read coeffs poilinom
	p=ellip;
	in>> *p ; p++; in>> *p; p++; in>> *p; p++; in>> *p; p++; in>> *p;
	p=triang;
	in>> *p ; p++; in>> *p; p++; in>> *p; p++; in>> *p; p++; in>> *p;
	p=bx_m;
	in>> *p ; p++; in>> *p; p++; in>> *p; p++; in>> *p; p++; in>> *p;
	p=bz_m;
	in>> *p ; p++; in>> *p; p++; in>> *p; p++; in>> *p; p++; in>> *p;
	p=dpsif;
	in2>> *p ; p++; in2>> *p; p++; in2>> *p; p++; in2>> *p; p++; in2>> *p; //dpsi/drho

	in2>>w_start; in2>>theta_start;

	in2>>zero;in2>>n1y0; in2>>n1z0;


	//if(alpha==empty){MessageBox(hwnd,"ERROR could not read alpha(mashtabn koef!!)!! Please press OK","Bad news",0);return 0;}
	if(w_start==empty){raylog<<"error w_start ne zadana!!!";return 0;}
	if(theta_start==empty){raylog<<"error theta ne zadana!!!";return 0;}
	if(n1y0==empty){raylog<<"error ny_start ne zadana!!!";return 0;}
	if(n1z0==empty){raylog<<"error nz_start ne zadana!!!";return 0;}

    while (in3 >> rho_temp && in3 >> n_temp)
    {
        rho_array.push_back(rho_temp);
        ne_array.push_back(n_temp*1e13);
    }

	//tk::spline ne_spline_deb(rho_array, ne_array, tk::spline::linear);
	ne_spline = tk::spline(rho_array, ne_array, tk::spline::linear);

	//in3>>zero;
	//in3>>plotn0;
	//in3>>plotn_a;
	//in3>>W_delt;
	//in3>>W_s;
	//if(plotn0==empty){raylog<<"plotn0 ne zadana!!!";return 0;}
	//if(plotn_a==empty){raylog<<"plotn_a ne zadana!!!";return 0;}
	//if(W_delt==empty){raylog<<"W_delt ne zadana!!!!";return 0;}
	//if(W_s==empty){raylog<<"W_s ne zadana!!!!";return 0;}

	double skip;

	initial>>alpha; //читаем масштабный коэффицент

	initial>>skip>>skip>>skip>>skip>>skip>>skip>>F>>moda>>inversT;

	F=2*pi*F;   //F=2*pi*f  f(Гц)

	if(F==empty){raylog<<"F ne zadana!!!!";return 0;}
	if(moda==empty){raylog<<"moda ne zadana!!!!";return 0;}
	if(inversT==empty){raylog<<"inversT ne zadano!!!!";return 0;}

	//

	btor>>Btor_0;
	btor >> R_center;

	cout <<"B_tor: " << Btor_0 << endl;
	cout << "R_center: " << R_center << endl;



	//________________________________________________________

	const double k0=F/c;

	//out << ncoeff <<"-ncoeff "<<rm<<"-rm (sm) "<<x0<<"-x0 (sm) "<<z0<<"-z0 (sm)"<< endl;  //write param

	out<<"2*pi*chastota = "<<F<<endl;
	out<<"k0 (1/sm) = "<<F/c<<endl;
	out<<"nomer modi = "<<moda<<endl;

	//out<<"plotn na granice "<<plotn_a<<endl;

	rm=(rm/100)*alpha;                //sm->m  alpha=1.217 - масштабный коэффицент
	x0=x0/100;
	z0=z0/100;

	cor.set_pol(shafr,ellip,triang,bx_m,bz_m,dpsif);      //set polinom and parametres
	cor.set_R0(rm,x0,z0);

	// начнём всё высчитывать для следущей точки
	//сначала считаем всё в текущ точке

	double W[3]={w_start,theta_start,0.0};     //{1.1913f,0.0f,0.0f}; x=0.62  стенка

	//W=[rho,theta,phi]

	double shf,ell,trin,bx,bz,dshf,dell,dtrin,dbx,dbz,ddshf,ddell,ddtrin,ddbx,ddbz;
	double x,y,z,x_real,y_real,z_real,R;
	double sinth,costh,sin2th,cos2th,sin4th,cos4th,sinph,cosph;
	double dx_W0,dx_W1,dx_W2,dy_W0,dy_W1,dy_W2,dz_W0,dz_W1,dz_W2,dR_W0,dR_W1;

	double ddx_W0W0,ddx_W1W1,ddx_W2W2,ddx_W0W1,ddx_W0W2,ddx_W1W2;
	double ddy_W0W0,ddy_W1W1,ddy_W2W2,ddy_W0W1,ddy_W0W2,ddy_W1W2;
	double ddz_W0W0,ddz_W1W1,ddz_W2W2,ddz_W0W1,ddz_W0W2,ddz_W1W2;

	double g00kov,g11kov,g22kov,g01kov,g02kov,g12kov,g;

	double g00,g11,g22,g01,g02,g12;
	double NX_0,NY_0,NZ_0;

	double dpsi,Bpol,B2,B,Btor,R0=R_center/cor.rm,plotn,Bx,By,Bz; // R0=Rcenter=0.36 пока не знаю почему!
	 //незабыть читать Btor_0 из файла!!
	double U,V;

	cout<< "R_0: " << R0 << endl;
	cout<< "cor.rm: " << cor.rm << endl;

	shf=cor.ret_shf(W[0]);
	ell=cor.ret_ell(W[0]);
	trin=cor.ret_trin(W[0]);
	bx=cor.ret_bx(W[0]);
	bz=cor.ret_bz(W[0]);

	dshf=cor.ret_dshf(W[0]);
	dell=cor.ret_dell(W[0]);
	dtrin=cor.ret_dtrin(W[0]);
	dbx=cor.ret_dbx(W[0]);
	dbz=cor.ret_dbz(W[0]);

	ddshf=cor.ret_ddshf(W[0]);
	ddell=cor.ret_ddell(W[0]);
	ddtrin=cor.ret_ddtrin(W[0]);
	ddbx=cor.ret_ddbx(W[0]);
	ddbz=cor.ret_ddbz(W[0]);

	sinth=sin(W[1]);  //Чтоб постоянно не считать посчитаем все синусы-косинусы
	sin2th=sin(2*W[1]);
	sin4th=sin(4*W[1]);
	costh=cos(W[1]);
	cos2th=cos(2*W[1]);
	cos4th=cos(4*W[1]);
	sinph=sin(W[2]);
	cosph=cos(W[2]);




	//x,y,z нормированные

	R= (cor.x0+            W[0]*costh  -shf
									 -trin*sinth*sinth
									 -bx*sin2th*sin2th*costh);


	x=R*cosph;




	y=R*sinph;


	z=cor.z0+            ell*sinth -bz  *sin2th    ;


	//x,y,z физические

			x_real= cor.rm*x;
			y_real= cor.rm*y;
			z_real= cor.rm*z;






	dR_W0=costh  -dshf
									 -dtrin*  sinth*sinth
									 -dbx   *sin2th*sin2th*costh   ;
	dR_W1=
									  -W[0]*sinth
									 -trin* sin2th
									 -bx*2.0f*sin4th*costh
									 +bx*sinth*sin2th*sin2th;






	//первые производные по координатам
	dx_W0 =        (    costh  -dshf
									 -dtrin*  sinth*sinth
									 -dbx   *sin2th*sin2th*costh
					) *cosph;

	dx_W1= (
									  -W[0]*sinth
									 -trin* sin2th
									 -bx*2.0*sin4th*costh
									 +bx*sinth*sin2th*sin2th
			)  *cosph;


	dx_W2= -(cor.x0+           W[0]*costh  -shf
									 -trin*  sinth*sinth
									 -bx   *sin2th*sin2th*costh

			)*sinph;

	dy_W0=  (     costh  -dshf
									 -dtrin*  sinth*sinth
									 -dbx   *sin2th*sin2th*costh
			) *sinph;


	dy_W1=(            -W[0]*sinth
									 -trin* sin2th
									 -bx*2.0*sin4th*costh
									 +bx*sinth*sin2th*sin2th
		  )  *sinph;


	dy_W2=(cor.x0+
									  W[0]*costh  -shf
									 -trin*  sinth*sinth
									 -bx   *sin2th*sin2th*costh

		  )*cosph;
	dz_W0= (      dell* sinth
				  -dbz  *sin2th
		   );

	dz_W1=(
									ell*costh
									-bz*2.0*cos2th
		   );



	dz_W2=0.0;

	//______________________________________________
	   //вторые производные по координатам


	ddx_W0W0= (    -ddshf
									 -ddtrin*  sinth*sinth
									 -ddbx   *sin2th*sin2th*costh
			  ) *cosph;

	ddx_W1W1=(                  -W[0]*costh
									 -2.0*trin* cos2th
									 -bx*(
														  +8*cos4th*costh
														  -4*sin4th*sinth
														  -sin2th*sin2th*costh
										 )
			)  *cosph;
	ddx_W2W2=-(cor.x0+            W[0]*costh  -shf
									 -trin*  sinth*sinth
									 -bx   *sin2th*sin2th*costh

			  )*cosph;
	ddx_W0W1=(                       -sinth
									 -dtrin* sin2th
									 -dbx*2.0*sin4th*costh
									 +dbx*sinth*sin2th*sin2th
			 )  *cosph;

	ddx_W0W2=  -(    costh  -dshf
									 -dtrin*  sinth*sinth
									 -dbx   *sin2th*sin2th*costh
				) *sinph;


	ddx_W1W2=-(
									  -W[0]*sinth
									 -trin* sin2th
									 -bx*2.0*sin4th*costh
									 +bx*sinth*sin2th*sin2th
			  )  *sinph;

	ddy_W0W0= (    -ddshf
									 -ddtrin*  sinth*sinth
									 -ddbx   *sin2th*sin2th*costh
			  ) *sinph;


	ddy_W1W1=(
									  -W[0]*costh
									 -2.0*trin* cos2th
									 -bx*(
														  +8.0*cos4th*costh
														  -4.0*sin4th*sinth
														  -sin2th*sin2th*costh


													)
			)  *sinph;

	ddy_W2W2=-(cor.x0+            W[0]*costh  -shf
									 -trin*  sinth*sinth
									 -bx   *sin2th*sin2th*costh

			 )*sinph;
	ddy_W0W1=(
									  -sinth
									 -dtrin* sin2th
									 -dbx*2.0*sin4th*costh
									 +dbx*sinth*sin2th*sin2th
			 )  *sinph;

	ddy_W0W2= (    costh  -dshf
									 -dtrin*  sinth*sinth
									 -dbx   *sin2th*sin2th*costh
			  ) *cosph;


	ddy_W1W2=(
									  -W[0]*sinth
									 -trin* sin2th
									 -bx*2.0*sin4th*costh
									 +bx*sinth*sin2th*sin2th
			 )  *cosph;

	ddz_W0W0= (        ddell* sinth
							   -ddbz  *sin2th
			  );

	ddz_W1W1=(
									-ell*sinth
									+bz*4.0*sin2th
			  );



	ddz_W2W2=0.0;
	ddz_W0W1=(
									dell*costh
									-dbz*2.0*cos2th
			 );





	ddz_W0W2=0.0;

	ddz_W1W2=0.0;



	//ковариантный метр. тензор

	g00kov=dx_W0*dx_W0+dy_W0*dy_W0+dz_W0*dz_W0;
	g11kov=dx_W1*dx_W1+dy_W1*dy_W1+dz_W1*dz_W1;
	g22kov=dx_W2*dx_W2+dy_W2*dy_W2+dz_W2*dz_W2;
	g01kov=dx_W0*dx_W1+dy_W0*dy_W1+dz_W0*dz_W1;
	g02kov=dx_W0*dx_W2+dy_W0*dy_W2+dz_W0*dz_W2; //=0
	g12kov=dx_W1*dx_W2+dy_W1*dy_W2+dz_W1*dz_W2; //=0
	g=g00kov*g11kov-g01kov*g01kov;

	//проверим чтобы небыло сингулярностей
	if(g==0.0){raylog<<"g=g00kov*g11kov-g01kov*g01kov=0!!!";return 0;}
	if(g22kov==0.0){raylog<<"g22kov=0!!!";return 0;}
	if(g11kov==0.0){raylog<<"g11kov=0!!!";return 0;}


	//контрвариантный метр. тензор
	g00=g11kov/g;
	g11=g00kov/g;
	g22=1.0/g22kov;
	g01=-g01kov/g;
	g02=0.0;
	g12=0.0;

	dpsi=cor.ret_dpsi(W[0]);

	Bpol=signBpol*10000*dpsi*(1.0/cor.rm)*(1.0/cor.rm)*sqrt(g00/g22kov);  //1)ф-лы в системе си!
	Btor=10000*Btor_0*R0/R;            //2)*10000 т.к нужно перевести в гаусы!


	out<<endl<<"(180/pi)*atan(Bpol/Btor)= "<<(180/pi)*atan(Bpol/Btor)<<endl;

	plotn=ne_spline(W[0]);

	out<<"plotn na granice "<<plotn<<endl;




	Bx=(dx_W1*Bpol/sqrt(g11kov))+(dx_W2*Btor/sqrt(g22kov));
	By=(dy_W1*Bpol/sqrt(g11kov))+(dy_W2*Btor/sqrt(g22kov));
	Bz=(dz_W1*Bpol/sqrt(g11kov))+(dz_W2*Btor/sqrt(g22kov));

	B2=Bx*Bx+By*By+Bz*Bz; //B^2
	B=sqrt(B2);

	//cout<<cor.z0*100<<cor.rm*100<<endl;
	raylog<<setiosflags(ios_base::scientific)<<setprecision(17);
	raylog<<"x= "<<x_real<<" y= "<<y_real<<" z= "<<z_real<<" Bpol= "<<Bpol<<" Btor= "<<Btor <<endl;
	raylog<<"Bx= "<<Bx<<" By= "<<By<<" Bz= " <<Bz<<endl;


///////////проверим найденные Bx By Bz, а затем Выберем моду

	double t0_n,t1_n,t2_n,t0_v,t1_v,t2_v,Bpoltst,Brtst,Btortst;
	double N[3];

	t0_n=dx_W0*Bx+dy_W0*By+dz_W0*Bz;
	t1_n=dx_W1*Bx+dy_W1*By+dz_W1*Bz;
	t2_n=dx_W2*Bx+dy_W2*By+dz_W2*Bz;
	t0_v=t0_n*g00+t1_n*g01+t2_n*g02;
	t1_v=t0_n*g01+t1_n*g11+t2_n*g12;
	t2_v=t0_n*g02+t1_n*g12+t2_n*g22;

	Bpoltst=t1_v*sqrt(g11kov);
	Brtst=t0_v*sqrt(g00kov);
	Btortst=t2_v*sqrt(g22kov);


	raylog<<"test! sdelaem obr preobr koord"<<endl<<"Br=  "<<Brtst<<"  Bpol= "<<Bpoltst<<"  Btor= "<<Btortst<<endl;

	//NX_0=0;
	NY_0=n1y0;
	NZ_0=n1z0;
	//NX_0=1,-1 в вакуме! в плазме должно слегка отличаться её и будем искать!
	//выберем его! Тем самым определим моду волны


	V=zar*zar*4.0*pi*plotn/(m_e*F*F);
	U=zar*zar*B*B/(m_e*m_e*c*c*F*F);

	raylog<<"B_start "<<B<<" plotn_start "<<plotn<<endl;
	raylog<<"U(B) "<<U<<" V(plotn) "<<V<<endl;

	//H=Nper^4*(1-u-v)+Nper^2*(2*U-U*V-2-2*V*V+4*V)+Nper^2*Npar^2*(2-2*u-2*v)+
	//Npar^4*(1-u-v+u*v)+Npar^2*(-2+2U+4V-2*U*V-2*V*V)+(1-U-3V+U*V+3*V*V-V*V*V)

	double K1,K2,K3,K4,K5,K6;
	double M1,M2,M3,M4,M5;
	double C0,C1,C2,C3,C4;

	K1=1.0-U-V;   //Nper^4
	K2=2*U-U*V-2.0-2.0*V*V+4.0*V;  //Nper^2
	K3=2.0-2.0*U-2.0*V+U*V; //Nper^2*Npar^2
	K4=1.0-U-V+U*V;  //Npar^4
	K5 =  -2.0+2.0*U+4.0*V-2.0*U*V-2.0*V*V;  //Npar^2
	K6=1.0-U-3.0*V+U*V+3.0*V*V-V*V*V;  // ^0

	M1=(NY_0*By+NZ_0*Bz)/B; //alpha

	M2=Bx/B;                 //betta
	M3=NY_0*NY_0+NZ_0*NZ_0;  //gamma^2

	M4=M1*M1;                //alpha^2
	M5=M2*M2;                //betta^2
	//итак ур-е которое нужно реш. содерж след коэф.

	C0=K1*(M3-M4)*(M3-M4)+K2*(M3-M4)+K3*(M3-M4)*M4+K4*(M4*M4)+K5*(M4)+K6;
	C1=K1*(-4*M1*M2*(M3-M4))+K2*(-2*M1*M2)+K3*(-2*M1*M1*M1*M2+2*M1*M2*(M3-M4))+K4*(4*M1*M1*M1*M2)+K5*(2*M1*M2);
	C2=K1*(4*M5*M4+2*(M3-M4)*(1-M5))+K2*(1-M5)+K3*(M4*(1-M5)+M5*(M3-M4)-4*M4*M5)+K4*(6*M4*M5)+K5*(M5);
	C3=K1*(-4*M1*M2*(1-M5))+K2*(0)+K3*(-2*M1*M2*M2*M2+2*M1*M2*(1-M5))+K4*(4*M1*M2*M2*M2)+K5*(0);
	C4=K1*((1-M5)*(1-M5))+K2*(0)+K3*(M5*(1-M5))+K4*(M5*M5)+K5*(0);


   //2) найдем корни ур-я с заданными коэфф-ми
	raylog << fixed << setprecision(17);

     int i;
	 const int M=4, MP1=M+1;
     const complex<DP> real1(1.0,0.0),imag1(0.0,1.0);
     const complex<DP> a_d[MP1]=
     {real1*C0,C1,C2,C3,C4};    //коээфиценты полинома 5 степени ax^5+bx^4+cx^3+dx^2+ex^1+fx^0
     bool polish;
     Vec_CPLX_DP a(a_d,MP1),roots(M);

     //  raylog << endl << "Roots of the polynomial x^4-1" << endl;


	out << fixed << setprecision(17);



    polish=false;
    NR::zroots(a,roots,polish);
    for (i=0;i<M;i++)roots[i]=(DP(1.0)+DP(0.01)*(i+1))*roots[i];
    polish=true;
    NR::zroots(a,roots,polish);
    out << endl << "Polished roots:" << endl;
    out << setw(14) << "root #" << setw(14) << "root:" << endl << endl;
    for (i=0;i<M;i++) {
    out << setw(11) << noshowpos << i;
    out << setw(25) << showpos << roots[i] << endl;
    }


		double H;
		double Npar,Npar2,Nper,Nper2;

		//задаем выбранную моду
		i=moda;
		//roots[1]=roots[3]/2+roots[2]/2;
        //roots[1]=1;
		Npar=(roots[i].real()*Bx+NY_0*By+NZ_0*Bz)/B;
		Nper2=(roots[i].real()*roots[i].real()+NY_0*NY_0+NZ_0*NZ_0-Npar*Npar);
        H=Nper2*Nper2*K1+Nper2*K2+Nper2*Npar*Npar*K3+Npar*Npar*Npar*Npar*K4+Npar*Npar*K5+K6;

		out<<"dlya tekyshego kornya na starte"<<endl;
		out<<scientific<<" Npar= "<<Npar<<" Nper= "<<sqrt(Nper2)<<" H= "<<H<<endl;

        NX_0=roots[i].real();

		out<<endl<<"Nx_start= "<<NX_0<<" Ny_start= "<<NY_0<<" Nz_start= "<<NZ_0;
        //exit(1);


//теперь Посчитаем коэфф преломления в новых коорд.

		t0_n=dx_W0*NX_0+dy_W0*NY_0+dz_W0*NZ_0; // нач норм волновой вектор в нижних коорд
		t1_n=dx_W1*NX_0+dy_W1*NY_0+dz_W1*NZ_0;
		t2_n=dx_W2*NX_0+dy_W2*NY_0+dz_W2*NZ_0;
		t0_v=t0_n*g00+t1_n*g01+t2_n*g02;         // нач норм волновой вектор в верхн коорд
		t1_v=t0_n*g01+t1_n*g11+t2_n*g12;
		t2_v=t0_n*g02+t1_n*g12+t2_n*g22;

		N[0]=t0_n;
		N[1]=t1_n;
		N[2]=t2_n;

		//И магнитное поле

		double B0_v,B1_v,B2_v;

		B0_v=0.0;
		B1_v=Bpol/sqrt(g11kov);
		B2_v=Btor/sqrt(g22kov);

		Npar=((N[1])*B1_v+N[2]*B2_v)/B;
		Npar2=Npar*Npar;
		Nper2=N[0]*N[0]*g00+2*N[0]*N[1]*g01+N[1]*N[1]*g11+N[2]*N[2]*g22-Npar2;
		Nper=sqrt(Nper2);
		H=Nper2*Nper2*K1+Nper2*K2+Nper2*Npar2*K3+Npar2*Npar2*K4+Npar2*K5+K6;
		raylog<<H;


///////////////////////////////////////////////////////
//Итак приступим к решению уравнений
		double s;

		Vec_DP P(6),dpdx(6),dpsav(6),Psav(6),pscal(6);
		DP eps,hdid,hnext,htry;


		P[0]=W[0];P[1]=W[1];P[2]=W[2];P[3]=N[0];P[4]=N[1];P[5]=N[2];


		s=0;
		eps=1e-15; //задаем точность
		htry=0.1;

		for (i=0;i<6;i++) pscal[i]=1.0;

		raylog<<endl<<"w"<<P[0];

		ray<<setprecision(10);
		ray<<left<<setw(14)<<fixed<<"rho"<<left<<setw(14)<<"X"<<left<<setw(14)<<"Y"<<left<<setw(14)<<"Z"<<left<<setw(14)<<"R"<<left<<setw(14)<<"Npar"<<left<<setw(14)<<"Nper"<<left<<setw(14)<<"Bpol/Btor"<<left<<setw(14)<<"Step"<<endl;

		derivs (s,P,dpdx); //Высчитываем для теста, и для вывода все правые части

		if(!ret) return 0;  // Всё ок?

		//exit(1);

		ray<<"-------------------------------------------------------------------------------------------"<<endl;

		out3<<endl<<setprecision(11)<<left<<setw(16)<<fixed<<"rho"<<left<<setw(16)<<"theta"<<left<<setw(16)<<"phi";
		out3<<endl<<setprecision(11)<<left<<setw(16)<<fixed<<P[0]<<left<<setw(16)<<P[1]<<left<<setw(16)<<P[2];

		i=0;

		double w_last,w_rem,th_rem,phi_rem,n0rem,n1rem,n2rem; //параметры в отсечке


		//flag3=1;
		w_last=W[0];

		derivs (s,P,dpdx);

		ne_log<<W[0]<<" "<<ne_spline(W[0])<<endl;

		if(!ret) return 0;

		while(P[0]<W[0]+0.00001){  //решаем пока не вернулись на то же rho
			i++;
			NR::rkqs(P,dpdx,s,htry,eps,pscal,hdid,hnext,derivs);
			//cout << "iteration: " << i << endl;
			//cout << P[0] << "<" << W[0] + 0.00001 << endl;
			flag=1;
			flag2=1;
			ne_log<<endl<<P[0]<<" "<<ne_spline(P[0]);
			//raylog<<endl<<"aa"<<endl;
			htry=hnext;
			//htry=0.1;
			raylog<<endl<<" rho = "<<P[0];
			derivs (s,P,dpdx);
			if(!ret) return 0;
			ray<<" "<<left<<setw(14)<<htry<<endl;

			if(flag3) if(P[0]-w_last>0){w_rem=P[0];th_rem=P[1],phi_rem=P[2],n0rem=P[3];n1rem=P[4];n2rem=P[5]; flag3=0;}

			w_last=P[0];

			out3<<endl<<setprecision(11)<<left<<setw(16)<<fixed<<P[0]<<left<<setw(16)<<P[1]<<left<<setw(16)<<P[2];

		}

		//смотрим как долго решали
		endtime=time(NULL);


		raylog<<endl<<endl<<"RunTime (full sec) = "<<difftime(endtime,starttime)<<endl;

		//const time_t timer = time(NULL);
		//raylog<<"DATA = " << std::ctime(&timer);

		raylog<<"rho min = "<<w_rem<<"plotn = " << ne_spline(W[0])<<endl;

		raylog<<endl<<endl<<"All ok!! Please press enter";


///****


	out.close();
	cutoff.close();
	out3.close();
	ray.close();
	initial.close();
	raylog.close();
	btor.close();
	ne_log.close();

	std::cout << "still alive" << std::endl;

	return 1;

}//END main


//обобщ метод использ алгоритма лаггера
//void NR::zroots(Vec_I_CPLX_DP &a, Vec_O_CPLX_DP &roots, const bool &polish)

//алгоритм лаггера для нахожд корней полиномов
//void NR::laguer(Vec_I_CPLX_DP &a, complex<DP> &x, int &its)

//void NR::rkqs(Vec_IO_DP &y, Vec_IO_DP &dydx, DP &x, const DP htry,
//	const DP eps, Vec_I_DP &yscal, DP &hdid, DP &hnext,
//	void derivs(const DP, Vec_I_DP &, Vec_O_DP &))

//void NR::rkck(Vec_I_DP &y, Vec_I_DP &dydx, const DP x,
//	const DP h, Vec_O_DP &yout, Vec_O_DP &yerr,
//	void derivs(const DP, Vec_I_DP &, Vec_O_DP &))

void derivs (const DP s, Vec_I_DP &P, Vec_O_DP &dydx) //Вычисляет правые части для рунге кутта
{
	//***************************************************
	//значение переменной по которой интегрируется
	//P вектор точки в который вычисляем правую часть
	//dydx значение правой части
	//***************************************************
	double W[3],N[3];

	 W[0]=P[0];
	 W[1]=P[1];
	 W[2]=P[2];
	 N[0]=P[3];
	 N[1]=P[4];
	 N[2]=P[5];



	double shf,ell,trin,bx,bz,dshf,dell,dtrin,dbx,dbz,ddshf,ddell,ddtrin,ddbx,ddbz;
    double x,y,z,x_real,y_real,z_real,R;
    double sinth,costh,sin2th,cos2th,sin4th,cos4th,sinph,cosph;
    double dx_W0,dx_W1,dx_W2,dy_W0,dy_W1,dy_W2,dz_W0,dz_W1,dz_W2,dR_W0,dR_W1;

    double ddx_W0W0,ddx_W1W1,ddx_W2W2,ddx_W0W1,ddx_W0W2,ddx_W1W2;
    double ddy_W0W0,ddy_W1W1,ddy_W2W2,ddy_W0W1,ddy_W0W2,ddy_W1W2;
    double ddz_W0W0,ddz_W1W1,ddz_W2W2,ddz_W0W1,ddz_W0W2,ddz_W1W2;



    double g00kov,g11kov,g22kov,g01kov,g02kov,g12kov,g;
    double dg00kov_W0,dg00kov_W1,dg00kov_W2,dg11kov_W0,dg11kov_W1,dg11kov_W2;
    double dg22kov_W0,dg22kov_W1,dg22kov_W2,dg01kov_W0,dg01kov_W1,dg01kov_W2;
    double dg02kov_W0,dg02kov_W1,dg02kov_W2,dg12kov_W0,dg12kov_W1,dg12kov_W2;
    double dg_W0,dg_W1,dg_W2,d1delg_W0,d1delg_W1,d1delg_W2;
    double dsq1delg11kov_W0,dsq1delg11kov_W1,dsq1delg22kov_W0,dsq1delg22kov_W1;
    double dsq1delg11kov_W2,dsq1delg22kov_W2;

    double g00,g11,g22,g01,g02,g12;
    double dg00_W0,dg00_W1,dg00_W2,dg11_W0,dg11_W1,dg11_W2,dg22_W0,dg22_W1,dg22_W2;
    double dg01_W0,dg01_W1,dg01_W2;

    double dpsi,ddpsi,Bpol,B2,B,Btor
	,R0=R_center/cor.rm,plotn,Bx,By,Bz; // R0=Rcenter=0.36 пока не знаю почему!
	//незабыть читать Btor_0 из файла!!
    double U,V;

    const double k0=F/c;

	double dN_W0,dB_W0,dB_W1,dB2_W0,dB2_W1,dBpol_W0,dBpol_W1,dBtor_W0,dBtor_W1;
    double dU_W0,dU_W1,dU_W2,dV_W0,dV_W1,dV_W2;
    double dK1_W0,dK1_W1,dK2_W0,dK2_W1,dK3_W0,dK3_W1,dK4_W0,dK4_W1,dK5_W0,dK5_W1,dK6_W0,dK6_W1;
	double K1,K2,K3,K4,K5,K6;

	double B0_v,B1_v,B2_v;
    double dH_Npar2,dH_Nper2;
    double dNpar_N0,dNpar_N1,dNpar_N2,dNpar_W0,dNpar_W1,dNpar_W2;
    double dNpar2_N0,dNpar2_N1,dNpar2_N2,dNpar2_W0,dNpar2_W1,dNpar2_W2;
    double dNper2_N0,dNper2_N1,dNper2_N2,dNper2_W0,dNper2_W1,dNper2_W2;

	double H;
	double Npar,Npar2,Nper,Nper2;


	double R_ekv; //эквивалентное R на магнитной оси, почти совпадает с экватариальной плоскостью. (для той же магнитной поверхности)


    shf=cor.ret_shf(W[0]);
	ell=cor.ret_ell(W[0]);
	trin=cor.ret_trin(W[0]);
	bx=cor.ret_bx(W[0]);
	bz=cor.ret_bz(W[0]);

	dshf=cor.ret_dshf(W[0]);
	dell=cor.ret_dell(W[0]);
	dtrin=cor.ret_dtrin(W[0]);
	dbx=cor.ret_dbx(W[0]);
	dbz=cor.ret_dbz(W[0]);

	ddshf=cor.ret_ddshf(W[0]);
	ddell=cor.ret_ddell(W[0]);
	ddtrin=cor.ret_ddtrin(W[0]);
	ddbx=cor.ret_ddbx(W[0]);
	ddbz=cor.ret_ddbz(W[0]);

	sinth=sin(W[1]);  //Чтоб постоянно не считать посчитаем все синусы-косинусы
	sin2th=sin(2.0*W[1]);
	sin4th=sin(4.0*W[1]);
	costh=cos(W[1]);
	cos2th=cos(2.0*W[1]);
	cos4th=cos(4.0*W[1]);
	sinph=sin(W[2]);
	cosph=cos(W[2]);

	 R= (cor.x0+            W[0]*costh  -shf
							     -trin*sinth*sinth
							     -bx*sin2th*sin2th*costh);



	 R_ekv=(cor.x0+            W[0]*1  -shf   );


	x=R*cosph;
	y=R*sinph;

	z=cor.z0+            ell*sinth -bz  *sin2th    ;


	//x,y,z физические

	x_real= cor.rm*x;
    y_real= cor.rm*y;
	z_real= cor.rm*z;

	dR_W0= (costh  -dshf
							     -dtrin*  sinth*sinth
							     -dbx   *sin2th*sin2th*costh)  ;
	dR_W1= (
			                      -W[0]*sinth
							     -trin* sin2th
							     -bx*2.0*sin4th*costh
								 +bx*sinth*sin2th*sin2th );






	//первые производные по координатам
	dx_W0 = dR_W0 *cosph;
	dx_W1= dR_W1 *cosph;
	dx_W2= -R * sinph;

	dy_W0=  dR_W0 *sinph;
	dy_W1= dR_W1*sinph;
	dy_W2= R*cosph;

	dz_W0= (      dell* sinth
				  -dbz  *sin2th
		   );

	dz_W1=(
									ell*costh
									-bz*2.0*cos2th
		   );



	dz_W2=0.0;

	//______________________________________________
	   //вторые производные по координатам


	ddx_W0W0= (    -ddshf
									 -ddtrin*  sinth*sinth
									 -ddbx   *sin2th*sin2th*costh
			  ) *cosph;

	ddx_W1W1=(                  -W[0]*costh
									 -2.0*trin* cos2th
									 -bx*(
														  +8.0*cos4th*costh
														  -4.0*sin4th*sinth
														  -sin2th*sin2th*costh
										 )
			)  *cosph;
	ddx_W2W2=-(cor.x0+            W[0]*costh  -shf
									 -trin*  sinth*sinth
									 -bx   *sin2th*sin2th*costh

			  )*cosph;
	ddx_W0W1=(                       -sinth
									 -dtrin* sin2th
									 -dbx*2.0*sin4th*costh
									 +dbx*sinth*sin2th*sin2th
			 )  *cosph;

	ddx_W0W2=  -(    costh  -dshf
									 -dtrin*  sinth*sinth
									 -dbx   *sin2th*sin2th*costh
				) *sinph;


	ddx_W1W2=-(
									  -W[0]*sinth
									 -trin* sin2th
									 -bx*2.0*sin4th*costh
									 +bx*sinth*sin2th*sin2th
			  )  *sinph;

	ddy_W0W0= (    -ddshf
									 -ddtrin*  sinth*sinth
									 -ddbx   *sin2th*sin2th*costh
			  ) *sinph;


	ddy_W1W1=(
									  -W[0]*costh
									 -2.0*trin* cos2th
									 -bx*(
														  +8.0*cos4th*costh
														  -4.0*sin4th*sinth
														  -sin2th*sin2th*costh


													)
			)  *sinph;

	ddy_W2W2=-(cor.x0+            W[0]*costh  -shf
									 -trin*  sinth*sinth
									 -bx   *sin2th*sin2th*costh

			 )*sinph;
	ddy_W0W1=(
									  -sinth
									 -dtrin* sin2th
									 -dbx*2.0*sin4th*costh
									 +dbx*sinth*sin2th*sin2th
			 )  *sinph;

	ddy_W0W2= (    costh  -dshf
									 -dtrin*  sinth*sinth
									 -dbx   *sin2th*sin2th*costh
			  ) *cosph;


	ddy_W1W2=(
									  -W[0]*sinth
									 -trin* sin2th
									 -bx*2.0*sin4th*costh
									 +bx*sinth*sin2th*sin2th
			 )  *cosph;

	ddz_W0W0= (        ddell* sinth
							   -ddbz  *sin2th
			  );

	ddz_W1W1=(
									-ell*sinth
									+bz*4.0*sin2th
			  );



	ddz_W2W2=0.0;
	ddz_W0W1=(
									dell*costh
									-dbz*2.0*cos2th
			 );





	ddz_W0W2=0.0;

	ddz_W1W2=0.0;



	//ковариантный метр. тензор

	g00kov=dx_W0*dx_W0+dy_W0*dy_W0+dz_W0*dz_W0;
	g11kov=dx_W1*dx_W1+dy_W1*dy_W1+dz_W1*dz_W1;
	g22kov=dx_W2*dx_W2+dy_W2*dy_W2+dz_W2*dz_W2;
	g01kov=dx_W0*dx_W1+dy_W0*dy_W1+dz_W0*dz_W1;
	//g02kov=dx_W0*dx_W2+dy_W0*dy_W2+dz_W0*dz_W2; //=0
	//g12kov=dx_W1*dx_W2+dy_W1*dy_W2+dz_W1*dz_W2; //=0
	g02kov=0;
	g12kov=0;
	g=g00kov*g11kov-g01kov*g01kov;

	//возможные ошибки и сингулярности?
	if(g==0.0){raylog<<"g=g00kov*g11kov-g01kov*g01kov=0!!!";ret=0;}
	if(g22kov==0.0){raylog<<"g22kov=0!!!";ret=0;}
	if(g11kov==0.0){raylog<<"g11kov=0!!!";ret=0;}

	//производные ковариантного метр. тензора

	dg00kov_W0=2.0*dx_W0*ddx_W0W0+2.0*dy_W0*ddy_W0W0+2.0*dz_W0*ddz_W0W0;
	dg00kov_W1=2.0*dx_W0*ddx_W0W1+2.0*dy_W0*ddy_W0W1+2.0*dz_W0*ddz_W0W1;
	//dg00kov_W2=2.0*dx_W0*ddx_W0W2+2.0*dy_W0*ddy_W0W2+2.0*dz_W0*ddz_W0W2; =0
	dg11kov_W0=2.0*dx_W1*ddx_W0W1+2.0*dy_W1*ddy_W0W1+2.0*dz_W1*ddz_W0W1;
	dg11kov_W1=2.0*dx_W1*ddx_W1W1+2.0*dy_W1*ddy_W1W1+2.0*dz_W1*ddz_W1W1;
	//dg11kov_W2=2.0*dx_W1*ddx_W1W2+2.0*dy_W1*ddy_W1W2+2.0*dz_W1*ddz_W1W2; =0
	dg22kov_W0=2.0*dx_W2*ddx_W0W2+2.0*dy_W2*ddy_W0W2+2.0*dz_W2*ddz_W0W2;
	dg22kov_W1=2.0*dx_W2*ddx_W1W2+2.0*dy_W2*ddy_W1W2+2.0*dz_W2*ddz_W1W2;
	//dg22kov_W2=2.0*dx_W2*ddx_W2W2+2.0*dy_W2*ddy_W2W2+2.0*dz_W2*ddz_W2W2; =0
	dg01kov_W0=dx_W0*ddx_W0W1+dx_W1*ddx_W0W0+dy_W0*ddy_W0W1+dy_W1*ddy_W0W0+dz_W0*ddz_W0W1+dz_W1*ddz_W0W0;
	dg01kov_W1=dx_W0*ddx_W1W1+dx_W1*ddx_W0W1+dy_W0*ddy_W1W1+dy_W1*ddy_W0W1+dz_W0*ddz_W1W1+dz_W1*ddz_W0W1;
	//dg01kov_W2=dx_W0*ddx_W1W2+dx_W1*ddx_W0W2+dy_W0*ddy_W1W2+dy_W1*ddy_W0W2+dz_W0*ddz_W1W2+dz_W1*ddz_W0W2; =0
	dg02kov_W0=dx_W0*ddx_W0W2+dx_W2*ddx_W0W0+dy_W0*ddy_W0W2+dy_W2*ddy_W0W0+dz_W0*ddz_W0W2+dz_W2*ddz_W0W0;
	dg02kov_W1=dx_W0*ddx_W1W2+dx_W2*ddx_W0W1+dy_W0*ddy_W1W2+dy_W2*ddy_W0W1+dz_W0*ddz_W1W2+dz_W2*ddz_W0W1;
	//dg02kov_W2=dx_W0*ddx_W2W2+dx_W2*ddx_W0W2+dy_W0*ddy_W2W2+dy_W2*ddy_W0W2+dz_W0*ddz_W2W2+dz_W2*ddz_W0W2; =0
	dg12kov_W0=dx_W1*ddx_W0W2+dx_W2*ddx_W0W1+dy_W1*ddy_W0W2+dy_W2*ddy_W0W1+dz_W1*ddz_W0W2+dz_W2*ddz_W0W1;
	dg12kov_W1=dx_W1*ddx_W1W2+dx_W2*ddx_W1W1+dy_W1*ddy_W1W2+dy_W2*ddy_W1W1+dz_W1*ddz_W1W2+dz_W2*ddz_W1W1;
	//dg12kov_W2=dx_W1*ddx_W2W2+dx_W2*ddx_W1W2+dy_W1*ddy_W2W2+dy_W2*ddy_W1W2+dz_W1*ddz_W2W2+dz_W2*ddz_W1W2; =0
	dg00kov_W2=0;
	dg11kov_W2=0;
	dg22kov_W2=0;
	dg01kov_W2=0;
	dg02kov_W2=0;
	dg12kov_W2=0;



	dg_W0=dg00kov_W0*g11kov+g00kov*dg11kov_W0-2.0*g01kov*dg01kov_W0;  //g=g00kov*g11kov-g01kov*g01kov
	dg_W1=dg00kov_W1*g11kov+g00kov*dg11kov_W1-2.0*g01kov*dg01kov_W1;
	//dg_W2=dg00kov_W2*g11kov+g00kov*dg11kov_W2-2.0*g01kov*dg01kov_W2;=0
	dg_W2=0;


	d1delg_W0=-dg_W0/(g*g);  //производные 1/g
	d1delg_W1=-dg_W1/(g*g);
	//d1delg_W2=-dg_W2/(g*g);=0
	d1delg_W2=0;

	dsq1delg11kov_W0=-0.5*(1.0/g11kov)*sqrt(1.0/g11kov)*dg11kov_W0;          //производные sqrt(1/gijcov)
	dsq1delg22kov_W0=-0.5*(1.0/g22kov)*sqrt(1.0/g22kov)*dg22kov_W0;
	dsq1delg11kov_W1=-0.5*(1.0/g11kov)*sqrt(1.0/g11kov)*dg11kov_W1;
	dsq1delg22kov_W1=-0.5*(1.0/g22kov)*sqrt(1.0/g22kov)*dg22kov_W1;
	//dsq1delg11kov_W2=-0.5*(1.0/g11kov)*sqrt(1.0/g11kov)*dg11kov_W2; =0
	//dsq1delg22kov_W2=-0.5*(1.0/g22kov)*sqrt(1.0/g22kov)*dg22kov_W2; =0

	dsq1delg11kov_W2=0;
	dsq1delg22kov_W2=0;


	//контрвариантный метр. тензор
	g00=g11kov/g;
	g11=g00kov/g;
	g22=1.0/g22kov;
	g01=-g01kov/g;
	g02=0.0;
	g12=0.0;

	//raylog<<endl<<"g02=0! "<<g12<<endl;
	//raylog<<endl;

	//производные контрвариантного метр. тензора

	dg00_W0=(dg11kov_W0/g)+g11kov*d1delg_W0;
	dg00_W1=(dg11kov_W1/g)+g11kov*d1delg_W1;
	//dg00_W2=(dg11kov_W2/g)+g11kov*d1delg_W2;=0
	dg00_W2=0;

	dg11_W0=(dg00kov_W0/g)+g00kov*d1delg_W0;
	dg11_W1=(dg00kov_W1/g)+g00kov*d1delg_W1;
	//dg11_W2=(dg00kov_W2/g)+g00kov*d1delg_W2;=0
	dg11_W2=0;

	dg22_W0=-dg22kov_W0/(g22kov*g22kov);
	dg22_W1=-dg22kov_W1/(g22kov*g22kov);
	//dg22_W2=-dg22kov_W2/(g22kov*g22kov);=0
	dg22_W2=0;

	dg01_W0=-(dg01kov_W0/g)-g01kov*d1delg_W0;
	dg01_W1=-(dg01kov_W1/g)-g01kov*d1delg_W1;
	//dg01_W2=-(dg01kov_W2/g)-g01kov*d1delg_W2;=0
	dg01_W2=0;

	dpsi=cor.ret_dpsi(W[0]);
	ddpsi=cor.ret_ddpsi(W[0]);

	Bpol=signBpol*10000*dpsi*(1/cor.rm)*(1/cor.rm)*sqrt(g00/g22kov);
	Btor=10000*Btor_0*R0/R;            //10000 т.к нужно перевести теслы в гаусы!

	//Bpol=pow(10.0f,-30.0f);
	//Btor=pow(10.0f,-30.0f);
	//plotn=1.0*pow(10.0f,+9.0f);
	plotn=ne_spline(W[0]);
	//plotn=1.0*pow(10.0,11.0);

	//plotn=0;

	Bx=(dx_W1*Bpol/sqrt(g11kov))+(dx_W2*Btor/sqrt(g22kov));
	By=(dy_W1*Bpol/sqrt(g11kov))+(dy_W2*Btor/sqrt(g22kov));
	Bz=(dz_W1*Bpol/sqrt(g11kov))+(dz_W2*Btor/sqrt(g22kov));

	//raylog<<endl<<endl<<Bz;exit(1);

	B2=Bx*Bx+By*By+Bz*Bz; //B^2
	B=sqrt(B2);

	//raylog<<endl<<B;exit(1);

	V=zar*zar*4.0*pi*plotn/(m_e*F*F);
	U=zar*zar*B2/(m_e*m_e*c*c*F*F);

	//dN_W0=   -plotn0*(1/W_s)*(1/cosh((W[0]-W_delt)/W_s)) *(1/cosh((W[0]-W_delt)/W_s))  ;
	dN_W0 = ne_spline.deriv(1, W[0]);
	//dN_W0=0;

	dBtor_W0=-dR_W0*10000*Btor_0*R0/(R*R);
	dBtor_W1=-dR_W1*10000*Btor_0*R0/(R*R);
	dBpol_W0=signBpol*  (10000*ddpsi*(1.0/cor.rm)*(1.0/cor.rm)*sqrt(g00/g22kov)+
			 10000*dpsi*(1.0/cor.rm)*(1.0/cor.rm)*(0.5)*sqrt(g22kov/g00)*(dg00_W0*g22kov-dg22kov_W0*g00)/(g22kov*g22kov)   );
	dBpol_W1=
			 signBpol*   (10000*dpsi*(1.0/cor.rm)*(1.0/cor.rm)*(0.5)*sqrt(g22kov/g00)*(dg00_W1*g22kov-dg22kov_W1*g00)/(g22kov*g22kov)  );

	dB2_W0=2.0*Bpol*dBpol_W0+2.0*Btor*dBtor_W0;
	dB2_W1=2.0*Bpol*dBpol_W1+2.0*Btor*dBtor_W1;
	dB_W0=0.5*(1.0/B)*dB2_W0;
	dB_W1=0.5*(1.0/B)*dB2_W1;

	dV_W0= zar*zar*4.0*pi*dN_W0/(m_e*F*F);
	dV_W1=0;dV_W2=0;
	dU_W0=zar*zar*dB2_W0/(m_e*m_e*c*c*F*F);
	dU_W1=zar*zar*dB2_W1/(m_e*m_e*c*c*F*F);
	dU_W2=0;

	K1=1.0-U-V;   //Nper^4
	K2=2*U-U*V-2.0-2.0*V*V+4.0*V;  //Nper^2
	K3=2.0-2.0f*U-2.0*V+U*V;  //Nper^2*Npar^2
	K4=1.0-U-V+U*V;  //Npar^4
	K5 =  -2.0+2.0f*U+4.0*V-2.0*U*V-2.0*V*V;  //Npar^2
	K6=1.0-U-3.0*V+U*V+3.0*V*V-V*V*V;  // ^0


	dK1_W0=-dU_W0-dV_W0;
	dK1_W1=-dU_W1;
	dK2_W0=2*dU_W0-dU_W0*V-U*dV_W0-4*V*dV_W0+4*dV_W0;
	dK2_W1=2*dU_W1-dU_W1*V;
	dK3_W0=-2.0*dU_W0-2.0*dV_W0+dU_W0*V+U*dV_W0;
	dK3_W1=-2.0*dU_W1+dU_W1*V;
	dK4_W0=-dU_W0-dV_W0+dU_W0*V+U*dV_W0;
	dK4_W1=-dU_W1+dU_W1*V;
	dK5_W0 =  2.0*dU_W0+4.0*dV_W0-2.0*dU_W0*V-2.0*U*dV_W0-4.0*V*dV_W0;
	dK5_W1 =  2.0*dU_W1-2.0*dU_W1*V;
	dK6_W0=-dU_W0-3.0*dV_W0+dU_W0*V+U*dV_W0+6.0*V*dV_W0-3*V*V*dV_W0;
	dK6_W1=-dU_W1+dU_W1*V;


	//raylog<<endl<<dK1_W1;
	//exit(1);

	B0_v=0.0;
	B1_v=Bpol/sqrt(g11kov);
	B2_v=Btor/sqrt(g22kov);
	//N[0]=t0_n;
	//N[1]=t1_n;
	//N[2]=t2_n;

	Npar=(N[1]*B1_v+N[2]*B2_v)/B;
	Npar2=Npar*Npar;
	Nper2=N[0]*N[0]*g00+2*N[0]*N[1]*g01+N[1]*N[1]*g11+N[2]*N[2]*g22-Npar2;
	Nper=sqrt(Nper2);
	H=Nper2*Nper2*K1+Nper2*K2+Nper2*Npar2*K3+Npar2*Npar2*K4+Npar2*K5+K6;
	if(flag2==1){ raylog<<scientific<<" H= "<<H; flag2=0;}


	if(flag33)if(!flag3)
	{

	out<<endl<<endl<<"rho min = "<<W[0]<<"  plotn*10e12(1/cm^3) = "<<(ne_spline(W[0]))/1e12;flag33=0;
	out<<endl<<"X(m) = "<<x_real<<" Y(m) = "<<y_real<<" Z(m) = "<<z_real<<" R(m) = "<<cor.rm*R;
	out<<endl<<"Npar = "<<Npar<<" Nper = "<<Nper<<" abs((Npar/Nper)*100)= "<<abs((Npar/Nper)*100);
	out<<endl<<"Kpar (1/sm) = "<<Npar*k0<<" Kper(1/sm) = "<<Nper*k0<<" угол относительно Kреш (в градусах) = "<<(189/pi)*atan(Npar/Nper);
	cutoff<<W[0]<<' '<<(ne_spline(W[0]))<<' '<<cor.rm*R<<' '<<z_real<<' '<<cor.rm*R_ekv<<' '<<cor.z0*cor.rm<<' '<<Npar*k0<<' '<<Nper*k0;
	cutoff<<endl<<"rho     |n_e (1/cm^3)  | R(m)  |  Z(m)   | R_ekv(m)|  Z_ekv(m) |Kpar (cm) |Kper(cm)";
	}


	dH_Npar2=Nper2*K3+2.0*Npar2*K4+K5;

	dH_Nper2=2*Nper2*K1+Npar2*K3+K2;
	dNpar2_N0=0;
	dNper2_N0=2*N[0]*g00+2*N[1]*g01;

	dNpar_N0=0;
	dNpar_N1=(Bpol/B)/sqrt(g11kov);
	dNpar_N2=(Btor/B)/sqrt(g22kov);
	dNpar_W0=     ( N[1]*dBpol_W0/sqrt(g11kov) + N[1]*Bpol*dsq1delg11kov_W0+
					N[2]*dBtor_W0/sqrt(g22kov) + N[2]*Btor*dsq1delg22kov_W0  )*(1.0/B)
				  -( N[1]*Bpol/sqrt(g11kov) + N[2]*Btor/sqrt(g22kov))*(1.0/B)*(1.0/B)*(dB_W0);

	dNpar_W1=     ( N[1]*dBpol_W1/sqrt(g11kov) + N[1]*Bpol*dsq1delg11kov_W1+
					N[2]*dBtor_W1/sqrt(g22kov) + N[2]*Btor*dsq1delg22kov_W1  )*(1.0/B)
				  -( N[1]*Bpol/sqrt(g11kov) + N[2]*Btor/sqrt(g22kov))*(1.0/B)*(1.0/B)*(dB_W1);

	//dNpar_W2=      (N[1]*Bpol*dsq1delg11kov_W2  + N[2]*Bpol*dsq1delg22kov_W2)*(1.0f/B);=0
	  dNpar_W2=0;



	dNpar2_N1=2*Npar*dNpar_N1;
	dNpar2_N2=2*Npar*dNpar_N2;

	dNpar2_W0=2*Npar*dNpar_W0;
	dNpar2_W1=2*Npar*dNpar_W1;
	dNpar2_W2=2*Npar*dNpar_W2;


	dNper2_N1=2*N[0]*g01+2*N[1]*g11-dNpar2_N1;
	dNper2_N2=2*N[2]*g22-dNpar2_N2;
	dNper2_W0=N[0]*N[0]*dg00_W0+2*N[0]*N[1]*dg01_W0+N[1]*N[1]*dg11_W0+N[2]*N[2]*dg22_W0-dNpar2_W0;
	dNper2_W1=N[0]*N[0]*dg00_W1+2*N[0]*N[1]*dg01_W1+N[1]*N[1]*dg11_W1+N[2]*N[2]*dg22_W1-dNpar2_W1;

	//dNper2_W2=N[0]*N[0]*dg00_W2+2*N[0]*N[1]*dg01_W2+N[1]*N[1]*dg11_W2+N[2]*N[2]*dg22_W2-dNpar2_W2;=0
	dNper2_W2=0;

	dydx[0]=dH_Npar2*dNpar2_N0+dH_Nper2*dNper2_N0;
	dydx[1]=dH_Npar2*dNpar2_N1+dH_Nper2*dNper2_N1;
	dydx[2]=dH_Npar2*dNpar2_N2+dH_Nper2*dNper2_N2;

	dydx[3]=dH_Npar2*dNpar2_W0+dH_Nper2*dNper2_W0
     +Nper2*Nper2*dK1_W0+Nper2*dK2_W0+Nper2*Npar2*dK3_W0+Npar2*Npar2*dK4_W0+Npar2*dK5_W0+dK6_W0;

	//m[i][3]=-m[i][3];  // "-" поскольку перед правой частью должен стоять минус для второй тройки ур-й

	dydx[4]=dH_Npar2*dNpar2_W1+dH_Nper2*dNper2_W1
     +Nper2*Nper2*dK1_W1+Nper2*dK2_W1+Nper2*Npar2*dK3_W1+Npar2*Npar2*dK4_W1+Npar2*dK5_W1+dK6_W1;

	//m[i][4]=-m[i][4];

	if(inversT==0)
	{

	dydx[3]=-dydx[3];
	dydx[4]=-dydx[4];
	//dydx[0]=-dydx[0];
	//dydx[1]=-dydx[1];
	//dydx[2]=-dydx[2];
	}

	else
	{
	dydx[3]=dydx[3];
	dydx[4]=dydx[4];
	dydx[0]=-dydx[0];
	dydx[1]=-dydx[1];
	dydx[2]=-dydx[2];
	}


	dydx[5]=0;


	if(flag==1)
	{
		ray<<setprecision(10);
		ray<<left<<setw(14)<<fixed <<W[0]<<setw(8)<<x_real<<' '<<y_real<<' '<<z_real<<' '<<cor.rm*R<<' '<<Npar<<' '<<Nper<<' '<<Bpol/Btor;
		flag=0;
	}

//return 0;

 }





