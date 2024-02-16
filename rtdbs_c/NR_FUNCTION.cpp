
#include "nr.h"
#include <limits>

//обобщ метод использ алгоритма лаггера
void NR::zroots(Vec_I_CPLX_DP &a, Vec_O_CPLX_DP &roots, const bool &polish)
{
	const DP EPS=1.0e-14;
	int i,its,j,jj;
	complex<DP> x,b,c;

	int m=a.size()-1;
	Vec_CPLX_DP ad(m+1);
	for (j=0;j<=m;j++) ad[j]=a[j];
	for (j=m-1;j>=0;j--) {
		x=0.0;
		Vec_CPLX_DP ad_v(j+2);
		for (jj=0;jj<j+2;jj++) ad_v[jj]=ad[jj];
		laguer(ad_v,x,its);
		if (fabs(imag(x)) <= 2.0*EPS*fabs(real(x)))
			x=complex<DP>(real(x),0.0);
		roots[j]=x;
		b=ad[j+1];
		for (jj=j;jj>=0;jj--) {
			c=ad[jj];
			ad[jj]=b;
			b=x*b+c;
		}
	}
	if (polish)
		for (j=0;j<m;j++)
			laguer(a,roots[j],its);
	for (j=1;j<m;j++) {
		x=roots[j];
		for (i=j-1;i>=0;i--) {
			if (real(roots[i]) <= real(x)) break;
			roots[i+1]=roots[i];
		}
		roots[i+1]=x;
	}
}

//алгоритм лаггера для нахожд корней полиномов
void NR::laguer(Vec_I_CPLX_DP &a, complex<DP> &x, int &its)
{
	const int MR=8,MT=10,MAXIT=MT*MR;
	const DP EPS=std::numeric_limits<DP>::epsilon();
	static const DP frac[MR+1]=
		{0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
	int iter,j;
	DP abx,abp,abm,err;
	complex<DP> dx,x1,b,d,f,g,h,sq,gp,gm,g2;

	int m=a.size()-1;
	for (iter=1;iter<=MAXIT;iter++) {
		its=iter;
		b=a[m];
		err=abs(b);
		d=f=0.0;
		abx=abs(x);
		for (j=m-1;j>=0;j--) {
			f=x*f+d;
			d=x*d+b;
			b=x*b+a[j];
			err=abs(b)+abx*err;
		}
		err *= EPS;
		if (abs(b) <= err) return;
		g=d/b;
		g2=g*g;
		h=g2-2.0*f/b;
		sq=sqrt(DP(m-1)*(DP(m)*h-g2));
		gp=g+sq;
		gm=g-sq;
		abp=abs(gp);
		abm=abs(gm);
		if (abp < abm) gp=gm;
		dx=MAX(abp,abm) > 0.0 ? DP(m)/gp : polar(1+abx,DP(iter));
		x1=x-dx;
		if (x == x1) return;
		if (iter % MT != 0) x=x1;
		else x -= frac[iter/MT]*dx;
	}
	nrerror("too many iterations in laguer");
	return;
}

//нахождение коэффицентов для полинома
void NR::polcoe(Vec_I_DP &x, Vec_I_DP &y, Vec_O_DP &cof)
{
	int k,j,i;
	DP phi,ff,b;

	int n=x.size();
	Vec_DP s(n);
	for (i=0;i<n;i++) s[i]=cof[i]=0.0;
	s[n-1]= -x[0];
	for (i=1;i<n;i++) {
		for (j=n-1-i;j<n-1;j++)
			s[j] -= x[i]*s[j+1];
		s[n-1] -= x[i];
	}
	for (j=0;j<n;j++) {
		phi=n;
		for (k=n-1;k>0;k--)
			phi=k*s[k]+x[j]*phi;
		ff=y[j]/phi;
		b=1.0;
		for (k=n-1;k>=0;k--) {
			cof[k] += b*ff;
			b=s[k]+x[j]*b;
		}
	}
}

// решение системы LU разложением

void NR::ludcmp(Mat_IO_DP &a, Vec_O_INT &indx, DP &d)
{
	const DP TINY=1.0e-20;
	int i,imax,j,k;
	DP big,dum,sum,temp;

	int n=a.nrows();
	Vec_DP vv(n);
	d=1.0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmp");
		vv[i]=1.0/big;
	}
	for (j=0;j<n;j++) {
		for (i=0;i<j;i++) {
			sum=a[i][j];
			for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<n;i++) {
			sum=a[i][j];
			for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ((dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=0;k<n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			d = -d;
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=TINY;
		if (j != n-1) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<n;i++) a[i][j] *= dum;
		}
	}
}
// решение линейной системы
void NR::lubksb(Mat_I_DP &a, Vec_I_INT &indx, Vec_IO_DP &b)
{
	int i,ii=0,ip,j;
	DP sum;

	int n=a.nrows();
	for (i=0;i<n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii != 0)
			for (j=ii-1;j<i;j++) sum -= a[i][j]*b[j];
		else if (sum != 0.0)
			ii=i+1;
		b[i]=sum;
	}
	for (i=n-1;i>=0;i--) {
		sum=b[i];
		for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}

void NR::rkqs(Vec_IO_DP &y, Vec_IO_DP &dydx, DP &x, const DP htry,
	const DP eps, Vec_I_DP &yscal, DP &hdid, DP &hnext,
	void derivs(const DP, Vec_I_DP &, Vec_O_DP &))
{
	const DP SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;
	int i;
	DP errmax,h,htemp,xnew;

	int n=y.size();
	h=htry;
	Vec_DP yerr(n),ytemp(n);
	for (;;) {
		rkck(y,dydx,x,h,ytemp,yerr,derivs);
		errmax=0.0;
		for (i=0;i<n;i++) errmax=MAX(errmax,fabs(yerr[i]/yscal[i]));
		errmax /= eps;
		if (errmax <= 1.0) break;
		htemp=SAFETY*h*pow(errmax,PSHRNK);
		h=(h >= 0.0 ? MAX(htemp,0.1*h) : MIN(htemp,0.1*h));
		xnew=x+h;
		if (xnew == x) nrerror("stepsize underflow in rkqs");
	}
	if (errmax > ERRCON) hnext=SAFETY*h*pow(errmax,PGROW);
	else hnext=5.0*h;
	x += (hdid=h);
	for (i=0;i<n;i++) y[i]=ytemp[i];
}

void NR::rkck(Vec_I_DP &y, Vec_I_DP &dydx, const DP x,
	const DP h, Vec_O_DP &yout, Vec_O_DP &yerr,
	void derivs(const DP, Vec_I_DP &, Vec_O_DP &))
{
	static const DP a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875,
		b21=0.2, b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42 = -0.9,
		b43=1.2, b51 = -11.0/54.0, b52=2.5, b53 = -70.0/27.0,
		b54=35.0/27.0, b61=1631.0/55296.0, b62=175.0/512.0,
		b63=575.0/13824.0, b64=44275.0/110592.0, b65=253.0/4096.0,
		c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0,
		dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0,
		dc4=c4-13525.0/55296.0, dc5 = -277.00/14336.0, dc6=c6-0.25;
	int i;

	int n=y.size();
	Vec_DP ak2(n),ak3(n),ak4(n),ak5(n),ak6(n),ytemp(n);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+b21*h*dydx[i];
	derivs(x+a2*h,ytemp,ak2);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
	derivs(x+a3*h,ytemp,ak3);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
	derivs(x+a4*h,ytemp,ak4);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
	derivs(x+a5*h,ytemp,ak5);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
	derivs(x+a6*h,ytemp,ak6);
	for (i=0;i<n;i++)
		yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
	for (i=0;i<n;i++)
		yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
}



