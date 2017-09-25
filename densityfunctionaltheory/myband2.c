#include<stdio.h>
#include<stdlib.h>
#include<math.h>

extern void dsyev( char* jobz, char* uplo, int* n, double* a, int* lda,
		                                double* w, double* work, int* lwork, int* info );
int main(){

	double a=10.26;
	int nk;
	double pi=3.14159265358979;
	double tpi=2.0*pi;
	double ecut=30.0;
	double h1[3],h2[3],h3[3];
	h1[0]=1.0;h1[1]=1.0;h1[2]=-1.0;
	h2[0]=1.0;h2[1]=-1.0;h2[2]=1.0;
	h3[0]=-1.0;h3[1]=1.0;h3[2]=1.0;

	int nm1,nm2,nm3;
	int nmax=30;
	double gx,gy,gz;
	nk=1;
	double kx[nk],ky[nk],kz[nk];
	kx[0]=0.125;ky[0]=0.125;kz[0]=0.0;
	double kgx,kgy,kgz,kg2;
        int npw=0;
	int i,n;double k[3];
	for(n=0;n<nk;n++){
		for(nm1=-nmax;nm1<=nmax;nm1++){
			for(nm2=-nmax;nm2<=nmax;nm2++){
				for(nm3=-nmax;nm3<=nmax;nm3++){
				//	kg2=0.0;
				//	for(i=0;i<3;i++){
						k[0]=kx[n];k[1]=ky[n];k[2]=kz[n];
						kgx=k[0]+nm1*h1[0]+nm2*h2[0]+nm3*h3[0];
						kgy=k[1]+nm1*h1[1]+nm2*h2[1]+nm3*h3[1];
						kgz=k[2]+nm1*h1[2]+nm2*h2[2]+nm3*h3[2];
						kg2=pow(kgx,2)+pow(kgy,2)+pow(kgz,2);//h1,h2 and h3 orthonormal na o hote pare
				//	}
					kg2=pow(tpi/a,2)*kg2;
					if(kg2<=ecut){
						npw++;
					}
				}
			}
		}
	}

	//printf("%i\n",npw);
	//double kgx,kgy,kgz;
	double *kgx1,*kgy1,*kgz1,*psix;
	kgx1=(double*) malloc(npw*sizeof(double));
	kgy1=(double*) malloc(npw*sizeof(double));
	kgz1=(double*) malloc(npw*sizeof(double));
	psix=(double*) malloc(npw*sizeof(double));
	int ij;
	ij=0;
	for(n=0;n<nk;n++){
		for(nm1=-nmax;nm1<=nmax;nm1++){
			for(nm2=-nmax;nm2<=nmax;nm2++){
				for(nm3=-nmax;nm3<=nmax;nm3++){
					//kg2=0.0;
					k[0]=kx[n];k[1]=ky[n];k[2]=kz[n];
					kgx=k[0]+nm1*h1[0]+nm2*h2[0]+nm3*h3[0];
					kgy=k[1]+nm1*h1[1]+nm2*h2[1]+nm3*h3[1];
					kgz=k[2]+nm1*h1[2]+nm2*h2[2]+nm3*h3[2];
					kg2=(pow(kgx,2)+pow(kgy,2)+pow(kgz,2));
					kg2=pow(tpi/a,2)*kg2;
					if(kg2<=ecut){
						kgx1[ij]=kgx;
						kgy1[ij]=kgy;
						kgz1[ij]=kgz;
						psix[ij]=kgx;
						ij++;
					}
				}
			}
		}
	}
       //printf("%i\n",ij);
        int t=0;
        double tau1,tau2,tau3;
	tau1=0.125;tau2=0.125;tau3=0.125;
	int j; double g2;
	for(i=0;i<npw;i++){
		for(j=0;j<npw;j++){
			gx=kgx1[i]-kgx1[j];
			gy=kgy1[i]-kgy1[j];
			gz=kgz1[i]-kgz1[j];
			g2=gx*gx+gy*gy+gz*gz;
		//	vsg[t]=cos(gx*tau1+gy*tau2+gz*tau3);
		}
	}
       double ff;
	double *vsg; double *gvec;double* ham;
	double *xx,*yy,*zz;
	vsg=(double*) malloc(npw*npw*sizeof(double));
	ham=(double*) malloc(npw*npw*sizeof(double));
	xx=(double*) malloc(npw*npw*sizeof(double));
	yy=(double*) malloc(npw*npw*sizeof(double));
	zz=(double*) malloc(npw*npw*sizeof(double));
	gvec=(double*) malloc(npw*npw*sizeof(double));
	int *info;
	for(i=0;i<npw;i++){
	//	ham[i]=0.0;
	//	xx[i]=0.0;
	//	yy[i]=0.0;
	//	zz[i]=0.0;
//	printf("%i %lf\n",i,psix[i]);
	}

	        for(i=0;i<npw;i++){
			for(j=0;j<npw;j++){
				gx=kgx1[i]-kgx1[j];
				gy=kgy1[i]-kgy1[j];
				gz=kgz1[i]-kgz1[j];
				g2=gx*gx+gy*gy+gz*gz;
				xx[t]=kgx1[i];yy[t]=kgy1[i];zz[t]=kgz1[i];
				if(fabs(g2-3)<1.0e-6){
					ff=-0.21;
				}
				else if(fabs(g2-4)<1.0e-6){
					ff=0.0;			               
				}
				else if(fabs(g2-8)<1.0e-6){
					ff=0.04;
				}
				else if(fabs(g2-11)<1.0e-6){
					ff=0.08;
				}
				else{
					ff=0.0;
				}
				if(i==j){
					ham[t]=pow(tpi/a,2)*(kgx1[i]*kgx1[i]+kgy1[i]*kgy1[i]+kgz1[i]*kgz1[i])+ff*cos(tpi*(gx*tau1+gy*tau2+gz*tau3));
//					xx[t]=gx;yy[t]=gy;zz[t]=gz;
				}
					else{

			       	ham[t]=ff*cos(tpi*(gx*tau1+gy*tau2+gz*tau3));
//				xx[t]=gx;yy[t]=gy;zz[t]=gz;
					}
				t++;
			}
		}

	
	double* e;double* work; int lwork;	
	e   = (double *) malloc ( npw * sizeof (double) );
	lwork = 3*npw;
	work= (double *) malloc ( lwork*sizeof (double) );

		dsyev( "V", "U", &npw, ham, &npw, e, work, &lwork, &info );
       int cal=0;
         FILE* enband;
	 enband=fopen("bandSi","w+");
	       for(j=0;j<8;j++){
	       fprintf(enband,"%lf\n",e[j]*13.6058);
       }
	       fclose(enband);


	return 0;
}

