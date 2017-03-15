#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#define ngrid 100 /*this is number of random particle*/
#define MIN(a,b) (((a)<(b)) ? (a):(b))

double particle(double x1,double x2, double y1, double y2,int rank);

int compare(const void *a,const void *b) {
double *x = (double *) a;
double *y = (double *) b;
if (*x < *y) return -1;
else if (*x > *y) return 1; return 0;
}



int main(){

	double y,delx,dely,x1,y1,x2,y2;
	double a;
	double b;
	int i,j,my_rank,comm_sz,source,ceil_count,floor_count;
	double maxx,maxy;
	double local_a,local_b,gmint,mint,xu,yu,xk,yk;
	double xk1,xk2,yk1,yk2;
	int local_n,count,nrow,ncol,nproc,x,rk,ck,rank;

	x1 = -1; x2 = 1;
	y1 = -1; y2 = 1;
	maxx=-1e32;
	maxy=maxx;


	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	local_n =(x2-x1)/(comm_sz);
	comm_sz=comm_sz;
	rank=my_rank;

	count=0;
	nproc=comm_sz;
	x=nproc;

	while(x%2==0){
		x=x/2;
		count++; 	
	}

	ceil_count=ceil((float)count/2);
	floor_count=floor((float)count/2);
	nrow=pow(2,ceil_count);
	ncol=pow(2,floor_count)*x;
	rk=floor((float)rank/ncol);
	ck=(rank)%ncol;

	xu=(x2-x1)/ncol;
	yu=(y2-y1)/nrow;
	/*xk and yk are starting points at each process*/

	xk=ck*xu+x1;
	yk=y2-rk*yu;


	/*xk1 and yk1 are ending points of each process */
	xk1=xk+xu;
	yk1=yk+yu;
	xk2=xk+xu;
	yk2=yk-yu; 




	mint=particle(xk,xk2,yk2,yk,my_rank);
	MPI_Reduce(&mint,&gmint,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);



	if (my_rank==0) printf("\nfrom all process minimum t is %f\n",gmint);

	MPI_Finalize();
	return 0;
}


double particle(double x1,double x2, double y1, double y2,int rank){
	
	double location[2*ngrid];
	double velx[ngrid];
	double vely[ngrid],delt[ngrid],x[ngrid],y[ngrid];
	double a,b,mint,randx,randy,mindx,denominator,denominator1,numerator;
	int i,j;
	char str2[15],str1[15];
	double min2=1e15, min1=1e15;


	FILE *fp;

	sprintf(str2, "%d", rank);
	strcpy (str1,"hwk4_");
	strcat (str1,str2);

	fp = fopen(str1, "w");
	a=1.0;
	b=0.5;
	srand(time(NULL)*2*rank);



	for (i=0;i<ngrid;i++){
	
		randx=((double)rand() /RAND_MAX)*(x2-x1+1);
		location[i]=x1+((double)rand()/RAND_MAX)*(x2-x1) ;
		if(location[i]>-0.09 && location[i]<0.09) location[i]=location[i]+0.01;
		randy=((double)rand() / (double)RAND_MAX);
		location[i+ngrid]=y1+randy*(y2-y1);
		if(location[i+ngrid]>-0.09 && location[i+ngrid]<0.09) location[i+ngrid]=location[i+ngrid]+0.01;
	
		denominator1 = (location[i] * location[i]) + (location[i+ngrid] * location[i+ngrid]);
		denominator  = denominator1 * denominator1;
		numerator    = -2*b*location[i]*location[i+ngrid];
		vely[i]= numerator/denominator;
		velx[i]= a + b*((location[i+ngrid]*location[i+ngrid])-(location[i]*location[i]))/denominator;
		x[i]=location[i];
		y[i]=location[i+ngrid];
		fprintf(fp,"%.8f %.8f %.8f %.8f\n", location[i],location[i+ngrid],velx[i],vely[i]);
		delt[i]=(sqrt(velx[i]*velx[i]+vely[i]*vely[i]));
/*		printf("\nmin1 =%f",delt[i]); */
	}
	
	qsort(x,ngrid,sizeof(double),compare);
	qsort(y,ngrid,sizeof(double),compare);
	qsort (delt, ngrid, sizeof(double), compare);

	mint=MIN((x[1]-x[0]),(y[1]-y[0]))/delt[0];
	fclose(fp);
	return mint;
}


