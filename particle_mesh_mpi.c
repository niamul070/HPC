#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#define ngrid 120 /*number of grid points in one direction*/
#define nparticle 10 /*number of particles*/
#define particleMass 1


int Gturn=0;

struct particle {
	double positionX;
	double positionY;
	double velocityX;
	double velocityY;
	double forceX;
	double forceY;
	int charge;
	double Ex;
	double Ey;
};


MPI_Datatype my_data_type;
double global_dt;
int totalElemPerProcess;
struct particle *chargeParticle;
struct particle *chargeParticleLeft;
struct particle *chargeParticleRight;
struct particle *chargeParticleUp;
struct particle *chargeParticleDown;
double *VnextX,*VnextY,*posNextX,*posNextY,*forceNextX,*forceNextY;
	
void init6Arrays(int size)
{
	VnextX = (double *) malloc (size * sizeof(double));
	VnextY = (double *) malloc (size * sizeof(double));
	posNextX = (double *) malloc (size * sizeof(double));
	posNextY = (double *) malloc (size * sizeof(double));
	forceNextX = (double *) malloc (size * sizeof(double));
	forceNextY = (double *) malloc (size * sizeof(double));
}

void swapParticle(int x, int y)
{
	struct particle temp = chargeParticle[x];
	chargeParticle[x] = chargeParticle[y];
	chargeParticle[y]= temp;
}
void initParticle(double xmax,double xmin,double ymax,double ymin,int rank){
	int i;
	for(i=0;i<totalElemPerProcess;i++){
			chargeParticle[i].positionX=xmin+((double)rand()/RAND_MAX)*(xmax-xmin) ;
			chargeParticle[i].positionY=ymin+((double)rand()/RAND_MAX)*(ymax-ymin) ;
			chargeParticle[i].velocityX=0.0;
			chargeParticle[i].velocityY=0.0;
			chargeParticle[i].forceX=0.0;
			chargeParticle[i].forceY=0.0;
			if(i%2==0) chargeParticle[i].charge=1;
			else 	chargeParticle[i].charge=-1;	
	}
}

void calcforces(double **delux,double **deluy,double xstart, double xend,double xincrement,double yincrement,int rank,double xmin,double xmax,double ymin,double ymax,int t,int up,int down,int left,int right,int xDim,int yDim,int istart,int nxlim,int jstart,int nylim){
	int forceI,forceJ;
	double kasi,neu;


	
	double dt,vmax;
	int i;
	if (t==1)
		dt=0.25; 
	else 
		dt=global_dt;
	
	dt=0.050;
	for(i=0;i<totalElemPerProcess;i++){
		forceJ=floor((chargeParticle[i].positionX-xmin)/xincrement);
		forceI=floor((chargeParticle[i].positionY-ymin)/yincrement);
		
		kasi=(chargeParticle[i].positionX-((double)forceJ*xincrement+xstart))/xincrement;
		neu=(chargeParticle[i].positionY-((double)forceI*yincrement+xstart))/yincrement;
		
		
		if(forceI<0 || (forceI+1)>=(yDim+2) || forceJ<0 || (forceJ+1)>=(xDim+2) ){
			swapParticle(i,totalElemPerProcess-1);
			totalElemPerProcess--;
			continue;
			
		}
		
		chargeParticle[i].Ex=(1-kasi)*(1-neu)*delux[forceI][forceJ]+kasi*(1-neu)*delux[forceI+1][forceJ]+kasi*neu*delux[forceI+1][forceJ+1]+(1-kasi)*neu*delux[forceI][forceJ+1];
		chargeParticle[i].Ey=(1-kasi)*(1-neu)*deluy[forceI][forceJ]+kasi*(1-neu)*deluy[forceI+1][forceJ]+kasi*neu*deluy[forceI+1][forceJ+1]+(1-kasi)*neu*deluy[forceI][forceJ+1];
		forceNextY[i]=(chargeParticle[i].charge)*(chargeParticle[i].Ey);
		posNextX[i]=chargeParticle[i].positionX+dt*chargeParticle[i].velocityX+0.5*dt*dt*forceNextX[i];
		posNextY[i]=chargeParticle[i].positionY+dt*chargeParticle[i].velocityY+0.5*dt*dt*forceNextY[i];		
		VnextX[i]=chargeParticle[i].velocityX+0.5*dt*((chargeParticle[i].forceX/particleMass)+forceNextX[i]);
		VnextY[i]=chargeParticle[i].velocityY+0.5*dt*((chargeParticle[i].forceY/particleMass)+forceNextY[i]);
		
	}
	
	
	
	vmax=sqrt(VnextX[0]*VnextX[0]+VnextY[0]*VnextY[0]);
	for(i=1;i<totalElemPerProcess;i++){
		if(vmax<sqrt(VnextX[i]*VnextX[i]+VnextY[i]*VnextY[i]))
			vmax=sqrt(VnextX[i]*VnextX[i]+VnextY[i]*VnextY[i]);
	}
	/*calculate dt */
	if(t !=1){
		dt=xincrement/vmax;
		global_dt=dt;
	}
	
	/*advance particles*/ 
	for(i=0;i<totalElemPerProcess;i++){
		chargeParticle[i].forceX=forceNextX[i];
		chargeParticle[i].forceY=forceNextY[i];
		chargeParticle[i].velocityX=VnextX[i];
		chargeParticle[i].velocityY=VnextY[i];
		chargeParticle[i].positionX=posNextX[i];
		chargeParticle[i].positionY=posNextY[i];
	}
	
	
	
		
	int varToSentRight;
	varToSentRight=0;
	int varToSentLeft;
	varToSentLeft=0;
	int varToSentUp;
	varToSentUp=0;
	int varToSentDown;
	varToSentDown=0;
	for (i=0;i<totalElemPerProcess;i++){
		
		
		if(chargeParticle[i].positionX < xstart || chargeParticle[i].positionX > xend ||
			chargeParticle[i].positionY < xstart || chargeParticle[i].positionY > xend)
		{
			swapParticle(i,totalElemPerProcess-1);
			totalElemPerProcess--;
			continue;
		}
		if(chargeParticle[i].positionX > xmax)
		{
			chargeParticleRight[varToSentRight] = chargeParticle[i];
			swapParticle(i,totalElemPerProcess-1);
			varToSentRight++;
			totalElemPerProcess--;
		}
		if(chargeParticle[i].positionX < xmin)
		{
			chargeParticleLeft[varToSentLeft] = chargeParticle[i];
			swapParticle(i,totalElemPerProcess-1);
			varToSentLeft++;
			totalElemPerProcess--;
		}
		if(chargeParticle[i].positionY > ymax)
		{
			chargeParticleUp[varToSentUp] = chargeParticle[i];
			swapParticle(i,totalElemPerProcess-1);
			varToSentUp++;
			totalElemPerProcess--;
		}
		if(chargeParticle[i].positionY < ymin)
		{
			chargeParticleDown[varToSentDown] = chargeParticle[i];
			swapParticle(i,totalElemPerProcess-1);
			varToSentDown++;
			totalElemPerProcess--;
		}
		
	}

		if(up>-1){
		
			int receiveFromUp;
			MPI_Recv(&receiveFromUp,1,MPI_INT,up,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			int tt;
			for(tt=0;tt<receiveFromUp;tt++)
			{
				struct particle tempP;
				MPI_Recv(&tempP,1,my_data_type,up,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				chargeParticle[totalElemPerProcess-1]=tempP;
				totalElemPerProcess++;
			}
			
			MPI_Send(&varToSentUp,1,MPI_INT,up,0, MPI_COMM_WORLD);
			for(tt=0;tt<varToSentUp;tt++)
			{
				struct particle tempP=chargeParticleUp[tt];
				MPI_Send(&tempP,1,my_data_type,up,0, MPI_COMM_WORLD);
			}
			
			
			
		}
	
	
		if(down>-1){
		
			int receiveFromDown;
			int tt;
			MPI_Send(&varToSentDown,1,MPI_INT,down,0, MPI_COMM_WORLD);
			for(tt=0;tt<varToSentDown;tt++)
			{
				struct particle tempP=chargeParticleDown[tt];
				MPI_Send(&tempP,1,my_data_type,down,0, MPI_COMM_WORLD);
			}
			
			
			MPI_Recv(&receiveFromDown,1,MPI_INT,down,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(tt=0;tt<receiveFromDown;tt++)
			{
				struct particle tempP;
				MPI_Recv(&tempP,1,my_data_type,down,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				chargeParticle[totalElemPerProcess-1]=tempP;
				totalElemPerProcess++;
			}
			
		
		}
	

	
		if(left>-1){
		
			int receiveFromLeft;
			MPI_Recv(&receiveFromLeft,1,MPI_INT,left,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			int tt;
			for(tt=0;tt<receiveFromLeft;tt++)
			{
				struct particle tempP;
				MPI_Recv(&tempP,1,my_data_type,left,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				chargeParticle[totalElemPerProcess-1]=tempP;
				totalElemPerProcess++;
			}
			
			MPI_Send(&varToSentLeft,1,MPI_INT,left,0, MPI_COMM_WORLD);
			for(tt=0;tt<varToSentLeft;tt++)
			{
				struct particle tempP=chargeParticleLeft[tt];
				MPI_Send(&tempP,1,my_data_type,left,0, MPI_COMM_WORLD);
			}
			
			
			
		}
	
		if(right>-1){
		
			int receiveFromRight;
			int tt;
			MPI_Send(&varToSentRight,1,MPI_INT,right,0, MPI_COMM_WORLD);
			for(tt=0;tt<varToSentRight;tt++)
			{
				struct particle tempP=chargeParticleRight[tt];
				MPI_Send(&tempP,1,my_data_type,right,0, MPI_COMM_WORLD);
			}
			
			
			MPI_Recv(&receiveFromRight,1,MPI_INT,right,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(tt=0;tt<receiveFromRight;tt++)
			{
				struct particle tempP;
				MPI_Recv(&tempP,1,my_data_type,right,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				chargeParticle[totalElemPerProcess-1]=tempP;
				totalElemPerProcess++;
			}
			
		
		}
	
	
	
	if(totalElemPerProcess<nparticle){
		for(i=totalElemPerProcess;i<nparticle;i++){
			chargeParticle[i].positionX=xmin+((double)rand()/RAND_MAX)*(xmax-xmin) ;
			chargeParticle[i].positionY=ymin+((double)rand()/RAND_MAX)*(ymax-ymin) ;
			chargeParticle[i].velocityX=0.0;
			chargeParticle[i].velocityY=0.0;
			chargeParticle[i].forceX=0.0;
			chargeParticle[i].forceY=0.0;
			if(i%2==0) chargeParticle[i].charge=1;
			else 	chargeParticle[i].charge=-1;
		
		}
		totalElemPerProcess=nparticle;
	}
	

	MPI_Barrier(MPI_COMM_WORLD);
	printf("\nTotal Turn=%d",Gturn++);
	
	
	if(t%10==0){
		char str1[15],str2[15],str3[15],str4[15];
		FILE *fp;
		sprintf(str2, "%d", rank);
		sprintf(str3,"%d",t);
		strcat(str3,str2);
		strcpy (str1,"pos_");
		strcat (str1,str3);
		fp = fopen(str1, "w");
	
		for(i=0;i<totalElemPerProcess;i++){
			fprintf(fp,"%0.4f	%0.4f\n",chargeParticle[i].positionX,chargeParticle[i].positionY);				
		
		}
		fclose(fp);
	}
	
	

	
}

int nodecalc(int istart,int nxlim,int jstart,int nylim,int rank,int rowNumber,int colNumber,int nrow,int ncol,double delx, double dely,double xstart,double xend,int t){
	
	
	
	
	int xDim = nxlim-istart;
	int yDim = nylim-jstart;
	int i,j,lenhight,lenwidth,left,right,up,down;
	i=istart;
	j=jstart;
	double u[xDim+2][yDim+2],u_new[xDim+2][yDim+2],f[xDim+2][yDim+2];
	double h,xmin,xmax,xi,eta,ymin,ymax;
	double sendLeft[yDim],recvLeft[yDim],sendRight[yDim],recvRight[yDim];
	double sendUp[xDim],recvUp[xDim],sendDown[xDim],recvDown[xDim];
	
	double **delux,**deluy;
	
	delux = (double **) malloc( (xDim+2)*sizeof(double *));
	deluy = (double **) malloc( (xDim+2)*sizeof(double *));
	
	int xxx;
	for(xxx= 0;xxx<xDim+2;xxx++)
	{
		delux[xxx]=(double *) malloc ((yDim+2)* sizeof(double));
		deluy[xxx]=(double *) malloc ((yDim+2)* sizeof(double));
	}
	
	
	int p;
	
	
	h=delx;
	
	left=findProcessNumberFromRowCol(rowNumber,colNumber-1,nrow,ncol);
	right=findProcessNumberFromRowCol(rowNumber,colNumber+1,nrow,ncol);
	up=findProcessNumberFromRowCol(rowNumber-1,colNumber,nrow,ncol);
	down=findProcessNumberFromRowCol(rowNumber+1,colNumber,nrow,ncol);

	/*array initialization */
	for(i=0;i<(nxlim+2-istart);i++){
		for(j=0;j<(nylim+2-jstart);j++){			
				u_new[i][j] = 0.0;
				u[i][j]=0.0;
				f[i][j]=0.0;
		}
	}
	
	
	
	
	/*physical boundary:each process will initialize its segment  of initial boundary */
	
	xmin=((double)colNumber/(double)ncol)*(xend-xstart)+xstart;
	xmax=((double)(colNumber+1)/(double)ncol)*(xend-xstart)+xstart;
	ymin=((double)rowNumber/(double)nrow)*(xend-xstart)+xstart; /*xend,xstart same for X,Y direction since domain is square */
	ymax=((double)(rowNumber+1)/(double)nrow)*(xend-xstart)+xstart;
	double xincrement,yincrement;
	xincrement=(xmax-xmin)/((double)(nxlim-istart));
	yincrement=(ymax-ymin)/((double)(nylim-jstart));
	
	
	/*calculate the density here. set flag based on t for random initialization or natural advancement */
	
	if (t==1){/*set random particles with zero velocities */
		initParticle(xmax,xmin,ymax,ymin,rank);
		
	}
	
	int particleI,particleJ,k;
	double chargeSum;
	
	/*loop over f for each i,j loop over particles and assign particle density to grid */
	for(i=0;i<(nxlim+2-istart);i++){
		for(j=0;j<(nylim+2-jstart);j++){
			chargeSum=0.0;
			for(k=0;k<totalElemPerProcess;k++){
				particleI=chargeParticle[k].positionX/xincrement;
				particleJ=chargeParticle[k].positionY/yincrement;
				if(particleI==i && particleJ==j){
					chargeSum=chargeSum+chargeParticle[k].charge;
					f[i][j]=chargeSum/(nparticle*xincrement*yincrement);
				}
				
			}
		}
	}

	if(istart==0) { /*define left boundary */
		for(j=1;j<(nylim+2-jstart-1);j++) {
			u_new[j][0+1]=-1.0;
			u[j][0+1]=-1.0;

		}	
	}
	
	if(nxlim==ngrid) { /*define right boundary */
		for(j=1;j<(nylim+2-jstart-1);j++) {
			u_new[j][nxlim-istart+1-1]=1.0;
			u[j][nxlim-istart+1-1]=1.0;

			
		}	
	}
	
	if(nylim==ngrid) { /*define top boundary */
		
		for(j=1;j<(nxlim+2-istart-1);j++) {
			xi=(j)*xincrement+xmin;
			eta=(xi-xstart)/(xend-xstart);
			u_new[0+1][j]=(1-eta)*(-1.0)+1.0*eta;  /*adding 1 for the purpose of cusion */

		}	
	}
	
	if(jstart==0) { /*define bottom boundary */
		
		
		for(j=1;j<(nxlim+2-istart);j++) {
			xi=(j)*xincrement+xmin;
			eta=(xi-xstart)/(xend-xstart);
			u_new[nylim-jstart+1-1][j]=(1-eta)*-1.0+1.0*eta;

		}	
	}
			
	/*end physical boundary*/
	
	int rowStart,colStart,rowEnd,colEnd;
	rowStart=0;
	colStart=0;
	rowEnd=0;
	colEnd=0;
	
	/*bottom left */
	if(istart==0 && jstart==0){rowEnd=1;colStart=1;}
	/*bottom right*/
	else if (nxlim==ngrid && jstart==0){rowEnd=1;colEnd=1;}
	/*top left */
	else if (istart==0 && nylim==ngrid){rowStart=1;colStart=1;}
	/*top right */
	else if(nxlim==ngrid && nylim==ngrid){rowStart=1;colEnd=1;}
	/*left */
	else if (istart==0){colStart=1;}
	/*right*/
	else if (nxlim==ngrid){colEnd=1;}
	/*top*/
	else if (nylim==ngrid){rowStart=1;}
	/*bottom*/
	else if (jstart==0){rowEnd=1;}
	else{
		
	}
	

	for(p=0;p<1000000;p++){
		/*calculation all interior nodes for each proc are done here */
		for(i=1+rowStart;i<(nxlim-istart+1-rowEnd);i++){
			for(j=1+colStart;j<(nylim-jstart+1-colEnd);j++){			
					u_new[i][j] = 0.25 * (u[i+1][j] + u[i-1][j] +u[i][j+1] + u[i][j-1] - h*h*f[i][j]);	
			}
		}
	
		/* check if process has physical boundary */
	
	
		/*define communication strategy:this is inside while loop*/	
	
	
		if(up>-1){
		
			MPI_Recv(recvUp,xDim,MPI_DOUBLE,up,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(i=0;i<xDim;i++){
				u_new[0][i+1]=recvUp[i];
			}
		
			for(i=0;i<xDim;i++){
				sendUp[i]=u_new[1][i+1];
			}
			MPI_Send(sendUp,xDim,MPI_DOUBLE,up,0, MPI_COMM_WORLD);
		}
	
		if(down>-1){
			for(i=0;i<xDim;i++){
				sendDown[i]=u_new[yDim][i+1];
			}		
			MPI_Send(sendDown,xDim,MPI_DOUBLE,down,0, MPI_COMM_WORLD);
		
			MPI_Recv(recvDown,xDim,MPI_DOUBLE,down,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(i=0;i<xDim;i++){
				u_new[yDim+1][i+1]=recvDown[i];
			}
		
		}
	
		if(left>-1){
		
			MPI_Recv(recvLeft,yDim,MPI_DOUBLE,left,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(i=0;i<yDim;i++){
				u_new[i+1][0]=recvLeft[i];
			}
			for(i=0;i<yDim;i++){
				sendLeft[i]=u_new[i+1][1];
			}
			MPI_Send(sendLeft,yDim,MPI_DOUBLE,left,0, MPI_COMM_WORLD);
		}
	
		if(right>-1){
		
			for(i=0;i<yDim;i++){

				sendRight[i]=u_new[i+1][xDim];
				
			}
			MPI_Send(sendRight,yDim,MPI_DOUBLE,right,0, MPI_COMM_WORLD);
			MPI_Recv(recvRight,yDim,MPI_DOUBLE,right,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(i=0;i<yDim;i++){
				u_new[i+1][xDim+1]=recvRight[i];

			}
		
		}
		
				

		for(i=0;i<(nxlim+2-istart);i++){
			for(j=0;j<(nylim+2-jstart);j++){
	   			   		u[i][j]=u_new[i][j];
	   		}
	   }
	}
	
	
	
	
	/* after field is calculated (u) calculate the gradient */
	for(i=1;i<(nxlim-istart+1);i++){
		for(j=1;j<(nylim-jstart+1);j++){			
			

				
				if(i==1 && j==1){ /*lower left corner*/
					delux[i][j]=(u[i][j+1]-u[i][j])/(xincrement);
					deluy[i][j]=(u[i+1][j]-u[i][j])/(xincrement);
				}
					
		
				else if(i==1 && j==(nylim-jstart)){ /*lower right corner*/
					delux[i][j]=(u[i][j]-u[i][j-1])/(xincrement);
					deluy[i][j]=(u[i+1][j]-u[i][j])/(xincrement);
				}
				else if (j==(nylim-jstart) && i==(nxlim-istart)){ /*top right*/
					delux[i][j]=(u[i][j]-u[i][j-1])/xincrement;
					deluy[i][j]=(u[i][j]-u[i-1][j])/xincrement;
				}
				else if (i==(nxlim-istart) && j==1){ /*top left*/
					delux[i][j]=(u[i][j+1]-u[i][j])/xincrement; 
					deluy[i][j]=(u[i][j]-u[i-1][j])/xincrement;
				}
				
				else if (j==1 && (i>1 || i<(nxlim-istart))){/*left side*/ 
					delux[i][j]=(u[i][j+1]-u[i][j])/xincrement;
					deluy[i][j]=(u[i+1][j]-u[i-1][j])/(2*xincrement);
				}
				else if (j==(nylim-jstart) && (i>1 || i<(nxlim-istart))){/*right side*/ 
					delux[i][j]=(u[i][j]-u[i][j-1])/xincrement;
					deluy[i][j]=(u[i+1][j]-u[i-1][j])/(2*xincrement);
				}
				else if (i==1 && (j>1 || j<(nylim-jstart))){/*bottom side */
					delux[i][j]=(u[i][j+1]-u[i][j-1])/(2*xincrement);
					deluy[i][j]=(u[i+1][j]-u[i][j])/xincrement;
				}
				else if (i==(nxlim-istart)&& (j>1 || j<(nylim-jstart))){/*top side*/ 
					delux[i][j]=(u[i][j+1]-u[i][j-1])/(2*xincrement);
					deluy[i][j]=(u[i][j]-u[i-1][j])/xincrement;
				}			
				
		
				
				else{
					delux[i][j]=(u[i][j+1]-u[i][j-1])/(2*xincrement);
					deluy[i][j]=(u[i+1][j]-u[i-1][j])/(2*xincrement);
				}
				
						
			}
	}
	
	
	
	if(up>-1){
		
			MPI_Recv(recvUp,xDim,MPI_DOUBLE,up,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(i=0;i<xDim;i++){
				delux[0][i+1]=recvUp[i];
			}
		
			for(i=0;i<xDim;i++){
				sendUp[i]=delux[1][i+1];
			}
			MPI_Send(sendUp,xDim,MPI_DOUBLE,up,0, MPI_COMM_WORLD);
		}
	
		if(down>-1){
			for(i=0;i<xDim;i++){
				sendDown[i]=delux[yDim][i+1];
			}		
			MPI_Send(sendDown,xDim,MPI_DOUBLE,down,0, MPI_COMM_WORLD);
		
			MPI_Recv(recvDown,xDim,MPI_DOUBLE,down,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(i=0;i<xDim;i++){
				delux[yDim+1][i+1]=recvDown[i];
			}
		
		}
	
		if(left>-1){
		
			MPI_Recv(recvLeft,yDim,MPI_DOUBLE,left,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(i=0;i<yDim;i++){
				delux[i+1][0]=recvLeft[i];
			}
			for(i=0;i<yDim;i++){
				sendLeft[i]=delux[i+1][1];
			}
			MPI_Send(sendLeft,yDim,MPI_DOUBLE,left,0, MPI_COMM_WORLD);
		}
	
		if(right>-1){
		
			for(i=0;i<yDim;i++){

				sendRight[i]=delux[i+1][xDim];
				
			}
			MPI_Send(sendRight,yDim,MPI_DOUBLE,right,0, MPI_COMM_WORLD);
			MPI_Recv(recvRight,yDim,MPI_DOUBLE,right,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(i=0;i<yDim;i++){
				delux[i+1][xDim+1]=recvRight[i];

			}
		
		}
		
		if(up>-1){
		
			MPI_Recv(recvUp,xDim,MPI_DOUBLE,up,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(i=0;i<xDim;i++){
				deluy[0][i+1]=recvUp[i];
			}
		
			for(i=0;i<xDim;i++){
				sendUp[i]=deluy[1][i+1];
			}
			MPI_Send(sendUp,xDim,MPI_DOUBLE,up,0, MPI_COMM_WORLD);
		}
	
		if(down>-1){
			for(i=0;i<xDim;i++){
				sendDown[i]=deluy[yDim][i+1];
			}		
			MPI_Send(sendDown,xDim,MPI_DOUBLE,down,0, MPI_COMM_WORLD);
		
			MPI_Recv(recvDown,xDim,MPI_DOUBLE,down,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(i=0;i<xDim;i++){
				deluy[yDim+1][i+1]=recvDown[i];
			}
		
		}
	
		if(left>-1){
		
			MPI_Recv(recvLeft,yDim,MPI_DOUBLE,left,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(i=0;i<yDim;i++){
				deluy[i+1][0]=recvLeft[i];
			}
			for(i=0;i<yDim;i++){
				sendLeft[i]=deluy[i+1][1];
			}
			MPI_Send(sendLeft,yDim,MPI_DOUBLE,left,0, MPI_COMM_WORLD);
		}
	
		if(right>-1){
		
			for(i=0;i<yDim;i++){

				sendRight[i]=deluy[i+1][xDim];
				
			}
			MPI_Send(sendRight,yDim,MPI_DOUBLE,right,0, MPI_COMM_WORLD);
			MPI_Recv(recvRight,yDim,MPI_DOUBLE,right,0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			for(i=0;i<yDim;i++){
				deluy[i+1][xDim+1]=recvRight[i];

			}
		
		}
	/*now we have the gradient. Next step will be to calculate the force on each particles */
	calcforces(delux,deluy,xstart,xend,xincrement,yincrement,rank,xmin,xmax,ymin,ymax,t,up,down,left,right,xDim,yDim,istart,nxlim,jstart,nylim);
	
	
	
return 0;
}



int findProcessNumberFromRowCol(int row,int col, int totalRow, int totalCol){
	if(row<0 || col <0 || row>=totalRow || col>=totalCol){
		return -1;
	}
	else {
		return row*totalCol+col;
	}
}


int main(){

	double y,delx,dely,x1,y1,x2,y2,xstart,xend,ystart,yend;
	double a;
	double b;
	int i,j,my_rank,comm_sz,source,ceil_count,floor_count;
	double maxx,maxy;
	double local_a,local_b,gmint,mint,xu,yu,xk,yk;
	double xk1,xk2,yk1,yk2;
	int count,nrow,ncol,nproc,x,rk,ck,rank,nxlim,nylim,sum;
	int up_proc,down_proc,left_proc,right_proc,t;

	xstart = -1; xend = 1;
	ystart = -1; yend = 1;
	maxx=-1e15;
	maxy=maxx;
	delx=(xend-xstart)/(ngrid);
	dely=(yend-ystart)/(ngrid);
	
	
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	rank=my_rank;
	
	int blocklengths[9] = {1,1,1,1,1,1,1,1,1};

	struct particle p;
 	MPI_Datatype types[9]={MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_INT,MPI_DOUBLE,MPI_DOUBLE};
	
    MPI_Aint  offsets[9];

    
	MPI_Aint mpiParticleDesc, mpiParticleDesc1;
    MPI_Get_address(&p, &mpiParticleDesc);

    MPI_Get_address(&p.positionX, &mpiParticleDesc1);
    offsets[0] = mpiParticleDesc1 - mpiParticleDesc;
    MPI_Get_address(&p.positionY, &mpiParticleDesc1);
    offsets[1] = mpiParticleDesc1 - mpiParticleDesc;

    MPI_Get_address(&p.velocityX, &mpiParticleDesc1);
    offsets[2] = mpiParticleDesc1 - mpiParticleDesc;
    MPI_Get_address(&p.velocityY, &mpiParticleDesc1);
    offsets[3] = mpiParticleDesc1 - mpiParticleDesc;

    MPI_Get_address(&p.forceX, &mpiParticleDesc1);
    offsets[4] = mpiParticleDesc1 - mpiParticleDesc;
    MPI_Get_address(&p.forceY, &mpiParticleDesc1);
    offsets[5] = mpiParticleDesc1 - mpiParticleDesc;

    MPI_Get_address(&p.charge, &mpiParticleDesc1);
    offsets[6] = mpiParticleDesc1 - mpiParticleDesc;
    MPI_Get_address(&p.Ex, &mpiParticleDesc1);
    offsets[7] = mpiParticleDesc1 - mpiParticleDesc;

    MPI_Get_address(&p.Ey, &mpiParticleDesc1);
    offsets[8] = mpiParticleDesc1 - mpiParticleDesc;
    

    MPI_Type_create_struct(9, blocklengths, offsets, types, &my_data_type);
    MPI_Type_commit(&my_data_type);
	
	
	
	
	
	

	count=0;
	nproc=comm_sz;
	x=nproc;
	chargeParticle=(struct particle *)malloc(nparticle*nproc*sizeof(struct particle));
	chargeParticleLeft=(struct particle *)malloc(nparticle*nproc*sizeof(struct particle));
	chargeParticleRight=(struct particle *)malloc(nparticle*nproc*sizeof(struct particle));
	chargeParticleUp=(struct particle *)malloc(nparticle*nproc*sizeof(struct particle));
	chargeParticleDown=(struct particle *)malloc(nparticle*nproc*sizeof(struct particle));
	totalElemPerProcess=nparticle;
	init6Arrays(nparticle*nproc);	
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


	xu=((double)ngrid)/ncol;
	yu=((double)ngrid)/nrow;
	
	/*xk and yk are starting points at each process*/
	x1=0;/*going right from grid*/

	y2=ngrid;/*going down from grid */

	xk=ck*xu+x1;
	yk=y2-rk*yu;


	/*xk1 and yk1 are ending points of each process */
	xk1=xk;
	yk1=yk;
	xk2=xk+xu;
	yk2=yk-yu; 

	xk1=floor(xk1);
	yk1=floor(yk1);
	xk2=floor(xk2);
	yk2=floor(yk2);

	

	i=(int)xk1;
	j=(int)yk2;
	nxlim=(int)xk2;
	nylim=(int)yk1;
	
	/*time loop should call poisson solver at each timestep */
	/*time loop will calculate charge density at each timestep */
	
	/*calculation of density, solution of poisson equation and calculation of gradient is done in 'nodecalc' */
	t=1; /*initial timestep*/
	for(t=1;t<100;t++){
		nodecalc(i,nxlim,j,nylim,rank,rk,ck,nrow,ncol,delx,dely,xstart,xend,t);
	}
	
	/* poisson solver will return the gradient array */
	/* with the help of gradient array interpolation of force to particle position will be calculated */
	/*from force advancement of particles will be done through verlet algorithm */
	/*minimum delta t will be calculated */
	
	/*finally cross boundary particles will be interchanged */ 
	
	/*go to next time step */



	MPI_Finalize();

	return 0;
}






