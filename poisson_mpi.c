#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#define ngrid 32 /*number of grid points in one direction*/



int nodecalc(int istart,int nxlim,int jstart,int nylim,int rank,int rowNumber,int colNumber,int nrow,int ncol,double delx, double dely,double xstart,double xend){
	
	int xDim = nxlim-istart;
	int yDim = nylim-jstart;
	int i,j,lenhight,lenwidth,left,right,up,down;
	i=istart;
	j=jstart;
	double u[xDim+2][yDim+2],u_new[xDim+2][yDim+2],f[xDim+2][yDim+2],delux[xDim+2][yDim+2],deluy[xDim+2][yDim+2];
	double h,xmin,xmax,xi,eta,minimum,diff1;
	double sendLeft[yDim],recvLeft[yDim],sendRight[yDim],recvRight[yDim];
	double sendUp[xDim],recvUp[xDim],sendDown[xDim],recvDown[xDim];
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
	
	double xincrement;
	xincrement=(xmax-xmin)/((double)(nxlim-istart));
	
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
			u_new[0+1][j]=(1-eta)*(-1.0)+1.0*eta;  		}	
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
	
	
	for(p=0;p<10000;p++){
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
	FILE *fp;
	char str1[15],str2[15];
	sprintf(str2, "%d", rank);
	strcpy (str1,"hwk5_");
	strcat (str1,str2);

	fp = fopen(str1, "w");
	
	for(i=1;i<(nxlim+1-istart);i++){
		fprintf(fp,"\n");
		for(j=1;j<(nylim+1-jstart);j++){			
				fprintf(fp,"%0.4f  ",u[i][j]);
				
		}
	}
	fclose(fp);
	
	delx=xincrement;
	dely=xincrement;
	
	/* after field is calculated (u) calculate the gradient */
	
	for(i=1;i<(nxlim-istart+2);i++){
		for(j=1;j<(nylim-jstart+2);j++){			
			

				
				if(i==1 && j==1){ //top left corner
					delux[i][j]=(u[i][j+1]-u[i][j])/(xincrement);
					deluy[i][j]=(u[i+1][j]-u[i][j])/(xincrement);
				}
					
		
				else if(i==(nxlim-istart) && j==(nylim-jstart)){ //lower right corner
					delux[i][j]=(u[i][j]-u[i][j-1])/(xincrement);
					deluy[i][j]=(u[i][j]-u[i-1][j])/(xincrement);
				}
				else if (i==1 && j==(nylim-jstart)){ //top right
					delux[i][j]=(u[i][j]-u[i][j-1])/xincrement;
					deluy[i][j]=(u[i+1][j]-u[i][j])/xincrement;
				}
				else if (i==(nxlim-istart) && j==1){ //bottom left
					delux[i][j]=(u[i][j+1]-u[i][j])/xincrement;
					deluy[i][j]=(u[i][j]-u[i-1][j])/xincrement;
				}
				
				else if (j==1 && i>1){//left side 
					delux[i][j]=(u[i][j+1]-u[i][j])/xincrement;
					deluy[i][j]=(u[i+1][j]-u[i-1][j])/(2*xincrement);
				}
				else if (j==nylim-jstart && i>1){//right side 
					delux[i][j]=(u[i][j]-u[i][j-1])/xincrement;
					deluy[i][j]=(u[i+1][j]-u[i-1][j])/(2*xincrement);
				}
				else if (i==1 && j>1){//top side 
					delux[i][j]=(u[i][j+1]-u[i][j-1])/(2*xincrement);
					deluy[i][j]=(u[i+1][j]-u[i][j])/xincrement;
				}
				else if (i==(nxlim-istart)&& j>1){//bottom side 
					delux[i][j]=(u[i][j+1]-u[i][j-1])/(2*xincrement);
					deluy[i][j]=(u[i+1][j]-u[i][j])/xincrement;
				}			
				
		
				
				else{
					delux[i][j]=(u[i][j+1]-u[i][j-1])/(2*xincrement);
					deluy[i][j]=(u[i+1][j]-u[i-1][j])/(2*xincrement);
				}
				
						
			}
	}
	
	
	sprintf(str2, "%d", rank);
	strcpy (str1,"deluy5_");
	strcat (str1,str2);

	fp = fopen(str1, "w");
	
	
	
	for(i=1;i<(nxlim-istart+2);i++){
		fprintf(fp,"\n");
		for(j=1;j<(nylim-jstart+2);j++){			
				fprintf(fp,"%0.4f	",deluy[i][j]);
				
		}
	}
	
	fclose(fp);
	
	
	sprintf(str2, "%d", rank);
	strcpy (str1,"delux5_");
	strcat (str1,str2);

	fp = fopen(str1, "w");
	
	
	
	for(i=1;i<(nxlim-istart+2);i++){
		fprintf(fp,"\n");
		for(j=1;j<(nylim-jstart+2);j++){			
				fprintf(fp,"%0.4f	",delux[i][j]);
				
		}
	}
	
	fclose(fp);
	
	
	
	
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
	double fi[ngrid][ngrid];
	int up_proc,down_proc,left_proc,right_proc;

	xstart = -1; xend = 1;
	ystart = -1; yend = 1;
	maxx=-1e15;
	maxy=maxx;
	delx=(xend-xstart)/(ngrid-1);
	dely=(yend-ystart)/(ngrid-1);
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
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
	
	
	



	nodecalc(i,nxlim,j,nylim,rank,rk,ck,nrow,ncol,delx,dely,xstart,xend);


	MPI_Barrier(MPI_COMM_WORLD);
	

	MPI_Finalize();

	return 0;
}






