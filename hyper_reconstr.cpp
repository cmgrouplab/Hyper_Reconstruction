#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

using namespace std;

#define SIZE 300 //system size
#define HS 150    //half size
#define N SIZE*SIZE   //total number of pixel in system
#define Range 250000000

int* x;		//store the coordinate of black pixel. Store the given target first for target S2 and then store for the evolution structure, 
int* y;		//so there is not any black pixel information stored of given structure when the evolution is going on

int BP;		//number of black pixel, declaration here but will be defined when read the given configuration. 
//int cluster[BP]; 
int sn[HS];
int bn[HS];
double s2[HS];
double s2T[HS];		//temperary s2 array to store the temp s2 in the evolution
int config[SIZE][SIZE] = {};	//store the phase information so that we could exchange black and white pixel every time.
double** distance_matrix;	//declaration here and define when read the given config. Used for store the distance segment for each pixel, convinient when change the black white pixel

double obj[HS] ={};	//target S2


int xind;	//the original coordinate for the selected black pixel
int yind;
int indexn;	//the index of the selected black pixel

double* energyE1;
double* energyE2;



double alpha = 0.95;
double beta = 0.85;

int TN = 450;
double T = 0.00001;

int Nevl = 50000;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////type 0 calculate the base SN distance and the type 1 calculate the black pixed distance//////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double pixeldistance(int ini, int inj, int type)
{
	double dx;
	double dy;
	double d;
	if(type == 0)
		{
		dx = fabs(ini%SIZE - inj%SIZE);   // In order to save the memory, dynamically calculate the distance instead of any global array, same complexity for base distance.
		if(dx>=HS)	dx = SIZE - dx;

		dy = fabs(ini/SIZE - inj/SIZE);
		if(dy>=HS)	dy = SIZE - dy;
		}
	else
		{
		dx = fabs(x[ini] - x[inj]);
		if(dx>=HS)	dx = SIZE - dx;
		
		dy = fabs(y[ini] - y[inj]);
		if(dy>=HS)	dy = SIZE - dy;
		}

	
	d = sqrt(dx*dx + dy*dy);
	return d;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////calculate the Sn and Bn so that we can calculate the S2//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void init_SN_BN ()
{
	for(int i=0; i<HS ; i++)
	{
		sn[i] = 0;
		bn[i] = 0;
	}

	double dis;
	int nd;
	double dd;

	printf("Initializing the total distance segement...\n");

	for(int i=0; i<N; i++)
	{
		dis = pixeldistance(i,0,0);
		
		nd = (int)floor(dis);
		dd = dis - nd;
		
		if(dd<0.5)
		{
			if(nd<HS)	sn[nd] += 1;
		}
		else
		{
			if(nd+1<HS)	sn[nd+1] += 1;
		}
	}

	for(int i=0; i<HS; i++)
	{
		sn[i] *= N;
	}
	
	sn[0] *= 2;

	printf("Initializing the pixel distance segement...\n");
	
	for(int i=0; i<BP; i++)
	for(int j=i; j<BP; j++)
	{	
		distance_matrix[i][j] = distance_matrix[j][i] = pixeldistance(i,j,1);
		
		nd = (int)floor(distance_matrix[i][j]);
		dd = distance_matrix[i][j] - nd;

		if(dd<0.5)
		{
			if(nd<HS)	bn[nd] += 2;
		}
		else
		{
			if(nd+1<HS)	bn[nd+1] += 2;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void get_config ()
{
	FILE* fp;
	if ((fp = fopen("Digital.2D.N344.Rho0.3.dat","r"))==NULL)
	{
	printf("Can not open input file!\n");
	exit(1);
	}

	int temp;

	for (int i=0; i<N; i++)
	{
		fscanf(fp,"%d",&temp);

		if(temp != 0)
		BP++;
	}

	fclose(fp);

	x = new int[BP];
	y = new int[BP];
	distance_matrix = new double*[BP];
	for(int i=0; i<BP; i++)
	distance_matrix[i] =new double[BP];
		

	if ((fp = fopen("Digital.2D.N344.Rho0.3.dat","r"))==NULL)
	{
	printf("Can not open input file!\n");
	exit(1);
	}

	
	int counter = -1;
	for(int i=0; i<BP; i++)
	{
		do
		{
		fscanf(fp,"%d",&temp);
		counter++;
		}while(temp == 0);
		

		y[i] = counter/SIZE;
		x[i] = counter%SIZE;

	}
	

	init_SN_BN ();
	printf("Compelete the statistical calculation for given configuration ...\n\n");
}												

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int touch(int index)
{
	if(config[x[index]][y[index]] == 1) return 1;
	else return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void init_config()
{
	for(int i=0; i<SIZE; i++)
	for(int j=0; j<SIZE; j++)
	config[i][j] = 0;

	for(int i=0; i<BP; i++)
	{
		do
		{
		x[i] = rand()%SIZE; 
		y[i] = rand()%SIZE;
		}while(touch(i)==1);

		config[x[i]][y[i]] = 1;
	}

	init_SN_BN();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void sample_S2(double S[HS])
{
	for(int i=0; i<HS; i++)
	{
		S[i] = (double)bn[i]/(double)sn[i];
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int movement()
{
	int tempid = rand()%(int)BP;
	
	int tempx = x[tempid];
	int tempy = y[tempid];

	int dx = rand()%(int)SIZE;
	int dy = rand()%(int)SIZE;

	x[tempid] += dx;
	if(x[tempid]>=SIZE)	x[tempid] -= SIZE;
	else if(x[tempid]<0)	{printf("1");cin.get();}	 

	y[tempid] += dy;
	if(y[tempid]>=SIZE)	y[tempid] -= SIZE;
	else if(y[tempid]<0)	{printf("2");cin.get();}	 

	if(touch(tempid) == 1)
	{
		x[tempid] = tempx;
		y[tempid] = tempy;
		
		return 0;
	}
	else
	{
		indexn = tempid;
		xind = tempx;
		yind = tempy;

		config[x[tempid]][y[tempid]] = 1;
		config[tempx][tempy] = 0;

		return 1;
	}

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void change_config()
{
	int s;
	do
	{
	s = movement();
	}while(s == 0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void resume_config()
{
	config[x[indexn]][y[indexn]] = 0;
	config[xind][yind] = 1;

	x[indexn] = xind;
	y[indexn] = yind;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double energy(double S[HS],double *e1, double *e2)
{
	double E1=0;
	double E2=0;

	for(int i=1; i<50; i++)
	E1 += (S[i]-obj[i])*(S[i]-obj[i]);

	*e1 = E1;

	for(int i=1; i<HS; i++)
	E2 += (S[i]-obj[0]*obj[0])*i;

	E2 = E2*E2;

	*e2 = E2;
	
	return E1+E2;
	
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double d_energy(double S2[HS], double ST2[HS])
{
	double d_E=0;
	d_E = energy(ST2,energyE1,energyE2) - energy(S2,energyE1,energyE2);
	return d_E;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double PE(double dE, double T)
{
	if(dE > 0) return exp(-dE/T);
	else return 1;
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main()
{

clock_t start, end;
struct timespec time1;
struct timespec time2;
clock_gettime(CLOCK_REALTIME,&time1);
start = clock();

srand(time(NULL));

double energyt = 0;
double energyb = Range;
double energyb1 = Range;
double energyb2 = Range;
double E1 = 5;
double E2 = 5;
energyE1 = &E1;
energyE2 = &E2;

get_config();
sample_S2(obj);

ofstream fout1("object_s2.dat");
for(int i=0; i<HS; i++)
fout1<<i<<"\t"<<obj[i]<<endl;

FILE* fp = fopen("Energy.txt","w");
fclose(fp);

for(int i=0; i<HS; i++)
{
	s2[i] = 0;
	s2T[i] = 0;
}

init_config();
sample_S2(s2);
ofstream fout2("ini_s2.dat");
for(int i=0; i<HS; i++)
fout2<<i<<"\t"<<s2[i]<<endl;


double distance_matrix_temp[BP];
for(int i =0; i<BP; i++)
{
	distance_matrix_temp[i] = 0;
}

int BNT[HS];
for(int i=0; i<HS; i++)
{	
	BNT[i] = 0;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
for(int q=0; q<TN; q++)
{
	T = alpha*T;

	int Nacc = 0;

	for(int i=0; i<Nevl; i++)
	{
		for(int r=0; r<HS; r++)
		{
			BNT[r] = 0;
		}

		change_config();
		
		for(int k=0; k<BP; k++)
		{
			distance_matrix_temp[k] = distance_matrix[indexn][k];
		}

		for(int k=0; k<BP; k++)
		{
			distance_matrix[indexn][k] = distance_matrix[k][indexn] = pixeldistance(indexn,k,1);
		}

		for(int k=0; k<BP; k++)
		{
			int dold = (int)floor(distance_matrix_temp[k]);
			double dd = distance_matrix_temp[k] - dold;

			if(dd<0.5)
			{
				if(dold<HS)	BNT[dold] -= 2;
			}
			else
			{
				if(dold+1<HS)	BNT[dold+1] -= 2;
			}
		}


		for(int k=0; k<BP; k++)
		{
			int dnew = (int)floor(distance_matrix[indexn][k]);
			double dd = distance_matrix[indexn][k] - dnew;
			if(dd<0.5)
			{
				if(dnew<HS)	BNT[dnew] += 2;
			}
			else
			{
				if(dnew+1<HS)	BNT[dnew+1] += 2;
			}
		}	

		for(int r=0; r<HS; r++)
		{
			bn[r] += BNT[r];
		}



		sample_S2(s2T);

		double P = (double)(rand()%Range)/(double)Range;

		if( P > PE(d_energy(s2,s2T),T))
		{
			resume_config();

			for(int r=0; r<HS; r++)
			{	
				bn[r] = bn[r] - BNT[r];
			}
			
			for(int k=0; k<BP; k++)
			{
				distance_matrix[indexn][k] = distance_matrix_temp[k];
				distance_matrix[k][indexn] = distance_matrix_temp[k];
			}
		}
		else
		{
			for(int r=0; r<HS; r++)
			{
				s2[r] = s2T[r];
			}
			Nacc++;
		}


		energyt = energy(s2,energyE1,energyE2);

		if(energyt < energyb)
		{
		energyb = energyt;
		energyb1 = *energyE1;
		energyb2 = *energyE2;
		}
	}



		printf("%d th change of temperature has finished... \n",q+1 );
   		cout<<"The acceptance rate: "<<(double)Nacc/(double)Nevl<<endl;
   		cout<<"The energy E = "<<energyb<<endl;


   		fp = fopen("Energy.txt","a");
   		fprintf(fp, "%1.12f\t %1.12f\t %1.12f \n", energyb,energyb1,energyb2);
   		fclose(fp);											


		printf("*************************************************\n");

		

}


//fp = fopen("Fconfig.txt","w");
//for(int it=0; it<BP; it++)
//  {
//    fprintf(fp, "%d \t %d \n", x[it], y[it]);
//  }
//fclose(fp);					


fp = fopen("S2.txt", "w");
for(int r=0; r<HS; r++)
  {
    fprintf(fp, "%d \t %f \n", r, s2[r]);
  }
fclose(fp);


fp = fopen("Fconfig.vtk","w");
for(int it=0; it<SIZE; it++)
for(int jt=0; jt<SIZE; jt++)
  {
    fprintf(fp, "%d\n", config[it][jt]);
  }
fclose(fp);					


clock_gettime(CLOCK_REALTIME,&time2);
end = clock();

cout<<"Real run time:"<<(time2.tv_sec-time1.tv_sec)<<"s"<<endl;
cout<<"CPU run time:"<<(double)(end - start)/CLOCKS_PER_SEC<<"s"<<endl;




return 0; 
}
