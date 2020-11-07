#include <stdio.h>
#include <math.h>
#include "ct.h"
#include "mri.h"


#define R11 64
#define C11 64
#define ROW 64
#define COL 32


int INPUT[64][64];
int INPUT1[64][64];
int coeff1[64][64],coeff2[64][64];
int rows=64;
int columns=64;
int Input[64][64];
int Even[64][28];
int Odd[64][32];
int Low[64][32];
int High[64][32];
int LEven[32][32];
int HOdd[32][32];
int LOdd[32][32];
int HEven[32][32];
int LL[32][32];
int LH[32][32];
int HL[32][32];
int HH[32][32];
int RLL[32][32];
int RHL[32][32];
int RHH[32][32];
int RLH[32][32];
int RL[64][32];
int RH[64][32];
int R[64][32];
int H[64][32];
int Output[64][64];
int Input1[64][64];
int Input2[64][64];

int rowS=32;
int colS=32;
int rowL=64;
int colL=64;

int wavedecode[64][64];
int waveencodeImage1[64][64];
int waveencodeImage2[64][64];
int FusedImage[64][64];
int aa,cc;
void integerdwt1();
void integerdwt2();
void reversedwt();








void integerdwt1()

{//begin of integer dwt

int i,j,k,a;
int columns1;
int rows1;
rows1=rows/2;
columns1=columns/2;
	
//////// one dimensional Even component////////

			for (j=0;j<rows;j++)
			{
				a=1;
 			   for (k=0;k<columns1;k++)
			   {
					Even[j][k]=Input1[j][a];
					a=a+2;
			   }
			}


	//////// one dimensional Odd component////////


			for (j=0;j<rows;j++)
			{
				a=0;
			   for (k=0;k<columns1;k++)
			   {
					Odd[j][k]=Input1[j][a];
					a=a+2;
			   }
			}


		
	////////////  comput L  AND  H pass component
			for (j=0;j<rows;j++)
			{
			   for (k=0;k<columns1;k++)
			   {

				 High[j][k]=Odd[j][k]-Even[j][k];
				 aa=Odd[j][k]-Even[j][k];
				 aa=aa/2;
				 cc=ceil(aa);
				 Low[j][k]=(Even[j][k]+cc);
			   }
			}

	///////////////////one dimensional L even

			for (j=0;j<columns1;j++)
			{
				a=1;
			   for (k=0;k<columns1;k++)
			   {
					LEven[k][j]=Low[a][j];
					HEven[k][j]=High[a][j];
					a=a+2;
			   }
			}
////////////////////////////////////////////////////////////////


	//////// one dimensional  L Odd component////////

			for (j=0;j<columns1;j++)
			{
				a=0;
			   for (k=0;k<columns1;k++)
			   {
					LOdd[k][j]=Low[a][j];
					HOdd[k][j]=High[a][j];

					a=a+2;
			   }
			}



/////////////////////////////////////////////////////////


			for (j=0;j<columns1;j++)
			{
			   for (k=0;k<columns1;k++)
			   {

				 LH[j][k]=LOdd[j][k]-LEven[j][k];
				 aa=LOdd[j][k]-LEven[j][k];
				 aa=aa/2;
				 cc=ceil(aa);
				 LL[j][k]=(LEven[j][k]+cc);
				 HH[j][k]=HOdd[j][k]-HEven[j][k];
				 aa=HOdd[j][k]-HEven[j][k];
				 aa=aa/2;
				 cc=ceil(aa);
				 HL[j][k]=(HEven[j][k]+cc);
			   }
			}


}// end of integer dwt



void integerdwt2()

{//begin of integer dwt

int i,j,k,a;
int columns1;
int rows1;
rows1=rows/2;
columns1=columns/2;



//////// one dimensional Even component////////

			for (j=0;j<rows;j++)
			{
				a=1;
 			   for (k=0;k<columns1;k++)
			   {
					Even[j][k]=Input2[j][a];
					a=a+2;
			   }
			}


	//////// one dimensional Odd component////////


			for (j=0;j<rows;j++)
			{
				a=0;
			   for (k=0;k<columns1;k++)
			   {
					Odd[j][k]=Input2[j][a];
					a=a+2;
			   }
			}


		
	////////////  comput L  AND  H pass component
			for (j=0;j<rows;j++)
			{
			   for (k=0;k<columns1;k++)
			   {

				 High[j][k]=Odd[j][k]-Even[j][k];
				 aa=Odd[j][k]-Even[j][k];
				 aa=aa/2;
				 cc=ceil(aa);
				 Low[j][k]=(Even[j][k]+cc);
			   }
			}

	///////////////////one dimensional L even

			for (j=0;j<columns1;j++)
			{
				a=1;
			   for (k=0;k<columns1;k++)
			   {
					LEven[k][j]=Low[a][j];
					HEven[k][j]=High[a][j];
					a=a+2;
			   }
			}
////////////////////////////////////////////////////////////////


	//////// one dimensional  L Odd component////////

			for (j=0;j<columns1;j++)
			{
				a=0;
			   for (k=0;k<columns1;k++)
			   {
					LOdd[k][j]=Low[a][j];
					HOdd[k][j]=High[a][j];

					a=a+2;
			   }
			}



/////////////////////////////////////////////////////////


			for (j=0;j<columns1;j++)
			{
			   for (k=0;k<columns1;k++)
			   {

				 LH[j][k]=LOdd[j][k]-LEven[j][k];
				 aa=LOdd[j][k]-LEven[j][k];
				 aa=aa/2;
				 cc=ceil(aa);
				 LL[j][k]=(LEven[j][k]+cc);
				 HH[j][k]=HOdd[j][k]-HEven[j][k];
				 aa=HOdd[j][k]-HEven[j][k];
				 aa=aa/2;
				 cc=ceil(aa);
				 HL[j][k]=(HEven[j][k]+cc);
			   }
			}


}// end of integer dwt


void reversedwt()
{//begin of reverse dwt



int i,j,k,a;
int columns1=32;
int Lenr =32;
int Lenc =32;
int rlen2r=64;
int valc;

/////init the RL RH arrays (4x2)

for(i=0;i<ROW;i++){
	for(j=0;j<COL;j++){
	RL[i][j]=0;
	RH[i][j]=0;
	}
}

	for (j=0;j<columns1;j++)
			{
			   for (k=0;k<columns1;k++)
			   {

				 aa=LH[j][k];
				 aa=aa/2;
				 cc=ceil(aa);

				 RLL[j][k]=(LL[j][k]-cc);

				 RLH[j][k]=RLL[j][k]+LH[j][k];


				 aa=HH[j][k];
				 aa=aa/2;
				 cc=ceil(aa);

				 RHL[j][k]=(HL[j][k]-cc);

				 RHH[j][k]=HH[j][k]+RHL[j][k];





			   }
			}



					k=0;
					for(i=1;i<ROW;i=i+2){
						for(j=0;j<COL;j++){
							RL[i][j]=RLL[k][j];
							RH[i][j]=RHL[k][j];

						}
						k++;
					}


						k=0;
					for(i=0;i<ROW;i=i+2){
						for(j=0;j<COL;j++){
							RL[i][j]=RLH[k][j];
							RH[i][j]=RHH[k][j];

						}
						k++;
					}


					for(i=0;i<ROW;i=i+1)
						{
						for(j=0;j<COL;j++)
						{
							aa=RH[i][j]/2;
							cc=ceil(aa);
							R[i][j]=(RL[i][j]-cc);
							H[i][j]=R[i][j]+RH[i][j];

						}

					}




					for(i=0;i<ROW;i++){
							k=0;
						for(j=1;j<ROW;j=j+2){
							Output[i][j]=R[i][k];
							k++;

						}

					}



					for(i=0;i<ROW;i++){
							k=0;
						for(j=0;j<ROW;j=j+2){
							Output[i][j]=H[i][k];
							k++;

						}

					}


}// end of reverse dwt





void main()

{
int i,j,m;
////////////////////////////////////////////////input

				   for( i=0;i<64;i++)
					{
						for( j=0;j<64;j++)
					{
				//	printf("%d \n",REDB[i][j]);
					Input1[i][j]=inpt[i][j];
					printf("%d \n",Input1[i][j]);
					}

					}//end of outer for
//									printf("Input2\n");
//
//
//////////////////////////// input2//////////////////////////////
//
				   for( i=0;i<R11;i++)
					{
						for( j=0;j<C11;j++)
					{
				//	printf("%d \n",REDB1[i][j]);
					Input2[i][j]=inpt1[i][j];
				printf("%d \n",Input2[i][j]);
					}

					}//end of outer for
									//	printf("value\n");

//
////////////////////////////////////////////////////////////////////////////
//


					integerdwt1();
				//	printf("matrix all: \n");

					for(i=0;i<rowL;i++){
							for(j=0;j<colL;j++){

								if(i<rowS){
									if(j<colS){
										waveencodeImage1[i][j]=LL[i][j];
									}
									else{
										waveencodeImage1[i][j]=LH[i][j-colS];
									}
								}
								else{
									if(j<colS){
										waveencodeImage1[i][j]=HL[i-rowS][j];
									}
									else{
										waveencodeImage1[i][j]=HH[i-rowS][j-colS];
									}

								}
								printf(" %d\n",waveencodeImage1[i][j]);
							}
							//printf("\n");
						}


									///////////////////


					integerdwt2();
				//	printf("matrix all: \n");

					for(i=0;i<rowL;i++){
							for(j=0;j<colL;j++){

								if(i<rowS){
									if(j<colS){
										waveencodeImage2[i][j]=LL[i][j];
									}
									else{
										waveencodeImage2[i][j]=LH[i][j-colS];
									}
								}
								else{
									if(j<colS){
										waveencodeImage2[i][j]=HL[i-rowS][j];
									}
									else{
										waveencodeImage2[i][j]=HH[i-rowS][j-colS];
									}

								}
							printf(" %d\n",waveencodeImage2[i][j]);
							}
							//printf("\n");
						}





					///////////////////////////////

					for(i=0;i<R11;i++)
						{
							for (j=0;j<C11;j++)
								{

								coeff1[i][j]=waveencodeImage1[i][j];
								coeff2[i][j]=waveencodeImage2[i][j];
								//if (coeff1[i][j] > coeff2[i][j])
								
								FusedImage[i][j]= (coeff1[i][j]+coeff2[i][j])/2;
								//FusedImage[i][j]= coeff1[i][j] ;
								
								//else
								//{
								//FusedImage[i][j]=coeff2[i][j];
								//




								}
						}



                         for(i=0;i<rowL;i++){
							for(j=0;j<colL;j++){

								if(i<rowS){
									if(j<colS){
										RLL[i][j]=FusedImage[i][j];
									}
									else{
										RLH[i][j-colS]=FusedImage[i][j];
									}
								}
								else{
									if(j<colS){
										RHL[i-rowS][j]=FusedImage[i][j];
									}
									else{
										RHH[i-rowS][j-colS]=FusedImage[i][j];
									}

								}
								//printf(" %d  ",all[i][j]);
							}
							//printf("\n");
						}


reversedwt();
//printf("exit\n");
for( i=0;i<64;i++)
					{
						for( j=0;j<64;j++)
					{
							printf("%d\n",Output[i][j]);
						//wavedecode[i][j]=Output[i][j];
					}

					}
				//	printf("exit\n");

//
//
//
//
//





}//end of main

