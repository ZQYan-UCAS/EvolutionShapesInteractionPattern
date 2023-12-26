//**************************************************************************************************************************************
//==================================Genetic Algorithm Optimization of Evolution Pathway for Real Protein================================
//=============================================Coded by Zhiqiang Yan at CIAC, 5/2018====================================================
//**************************************************************************************************************************************
#include "time.h"
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#define IT 5000	 			//number of optimization step
#define COMNUM 4
#define RMSD_CUTOFF 3.0
#define SIZE 500
#define PER 0.05
#define NM 1
#define NC 0
#define PM	1.0
#define PC 	0.0
#define SEQBIN 100
#define CONFBIN 100
#define BIN 10
//#define RUNSEQ 4

static long i,run_num,run_num1,run_num2,min_conformation[IT+1],fixconf;
static long domain_min_index,domain_pop_conf[SIZE+1],domain_pre_pop_conf[SIZE+1],all_conf[100000],domain_temp_pop_conf[SIZE+1];
static long complex_min_index,complex_pop_conf[SIZE+1],complex_pre_pop_conf[SIZE+1],complex_temp_pop_conf[SIZE+1];

static int flag_domain_design,flag_complex_design,hist_survive[IT+1],hist_die[IT+1];
static int hist_new_seq[IT+1],hist_lost_seq[IT+1],hist_new_conf[IT+1],hist_lost_conf[IT+1],new_pop_num;
static int sequence1[600],sequence2[600],select_seq1,select_seq2,sequence_native[COMNUM+1][600];
static int domain_conformation[3000][400][2],complex_conformation[3000][400][2],sequence[600],domain_num_contact[3000],complex_num_contact[3000];
static int domain_pre_pop_seq[SIZE+1][600],domain_pop_seq[SIZE+1][600],pre_pop_seq[SIZE+1][600],pop_seq[SIZE+1][600],start_seq[SIZE+1][600],all_seq[10000][600],all_seq_num,all_conf_num;
static int domain_temp_pop_seq[SIZE+1][600],complex_temp_pop_seq[SIZE+1][600],temp_pop_seq[SIZE+1][600];
static int length[COMNUM+1],complex_pre_pop_seq[SIZE+1][600],complex_pop_seq[SIZE+1][600];
static int partner_one_length[COMNUM+1],partner_two_length[COMNUM+1],complex_decoy_num[COMNUM+1],domain_one_decoy_num[COMNUM+1],domain_two_decoy_num[COMNUM+1],contact_num[COMNUM+1][3000];

static char filename[200],path[200],residue[20][4],residue_large[20][4],residue_simple[20][2];
static char complex_name[COMNUM+1][20],dock_partner_one[COMNUM+1][10],dock_partner_two[COMNUM+1][10];

static double aver_e[50],square_aver[50],contactorder[3000],lrmsd[COMNUM+1][3000],irmsd[COMNUM+1][3000],fnat[COMNUM+1][3000],score[COMNUM+1][3000],gap_domain,gap_complex;
static double domain_gap,domain_average_gap,domain_roughness,domain_isr,domain_min_e,domain_delta_g;
static double complex_gap,complex_average_gap,complex_roughness,complex_isr,complex_min_e,complex_delta_g;
static double domain_pre_fitfunction[SIZE+1],domain_fitfunction[SIZE+1],domain_temp_fitfunction[SIZE+1],domain_temp_parameter[SIZE+1][6];
static double complex_pre_fitfunction[SIZE+1],complex_fitfunction[SIZE+1],complex_temp_fitfunction[SIZE+1],complex_temp_parameter[SIZE+1][6];
static double mj_matrix[20][20],domain_energy[3000],complex_energy[3000],domain_parameter[SIZE+1][6],complex_parameter[SIZE+1][6],domain_average_parameter[6],complex_average_parameter[6],beta,temperature;

FILE *fp;
FILE *fp_fit;
FILE *fp_seq;
FILE *fp_e;
FILE *fp_isr;
FILE *fp_gap;
FILE *fp_average_gap;
FILE *fp_delta_g;
FILE *fp_contactorder;
FILE *fp_rough;
FILE *fp_conf;
FILE *fp_stat_seq;
FILE *fp_stat_conf;
FILE *fp_stat_parameter;
//-------------------------------------------------Generate random number--------------------------------------------------------
#include <iostream>
 
const double InvMaxRand=2.3283064359965952029459655278022e-10;
const int NBIT_FORRAND=32;

int irdum1[16384], irdum2;

int i_rand()
{
	irdum2=((irdum2+1)&16383);
	irdum1[irdum2]=(((irdum1[((irdum2-157)&(16383))])^(irdum1[((irdum2- 314)&(16383))]))^((irdum1[((irdum2-471)&(16383))])^(irdum1[((irdum2-9689)&(16383))])));
	
	return(irdum1[irdum2]);
}

double f_rand()
{
	return((i_rand()+0.5)*InvMaxRand+0.5);
}

void initrand(unsigned int nseed)
{
	using namespace std;
 	int ibm,idum;
	int i;
	
	ibm=static_cast<int>(nseed+nseed+1);

	for(irdum2=0;irdum2<16384;irdum2++)
	{
		idum=0;
		for(i=0;i<NBIT_FORRAND;i++)
		{
			idum+=idum;
			ibm*=16807;
			if(ibm<=0) 
				idum++;
		}
		irdum1[irdum2]=idum;
	}
	
	irdum2=0;
	
	for(i=0;i<200000;i++)     //warm up the generator
		idum=i_rand();
}
//-----------------------------------Read the index of the amino acid in MJ matrix-----------------------------------------------
void ReadAAIndex()
{
	int m;

//	sprintf(filename,"%s\\RunEvolution\\AAIndex.dat",path);
	sprintf(filename,"%s/RunEvolution/AAIndex.dat",path);
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}

	for(m=0;m<20;m++)
	{
		fscanf(fp,"%s",residue[m]);
		fscanf(fp,"%s",residue_large[m]);
		fscanf(fp,"%s",residue_simple[m]);
	}

//	for(m=0;m<20;m++)
//		printf("%s %s\n",residue[m],residue_simple[m]);

	fclose(fp);
}
//-----------------------------------Read the MJ matrix-------------------------------------------------------------------------
void ReadMJPotential()
{
	int m,n;

//	sprintf(filename,"%s\\RunEvolution\\MJMatrix.dat",path);
	sprintf(filename,"%s/RunEvolution/MJMatrix.dat",path);
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("can't find the file MJMatrix.dat!\n");
		exit(0);
	}

	for(m=0;m<20;m++)
	{
		for(n=m;n<20;n++)
		{
			fscanf(fp,"%lf",&mj_matrix[m][n]);
			mj_matrix[n][m]=mj_matrix[m][n];
//			printf("%8.3lf",mj_matrix[m][n]);
		}
		fscanf(fp,"%*[^\n]");
//		printf("\n");
	}

	fclose(fp);
}
//------------------------------------------Read the information for the domain family-------------------------------------------
void ReadComplexSet()
{
	int m,n;

//	sprintf(filename,"%s\\RunEvolution\\3DGE_ComplexSet.dat",path);
	sprintf(filename,"%s/RunEvolution/3DGE_ComplexSet.dat",path);
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("can't find the file WW_Complex_Final.dat!\n");
		exit(0);
	}

	m=0,n=0;
	while(feof(fp)==0)
	{
		fscanf(fp,"%d",&m);
		fscanf(fp,"%s",complex_name[m]);
		fscanf(fp,"%s",dock_partner_one[m]);
		fscanf(fp,"%s",dock_partner_two[m]);
		fscanf(fp,"%ld",&complex_decoy_num[m]);
		fscanf(fp,"%ld",&domain_one_decoy_num[m]);
		fscanf(fp,"%ld",&domain_two_decoy_num[m]);
		fscanf(fp,"%*[^\n]");

		if(m==COMNUM)
			break;
	}
	fclose(fp);
}
//---------------------------------------Read Native Sequence for Each Domain----------------------------------------------------
void ReadNativeComplexSequence()
{
	int tempd,m,n;
	char temps[600];

		
	for(m=1;m<=COMNUM;m++)
	{	
//		sprintf(filename,"%s\\CalContact\\NativeSequence\\%s_NativeSequence.dat",path,complex_name[m]);
		sprintf(filename,"%s/CalContact/NativeSequence/%s_NativeSequence.dat",path,complex_name[m]);
		if((fp=fopen(filename,"r"))==NULL)
		{
			printf("can't find the file %s!\n",filename);
			exit(0);
		}
	
		fscanf(fp,"%d",&tempd);
		fscanf(fp,"%s",temps);

		for(n=1;n<=strlen(dock_partner_one[m]);n++)
		{
			fscanf(fp,"%s",temps);
			partner_one_length[m]+=strlen(temps);
		}
		for(n=1;n<=strlen(dock_partner_two[m]);n++)
		{
			fscanf(fp,"%s",temps);
			partner_two_length[m]=strlen(temps);
		}

		for(n=1;n<=partner_one_length[m];n++)
			fscanf(fp,"%d",&sequence_native[m][n]);

		for(n=partner_one_length[m]+1;n<=partner_one_length[m]+partner_two_length[m];n++)
			fscanf(fp,"%d",&sequence_native[m][n]);
			
		fscanf(fp,"%*[^\n]");

		length[m]=partner_one_length[m]+partner_two_length[m];

		fclose(fp);
	}

//	printf("one=%4d, two=%4d\n",partner_one_length[1],partner_two_length[1]);
//	for(n=1;n<=partner_one_length[1]+partner_two_length[1];n++)
//		printf("%3d",sequence_native[1][n]);
//	printf("\n");
}
//---------------------------------Read the number of the starting sequences-----------------------------------------------------
void ReadRunNum()
{
	int tempd;	

	sprintf(filename,"RunNum.dat");
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("can't find the file %s\n",filename);
		exit(0);
	}

	fscanf(fp,"%d",&tempd);
	run_num=tempd;
	fscanf(fp,"%lf",&gap_domain);
	fscanf(fp,"%lf",&gap_complex);
	fscanf(fp,"%d",&tempd);
	run_num1=SIZE*(tempd-1)+1,run_num2=SIZE*tempd;
	fclose(fp);

//	printf("%4d %4d gap_domain=%8.3lf, gap_complex=%8.3lf\n",run_num,run_num1,gap_domain,gap_complex);
//	run_num=21;
//	run_num1=1;
//	run_num2=500;
}
//-----------------------------------Read all the conformations------------------------------------------------------------------
void ReadDomainNativeContactMap()
{
	int m,n;

//	sprintf(filename,"%s\\ContactMap\\Native\\%s\\RR_%s_%s_ContactMap_Domain_HA_Native_5.0.dat",path,complex_name[run_num],complex_name[run_num],dock_partner_two[run_num]);
	sprintf(filename,"%s/ContactMap/Native/%s/RR_%s_%s_ContactMap_Domain_HA_Native_5.0.dat",path,complex_name[run_num],complex_name[run_num],dock_partner_two[run_num]);
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}
	
	m=0,n=0;				//here m=0 means the native is set to be 0
	while(feof(fp)==0)
	{
		n++;
			
		fscanf(fp,"%d",&domain_conformation[m][n][0]);
		fscanf(fp,"%d",&domain_conformation[m][n][1]);
			
		fscanf(fp,"%*[^\n]");
	}
	fclose(fp);

//	printf("%4d %4s is reading ...\n",run_num,complex_name[run_num]);

	domain_num_contact[m]=n-1;
//	for(n=1;n<=domain_num_contact[m];n++)
//		printf("%4d %4d\n",domain_conformation[m][n][0],domain_conformation[m][n][1]);
}
//-----------------------------------Read all the conformations------------------------------------------------------------------
void ReadDomainDecoyContactMap()
{
	int m,n;

	for(m=1;m<=domain_two_decoy_num[run_num];m++)
	{
//		sprintf(filename,"%s\\ContactMap\\Decoy\\%s\\RR_%s_ContactMap_Domain_HA_Decoy_5.0_%d.dat",path,complex_name[run_num],complex_name[run_num],m);
		sprintf(filename,"%s/ContactMap/Decoy/%s/RR_%s_ContactMap_Domain_HA_Decoy_5.0_%d.dat",path,complex_name[run_num],complex_name[run_num],m);
		if((fp=fopen(filename,"r"))==NULL)
		{
			printf("can't find the file %s!\n",filename);
			exit(0);
		}
	
		n=0;
		while(feof(fp)==0)
		{
			n++;
				
			fscanf(fp,"%d",&domain_conformation[m][n][0]);
			fscanf(fp,"%d",&domain_conformation[m][n][1]);
				
			fscanf(fp,"%*[^\n]");
		}
		fclose(fp);

//		printf("%4d %4s is reading ...%4d\n",m,complex_name[run_num],n-1);

		domain_num_contact[m]=n-1;
//		for(n=1;n<=domain_num_contact[m];n++)
//			printf("%4d %4d\n",domain_conformation[m][n][0],domain_conformation[m][n][1]);
	}
}
//-----------------------------------Read all the conformations------------------------------------------------------------------
void ReadComplexNativeContactMap()
{
	int m,n;

//	sprintf(filename,"%s\\ContactMap\\Native\\%s\\%s_ContactMap_Complex_HA_Native_5.0.dat",path,complex_name[run_num],complex_name[run_num]);
	sprintf(filename,"%s/ContactMap/Native/%s/%s_ContactMap_Complex_HA_Native_5.0.dat",path,complex_name[run_num],complex_name[run_num]);
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}
	
	m=0,n=0;
	while(feof(fp)==0)
	{
		n++;
			
		fscanf(fp,"%d",&complex_conformation[m][n][0]);
		fscanf(fp,"%d",&complex_conformation[m][n][1]);
			
		fscanf(fp,"%*[^\n]");
	}
	fclose(fp);

//	printf("%4d %4s is reading ...\n",m,pdbid[m]);

	complex_num_contact[m]=n-1;
//	for(n=1;n<=complex_num_contact[m];n++)
//		printf("%4d %4d\n",complex_conformation[m][n][0],complex_conformation[m][n][1]);
}
//-----------------------------------Read all the conformations------------------------------------------------------------------
void ReadComplexDecoyContactMap()
{
	int m,n;

	for(m=1;m<=complex_decoy_num[run_num];m++)
	{
//		sprintf(filename,"%s\\ContactMap\\Decoy\\%s\\%s_ContactMap_Complex_HA_Decoy_5.0_%d.dat",path,complex_name[run_num],complex_name[run_num],m);
		sprintf(filename,"%s/ContactMap/Decoy/%s/%s_ContactMap_Complex_HA_Decoy_5.0_%d.dat",path,complex_name[run_num],complex_name[run_num],m);
		if((fp=fopen(filename,"r"))==NULL)
		{
			printf("can't find the file %s!\n",filename);
			exit(0);
		}
	
		n=0;
		while(feof(fp)==0)
		{
			n++;
				
			fscanf(fp,"%d",&complex_conformation[m][n][0]);
			fscanf(fp,"%d",&complex_conformation[m][n][1]);
				
			fscanf(fp,"%*[^\n]");
		}
		fclose(fp);

//		printf("%4d %4s is reading ...\n",m,pdbid[m]);

		complex_num_contact[m]=n-1;
//		for(n=1;n<=num_contact[m];n++)
//			printf("%4d %4d\n",complex_conformation[m][n][0],complex_conformation[m][n][1]);
	}
}
//-------------------------------------Read the RMSD for the binding decoys------------------------------------------------------
void ReadComplexInformation()
{
	int m,tempd;

//	sprintf(filename,"%s\\CalContact\\ContactNumber\\%s_ContactNumber.dat",path,complex_name[run_num]);
	sprintf(filename,"%s/CalContact/ContactNumber/%s_ContactNumber.dat",path,complex_name[run_num]);
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}

	for(m=0;m<=complex_decoy_num[run_num];m++)
	{
		fscanf(fp,"%d",&tempd);
		fscanf(fp,"%d",&contact_num[run_num][tempd]);
		fscanf(fp,"%lf",&lrmsd[run_num][tempd]);
		fscanf(fp,"%lf",&fnat[run_num][tempd]);
		fscanf(fp,"%lf",&score[run_num][tempd]);
		fscanf(fp,"%lf",&irmsd[run_num][tempd]);
		fscanf(fp,"%*[^\n]");

//		printf("%4d %8.3lf\n",tempd,lrmsd[run_num][tempd]);
	}
	fclose(fp);
}
//-----------------------------------Calculate the contact order of each conforamtion--------------------------------------------
void CalContactOrder()
{
	long o,p,u,v;
	double dist;

	for(o=0;o<=domain_two_decoy_num[run_num];o++)
	{	
		dist=0;
		for(p=1;p<=domain_num_contact[o];p++)
		{
			u=domain_conformation[o][p][0],v=domain_conformation[o][p][1];
				dist=dist+abs(u-v);
		}
		contactorder[o]=dist/domain_num_contact[o];
//		printf("Contact Order is %5d %8.3lf\n",o,contactorder[o]);
	}
}
//---------------------------------Read the generated population-----------------------------------------------------------------
void ReadPopulation()
{
	int m,n,o,tempd;

//	sprintf(filename,"%s\\GeneratePopulation\\%s\\%s_Population_HA_%2.1lf_%2.1lf.dat",path,complex_name[run_num],complex_name[run_num],gap_domain,gap_complex);
	sprintf(filename,"%s/GeneratePopulation/%s/%s_Population_HA_%2.1lf_%2.1lf.dat",path,complex_name[run_num],complex_name[run_num],gap_domain,gap_complex);
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("can't find the file %s\n",filename);
		exit(0);
	}

	o=0;
	for(m=1;m<=100000;m++)
	{
		fscanf(fp,"%d",&tempd);

		if(tempd>=run_num1&&tempd<=run_num2)
		{
			o++;
			for(n=partner_one_length[run_num]+1;n<=length[run_num];n++)
				fscanf(fp,"%d",&start_seq[o][n]);
			fscanf(fp,"%*[^\n]");

			for(n=1;n<=partner_one_length[run_num];n++)								//here take the native peptide sequence
				start_seq[o][n]=sequence_native[run_num][n];
		}
		fscanf(fp,"%*[^\n]");

		if(tempd==run_num2)
			break;
	}
	fclose(fp);

//	for(m=1;m<=SIZE;m++)
//	{
//		for(n=1;n<=NG;n++)
//			start_seq[m][n]=temp_start_seq[SIZE*(run_num2-1)+m][n];
//	}

//	for(n=partner_one_length[run_num]+1;n<=length[run_num];n++)
//		printf("%4d ",start_seq[1][n]);
//	printf("\n");
}
//---------------------------------Generate the sequence randomly as the start point---------------------------------------------
void RandomSequence(int x)
{
	int m,n;

	for(n=partner_one_length[run_num]+1;n<=length[run_num];n++)
	{
//		m=int(20*f_rand());
//		sequence[n]=m;

		sequence[n]=start_seq[x][n];
//		sequence[n]=start_seq[RUNSEQ][n];				//Start with the same sequence/conformation in the population
	}

//	for(n=partner_one_length[run_num]+1;n<=length[run_num];n++)
//		printf("%4d ",sequence[n]);
//	printf("\n");
}
//------------------------------------------Select the designed sequence---------------------------------------------------------
int SeekDomainNative()
{
	int p,u,v,flag[3000];
	long o;
	double second_energy,min_energy;
	double average,square,squareaver,aversquare,lamada;
	double partition,prob;

	min_energy=10000,flag_domain_design=0;
	average=0,square=0,squareaver=0,aversquare=0,lamada=0;
	for(o=0;o<=domain_two_decoy_num[run_num];o++)
	{
		domain_energy[o]=0.0;
		flag[o]=0;

		for(p=1;p<=domain_num_contact[o];p++)
		{
			u=partner_one_length[run_num]+domain_conformation[o][p][0],v=partner_one_length[run_num]+domain_conformation[o][p][1];
			domain_energy[o]+=0.1*mj_matrix[sequence[u]][sequence[v]];			//0.1 is for delta_g calculation 
		}

		if(min_energy>=domain_energy[o])					//note that the min_energy is not the last conformation
			if(min_energy>domain_energy[o])
			{
				min_energy=domain_energy[o];
				domain_min_index=o;
			}
			else 
			{
				flag[domain_min_index]=1;
				flag[o]=1;
			}

		average+=domain_energy[o];
		square+=domain_energy[o]*domain_energy[o];
	}

	aversquare=square/(domain_two_decoy_num[run_num]+1);
	squareaver=(average/(domain_two_decoy_num[run_num]+1))*(average/(domain_two_decoy_num[run_num]+1));	
	lamada=(average/domain_two_decoy_num[run_num]-min_energy)/sqrt(aversquare-squareaver);

	if(flag[domain_min_index]==0)
	{
		second_energy=1000;
		for(o=0;o<=domain_two_decoy_num[run_num];o++)
			if(domain_energy[o]>min_energy&&domain_energy[o]<second_energy)
				second_energy=domain_energy[o];
	}
	else
		second_energy=min_energy;

//	if(flag[domain_min_index]==0)
//	if(flag[domain_min_index]==0&&lamada>EGAP)
	if(flag[domain_min_index]==0&&(second_energy-min_energy)>gap_domain)
		flag_domain_design=1;

	domain_min_e=min_energy;
	domain_gap=second_energy-min_energy;
	domain_average_gap=average/(domain_two_decoy_num[run_num]+1)-min_energy;
	domain_roughness=sqrt(aversquare-squareaver);
	domain_isr=lamada;

	partition=0,prob=0,domain_delta_g=0;
	if(flag[domain_min_index]==0)				//Calculate the free energy (stability, DeltaG) of the native conformation
	{
		for(o=0;o<=domain_two_decoy_num[run_num];o++)
			partition=partition+exp(-(domain_energy[o]-domain_min_e));
		prob=exp(domain_min_e-domain_min_e);
		prob=prob/partition;
		domain_delta_g=-log(prob/(1-prob));
	}
//	printf("%4d %8.3lf %8.3lf %10.3lf %10.3lf %8.3lf %8.3lf\n",domain_min_index,domain_isr,domain_delta_g,prob,partition,domain_gap,domain_average_gap);
//	printf("...\n");

	return(flag_domain_design);
}
//------------------------------------------Select the designed sequence---------------------------------------------------------
int SeekComplexNative()
{
	int p,u,v,flag[3000];
	long o;
	double second_energy,min_energy;
	double average,square,squareaver,aversquare,lamada;
	double partition,prob;

	min_energy=10000,flag_complex_design=0;
	average=0,square=0,squareaver=0,aversquare=0,lamada=0;
	for(o=0;o<=complex_decoy_num[run_num];o++)
	{
		complex_energy[o]=0;
		flag[o]=0;

		for(p=1;p<=complex_num_contact[o];p++)
		{
			u=complex_conformation[o][p][0],v=complex_conformation[o][p][1];
			complex_energy[o]+=0.1*mj_matrix[sequence[u]][sequence[v]];			//0.1 is for delta_g calculation 
		}

		if(min_energy>=complex_energy[o])					//note that the min_energy is not the last conformation
			if(min_energy>complex_energy[o])
			{
				min_energy=complex_energy[o];
				complex_min_index=o;
			}
			else 
			{
				flag[complex_min_index]=1;
				flag[o]=1;
			}

		average+=complex_energy[o];
		square+=complex_energy[o]*complex_energy[o];

//		printf("complex energy =%4d %8.3lf\n",o,complex_energy[o]);
	}

	aversquare=square/(complex_decoy_num[run_num]+1);
	squareaver=(average/(complex_decoy_num[run_num]+1))*(average/(complex_decoy_num[run_num]+1));	
	lamada=(average/(complex_decoy_num[run_num]+1)-min_energy)/sqrt(aversquare-squareaver);

	if(flag[complex_min_index]==0)
	{
		second_energy=1000;
		for(o=0;o<=complex_decoy_num[run_num];o++)
			if(complex_energy[o]>min_energy&&complex_energy[o]<second_energy)
				second_energy=complex_energy[o];
	}
	else
		second_energy=min_energy;

//	if(flag[min_index]==0)
//	if(flag[min_index]==0&&lamada>EGAP)
	if(flag[complex_min_index]==0&&(second_energy-min_energy)>gap_complex)
		flag_complex_design=1;

	complex_min_e=min_energy;
	complex_gap=second_energy-min_energy;
	complex_average_gap=average/(complex_decoy_num[run_num]+1)-min_energy;
	complex_roughness=sqrt(aversquare-squareaver);
	complex_isr=lamada;

	partition=0,prob=0,complex_delta_g=0;
	if(flag[complex_min_index]==0)				//Calculate the free energy (stability, DeltaG) of the native conformation
	{
		for(o=0;o<=complex_decoy_num[run_num];o++)
			partition=partition+exp(-(complex_energy[o]-complex_min_e));
		prob=exp(complex_min_e-complex_min_e);
		prob=prob/partition;
		complex_delta_g=-log(prob/(1-prob));
	}
//	printf("\n%4d %8.3lf %8.3lf %10.3lf %10.3lf %8.3lf %8.3lf\n",complex_min_index,complex_isr,complex_delta_g,prob,partition,complex_gap,complex_average_gap);
//	printf("...\n");

	return(flag_complex_design);
}
//-------------------------------Generate the initial population of the sequences------------------------------------------------
void GeneratePopulation()
{
	int m,n;
	long o;

	m=0;
	for(n=1;n<=partner_one_length[run_num];n++)
		sequence[n]=sequence_native[run_num][n];

	for(;;)
	{
		RandomSequence(m+1);
		if(SeekDomainNative()==1&&SeekComplexNative()==1&&domain_min_index==0&&lrmsd[run_num][complex_min_index]<=RMSD_CUTOFF)				//Here settings determines the types of how evolving in the sequence and structure space
		{
			m++;

			domain_fitfunction[m]=domain_min_e*domain_isr;
			domain_parameter[m][0]=domain_min_e;
			domain_parameter[m][1]=domain_gap;
			domain_parameter[m][2]=domain_average_gap;
			domain_parameter[m][3]=domain_roughness;
			domain_parameter[m][4]=domain_isr;
			domain_parameter[m][5]=domain_delta_g;
			domain_pop_conf[m]=domain_min_index;

			for(n=partner_one_length[run_num]+1;n<=length[run_num];n++)
				domain_pop_seq[m][n]=sequence[n];

			
			complex_fitfunction[m]=complex_min_e*complex_isr;
			complex_parameter[m][0]=complex_min_e;
			complex_parameter[m][1]=complex_gap;
			complex_parameter[m][2]=complex_average_gap;
			complex_parameter[m][3]=complex_roughness;
			complex_parameter[m][4]=complex_isr;
			complex_parameter[m][5]=complex_delta_g;
			complex_pop_conf[m]=complex_min_index;

			for(n=1;n<=partner_one_length[run_num];n++)
				complex_pop_seq[m][n]=sequence[n];

			for(n=1;n<=length[run_num];n++)
				pop_seq[m][n]=sequence[n];

//			printf("%4d is generated......\n",m);
		}

		if(m%100==0)
		{
			sprintf(filename,"Domain_EnergyDistribution_%d.dat",m);
			if((fp=fopen(filename,"w"))==NULL)
			{
				printf("can't find the file %s\n",filename);
				exit(0);
			}

			o=0;
			fprintf(fp,"%5d %8.3lf\n",o,domain_min_e);
			for(o=0;o<=domain_two_decoy_num[run_num];o++)
				fprintf(fp,"%5d %8.3lf\n",o,domain_energy[o]);
			fclose(fp);
		}

		if(m%100==0)
		{
			sprintf(filename,"Complex_EnergyDistribution_%d.dat",m);
			if((fp=fopen(filename,"w"))==NULL)
			{
				printf("can't find the file %s\n",filename);
				exit(0);
			}

			o=0;
			fprintf(fp,"%5d %8.3lf\n",o,complex_min_e);
			for(o=0;o<=complex_decoy_num[run_num];o++)
				fprintf(fp,"%5d %8.3lf\n",o,complex_energy[o]);
			fclose(fp);
		}
				
		if(m==SIZE)
			break;
	}

//	fixconf=min_index;				//set the fixed conformation
	
//	for(n=1;n<=length[run_num];n++)
//		printf("%4d ",pop_seq[SIZE][n]);
//	printf("\n");
}
//----------------------------------------------Update the population with evolved sequence--------------------------------------
void UpdatePopulation(int x)
{
	int y;

	for(x=1;x<=SIZE;x++)
	{
		domain_fitfunction[x]=domain_temp_fitfunction[x];
		domain_parameter[x][0]=domain_temp_parameter[x][0];
		domain_parameter[x][1]=domain_temp_parameter[x][1];
		domain_parameter[x][2]=domain_temp_parameter[x][2];
		domain_parameter[x][3]=domain_temp_parameter[x][3];
		domain_parameter[x][4]=domain_temp_parameter[x][4];
		domain_parameter[x][5]=domain_temp_parameter[x][5];
		domain_pop_conf[x]=domain_temp_pop_conf[x];

		for(y=partner_one_length[run_num]+1;y<=length[run_num];y++)
			domain_pop_seq[x][y]=domain_temp_pop_seq[x][y];
//			domain_pop_seq[x][y]=sequence[y];

		complex_fitfunction[x]=complex_temp_fitfunction[x];
		complex_parameter[x][0]=complex_temp_parameter[x][0];
		complex_parameter[x][1]=complex_temp_parameter[x][1];
		complex_parameter[x][2]=complex_temp_parameter[x][2];
		complex_parameter[x][3]=complex_temp_parameter[x][3];
		complex_parameter[x][4]=complex_temp_parameter[x][4];
		complex_parameter[x][5]=complex_temp_parameter[x][5];
		complex_pop_conf[x]=complex_temp_pop_conf[x];

		for(y=1;y<=partner_one_length[run_num];y++)
			complex_pop_seq[x][y]=complex_temp_pop_seq[x][y];
//			complex_pop_seq[x][y]=sequence[y];

		for(y=1;y<=length[run_num];y++)
			pop_seq[x][y]=temp_pop_seq[x][y];
//			pop_seq[x][y]=sequence[y];
	}
}
//---------------------------------Use the ISR as the Roulette wheel selection criteria------------------------------------------
void SelectSequence()
{
	int m,n,rank_isr[2][SIZE+1],rank_stability[2][SIZE+1],rank_isr_index[2][SIZE+1],rank_stability_index[2][SIZE+1],rank[SIZE+1],rank_rank[SIZE+1],rank_index[SIZE+1];
	double prob,select_fit[SIZE+1];
//-----------------Rank the Delta_G and ISR for the domain-----------------
	for(m=1;m<=SIZE;m++)
	{
		rank_stability[0][m]=1;
		for(n=1;n<=SIZE;n++)
			if(-domain_parameter[m][5]<-domain_parameter[n][5])
				rank_stability[0][m]++;
			else if(domain_parameter[m][5]==domain_parameter[n][5]&&m!=n)
			{
				if(m<n)
					rank_stability[0][m]++;
			}
		rank_stability_index[0][rank_stability[0][m]]=m;
	}

	for(m=1;m<=SIZE;m++)
	{
		rank_isr[0][m]=1;
		for(n=1;n<=SIZE;n++)
			if(domain_parameter[m][4]<domain_parameter[n][4])
				rank_isr[0][m]++;
			else if(domain_parameter[m][4]==domain_parameter[n][4]&&m!=n)
			{
				if(m<n)
					rank_isr[0][m]++;
			}
		rank_isr_index[0][rank_isr[0][m]]=m;
	}
//----------------Rank the Delta_G and ISR for the complex-----------------
	for(m=1;m<=SIZE;m++)
	{
		rank_stability[1][m]=1;
		for(n=1;n<=SIZE;n++)
			if(-complex_parameter[m][5]<-complex_parameter[n][5])
				rank_stability[1][m]++;
			else if(complex_parameter[m][5]==complex_parameter[n][5]&&m!=n)
			{
				if(m<n)
					rank_stability[1][m]++;
			}
		rank_stability_index[1][rank_stability[1][m]]=m;
	}

	for(m=1;m<=SIZE;m++)
	{
		rank_isr[1][m]=1;
		for(n=1;n<=SIZE;n++)
			if(complex_parameter[m][4]<complex_parameter[n][4])
				rank_isr[1][m]++;
			else if(complex_parameter[m][4]==complex_parameter[n][4]&&m!=n)
			{
				if(m<n)
					rank_isr[1][m]++;
			}
		rank_isr_index[1][rank_isr[1][m]]=m;
	}
//----------------------------------------------------------------------------
	for(m=1;m<=SIZE;m++)
		rank[m]=rank_isr[0][m]+rank_stability[0][m]+rank_isr[1][m]+rank_stability[1][m];

	for(m=1;m<=SIZE;m++)
	{
		rank_rank[m]=1;
		for(n=1;n<=SIZE;n++)
			if(rank[m]>rank[n])
				rank_rank[m]++;
			else if(rank[m]==rank[n]&&m!=n)
			{
				if(m<n)
					rank_rank[m]++;
			}
		rank_index[rank_rank[m]]=m;
	}

	for(m=0;m<=SIZE;m++)
		select_fit[m]=0.0;

	select_fit[1]=PER;
	for(m=2;m<=SIZE;m++)
		select_fit[m]=select_fit[m-1]*(1.0-PER);

	for(m=1;m<=SIZE;m++)
		select_fit[m]=select_fit[m-1]+select_fit[m];

	for(m=1;m<=SIZE;m++)
	{
		select_fit[m]=select_fit[m]/select_fit[SIZE];
//		printf("%4d %8.3lf\n",m,select_fit[m]);
	}

	prob=f_rand();
	for(m=1;m<=SIZE;m++)
		if(prob>=select_fit[m-1]&&prob<select_fit[m])
			select_seq1=rank_index[m];

	prob=f_rand();
	for(m=1;m<=SIZE;m++)
		if(prob>=select_fit[m-1]&&prob<select_fit[m])
			select_seq2=rank_index[m];
}
//---------------------------------Random mutation on the sequence---------------------------------------------------------------
void RandomMutation(int x)
{
	int m,n;

	for(m=1;m<=length[run_num];m++)
		sequence[m]=pop_seq[x][m];

//	m=1+int(length[run_num]*f_rand());								//Free domain and peptide sequence
//	m=domain_length[run_num]+1+int(pep_length[run_num]*f_rand());	//Just free peptide sequence
	m=partner_one_length[run_num]+1+int(partner_two_length[run_num]*f_rand());						//just free domain sequence
//	printf("%4d\n",m);
	for(;;)
	{
		n=int(20*f_rand());
//		printf("%4d.....",n); 
		if(sequence[m]!=n)
		{
			sequence[m]=n;
			break;
		}
		else
			continue;
	}	
}
//-------------------------------Crossover between two selected sequences--------------------------------------------------------
void Crossover(int m, int n)
{
	int p,q;

	for(q=1;q<=length[run_num];q++)
	{
		sequence1[q]=pop_seq[m][q];
		sequence2[q]=pop_seq[n][q];
	}

	p=1+int((length[run_num])*f_rand());
//	sequence1[p]=pop_seq[n][p];
//	sequence2[p]=pop_seq[m][p];

	for(q=1;q<=p;q++)
	{
		sequence1[q]=pop_seq[n][q];
		sequence2[q]=pop_seq[m][q];
	}
}
//---------------------------Check whether there is the same sequence for the survival or death sequence-------------------------
int CheckSequence(int m)
{
	int n,o,p,q;

	q=0;
	for(n=1;n<=SIZE;n++)
		if(n!=m)
		{
			p=0;
			for(o=1;o<=length[run_num];o++)
				if(temp_pop_seq[m][o]==pop_seq[n][o])
					p++;

			if(p==length[run_num])
				q=1;
		}
	return(q);
}
//------------------------------Evolve the sequence with mutation or crossover---------------------------------------------------
void Evolution(int m)
{
	int n,o;

	if(SeekDomainNative()==1&&SeekComplexNative()==1&&domain_min_index==0&&lrmsd[run_num][complex_min_index]<=RMSD_CUTOFF)		//Here settings determines the types of how evolving in the sequence and structure space
//	if(SeekDomainNative()==1)
//	if(RMSD_CUTOFF>1)
	{
		new_pop_num++;
//		new_pop_num=m;
		hist_survive[i]++;

		domain_temp_fitfunction[new_pop_num]=domain_min_e*domain_isr;
		domain_temp_parameter[new_pop_num][0]=domain_min_e;
		domain_temp_parameter[new_pop_num][1]=domain_gap;
		domain_temp_parameter[new_pop_num][2]=domain_average_gap;
		domain_temp_parameter[new_pop_num][3]=domain_roughness;
		domain_temp_parameter[new_pop_num][4]=domain_isr;
		domain_temp_parameter[new_pop_num][5]=domain_delta_g;
		domain_temp_pop_conf[new_pop_num]=domain_min_index;

		for(n=partner_one_length[run_num]+1;n<=length[run_num];n++)
			domain_temp_pop_seq[new_pop_num][n]=sequence[n];

		complex_temp_fitfunction[new_pop_num]=complex_min_e*complex_isr;
		complex_temp_parameter[new_pop_num][0]=complex_min_e;
		complex_temp_parameter[new_pop_num][1]=complex_gap;
		complex_temp_parameter[new_pop_num][2]=complex_average_gap;
		complex_temp_parameter[new_pop_num][3]=complex_roughness;
		complex_temp_parameter[new_pop_num][4]=complex_isr;
		complex_temp_parameter[new_pop_num][5]=complex_delta_g;
		complex_temp_pop_conf[new_pop_num]=complex_min_index;

		for(n=1;n<=partner_one_length[run_num];n++)
			complex_temp_pop_seq[new_pop_num][n]=sequence[n];

		for(n=1;n<=length[run_num];n++)
			temp_pop_seq[new_pop_num][n]=sequence[n];

//		if(CheckSequence(new_pop_num)==0)					//Check every evolution step, rather than every mutation step
//			hist_new_seq[i]++;

/*		hist_survive[i]++;
		fitfunction[m]=min_e*isr;
		parameter[m][0]=min_e;
		parameter[m][1]=gap;
		parameter[m][2]=average_gap;
		parameter[m][3]=roughness;
		parameter[m][4]=isr;
		parameter[m][5]=delta_g;
		pop_conf[m]=min_index;

		for(n=1;n<=length[run_num];n++)
			pop_seq[m][n]=sequence[n];

		if(CheckSequence(m)==0)					//Check every evolution step, rather than every mutation step
			hist_new_seq[i]++;*/
	}
	else
	{
//		if(CheckSequence(m)==0)
//			hist_lost_seq[i]++;

		hist_die[i]++;

		for(;;)
		{
			SelectSequence();
			o=select_seq1;
			if(o!=m)
				break;
		}

		new_pop_num++;
//		new_pop_num=m;

		domain_temp_fitfunction[new_pop_num]=domain_fitfunction[o];
		domain_temp_parameter[new_pop_num][0]=domain_parameter[o][0];
		domain_temp_parameter[new_pop_num][1]=domain_parameter[o][1];
		domain_temp_parameter[new_pop_num][2]=domain_parameter[o][2];
		domain_temp_parameter[new_pop_num][3]=domain_parameter[o][3];
		domain_temp_parameter[new_pop_num][4]=domain_parameter[o][4];
		domain_temp_parameter[new_pop_num][5]=domain_parameter[o][5];
		domain_temp_pop_conf[new_pop_num]=domain_pop_conf[o];

		for(n=partner_one_length[run_num]+1;n<=length[run_num];n++)
			domain_temp_pop_seq[new_pop_num][n]=pop_seq[o][n];

		complex_temp_fitfunction[new_pop_num]=complex_fitfunction[o];
		complex_temp_parameter[new_pop_num][0]=complex_parameter[o][0];
		complex_temp_parameter[new_pop_num][1]=complex_parameter[o][1];
		complex_temp_parameter[new_pop_num][2]=complex_parameter[o][2];
		complex_temp_parameter[new_pop_num][3]=complex_parameter[o][3];
		complex_temp_parameter[new_pop_num][4]=complex_parameter[o][4];
		complex_temp_parameter[new_pop_num][5]=complex_parameter[o][5];
		complex_temp_pop_conf[new_pop_num]=complex_pop_conf[o];

		for(n=1;n<=partner_one_length[run_num];n++)
			complex_temp_pop_seq[new_pop_num][n]=pop_seq[o][n];

		for(n=1;n<=length[run_num];n++)
			temp_pop_seq[new_pop_num][n]=pop_seq[o][n];

/*		for(;;)
		{
			o=1+int(SIZE*f_rand());
			if(o!=m)
				break;
		}

		for(n=1;n<=length[run_num];n++)
			pop_seq[m][n]=pop_seq[o][n];
		fitfunction[m]=fitfunction[o];
		parameter[m][0]=parameter[o][0];
		parameter[m][1]=parameter[o][1];
		parameter[m][2]=parameter[o][2];
		parameter[m][3]=parameter[o][3];
		parameter[m][4]=parameter[o][4];
		parameter[m][5]=parameter[o][5];
		pop_conf[m]=pop_conf[o];*/
	}
}
//-----------------------------Open the files for output the data----------------------------------------------------------------
void OpenFiles()
{
	sprintf(filename,"Evolution_Seq.dat");
	if((fp_seq=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s\n",filename);
		exit(0);
	}

	sprintf(filename,"Evolution_Conf.dat");
	if((fp_conf=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s\n",filename);
		exit(0);
	}

	sprintf(filename,"Evolution_Fitness.dat");
	if((fp_fit=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s\n",filename);
		exit(0);
	}

	sprintf(filename,"Evolution_E.dat");
	if((fp_e=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s\n",filename);
		exit(0);
	}

	sprintf(filename,"Evolution_Gap.dat");
	if((fp_gap=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s\n",filename);
		exit(0);
	}

	sprintf(filename,"Evolution_Average_Gap.dat");
	if((fp_average_gap=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s\n",filename);
		exit(0);
	}

	sprintf(filename,"Evolution_Rough.dat");
	if((fp_rough=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s\n",filename);
		exit(0);
	}

	sprintf(filename,"Evolution_ISR.dat");
	if((fp_isr=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s\n",filename);
		exit(0);
	}

	sprintf(filename,"Evolution_DeltaG.dat");
	if((fp_delta_g=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s\n",filename);
		exit(0);
	}

	sprintf(filename,"Evolution_ContactOrder.dat");
	if((fp_contactorder=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s\n",filename);
		exit(0);
	}

	sprintf(filename,"Evolution_Stat_Seq.dat");
	if((fp_stat_seq=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s\n",filename);
		exit(0);
	}

	sprintf(filename,"Evolution_Stat_Conf.dat");
	if((fp_stat_conf=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s\n",filename);
		exit(0);
	}

	sprintf(filename,"Evolution_Stat_Parameter.dat");
	if((fp_stat_parameter=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s\n",filename);
		exit(0);
	}
}
//-----------------------------Calculate the similarity between the sequences----------------------------------------------------
double SeqSimilarity(int m, int n)
{
	int x,y;

	y=0;
	for(x=1;x<=length[run_num];x++)
		if(pop_seq[m][x]==pop_seq[n][x])
			y++;

	return(double(y)/length[run_num]);
}
//-----------------------------Calculate the similarity between the sequences----------------------------------------------------
double DomainSeqSimilarity(int m, int n)
{
	int x,y;

	y=0;
	for(x=partner_one_length[run_num]+1;x<=length[run_num];x++)
		if(domain_pop_seq[m][x]==domain_pop_seq[n][x])
			y++;

	return(double(y)/partner_two_length[run_num]);
}
//-----------------------------Calculate the similarity between the sequences----------------------------------------------------
double PepSeqSimilarity(int m, int n)
{
	int x,y;

	y=0;
	for(x=1;x<=partner_one_length[run_num];x++)
		if(complex_pop_seq[m][x]==complex_pop_seq[n][x])
			y++;

	return(double(y)/partner_one_length[run_num]);
}/*
//-----------------------------Calculate the similarity between the sequences by considering the symmetry------------------------
double SeqSimilaritySym(int m, int n)
{
	int x,y,z;

	y=0;
	for(x=1;x<=length[run_num];x++)
		if(pop_seq[m][x]==pop_seq[n][x])
			y++;

	z=0;
	for(x=1;x<=length[run_num];x++)
		if(pop_seq[m][x]==pop_seq[n][length[run_num]-x+1])
			z++;

	if(y>=z)
		return(double(y)/length[run_num]);
	else
		return(double(z)/length[run_num]);
}*/
//---------------------------Check the structural similarity between two consecutive sequences-----------------------------------
double DomainConfSimilarity(int x, int y)
{
	long m,n,p;
	int firstconf[200][200],secondconf[200][200];
	double sum=0;

	for(m=1;m<=partner_two_length[run_num];m++)
		for(n=1;n<=partner_two_length[run_num];n++)
		{
			firstconf[m][n]=0;
			secondconf[m][n]=0;
		}

	for(p=1;p<=domain_num_contact[domain_pop_conf[x]];p++)
	{
		m=domain_conformation[domain_pop_conf[x]][p][0],n=domain_conformation[domain_pop_conf[x]][p][1];
		firstconf[m][n]=1;
		firstconf[n][m]=1;

//		printf("%4d %4d\n",m,n);
	}

	for(p=1;p<=domain_num_contact[domain_pop_conf[y]];p++)
	{
		m=domain_conformation[domain_pop_conf[y]][p][0],n=domain_conformation[domain_pop_conf[y]][p][1];
		secondconf[m][n]=1;
		secondconf[n][m]=1;
	}

	sum=0;
	for(m=1;m<=partner_two_length[run_num]-1;m++)
		for(n=m+1;n<=partner_two_length[run_num];n++)
			if(firstconf[m][n]!=0&&secondconf[m][n]!=0)
					sum++;
	if(domain_num_contact[domain_pop_conf[x]]>=domain_num_contact[domain_pop_conf[y]])
		sum=sum/domain_num_contact[domain_pop_conf[y]];
	else
		sum=sum/domain_num_contact[domain_pop_conf[x]];

	return(sum);
}
//------------------------------Output the results at each evolving step---------------------------------------------------------
void OutputStep()
{
	int m,n,o,x,y,z,sum_seq,sum_conf,flag_1[SIZE+1],flag_2[SIZE+1],seq_sim_index[SIZE+1],conf_sim_index[SIZE+1];
	double sum,sum2;
//---------------------------------------------------------------------------
	if(i%SEQBIN==0)
	{
		for(m=1;m<=SIZE;m++)
		{
			fprintf(fp_seq,"%5d%4d  ",i,m);
			for(n=1;n<=length[run_num];n++)
			{
				fprintf(fp_seq,"%s",residue_simple[pop_seq[m][n]]);
				if(n==partner_one_length[run_num])
					fprintf(fp_seq,"   ");
			}
			fprintf(fp_seq,"   ");
			for(n=1;n<=length[run_num];n++)
			{
				fprintf(fp_seq,"%3d",pop_seq[m][n]);
				if(n==partner_one_length[run_num])
					fprintf(fp_seq,"   ");
			}
			fprintf(fp_seq,"\n");
		}
	}
//----------------Below is for the complex sequence------------------------
	for(m=1;m<=SIZE;m++)
		flag_1[m]=0;
	n=0;
	for(m=1;m<=SIZE;m++)
	{
		for(x=m+1;x<=SIZE;x++)
		{
			z=0;
			for(y=1;y<=length[run_num];y++)
				if(pop_seq[m][y]==pop_seq[x][y]&&flag_1[x]==0)
					z++;

			if(z==length[run_num])
			{
				n++;
				flag_1[x]=1;
			}
		}
	}
	fprintf(fp_stat_seq,"%6d    %5d",i,SIZE-n);
//--------------------------------------------------------------------------
	sum_seq=0;
	for(m=1;m<=SIZE;m++)
		if(flag_1[m]==0)
		{
			sum_seq++;
			seq_sim_index[sum_seq]=m;
		}

	for(m=1;m<=SIZE;m++)
	{
		flag_2[m]=0;
		for(n=1;n<=SIZE;n++)
		{
			y=0;
			for(x=1;x<=length[run_num];x++)
				if(pop_seq[m][x]==pre_pop_seq[n][x]&&flag_1[m]==0)
					y++;
			
			if(y==length[run_num])
				flag_2[m]=1;
		}
	}
	o=0;
	for(m=1;m<=SIZE;m++)
		if(flag_2[m]==0&&flag_1[m]==0)
			o++;
	fprintf(fp_stat_seq,"%5d",o);
//-------------------------------------------------------------------------
	for(m=1;m<=SIZE;m++)
		flag_1[m]=0;
	n=0;
	for(m=1;m<=SIZE;m++)
	{
		for(x=m+1;x<=SIZE;x++)
		{
			z=0;
			for(y=1;y<=length[run_num];y++)
				if(pre_pop_seq[m][y]==pre_pop_seq[x][y]&&flag_1[x]==0)
					z++;

			if(z==length[run_num])
			{
				n++;
				flag_1[x]=1;
			}
		}
	}

	for(m=1;m<=SIZE;m++)
	{
		flag_2[m]=0;
		for(n=1;n<=SIZE;n++)
		{
			y=0;
			for(x=1;x<=length[run_num];x++)
				if(pre_pop_seq[m][x]==pop_seq[n][x]&&flag_1[m]==0)
					y++;
			
			if(y==length[run_num])
				flag_2[m]=1;
		}
	}
	o=0;
	for(m=1;m<=SIZE;m++)
		if(flag_2[m]==0&&flag_1[m]==0)
			o++;
	if(i==0)
		o=0;
	fprintf(fp_stat_seq,"%5d",o);
//-------------------------------------------------------------------------
	for(m=1;m<=SIZE;m++)
		for(n=1;n<=length[run_num];n++)
			pre_pop_seq[m][n]=pop_seq[m][n];			

	y=0,sum=0,sum2=0;
	for(m=1;m<SIZE;m++)
		for(n=m+1;n<=SIZE;n++)
		{
			y++;
			sum+=SeqSimilarity(m,n);
		}
	fprintf(fp_stat_seq,"%8.3lf  ",sum/y);

	y=0,sum=0,sum2=0;
	for(m=1;m<sum_seq;m++)
		for(n=m+1;n<=sum_seq;n++)
		{
			y++;
			sum+=SeqSimilarity(seq_sim_index[m],seq_sim_index[n]);
		}
	if(sum_seq==1)							//in case there is only one different sequence
		sum=1,sum2=1,y=1;
	fprintf(fp_stat_seq,"%8.3lf     ",sum/y);
//----------------Below is for the domain sequence--------------------------
	for(m=1;m<=SIZE;m++)
		flag_1[m]=0;
	n=0;
	for(m=1;m<=SIZE;m++)
	{
		for(x=m+1;x<=SIZE;x++)
		{
			z=0;
			for(y=partner_one_length[run_num]+1;y<=length[run_num];y++)
				if(domain_pop_seq[m][y]==domain_pop_seq[x][y]&&flag_1[x]==0)
					z++;

			if(z==partner_two_length[run_num])
			{
				n++;
				flag_1[x]=1;
			}
		}
	}
	fprintf(fp_stat_seq,"%5d",SIZE-n);
//--------------------------------------------------------------------------
	sum_seq=0;
	for(m=1;m<=SIZE;m++)
		if(flag_1[m]==0)
		{
			sum_seq++;
			seq_sim_index[sum_seq]=m;
		}

	for(m=1;m<=SIZE;m++)
	{
		flag_2[m]=0;
		for(n=1;n<=SIZE;n++)
		{
			y=0;
			for(x=partner_one_length[run_num]+1;x<=length[run_num];x++)
				if(domain_pop_seq[m][x]==domain_pre_pop_seq[n][x]&&flag_1[m]==0)
					y++;
			
			if(y==partner_two_length[run_num])
				flag_2[m]=1;
		}
	}
	o=0;
	for(m=1;m<=SIZE;m++)
		if(flag_2[m]==0&&flag_1[m]==0)
			o++;
	fprintf(fp_stat_seq,"%5d",o);
//-------------------------------------------------------------------------
	for(m=1;m<=SIZE;m++)
		flag_1[m]=0;
	n=0;
	for(m=1;m<=SIZE;m++)
	{
		for(x=m+1;x<=SIZE;x++)
		{
			z=0;
			for(y=partner_one_length[run_num]+1;y<=length[run_num];y++)
				if(domain_pre_pop_seq[m][y]==domain_pre_pop_seq[x][y]&&flag_1[x]==0)
					z++;

			if(z==partner_two_length[run_num])
			{
				n++;
				flag_1[x]=1;
			}
		}
	}

	for(m=1;m<=SIZE;m++)
	{
		flag_2[m]=0;
		for(n=1;n<=SIZE;n++)
		{
			y=0;
			for(x=partner_one_length[run_num]+1;x<=length[run_num];x++)
				if(domain_pre_pop_seq[m][x]==domain_pop_seq[n][x]&&flag_1[m]==0)
					y++;
			
			if(y==partner_two_length[run_num])
				flag_2[m]=1;
		}
	}
	o=0;
	for(m=1;m<=SIZE;m++)
		if(flag_2[m]==0&&flag_1[m]==0)
			o++;
	if(i==0)
		o=0;
	fprintf(fp_stat_seq,"%5d",o);
//-------------------------------------------------------------------------
	for(m=1;m<=SIZE;m++)
		for(n=partner_one_length[run_num]+1;n<=length[run_num];n++)
			domain_pre_pop_seq[m][n]=domain_pop_seq[m][n];			

	y=0,sum=0,sum2=0;
	for(m=1;m<SIZE;m++)
		for(n=m+1;n<=SIZE;n++)
		{
			y++;
			sum+=DomainSeqSimilarity(m,n);
		}
	fprintf(fp_stat_seq,"%8.3lf  ",sum/y);

	y=0,sum=0,sum2=0;
	for(m=1;m<sum_seq;m++)
		for(n=m+1;n<=sum_seq;n++)
		{
			y++;
			sum+=DomainSeqSimilarity(seq_sim_index[m],seq_sim_index[n]);
		}
	if(sum_seq==1)							//in case there is only one different sequence
		sum=1,sum2=1,y=1;
	fprintf(fp_stat_seq,"%8.3lf     ",sum/y);
//----------------Below is for the peptide sequence--------------------------
	for(m=1;m<=SIZE;m++)
		flag_1[m]=0;
	n=0;
	for(m=1;m<=SIZE;m++)
	{
		for(x=m+1;x<=SIZE;x++)
		{
			z=0;
			for(y=1;y<=partner_one_length[run_num];y++)
				if(complex_pop_seq[m][y]==complex_pop_seq[x][y]&&flag_1[x]==0)
					z++;

			if(z==partner_one_length[run_num])
			{
				n++;
				flag_1[x]=1;
			}
		}
	}
	fprintf(fp_stat_seq,"%5d",SIZE-n);
//--------------------------------------------------------------------------
	sum_seq=0;
	for(m=1;m<=SIZE;m++)
		if(flag_1[m]==0)
		{
			sum_seq++;
			seq_sim_index[sum_seq]=m;
		}

	for(m=1;m<=SIZE;m++)
	{
		flag_2[m]=0;
		for(n=1;n<=SIZE;n++)
		{
			y=0;
			for(x=1;x<=partner_one_length[run_num];x++)
				if(complex_pop_seq[m][x]==complex_pre_pop_seq[n][x]&&flag_1[m]==0)
					y++;
			
			if(y==partner_one_length[run_num])
				flag_2[m]=1;
		}
	}
	o=0;
	for(m=1;m<=SIZE;m++)
		if(flag_2[m]==0&&flag_1[m]==0)
			o++;
	fprintf(fp_stat_seq,"%5d",o);
//-------------------------------------------------------------------------
	for(m=1;m<=SIZE;m++)
		flag_1[m]=0;
	n=0;
	for(m=1;m<=SIZE;m++)
	{
		for(x=m+1;x<=SIZE;x++)
		{
			z=0;
			for(y=1;y<=partner_one_length[run_num];y++)
				if(complex_pre_pop_seq[m][y]==complex_pre_pop_seq[x][y]&&flag_1[x]==0)
					z++;

			if(z==partner_one_length[run_num])
			{
				n++;
				flag_1[x]=1;
			}
		}
	}

	for(m=1;m<=SIZE;m++)
	{
		flag_2[m]=0;
		for(n=1;n<=SIZE;n++)
		{
			y=0;
			for(x=1;x<=partner_one_length[run_num];x++)
				if(complex_pre_pop_seq[m][x]==complex_pop_seq[n][x]&&flag_1[m]==0)
					y++;
			
			if(y==partner_one_length[run_num])
				flag_2[m]=1;
		}
	}
	o=0;
	for(m=1;m<=SIZE;m++)
		if(flag_2[m]==0&&flag_1[m]==0)
			o++;
	if(i==0)
		o=0;
	fprintf(fp_stat_seq,"%5d",o);
//-------------------------------------------------------------------------
	for(m=1;m<=SIZE;m++)
		for(n=1;n<=partner_one_length[run_num];n++)
			complex_pre_pop_seq[m][n]=complex_pop_seq[m][n];			

	y=0,sum=0,sum2=0;
	for(m=1;m<SIZE;m++)
		for(n=m+1;n<=SIZE;n++)
		{
			y++;
			sum+=PepSeqSimilarity(m,n);
		}
	fprintf(fp_stat_seq,"%8.3lf  ",sum/y);

	y=0,sum=0,sum2=0;
	for(m=1;m<sum_seq;m++)
		for(n=m+1;n<=sum_seq;n++)
		{
			y++;
			sum+=PepSeqSimilarity(seq_sim_index[m],seq_sim_index[n]);
		}
	if(sum_seq==1)							//in case there is only one different sequence
		sum=1,sum2=1,y=1;
	fprintf(fp_stat_seq,"%8.3lf     ",sum/y);
//--------------------------------------------------------------------------
	if(i!=0)
		fprintf(fp_stat_seq," %5d %5d %5d %5d\n",hist_survive[i],hist_die[i],hist_new_seq[i],hist_lost_seq[i]);
	else
		fprintf(fp_stat_seq," %5d %5d %5d %5d\n",hist_survive[i],hist_die[i],hist_new_seq[i],hist_lost_seq[i]);
//-----------------Below is for the conformation of the domain-------------
	if(i%CONFBIN==0)
	{
		fprintf(fp_conf,"%5d	",i);
		for(x=1;x<=SIZE;x++)
		{
			fprintf(fp_conf,"%4d ",domain_pop_conf[x]);
//			if(x%500==0)
//				fprintf(fp_conf,"\n");
		}
		fprintf(fp_conf,"\n");

		fprintf(fp_conf,"%5d	",i);
		for(x=1;x<=SIZE;x++)
		{
			fprintf(fp_conf,"%4d ",complex_pop_conf[x]);
//			if(x%500==0)
//				fprintf(fp_conf,"\n");
		}
		fprintf(fp_conf,"\n");
	}
//-------------------------------------------------------------------------
	for(m=1;m<=SIZE;m++)
		flag_1[m]=0;
	n=0;

	for(x=1;x<=SIZE;x++)
	{
		for(m=x+1;m<=SIZE;m++)
			if(domain_pop_conf[x]==domain_pop_conf[m]&&flag_1[m]==0)
			{
				n++;
				flag_1[m]=1;
			}
	}
	fprintf(fp_stat_conf,"%6d    %5d",i,SIZE-n);

	sum_conf=0;
	for(m=1;m<=SIZE;m++)
		if(flag_1[m]==0)
		{
			sum_conf++;
			conf_sim_index[sum_conf]=m;
		}

	y=0;
	for(m=1;m<=SIZE;m++)
	{
		flag_2[m]=0;
		for(n=1;n<=SIZE;n++)
			if(domain_pop_conf[m]==domain_pre_pop_conf[n]&&flag_1[m]==0)
			{
					y++;
					flag_2[m]=1;
			}
	}
	o=0;
	for(m=1;m<=SIZE;m++)
		if(flag_2[m]==0&&flag_1[m]==0)
			o++;
	fprintf(fp_stat_conf,"%5d",o);
//------------------------------------------------------------------------
	for(m=1;m<=SIZE;m++)
		flag_1[m]=0;
	n=0;
	for(x=1;x<=SIZE;x++)
	{
		for(m=x+1;m<=SIZE;m++)
			if(domain_pre_pop_conf[x]==domain_pre_pop_conf[m]&&flag_1[m]==0)
			{
				n++;
				flag_1[m]=1;
			}
	}

	y=0;
	for(m=1;m<=SIZE;m++)
	{
		flag_2[m]=0;
		for(n=1;n<=SIZE;n++)
			if(domain_pre_pop_conf[m]==domain_pop_conf[n]&&flag_1[m]==0)
			{
					y++;
					flag_2[m]=1;
			}
	}
	o=0;
	for(m=1;m<=SIZE;m++)
		if(flag_2[m]==0&&flag_1[m]==0)
			o++;
	if(i==0)
		o=0;
	fprintf(fp_stat_conf,"%5d",o);
//------------------------------------------------------------------------
	for(m=1;m<=SIZE;m++)
		domain_pre_pop_conf[m]=domain_pop_conf[m];

	y=0,sum=0;
	for(m=1;m<SIZE;m++)
		for(n=m+1;n<=SIZE;n++)
		{
			y++;
			sum+=DomainConfSimilarity(m,n);
		}
	fprintf(fp_stat_conf,"%8.3lf ",sum/y);

	y=0,sum=0;
	for(m=1;m<sum_conf;m++)
		for(n=m+1;n<=sum_conf;n++)
		{
			y++;
			sum+=DomainConfSimilarity(conf_sim_index[m],conf_sim_index[n]);
		}
	if(sum_conf==1)							//in case there is only one different confomation
		sum=1,y=1;
	fprintf(fp_stat_conf,"%8.3lf     ",sum/y);
//---------------Below is the conformation for the complex-----------------
	for(m=1;m<=SIZE;m++)
		flag_1[m]=0;
	n=0;

	for(x=1;x<=SIZE;x++)
	{
		for(m=x+1;m<=SIZE;m++)
			if(complex_pop_conf[x]==complex_pop_conf[m]&&flag_1[m]==0)
			{
				n++;
				flag_1[m]=1;
			}
	}
	fprintf(fp_stat_conf,"     %6d    %5d",i,SIZE-n);

	sum_conf=0;
	for(m=1;m<=SIZE;m++)
		if(flag_1[m]==0)
		{
			sum_conf++;
			conf_sim_index[sum_conf]=m;
		}

	y=0;
	for(m=1;m<=SIZE;m++)
	{
		flag_2[m]=0;
		for(n=1;n<=SIZE;n++)
			if(complex_pop_conf[m]==complex_pre_pop_conf[n]&&flag_1[m]==0)
			{
					y++;
					flag_2[m]=1;
			}
	}
	o=0;
	for(m=1;m<=SIZE;m++)
		if(flag_2[m]==0&&flag_1[m]==0)
			o++;
	fprintf(fp_stat_conf,"%5d",o);
//------------------------------------------------------------------------
	for(m=1;m<=SIZE;m++)
		flag_1[m]=0;
	n=0;
	for(x=1;x<=SIZE;x++)
	{
		for(m=x+1;m<=SIZE;m++)
			if(complex_pre_pop_conf[x]==complex_pre_pop_conf[m]&&flag_1[m]==0)
			{
				n++;
				flag_1[m]=1;
			}
	}

	y=0;
	for(m=1;m<=SIZE;m++)
	{
		flag_2[m]=0;
		for(n=1;n<=SIZE;n++)
			if(complex_pre_pop_conf[m]==complex_pop_conf[n]&&flag_1[m]==0)
			{
					y++;
					flag_2[m]=1;
			}
	}
	o=0;
	for(m=1;m<=SIZE;m++)
		if(flag_2[m]==0&&flag_1[m]==0)
			o++;
	if(i==0)
		o=0;
	fprintf(fp_stat_conf,"%5d\n",o);
//------------------Below is for the parameters----------------------------
	sum=0;
	fprintf(fp_contactorder,"%5d",i);
	for(x=1;x<=SIZE;x++)
	{
		sum=sum+contactorder[domain_pop_conf[x]];
		fprintf(fp_contactorder,"%8.2lf",contactorder[domain_pop_conf[x]]);
	}
	fprintf(fp_contactorder,"	%10.3lf\n",sum/SIZE);
	fprintf(fp_stat_parameter,"%6d %10.3lf",i,sum/SIZE);

//-------------------------------------------------------------------------
	sum=0;
	fprintf(fp_e,"%5d",i);
	for(x=1;x<=SIZE;x++)
	{
		sum+=domain_parameter[x][0];
		fprintf(fp_e,"%8.2lf",domain_parameter[x][0]);
	}
	fprintf(fp_e,"	%10.3lf\n",sum/SIZE);
	fprintf(fp_stat_parameter," %10.3lf",sum/SIZE);
	domain_average_parameter[0]=sum/SIZE;

	sum=0;
	fprintf(fp_gap,"%5d",i);
	for(x=1;x<=SIZE;x++)
	{
		sum+=domain_parameter[x][1];
		fprintf(fp_gap,"%8.2lf",domain_parameter[x][1]);
	}
	fprintf(fp_gap,"	%10.3lf\n",sum/SIZE);
	fprintf(fp_stat_parameter," %10.3lf",sum/SIZE);
	domain_average_parameter[1]=sum/SIZE;

	sum=0;
	fprintf(fp_average_gap,"%5d",i);
	for(x=1;x<=SIZE;x++)
	{
		sum+=domain_parameter[x][2];
		fprintf(fp_average_gap,"%8.2lf",domain_parameter[x][2]);
	}
	fprintf(fp_average_gap,"	%10.3lf\n",sum/SIZE);
	fprintf(fp_stat_parameter," %10.3lf",sum/SIZE);
	domain_average_parameter[2]=sum/SIZE;

	sum=0;
	fprintf(fp_rough,"%5d",i);
	for(x=1;x<=SIZE;x++)
	{
		sum+=domain_parameter[x][3];
		fprintf(fp_rough,"%8.2lf",domain_parameter[x][3]);
	}
	fprintf(fp_rough,"	%10.3lf\n",sum/SIZE);
	fprintf(fp_stat_parameter," %10.3lf",sum/SIZE);
	domain_average_parameter[3]=sum/SIZE;

	sum=0;
	fprintf(fp_isr,"%5d",i);
	for(x=1;x<=SIZE;x++)
	{
		sum+=domain_parameter[x][4];
		fprintf(fp_isr,"%8.2lf",domain_parameter[x][4]);
	}
	fprintf(fp_isr,"	%10.3lf\n",sum/SIZE);
	fprintf(fp_stat_parameter," %10.3lf",sum/SIZE);
	domain_average_parameter[4]=sum/SIZE;

	sum=0;
	fprintf(fp_delta_g,"%5d",i);
	for(x=1;x<=SIZE;x++)
	{
		sum+=domain_parameter[x][5];
		fprintf(fp_delta_g,"%8.2lf",domain_parameter[x][5]);
	}
	fprintf(fp_delta_g,"	%10.3lf\n",sum/SIZE);
	fprintf(fp_stat_parameter," %10.3lf",sum/SIZE);
	domain_average_parameter[5]=sum/SIZE;

	sum=0;
	fprintf(fp_fit,"%5d",i);
	for(x=1;x<=SIZE;x++)
	{
		sum+=domain_fitfunction[x];
		fprintf(fp_fit,"%8.2lf",domain_fitfunction[x]);
	}
	fprintf(fp_fit,"	%10.3lf\n",sum/SIZE);
	fprintf(fp_stat_parameter,"%12.3lf",sum/SIZE);

//-------------------------------------------------------------------------
	sum=0;
	fprintf(fp_e,"%5d",i);
	for(x=1;x<=SIZE;x++)
	{
		sum+=complex_parameter[x][0];
		fprintf(fp_e,"%8.2lf",complex_parameter[x][0]);
	}
	fprintf(fp_e,"	%10.3lf\n",sum/SIZE);
	fprintf(fp_stat_parameter," %10.3lf",sum/SIZE);
	complex_average_parameter[0]=sum/SIZE;

	sum=0;
	fprintf(fp_gap,"%5d",i);
	for(x=1;x<=SIZE;x++)
	{
		sum+=complex_parameter[x][1];
		fprintf(fp_gap,"%8.2lf",complex_parameter[x][1]);
	}
	fprintf(fp_gap,"	%10.3lf\n",sum/SIZE);
	fprintf(fp_stat_parameter," %10.3lf",sum/SIZE);
	complex_average_parameter[1]=sum/SIZE;

	sum=0;
	fprintf(fp_average_gap,"%5d",i);
	for(x=1;x<=SIZE;x++)
	{
		sum+=complex_parameter[x][2];
		fprintf(fp_average_gap,"%8.2lf",complex_parameter[x][2]);
	}
	fprintf(fp_average_gap,"	%10.3lf\n",sum/SIZE);
	fprintf(fp_stat_parameter," %10.3lf",sum/SIZE);
	complex_average_parameter[2]=sum/SIZE;

	sum=0;
	fprintf(fp_rough,"%5d",i);
	for(x=1;x<=SIZE;x++)
	{
		sum+=complex_parameter[x][3];
		fprintf(fp_rough,"%8.2lf",complex_parameter[x][3]);
	}
	fprintf(fp_rough,"	%10.3lf\n",sum/SIZE);
	fprintf(fp_stat_parameter," %10.3lf",sum/SIZE);
	complex_average_parameter[3]=sum/SIZE;

	sum=0;
	fprintf(fp_isr,"%5d",i);
	for(x=1;x<=SIZE;x++)
	{
		sum+=complex_parameter[x][4];
		fprintf(fp_isr,"%8.2lf",complex_parameter[x][4]);
	}
	fprintf(fp_isr,"	%10.3lf\n",sum/SIZE);
	fprintf(fp_stat_parameter," %10.3lf",sum/SIZE);
	complex_average_parameter[4]=sum/SIZE;

	sum=0;
	fprintf(fp_delta_g,"%5d",i);
	for(x=1;x<=SIZE;x++)
	{
		sum+=complex_parameter[x][5];
		fprintf(fp_delta_g,"%8.2lf",complex_parameter[x][5]);
	}
	fprintf(fp_delta_g,"	%10.3lf\n",sum/SIZE);
	fprintf(fp_stat_parameter," %10.3lf",sum/SIZE);
	complex_average_parameter[5]=sum/SIZE;

	sum=0;
	fprintf(fp_fit,"%5d",i);
	for(x=1;x<=SIZE;x++)
	{
		sum+=complex_fitfunction[x];
		fprintf(fp_fit,"%8.2lf",complex_fitfunction[x]);
	}
	fprintf(fp_fit,"	%10.3lf\n",sum/SIZE);
	fprintf(fp_stat_parameter,"%12.3lf\n",sum/SIZE);
}
//----------------------------Main loop for the optimization of the simulated annealing for the evolotion path-------------------
int main()
{
	int j,k,l,x,y,z,num_mutate,new_pop[SIZE+1];
	double prob;

//	strcpy(path,"C:\\Work\\Program\\Evolution\\Coevolution");
	strcpy(path,"/public/home/Skleac_wj_02/Work/Evolution/Coevolution");

	initrand((unsigned)time(NULL));

	OpenFiles();
	ReadAAIndex();
	ReadMJPotential();
	ReadComplexSet();
	ReadNativeComplexSequence();
	ReadRunNum();

	ReadDomainNativeContactMap();
	ReadDomainDecoyContactMap();
	ReadComplexNativeContactMap();
	ReadComplexDecoyContactMap();
	ReadComplexInformation();

	CalContactOrder();
	ReadPopulation();						//Read the already generated population in order to evolve from the same population
	GeneratePopulation();
	OutputStep();

	for(i=1;i<=IT;i++)						//Test how long IT is at first, see when it converges
	{
		printf("%5d...\n",i);

		new_pop_num=0;
		for(x=1;x<=SIZE;x++)
			new_pop[x]=0;

//		for(x=1;x<=SIZE/2;x++)
//			printf("%4d %4d\n",x,new_pop[x]);

		num_mutate=0;
		for(;;)
		{
				
//			prob=f_rand();
//			if(prob<PM)
//			{		
//				j=1+int(SIZE*f_rand());
				num_mutate++;
				j=num_mutate;
				RandomMutation(j);
				Evolution(j);

				if(new_pop_num==SIZE)
//				if(num_mutate==SIZE)
				{
					UpdatePopulation(j);
					break;
				}
//			}
/*			else
			{				
				for(l=1;l<=length[run_num];l++)
					sequence[l]=pop_seq[j][l];

				Evolution(j);

				if(new_pop_num==SIZE)
				{
					UpdatePopulation();
					break;
				}
			}*/

//			printf("Seq %4d ....\n",y);

		}
	
		if(i%BIN==0)
			OutputStep();

//		printf("%5d %5d \n",hist_new_seq[i],hist_lost_seq[i]);
	}	

	fclose(fp_fit);
	fclose(fp_conf);
	fclose(fp_e);
	fclose(fp_gap);
	fclose(fp_average_gap);
	fclose(fp_delta_g);
	fclose(fp_contactorder);
	fclose(fp_isr);	
	fclose(fp_seq);
	fclose(fp_stat_seq);
	fclose(fp_stat_conf);
	fclose(fp_stat_parameter);
}
