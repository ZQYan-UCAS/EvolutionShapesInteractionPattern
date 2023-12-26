//**************************************************************************************************************************************
//==================================    Statistical Coupling Analysis of Aligned Seqeuence==============================================
//=============================================Coded by Zhiqiang Yan at CIAC, 01/2020===================================================
//**************************************************************************************************************************************
#include "time.h"
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"


#define COMNUM 1
#define NA 20
#define POS_CUT 0.5		//The value to remove the position, fraction of gaps (. or-) in the postion
#define SEQ_CUT 0.2			//The value to remove the sequence, fraction of gaps in the sequence

static int i,j,flag_seq[20000],flag_aa[200],group_type[NA+1][3],group_number[3][NA+1],type_number[3],couple_position[15000][2],rank_couple[15000],rank_couple_index[15000];
static int complex_length[COMNUM+1],protein_one_length[COMNUM+1],protein_two_length[COMNUM+1],seqnum[COMNUM+1];
static int hist_position_aa[200][NA+1],hist_pair_aa[200][200][NA+1][NA+1],sequence_length[COMNUM+1],sequence_num[COMNUM+1];
static char filename[120],residue[NA+1][4],residue_simple[NA+1][2],sequence_name[10000][300],sequence_id[10000][200],sequence_origin[10000][200];
static char complex_name[COMNUM+1][20],chainid[COMNUM+1][2];
static double freq_position_aa[200][NA+1],sum_freq_position[200],sum_freq_pair[200][200],position_conservation[200],freq_pair_aa[200][200][NA+1][NA+1];
static double freq_position_group[3][200][NA+1],freq_pair_group[200][200][NA+1][NA+1];
static double freq_coupling[200][200][NA+1][NA+1],weight_freq_coupling[200][200][NA+1][NA+1],sum_freq_coupling[200][200];
static double mean_freq[NA+1],group_freq[3][NA+1],phi[200][NA+1],couple_value[15000];

FILE *fp;
FILE *fp_out;
/*
//-----------------------------------Read the mean fequency of 20 amino acids----------------------------------------------------
void ReadAAMeanFrequency()
{
	int m;
	char temps[30];

	if((fp=fopen("AAMeanFrequency.dat","r"))==NULL)
	{
		printf("can't find the file AAIndex.dat!\n");
		exit(0);
	}

	for(m=0;m<20;m++)
	{
		fscanf(fp,"%s",residue[m]);
		fscanf(fp,"%s",temps);
		fscanf(fp,"%s",residue_simple[m]);
		fscanf(fp,"%lf",&mean_freq[m]);
		fscanf(fp,"%*[^\n]");
	}
	fclose(fp);

	for(m=0;m<20;m++)
		printf("%4d %s %s %8.3lf\n",m,residue[m],residue_simple[m],mean_freq[m]);
	
}*/
//-----------------------------------Read the mean fequency of 20 amino acids----------------------------------------------------
void ReadAAGroup()
{
	int m,n;
	char temps[30];

	if((fp=fopen("AAGroup.dat","r"))==NULL)
	{
		printf("can't find the file AAIndex.dat!\n");
		exit(0);
	}

	for(m=0;m<20;m++)
	{
		fscanf(fp,"%s",residue[m]);
		fscanf(fp,"%s",temps);
		fscanf(fp,"%s",residue_simple[m]);
		fscanf(fp,"%lf",&mean_freq[m]);

		for(n=0;n<=2;n++)
			fscanf(fp,"%d",&group_type[m][n]);

		fscanf(fp,"%*[^\n]");
	}
	fclose(fp);

	for(n=0;n<=2;n++)
	{
		for(m=0;m<20;m++)
		{
			group_number[n][group_type[m][n]]++;
			group_freq[n][group_type[m][n]]+=mean_freq[m];
		}
	}

	for(n=0;n<=2;n++)
		for(m=0;m<20;m++)
			if(group_number[n][m]==0)
			{
				type_number[n]=m;
				break;
			}


	for(m=0;m<20;m++)
		printf("%4d %s %s %8.3lf\n",m,residue[m],residue_simple[m],mean_freq[m]);

	for(m=0;m<20;m++)
	{
		printf("%4d	",m);
		for(n=0;n<=2;n++)
			printf("%4d ",group_type[m][n]);
		printf("\n");
	}

	for(n=0;n<=2;n++)
		for(m=0;m<20;m++)
			printf("%4d %4d %4d %8.3lf\n",n,m,group_number[n][m],group_freq[n][m]);

	for(n=0;n<=2;n++)
		printf("%4d %4d\n",n,type_number[n]);
	
}
//------------------------------------------Read the information for the domain family-------------------------------------------
void ReadCoevolutionComplex()
{
	int m,n;

	if((fp=fopen("CoevolutionComplex.dat","r"))==NULL)
	{
		printf("can't find the file!\n");
		exit(0);
	}

	m=0,n=0;
	while(feof(fp)==0)
	{
		m++;
		fscanf(fp,"%d",&n);
		fscanf(fp,"%s",complex_name[m]);
		fscanf(fp,"%d",&protein_one_length[m]);
		fscanf(fp,"%d",&protein_two_length[m]);
		fscanf(fp,"%d",&complex_length[m]);
		fscanf(fp,"%*[^\n]");
									
		printf("%4d %10s %4d %4d\n",m,complex_name[m],protein_one_length[m],protein_two_length[m]);

		if(m==COMNUM)
			break;	
	}
	fclose(fp);
}
//-----------------------------------Read native sequences from pfam-------------------------------------------------------------
void ReadNativeSequence(int x)
{
	int n;
	char temps1[100],temps2[100],temps3[100];

	sprintf(filename,"D:\\Work\\Program\\Evolution\\Coevolution\\OriginalData\\Sequence\\TCS\\StandardDataset\\Standard_%s_dataset.fasta",complex_name[x]);
	if((fp=fopen(filename,"r"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}

	n=0;
	while(feof(fp)==0)
	{
		n++;
		fscanf(fp,"%s",sequence_name[n]);
		fscanf(fp,"%s",temps1);
		fscanf(fp,"%s",temps2);
		fscanf(fp,"%s",temps3);
		strcpy(sequence_id[n],temps1);
		strcat(sequence_id[n],temps2);
		strcat(sequence_id[n],temps3);

		fscanf(fp,"%*[^\n]");
	}
	fclose(fp);
	seqnum[x]=n-1;

	printf("The sequence number is %4d\n %s\n %s\n",n-1,sequence_name[n-1],sequence_id[n-1]);

}
//-------------------------------------------------Preprocess the sequences------------------------------------------------------
void PreprocessNativeSequence(int x)
{
	int m,n,o,p,q,flag[200];
	double r;

	sprintf(filename,"%s_StandardDataset_Processed.dat",complex_name[x]);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}

	for(m=0;m<200;m++)
		for(p=0;p<NA;p++)
			hist_position_aa[m][p]=0;

	o=strlen(sequence_id[seqnum[x]]);
	q=0;
	for(n=1;n<=seqnum[x];n++)
	{
		q=0;
		for(m=0;m<o;m++)
		{
			flag[m]=0;
			for(p=0;p<NA;p++)
				if(sequence_id[n][m]==residue_simple[p][0])
				{
					flag[m]=1;
					sequence_origin[n][m]=p;
					q++;

					hist_position_aa[m][p]++;
				}

			if(flag[m]==0)
			{
				sequence_origin[n][m]=NA;
				hist_position_aa[m][NA]++;
			}

		}
	}

	sequence_length[x]=0,sequence_num[x]=0;
	for(m=0;m<o;m++)
	{
		flag_aa[m]=0;
		r=double(hist_position_aa[m][NA])/double(seqnum[x]);
		if(r<POS_CUT)					//The positition with residue frequence larger than POS_CUT=0.4 are considered
		{
			sequence_length[x]++;
			flag_aa[m]=1;
		}
		else
			printf("%4d position is low frequent...%8.3lf\n",m,r);
	}

	for(n=1;n<=seqnum[x];n++)
	{
		q=0,flag_seq[n]=0;
		for(m=0;m<o;m++)
			for(p=0;p<NA;p++)
				if(sequence_id[n][m]==residue_simple[p][0]&&flag_aa[m]==1)
					q++;

		if(double(q)/sequence_length[x]>=1-SEQ_CUT)			//24/30=0.8 is the cutoff for removing sequences with gaps>30-24
		{
			flag_seq[n]=1;
			sequence_num[x]++;
		}			
		else
			printf("%4d sequence is removed...%8.3lf %s\n",n,double(q)/sequence_length[x],sequence_id[n]);
	}

	for(n=1;n<=seqnum[x];n++)
		for(m=n+1;m<=seqnum[x];m++)
			if(flag_seq[n]==1&&flag_seq[m]==1)
			{
				if(strcmp(sequence_id[n],sequence_id[m])==0)			//Remove the repeated sequences
				{
					flag_seq[n]=0;
					sequence_num[x]--;
					printf("%4d sequence is repeated...%s\n",n,sequence_id[n]);
					printf("%4d sequence is repeated...%s\n\n",m,sequence_id[m]);
				}
			}

	for(n=1;n<=seqnum[x];n++)
		if(flag_seq[n]==1)
			fprintf(fp,"%260s     %s\n",sequence_name[n],sequence_id[n]);
	fclose(fp);
//-----------------------------------------
	sprintf(filename,"%s_StandardDataset_Processed_Fasta.dat",complex_name[x]);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}
	for(n=1;n<=seqnum[x];n++)
		if(flag_seq[n]==1)
			fprintf(fp,"%s\n%s\n",sequence_name[n],sequence_id[n]);
	fclose(fp);

	printf("%4d %4d %4d %4d\n",seqnum[x],sequence_num[x],complex_length[x],sequence_length[x]);
}
//-----------------------------------------Calculate residue frequency-----------------------------------------------------------
void CalAAFrequency(int x)
{
	int m,n,o,p;
	
	sprintf(filename,"%s_StandardDataset_AAFrequency.dat",complex_name[x]);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}

	for(m=0;m<200;m++)
		for(p=0;p<NA;p++)
			hist_position_aa[m][p]=0;

	o=strlen(sequence_id[seqnum[x]]);
	for(n=1;n<=seqnum[x];n++)
		if(flag_seq[n]==1)
		{
			for(m=0;m<o;m++)
				if(flag_aa[m]==1)
				{
					for(p=0;p<NA;p++)
						if(sequence_id[n][m]==residue_simple[p][0])
							hist_position_aa[m][p]++;
				}
		}

	for(m=0;m<o;m++)
		if(flag_aa[m]==1)
		{
			sum_freq_position[m]=0;
			for(p=0;p<NA;p++)
			{
				freq_position_aa[m][p]=double(hist_position_aa[m][p])/sequence_num[x];
				sum_freq_position[m]+=freq_position_aa[m][p];
			}
		}

//	for(m=0;m<o;m++)
//		if(flag_aa[m]==1)
//			for(p=0;p<NA;p++)
//				freq_position_aa[m][p]=freq_position_aa[m][p]/sum_freq_position[m];   //Normalize the frequences by considering the gaps

	for(p=0;p<NA;p++)
	{
		fprintf(fp,"%4d %s ",p,residue_simple[p]);
		for(m=0;m<o;m++)
			if(flag_aa[m]==1)
				fprintf(fp,"%9.4lf",freq_position_aa[m][p]);
		fprintf(fp,"\n");
	}
	fprintf(fp,"Total  ");
	for(m=0;m<o;m++)
		if(flag_aa[m]==1)
			fprintf(fp,"%9.5lf",sum_freq_position[m]);
	fclose(fp);
//----------------------------------------
	sprintf(filename,"%s_StandardDataset_AAFrequency_RR.dat",complex_name[x]);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}
	for(p=0;p<NA;p++)
	{
		fprintf(fp,"%4d %s ",p,residue_simple[p]);
		for(m=64;m<o;m++)
			if(flag_aa[m]==1)
				fprintf(fp,"%9.4lf",freq_position_aa[m][p]);
		fprintf(fp,"\n");
	}
	fprintf(fp,"Total  ");
	for(m=64;m<o;m++)
		if(flag_aa[m]==1)
			fprintf(fp,"%9.5lf",sum_freq_position[m]);
	fclose(fp);
}
//-------------------------------------Calculate the frequencies for AA Groups---------------------------------------------------
void CalGroupFrequency(int x)
{
	int m,n,o,p;
	
	sprintf(filename,"%s_StandardDataset_GroupFrequency.dat",complex_name[x]);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}

	o=strlen(sequence_id[seqnum[x]]);
	for(n=0;n<=2;n++)
	{
		for(m=0;m<o;m++)
			if(flag_aa[m]==1)
			{
				sum_freq_position[m]=0;
				for(p=0;p<NA;p++)
					freq_position_group[n][m][group_type[p][n]]+=freq_position_aa[m][p];

				for(p=0;p<type_number[n];p++)
					sum_freq_position[m]+=freq_position_group[n][m][p];
			}

		for(p=0;p<type_number[n];p++)
		{
			fprintf(fp,"%4d	",p);
			for(m=0;m<o;m++)
				if(flag_aa[m]==1)
					fprintf(fp,"%9.4lf",freq_position_group[n][m][p]);
			fprintf(fp,"\n");
		}
		fprintf(fp,"Total  ");
		for(m=0;m<o;m++)
			if(flag_aa[m]==1)
				fprintf(fp,"%9.5lf",sum_freq_position[m]);
		fprintf(fp,"\n");
	}
	fclose(fp);
//--------------------------------------------
	sprintf(filename,"%s_StandardDataset_GroupFrequency_RR.dat",complex_name[x]);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}

	o=strlen(sequence_id[seqnum[x]]);
	for(n=0;n<=2;n++)
	{
		for(p=0;p<type_number[n];p++)
		{
			fprintf(fp,"%4d	",p);
			for(m=64;m<o;m++)
				if(flag_aa[m]==1)
					fprintf(fp,"%9.4lf",freq_position_group[n][m][p]);
			fprintf(fp,"\n");
		}
		fprintf(fp,"Total  ");
		for(m=64;m<o;m++)
			if(flag_aa[m]==1)
				fprintf(fp,"%9.5lf",sum_freq_position[m]);
		fprintf(fp,"\n");
	}
	fclose(fp);
}
//-----------------------------------------Calculate position (first-order) conservation-----------------------------------------
void CalPositionConservation(int x)
{
	int m,n,o,p,q;
	double entropy;

	sprintf(filename,"%s_StandardDataset_PositionConservation.dat",complex_name[x]);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}

	q=0,o=strlen(sequence_id[seqnum[x]]);
	for(m=0;m<o;m++)
	{
		if(flag_aa[m]==1)
		{	
			q++;
			position_conservation[m]=0;
			for(p=0;p<NA;p++)
			{
				if(freq_position_aa[m][p]==0)						//sometimes some residues are not occurred in the some positions
					entropy=0;
				else
					entropy=freq_position_aa[m][p]*log(freq_position_aa[m][p]/mean_freq[p]);

				position_conservation[m]+=entropy;
			}			
			fprintf(fp,"%4d %4d %8.3lf",q,m+1,position_conservation[m]);

			for(n=0;n<=2;n++)
			{
				position_conservation[m]=0;
				for(p=0;p<type_number[n];p++)
				{
					if(freq_position_group[n][m][p]==0)						//sometimes some residues are not occurred in the some positions
						entropy=0;
					else
						entropy=freq_position_group[n][m][p]*log(freq_position_group[n][m][p]/group_freq[n][p]);

					position_conservation[m]+=entropy;
				}
				fprintf(fp,"%8.3lf",position_conservation[m]);
			}	
			fprintf(fp,"\n");	
		}	
	}
	fclose(fp);
//------------------------------------------------------------
	sprintf(filename,"%s_StandardDataset_PositionConservation_RR.dat",complex_name[x]);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}

	q=0,o=strlen(sequence_id[seqnum[x]]);
	for(m=64;m<o;m++)
	{
		if(flag_aa[m]==1)
		{	
			q++;
			position_conservation[m]=0;
			for(p=0;p<NA;p++)
			{
				if(freq_position_aa[m][p]==0)						//sometimes some residues are not occurred in some positions
					entropy=0;
				else
					entropy=freq_position_aa[m][p]*log(freq_position_aa[m][p]/mean_freq[p]);

				position_conservation[m]+=entropy;
			}			
			fprintf(fp,"%4d %4d ",q,m+1);
			if(q+3<70)
				fprintf(fp,"%4d %4d %8.3lf",q+3,q+4,position_conservation[m]);
			else 
				fprintf(fp,"%4d %4d %8.3lf",q+5,q+6,position_conservation[m]);

			for(n=0;n<=2;n++)
			{
				position_conservation[m]=0;
				for(p=0;p<type_number[n];p++)
				{
					if(freq_position_group[n][m][p]==0)						//sometimes some residues are not occurred in the some positions
						entropy=0;
					else
						entropy=freq_position_group[n][m][p]*log(freq_position_group[n][m][p]/group_freq[n][p]);

					position_conservation[m]+=entropy;
				}
				fprintf(fp,"%8.3lf",position_conservation[m]);
			}	
			fprintf(fp,"\n");	
		}	
	}
	fclose(fp);
}
//-----------------------------------------Calculate residue pair frequency------------------------------------------------------
void CalPairFrequency(int x)
{
	int m,n,o,p,q,r;
	
	sprintf(filename,"%s_StandardDataset_SumPairFrequency.dat",complex_name[x]);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}

	for(m=0;m<200;m++)
		for(n=0;n<200;n++)
			for(o=0;o<NA;o++)
				for(p=0;p<NA;p++)
					hist_pair_aa[m][n][o][p]=0;

	o=strlen(sequence_id[seqnum[x]]);
	for(n=1;n<=seqnum[x];n++)
		if(flag_seq[n]==1)
		{
			for(m=0;m<o;m++)
				if(flag_aa[m]==1)
					for(r=0;r<o;r++)
						if(flag_aa[r]==1)
						{
							for(p=0;p<NA;p++)
								if(sequence_id[n][m]==residue_simple[p][0])
									for(q=0;q<NA;q++)
										if(sequence_id[n][r]==residue_simple[q][0])
											hist_pair_aa[m][r][p][q]++;
						}
		}

	for(m=0;m<o;m++)
		if(flag_aa[m]==1)
			for(r=0;r<o;r++)
				if(flag_aa[r]==1)
				{
					sum_freq_pair[m][r]=0;
					for(p=0;p<NA;p++)
						for(q=0;q<NA;q++)
						{
							freq_pair_aa[m][r][p][q]=double(hist_pair_aa[m][r][p][q])/sequence_num[x];
							sum_freq_pair[m][r]+=freq_pair_aa[m][r][p][q];
						}
				}
	
	r=0;
	for(m=0;m<o;m++)
		if(flag_aa[m]==1)
		{
			r++;
			fprintf(fp,"%4d %4d",r,m+1);
			for(n=0;n<o;n++)
				if(flag_aa[n]==1)	
					fprintf(fp,"%9.4lf",sum_freq_pair[m][n]);
			fprintf(fp,"\n");
		}
	fclose(fp);
}
//-----------------------------------Calculate the coupling (second-order) conservation------------------------------------------
void CalCouplingConservation(int x)
{
	int m,n,o,p,q,r;

	sprintf(filename,"%s_StandardDataset_CouplingConservation.dat",complex_name[x]);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}

	o=strlen(sequence_id[seqnum[x]]);
	for(m=0;m<o;m++)
		if(flag_aa[m]==1)
			for(p=0;p<NA;p++)
			{
				if(freq_position_aa[m][p]==0||freq_position_aa[m][p]==1)						//sometimes some residues are not occurred in the some positions
					phi[m][p]=0;
				else
					phi[m][p]=log((freq_position_aa[m][p]*(1-mean_freq[p]))/((1-freq_position_aa[m][p])*mean_freq[p]));
			}

	for(m=0;m<o;m++)
		if(flag_aa[m]==1)
			for(r=0;r<o;r++)
				if(flag_aa[r]==1)
					for(p=0;p<NA;p++)
						for(q=0;q<NA;q++)
						{
							freq_coupling[m][r][p][q]=freq_pair_aa[m][r][p][q]-freq_position_aa[m][p]*freq_position_aa[r][q];
							weight_freq_coupling[m][r][p][q]=phi[m][p]*phi[r][q]*freq_coupling[m][r][p][q];
							sum_freq_coupling[m][r]+=weight_freq_coupling[m][r][p][q]*weight_freq_coupling[m][r][p][q];
						}

	r=0;
	for(m=0;m<o;m++)
		if(flag_aa[m]==1)
		{
			r++;
			fprintf(fp,"%4d %4d",r,m+1);
			for(n=0;n<o;n++)
				if(flag_aa[n]==1)	
					fprintf(fp,"%9.4lf",sqrt(sum_freq_coupling[m][n]));
			fprintf(fp,"\n");
		}
	fclose(fp);
	
	sprintf(filename,"%s_StandardDataset_CouplingConservation_Matrix.dat",complex_name[x]);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}
	r=0;
	for(m=0;m<o;m++)
		if(flag_aa[m]==1)
		{
			r++;
			q=0;
			for(n=0;n<o;n++)
				if(flag_aa[n]==1)	
				{
					q++;
					fprintf(fp,"%4d %4d %4d %4d %9.4lf\n",r,q,m+1,n+1,sqrt(sum_freq_coupling[m][n]));
				}
		}
	fclose(fp);
//-----------------------------------------
	sprintf(filename,"%s_StandardDataset_CouplingConservation_RR.dat",complex_name[x]);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}
	r=0;
	for(m=64;m<o;m++)					//64 is the start of the RR domain
		if(flag_aa[m]==1)
		{
			r++;
			fprintf(fp,"%4d %4d",r,m+1);
			for(n=64;n<o;n++)
				if(flag_aa[n]==1)	
					fprintf(fp,"%9.4lf",sqrt(sum_freq_coupling[m][n]));
			fprintf(fp,"\n");
		}
	fclose(fp);

	sprintf(filename,"%s_StandardDataset_CouplingConservation_Matrix_RR.dat",complex_name[x]);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}
	r=0;
	for(m=64;m<o;m++)
		if(flag_aa[m]==1)
		{
			r++;
			q=0;
			for(n=64;n<o;n++)
				if(flag_aa[n]==1)	
				{
					q++;
					if(r+3<70&&q+3<70)
						fprintf(fp,"%4d %4d %4d %4d %4d %4d %4d %4d %9.4lf\n",r,q,m+1,n+1,r+3,q+3,r+4,q+4,sqrt(sum_freq_coupling[m][n]));
					else if(r+3<70&&q+3>=70)
						fprintf(fp,"%4d %4d %4d %4d %4d %4d %4d %4d %9.4lf\n",r,q,m+1,n+1,r+3,q+5,r+4,q+6,sqrt(sum_freq_coupling[m][n]));
					else if(r+3>=70&&q+3<70)
						fprintf(fp,"%4d %4d %4d %4d %4d %4d %4d %4d %9.4lf\n",r,q,m+1,n+1,r+5,q+3,r+6,q+4,sqrt(sum_freq_coupling[m][n]));
					else
						fprintf(fp,"%4d %4d %4d %4d %4d %4d %4d %4d %9.4lf\n",r,q,m+1,n+1,r+5,q+5,r+6,q+6,sqrt(sum_freq_coupling[m][n]));					
				}
		}
	fclose(fp);
}
//------------------------------------------Separate strong and weak coupling----------------------------------------------------
void SeparateStrongWeakCoupling(int x)
{
	int m,n,o,p,q,r,s;

	sprintf(filename,"%s_StandardDataset_CouplingConservation_Matrix_Updiagonal.dat",complex_name[x]);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}

	o=strlen(sequence_id[seqnum[x]]);
	printf("The length is %4d\n",o);
	r=0,s=0;
	for(m=0;m<o;m++)
		if(flag_aa[m]==1)
		{
			r++;
			q=r;
			for(n=m+1;n<o;n++)
				if(flag_aa[n]==1)	
				{
					q++;
					fprintf(fp,"%4d %4d %4d %4d %9.4lf\n",r,q,m+1,n+1,sqrt(sum_freq_coupling[m][n]));

					s++;
					couple_value[s]=sqrt(sum_freq_coupling[m][n]);
					couple_position[s][0]=r;
					couple_position[s][1]=q;
				}
		}
	fclose(fp);

	sprintf(filename,"%s_StandardDataset_CouplingConservation_Matrix_RankCoupling.dat",complex_name[x]);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}
	for(m=1;m<=s;m++)
	{
		rank_couple[m]=1;
		for(n=1;n<=s;n++)
			if(couple_value[n]>couple_value[m])
				rank_couple[m]++;
			else if(couple_value[n]==couple_value[m]&&m!=n)
			{
				if(m<n)
					rank_couple[m]++;
			}
		rank_couple_index[rank_couple[m]]=m;
	}

	for(m=1;m<=s;m++)
		fprintf(fp,"%4d %4d %4d %8.3lf\n",m,couple_position[rank_couple_index[m]][0],couple_position[rank_couple_index[m]][1],couple_value[rank_couple_index[m]]);
	fclose(fp);
//-----------------------------------------------------
	sprintf(filename,"%s_StandardDataset_CouplingConservation_Matrix_RR_Updiagonal.dat",complex_name[x]);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}

	for(m=0;m<15000;m++)
	{
		couple_value[m]=0;
		couple_position[m][0]=0;
		couple_position[m][1]=0;
		rank_couple[m]=0;
		rank_couple_index[m]=0;
	}

	r=0,s=0;
	for(m=64;m<o;m++)
		if(flag_aa[m]==1)
		{
			r++;
			q=r;
			for(n=m+1;n<o;n++)
				if(flag_aa[n]==1)	
				{
					q++;
					if(r+3<70&&q+3<70)
						fprintf(fp,"%4d %4d %4d %4d %4d %4d %4d %4d %9.4lf\n",r,q,m+1,n+1,r+3,q+3,r+4,q+4,sqrt(sum_freq_coupling[m][n]));
					else if(r+3<70&&q+3>=70)
						fprintf(fp,"%4d %4d %4d %4d %4d %4d %4d %4d %9.4lf\n",r,q,m+1,n+1,r+3,q+5,r+4,q+6,sqrt(sum_freq_coupling[m][n]));
					else if(r+3>=70&&q+3<70)
						fprintf(fp,"%4d %4d %4d %4d %4d %4d %4d %4d %9.4lf\n",r,q,m+1,n+1,r+5,q+3,r+6,q+4,sqrt(sum_freq_coupling[m][n]));
					else
						fprintf(fp,"%4d %4d %4d %4d %4d %4d %4d %4d %9.4lf\n",r,q,m+1,n+1,r+5,q+5,r+6,q+6,sqrt(sum_freq_coupling[m][n]));	

					s++;
					couple_value[s]=sqrt(sum_freq_coupling[m][n]);
					couple_position[s][0]=r;
					couple_position[s][1]=q;
				}
		}
	fclose(fp);

	sprintf(filename,"%s_StandardDataset_CouplingConservation_Matrix_RR_RankCoupling.dat",complex_name[x]);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}
	for(m=1;m<=s;m++)
	{
		rank_couple[m]=1;
		for(n=1;n<=s;n++)
			if(couple_value[n]>couple_value[m])
				rank_couple[m]++;
			else if(couple_value[n]==couple_value[m]&&m!=n)
			{
				if(m<n)
					rank_couple[m]++;
			}
		rank_couple_index[rank_couple[m]]=m;
	}

	for(m=1;m<=s;m++)
	{
		if(couple_position[rank_couple_index[m]][0]+3<70&&couple_position[rank_couple_index[m]][1]+3<70)
			fprintf(fp,"%4d %4d %4d %4d %4d %4d %4d %8.3lf\n",m,couple_position[rank_couple_index[m]][0],couple_position[rank_couple_index[m]][1],
			couple_position[rank_couple_index[m]][0]+3,couple_position[rank_couple_index[m]][1]+3,couple_position[rank_couple_index[m]][0]+4,
			couple_position[rank_couple_index[m]][1]+4,couple_value[rank_couple_index[m]]);
		else if(couple_position[rank_couple_index[m]][0]+3<70&&couple_position[rank_couple_index[m]][1]+3>=70)
			fprintf(fp,"%4d %4d %4d %4d %4d %4d %4d %8.3lf\n",m,couple_position[rank_couple_index[m]][0],couple_position[rank_couple_index[m]][1],
			couple_position[rank_couple_index[m]][0]+3,couple_position[rank_couple_index[m]][1]+5,couple_position[rank_couple_index[m]][0]+4,
			couple_position[rank_couple_index[m]][1]+6,couple_value[rank_couple_index[m]]);
		else if(couple_position[rank_couple_index[m]][0]+3>=70&&couple_position[rank_couple_index[m]][1]+3<70)
			fprintf(fp,"%4d %4d %4d %4d %4d %4d %4d %8.3lf\n",m,couple_position[rank_couple_index[m]][0],couple_position[rank_couple_index[m]][1],
			couple_position[rank_couple_index[m]][0]+5,couple_position[rank_couple_index[m]][1]+3,couple_position[rank_couple_index[m]][0]+6,
			couple_position[rank_couple_index[m]][1]+4,couple_value[rank_couple_index[m]]);
		else
			fprintf(fp,"%4d %4d %4d %4d %4d %4d %4d %8.3lf\n",m,couple_position[rank_couple_index[m]][0],couple_position[rank_couple_index[m]][1],
			couple_position[rank_couple_index[m]][0]+5,couple_position[rank_couple_index[m]][1]+5,couple_position[rank_couple_index[m]][0]+6,
			couple_position[rank_couple_index[m]][1]+6,couple_value[rank_couple_index[m]]);
	}
	fclose(fp);
}
//---------------------------------------------Calculate group pair frequency----------------------------------------------------
void CalGroupPairFrequency(int x, int y)
{
	int m,n,o,p,q,r;

	sprintf(filename,"%s_StandardDataset_SumPairFrequency_%d.dat",complex_name[x],y);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}

	for(m=0;m<200;m++)
		for(n=0;n<200;n++)
			for(o=0;o<NA;o++)
				for(p=0;p<NA;p++)
					freq_pair_group[m][n][o][p]=0;

	o=strlen(sequence_id[seqnum[x]]);
	for(m=0;m<o;m++)
		if(flag_aa[m]==1)
			for(r=0;r<o;r++)
				if(flag_aa[r]==1)
				{
					sum_freq_pair[m][r]=0;
					for(p=0;p<NA;p++)
						for(q=0;q<NA;q++)
							freq_pair_group[m][r][group_type[p][y]][group_type[q][y]]+=freq_pair_aa[m][r][p][q];
					
					for(p=0;p<type_number[y];p++)
						for(q=0;q<type_number[y];q++)
							sum_freq_pair[m][r]+=freq_pair_group[m][r][p][q];
				}
	
	r=0;
	for(m=0;m<o;m++)
		if(flag_aa[m]==1)
		{
			r++;
			fprintf(fp,"%4d %4d",r,m+1);
			for(n=0;n<o;n++)
				if(flag_aa[n]==1)	
					fprintf(fp,"%9.4lf",sum_freq_pair[m][n]);
			fprintf(fp,"\n");
		}

	fclose(fp);
}
//-----------------------------------Calculate group coupling (second-order) conservation----------------------------------------
void CalGroupCouplingConservation(int x, int y)
{
	int m,n,o,p,q,r;

	sprintf(filename,"%s_StandardDataset_CouplingConservation_%d.dat",complex_name[x],y);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}

	for(m=0;m<200;m++)
		for(r=0;r<200;r++)
			sum_freq_coupling[m][r]=0;

	o=strlen(sequence_id[seqnum[x]]);
	for(m=0;m<o;m++)
		if(flag_aa[m]==1)
			for(p=0;p<type_number[y];p++)
			{
				if(freq_position_group[y][m][p]==0||freq_position_group[y][m][p]==1)						//sometimes some residues are not occurred in some positions
					phi[m][p]=0;
				else
					phi[m][p]=log((freq_position_group[y][m][p]*(1-group_freq[y][p]))/((1-freq_position_group[y][m][p])*group_freq[y][p]));
			}

	for(m=0;m<o;m++)
		if(flag_aa[m]==1)
			for(r=0;r<o;r++)
				if(flag_aa[r]==1)
					for(p=0;p<type_number[y];p++)
						for(q=0;q<type_number[y];q++)
						{
							freq_coupling[m][r][p][q]=freq_pair_group[m][r][p][q]-freq_position_group[y][m][p]*freq_position_group[y][r][q];
							weight_freq_coupling[m][r][p][q]=phi[m][p]*phi[r][q]*freq_coupling[m][r][p][q];
							sum_freq_coupling[m][r]+=weight_freq_coupling[m][r][p][q]*weight_freq_coupling[m][r][p][q];
						}

	r=0;
	for(m=0;m<o;m++)
		if(flag_aa[m]==1)
		{
			r++;
			fprintf(fp,"%4d %4d",r,m+1);
			for(n=0;n<o;n++)
				if(flag_aa[n]==1)	
					fprintf(fp,"%9.4lf",sqrt(sum_freq_coupling[m][n]));
			fprintf(fp,"\n");
		}
	fclose(fp);
	
	sprintf(filename,"%s_StandardDataset_CouplingConservation_Matrix_%d.dat",complex_name[x],y);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}
	r=0;
	for(m=0;m<o;m++)
		if(flag_aa[m]==1)
		{
			r++;
			q=0;
			for(n=0;n<o;n++)
				if(flag_aa[n]==1)	
				{
					q++;
					fprintf(fp,"%4d %4d %4d %4d %9.4lf\n",r,q,m+1,n+1,sqrt(sum_freq_coupling[m][n]));
				}
		}
	fclose(fp);
//-------------------------------------------------------	
	sprintf(filename,"%s_StandardDataset_CouplingConservation_%d_RR.dat",complex_name[x],y);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}
	r=0;
	for(m=64;m<o;m++)
		if(flag_aa[m]==1)
		{
			r++;
			fprintf(fp,"%4d %4d",r,m+1);
			for(n=64;n<o;n++)
				if(flag_aa[n]==1)	
					fprintf(fp,"%9.4lf",sqrt(sum_freq_coupling[m][n]));
			fprintf(fp,"\n");
		}
	fclose(fp);

	sprintf(filename,"%s_StandardDataset_CouplingConservation_Matrix_%d_RR.dat",complex_name[x],y);
	if((fp=fopen(filename,"w"))==NULL)
	{
		printf("can't find the file %s!\n",filename);
		exit(0);
	}
	r=0;
	for(m=64;m<o;m++)
		if(flag_aa[m]==1)
		{
			r++;
			q=0;
			for(n=64;n<o;n++)
				if(flag_aa[n]==1)	
				{
					q++;
					fprintf(fp,"%4d %4d %4d %4d %9.4lf\n",r,q,m+1,n+1,sqrt(sum_freq_coupling[m][n]));
				}
		}
	fclose(fp);
}
//-------------------------------------------------Main loop---------------------------------------------------------------------
void main()
{
//	ReadAAMeanFrequency();
	ReadAAGroup();
	ReadCoevolutionComplex();

	for(i=1;i<=1;i++)
	{
		ReadNativeSequence(i);
		PreprocessNativeSequence(i);

		CalAAFrequency(i);
		CalGroupFrequency(i);
		CalPositionConservation(i);


		CalPairFrequency(i);
		CalCouplingConservation(i);
		SeparateStrongWeakCoupling(i);
	
		for(j=0;j<=2;j++)
		{
			CalGroupPairFrequency(i,j);
			CalGroupCouplingConservation(i,j);
		}
	}
}