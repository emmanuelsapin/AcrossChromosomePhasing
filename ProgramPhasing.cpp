/***************************************************************************************************************************************************
*
*                                 This program will read segment and output condense files
*                                 
****************************************************************************************************************************************************/
#include <time.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <errno.h>
#include <sys/ipc.h> 
#include <sys/shm.h> 
#include <sys/sysinfo.h>
#include <unistd.h>
#include <sched.h>
#include <sys/syscall.h>
#include <omp.h>

#include "readinteger.h"
#include "readreal.h"
#include "readnegativereal.h"

 
#define NMAX 664 
#define NCHR 23
#define NSNPPERCHR 26230
#define MAXSNP 500005
#define MAXPOP 435188
#define NBGROUPPERPERSON NMAX+1
#define MAXSEGMENTPAIR 200
#define MAXSUMLENPAIR 35000
#define MAXCLOSERELAT175 310000
#define MAXCLOSERELAT 6000
#define MAXCLOSERELATTEMP 450000
#define MAXID 6026457
#define MINPIHAT 0
#define MAXPIHAT 1
#define NBINMATRIX1 100
#define NBINMATRIX2 200
#define NBINMATRIX3 500
#define NBINMATRIX4 1000
#define NBINMATRIX5 2000
#define NBINMATRIX6 5000
#define NBINMATRIX7 10000
#define nbbyteforsegment 7
#define MAXSEGINDIVTEMP 550000
#define MAXSEGINDIV 10000
#define MAXNBPAIR  152252 
#define MAXNBPO 3719
#define MAXNBDIVISOR 25
#define NBINDIVMAX 100000
#define NBINDIV 100000

#define MAXNBTRIO 978
#define NBCLOSER 4351
#define MAXFILE 6
#define NBINDIVEA 1
#define NBTOURNAMENT 0
#define MAXGEN 1 
#define MAXVARIABLE 28
#define INCREMENTLOOP 1
#define INCREMENTLOOPMISS 3
#define INCREMENTLOOPNB (3/2)

#define X1SIZE 18
#define Y1SIZE 25
const int XSTART[10] = {648, 71, 136, 191, 263,327,390,448,519,576}; 
#define Y1START 30

int bestpihatagainstallID[100];
int IDbestpihat;
int IDbestpihat2;

float pihatagainstall[MAXPOP];
float pihatagainstall2[MAXPOP];
uint64_t totseuil=0; 
int nbchrdivider[51][23];

int IDjob=0;
int sumtotindivcount=0;
int placefirttoconsider=0;

typedef struct 
{	int start;
	int end;
	int phasing;
	int nbright;
	int nbwrong;
	int segment;
} typechrdivider;

typechrdivider chrdivider[50][23][20];
	
unsigned int focaltotlenghtpos=0;
unsigned int focaltotseg=0;
float seuilpihat[3];
int result[4][6];
float bestpihat[7];
float bestcor=20000;
int printdetail;
int UKBID[MAXID];
int UKBIDback[MAXPOP];
int nbsnpperchr[23];
int nbsnpperchrinfile[23];

int MAF[NSNPPERCHR][23]={0};
	
int genomeoffpss[3][NSNPPERCHR][23];	

int NbIndiv=0;

//int isoffspring[MAXPOP]; //-1 not offspring;>-1 trio number 
//int isparent[MAXPOP]; //-1 not offspring;>-1 trio number 
//int nblistrelat[MAXPOP];
//int listrelat[MAXPOP][13];

typedef struct
{	int chi;
	unsigned char haptohap;
} resultphaing;

typedef struct
{	int ID;
	unsigned char coef;
} indivclose;

typedef struct
{	int ID;
	int PIHAT10000;
} structrelatif;
	
#define NBMAXRELATIVE 450
	
unsigned char * genomes[23];

typedef struct 
{	int chr;
	int snpstart;
	int snpend;
	float cor1;
	float cor2;	
	int nblastsplit;	
} structphaseerror;


typedef struct 
{	int IDoffspring;
	int IDp1;
	int IDp2;
	float score;
	int nbbad;
	int nbkeep;
} structtrio;

structtrio tabtrio[979];

typedef struct
{	double cor;
	int nbgroup1;
	int nbgroup2;
	
} pointdecision;

pointdecision tappointdec[MAXNBDIVISOR*22];
		
typedef struct 
{	int IDoffspring;
	int IDp1;
	float score;
	int nbphaseerror;
} structPO;

structPO tabPO[MAXNBPO];


int nbbreak=0;
#define MAXBREAK 1000

typedef struct 
{	int ID;
	int start;
	int end;
	int chr;
	int hapfocal;
	int haprelat;
	int averageseg;
	int nbhet;
} structseg;
//int seuil[23][23];

//unsigned char * genomes;

uint64_t * tempload;

//unsigned int proportion[MAXNBTRIO][23];
	
uint64_t sumlenseg=0;
uint64_t sumlenseghet=0;
int nbsnpwoseg=0;
	
int distrigametic[12];	
int distrigametickeep[23];	
/*
typedef struct 
{	unsigned long long pos;
//	int cm;
} structSNP;

structSNP map[NSNPPERCHR];*/


//snpinteger genomeinteger[MAXPOP][5230];//NSNPPERCHR
//5227
//int listrelativestart[MAXPOP];
//int listrelativeend[MAXPOP];
//int listrelativeP1[MAXPOP];
//int listrelativeP2[MAXPOP];
	
//structSNP map[NCHR][NSNPPERCHR];

void delay(int number_of_seconds)  
{ 
    // Converting time into milli_seconds 
    int milli_seconds = 1000 * number_of_seconds; 
  
    // Stroing start time 
    clock_t start_time = clock(); 
  
    // looping till required time is not acheived 
    while (clock() < start_time + milli_seconds); 
} 

int readID(char pathID[])
{	FILE * fileID;
	if ((fileID = fopen(pathID, "r")) == NULL) 
	{	if (printdetail) printf("file %s is not found\n",pathID);		
		return (1);
	}; 
	if (printdetail) printf("file %s opened\n",pathID);
	char first; 
	//do {first=getc(fileID);} while (first!='\n'); 
	//do {first=getc(fileID);} while (first!='\n'); 
	int place=1;
	int ID1;
	for(int ID1=0;ID1<MAXID;ID1++)
	{	UKBID[ID1]=-1;
	};
	for(int ID1=0;ID1<MAXPOP;ID1++)
	{	UKBIDback[ID1]=0;
	};
	
	do  
	{	ID1=readinteger(fileID);
		if (ID1==0)
		{	if (printdetail) printf("ID1=0 at %d",place);
			exit(1);
		};
		if ((ID1%10000)==0) printf("%d\n",ID1);
		if (ID1<INT32_MAX) 
		{	UKBID[ID1]=place;
			UKBIDback[place]=ID1;
			place++;
	
			do {first=getc(fileID);} while (first!='\n'); 
		}
	} while (place<MAXPOP && ID1<INT32_MAX);
	fclose(fileID); 
	return 0; 
};


unsigned char * buffer;

int readtrio()
{	for(int nbtrio=0;nbtrio<MAXPOP;nbtrio++)
	{	//isoffspring[nbtrio]=-1;
	};
	FILE * triofile;
	if ((triofile = fopen("/pl/active/KellerLab/Emmanuel/buildtree/true_trio_dat.csv", "r")) == NULL) 
	{	if (printdetail) printf("file triofile is not found\n");		
		return (1);
	};
	if (printdetail) printf("file triofile opened\n");		
	char first;
	do {first=getc(triofile);} while (first!='\n');
	for(int nbtrio=0;nbtrio<MAXNBTRIO;nbtrio++)
	{	tabtrio[nbtrio].IDoffspring=readinteger(triofile);
		tabtrio[nbtrio].IDp1=readinteger(triofile);
		do {first=getc(triofile);} while (first!=','); 
		do {first=getc(triofile);} while (first!=',');
		tabtrio[nbtrio].IDp2=readinteger(triofile);
		do {first=getc(triofile);} while (first!='\n'); 
//		if (printdetail)
		printf("%d %d %d\n",tabtrio[nbtrio].IDoffspring,tabtrio[nbtrio].IDp1,tabtrio[nbtrio].IDp2 );
		printf("%d %d %d\n",UKBID[tabtrio[nbtrio].IDoffspring]-1,UKBID[tabtrio[nbtrio].IDp1]-1,UKBID[tabtrio[nbtrio].IDp2]-1 );
		//isoffspring[UKBID[tabtrio[nbtrio].IDoffspring]-1]=nbtrio;
	}
	fclose(triofile);
	return (0);
}

int readgenomelocal(char pathfile[],int chr,int run,int step,unsigned char * gentomodify)
{	FILE * endfile;
	char number[100];
	char pathfilechr[300]; 
	//strcpy(pathfile,"/pl/active/KellerLab/Emmanuel/correctUKB/run1step2/error");
	strcpy(pathfilechr,pathfile);
	sprintf(number, "%d", chr);
	strcat(pathfile,number);
	strcat(pathfile,".hap");
	printf("%s\n",pathfile);
	if ((endfile = fopen(pathfile, "r")) == NULL) 
	{	printf("file end is not found\n");		 
		return (1);
	};
	char first;
	int snp=0;
	//for(snp<nbsnpperchr[chr];snp++)// if (ID<10)26229
	do 
	{//	printf("%d\n",snp);
		first=getc(endfile);
		if (first!=EOF)
		{	if (snp<2) printf("%c ",first);
			first=getc(endfile);
			if (snp<2) printf("%c ",first);
			int i=readinteger(endfile);
			if (snp<2) printf("%d ",i);
			
			/*int stringrank=-1;
			char name[50];
			do 
			{	stringrank++;
				first=getc(endfile);
				if (first!='\0' && first!='\n') name[stringrank]=first;
				
			} while (first!='\0' && first!='\n')
			name[stringrank]='\0';	
			if (strcmp(name,tabsnp[snp].name*/
			do {first=getc(endfile);} while (first!=32);
			i=readinteger(endfile);
			if (snp<2) printf("%d ",i);
			first=getc(endfile);
			
			if (snp<2) printf("%c ",first);
			getc(endfile);
			getc(endfile);
			for (int relat=0;relat<NbIndiv;relat++)//NBIRISH
			{	first=getc(endfile);
				int geno=0;
				first=getc(endfile);
				if (snp<2) printf("%c ",first);
				if (first==49) geno=geno+2;
				first=getc(endfile);
				//if (snp==0) printf("%c ",first);
				first=getc(endfile);
				first=getc(endfile);
				
				
				//if (snp==0) printf("%c ",first);
				first=getc(endfile);
				if (first==49) geno=geno+1;
				*( gentomodify+(unsigned long long) (((((unsigned long long) relat)*(nbsnpperchr[chr]/4+((nbsnpperchr[chr]%4)>0)) )+snp/4)))=
										(*( gentomodify+(unsigned long long) (((((unsigned long long) relat)*(nbsnpperchr[chr]/4+((nbsnpperchr[chr]%4)>0)) )+snp/4))) & 
											(~(3<<((snp%4)*2))))
												|	((geno)<<((snp%4)*2));
				getc(endfile);
				getc(endfile);
			};
		};
		snp++;
	} while (first!=EOF);
	nbsnpperchrinfile[chr]=snp;
	fclose(endfile);
}


int checkrelat[MAXPOP];

int segstart[24];

int nbrelatpihat=0;

int loadsegment(int ID,int numtrio,int IDp1loop,int IDp2loop,int lenminseg,int version,int gentostart,char pathresult[],char PathMAF[])
{	printf("IDjob %d\n",IDjob);
	int tabseuil[10]; 
	tabseuil[0]=1;tabseuil[1]=2;tabseuil[2]=3;tabseuil[3]=4;tabseuil[4]=5;tabseuil[5]=6;tabseuil[6]=8;tabseuil[7]=11;tabseuil[8]=15;tabseuil[9]=20;
	int seuilacrosschr=0;
	int seuilinchr=3;
	int hetsquare=0;
	int prodadd=0;
	int prodaddpersnp=1;
	int keepone=0;	
	float seuilpihat=0.33;
	int LD=0;
	int takelenasweight=0;
	int howcalculfreq=2;
	int MAFintoaccount=1;
	int limitweight=1;
	int limitweightproduct=13000; 
	int removelastchrsassoc=3;
	int addiftwicelimit=0;
	int hetsnp=0;
	int relatpihatP1[MAXCLOSERELATTEMP];
	int relatpihatP2[MAXCLOSERELAT];
	int relatpihatID[MAXPOP];
	int relatpihatchr[MAXCLOSERELAT][23][2];
	int seuilscore=1;
	int powersegment=1;	
	printf("IDjob %d\n",IDjob);
	for(int IDrun=0;IDrun<MAXPOP;IDrun++)// if (ID<10)
	{	relatpihatID[IDrun]=-1;	
	};
	int64_t segnum=0;
	
	FILE * MAFfile;
	if ((MAFfile = fopen(PathMAF, "r")) == NULL) 
	{	if (printdetail) printf("file triofile is not found\n");		
		return (1);
	};
	int chr=1;
	do 
	{	chr=readinteger(MAFfile);	
		if (chr<23)
		{	int snp=readinteger(MAFfile);	
			MAF[snp][chr]=readinteger(MAFfile);	
			
		};
	} while (chr<23);
	fclose(MAFfile);
	
	for(int relat=0;relat<MAXCLOSERELATTEMP;relat++)// if (ID<10)
	{	relatpihatP1[relat]=0;
	//	relatpihatP2[relat]=0;
	};	
	printf("IDjob %d\n",IDjob);
	FILE * endfile;
	uint64_t averageofaverage[23];
	char number[100];
	char pathfile[300]; 								
	int nbseg=0;
	int segstarttemp[24];
	printf("IDjob %d\n",IDjob);
	float bestpihatagainstall[100];
	for(int relat=0;relat<MAXPOP;relat++)// if (ID<10)
	{	pihatagainstall[relat]=0;
			
	};
	FILE * filepihat;
	char pathpihat[300]; 	
	strcat(pathpihat,"/pl/active/KellerLab/Emmanuel/correctphaseerrorbasedonrelative/allpihat");///0.TXT
	sprintf(number, "%d",ID/100);
	strcat(pathpihat,number);
	strcat(pathpihat,".TXT");
	printf("%s\n",pathpihat);
	if ((filepihat = fopen(pathpihat, "r")) == NULL) 
	{	if (printdetail) printf("file pihat is not found\n");		
		return (1); 
	};
	char first;
	for(int relat=0;relat<(ID%100);relat++)// if (ID<10)
	{	do {first=getc(filepihat);} while (first!='\n'); 
	};
	for(int relat=0;relat<100;relat++)// if (ID<10)
	{	readinteger(filepihat);	
		int IDrelat=readinteger(filepihat);
		float pihat=readnegativereal(filepihat);
		float pihat2=readnegativereal(filepihat);
		bestpihatagainstallID[relat]=IDrelat;
		bestpihatagainstall[relat]=pihat;
		pihatagainstall[IDrelat]=pihat;
		pihatagainstall2[IDrelat]=pihat2;
		printf("%d %f %d %f %f\n",relat,bestpihatagainstall[relat],bestpihatagainstallID[relat],pihatagainstall[IDrelat], pihatagainstall2[IDrelat]);
	//	fprintf(filepihat,"%d %f %d %d\n",relat,bestpihatagainstall[relat],bestpihatagainstallID[relat],relatpihatP1[bestpihatagainstallID[relat]] );
	};
	fclose(filepihat);
	int compteur=0;
	for(int comp=0;comp<7;comp++) bestpihat[comp]=0;
	bestpihat[0]=bestpihatagainstall[compteur];
	int nbtimeIDchange=0;
		
	do {compteur++;} while (bestpihatagainstall[compteur]>seuilpihat);
	placefirttoconsider=compteur;
	bestpihat[1]=bestpihatagainstall[compteur];
	IDbestpihat=bestpihatagainstallID[compteur];
	IDbestpihat2=bestpihatagainstallID[compteur+1];
	printf("best ID pihat = %d\n",IDbestpihat);
	for(int comp=compteur;comp<compteur+2;comp++)	bestpihat[2]=bestpihat[2]+bestpihatagainstall[comp];
	for(int comp=compteur;comp<compteur+5;comp++)	bestpihat[3]=bestpihat[3]+bestpihatagainstall[comp];
	for(int comp=compteur;comp<compteur+10;comp++)	bestpihat[4]=bestpihat[4]+bestpihatagainstall[comp];
	for(int comp=compteur;comp<compteur+20;comp++)	bestpihat[5]=bestpihat[5]+bestpihatagainstall[comp];
	for(int comp=compteur;comp<compteur+50;comp++)	bestpihat[6]=bestpihat[6]+bestpihatagainstall[comp];	
	segstarttemp[0]=0;
	
	float pihatagainstallchrMP[MAXPOP][2]={0};
	for(int relat=0;relat<MAXPOP;relat++) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
	{	pihatagainstallchrMP[relat][0]=0;pihatagainstallchrMP[relat][1]=0;
	};
	printf("%d\n",IDjob);
	return (1);
}

int compareressultindivallcombinV15D4(int ID,int chr1, int chr2,int numtrio,float limitnbsnp,int run,int IDp1loop,int IDp2loop,int lenminseg,int version,int gentostart,char pathresult[])
{	printf("%d\n",IDjob);
	seuilpihat[0]=0.33;
	seuilpihat[1]=0.33;
	seuilpihat[2]=0.33;
	int increaseincrementwhenneg=1;
	int takesumintoaccout=1;
	int limitweight=105;
	int howtocalculatenotimprove=0;
	int howtocalculatenotimprove2=0;
	int nbmutdependonnotimprove=0;
	int mutateallhapwhennoimprove=0;
	int changeonlyifneg=2;
	int crossclever=0;
	int numgentowait=17;
	int howtocalculaterelattochr=0; 
	int nbmutvari=1;
	int64_t tabresmut[50][50];
	int nbtimetabresmut[50][50];
	for(int mut1=0;mut1<50;mut1++)
	{	for(int mut2=0;mut2<50;mut2++)
		{	tabresmut[mut1][mut2]=0;
			nbtimetabresmut[mut1][mut2]=1;
		};
	};
	int mostrelatinscore=0;
	int howcalculmostrelaetd=0;
	int howcalculmostrelaetd2=0;
	int addoneandzero=1;
	int localsearch=1;
	int randmuthap=0;
	unsigned char relatsuperpose[MAXCLOSERELAT][MAXCLOSERELAT];
	int64_t relatpihatchr[MAXCLOSERELAT][23][2];
	int seuilscore=1;
	int powersegment=1;		
	for(int relat=0;relat<MAXCLOSERELAT;relat++)// if (ID<10)
	{	for(int  chrtemp1=chr1;chrtemp1<chr2+1;chrtemp1++)		
		{	relatpihatchr[relat][chrtemp1][0]=0;relatpihatchr[relat][chrtemp1][1]=0;
		};
		for(int relat2=relat+1;relat2<MAXCLOSERELAT;relat2++)// if (ID<10)
		{	relatsuperpose[relat][relat2]=127;
			relatsuperpose[relat2][relat]=0;
		};
	};
	int64_t segnum=0;
	unsigned char relatbreak[MAXBREAK][MAXCLOSERELAT];
	for(int relat=0;relat<nbrelatpihat;relat++)// if (ID<10)
	{	for(int  break1=0;break1<nbbreak;break1++)		
		{	relatbreak[break1][relat]=0;
		};
	};
	for(int relat=0;relat<nbrelatpihat;relat++)// if (ID<10)
	{	int highestnbhet=0;
		int highestnbhetbreak=0;
		for(int  break1=0;break1<nbbreak;break1++)		
		{	if (highestnbhet<relatbreak[break1][relat]) 
			{	highestnbhet=relatbreak[break1][relat];
				highestnbhetbreak=break1;
			};
			if (relatbreak[break1][relat]) relatbreak[break1][relat]=1;
		};
		for(int  break1=0;break1<nbbreak;break1++)		
		{	relatbreak[break1][relat]=(highestnbhetbreak==break1);
			
		};
	};	
	for(int relat=0;relat<nbrelatpihat;relat++)// if (ID<10)
	{	printf("%d ",relat);
		for(int  break1=0;break1<nbbreak;break1++)		
		{	printf("%d ",relatbreak[break1][relat]);
		};
		
	};
	int relattochr[MAXCLOSERELAT];
	for(int relat=0;relat<nbrelatpihat;relat++)// if (ID<10)
	{	int bestchr=1;
		for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
		{	int hapmost=(relatpihatchr[relat][chrtemp1][0]>relatpihatchr[relat][chrtemp1][1])?0:1;
			int pihatnorm=relatpihatchr[relat][chrtemp1][hapmost]*nbsnpperchr[1]/nbsnpperchr[chrtemp1];
			if (chrtemp1==1 || relattochr[relat]<pihatnorm)  
			{	relattochr[relat]=pihatnorm;
				bestchr=chrtemp1+22*hapmost;
			};
		};
		relattochr[relat]=bestchr;
	};
	for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
	{	for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)		
		{	genomeoffpss[0][snp][chrtemp1]=(*((genomes[chrtemp1]+(unsigned long long) ID*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			genomeoffpss[1][snp][chrtemp1]=(*((genomes[chrtemp1]+(unsigned long long) IDp1loop*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			genomeoffpss[2][snp][chrtemp1]=(*((genomes[chrtemp1]+(unsigned long long) IDp2loop*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
		};	
	};
	int maxpihat=0;
	int64_t nbseghap[2][23][26229];
	int64_t nbsegoverchr[23][2];
	int ratio=2;	
	printf("maxpihat %d\n",maxpihat);
	FILE * endfile;
	uint64_t averageofaverage[23];
	char number[100];
	char pathfile[300]; 								
	uint64_t sumsumnincrementplus=0;
	uint64_t sumsumnincrementmoins=0;
	sumsumnincrementplus=0;
	sumsumnincrementmoins=0;
	float plusmoins=1.0*sumsumnincrementplus/sumsumnincrementmoins;
	float moinsplus=1.0*sumsumnincrementmoins/sumsumnincrementplus;
	printf("%f %f\n",plusmoins,moinsplus);
	sumsumnincrementplus=0;
	sumsumnincrementmoins=0;
	int64_t sumallpihatpos=0;
	int64_t sumallpihatneg=0;
	int nbpos=0;
	int nbneg=0;
	int64_t minpihat=0;
	int minpihatID1=0;
	int minpihatID2=0;
	unsigned char indivEA[NBINDIVEA][MAXCLOSERELAT];
	unsigned char bestindiv[MAXCLOSERELAT];
	for(int relat=0;relat<nbrelatpihat;relat++)// if (ID<10)
	{	bestindiv[relat]=0;
	};
	int64_t scoreEA[NBINDIVEA];	
	int64_t bestscore=-1000000;
	int64_t nbsumerrorbestscore=0;
	int64_t bestscoreofbestfitness=0;
	int64_t bestphasing=0;
	int same;
	int gen=gentostart;	
	int lastgen=100;	
	int nuincreasenewtop=50;
	int64_t averageparentold5=-1000000;
	double averageparent[MAXGEN];
	int print=1;
	int lastgenimprove=gen;
	int genintiindiv=gen;
	int64_t bestscoresinceinit=bestscore;
	int64_t bestscorefitness=-100000;
	int64_t bestscoreclustering=-100000;
	
	do
	{	int64_t scorerelattab[NBINDIVEA][MAXCLOSERELAT];
		int64_t scorerelattaballchr[NBINDIVEA][MAXCLOSERELAT][MAXBREAK];
		double scorehap[NBINDIVEA][MAXBREAK];
		int bestmut[NBINDIVEA];
	 	for(int IDEA=0;IDEA<NBINDIVEA;IDEA++)// if (ID<10)
		{	int groumdtruthdone=0;
			int bestwellphased=0;
			int bestunwellphased=0;
			int64_t score1=0;
			int64_t score2=0;
			int64_t score1tab[MAXBREAK];
			int64_t score2tab[MAXBREAK];			
			for(int  chrtemp1=0;chrtemp1<MAXBREAK;chrtemp1++)		
			{	score1tab[chrtemp1]=0;
				score2tab[chrtemp1]=0;
			}
				
			int nbonebalanced=0;
		/*	for(int relat=0;relat<nbrelatpihat;relat++)// if (ID<10)
			{	if (indivEA[IDEA][relat]) nbone++;
			};	*/
			int64_t score1averageneg=0;
			int score1nbneg=1;
			int64_t score2averageneg=0;
			int score2nbneg=1;
			
			bestmut[IDEA]=-1;
		
			int chagemade=1;
			int64_t iter=0;
			int iterbreak=0;
			
			int64_t maxscorerelattab=0;
			int64_t plus[MAXCLOSERELAT]={0};
			int64_t minus[MAXCLOSERELAT]={0};
			maxscorerelattab=0;
			int64_t maxscorerelattabID=0;
			
			
			int64_t score1temp=0;
			 int64_t  score2temp=0;
			printf("%d %d %d %d\n",score1,score2,score1nbneg,score1nbneg);
			int64_t score;
			int64_t averagescore1=1;//((nbonebalanced+1)!=0)?score1/(nbonebalanced+1):score1;
			int64_t averagescore2=1;//((1+nbrelatpihat-nbonebalanced)!=0)?score2/(1+nbrelatpihat-nbonebalanced):score2;
			if (nbonebalanced*2>nbrelatpihat) nbonebalanced=1+nbonebalanced*2-nbrelatpihat; else nbonebalanced=1+nbrelatpihat-nbonebalanced*2;
			if (bestscoresinceinit<score) 
			{//	printf("improve %d %d %lld %lld\n",lastgenimprove,gen,bestscoresinceinit,score);
				lastgenimprove=gen;
				bestscoresinceinit=score;
			};
			int64_t highest=score;
			scoreEA[IDEA]=score/bestcor;
			
			{	//bestscore=highest;
				bestscoreclustering=score;
				if (bestscoreofbestfitness==0) bestscoreofbestfitness=1;
				char number[100];
				char pathfile[300]; 
		
				int64_t scorechr1[23];
				int64_t scorechr2[23];
				int64_t scorechr3[23];
				int64_t scorechr4[23];
				int keep[23][4];
				int64_t nbonefromhap[2][23][26229];
				int64_t nbzerofromhap[2][23][26229];
				
				int64_t nbonequalified[2];
				int64_t nbzeroqualified[2];
				nbonequalified[0]=0;
				nbonequalified[1]=0;
				nbzeroqualified[0]=0;
				nbzeroqualified[1]=0;
				int64_t sumdiff[23];
				int64_t nbcontracditperchr[23];
				
				for(int nbkickout=1;nbkickout<23;nbkickout++)
				{	keep[nbkickout][ratio]=1;
				};
				{	nbsumerrorbestscore=1;
					for(int relat=0;relat<nbrelatpihat;relat++)// if (ID<10)
					{	bestindiv[relat]=indivEA[IDEA][relat];
					};	
					int phase1[23][1][1];
					int guessphase1[23][1][1];
					int phase1guessphase1[23][1][1];
					int phase1guessphase2[23][1][1];
					int phase2[23][1][1];
					int guessphase2[23][1][1];
					int phase2guessphase2[23][1][1];
					int phase2guessphase1[23][1][1];
					int phaseguessright[23][1][1];
					int phaseguesswrong[23][1][1];
					int phasenotguessed[23][1][1];
					int nbsum=1;
				
					int wrong[23];
					int right[23];
					unsigned int x, y, width, height;
							
					for(int prod=0;prod<1;prod++)
					{	for(int chr=0;chr<23;chr++)
						{	for(int ratio=0;ratio<1;ratio++)
							{	phaseguessright[chr][prod][ratio]=0;
								phaseguesswrong[chr][prod][ratio]=0;
								phasenotguessed[chr][prod][ratio]=0;
								phase1[chr][prod][ratio]=0;
								guessphase1[chr][prod][ratio]=0;
								phase1guessphase1[chr][prod][ratio]=0;
								phase1guessphase2[chr][prod][ratio]=0;
								phase2[chr][prod][ratio]=0;
								guessphase2[chr][prod][ratio]=0;
								phase2guessphase2[chr][prod][ratio]=0;
								phase2guessphase1[chr][prod][ratio]=0;
							};	
						};
					};
					int lenshowchr=height*0.9;
					int widhtshowchr=5;
					int gapshowchry=63;
					int widthshowchry=9;
					int gapshowchrx=5;
					int64_t sumdifftot=0;
					int phaseparent1=-1;
					int phaseparent1tap[23][26229];
					int phaseguesedparent1tap[23][26229];
					for(int  chrtemp1=1;chrtemp1<22+1;chrtemp1++)		
					{	for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)
						{	phaseparent1tap[chrtemp1][snp]=0;
							phaseguesedparent1tap[chrtemp1][snp]=0;
						};
					};
					int phaseunknow=0;	
					int bestpur=0;
					double dataperchr[23][110];		
					int nbdivision=2;
					int powerpihat=0;
					int exponentcorp=0;
					double sumoftabcor=0;
					int firstexponent;
					int thirdindice;
					double pihatagainstallchrMPphaseerror[23][MAXPOP][200][2];//={0};
					int purcentage;
					int powerpihatDEGREE4=0;
					for(int size=0;size<51;size++)
					{	for(int chr=1;chr<23;chr++)
						{	for(int window=0;window<1;window++)
							{	chrdivider[size][chr][window].segment=0;
							};
						};
					};

					//MIN.CM = 25 
					chrdivider[25][1][0].end=2203;
					chrdivider[25][1][1].start=2203;
					chrdivider[25][1][1].end=7136;
					chrdivider[25][1][2].start=7136;
					chrdivider[25][1][2].end=14606;
					chrdivider[25][1][3].start=14606;
					chrdivider[25][1][3].end=18289;
					chrdivider[25][1][4].start=18289;
					chrdivider[25][1][4].end=nbsnpperchr[1];
					nbchrdivider[25][1]=5;
					nbchrdivider[25][0]=nbchrdivider[25][1];
					chrdivider[25][2][0].end=2097;
					chrdivider[25][2][1].start=2097;
					chrdivider[25][2][1].end=6555;
					chrdivider[25][2][2].start=6555;
					chrdivider[25][2][2].end=12181;
					chrdivider[25][2][3].start=12181;
					chrdivider[25][2][3].end=18340;
					chrdivider[25][2][4].start=18340;
					chrdivider[25][2][4].end=nbsnpperchr[2];
					nbchrdivider[25][2]=5;
					if (nbchrdivider[25][2]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][2];
					chrdivider[25][3][0].end=4737;
					chrdivider[25][3][1].start=4737;
					chrdivider[25][3][1].end=9553;
					chrdivider[25][3][2].start=9553;
					chrdivider[25][3][2].end=14838;
					chrdivider[25][3][3].start=14838;
					chrdivider[25][3][3].end=18519;
					chrdivider[25][3][4].start=18519;
					chrdivider[25][3][4].end=nbsnpperchr[3];
					nbchrdivider[25][3]=5;
					if (nbchrdivider[25][3]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][3];
					chrdivider[25][4][0].end=3754;
					chrdivider[25][4][1].start=3754;
					chrdivider[25][4][1].end=6615;
					chrdivider[25][4][2].start=6615;
					chrdivider[25][4][2].end=9749;
					chrdivider[25][4][3].start=9749;
					chrdivider[25][4][3].end=16124;
					chrdivider[25][4][4].start=16124;
					chrdivider[25][4][4].end=nbsnpperchr[4];
					nbchrdivider[25][4]=5;
					if (nbchrdivider[25][4]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][4];
					chrdivider[25][5][0].end=6946;
					chrdivider[25][5][1].start=6946;
					chrdivider[25][5][1].end=11357;
					chrdivider[25][5][2].start=11357;
					chrdivider[25][5][2].end=15071;
					chrdivider[25][5][3].start=15071;
					chrdivider[25][5][3].end=nbsnpperchr[5];
					nbchrdivider[25][5]=4;
					if (nbchrdivider[25][5]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][5];
					chrdivider[25][6][0].end=3976;
					chrdivider[25][6][1].start=3976;
					chrdivider[25][6][1].end=11907;
					chrdivider[25][6][2].start=11907;
					chrdivider[25][6][2].end=nbsnpperchr[6];
					nbchrdivider[25][6]=3;
					if (nbchrdivider[25][6]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][6];
					chrdivider[25][7][0].end=3082;
					chrdivider[25][7][1].start=3082;
					chrdivider[25][7][1].end=7188;
					chrdivider[25][7][2].start=7188;
					chrdivider[25][7][2].end=12684;
					chrdivider[25][7][3].start=12684;
					chrdivider[25][7][3].end=nbsnpperchr[7];
					nbchrdivider[25][7]=4;
					if (nbchrdivider[25][7]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][7];
					chrdivider[25][8][0].end=7271;
					chrdivider[25][8][1].start=7271;
					chrdivider[25][8][1].end=12465;
					chrdivider[25][8][2].start=12465;
					chrdivider[25][8][2].end=nbsnpperchr[8];
					nbchrdivider[25][8]=3;
					if (nbchrdivider[25][8]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][8];
					chrdivider[25][9][0].end=3384;
					chrdivider[25][9][1].start=3384;
					chrdivider[25][9][1].end=7020;
					chrdivider[25][9][2].start=7020;
					chrdivider[25][9][2].end=9803;
					chrdivider[25][9][3].start=9803;
					chrdivider[25][9][3].end=nbsnpperchr[9];
					nbchrdivider[25][9]=4;
					if (nbchrdivider[25][9]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][9];
					chrdivider[25][10][0].end=2407;
					chrdivider[25][10][1].start=2407;
					chrdivider[25][10][1].end=8424;
					chrdivider[25][10][2].start=8424;
					chrdivider[25][10][2].end=12288;
					chrdivider[25][10][3].start=12288;
					chrdivider[25][10][3].end=nbsnpperchr[10];
					nbchrdivider[25][10]=4;
					if (nbchrdivider[25][10]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][10];
					chrdivider[25][11][0].end=2953;
					chrdivider[25][11][1].start=2953;
					chrdivider[25][11][1].end=5344;
					chrdivider[25][11][2].start=5344;
					chrdivider[25][11][2].end=nbsnpperchr[11];
					nbchrdivider[25][11]=3;
					if (nbchrdivider[25][11]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][11];
					chrdivider[25][12][0].end=5937;
					chrdivider[25][12][1].start=5937;
					chrdivider[25][12][1].end=9107;
					chrdivider[25][12][2].start=9107;
					chrdivider[25][12][2].end=14609;
					chrdivider[25][12][3].start=14609;
					chrdivider[25][12][3].end=nbsnpperchr[12];
					nbchrdivider[25][12]=4;
					if (nbchrdivider[25][12]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][12];
					chrdivider[25][13][0].end=2524;
					chrdivider[25][13][1].start=2524;
					chrdivider[25][13][1].end=6229;
					chrdivider[25][13][2].start=6229;
					chrdivider[25][13][2].end=nbsnpperchr[13];
					nbchrdivider[25][13]=3;
					if (nbchrdivider[25][13]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][13];
					chrdivider[25][14][0].end=3228;
					chrdivider[25][14][1].start=3228;
					chrdivider[25][14][1].end=6437;
					chrdivider[25][14][2].start=6437;
					chrdivider[25][14][2].end=nbsnpperchr[14];
					nbchrdivider[25][14]=3;
					if (nbchrdivider[25][14]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][14];
					chrdivider[25][15][0].end=3288;
					chrdivider[25][15][1].start=3288;
					chrdivider[25][15][1].end=6043;
					chrdivider[25][15][2].start=6043;
					chrdivider[25][15][2].end=8861;
					chrdivider[25][15][3].start=8861;
					chrdivider[25][15][3].end=nbsnpperchr[15];
					nbchrdivider[25][15]=4;
					if (nbchrdivider[25][15]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][15];
					chrdivider[25][16][0].end=3860;
					chrdivider[25][16][1].start=3860;
					chrdivider[25][16][1].end=nbsnpperchr[16];
					nbchrdivider[25][16]=2;
					if (nbchrdivider[25][16]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][16];
					chrdivider[25][17][0].end=2196;
					chrdivider[25][17][1].start=2196;
					chrdivider[25][17][1].end=4676;
					chrdivider[25][17][2].start=4676;
					chrdivider[25][17][2].end=8643;
					chrdivider[25][17][3].start=8643;
					chrdivider[25][17][3].end=nbsnpperchr[17];
					nbchrdivider[25][17]=4;
					if (nbchrdivider[25][17]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][17];
					chrdivider[25][18][0].end=3734;
					chrdivider[25][18][1].start=3734;
					chrdivider[25][18][1].end=7591;
					chrdivider[25][18][2].start=7591;
					chrdivider[25][18][2].end=nbsnpperchr[18];
					nbchrdivider[25][18]=3;
					if (nbchrdivider[25][18]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][18];
					chrdivider[25][19][0].end=2982;
					chrdivider[25][19][1].start=2982;
					chrdivider[25][19][1].end=6035;
					chrdivider[25][19][2].start=6035;
					chrdivider[25][19][2].end=nbsnpperchr[19];
					nbchrdivider[25][19]=3;
					if (nbchrdivider[25][19]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][19];
					chrdivider[25][20][0].end=2054;
					chrdivider[25][20][1].start=2054;
					chrdivider[25][20][1].end=5072;
					chrdivider[25][20][2].start=5072;
					chrdivider[25][20][2].end=nbsnpperchr[20];
					nbchrdivider[25][20]=3;
					if (nbchrdivider[25][20]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][20];
					chrdivider[25][21][0].end=3025;
					chrdivider[25][21][1].start=3025;
					chrdivider[25][21][1].end=nbsnpperchr[21];
					nbchrdivider[25][21]=2;
					if (nbchrdivider[25][21]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][21];
					chrdivider[25][22][0].end=1615;
					chrdivider[25][22][1].start=1615;
					chrdivider[25][22][1].end=nbsnpperchr[22];
					nbchrdivider[25][22]=2;
					if (nbchrdivider[25][22]>nbchrdivider[25][0]) nbchrdivider[25][0]=nbchrdivider[25][22];



					int nbdivisionDEGREE4=0; //for(int powerpihatDEGREE4=0;powerpihatDEGREE4<3;powerpihatDEGREE4++)// for(int exponentcorp=0;exponentcorp<3;exponentcorp++) 
				//	for(int thirdindiceloop=0;thirdindiceloop<2;thirdindiceloop++)			
					{	char segwithav[23][26229];
						double seuilpihatcorrection=0.022;
						printf("%d %f\n",placefirttoconsider,pihatagainstall[placefirttoconsider]);
						for(int  relattocompare=0;relattocompare<20;relattocompare++)	
						{	int countsnp=0;
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbsnpperchr[chrtemp1];chrdividerrun++)
								{	segwithav[chrtemp1][chrdividerrun]=0; 
								};
							};
							
							int IDrelattocompare=bestpihatagainstallID[placefirttoconsider+relattocompare];
							if (pihatagainstall[IDrelattocompare]>seuilpihatcorrection	)
							{	
								#pragma omp parallel for 	
								for(int  chrtemp1=22;chrtemp1>0;chrtemp1--)		
								{	printf("Start correctin chr %d\n",chrtemp1);
									int nbindivmatchseg[MAXPOP][4];
									for(int64_t relat=0;relat<(int64_t) MAXPOP;relat++)// if (pihatagainstall[relat]<seuilpihat[0] )
									{	nbindivmatchseg[relat][0]=0;
										nbindivmatchseg[relat][1]=0;
										nbindivmatchseg[relat][2]=0;
										nbindivmatchseg[relat][3]=0;
									};
									int nbsegmentofthislength[NSNPPERCHR];
									for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)
									{	nbsegmentofthislength[snp]=0;
									};
									int typesegment=-1;
									int lasttypesegment=-1;
									int endlastsegment[4]={-1};
									double ratetoconsiderseg=0.0001;
									double ratetoconsidersegPE=ratetoconsiderseg*2; 
									ratetoconsiderseg=0.0002;   
									ratetoconsidersegPE=2*ratetoconsiderseg; 
									int phaseerrorpossible[4]={0};
									int breaknubercm=25;
									for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)
									{	//printf("%d %d %d %d",chrdividerrun,lastmatchend	,chrdividerrun,lastmatchtype); 
										int snpvalue0=genomeoffpss[0][snp][chrtemp1];		
										int nbsnpperchrby4=(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0));
										int snpmod4=(((snp%4)*2));
										int snpdiv4=snp/4;
										//#pragma omp parallel for 			
										for(int64_t relat=0;relat<(int64_t) MAXPOP;relat++) if (relat!=ID && pihatagainstall[relat]<seuilpihat[0] )
										{	int snpvalue1=(*((genomes[chrtemp1]+(unsigned long long) relat*nbsnpperchrby4 )+snpdiv4)>>snpmod4)&3;
											if ((snpvalue0&1)==(snpvalue1&1))
											{	if (nbindivmatchseg[relat][0]<10) nbindivmatchseg[relat][0]++; else 
												nbsegmentofthislength[++nbindivmatchseg[relat][0]]++;
											} else 
											{	if (nbindivmatchseg[relat][0]<10) nbindivmatchseg[relat][0]=0; else 
												{	for(int lenght=10;lenght<nbindivmatchseg[relat][0]+1;lenght++) nbsegmentofthislength[lenght]--;
													nbindivmatchseg[relat][0]=0;
												}
											};
											if ((snpvalue0&1)==(snpvalue1>>1))
											{	if (nbindivmatchseg[relat][1]<10) nbindivmatchseg[relat][1]++;
												else nbsegmentofthislength[++nbindivmatchseg[relat][1]]++;
											} else 
											{	if (nbindivmatchseg[relat][1]<10) nbindivmatchseg[relat][1]=0; else 
												{	for(int lenght=1;lenght<nbindivmatchseg[relat][1]+1;lenght++) nbsegmentofthislength[lenght]--;
													nbindivmatchseg[relat][1]=0;
												};
											};
											if ((snpvalue0>>1)==(snpvalue1&1))
											{	if (nbindivmatchseg[relat][2]<10) nbindivmatchseg[relat][2]++;
												else nbsegmentofthislength[++nbindivmatchseg[relat][2]]++;
											} else 
											{	if (nbindivmatchseg[relat][2]<10) nbindivmatchseg[relat][2]=0;
												else 
												{	for(int lenght=1;lenght<nbindivmatchseg[relat][2]+1;lenght++) nbsegmentofthislength[lenght]--;
													nbindivmatchseg[relat][2]=0;
												}
											};
											if ((snpvalue0>>1)==(snpvalue1>>1))
											{	if (nbindivmatchseg[relat][3]<10) nbindivmatchseg[relat][3]++;
												else nbsegmentofthislength[++nbindivmatchseg[relat][3]]++;
											} else 
											{	if (nbindivmatchseg[relat][3]<10) nbindivmatchseg[relat][3]=0;	
												else 
												{	for(int lenght=1;lenght<nbindivmatchseg[relat][3]+1;lenght++) nbsegmentofthislength[lenght]--;
													nbindivmatchseg[relat][3]=0;	
												};
											};
										};
										
						//				for(int snprun=0;snprun<100;snprun++)
						//				{	printf("%d %d\n",snprun,nbsegmentofthislength[snprun]);
						//				};
										
										int relat=IDrelattocompare;
										if (nbindivmatchseg[relat][0]==0) phaseerrorpossible[0]=0;
										if (nbindivmatchseg[relat][1]==0) phaseerrorpossible[1]=0;
										if (nbindivmatchseg[relat][2]==0) phaseerrorpossible[2]=0;
										if (nbindivmatchseg[relat][3]==0) phaseerrorpossible[3]=0;
										
										int snpvalue1=(*((genomes[chrtemp1]+(unsigned long long) relat*nbsnpperchrby4 )+snpdiv4)>>snpmod4)&3;
										int lenght;
										if (nbindivmatchseg[relat][0]>10 && 
												nbsegmentofthislength[nbindivmatchseg[relat][0]]<(phaseerrorpossible[0]?MAXPOP*ratetoconsidersegPE:MAXPOP*ratetoconsiderseg) && 
												nbindivmatchseg[relat][0]>nbindivmatchseg[relat][1] &&
												nbindivmatchseg[relat][0]>nbindivmatchseg[relat][2] &&
												nbindivmatchseg[relat][0]>nbindivmatchseg[relat][3] 
												)
										{	printf("chr %d segment type 1 at snp %d of lenght %d as only %d indivs has it\n",chrtemp1,snp,nbindivmatchseg[relat][0], nbsegmentofthislength[nbindivmatchseg[relat][0]]);
											lenght=nbindivmatchseg[relat][0];
											typesegment=1;
											for(int spnrun=snp-nbindivmatchseg[relat][0];spnrun<snp+1;spnrun++) segwithav[chrtemp1][spnrun]=1;
											endlastsegment[0]=snp;	
										//	chrdivider[breaknubercm][chrtemp1][window].segment
										} else if (nbindivmatchseg[relat][1]>10 && nbsegmentofthislength[nbindivmatchseg[relat][1]]<(phaseerrorpossible[1]?MAXPOP*ratetoconsidersegPE:MAXPOP*ratetoconsiderseg) &&
													nbindivmatchseg[relat][1]>nbindivmatchseg[relat][2] &&
													nbindivmatchseg[relat][1]>nbindivmatchseg[relat][3] 
												)
										{	printf("chr %d segment type 3 at snp %d of lenght %d as only %d indivs has it\n",chrtemp1,snp,nbindivmatchseg[relat][1], nbsegmentofthislength[nbindivmatchseg[relat][1]]);
											lenght=nbindivmatchseg[relat][1];
											typesegment=3;
											endlastsegment[1]=snp;
											for(int spnrun=snp-nbindivmatchseg[relat][1];spnrun<snp+1;spnrun++) segwithav[chrtemp1][spnrun]=1;
										} else if (nbindivmatchseg[relat][2]>10 && nbsegmentofthislength[nbindivmatchseg[relat][2]]<(phaseerrorpossible[2]?MAXPOP*ratetoconsidersegPE:MAXPOP*ratetoconsiderseg) &&
													nbindivmatchseg[relat][2]>nbindivmatchseg[relat][3] 
												)
										{	printf("chr %d segment type 2 at snp %d of lenght %d as only %d indivs has it\n",chrtemp1,snp,nbindivmatchseg[relat][2], nbsegmentofthislength[nbindivmatchseg[relat][2]]);
											lenght=nbindivmatchseg[relat][2];
											typesegment=2;
											endlastsegment[2]=snp;
											for(int spnrun=snp-nbindivmatchseg[relat][2];spnrun<snp+1;spnrun++) segwithav[chrtemp1][spnrun]=2;
											
										} else if (nbindivmatchseg[relat][3]>10 && nbsegmentofthislength[nbindivmatchseg[relat][3]]<(phaseerrorpossible[3]?MAXPOP*ratetoconsidersegPE:MAXPOP*ratetoconsiderseg)
												)
										{	printf("chr %d segment type 4 at snp %d of lenght %d as only %d indivs has it\n",chrtemp1,snp,nbindivmatchseg[relat][3], nbsegmentofthislength[nbindivmatchseg[relat][3]]);
											lenght=nbindivmatchseg[relat][3];
											typesegment=4;
											endlastsegment[3]=snp;
											for(int spnrun=snp-nbindivmatchseg[relat][3];spnrun<snp+1;spnrun++) segwithav[chrtemp1][spnrun]=2;
										};
										if (endlastsegment[0]==snp-1 || endlastsegment[1]==snp-1 || endlastsegment[2]==snp-1 || endlastsegment[3]==snp-1)
										{	if ((snpvalue1&1)!=(snpvalue1>>1))
											{	phaseerrorpossible[0]=1;
												phaseerrorpossible[1]=1;
												phaseerrorpossible[2]=1;
												phaseerrorpossible[3]=1;
										
											}
										}
										if (typesegment>-1 && lasttypesegment>-1 && (typesegment&1)!=(lasttypesegment&1))
										{	int end;
											if (typesegment==1) end=endlastsegment[0];
											if (typesegment==2) end=endlastsegment[2];
											if (typesegment==3) end=endlastsegment[1];
											if (typesegment==4) end=endlastsegment[3];
											printf("change from %d to %d\n",(end+snp-lenght)/2,nbsnpperchr[chrtemp1]);
										//	nbphasecorrect++;
											for(int snprun=(snp-lenght);snprun<nbsnpperchr[chrtemp1];snprun++)		
											{//	nbsnpcorrect++;
												genomeoffpss[0][snprun][chrtemp1]=((genomeoffpss[0][snprun][chrtemp1]&1)<<1)+(genomeoffpss[0][snprun][chrtemp1]>>1);
											};
											snp=end;
											for(int snprun=end;snprun<nbsnpperchr[chrtemp1];snprun++)		
											{	segwithav[chrtemp1][snprun]=0;
											};
											
											
											typesegment=-1;
											//lasttypesegment=-1;
										//	endlastsegment=-1;
											for(int64_t relat=0;relat<(int64_t) MAXPOP;relat++) if (pihatagainstall[relat]<seuilpihat[0] )
											{	nbindivmatchseg[relat][0]=0;
												nbindivmatchseg[relat][1]=0;
												nbindivmatchseg[relat][2]=0;
												nbindivmatchseg[relat][3]=0;
											};
											for(int snprun=0;snprun<nbsnpperchr[chrtemp1];snprun++)
											{	nbsegmentofthislength[snprun]=0;
											};
										//	snp=snp-lenght;
										} else if (typesegment>-1)
										{	lasttypesegment=typesegment;
										};									
										
									
									};		
								};
							};
						};
						int breaknubercm=25;
						int breaknubercmloop=25;
						{	for(int prod=0;prod<1;prod++)
							{	for(int chr=0;chr<23;chr++)
								{	for(int ratio=0;ratio<1;ratio++)
									{	phaseguessright[chr][prod][ratio]=0;
										phaseguesswrong[chr][prod][ratio]=0;
										phasenotguessed[chr][prod][ratio]=0;
										phase1[chr][prod][ratio]=0;
										guessphase1[chr][prod][ratio]=0;
										phase1guessphase1[chr][prod][ratio]=0;
										phase1guessphase2[chr][prod][ratio]=0;
										phase2[chr][prod][ratio]=0;
										guessphase2[chr][prod][ratio]=0;
										phase2guessphase2[chr][prod][ratio]=0;
										phase2guessphase1[chr][prod][ratio]=0;
									};	
								};
							};
							
							for(int  chrtemp1=1;chrtemp1<22+1;chrtemp1++)		
							{	for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)
								{	phaseparent1tap[chrtemp1][snp]=0;
									phaseguesedparent1tap[chrtemp1][snp]=0;
								};
							};							
							
						//	if (breaknubercmloop==23) breaknubercm=24;
						//	else if (breaknubercmloop==24) breaknubercm=25;
						//	else if (breaknubercmloop==25) breaknubercm=24;
							char segwithav[23][26229];
							int countsnp=0;
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbsnpperchr[chrtemp1];chrdividerrun++)
								{	segwithav[chrtemp1][chrdividerrun]=0; 
								};
							};
							for(int  chrtemp1=1;chrtemp1<22+1;chrtemp1++)		
							{	phaseguessright[chrtemp1][0][0]=0;
								phaseguesswrong[chrtemp1][0][0]=0;	
								wrong[chrtemp1]=0;
								right[chrtemp1]=0;
								int phaseparent1=-1;
								int phaseparent1temp=-1;
								int nbsnperror=0;
								int64_t nbonesnp=1;//rand();
								int64_t nbzerosnp=0;//rand();
							//	printf("%d %d %d\n",6,phaseguessright[6][0][0],phaseguesswrong[6][0][0]);
								for(int snp=0;snp<nbsnpperchr[chrtemp1];snp++)
								{	int marker=genomeoffpss[0][snp][chrtemp1];
									int markerp1=genomeoffpss[1][snp][chrtemp1];
								
									int snperrofind=0;
									if (marker==0 || marker==3)
									{	if (markerp1==3-marker) snperrofind=1;
									
									};
									if (snperrofind==0)
									{	int phase=-1;
										if (marker!=0 && marker!=3)
										{	if (markerp1==0) phase=marker;
											else if (markerp1==3) phase=3-marker;	
										};
									//	 printf("snp %d: %d %d %d\n",snp,phaseparent1temp,phaseparent1,phase);
										if (phase>-1)
										{	if (phaseparent1==-1) 
											{	phaseparent1=phase;
												phaseparent1temp=phase;
												printf("phase %d detected at snp %d\n",phase,snp);
											} else if (phaseparent1!=phase) 
											{	if (phaseparent1temp==phase)
												{	printf("phase change at snp %d: %d %d %d\n",snp,phaseparent1temp,phaseparent1,phase);
													phaseparent1=phase;
													phaseparent1temp=phase;
												} else 
												{	phaseparent1temp=phase;
												};
											} else 
											{	phaseparent1temp=phaseparent1;
											}
										};
									} else 
									{	nbsnperror++;
									};
									
									phaseparent1tap[chrtemp1][snp]=phaseparent1;								
									if (phaseparent1>-1)
									{	if (nbonesnp>1*nbzerosnp && phaseparent1==1) // && nbone[snp]>sum+((prod+1)*nbzero[snp])
										{	phase1[ratio][0][0]++;
											if (keep[chrtemp1][ratio]) guessphase1[ratio][0][0]++; else guessphase2[ratio][0][0]++;
											if (keep[chrtemp1][ratio]) phase1guessphase1[ratio][0][0]++; else  phase1guessphase2[ratio][0][0]++;
											for(int ratio=2;ratio<3;ratio++)
											{	if (keep[chrtemp1][ratio]) 
												{	phaseguessright[chrtemp1][0][0]++;
													phaseguesedparent1tap[chrtemp1][snp]=1;
												}
												else 
												{	phaseguesswrong[chrtemp1][0][0]++;//&& 2*nbone>3*nbzero+5 
													phaseguesedparent1tap[chrtemp1][snp]=0;
												};
											};											
											right[chrtemp1]++;
									//		printf(" R");
										}
										else if (nbonesnp>1*nbzerosnp && phaseparent1==2)//&& nbone[snp]>sum+((prod+1)*nbzero[snp]) 
										{	phase2[ratio][0][0]++;
											if (keep[chrtemp1][ratio]) guessphase1[ratio][0][0]++; else guessphase2[ratio][0][0]++;
											if (keep[chrtemp1][ratio]) phase2guessphase1[ratio][0][0]++; else phase2guessphase2[ratio][0][0]++;
											for(int ratio=2;ratio<3;ratio++)
											{	if (keep[chrtemp1][ratio]) 
												{	phaseguesswrong[chrtemp1][0][0]++;//&& 2*nbone>3*nbzero+5 
													phaseguesedparent1tap[chrtemp1][snp]=1;
												}
												else 
												{	phaseguessright[chrtemp1][0][0]++;
													phaseguesedparent1tap[chrtemp1][snp]=0;
												};
											};
											wrong[chrtemp1]++;
									//		printf(" W");
										}
										else if (nbonesnp*1<nbzerosnp && phaseparent1==2)// && sum+((0+1)*nbone[snp])<nbzero[snp]
										{	phase2[ratio][0][0]++;
											if (keep[chrtemp1][ratio]) guessphase2[ratio][0][0]++; else guessphase1[ratio][0][0]++;
											if (keep[chrtemp1][ratio]) phase2guessphase2[ratio][0][0]++; else phase2guessphase1[ratio][0][0]++;
											for(int ratio=2;ratio<3;ratio++)
											{	if (keep[chrtemp1][ratio]) 
												{	phaseguessright[chrtemp1][0][0]++;//&& 3*nbone+5<nbzero*2
													phaseguesedparent1tap[chrtemp1][snp]=0;
												}	
												else 
												{	phaseguesswrong[chrtemp1][0][0]++;
													phaseguesedparent1tap[chrtemp1][snp]=1;
												}
											};
											right[chrtemp1]++;
									//		printf(" r");
										}
										else if (nbonesnp*1<nbzerosnp && phaseparent1==1) //&& sum+((prod+1)*nbone[snp])<nbzero[snp]
										{	phase1[ratio][0][0]++;
											if (keep[chrtemp1][ratio]) guessphase2[ratio][0][0]++; else guessphase1[ratio][0][0]++;
											if (keep[chrtemp1][ratio]) phase1guessphase2[ratio][0][0]++; else phase1guessphase1[ratio][0][0]++;
											for(int ratio=2;ratio<3;ratio++)
											{	if (keep[chrtemp1][ratio]) 
												{	phaseguesswrong[chrtemp1][0][0]++;	//&& 3*nbone+5<nbzero*2
													phaseguesedparent1tap[chrtemp1][snp]=0;
												}
												else 
												{	phaseguessright[chrtemp1][0][0]++;
													phaseguesedparent1tap[chrtemp1][snp]=1;
												};	
											};
											wrong[chrtemp1]++;
									//		printf(" w");
										}
										else 
										{	for(int ratio=2;ratio<3;ratio++)
											{	if (keep[chrtemp1][ratio]) phasenotguessed[ratio][0][0]++;
												else phasenotguessed[ratio][0][0]++;
											};
										};			
									} else 
									{	phaseunknow++;
									};
									
									//printf("%d s %d se %d 1 %d 0 %d P %d r %d w %d\n",chrtemp1, snp,nbsnperror,nbone,nbzero,phaseparent1,right[chrtemp1],wrong[chrtemp1]);
								};
								bestpur=bestpur+(wrong[chrtemp1]>right[chrtemp1]?wrong[chrtemp1]:right[chrtemp1]);
							//	printf("chr %d right %d wrong %d nbsnperror %d\n",chrtemp1,right[chrtemp1],wrong[chrtemp1],nbsnperror);
							//	sumnotguessed=sumnotguessed+phasenotguessed;
								//sumright=sumright+(phaseguessright,phaseguesswrong
						//		printf("ID %d chr %d diff %lld diffbysnp %lld\n",IDEA,chrtemp1,sumdiff[chrtemp1],sumdiff[chrtemp1]/nbsnpperchr[chrtemp1]);
								 
								sumdifftot=sumdifftot+sumdiff[chrtemp1];
								printf("%d %d %d\n",chrtemp1,phaseguessright[chrtemp1][0][0],phaseguesswrong[chrtemp1][0][0]);
								dataperchr[chrtemp1][9]=(phaseguessright[chrtemp1][0][0]>phaseguesswrong[chrtemp1][0][0]?phaseguessright[chrtemp1][0][0]:phaseguesswrong[chrtemp1][0][0])/(phaseguessright[chrtemp1][0][0]+phaseguesswrong[chrtemp1][0][0]);
								printf("%d %d %d\n",6,phaseguessright[6][0][0],phaseguesswrong[6][0][0]);
							};
							int64_t purcentage=0;
							ratio=2;
							printf("%d %d %d\n",6,phaseguessright[6][0][0],phaseguesswrong[6][0][0]);
							purcentage=(int64_t) 10000*
													(1+(phaseguessright[ratio][0][0]>phaseguesswrong[ratio][0][0]?
															phaseguessright[ratio][0][0]:
															phaseguesswrong[ratio][0][0]))/ 
													(phaseguesswrong[ratio][0][0]+phaseguessright[ratio][0][0]+1.0);
							printf("%d pourcentage before %d %f %f\n",IDEA,purcentage,(float) bestpur/330005); 
						
							if (purcentage>bestphasing) bestphasing=purcentage;
							
							int nbphaseerror=0;
							double pihatthrehjold[23];
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	pihatthrehjold[chrtemp1]=0.10;//-exponentcorp*0.05 ;
							};
							int nbphaseright=0;
							int nbphasewrong=0;
							int tabnbphaseright[23][25];
							int tabnbphasewrong[23][25];
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	//for(int  snpphaseerror=0;snpphaseerror<nbsnpperchr[chrtemp1];snpphaseerror=snpphaseerror+100)		
								for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	tabnbphaseright[chrtemp1][chrdividerrun]=0;
									tabnbphasewrong[chrtemp1][chrdividerrun]=0;
								};
							};
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	//for(int  snpphaseerror=0;snpphaseerror<nbsnpperchr[chrtemp1];snpphaseerror=snpphaseerror+100)		
								for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	int nbphaseright=0;
									int nbphasewrong=0;
									for(int snp=chrdivider[breaknubercm][chrtemp1][chrdividerrun].start;snp<chrdivider[breaknubercm][chrtemp1][chrdividerrun].end;snp++)		
									{	if (1==phaseparent1tap[chrtemp1][snp]-1) nbphaseright++; else nbphasewrong++;
									};
									printf("%d %d %d %d\n",chrtemp1,chrdividerrun,nbphaseright,nbphasewrong);
									tabnbphaseright[chrtemp1][chrdividerrun]=nbphaseright;
									tabnbphasewrong[chrtemp1][chrdividerrun]=nbphasewrong;
								};	
							}			
							double pihatagainstallchrMP[MAXPOP][2]={0};
							 #pragma omp parallel for
							for(int relat=0;relat<MAXPOP;relat++)// if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
							{	pihatagainstallchrMP[relat][0]=0;
								pihatagainstallchrMP[relat][1]=0;
							}
							
							int sizebin=450;
							
							#pragma omp parallel for					
							for(int relat=0;relat<MAXPOP;relat++)// if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
							{	pihatagainstallchrMPphaseerror[0][relat][0][0]=0;
								pihatagainstallchrMPphaseerror[0][relat][0][1]=0;
								for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][0]*4;chrdividerrun++)
								{	pihatagainstallchrMPphaseerror[0][relat][chrdividerrun][0]=0;
									pihatagainstallchrMPphaseerror[0][relat][chrdividerrun][1]=0;
								};
							}
							int hetsnp=0;
							double exponentinterchr=1;
							for(int relat=0;relat<MAXPOP;relat++)// if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
							{	//printf("%d %f %f\n",relat,pihatagainstall[relat],seuilpihat[powerpihat]);
							}
							for(int relat=0;relat<MAXPOP;relat++)// if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
							{	if (pihatagainstall[relat]>0)  printf("PIHAT %d %f %f\n",relat,pihatagainstall[relat],seuilpihat[0]);
							}
							float firstexponent=2;//2+nbdivisionDEGREE4 ;//4.5+.5*nbdivision;//+exponentcorp*.1; DONE
							int powersnppihat=4;
							int proprotiontokeep=3025;
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	//for(int  snpphaseerror=0;snpphaseerror<nbsnpperchr[chrtemp1];snpphaseerror=snpphaseerror+100)		
								{	hetsnp=0;
								//chrdivider=nbsnpperchr[chrtemp1]/sizebin;
									#pragma omp parallel for 				
									for(int relat=0;relat<MAXPOP;relat++)// if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
									{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1]*4;chrdividerrun++)
										{	pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun][0]=0;
											pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun][1]=0;
										};
									//	pihatagainstallchrMPphaseerror[chrtemp1][1][relat][0]=0;
									};
								//	int  chrtemp2=chrtemp1;
								
									for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
									{	double temp[MAXPOP][4];
										for(int64_t relat=0;relat<(int64_t) MAXPOP;relat++) 
										{	temp[relat][0]=0;
											temp[relat][1]=0;
											temp[relat][2]=0;
											temp[relat][3]=0;
										};
										for(int snp=chrdivider[breaknubercm][chrtemp1][chrdividerrun].start;snp<chrdivider[breaknubercm][chrtemp1][chrdividerrun].end;snp++)		
										{	int snpvalue0=genomeoffpss[0][snp][chrtemp1];
											if (snpvalue0!=0 && snpvalue0!=3 )//&& phaseparent1tap[chrtemp1][snp]!=-1)
											{//	if (snp<snpphaseerror) snpvalue0=3-snpvalue0;
												hetsnp++;
												double maffloat=1.0*MAF[snp][chrtemp1]/MAXPOP/2;
												//if (phaseguesedparent1tap[chrtemp1][snp]==0) printf("STOP\n");
												int parent0indiv0=(phaseguesedparent1tap[chrtemp1][snp]==1)?(snpvalue0>>1):(snpvalue0&1);
												int parent1indiv0=(phaseguesedparent1tap[chrtemp1][snp]==1)?(snpvalue0&1):(snpvalue0>>1);
											//	printf("%d %d\n",phaseguesedparent1tap[chrtemp1][snp],snpvalue0);
												double pcontribu10=((parent0indiv0)-1.0*maffloat);
												double pcontribu11=((parent1indiv0)-1.0*maffloat);
												double maffloatdiviseur=(maffloat)*(2-maffloat*2)/3*2;
												double pconttab[4][2];
												pconttab[0][0]=pow(((0)-maffloat)/maffloatdiviseur*pcontribu10,1.0/powersnppihat);
												pconttab[0][1]=pow(((0)-maffloat)/maffloatdiviseur*pcontribu11,1.0/powersnppihat);
												pconttab[1][0]=pow(((1)-maffloat)/maffloatdiviseur*pcontribu10,1.0/powersnppihat);
												pconttab[1][1]=pow(((1)-maffloat)/maffloatdiviseur*pcontribu11,1.0/powersnppihat); 
												int nbsnpperchrby4=(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0));
												int snpmod4=(((snp%4)*2));
												int snpdiv4=snp/4;
												#pragma omp parallel for
												for(int64_t relat=0;relat<(int64_t) MAXPOP;relat++) 
												if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 && (pihatagainstall[relat]>(1.0/64+1.0/128)/2 || (relat%proprotiontokeep)==0 ) )
												{	int snpvalue1=(*((genomes[chrtemp1]+(unsigned long long) relat*nbsnpperchrby4 )+snpdiv4)>>snpmod4)&3;
													if (pconttab[snpvalue1&1][0]>0) temp[relat][0]=temp[relat][0]+(pconttab[snpvalue1&1][0]);
													else temp[relat][0]=0;
													if (temp[relat][0]<0) temp[relat][0]=0;
													else if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]<temp[relat][0]) pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]=temp[relat][0];
													
												//	pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+1][1]=pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+1][1]+pconttab[snpvalue1>>1][0];
													if (pconttab[snpvalue1>>1][0]>0) temp[relat][1]=temp[relat][1]+(pconttab[snpvalue1>>1][0]);
													else temp[relat][1]=0;
													if (temp[relat][1]<0) temp[relat][1]=0;
													else if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+1][0]<temp[relat][1]) pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+1][0]=temp[relat][1];
													
												//	pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][1]=pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][1]+pconttab[snpvalue1&1][1];
													if (pconttab[snpvalue1&1][1]>0) temp[relat][2]=temp[relat][2]+(pconttab[snpvalue1&1][1]);
													else temp[relat][2]=0;
													if (temp[relat][2]<0) temp[relat][2]=0;
													else if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]<temp[relat][2]) pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]=temp[relat][2];
													
													//pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+3][1]=pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+3][1]+pconttab[snpvalue1>>1][1];
													if (pconttab[snpvalue1>>1][1]>0) temp[relat][3]=temp[relat][3]+(pconttab[snpvalue1>>1][1]);
													else temp[relat][3]=0;
													if (temp[relat][3]<0) temp[relat][3]=0;
													else if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+3][0]<temp[relat][3]) pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+3][0]=temp[relat][3];
														
												};
											};
										};
										for(int relat=0;relat<MAXPOP;relat++) 
											if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 && (pihatagainstall[relat]>(1.0/64+1.0/128)/2 || (relat%proprotiontokeep)==0 ) )
										{	int value=0;
											//for(int value=0;value<2;value++)
											{	pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]=//(1+pihatagainstall[relat])*
																				(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]>pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+1][value]?
																												pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]:
																												pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+1][value]);
												pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]=//(1+pihatagainstall[relat])*
																				(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]>pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+3][value]?
																												pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]:
																												pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+3][value]);
											};
										};
										printf("%d %d %d %f %f \n",chrtemp1,0,chrdividerrun,
																	pihatagainstallchrMPphaseerror[chrtemp1][0][chrdividerrun*4][0],
																	pihatagainstallchrMPphaseerror[chrtemp1][0][chrdividerrun*4+2][0]);
									};
								};	
							};
							float thirdindice=3;//2+nbdivisionDEGREE4; DONE
							double corglobmax=0;						
							int chrmax1;						
							int chrmax2;						
							int dividmax1;						
							int dividmax2;						
							double sumall[23][25][2][2];
							double sumallsquare[23][25][2][2];
							#pragma omp parallel for 
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	for(int value=0;value<2;value++)
									{	sumall[chrtemp1][chrdividerrun][0][value]=0;
										sumall[chrtemp1][chrdividerrun][1][value]=0;
										sumallsquare[chrtemp1][chrdividerrun][0][value]=0;
										sumallsquare[chrtemp1][chrdividerrun][1][value]=0;
										for(int relat=0;relat<MAXPOP;relat++) 
											if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 && (pihatagainstall[relat]>(1.0/64+1.0/128)/2 || (relat%proprotiontokeep)==0 ) ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
										{	sumall[chrtemp1][chrdividerrun][0][value]=sumall[chrtemp1][chrdividerrun][0][value]+
												pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]),firstexponent)*
													(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]>0?1:-1);
											sumall[chrtemp1][chrdividerrun][1][value]=sumall[chrtemp1][chrdividerrun][1][value]+
												pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]),firstexponent)*
													(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]>0?1:-1);
											sumallsquare[chrtemp1][chrdividerrun][0][value]=sumallsquare[chrtemp1][chrdividerrun][0][value]+
												pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][value]),firstexponent*2);
											sumallsquare[chrtemp1][chrdividerrun][1][value]=sumallsquare[chrtemp1][chrdividerrun][1][value]+
												pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][value]),firstexponent*2);
										};
									};
								};
							};
							double allcor[23][MAXNBDIVISOR][23][MAXNBDIVISOR];
							double allseg[23][MAXNBDIVISOR][23][MAXNBDIVISOR];
							#pragma omp parallel for		
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	for(int  chrtemp2=chrtemp1;chrtemp2<23;chrtemp2++)		
									{	for(int chrdividerrun2=0;chrdividerrun2<nbchrdivider[breaknubercm][chrtemp2];chrdividerrun2++) 
										{	allseg[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2]=0;
										};
									};
								};
							};			
							int countseg[23][MAXNBDIVISOR];
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	countseg[chrtemp1][chrdividerrun]=1;
								};
							};
							
							float seuil1=0.5;
							#pragma omp parallel for		
							for(int relat=0;relat<MAXPOP;relat++) if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 && (pihatagainstall[relat]>(1.0/64+1.0/128)/2 || (relat%proprotiontokeep)==0 ) )
							{	
								#pragma omp parallel for		
								for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
								{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
									{	//double thresholdseg=seuil1/nbchrdivider[breaknubercm][chrtemp1]*nbsnpperchr[chrtemp1];				
										double thresholdseg=seuil1*(chrdivider[breaknubercm][chrtemp1][chrdividerrun].end-chrdivider[breaknubercm][chrtemp1][chrdividerrun].start);				
										if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]>thresholdseg || pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]>thresholdseg)
										{	if (pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]<thresholdseg || pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]<thresholdseg)
											{	countseg[chrtemp1][chrdividerrun]++;
												
											};
										};
									};
								};
							};
							double sumallseg=0;
							float penalty=0.7;
							double coefforcor=1.0/4;//nbdivisionDEGREE4
							#pragma omp parallel for		
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	for(int  chrtemp2=chrtemp1;chrtemp2<23;chrtemp2++)		
									{	for(int chrdividerrun2=0;chrdividerrun2<nbchrdivider[breaknubercm][chrtemp2];chrdividerrun2++) if (chrtemp1!=chrtemp2 || chrdividerrun<chrdividerrun2)
										{	double sumproduct; 
											sumproduct=0;
											double nbelem=0;
											for(int relat=0;relat<MAXPOP;relat++) 
												if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 && (pihatagainstall[relat]>(1.0/64+1.0/128)/2 || (relat%proprotiontokeep)==0 ) ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
											{	nbelem++;
												sumproduct=sumproduct+
															pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]),firstexponent)*
															pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4][0]),firstexponent);
											};
											double  cor1=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])/
														sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][0][0]-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp1][chrdividerrun][0][0])*
															((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][0][0]-sumall[chrtemp2][chrdividerrun2][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])));
											//	if (fabs(cor1)<fabs(cor11)*coefforcor) cor1=cor11;
											double sumproduct1=sumproduct;
											sumproduct=0;
											nbelem=0;
											for(int relat=0;relat<MAXPOP;relat++) 
												if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 && (pihatagainstall[relat]>(1.0/64+1.0/128)/2 || (relat%proprotiontokeep)==0 ) )//if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
											{	nbelem++;
												sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]),firstexponent)*
																	  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4+2][0]),firstexponent);
											};
											double  cor2=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp2][chrdividerrun2][1][0])/
														sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][0][0]-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp1][chrdividerrun][0][0])*
															((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][1][0]-sumall[chrtemp2][chrdividerrun2][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])));
									//		if (fabs(cor2)<fabs(cor22)*coefforcor) cor2=cor22;
											sumproduct=0;
											 nbelem=0;
											for(int relat=0;relat<MAXPOP;relat++) 
												if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 && (pihatagainstall[relat]>(1.0/64+1.0/128)/2 || (relat%proprotiontokeep)==0 ) ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
											{	nbelem++;
												sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]),firstexponent)*
																	  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4+2][0]),firstexponent);
											};
											double  cor3=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])/
														sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][1][0]-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp1][chrdividerrun][1][0])*
															((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][1][0]-sumall[chrtemp2][chrdividerrun2][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])));
										
											//		if (fabs(cor3)<fabs(cor33)*coefforcor) cor3=cor33;				
											sumproduct=0;
											 nbelem=0;
											for(int relat=0;relat<MAXPOP;relat++) 
												if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 && (pihatagainstall[relat]>(1.0/64+1.0/128)/2 || (relat%proprotiontokeep)==0 ) ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
											{	nbelem++;
												sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]),firstexponent)*	
																	  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4][0]),firstexponent);
											};
											double  cor4=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp2][chrdividerrun2][0][0])/
														sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][1][0]-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp1][chrdividerrun][1][0])*
															((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][0][0]-sumall[chrtemp2][chrdividerrun2][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])));
											//		if (fabs(cor4)<fabs(cor44)*coefforcor) cor4=cor44;				
											
											double corglob=(pow(fabs(cor1),thirdindice-1)*cor1 -
																							pow(fabs(cor2),thirdindice-1)*cor2+ 
																							pow(fabs(cor3),thirdindice-1)*cor3-
																							pow(fabs(cor4),thirdindice-1)*cor4);///
																						//	pow(countseg[chrtemp1][chrdividerrun]*countseg[chrtemp2][chrdividerrun2],1/(1));//
																							//pow(allseg[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2],1.0*(1+nbdivisionDEGREE4));[1] "35 92.513717948718" ->0.5
																														// [1] "35 92.4825" -> 0.33
																														
											if (corglob<0 && chrtemp1==chrtemp2) corglob=corglob*penalty;
											allcor[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2]=corglob;
											allcor[chrtemp2][chrdividerrun2][chrtemp1][chrdividerrun]=corglob;
											//printf("%d %d %d %d %f\n",chrtemp1,chrtemp2,chrdividerrun,chrdividerrun2,corglob
											
											if (fabs(corglob)>fabs(corglobmax))
											{	corglobmax=corglob;
												chrmax1=chrtemp1;						
												chrmax2=chrtemp2;						
												dividmax1=chrdividerrun;						
												dividmax2=chrdividerrun2;
												printf("%d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f  \n",chrtemp1,chrtemp2,chrdividerrun,chrdividerrun2,cor1,cor2,cor3,cor4,
																							sumproduct,nbelem,
																							sumallsquare[chrtemp1][chrdividerrun][0][0],sumallsquare[chrtemp2][chrdividerrun2][0][0],
																							sumall[chrtemp1][chrdividerrun][0][0],sumall[chrtemp2][chrdividerrun2][0][0],
																							sumproduct1,allseg[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2],corglobmax);
											};											
										};
									};
								};
							};
							
							printf("MAX %d %d %d %d %f\n",chrmax1,dividmax1,chrmax2,dividmax2,corglobmax);
							int group[23][MAXNBDIVISOR];
							
							int havemerged[23][MAXNBDIVISOR];
							int nbingroup[23][MAXNBDIVISOR];
							#pragma omp parallel for 	
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	havemerged[chrtemp1][chrdividerrun]=0;
									nbingroup[chrtemp1][chrdividerrun]=1;
									group[chrtemp1][chrdividerrun]=chrtemp1+chrdividerrun*23; 
								};
							};	
							int nbmerge=0;
							float mergingexponent[2];//nbdivisionDEGREE4
							mergingexponent[0]=2;
							mergingexponent[1]=5; 
							double sumcorglobal=0;
							int summingroup=0;
							int summaxgroup=0;
							int sumgroup=0;
							int nbpair=0;
							do 
							{	group[chrmax2][dividmax2]=chrmax1+dividmax1*23;
								printf("MERGE chr %d chunk %d to chr %d chnk %d :%d %d %d %d corglobmax %f \n",chrmax2,dividmax2,chrmax1,dividmax1,tabnbphaseright[chrmax1][dividmax1],
											tabnbphasewrong[chrmax1][dividmax1],tabnbphaseright[chrmax2][dividmax2],tabnbphasewrong[chrmax2][dividmax2],corglobmax);
								
								tappointdec[nbmerge].nbgroup1=nbingroup[chrmax1][dividmax1];
								tappointdec[nbmerge].nbgroup2=nbingroup[chrmax2][dividmax2];
								tappointdec[nbmerge].cor=fabs(corglobmax);
								
								nbmerge++;
								
								sumcorglobal=sumcorglobal+fabs(corglobmax);
								sumgroup=sumgroup+nbingroup[chrmax1][dividmax1]+nbingroup[chrmax2][dividmax2];
								summingroup=summingroup+((nbingroup[chrmax1][dividmax1]<nbingroup[chrmax2][dividmax2])?nbingroup[chrmax1][dividmax1]:nbingroup[chrmax2][dividmax2]);
								summaxgroup=summaxgroup+((nbingroup[chrmax1][dividmax1]>nbingroup[chrmax2][dividmax2])?nbingroup[chrmax1][dividmax1]:nbingroup[chrmax2][dividmax2]);
								nbingroup[chrmax1][dividmax1]=nbingroup[chrmax1][dividmax1]+nbingroup[chrmax2][dividmax2];
								if (corglobmax>0)
								{	havemerged[chrmax2][dividmax2]=1;
									for(int  chrtemp3=1;chrtemp3<23;chrtemp3++)	
									{	for(int chrdividerrun3=0;chrdividerrun3<nbchrdivider[breaknubercm][chrtemp3];chrdividerrun3++) 
										{	//allseg[chrtemp3][chrdividerrun3][chrmax1][dividmax1]=(allseg[chrtemp3][chrdividerrun3][chrmax1][dividmax1]+allseg[chrtemp3][chrdividerrun3][chrmax2][dividmax2])/2;
											//allseg[chrmax1][dividmax1][chrtemp3][chrdividerrun3]=allseg[chrtemp3][chrdividerrun3][chrmax1][dividmax1];
										};
									};
									if (nbdivisionDEGREE4==1) // "35 92.0095512820513"
										countseg[chrmax1][dividmax1]=countseg[chrmax1][dividmax1]+countseg[chrmax2][dividmax2];
									if (nbdivisionDEGREE4==0) countseg[chrmax1][dividmax1]=(countseg[chrmax1][dividmax1]>countseg[chrmax2][dividmax2])?//[1] "35 92.0783333333333"
																							countseg[chrmax1][dividmax1]:countseg[chrmax2][dividmax2];
									if (nbdivisionDEGREE4==2) countseg[chrmax1][dividmax1]=countseg[chrmax1][dividmax1]+countseg[chrmax2][dividmax2]+//[1] "35 86.1746794871795"
																							((countseg[chrmax1][dividmax1]>countseg[chrmax2][dividmax2])?
																							countseg[chrmax1][dividmax1]:countseg[chrmax2][dividmax2]);
									
									#pragma omp parallel for		
									for(int relat=0;relat<MAXPOP;relat++) 
										if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 && (pihatagainstall[relat]>(1.0/64+1.0/128)/2 || (relat%proprotiontokeep)==0 ) )										//if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
									{	//printf("%f %f %f %f\n",pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4],pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2],
										//		pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4],pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4+2]);
										for(int value=0;value<2;value++)
										{	double p1=pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]>0?1:-1);
											double p2=pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]>0?1:-1);
											double p3=pow(fabs(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4][value]>0?1:-1);
											double p4=pow(fabs(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4+2][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4+2][value]>0?1:-1);
											
											pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]=pow(fabs(p1+p3),1.0/mergingexponent[value])*(p1+p3>0?1:-1);
											pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]=pow(fabs(p2+p4),1.0/mergingexponent[value])*(p2+p4>0?1:-1);
										//	printf("%f %f\n",pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4],pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2]);
										};
									}	//866 7561957
									
								} else
								{	havemerged[chrmax2][dividmax2]=-1;
									#pragma omp parallel for
									for(int  chrtemp3=1;chrtemp3<23;chrtemp3++)		
									{	for(int chrdividerrun3=0;chrdividerrun3<nbchrdivider[breaknubercm][chrtemp3];chrdividerrun3++)
										{//	allseg[chrtemp3][chrdividerrun3][chrmax1][dividmax1]=(allseg[chrtemp3][chrdividerrun3][chrmax1][dividmax1]+allseg[chrtemp3][chrdividerrun3][chrmax2][dividmax2])/2;
											//allseg[chrmax1][dividmax1][chrtemp3][chrdividerrun3]=allseg[chrtemp3][chrdividerrun3][chrmax1][dividmax1];
										};
									};
									if (nbdivisionDEGREE4==2) // "35 92.0095512820513"
										countseg[chrmax1][dividmax1]=countseg[chrmax1][dividmax1]+countseg[chrmax2][dividmax2];
									if (nbdivisionDEGREE4==0) countseg[chrmax1][dividmax1]=(countseg[chrmax1][dividmax1]>countseg[chrmax2][dividmax2])?//[1] "35 92.0575"
																							countseg[chrmax1][dividmax1]:countseg[chrmax2][dividmax2];
									if (nbdivisionDEGREE4==1) countseg[chrmax1][dividmax1]=countseg[chrmax1][dividmax1]+countseg[chrmax2][dividmax2]+//[1] "35 86.1746794871795"
																							((countseg[chrmax1][dividmax1]>countseg[chrmax2][dividmax2])?
																							countseg[chrmax1][dividmax1]:countseg[chrmax2][dividmax2]);
								
									#pragma omp parallel for		
									for(int relat=0;relat<MAXPOP;relat++)
										if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 && (pihatagainstall[relat]>(1.0/64+1.0/128)/2 || (relat%proprotiontokeep)==0 ) ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
									{	for(int value=0;value<2;value++)
										{
											//	printf("%f %f %f %f\n",pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4],pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2],
											//		pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4],pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4+2]);
											double p1=pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]>0?1:-1);
											double p2=pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]>0?1:-1);
											double p3=pow(fabs(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4][value]>0?1:-1);
											double p4=pow(fabs(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4+2][value]),mergingexponent[value])*(pihatagainstallchrMPphaseerror[chrmax2][relat][dividmax2*4+2][value]>0?1:-1);
											pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]=pow(fabs(p1+p4),1.0/mergingexponent[value])*(p1+p4>0?1:-1);
											pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]=pow(fabs(p2+p3),1.0/mergingexponent[value])*(p2+p3>0?1:-1);
										//	printf("%f %f\n",pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4],pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2]);
										};
									};					
								};
								sumall[chrmax1][dividmax1][0][0]=0;
								sumall[chrmax1][dividmax1][1][0]=0;
								sumallsquare[chrmax1][dividmax1][0][0]=0;
								sumallsquare[chrmax1][dividmax1][1][0]=0;
								sumall[chrmax1][dividmax1][0][1]=0;
								sumall[chrmax1][dividmax1][1][1]=0;
								sumallsquare[chrmax1][dividmax1][0][1]=0;
								sumallsquare[chrmax1][dividmax1][1][1]=0;
								for(int relat=0;relat<MAXPOP;relat++) 
									if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 && (pihatagainstall[relat]>(1.0/64+1.0/128)/2 || (relat%proprotiontokeep)==0 ) ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
								{	for(int value=0;value<2;value++)
									{	sumall[chrmax1][dividmax1][0][value]=sumall[chrmax1][dividmax1][0][value]+
											pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]),firstexponent*(0.3+0.7))*
												(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]>0?1:-1);
										sumall[chrmax1][dividmax1][1][value]=sumall[chrmax1][dividmax1][1][value]+
											pow(fabs(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]),firstexponent*(0.3+0.7))*
											(pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]>0?1:-1);
										sumallsquare[chrmax1][dividmax1][0][value]=sumallsquare[chrmax1][dividmax1][0][value]+ pow((pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4][value]),firstexponent*2);
										sumallsquare[chrmax1][dividmax1][1][value]=sumallsquare[chrmax1][dividmax1][1][value]+ pow((pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4+2][value]),firstexponent*2);
									//	printf("%f %f\n",sumall[chrmax1][dividmax1][0],pihatagainstallchrMPphaseerror[chrmax1][relat][dividmax1*4]);
									};
								};
								int chrtemp1=chrmax1;
								int chrdividerrun=dividmax1;
								corglobmax=0;	
								
								double highestnew=0;		
								
								for(int  chrtemp2=1;chrtemp2<23;chrtemp2++)		
								{	for(int chrdividerrun2=0;chrdividerrun2<nbchrdivider[breaknubercm][chrtemp2];chrdividerrun2++) if (havemerged[chrtemp2][chrdividerrun2]==0 && (chrtemp1!=chrtemp2 || chrdividerrun!=chrdividerrun2))
									{	double sumproduct; 
										sumproduct=0;
										int nbelem=0;
										for(int relat=0;relat<MAXPOP;relat++) 
											if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 && (pihatagainstall[relat]>(1.0/64+1.0/128)/2 || (relat%proprotiontokeep)==0 ) )											//if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
										{	nbelem++;
											sumproduct=sumproduct+
														pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]),firstexponent)*
														pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4][0]),firstexponent);
										};
										double  cor1=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])/
													sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][0][0]-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp1][chrdividerrun][0][0])*
														((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][0][0]-sumall[chrtemp2][chrdividerrun2][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])));
										sumproduct=0;
												//	if (fabs(cor1)<fabs(cor11)*coefforcor) cor1=cor11;
										sumproduct=0;
										nbelem=0;
										for(int relat=0;relat<MAXPOP;relat++) 
											if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 && (pihatagainstall[relat]>(1.0/64+1.0/128)/2 || (relat%proprotiontokeep)==0 ) ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
										{	nbelem++;
											sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4][0]),firstexponent*(0.3+0.7))*
																  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4+2][0]),firstexponent*(0.3+0.7));
										};
										double  cor2=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp2][chrdividerrun2][1][0])/
													sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][0][0]-sumall[chrtemp1][chrdividerrun][0][0]*sumall[chrtemp1][chrdividerrun][0][0])*
														((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][1][0]-sumall[chrtemp2][chrdividerrun2][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])));
										sumproduct=0;
									//		if (fabs(cor2)<fabs(cor22)*coefforcor) cor2=cor22;
										
										sumproduct=0;
										 nbelem=0;
										for(int relat=0;relat<MAXPOP;relat++) 
											if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 && (pihatagainstall[relat]>(1.0/64+1.0/128)/2 || (relat%proprotiontokeep)==0 ) ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
										{	nbelem++;
											sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]),firstexponent*(0.3+0.7))*
																  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4+2][0]),firstexponent*(0.3+0.7));
										};
										double  cor3=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])/
													sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][1][0]-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp1][chrdividerrun][1][0])*
														((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][1][0]-sumall[chrtemp2][chrdividerrun2][1][0]*sumall[chrtemp2][chrdividerrun2][1][0])));
									
										sumproduct=0;
															
								//		if (fabs(cor3)<fabs(cor33)*coefforcor) cor3=cor33;				
										sumproduct=0;
										 nbelem=0;
										for(int relat=0;relat<MAXPOP;relat++) 
											if (pihatagainstall[relat]<seuilpihat[0] && pihatagainstall2[relat]<0.011 && (pihatagainstall[relat]>(1.0/64+1.0/128)/2 || (relat%proprotiontokeep)==0 ) ) //if (ID!=relat && relat!=IDp1loop && relat!=IDp2loop)
										{	nbelem++;
											sumproduct=sumproduct+pow(fabs(pihatagainstallchrMPphaseerror[chrtemp1][relat][chrdividerrun*4+2][0]),firstexponent*(0.3+0.7))*
																  pow(fabs(pihatagainstallchrMPphaseerror[chrtemp2][relat][chrdividerrun2*4][0]),firstexponent*(0.3+0.7));
										};
										double  cor4=(((nbelem)*sumproduct-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp2][chrdividerrun2][0][0])/
													sqrt(((nbelem)*sumallsquare[chrtemp1][chrdividerrun][1][0]-sumall[chrtemp1][chrdividerrun][1][0]*sumall[chrtemp1][chrdividerrun][1][0])*
														((nbelem)*sumallsquare[chrtemp2][chrdividerrun2][0][0]-sumall[chrtemp2][chrdividerrun2][0][0]*sumall[chrtemp2][chrdividerrun2][0][0])));
										sumproduct=0;
									//		if (fabs(cor4)<fabs(cor44)*coefforcor) cor4=cor44;				
										float exponantnumer=1;
										double corglob=(pow(fabs(cor1),thirdindice-1)*cor1 -
																						pow(fabs(cor2),thirdindice-1)*cor2+ 
																						pow(fabs(cor3),thirdindice-1)*cor3-
																						pow(fabs(cor4),thirdindice-1)*cor4)
																				*pow((nbingroup[chrtemp1][chrdividerrun]+nbingroup[chrtemp2][chrdividerrun2]),exponantnumer);///
																					//	pow(countseg[chrtemp1][chrdividerrun]*countseg[chrtemp2][chrdividerrun2],1/(1));//
																						//pow(allseg[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2],1.0*(1+nbdivisionDEGREE4));[1] "35 92.513717948718" ->0.5
																													// [1] "35 92.4825" -> 0.33
										
																													
										if (highestnew<fabs(corglob)) 
										{	highestnew=fabs(corglob); 
											printf("%f",corglob);
											printf("NB1 %d MB2 %d \n",nbingroup[chrtemp1][chrdividerrun],nbingroup[chrtemp2][chrdividerrun2]);  
										};
										allcor[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2]=corglob;										
										allcor[chrtemp2][chrdividerrun2][chrtemp1][chrdividerrun]=corglob;										
									};
								};
							//	printf("highestnew=%f\n",highestnew);
								nbpair=0;
								for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
								{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++) if (havemerged[chrtemp1][chrdividerrun]==0)
									{	for(int  chrtemp2=chrtemp1;chrtemp2<23;chrtemp2++)		
										{	for(int chrdividerrun2=0;chrdividerrun2<nbchrdivider[breaknubercm][chrtemp2];chrdividerrun2++) if (havemerged[chrtemp2][chrdividerrun2]==0 && (chrtemp1!=chrtemp2 || chrdividerrun<chrdividerrun2))
											{	nbpair++;
												double corglob=allcor[chrtemp1][chrdividerrun][chrtemp2][chrdividerrun2];
												if (fabs(corglob)>fabs(corglobmax))
												{//	printf("%d %d %d %d %f\n",chrtemp1,chrtemp2,chrdividerrun,chrdividerrun2,corglob	);
													corglobmax=corglob;
													chrmax1=chrtemp1;						
													chrmax2=chrtemp2;						
													dividmax1=chrdividerrun;						
													dividmax2=chrdividerrun2;
													printf("NB1 %d MB2 %d c %f \n",nbingroup[chrtemp1][chrdividerrun],nbingroup[chrtemp2][chrdividerrun2],corglob);  
													if (nbingroup[chrtemp1][chrdividerrun]==1 || nbingroup[chrtemp2][chrdividerrun2]==1)
													{	int chrone=(nbingroup[chrtemp1][chrdividerrun]==1)?chrtemp2:chrtemp1;
														int chrdividerrunone=(nbingroup[chrtemp1][chrdividerrun]!=1)?chrdividerrun2:chrdividerrun;
														int chrnotone=(nbingroup[chrtemp1][chrdividerrun]==1)?chrtemp2:chrtemp1;
														int chrdividerrunnotone=(nbingroup[chrtemp1][chrdividerrun]!=1)?chrdividerrun2:chrdividerrun;
														
														printf("1 is %d %d\n",chrnotone,chrdividerrunnotone);
														int numberattached=0;
														for(int chrdividerrunrun=0;chrdividerrunrun<nbchrdivider[breaknubercm][chrone];chrdividerrunrun++) 
														{//	printf("checking window %d\n",chrdividerrunrun);
															
															if (chrdividerrunrun==chrdividerrunone)
															{	//printf("window n %d is consider window\n", chrdividerrunrun);
															} else
															{	if (havemerged[chrone][chrdividerrunrun]==0)
																{//	printf("window n%d hasn t merge\n", chrdividerrunrun);
																} else 
																{	int hrdividerloop=chrdividerrunrun;	
																	int chrloop=chrone;
																	int phaseknown=havemerged[chrloop][hrdividerloop];
																	while (iter<1000 && havemerged[chrloop][hrdividerloop]!=0)
																	{	//printf("%d %d %d %d\n", phaseknown,chr,chrdivi,havemerged[chr][chrdivi]);	
																		
																		iter++;
																		phaseknown=phaseknown*(havemerged[chrloop][hrdividerloop]!=0?havemerged[chrloop][hrdividerloop]:1);
																		int savechr=chrloop;
																		chrloop=group[chrloop][hrdividerloop]%23;
																		hrdividerloop=group[savechr][hrdividerloop]/23;
																		
																	};
																	if 	(chrloop==chrnotone && hrdividerloop==chrdividerrunnotone)
																	{	printf("window n%d is attached to %d %d with phaseknown %d\n", chrdividerrunrun,chrnotone,chrdividerrunnotone,phaseknown);
																		numberattached=numberattached+phaseknown;
																	} else 
																	{	printf("window n%d is not attached to %d %d\n", chrdividerrunrun,chrnotone,chrdividerrunnotone);
																	};																		
																};
															};
														};
														if (corglobmax*numberattached<0) corglobmax=corglobmax*penalty;
													};
												};	
											};
										};
									};
								};
								printf("Nb pairs: %d\n",nbpair);
							} while (nbpair>0 && nbmerge<22*19-1);
							
							 nbphaseright=0;
							 nbphasewrong=0;
							
							
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		 
							{	for(int chrdividerrun=0;chrdividerrun<nbchrdivider[breaknubercm][chrtemp1];chrdividerrun++)
								{	int chr=chrtemp1;
									int chrdivi=chrdividerrun;
									int phaseknown=1;
									int iter=0;
									while (iter<1000 && havemerged[chr][chrdivi]!=0)
									{	//printf("%d %d %d %d\n", phaseknown,chr,chrdivi,havemerged[chr][chrdivi]);	
										iter++;
										phaseknown=phaseknown*(havemerged[chr][chrdivi]!=0?havemerged[chr][chrdivi]:1);
										int savechr=chr;
										chr=group[chr][chrdivi]%23;
										chrdivi=group[savechr][chrdivi]/23;
										
									} 
									if (iter>900) exit(0);
									printf("pahse chr %d div %d is %d\n",chrtemp1,chrdividerrun,phaseknown);
									chrdivider[breaknubercm][chrtemp1][chrdividerrun].phasing=phaseknown;
									for(int snp=chrdivider[breaknubercm][chrtemp1][chrdividerrun].start;snp<chrdivider[breaknubercm][chrtemp1][chrdividerrun].end;snp++)		
									{	if ((phaseknown+3)/2==2) 
										{	if (genomeoffpss[0][snp][chrtemp1]==1) 
											{	genomeoffpss[0][snp][chrtemp1]=2;
												*( genomes[chrtemp1]+(unsigned long long) (((((unsigned long long) ID+NBINDIV)*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)))=
													(*( genomes[chrtemp1]+(unsigned long long) (((((unsigned long long) ID+NBINDIV)*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4))) & 
														(~(3<<((snp%4)*2))))
															|	((2)<<((snp%4)*2));
											}
											else if (genomeoffpss[0][snp][chrtemp1]==2) 
											{	genomeoffpss[0][snp][chrtemp1]=1;
												*( genomes[chrtemp1]+(unsigned long long) (((((unsigned long long) ID+NBINDIV)*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4)))=
													(*( genomes[chrtemp1]+(unsigned long long) (((((unsigned long long) ID+NBINDIV)*(nbsnpperchr[chrtemp1]/4+((nbsnpperchr[chrtemp1]%4)>0)) )+snp/4))) & 
														(~(3<<((snp%4)*2))))
															|	((1)<<((snp%4)*2));
											};
										};
										if ((phaseknown+3)/2==phaseparent1tap[chrtemp1][snp]) 
										{	nbphaseright++;
											chrdivider[breaknubercm][0][0].nbright++;
											chrdivider[breaknubercm][chrtemp1][chrdividerrun].nbright++;
										} else 
										{	nbphasewrong++;
											chrdivider[breaknubercm][0][0].nbwrong++;
											chrdivider[breaknubercm][chrtemp1][chrdividerrun].nbwrong++;
										};
									
										//printf("%d %d %d %d\n",phaseknown,phaseparent1tap[chrtemp1][snp],nbphaseright,nbphasewrong);
									};
									printf("chr %d div %d nbphaseright %d nbphasewrong %d\n",chrtemp1,chrdividerrun,nbphaseright,nbphasewrong);
								};
								printf("chr %d nbphaseright %d nbphasewrong %d ratio %f \n",chrtemp1,nbphaseright,nbphasewrong,
																							(nbphaseright>nbphasewrong?nbphaseright:nbphasewrong)/(1.0+nbphaseright+nbphasewrong));
							};
							
							purcentage=(int64_t) 10000*
											(nbphaseright>nbphasewrong?nbphaseright:nbphasewrong)/(nbphaseright+nbphasewrong>0?nbphaseright+nbphasewrong:1);						
							for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)	
							{	dataperchr[chrtemp1][52+nbdivisionDEGREE4-18+3*powerpihatDEGREE4+breaknubercmloop-23]=purcentage;		
							};
						}
					};
				};					
			};	
		};
	} while (gen<MAXGEN);	
	return (0);
}

int GetCPUCount()
{ 	cpu_set_t cs;
	CPU_ZERO(&cs);
	sched_getaffinity(0, sizeof(cs), &cs);
	int count = 0;
	for (int i = 0; i < 20; i++)
	{  if (CPU_ISSET(i, &cs))
		count++;
	};
	return count;
}

int writeoutput(int chr, char pathoutput[] )
{	FILE *fp; 
    char filename[200]; 
    strcpy(filename,pathoutput);
	char number[100];
    sprintf(number, "%d", chr);
	strcat(filename,number);
	strcat(filename,".ped");
	printf("The file to write the data is %s\n", filename); 
	fp = fopen(filename, "w"); 
    if (fp == NULL) 
    {   printf("Could not open file %s\n", filename); 
        return (0); 
    };
	int32_t person;
	for(int snp=0;snp<nbsnpperchrinfile[chr];snp++)
	{	for(person=0;person<NbIndiv;person++)
		{	fprintf(fp,	"%d %d 0 0 -9 ",chr,person+1);
			int geno=(*((genomes[chr]+(unsigned long long) person*(nbsnpperchr[chr]/4+((nbsnpperchr[chr]%4)>0)) )+snp/4)>>(((snp%4)*2)))&3;
			fprintf(fp, "%d.0 %d.0 ", geno%2, geno/2);
		};
		fprintf(fp,"\n");
	};
	fclose (fp);
	return (0);
}


//
//  main program
//	


int main(int argc, char *argv[])
{	clock_t begin0 = clock();
	printdetail=1;
	int chrtodo1;
	char PathInput[200]; 
	char PathOutput[200]; 
	char PathMAF[200];
	for(int input=1;input<argc;input++)
	{/*	if( strncmp(argv[input], "-off", strlen("-off")) == 0 && input < argc-1) 
		{	off = atoi(argv[++input]);
		
		} 
		else if( strncmp(argv[input], "-parent1", strlen("-parent1")) == 0 && input < argc-1) 
		{	parent1 = atoi(argv[++input]);
		} 
		else if( strncmp(argv[input], "-parent2", strlen("-parent2")) == 0 && input < argc-1) 
		{	parent2 = atoi(argv[++input]);
		} */
		if( strncmp(argv[input], "-NbIndiv", strlen("-NbIndiv")) == 0 && input < argc-1) NbIndiv = atoi(argv[++input]); 
		else if( strncmp(argv[input], "-PathInput", strlen("-PathInput")) == 0 && input < argc-1) strcpy(PathInput,argv[++input]);
		else if( strncmp(argv[input], "-PathOutput", strlen("-PathOutput")) == 0 && input < argc-1) strcpy(PathOutput,argv[++input]);
		else if( strncmp(argv[input], "-PathMAF", strlen("-PathMAF")) == 0 && input < argc-1) strcpy(PathMAF,argv[++input]);
	};
	if (NbIndiv==0)
	{	printf("ERROR: Number of indivudals is zero or undefined\n");
		exit(0);
	}
	
	if (NbIndiv>NBINDIVMAX)
	{	printf("ERROR: Number of indivudals is higher than %d\n",NBINDIVMAX);
		exit(0);
	}
	
	
	nbsnpperchr[0]=330005*2;
	nbsnpperchr[1]=26229*2;
	
	nbsnpperchr[2]=26210*2;nbsnpperchr[3]=22209*2;nbsnpperchr[4]=20690*2;nbsnpperchr[5]=19027*2;nbsnpperchr[6]=18418*2;nbsnpperchr[7]=18367*2;nbsnpperchr[8]=16283*2;nbsnpperchr[9]=14990*2;
	nbsnpperchr[10]=16494*2;nbsnpperchr[11]=15818*2;nbsnpperchr[12]=16008*2;nbsnpperchr[13]=11510*2;nbsnpperchr[14]=10804*2;nbsnpperchr[15]=10884*2;nbsnpperchr[16]=12195*2;
	nbsnpperchr[17]=11486*2;nbsnpperchr[18]=10222*2;nbsnpperchr[19]=9806*2;nbsnpperchr[20]=8985*2;nbsnpperchr[21]=5227*2;nbsnpperchr[22]=5882*2;
	
	
	
	for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
	{	chrdivider[25][chrtemp1][0].start=0;
		chrdivider[25][chrtemp1][0].end=nbsnpperchr[chrtemp1];
		nbchrdivider[25][chrtemp1]=1;
	};
	srand(0);
	for(int method=0;method<4;method++)// 23
	{	for(int res=0;res<6;res++)// 23
		{	result[method][res]=0;
		};
	};
	FILE * fileresult;
	char number[100];
	char pathfile[200]; 
//	readID(pathID); 
//	readtrio();
//	readPO();
//	readgroundtruth(chrtodo1,chr2+1);
	if (printdetail) printf("read trio done \n");
	if (printdetail) printf("read ID done \n");
	//genomes = (unsigned char*) calloc ((unsigned long long) (1+nbsnpperchr[chrtodo1]/4)*MAXPOP, sizeof(char));
	printf("%d\n",sizeof(genomes));
	clock_t step1 = clock(); 
	float elapsed_secs1 = (float)(step1-begin0);
    printf("\nTotal time 1:%f cpu click so %f seconds\n",elapsed_secs1,elapsed_secs1/CLOCKS_PER_SEC );
	
	
	clock_t step11 = clock();
	float elapsed_secs11 = (float)(step11-step1);
    printf("\nTotal time 1:%f cpu click so %f seconds\n",elapsed_secs11,elapsed_secs11/CLOCKS_PER_SEC );
	for(int IDrun=0;IDrun<12;IDrun++)// if (ID<10) MAXNBTRIO
	{	distrigametic[IDrun]=0; 
	}
	for(int IDrun=0;IDrun<23;IDrun++)// if (ID<10) MAXNBTRIO
	{	distrigametickeep[IDrun]=0;
	}
	//
	//for(int ID=strattrio;ID<(trioend>MAXNBTRIO?MAXNBTRIO:trioend);ID++)// if (ID<10) MAXNBTRIO 
	
	//for(int parameter=0;parameter<1;parameter++)// if (ID<10) MAXNBTRIO 
		
	//caluclatenumbersegment(pathresult);
	for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)
	{	genomes[chrtemp1] = (unsigned char*) calloc ((unsigned long long) (1+nbsnpperchr[chrtemp1]/4)*NbIndiv, sizeof(char));
		readgenomelocal(PathInput,chrtemp1,100+chrtemp1,105,genomes[chrtemp1]);
		
	};
	return 0;

	for (int indiv=0;indiv<NbIndiv;indiv++)
	{	int strattrio=indiv;
		int ID=strattrio;
		printf("start ID %d\n",ID);
		loadsegment( indiv, 
						indiv,
						indiv,
						indiv,
						0,
						2,
						0,
						PathInput,
						PathMAF);
		compareressultindivallcombinV15D4( indiv,
						chrtodo1,
						0,
						indiv,
						0,
						0,
						indiv,
						indiv,
						0,
						2,
						0,
						PathOutput); 
		for(int  chrtemp1=1;chrtemp1<23;chrtemp1++)		
		{	writeoutput(chrtemp1,PathOutput);
		};
	}
	
	//int compareressultindivallcombinV2(int ID,int chr1, int chr2,int numtrio,int limitnbsnp,int run,int IDp1loop,int IDp2loop,int lenminseg)
	return 0;
}

