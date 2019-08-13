
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

void ETable_sort(int size, double *ETable, int *ETable_seq)
{
    int i,j;
    int ESeq_Temp;
    double E_Temp;
    
    for(i=0;i<size-1;i++)
    {
        for(j=0;j<size-1-i;j++)
        {
            if(ETable[j]>ETable[j+1])
            {
                E_Temp=ETable[j];
                ETable[j]=ETable[j+1];
                ETable[j+1]=E_Temp;
                
                ESeq_Temp=ETable_seq[j];
                ETable_seq[j]=ETable_seq[j+1];
                ETable_seq[j+1]=ESeq_Temp;
            }
        }
    }
}

double PTD(int RouteL, int *Route1, double *DistVec, double *RouteProb, double D_inf, double **SpotsPos, double *PosCab, double *PosDest)
{
    int i,j;
    double CumuRouteLength,CumuProb;
    double PTD_E;
    
    DistVec[0]=pow(pow(PosCab[0]-SpotsPos[Route1[0]][0],2)+pow(PosCab[1]-SpotsPos[Route1[0]][1],2),0.5);
    for(i=1;i<RouteL;i++){
        DistVec[i]=pow(pow(SpotsPos[Route1[i-1]][0]-SpotsPos[Route1[i]][0],2)+pow(SpotsPos[Route1[i-1]][1]-SpotsPos[Route1[i]][1],2),0.5);
    }
    for(i=0;i<RouteL;i++){
        RouteProb[i]=SpotsPos[Route1[i]][2];
    }
    
    PTD_E=0;
    CumuRouteLength=0;
    CumuProb=1;
    CumuRouteLength+=DistVec[0];
    PTD_E+=CumuRouteLength*RouteProb[0];
    CumuProb*=1-RouteProb[0];
    
    for(i=1;i<RouteL;i++){
        CumuRouteLength+=DistVec[i];
        PTD_E+=CumuRouteLength*CumuProb*RouteProb[i];
        CumuProb*=1-RouteProb[i];
    }
    //D_inf=CumuRouteLength+pow(pow(PosDest[0]-SpotsPos[Route1[RouteL-1]][0],2)+pow(PosDest[1]-SpotsPos[Route1[RouteL-1]][1],2),0.5);
    PTD_E+=D_inf*CumuProb;
    
    return PTD_E;
}

int MyMin(int a, int b){
	if(a<b)
		return a;
	else
		return b;
}

int MyMax(int a, int b){
	if(a<b)
		return b;
	else
		return a;
}

void DisturbX_RL(int *Route_X0, int *Route_X1, int RouteL, int SpotsNum, int ***GridIndex, int **GridInNum, int **SpotsPosIndex, int MVLength, int XInterNum, int YInterNum)
{
    int i,j,ChangedPosition,NewCityNO,level,CountCity,XIndex,YIndex,NewCityNO_pre,NewCityLevel,CountCity_pre,flag_chosen;
    int level_array[100],BoundGrid[4][2];
    ChangedPosition=rand()%RouteL;
    
    //RL begins
    XIndex=SpotsPosIndex[Route_X0[ChangedPosition]][0];
    YIndex=SpotsPosIndex[Route_X0[ChangedPosition]][1];
    CountCity=GridInNum[XIndex][YIndex];

    level=0;

    BoundGrid[0][1]=XIndex;
    BoundGrid[1][1]=XIndex;
    BoundGrid[2][1]=YIndex;
    BoundGrid[3][1]=YIndex;

    level_array[level]=CountCity;
    
    BoundGrid[0][0]=BoundGrid[0][1];
    BoundGrid[1][0]=BoundGrid[1][1];
    BoundGrid[2][0]=BoundGrid[2][1];
    BoundGrid[3][0]=BoundGrid[3][1];

    while(CountCity<MVLength){
        level++;
        if(XIndex-level>=0){
            for(i=0;BoundGrid[2][0]+i<=BoundGrid[3][0];i++){
                CountCity+=GridInNum[XIndex-level][BoundGrid[2][0]+i];
            }
            BoundGrid[0][1]=XIndex-level;
        }
        if(XIndex+level<=XInterNum-1){
            for(i=0;BoundGrid[2][0]+i<=BoundGrid[3][0];i++){
                CountCity+=GridInNum[XIndex+level][BoundGrid[2][0]+i];
            }
            BoundGrid[1][1]=XIndex+level;
        }
        if(YIndex-level>=0){
            for(i=0;BoundGrid[0][0]+i<=BoundGrid[1][0];i++){
                CountCity+=GridInNum[BoundGrid[0][0]+i][YIndex-level];
            }
            BoundGrid[2][1]=YIndex-level;
        }
        if(YIndex+level<=YInterNum-1){
            for(i=0;BoundGrid[0][0]+i<=BoundGrid[1][0];i++){
                CountCity+=GridInNum[BoundGrid[0][0]+i][YIndex+level];
            }
            BoundGrid[3][1]=YIndex+level;
        }
        BoundGrid[0][0]=BoundGrid[0][1];
        BoundGrid[1][0]=BoundGrid[1][1];
        BoundGrid[2][0]=BoundGrid[2][1];
        BoundGrid[3][0]=BoundGrid[3][1];
        //? && or &
        if(XIndex-level>=0&&YIndex-level>=0){
            CountCity+=GridInNum[XIndex-level][YIndex-level];
        }
        if(XIndex-level>=0&&YIndex+level<=YInterNum-1){
            CountCity+=GridInNum[XIndex-level][YIndex+level];
        }
        if(XIndex+level<=XInterNum-1&&YIndex-level>=0){
            CountCity+=GridInNum[XIndex+level][YIndex-level];
        }
        if(XIndex+level<=XInterNum-1&&YIndex+level<=YInterNum-1){
            CountCity+=GridInNum[XIndex+level][YIndex+level];
        }
        level_array[level]=CountCity;
    }
    
    NewCityNO_pre=rand()%CountCity+1;
    for(i=0;i<=level;i++){
        if(NewCityNO_pre<=level_array[i])
            break;
    }
    NewCityLevel=i;
    
    
    level=NewCityLevel;
    
    if(level==0){
        NewCityNO=GridIndex[XIndex][YIndex][NewCityNO_pre-1];
    }else
    {
        flag_chosen=0;
        
    		BoundGrid[0][0]=MyMax(XIndex-(level-1),0);
    		BoundGrid[1][0]=MyMin(XIndex+level-1,XInterNum-1);
    		BoundGrid[2][0]=MyMax(YIndex-(level-1),0);
    		BoundGrid[3][0]=MyMin(YIndex+level-1,YInterNum-1);
    		
        CountCity=level_array[level-1];

        if(XIndex-level>=0){
            for(i=0;BoundGrid[2][0]+i<=BoundGrid[3][0];i++){
                CountCity_pre=CountCity;
                CountCity+=GridInNum[XIndex-level][BoundGrid[2][0]+i];
                if(NewCityNO_pre<=CountCity){
                    NewCityNO=GridIndex[XIndex-level][BoundGrid[2][0]+i][NewCityNO_pre-CountCity_pre-1];
                    flag_chosen=1;
                    break;
                }
            }
        }
        
        if(XIndex+level<=XInterNum-1&&flag_chosen==0){
            for(i=0;BoundGrid[2][0]+i<=BoundGrid[3][0];i++){
                CountCity_pre=CountCity;
                CountCity+=GridInNum[XIndex+level][BoundGrid[2][0]+i];
                if(NewCityNO_pre<=CountCity){
                    NewCityNO=GridIndex[XIndex+level][BoundGrid[2][0]+i][NewCityNO_pre-CountCity_pre-1];
                    flag_chosen=1;
                    break;
                }
            }
        }
        
        if(YIndex-level>=0&&flag_chosen==0){
            for(i=0;BoundGrid[0][0]+i<=BoundGrid[1][0];i++){
                CountCity_pre=CountCity;
                CountCity+=GridInNum[BoundGrid[0][0]+i][YIndex-level];
                if(NewCityNO_pre<=CountCity){
                    NewCityNO=GridIndex[BoundGrid[0][0]+i][YIndex-level][NewCityNO_pre-CountCity_pre-1];
                    flag_chosen=1;
                    break;
                }
            }
        }
        
        if(YIndex+level<=YInterNum-1&&flag_chosen==0){
            for(i=0;BoundGrid[0][0]+i<=BoundGrid[1][0];i++){
                CountCity_pre=CountCity;
                CountCity+=GridInNum[BoundGrid[0][0]+i][YIndex+level];
                if(NewCityNO_pre<=CountCity){
                    NewCityNO=GridIndex[BoundGrid[0][0]+i][YIndex+level][NewCityNO_pre-CountCity_pre-1];
                    flag_chosen=1;
                    break;
                }
            }
        }
        
        if(XIndex-level>=0&&YIndex-level>=0&&flag_chosen==0){
            CountCity_pre=CountCity;
            CountCity+=GridInNum[XIndex-level][YIndex-level];
            if(NewCityNO_pre<=CountCity){
                NewCityNO=GridIndex[XIndex-level][YIndex-level][NewCityNO_pre-CountCity_pre-1];
                flag_chosen=1;
            }
        }
        
        if(XIndex-level>=0&&YIndex+level<=YInterNum-1&&flag_chosen==0){
            CountCity_pre=CountCity;
            CountCity+=GridInNum[XIndex-level][YIndex+level];
            if(NewCityNO_pre<=CountCity){
                NewCityNO=GridIndex[XIndex-level][YIndex+level][NewCityNO_pre-CountCity_pre-1];
                flag_chosen=1;
            }
        }
        
        if(XIndex+level<=XInterNum-1&&YIndex-level>=0&&flag_chosen==0){
            CountCity_pre=CountCity;
            CountCity+=GridInNum[XIndex+level][YIndex-level];
            if(NewCityNO_pre<=CountCity){
                NewCityNO=GridIndex[XIndex+level][YIndex-level][NewCityNO_pre-CountCity_pre-1];
                flag_chosen=1;
            }
        }
        
        if(XIndex+level<=XInterNum-1&&YIndex+level<=YInterNum-1&&flag_chosen==0){
            CountCity_pre=CountCity;
            CountCity+=GridInNum[XIndex+level][YIndex+level];
            if(NewCityNO_pre<=CountCity){
                NewCityNO=GridIndex[XIndex+level][YIndex+level][NewCityNO_pre-CountCity_pre-1];
                flag_chosen=1;
            }
        }
        
        if(NewCityNO<0||NewCityNO>=SpotsNum){
            NewCityNO=0;
            printf("NewCityNO: %d.\n Error! NewCityNO is chosen as a wrong number.\n",NewCityNO);
        }
    }

    
    //RL ends
    

    for(i=0;i<RouteL;i++){
        Route_X1[i]=Route_X0[i];
    }
    Route_X1[ChangedPosition]=NewCityNO;
    for(i=0;i<RouteL;i++){
        if(NewCityNO==Route_X0[i])
            Route_X1[i]=Route_X0[ChangedPosition];
    }
}

void DisturbX(int *Route_X0, int *Route_X1, int RouteL, int SpotsNum, int **XYorder, int **XYindex, int MVLength)
{
    int i,j,ChangedIndex,ChangedPosition,ChangedToCityNum,ChangedToCityIndex,flag,LeftIndex,RightIndex,XOrY;
    ChangedPosition=rand()%RouteL;
    
    //NWV begins
    XOrY=rand()%2;
    LeftIndex=XYindex[XOrY][ChangedPosition]-MVLength/2;
    RightIndex=XYindex[XOrY][ChangedPosition]+MVLength/2;
    if(LeftIndex<0)
        LeftIndex=0;
    if(RightIndex>SpotsNum-1)
        RightIndex=SpotsNum-1;
    
    ChangedToCityIndex=LeftIndex+rand()%(RightIndex-LeftIndex+1);//It can reach [LeftIndex,RightIndex]
    ChangedToCityNum=XYorder[XOrY][ChangedToCityIndex];
    //NWV ends
    
    flag=1;
    for(i=0;i<RouteL;i++){
        Route_X1[i]=Route_X0[i];
    }
    Route_X1[ChangedPosition]=ChangedToCityNum;
    for(i=0;i<RouteL;i++){
        if(ChangedToCityNum==Route_X0[i])
            Route_X1[i]=Route_X0[ChangedPosition];
    }
}

void Initial_RL(int *Route_X0, int RouteL, int ***GridIndex, int **GridInNum, int XInterNum, int YInterNum, double *PosCab, double GridXInterL, double GridYInterL, double Xmin1, double Ymin1)
{
    int i,j,level,CountCity,XIndex,YIndex,CountCity_pre;
    int BoundGrid[4][2];
    
    //RL begins
    XIndex=floor((PosCab[0]-Xmin1)/GridXInterL);
    YIndex=floor((PosCab[1]-Ymin1)/GridYInterL);
    CountCity_pre=0;
    CountCity=GridInNum[XIndex][YIndex];


    level=0;

    BoundGrid[0][1]=XIndex;
    BoundGrid[1][1]=XIndex;
    BoundGrid[2][1]=YIndex;
    BoundGrid[3][1]=YIndex;

    for(i=0;i<MyMin(CountCity,RouteL);i++){
    	Route_X0[i]=GridIndex[XIndex][YIndex][i];
    }


    BoundGrid[0][0]=BoundGrid[0][1];
    BoundGrid[1][0]=BoundGrid[1][1];
    BoundGrid[2][0]=BoundGrid[2][1];
    BoundGrid[3][0]=BoundGrid[3][1];

    while(CountCity<RouteL){
        level++;
        
        if(XIndex-level>=0){
            for(i=0;BoundGrid[2][0]+i<=BoundGrid[3][0];i++){
            		CountCity_pre=CountCity;
                CountCity+=GridInNum[XIndex-level][BoundGrid[2][0]+i];
                for(j=CountCity_pre;j<MyMin(CountCity,RouteL);j++){
                	Route_X0[j]=GridIndex[XIndex-level][BoundGrid[2][0]+i][j-CountCity_pre];
                }
            }
            BoundGrid[0][1]=XIndex-level;
        }
        if(CountCity>=RouteL)
        		break;
        		
        if(XIndex+level<=XInterNum-1){
            for(i=0;BoundGrid[2][0]+i<=BoundGrid[3][0];i++){
            		CountCity_pre=CountCity;
                CountCity+=GridInNum[XIndex+level][BoundGrid[2][0]+i];
                for(j=CountCity_pre;j<MyMin(CountCity,RouteL);j++){
                	Route_X0[j]=GridIndex[XIndex+level][BoundGrid[2][0]+i][j-CountCity_pre];
                }
            }
            BoundGrid[1][1]=XIndex+level;
        }
        if(CountCity>=RouteL)
        		break;
        		
        if(YIndex-level>=0){
            for(i=0;BoundGrid[0][0]+i<=BoundGrid[1][0];i++){
            		CountCity_pre=CountCity;
                CountCity+=GridInNum[BoundGrid[0][0]+i][YIndex-level];
                for(j=CountCity_pre;j<MyMin(CountCity,RouteL);j++){
                	Route_X0[j]=GridIndex[BoundGrid[0][0]+i][YIndex-level][j-CountCity_pre];
                }
            }
            BoundGrid[2][1]=YIndex-level;
        }
        if(CountCity>=RouteL)
        		break;
        		
        if(YIndex+level<=YInterNum-1){
            for(i=0;BoundGrid[0][0]+i<=BoundGrid[1][0];i++){
            		CountCity_pre=CountCity;
                CountCity+=GridInNum[BoundGrid[0][0]+i][YIndex+level];
                for(j=CountCity_pre;j<MyMin(CountCity,RouteL);j++){
                	Route_X0[j]=GridIndex[BoundGrid[0][0]+i][YIndex+level][j-CountCity_pre];
                }
            }
            BoundGrid[3][1]=YIndex+level;
        }
        if(CountCity>=RouteL)
        		break;
        
        BoundGrid[0][0]=BoundGrid[0][1];
        BoundGrid[1][0]=BoundGrid[1][1];
        BoundGrid[2][0]=BoundGrid[2][1];
        BoundGrid[3][0]=BoundGrid[3][1];
        //? && or &
        if(XIndex-level>=0&&YIndex-level>=0){
        		CountCity_pre=CountCity;
            CountCity+=GridInNum[XIndex-level][YIndex-level];
            for(j=CountCity_pre;j<MyMin(CountCity,RouteL);j++){
                	Route_X0[j]=GridIndex[XIndex-level][YIndex-level][j-CountCity_pre];
            }
        }
        if(XIndex-level>=0&&YIndex+level<=YInterNum-1){
        		CountCity_pre=CountCity;
            CountCity+=GridInNum[XIndex-level][YIndex+level];
            for(j=CountCity_pre;j<MyMin(CountCity,RouteL);j++){
                	Route_X0[j]=GridIndex[XIndex-level][YIndex+level][j-CountCity_pre];
            }
        }
        if(XIndex+level<=XInterNum-1&&YIndex-level>=0){
        		CountCity_pre=CountCity;
            CountCity+=GridInNum[XIndex+level][YIndex-level];
            for(j=CountCity_pre;j<MyMin(CountCity,RouteL);j++){
                	Route_X0[j]=GridIndex[XIndex+level][YIndex-level][j-CountCity_pre];
            }
        }
        if(XIndex+level<=XInterNum-1&&YIndex+level<=YInterNum-1){
        		CountCity_pre=CountCity;
            CountCity+=GridInNum[XIndex+level][YIndex+level];
            for(j=CountCity_pre;j<MyMin(CountCity,RouteL);j++){
                	Route_X0[j]=GridIndex[XIndex+level][YIndex+level][j-CountCity_pre];
            }
        }
 
    }

}


int AcceptOrNot(double PTD_E0, double PTD_E1, double Temperature)
{
    int AcceptFlag;
    double AccProb;
    
    if(PTD_E1<PTD_E0)
        AcceptFlag=1;
    else{
        if(Temperature<0.000000000001)
            AcceptFlag=0;
        else if((PTD_E1-PTD_E0)/Temperature>30)
            AcceptFlag=0;
        else{
            AccProb=exp((PTD_E0-PTD_E1)/Temperature);
            if(rand()/RAND_MAX<AccProb)
                AcceptFlag=1;
            else
                AcceptFlag=0;
        }
    }
    return AcceptFlag;
}

int main(int argc, char** argv)
{
    int i,j,k,RepeatFlag;
    int SpotsNum,RouteL,RandomSeed,*Route_X0,*Route_X1;
    double PTD_E0,PTD_E1,Temperature,CoolAlpha,PosCab[2],PosDest[2],PTD_Record;
    int StepNum,AcceptFlag,StopCount,StopNum;
    int test1,test2,rank,size,AccNum;
    double temp_input,D_inf,LogarC;
    double **SpotsPos,*DistVec,*RouteProb;
    
    int RepeatedExpNum,OneExpNum,GroupAmount,Grank,Granks[1000],*ETable_seq,MixTimes,GroupNO,ParaRecord,FixedStopN;
    double x_time,x_time2,x_time3,*ETable,*PosRecord_temp,PTD_Stop_temp,ParaR_time,Single_time,FixedIniTemp;
    int rand_a,rand_b,rand_c;
    int ParaM[]={1,2,4,8,16,32,64,128};//ParaMax
    //double CoolBeta[]={100.0,2000.0,1800.0,1600.0,1400.0,1200.0,1000.0,800.0,600.0,400.0,200.0,100.0};//TestDiffBeta
    double CoolBeta[]={10000.0,1000.0,100.0,10.0,1.0,0.1,0.01,0.001};//TestDiffBeta
    double FinishedTimePercentage[]={2.0/40.0,4.0/40.0,6.0/40.0,8.0/40.0,10.0/40.0,12.0/40.0,14.0/40.0,16.0/40.0,18.0/40.0,20.0/40.0,22.0/40.0,24.0/40.0,26.0/40.0,28.0/40.0,30.0/40.0,32.0/40.0,34.0/40.0,36.0/40.0,38.0/40.0,40.0/40.0,};
    int ParaRecordMax,RepeatedExpIndex;
    double PTD_min;
    int MVLength,MVLength_Ini,**XYorder,**XYindex,TemperatureNo;//NMV
    double *Xpos,*Ypos;//NMV
    struct {
    	double PTD_val;
    	int rank_index;
    }Comp_E_in,Comp_E_out;
		//RL variables:
    double Xmax,Xmin,Ymax,Ymin,Xmax1,Xmin1,Ymax1,Ymin1,GridXInterL,GridYInterL;
    int XInterNum,YInterNum,**GridInNum,***GridIndex,**GridIndexCount,XIndexTemp,YIndexTemp,**SpotsPosIndex;
    int MV_K,MV_Acc,MV_Total,MoveMode,MinMVLength,MaxMVLength,SynData_RealData,main_XIndex,main_YIndex,main_CountCity,AccNum_pre,proc_PrintInter;

    FILE *SpotsPosFile,*accE,*PositionFile,*FinalFile;
    
    MPI_Comm Comm_World,Comm_OneExp;
    MPI_Group Group_World,Group_OneExp;
    
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    
    Comm_World=MPI_COMM_WORLD;
    MPI_Comm_group(Comm_World,&Group_World);
    
    
    
    
    
    
    
    
    //Parameter Settings:
    OneExpNum=1;
    RepeatedExpNum=10;
    ParaRecordMax=20;
    
    //RL:
    SynData_RealData=3;// 1 is syn, 2 is SF real data, and 3 is BJ real data.
    XInterNum=100;
    YInterNum=100;
    
    RouteL=10;
    FixedIniTemp=10.0;
    FixedStopN=1000000;
    
    
    
    
    
    

    
    
    

    //TemperatureNo=5;
    
    GroupNO=(int)(rank/OneExpNum);//If I want to submit the code 2 times, then there are 2 groups. The same as exp num.
    //ParaRecord=(int)(GroupNO/RepeatedExpNum);
    
    
    for(i=0;i<OneExpNum;i++){
        Granks[i]=GroupNO*OneExpNum+i;
    }
    MPI_Group_incl(Group_World,OneExpNum,Granks,&Group_OneExp);
    MPI_Comm_create(Comm_World,Group_OneExp,&Comm_OneExp);
    MPI_Comm_rank(Comm_OneExp,&Grank);
    GroupAmount=size/OneExpNum;
    ETable = (double*)malloc(OneExpNum*sizeof(double));
    ETable_seq = (int*)malloc(OneExpNum*sizeof(int));
    
    x_time=MPI_Wtime();
    x_time2=x_time-(int)x_time;
    x_time3=(10000*x_time2-(int)(10000*x_time2))*100;
    rand_a=1.0;//(getpid()+1);
    rand_b=(int)x_time3;
    rand_c=rank+1;
    RandomSeed=rand_a*rand_b*rand_c; //TH
    //RandomSeed=((int)getpid()+1)*(x_time+0.01)*100*(rank+1); //SW
    srand(RandomSeed);
    printf("%d %d %d %d\n",RandomSeed,rand_a,rand_b,rand_c);//TH
    
    //printf("%d %d %d %f\n",RandomSeed,getpid(),(int)getpid(),x_time);//SW
    
    if(SynData_RealData==1){
        if(rank==0){
            SpotsPosFile=fopen("SD_6467.txt","r");
            fscanf(SpotsPosFile,"%d",&SpotsNum);
        }
        D_inf=sqrt(2.0);//(double)RouteL*sqrt(2);
    }else if(SynData_RealData==2){
        if(rank==0){
            SpotsPosFile=fopen("RD_11892.txt","r");
            fscanf(SpotsPosFile,"%d",&SpotsNum);
        }
        D_inf=1;//(double)RouteL*sqrt(2);
    }else if(SynData_RealData==3){
        if(rank==0){
            SpotsPosFile=fopen("RD_BJ_21824.txt","r");
            fscanf(SpotsPosFile,"%d",&SpotsNum);
        }
        D_inf=1;//(double)RouteL*sqrt(2);
    }
    
    MPI_Bcast(&SpotsNum,1,MPI_INT,0,MPI_COMM_WORLD);
    
    PosRecord_temp=(double*)malloc(3*SpotsNum*sizeof(double));
    SpotsPos=(double**)malloc(SpotsNum*sizeof(double*));
    for(i=0;i<SpotsNum;i++){
        SpotsPos[i]=(double*)malloc(3*sizeof(double));
    }
    Xpos=(double*)malloc(SpotsNum*sizeof(double));
    Ypos=(double*)malloc(SpotsNum*sizeof(double));
    XYorder=(int**)malloc(2*sizeof(int*));
    XYorder[0]=(int*)malloc(SpotsNum*sizeof(int));
    XYorder[1]=(int*)malloc(SpotsNum*sizeof(int));
    XYindex=(int**)malloc(SpotsNum*sizeof(int*));
    XYindex[0]=(int*)malloc(SpotsNum*sizeof(int));
    XYindex[1]=(int*)malloc(SpotsNum*sizeof(int));
    if(rank==0){
        for(i=0;i<3*SpotsNum;i++){
            fscanf(SpotsPosFile,"%lf",&temp_input);
            PosRecord_temp[i]=temp_input;
        }
        fclose(SpotsPosFile);
    }
    MPI_Bcast(PosRecord_temp,3*SpotsNum,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for(i=0;i<SpotsNum;i++){
        SpotsPos[i][0]=PosRecord_temp[3*i+0];
        SpotsPos[i][1]=PosRecord_temp[3*i+1];
        SpotsPos[i][2]=PosRecord_temp[3*i+2];
        Xpos[i]=SpotsPos[i][0];
        Ypos[i]=SpotsPos[i][1];
        XYorder[0][i]=i;
        XYorder[1][i]=i;
    }
    
    //My Real Local move:
    
    GridInNum=(int**)malloc(XInterNum*sizeof(int*));
    for(i=0;i<XInterNum;i++){
        GridInNum[i]=(int*)malloc(YInterNum*sizeof(int));
    }
    for(i=0;i<XInterNum;i++){
        for(j=0;j<YInterNum;j++){
            GridInNum[i][j]=0;
        }
    }

    GridIndexCount=(int**)malloc(XInterNum*sizeof(int*));
    for(i=0;i<XInterNum;i++){
        GridIndexCount[i]=(int*)malloc(YInterNum*sizeof(int));
    }
    for(i=0;i<XInterNum;i++){
        for(j=0;j<YInterNum;j++){
            GridIndexCount[i][j]=0;
        }
    }
    
    SpotsPosIndex=(int**)malloc(SpotsNum*sizeof(int*));
    for(i=0;i<SpotsNum;i++){
        SpotsPosIndex[i]=(int*)malloc(2*sizeof(int));
    }
    
    Xmax=-0.1;
    Xmin=1.1;
    Ymax=-0.1;
    Ymin=1.1;
    for(i=0;i<SpotsNum;i++){
        if(SpotsPos[i][0]<Xmin)
            Xmin=SpotsPos[i][0];
        if(SpotsPos[i][0]>Xmax)
            Xmax=SpotsPos[i][0];
        if(SpotsPos[i][1]<Ymin)
            Ymin=SpotsPos[i][1];
        if(SpotsPos[i][1]>Ymax)
            Ymax=SpotsPos[i][1];
    }
    Xmax1=Xmax+0.1*(Xmax-Xmin)/(double)XInterNum;
    Xmin1=Xmin-0.1*(Xmax-Xmin)/(double)XInterNum;
    Ymax1=Ymax+0.1*(Ymax-Ymin)/(double)YInterNum;
    Ymin1=Ymin-0.1*(Ymax-Ymin)/(double)YInterNum;
    
    GridXInterL=(Xmax1-Xmin1)/(double)XInterNum;
    GridYInterL=(Ymax1-Ymin1)/(double)YInterNum;
    
    for(i=0;i<SpotsNum;i++){
        SpotsPosIndex[i][0]=floor((SpotsPos[i][0]-Xmin1)/GridXInterL);
        SpotsPosIndex[i][1]=floor((SpotsPos[i][1]-Ymin1)/GridYInterL);
        GridInNum[SpotsPosIndex[i][0]][SpotsPosIndex[i][1]]++;
    }
    
    GridIndex=(int***)malloc(XInterNum*sizeof(int**));
    for(i=0;i<XInterNum;i++){
        GridIndex[i]=(int**)malloc(YInterNum*sizeof(int*));
        for(j=0;j<YInterNum;j++){
            GridIndex[i][j]=(int*)malloc(GridInNum[i][j]*sizeof(int));
        }
    }

    for(i=0;i<SpotsNum;i++){
        XIndexTemp=SpotsPosIndex[i][0];
        YIndexTemp=SpotsPosIndex[i][1];
        GridIndex[XIndexTemp][YIndexTemp][GridIndexCount[XIndexTemp][YIndexTemp]++]=i;//Make sure it is as expected.
    }
    
    //End
    
    
    ETable_sort(SpotsNum,Xpos,XYorder[0]);
    ETable_sort(SpotsNum,Ypos,XYorder[1]);
    for(i=0;i<SpotsNum;i++){
        XYindex[0][XYorder[0][i]]=i;
        XYindex[1][XYorder[1][i]]=i;
    }
    
    Route_X0=(int*)malloc(RouteL*sizeof(int));
    Route_X1=(int*)malloc(RouteL*sizeof(int));
    RouteProb=(double*)malloc(RouteL*sizeof(double));
    DistVec=(double*)malloc(RouteL*sizeof(double));
    
    accE=fopen("accE.txt","a");
    FinalFile=fopen("FinalFile.txt","a");

    for(ParaRecord=0;ParaRecord<ParaRecordMax;ParaRecord++)
    {
        ParaR_time=MPI_Wtime();
        for(RepeatedExpIndex=0;RepeatedExpIndex<RepeatedExpNum;RepeatedExpIndex++)
        {
            Single_time=MPI_Wtime();
            
            
            
            
            
            //Parameter Settings:

            //MV:
            MVLength_Ini=0.1*SpotsNum;//SpotsNum;//40+300*ParaRecord;//floor(0.5*((double)RouteL)*log((double)SpotsNum));
            MV_K=10;//4010;//10+100*ParaRecord;
            MV_Total=100;
            MinMVLength=floor(0.01*(1+0)*SpotsNum);
            MaxMVLength=0.1*SpotsNum;
            CoolAlpha=pow((1.0-CoolBeta[6]/10000.0),1.0/FinishedTimePercentage[ParaRecord]);//search for FinishedTimePercentage and make adjustments
            MoveMode=0;// 0 is fixed move and 1 is dynamic
            
        
            
            
            
            if(SynData_RealData==1){
                PosCab[0]=0.5;
                PosCab[1]=0.5;
            }else if(SynData_RealData==2){
                PosCab[0]=0.4;
                PosCab[1]=0.7;
            }else if(SynData_RealData==3){
                PosCab[0]=0.4;
                PosCab[1]=0.7;
            }
            
            MVLength=MVLength_Ini;
            //MVLength=MyMin(MVLength_Ini,MaxMVLength);
            //MVLength=MyMax(MVLength,MinMVLength);
            PosDest[0]=750.0;
            PosDest[1]=750.0;
            //The core part begins:
            //StopNum=(int)(FixedStopN/OneExpNum);
            StopNum=FixedStopN;
            //StopNum=(int)(((int)pow(10,5+ParaRecord))/OneExpNum);
            PTD_Record=0;
            PTD_Stop_temp=0;
            StopCount=0;
            AccNum=0;
            proc_PrintInter=100000;
            AccNum_pre=-proc_PrintInter+1;
            test1=0;
            test2=0;
            
            MixTimes=0;
            
            //Temperature=pow(10,0+ParaRecord);//1000000;
            Temperature=FixedIniTemp;

            //CoolAlpha=1.0-CoolBeta[ParaRecord]/10000.0;
            //CoolAlpha=1.0-CoolBeta[TemperatureNo]/10000.0;
            /*
            for(i=0;i<RouteL;i++){
                RepeatFlag=0;
                while(RepeatFlag==0){
                    RepeatFlag=1;
                    Route_X0[i]=floor((double)SpotsNum*rand()/RAND_MAX);
                    if(Route_X0[i]==SpotsNum)
                        RepeatFlag=0;
                    for(j=0;j<i;j++)
                        if(Route_X0[j]==Route_X0[i])
                            RepeatFlag=0;
                }
            }*/
            Initial_RL(Route_X0,RouteL,GridIndex,GridInNum,XInterNum,YInterNum,PosCab,GridXInterL,GridYInterL,Xmin1,Ymin1);
            //MPI_Bcast(Route_X0,RouteL,MPI_INT,0,Comm_OneExp);
            PTD_E0=PTD(RouteL,Route_X0,DistVec,RouteProb,D_inf,SpotsPos,PosCab,PosDest);
            
            PTD_min=1000000.0;
            StepNum=0;
            MV_Acc=0;
            while(1){
				StepNum++;
                //DisturbX(Route_X0,Route_X1,RouteL,SpotsNum,XYorder,XYindex,MVLength);
                DisturbX_RL(Route_X0,Route_X1,RouteL,SpotsNum,GridIndex,GridInNum,SpotsPosIndex,MVLength,XInterNum,YInterNum);
                PTD_E1=PTD(RouteL,Route_X1,DistVec,RouteProb,D_inf,SpotsPos,PosCab,PosDest);
                AcceptFlag=AcceptOrNot(PTD_E0,PTD_E1,Temperature);
                if(AcceptFlag==1){
                    for(i=0;i<RouteL;i++)
                        Route_X0[i]=Route_X1[i];
                    PTD_E0=PTD_E1;
                    AccNum++;
                    MV_Acc++;
                }
                
                if(StepNum%proc_PrintInter==1&rank==0){
                    fprintf(accE,"%f %d %f %d %d %d %d\n",PTD_E0,StepNum,((double)(AccNum-AccNum_pre))/(double)proc_PrintInter,GroupNO,ParaRecord,RepeatedExpIndex,MVLength);
                    AccNum_pre=AccNum;
                    fflush(accE);
                }
                //if(StepNum==27)
                //if(StepNum%ParaM[ParaRecord]==0)//Seq has to go through this
                if(StepNum%ParaM[0]==0)
                {
                    /*for(j=0;j<OneExpNum;j++){
                        ETable_seq[j]=j;
                    }
                    MPI_Gather(&PTD_E0,1,MPI_DOUBLE,ETable,1,MPI_DOUBLE,0,Comm_OneExp);
                    if(Grank==0){
                        ETable_sort(OneExpNum,ETable,ETable_seq);
                    }
                    MPI_Bcast(ETable,OneExpNum,MPI_DOUBLE,0,Comm_OneExp);
                    MPI_Bcast(ETable_seq,OneExpNum,MPI_INT,0,Comm_OneExp);
                    MPI_Bcast(Route_X0,RouteL,MPI_INT,ETable_seq[0],Comm_OneExp);
                    PTD_E0=ETable[0];*/
										//printf("%f %d %d %d %d %d %d\n",PTD_E0,Grank,GroupNO,Route_X0[0],Route_X0[1],Route_X0[2],Route_X0[3]);
                    /*Comp_E_in.PTD_val=PTD_E0;
                    Comp_E_in.rank_index=Grank;
                    MPI_Allreduce(&Comp_E_in,&Comp_E_out,1,MPI_DOUBLE_INT,MPI_MINLOC,Comm_OneExp);                   
                    MPI_Bcast(Route_X0,RouteL,MPI_INT,Comp_E_out.rank_index,Comm_OneExp);
                    PTD_E0=Comp_E_out.PTD_val;*/
                    //printf("%f %d %d %d %d %d %d %d\n",PTD_E0,Grank,GroupNO,Comp_E_out.rank_index,Route_X0[0],Route_X0[1],Route_X0[2],Route_X0[3]);
		                
		                if(PTD_min>PTD_E0)
		                    PTD_min=PTD_E0;
		                if(PTD_E0==PTD_Record)
		                    StopCount++;
		                else{
		                    StopCount=0;
		                    PTD_Record=PTD_E0;
		                }
		                if(StopCount*ParaM[0]>StopNum)
		                    break;
                    
                }

                
                for(i=0;i<OneExpNum;i++)//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                		Temperature=Temperature*CoolAlpha;
                
                
                
                if(StepNum%MV_Total==0){
                    if(MoveMode==1){
                        MVLength=MyMax(MVLength+floor(MV_K*(((double)MV_Acc)/(double)MV_Total-0.44)),MinMVLength);
                        MVLength=MyMin(MVLength,MaxMVLength);
                    }
                	MV_Acc=0;
              	}
                
            }
            if(Grank==0){
                
                main_XIndex=floor((PosCab[0]-Xmin1)/GridXInterL);
                main_YIndex=floor((PosCab[1]-Ymin1)/GridYInterL);
                main_CountCity=GridInNum[main_XIndex][main_YIndex];
                
                fprintf(FinalFile,"%f %f %d %d %d %d %f %d %d %f %f %f %d %d %d ",PTD_min,PTD_E0,StepNum,ParaRecord,GroupNO,RepeatedExpIndex,Temperature,AccNum,StopNum,(double)AccNum/(double)StepNum,MPI_Wtime()-Single_time,CoolAlpha,MVLength,MinMVLength,main_CountCity);
                for(i=0;i<RouteL;i++)
                    fprintf(FinalFile,"%d ",Route_X0[i]);
                fprintf(FinalFile,"\n");
                fflush(FinalFile);
            }
        }
        //printf("ParaRecord Time: %f %d\n",(MPI_Wtime()-ParaR_time)/(double)RepeatedExpNum,ParaRecord);
    }

    for(i=0;i<size;i++)
        MPI_Bcast(&size,1,MPI_INT,i,MPI_COMM_WORLD);
    if(rank==0)
        printf("%f %f\n",MPI_Wtime()-x_time,x_time);//TH
    
    fclose(FinalFile);
    fclose(accE);

		//RL:
    for(i=0;i<XInterNum;i++){
        free(GridInNum[i]);
    }
    free(GridInNum);

    for(i=0;i<XInterNum;i++){
        free(GridIndexCount[i]);
    }
    free(GridIndexCount);
    
    for(i=0;i<SpotsNum;i++)
        free(SpotsPosIndex[i]);
    free(SpotsPosIndex);

    for(i=0;i<XInterNum;i++){
        for(j=0;j<YInterNum;j++)
            free(GridIndex[i][j]);
        free(GridIndex[i]);
    }
    free(GridIndex);


    //RL ends.

    for(i=0;i<SpotsNum;i++)
        free(SpotsPos[i]);
    free(SpotsPos);
    free(Xpos);
    free(Ypos);
    free(XYorder[0]);
    free(XYorder[1]);
    free(XYorder);
    free(XYindex[0]);
    free(XYindex[1]);
    free(XYindex);
    free(DistVec);
    free(Route_X0);
    free(Route_X1);
    free(RouteProb);
    MPI_Finalize();
    return 0;
}
