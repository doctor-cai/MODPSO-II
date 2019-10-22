/*==================================================================================================*/
// Please find details of the method from:
//
//  Cai, Qing, et al. "Network structural balance based on evolutionary multiobjective 
//  optimization: A two-step approach." IEEE Transactions on Evolutionary Computation,
//  19(6): 903-916, 2015.
/*==================================================================================================*/
// The source codes are free for research work. If you have any problem with the source codes, please 
// contact with:

//	Maoguo Gong, 
//		Key Laboratory of Intelligent Perception and Image Understanding of Ministry of Education,
//		University of Xidian,
//		Xi'an, Shaanxi Province 710071, China.
//		http://web.xidian.edu.cn/mggong/index.html
//		Email: gong@ieee.org or mggong@mail.xidian.edu.cn

//    Qing Cai
//		Key Laboratory of Intelligent Perception and Image Understanding of Ministry of Education,
//		University of Xidian,
//		Xi'an, Shaanxi Province 710071, China.
//		Email: 506183509@qq.com
/*==================================================================================================*/
// Programmer:		
//		Qing Cai
// Last Update:
//		Nov. 11, 2012
//		Apr. 10, 2016
//		May. 11, 2018
/*==================================================================================================*/
// The package is passed compiling under: MS Visual C++ in Windows XP. 
// The latest package is passed compiling under: MS Visual Studio 2015 in Windows 10. 
// After running the program, please find the final results in the corresponding document archives.
// All the network topology structures shown in the paper are drawn with the software Pajek.
/*==================================================================================================*/



#include "moead.h"
	

void main()
{

	strcpy(strFunctionType,"_TCH1"); //这个方式是最好的
//	strcpy(strFunctionType,"_TCH2"); //这个方式得到的结果太差了
//	strcpy(strFunctionType,"_PBI");  //这个方式的结果还可以，但缺乏多样性

/*-----------------实验发现就以下几个参数比较好，其他的都不行-----------------------*/
	max_gen         = 100;       // maximal number of generations
	
//	int  niche      = 40;        // neighborhood size //这个参数40 perfect
//	Tr              = 10;        // 0.1*populationsize //根据tradeoff between the convergence and diversity 
	niche           = 10;								  // Tr有待调整
	nr              = 2;
	int  rn         = 30;
/*----------------------------------------------------------------------------------*/

	popsize = 100;
	numObjectives = 2;
	int runtimes = rn;

	char *SignedData[] = {  "inter_relation5.txt", "slovene1.txt", \
							"gahuku2.txt", "data1.txt",            \
							"data2.txt",   "EGFR.mat",             \
							"Macrophage.txt", "yeast.txt", 	       \
							 "ecoli.txt",  "wiki_data.txt",             \
							 "slashdot_data.txt", "epinions_data.txt" };

	//remove是去掉度为0的,remove_remove是去掉度为1的，remove-remove-remove是去掉只有负边的节点的
  
	char *SignedDataLabel[] = { "sloeven1label.txt", "gahuku2label.txt", "datalabel.txt",\
							  "datalabel.txt",     "datalabel.txt", "", "", "", "", "", "", "", "" };

	int nnodes[] = {6, 10, 16, 28, 28, 329, 678, 690, 1461, 7114, 77357, 131828};

/*----------------------------------the signed LFR benchmark networks------------------------------*/
	char *LFRData[] = {""};
	
/*-----------------------the real labels of the LFR benchmark networks-------------------------*/
	char *LFRDataLabel[] ={ "LFRlabel\\community0.1.txt", "LFRlabel\\community0.2.txt",\
							"LFRlabel\\community0.3.txt",	"LFRlabel\\community0.4.txt",\
							"LFRlabel\\community0.5.txt",	"LFRlabel\\community0.6.txt"};


/*-------------------------------------------------------------------*/
/***********************testing for LFR networks**********************/
/*-------------------------------------------------------------------*/
/*  int LFRDataIndex = 35;
  optimization = 0;
  numVariables  = 500;
  mutate_posibility = 0.9;
  SignedFlag = 1;  // signed network
  AdjacentMatrix = new int *[numVariables];
  for ( int i = 0;i < numVariables;i++ )
  {
	AdjacentMatrix[i] = new int [numVariables];
  }
  node = new network[numVariables];
  
  strcpy(strGNExtend,LFRData[LFRDataIndex]);
  ReadFile(strGNExtend,AdjacentMatrix,numVariables,numVariables); 
  strcpy(strlabel,LFRDataLabel[5]);*/

/*-------------------------------------------------------------------*/
/***********************testing for signed networks*******************/
/*-------------------------------------------------------------------*/
	int RealIndex = 1;
	optimization = 0;
	numVariables = nnodes[RealIndex];
	SignedFlag = 1;  // signed network pm=0.9 perfect
	mutate_posibility = 0.9;
	AdjacentMatrix = (char**)malloc(sizeof(char*)*numVariables);
	for (int i = 0; i < numVariables; i++)
		AdjacentMatrix[i] = (char*)malloc(sizeof(char)*numVariables);
	node = new network[numVariables];
	strcpy(FilePath, "RealWorld\\");
	strcpy(NetworkName, SignedData[RealIndex]);
	strcpy(NetworkLabel, SignedDataLabel[RealIndex]);
	strcpy(FileName, FilePath);
	strcat(FileName, NetworkName);
	ReadFile(FileName, AdjacentMatrix, numVariables, numVariables);
	strcpy(LabelName, FilePath);
	strcat(LabelName, NetworkLabel);
//--------------------------------------------------------------//
	
	

//--------------------!!! warning !!!-------------------------------------------//
/*------------------------------------------------------------------------------*/
/*-------- do not change the following codes unless u really know them ---------*/
/*------------------------------------------------------------------------------*/

	NodeInformation();	

	seed = (seed + 111)%1235;	
	rnd_uni_init = -(long)seed;	

	TMOEAD  MOEAD;      
	
	
	for ( mt = 0; mt < runtimes; mt++ )
	{
		
		cout << "The " << mt+1 << "th" << " run"<< endl;
		MOEAD.run(max_gen, rn, RealIndex + 1, mt + 1);
		
	}


	for (int j = 0; j < runtimes; j++)
	{
		cout<<"Q = "<<Thirty_Run_modularity[j]<< setprecision(5) << "  ";
		if (check_label())	
		{
			cout<<"NMI = "<<Thirty_Run_NMI[j] << setprecision(5);
		}
		cout << endl;
	}


	if (check_label())	NMImax = *max_element(Thirty_Run_NMI.begin(),Thirty_Run_NMI.end());
	Qmax = *max_element(Thirty_Run_modularity.begin(),Thirty_Run_modularity.end());

	for (int k = 0; k < runtimes; k++)
	{
		if (check_label())	NMIavg += Thirty_Run_NMI[k];
		Qavg += Thirty_Run_modularity[k];
	}

	Qavg = Qavg/Thirty_Run_modularity.size();

	if (check_label())
	{
		NMIavg = NMIavg/Thirty_Run_NMI.size();
		cout << "max NMI = " << NMImax << "   " << "avg NMI = " << NMIavg << endl;
	}
	cout<<"max Q = "<<Qmax<<"    "<<"avg Q = "<<Qavg<<endl<<endl;
	cout<<"The end of the algorithm! Find the detailed results in the document archives"<<endl;
	
	/*for (int j = 0; j < runtimes; j++)
	{
		cout<<"Q = "<<Thirty_Run_modularity[j]<<"  ";
		cout<<"NMI = "<<Thirty_Run_NMI[j]<<endl;
	}
	
	NMImax = *max_element( Thirty_Run_NMI.begin(), Thirty_Run_NMI.end() );
	Qmax = *max_element( Thirty_Run_modularity.begin(), Thirty_Run_modularity.end() );

	for (int k = 0; k < runtimes; k++)
	{
		NMIavg += Thirty_Run_NMI[k];
		Qavg += Thirty_Run_modularity[k];
	}

	
	NMIavg = NMIavg/Thirty_Run_NMI.size();
	Qavg = Qavg/Thirty_Run_modularity.size();
	cout<<"max NMI = "<<NMImax<<"   "<<"avg NMI = "<<NMIavg<<endl;
	cout<<"max Q = "<<Qmax<<"    "<<"avg Q = "<<Qavg<<endl;*/
	
	/*
	for (mutate_posibility = 0.01; mutate_posibility<=0.9;mutate_posibility+=0.01)
	{
		cout<<"pm = "<<mutate_posibility<<"    ";
		MOEAD.run(sd, niche, max_gen, rn); 
	}//*/


	for (int i = 0; i < numVariables; i++ )
	{
		delete[] AdjacentMatrix[i];
	}
	delete[] AdjacentMatrix;

	Thirty_Run_modularity.clear();
	Thirty_Run_NMI.clear();
	pm_modularity.clear();
	chrom.clear();    
	NMI.clear();
	pm_NMI.clear();
	clusters.clear();
}
