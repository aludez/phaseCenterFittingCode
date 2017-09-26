#include "AnitaConventions.h"
#include "PrettyAnitaEvent.h"
#include "UsefulAdu5Pat.h"
#include "AnitaGeomTool.h"
#include "TFile.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TString.h"
#include "TTree.h"
#include "TMath.h"
#include "TChain.h"
#include "RawAnitaHeader.h"
#include "CorrelationSummaryAnita4.h"
#include "TTreeIndex.h"
#include <iostream>
#include <fstream>

#define MAX_ANTENNAS 48

void fillArrays(double *eventNumberIndex, double *thetaWaveIndex, double *phiWaveIndex, int *antIndex1, int *antIndex2,  double *maxCorrTimeIndex, int *corrInds, int *polIndex, int setPol, bool* adjacent);

double fitObject(double *par, double *eventNumberIndex, double *thetaWaveIndex, double *phiWaveIndex, int *antIndex1, int *antIndex2, double *maxCorrTimeIndex, int *corrInds, int *polIndex, bool* adjacent);

double findSlope(double **x, double **y, int *count);
double getMean(double **x, int *count, int** corrs);
double getRMS(double **x, int *count);
double getDeltaTExpected(int ant1, int ant2, double thetaWave, double phiWave, double *deltaR, double *deltaZ, double *deltaPhi);

const int arrayMax = 5000000;
double eventNumberIndex[arrayMax];
double thetaWaveIndex[arrayMax];
double phiWaveIndex[arrayMax];
int antIndex1[arrayMax];
int antIndex2[arrayMax];
double maxCorrTimeIndex[arrayMax];
int corrInds[arrayMax];
int polIndex[arrayMax];
bool adjacent[arrayMax];

void fitFCN(Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t flag);

void fitLindaSplitCost(int which)
{
	int setPol = which; //pol H=1 V=0
	
	double stepSize = 0.0001;

	gSystem->Load("libMathMore.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("libGeom.so");
	gSystem->Load("libMinuit.so");

	AnitaGeomTool *agt = AnitaGeomTool::Instance();
	agt->usePhotogrammetryNumbers(1);

	fillArrays(eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, corrInds, polIndex, setPol, adjacent);
	
	TMinuit * myMin = new TMinuit(4*MAX_ANTENNAS);
	myMin->SetPrintLevel(-1);
	//myMin->SetObjectFit(fitObject);
	myMin->SetFCN(fitFCN);

	Double_t deltaCableDelays[MAX_ANTENNAS], deltaCableDelaysErr[MAX_ANTENNAS];
	Double_t deltaR[MAX_ANTENNAS], deltaRErr[MAX_ANTENNAS];
	Double_t deltaZ[MAX_ANTENNAS], deltaZErr[MAX_ANTENNAS];
	Double_t deltaPhi[MAX_ANTENNAS], deltaPhiErr[MAX_ANTENNAS];

	for(int j = 0; j < MAX_ANTENNAS; j++)
	{
		deltaR[j] = 0.;
		deltaZ[j] = 0.;
		deltaPhi[j] = 0.;
		deltaCableDelays[j] = 0.;
	
		char name[30];
		sprintf(name, "r%d", j);
		myMin->DefineParameter(j, name, deltaR[j], stepSize, 0., 0.);
		sprintf(name, "z%d", j);
		myMin->DefineParameter(j+MAX_ANTENNAS,name, deltaZ[j], stepSize, 0.,0.);
		sprintf(name, "phi%d", j);
		myMin->DefineParameter(j+MAX_ANTENNAS*2, name, deltaPhi[j], stepSize, 0., 0.);
		sprintf(name, "cable%d", j);
		myMin->DefineParameter(j+MAX_ANTENNAS*3, name, deltaCableDelays[j], stepSize, 0., 0.);
	}

	//myMin->FixParameter(MAX_ANTENNAS*3); //fixt0
	//fix things that we arent actively fitting for
	for (int y = 0; y < MAX_ANTENNAS; y++)
	{
		myMin->FixParameter(y);
		myMin->FixParameter(y+MAX_ANTENNAS);
		myMin->FixParameter(y+MAX_ANTENNAS*2);
	}

	myMin->Migrad();

	for (int y = 0; y < MAX_ANTENNAS; y++) myMin->GetParameter(y+MAX_ANTENNAS*3, deltaCableDelays[y], deltaCableDelaysErr[y]);

	//first step done
	TMinuit* myMin2 = new TMinuit(4*MAX_ANTENNAS);
	myMin2->SetPrintLevel(-1);
	//myMin->SetObjectFit(fitObject);
	myMin2->SetFCN(fitFCN);

	for(int j = 0; j < MAX_ANTENNAS; j++)
	{
		char name[30];
		sprintf(name, "r%d", j);
		deltaR[j] = -.15;
		myMin2->DefineParameter(j, name, deltaR[j], stepSize, 0., 0.);
		sprintf(name, "z%d", j);
		myMin2->DefineParameter(j+MAX_ANTENNAS,name, deltaZ[j], stepSize, 0.,0.);
		sprintf(name, "phi%d", j);
		myMin2->DefineParameter(j+MAX_ANTENNAS*2, name, deltaPhi[j], stepSize, 0., 0.);
		sprintf(name, "cable%d", j);
		myMin2->DefineParameter(j+MAX_ANTENNAS*3, name, deltaCableDelays[j], stepSize, 0., 0.);
	}

	for (int y = 0; y < MAX_ANTENNAS; y++)
	{
		//myMin2->FixParameter(y);
		myMin2->FixParameter(y+MAX_ANTENNAS);
		myMin2->FixParameter(y+MAX_ANTENNAS*2);
		myMin2->FixParameter(y+MAX_ANTENNAS*3);
	}

	myMin2->Migrad();
	
	for (int y = 0; y < MAX_ANTENNAS; y++) myMin2->GetParameter(y, deltaR[y], deltaRErr[y]);
	
	//second step done
	TMinuit* myMin3 = new TMinuit(4*MAX_ANTENNAS);
	myMin3->SetPrintLevel(-1);
	//myMin->SetObjectFit(fitObject);
	myMin3->SetFCN(fitFCN);

	for(int j = 0; j < MAX_ANTENNAS; j++)
	{
		char name[30];
		sprintf(name, "r%d", j);
		myMin3->DefineParameter(j, name, deltaR[j], stepSize, 0., 0.);
		sprintf(name, "z%d", j);
		myMin3->DefineParameter(j+MAX_ANTENNAS,name, deltaZ[j], stepSize, 0.,0.);
		sprintf(name, "phi%d", j);
		myMin3->DefineParameter(j+MAX_ANTENNAS*2, name, deltaPhi[j], stepSize, -1., 1.);
		sprintf(name, "cable%d", j);
		myMin3->DefineParameter(j+MAX_ANTENNAS*3, name, deltaCableDelays[j], stepSize, 0., 0.);
	}

	for (int y = 0; y < MAX_ANTENNAS; y++)
	{
		myMin3->FixParameter(y);
		myMin3->FixParameter(y+MAX_ANTENNAS);
		//myMin3->FixParameter(y+MAX_ANTENNAS*2);
		myMin3->FixParameter(y+MAX_ANTENNAS*3);
	}

	myMin3->Migrad();
	
	for (int y = 0; y < MAX_ANTENNAS; y++){
		myMin3->GetParameter(y+MAX_ANTENNAS*2, deltaPhi[y], deltaPhiErr[y]);
		//printf("%g\n", deltaPhi[y]);
	}

	//third step done
	TMinuit* myMin4 = new TMinuit(4*MAX_ANTENNAS);
	myMin4->SetPrintLevel(-1);
	//myMin->SetObjectFit(fitObject);
	myMin4->SetFCN(fitFCN);

	for(int j = 0; j < MAX_ANTENNAS; j++)
	{
		char name[30];
		sprintf(name, "r%d", j);
		myMin4->DefineParameter(j, name, deltaR[j], stepSize, 0., 0.);
		sprintf(name, "z%d", j);
		myMin4->DefineParameter(j+MAX_ANTENNAS,name, deltaZ[j], stepSize, 0.,0.);
		sprintf(name, "phi%d", j);
		myMin4->DefineParameter(j+MAX_ANTENNAS*2, name, deltaPhi[j], stepSize, 0., 0.);
		sprintf(name, "cable%d", j);
		myMin4->DefineParameter(j+MAX_ANTENNAS*3, name, deltaCableDelays[j], stepSize, 0., 0.);
	}

	for (int y = 0; y < MAX_ANTENNAS; y++)
	{
		myMin4->FixParameter(y);
		//myMin4->FixParameter(y+MAX_ANTENNAS);
		myMin4->FixParameter(y+MAX_ANTENNAS*2);
		myMin4->FixParameter(y+MAX_ANTENNAS*3);
	}

	myMin4->Migrad();
	
	for (int y = 0; y < MAX_ANTENNAS; y++){
		myMin4->GetParameter(y+MAX_ANTENNAS, deltaZ[y], deltaZErr[y]);
	}
	TMinuit* myMin5 = new TMinuit(4*MAX_ANTENNAS);
	myMin5->SetPrintLevel(-1);
	//myMin->SetObjectFit(fitObject);
	myMin5->SetFCN(fitFCN);

	for(int j = 0; j < MAX_ANTENNAS; j++)
	{
		char name[30];
		sprintf(name, "r%d", j);
		myMin5->DefineParameter(j, name, deltaR[j], stepSize, 0., 0.);
		sprintf(name, "z%d", j);
		myMin5->DefineParameter(j+MAX_ANTENNAS,name, deltaZ[j], stepSize, 0.,0.);
		sprintf(name, "phi%d", j);
		myMin5->DefineParameter(j+MAX_ANTENNAS*2, name, deltaPhi[j], stepSize, 0., 0.);
		sprintf(name, "cable%d", j);
		myMin5->DefineParameter(j+MAX_ANTENNAS*3, name, deltaCableDelays[j], stepSize, 0., 0.);
	}

	for (int y = 0; y < MAX_ANTENNAS; y++)
	{
		//myMin4->FixParameter(y);
		myMin5->FixParameter(y+MAX_ANTENNAS);
		myMin5->FixParameter(y+MAX_ANTENNAS*2);
		myMin5->FixParameter(y+MAX_ANTENNAS*3);
	}

	myMin5->Migrad();
	
	std::time_t now = std::time(NULL);
	std::tm * ptm = std::localtime(&now);
	char timeBuffer[32];
	std::strftime(timeBuffer, 32, "%Y_%m_%d_time_%H_%M_%S", ptm);

	TString polStr;
	if (setPol == 1) polStr = "H";
	if (setPol == 0) polStr = "V";
	
	ofstream newfile(Form("TEST/splitCost_%s.txt", polStr.Data())); 
	
	int currPol = 0;
	if (setPol == 0) currPol = 1;
	newfile << "Antenna POL deltaR(m) deltaPhi(deg) deltaZ(m)" << std::endl;
	for(int j = 0; j < MAX_ANTENNAS; j++)
	{
		myMin5->GetParameter(j, deltaR[j], deltaRErr[j]);
		std::cout << " deltaR[" << j << "] = " << deltaR[j] << " +/- " << deltaRErr[j] << std::endl; 

		myMin5->GetParameter(j+MAX_ANTENNAS*2, deltaPhi[j], deltaPhiErr[j]);
		std::cout << " deltaPhi[" << j << "] = " << deltaPhi[j] << " +/- " << deltaPhiErr[j] << std::endl; 
		
		myMin5->GetParameter(j+MAX_ANTENNAS, deltaZ[j], deltaZErr[j]);
		std::cout << " deltaZ[" << j << "] = " << deltaZ[j] << " +/- " << deltaZErr[j] << std::endl;

		newfile << j << "	" << currPol << " " << deltaR[j] << " " << deltaPhi[j]*TMath::RadToDeg() << " " << deltaZ[j] << " " << deltaCableDelays[j] << std::endl;

	}
}


void fitFCN(Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t flag)
{
	double diffErr = fitObject(par, eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, corrInds, polIndex, adjacent);
	f = diffErr;
}


void fillArrays(double* eventNumberIndex, double* thetaWaveIndex, double* phiWaveIndex, int* antIndex1, int* antIndex2,  double* maxCorrTimeIndex, int* corrInds, int* polIndex, int setPol, bool* adjacent)
{
	char patName[FILENAME_MAX];
	char corrName[FILENAME_MAX];
	char headName[FILENAME_MAX];
	AnitaGeomTool * agt = AnitaGeomTool::Instance();
	agt->usePhotogrammetryNumbers(1);

	Adu5Pat * pat = 0;
	CorrelationSummaryAnita4 * corr = 0;
	RawAnitaHeader * header = 0;
	int pol = 0;
	TChain * gpsChain = new TChain("adu5PatTree");
	TChain * corrChain = new TChain("corrTree");
	TChain * headChain = new TChain("headTree");
	
	
	//change runs
	
	for (int run = 123; run < 153; run++)
	{
		//if (run==135 || run==136 || run==137 || run==138) continue; //low snr cut
		if (run==129 || run==130 || run==131 || run==132) continue; //45degrees
		//if (run!=129 && run!=130 && run!=131 && run!=148 && run!=149 && run!=150 && run!=151) continue; //not 45degrees
		sprintf(patName,  "/project/kicp/avieregg/anitaIV/flight1617/root/run%d/gpsEvent%d.root",run,run);
		sprintf(headName, "/project/kicp/avieregg/anitaIV/flight1617/root/run%d/headFile%d.root",run,run);
		sprintf(corrName, "corrTrees/run%dCorrTree.root",run);
		headChain->Add(headName);	
		gpsChain->Add(patName);
		corrChain->Add(corrName);
	}

	gpsChain->SetBranchAddress("pat", &pat);
	corrChain->SetBranchAddress("corr", &corr);
	corrChain->SetBranchAddress("pol", &pol);
	headChain->SetBranchAddress("header", &header);
	headChain->BuildIndex("eventNumber");
	TTreeIndex *ind = (TTreeIndex*) headChain->GetTreeIndex();

	int maxEntry = corrChain->GetEntries();

	Double_t sourceLat, sourceLon, sourceAlt, timeOffset;
	
	Double_t phiWave, thetaWave, lower, upper;
	Double_t maxCorrTime, deltaTExpected;
	int ant1, ant2;

	int countIndex = 0;
	int countNotsave = 0;
	double antPhi[MAX_ANTENNAS] = {0};
	AnitaPol::AnitaPol_t HPOL = AnitaPol::kHorizontal;//pol
	AnitaPol::AnitaPol_t VPOL = AnitaPol::kVertical;
	if(setPol == 0) for (int ant = 0; ant < MAX_ANTENNAS; ant++) antPhi[ant] = agt->getAntPhiPositionRelToAftFore(ant, VPOL);
	if(setPol == 1) for (int ant = 0; ant < MAX_ANTENNAS; ant++) antPhi[ant] = agt->getAntPhiPositionRelToAftFore(ant, HPOL);
	Double_t additionalPhi = 22.5*TMath::DegToRad();
	Double_t twoPi = 2*TMath::Pi();

	double passed = 0;
	double tot = 0;
	printf("total events = %d\n", maxEntry);
	for(Long64_t entry=0; entry<maxEntry; entry++)
	{
		corrChain->GetEntry(entry);
		if(setPol == 1){
			if (pol == 0) continue;
		}
		if(setPol == 0){
			if (pol == 1) continue;
		}
		Long64_t gpsEntry = ind->GetEntryNumberWithIndex(corr->eventNumber, 0);
		if(gpsEntry < 0) continue;
		gpsChain->GetEntry(gpsEntry);

		UsefulAdu5Pat usefulPat(pat);
		usefulPat.getThetaAndPhiWave(AnitaLocations::LONGITUDE_WAIS_A4, AnitaLocations::LATITUDE_WAIS_A4, AnitaLocations::ALTITUDE_WAIS_A4, thetaWave, phiWave);

		for (int corrInd=0; corrInd < NUM_CORRELATIONS_ANITA4; corrInd++)
		{
			//pick correlations
			//this line ensures only close correlations are considered
			//if (corrInd > 11 && corrInd!=37 && corrInd!=38 && corrInd!=39 && corrInd!=49 && corrInd!=50 && corrInd!=51) continue;
			if (corrInd > 11) continue; //only nearest neighbors
			//if (corrInd!=1 && corrInd!=4 && corrInd!=6 && corrInd!=7 && corrInd!=8 && corrInd!=9 && corrInd!=10 && corrInd!=11) continue; //only using corrs with center ants
			//if (corrInd!=1 && corrInd!=4) continue; //only using corrs with only center ants

			maxCorrTime = corr->maxCorTimes[corrInd];
			ant1 = corr->firstAnt[corrInd];
			ant2 = corr->secondAnt[corrInd];
			if(setPol == 0){
				if(ant1 == 45 || ant2 == 45) continue; //BV14 was broken
			}

			deltaTExpected = usefulPat.getDeltaTExpected(ant1, ant2, AnitaLocations::LONGITUDE_WAIS_A4, AnitaLocations::LATITUDE_WAIS_A4, AnitaLocations::ALTITUDE_WAIS_A4);

			lower = antPhi[ant1] - additionalPhi;
			upper = antPhi[ant2] + additionalPhi;
			if (lower<0) lower+=twoPi;
			if (upper>twoPi) upper-=twoPi;

			if (lower>upper)
			{
				if (phiWave<TMath::Pi()) lower -=twoPi;
				else upper +=twoPi;
			}
			tot++;
			if((maxCorrTime-deltaTExpected)*(maxCorrTime-deltaTExpected) < 1) passed++;
			if (phiWave > lower && phiWave < upper && (maxCorrTime - deltaTExpected) * (maxCorrTime - deltaTExpected) < 1)
			{
				//printf("delta t = %g\n", maxCorrTime - deltaTExpected);
				thetaWaveIndex[countIndex] = thetaWave;
				phiWaveIndex[countIndex] = phiWave;
				maxCorrTimeIndex[countIndex] = maxCorrTime;
				eventNumberIndex[countIndex] = corr->eventNumber;
				antIndex1[countIndex] = ant1;
				antIndex2[countIndex] = ant2;
				polIndex[countIndex] = pol;
				corrInds[countIndex] = corrInd;
				if((corrInd>5 && corrInd<12) || corrInd>48) adjacent[countIndex] = true;
				else
				{
					adjacent[countIndex] = false;
				}
				countIndex++;
			}
		}
	}
	delete pat;
	delete corr;
	delete gpsChain;
	delete corrChain;

	printf("%g pass\n", passed/tot);

	antIndex1[countIndex] = -999;
}

double fitObject(double *par, double *eventNumberIndex, double *thetaWaveIndex, double *phiWaveIndex, int *antIndex1, int *antIndex2, double *maxCorrTimeIndex, int *corrInds, int *polIndex, bool* adjacent)
{
	int arrayNum = 500000;
	AnitaGeomTool *agt = AnitaGeomTool::Instance();
	agt->usePhotogrammetryNumbers(1);
	
	Double_t deltaR[MAX_ANTENNAS] = {0};
	Double_t deltaZ[MAX_ANTENNAS] = {0};
	Double_t deltaPhi[MAX_ANTENNAS] = {0};
	Double_t deltaCableDelays[MAX_ANTENNAS] = {0};

	for(int i = 0; i < MAX_ANTENNAS; i++)
	{
		deltaR[i] = par[i];
		deltaZ[i] = par[i + MAX_ANTENNAS];
		deltaPhi[i] = par[i + 2*MAX_ANTENNAS];
		deltaCableDelays[i] = par[i + 3*MAX_ANTENNAS];
	}

	//printf("deltaphi = %g, deltaR = %g, deltaZ = %g\n", deltaPhi[0], deltaR[0], deltaZ[0]);

	double **corrDT;
	double **corrPhiWave;
	int **corrs;
	double** adjacentAntDeltaT;
	double** verticalAntDeltaT;
	double** adjacentAntPhiWave;
	double** verticalAntPhiWave;

	corrDT = new double*[48];
	corrPhiWave = new double*[48];
	corrs = new int*[48];
	adjacentAntDeltaT = new double*[48];
	verticalAntDeltaT = new double*[48];
	adjacentAntPhiWave = new double*[48];
	verticalAntPhiWave = new double*[48];
	for (int i = 0; i < 48; i++)
	{
		corrDT[i] = new double[arrayNum];
		corrPhiWave[i] = new double[arrayNum];
		corrs[i] = new int[arrayNum];
		adjacentAntDeltaT[i] = new double[arrayNum];
		verticalAntDeltaT[i] = new double[arrayNum];
		adjacentAntPhiWave[i] = new double[arrayNum];
		verticalAntPhiWave[i] = new double[arrayNum];
	}

	Int_t countCorrs[48] = {0};
	Int_t countAdjacent[48] = {0};
	Int_t countVertical[48] = {0};

	Double_t additionalPhi = 22.5 * TMath::DegToRad();

	Double_t phiWave, thetaWave, deltaTExpected, maxCorrTime;

	Double_t lower, upper;
	Double_t twoPi = TMath::TwoPi();

	int ant1, ant2, vert;
	int arrayCount = 0;
	Int_t entry = 0;

	double sigmaNorm = 0;
	while(antIndex1[entry]!=-999)
	{
		thetaWave = thetaWaveIndex[entry];
		phiWave = phiWaveIndex[entry];

		maxCorrTime = maxCorrTimeIndex[entry];
		ant1 = antIndex1[entry];
		ant2 = antIndex2[entry];

		//printf("ant1 = %d, ant2 = %d, theta = %g, phi = %g, deltaR = %g, deltaZ = %g, deltaPhi = %g\n", ant1, ant2, thetaWave, phiWave, deltaR[ant1], deltaZ[ant1], deltaPhi[ant1]);
		deltaTExpected = getDeltaTExpected(ant1, ant2, thetaWave, phiWave, deltaR, deltaZ, deltaPhi) + deltaCableDelays[ant1] - deltaCableDelays[ant2];
		//printf("deltaT = %g\n", maxCorrTime - deltaTExpected);
		sigmaNorm += maxCorrTime;
		corrDT[ant1][countCorrs[ant1]] = maxCorrTime - deltaTExpected;
		corrPhiWave[ant1][countCorrs[ant1]] = phiWave;
		corrs[ant1][countCorrs[ant1]] = corrInds[entry];
		countCorrs[ant1]++;

		if(adjacent[entry])
		{
			adjacentAntDeltaT[ant1][countAdjacent[ant1]] = maxCorrTime - deltaTExpected;
			adjacentAntPhiWave[ant1][countAdjacent[ant1]] = phiWave;
			countAdjacent[ant1]++;
		} else {
			if(ant1 < 16)
			{
				if(ant2 < 32) vert = ant1;
				else vert = ant2;
			}
			else vert = ant1;
			verticalAntDeltaT[vert][countVertical[vert]] = maxCorrTime - deltaTExpected;
			verticalAntPhiWave[vert][countVertical[vert]] = phiWave;
			countVertical[vert]++;
		}


		entry++;
	}

	Double_t sumMeanAdj = getMean(corrDT, countCorrs, corrs);
	//Double_t sumRMSAdj = getRMS(corrDT,  countCorrs);
	Double_t sumGradAdj = findSlope(adjacentAntPhiWave, adjacentAntDeltaT, countAdjacent);
	Double_t sumGradVert = findSlope(verticalAntPhiWave, verticalAntDeltaT, countVertical);

	for (int ants = 0; ants < MAX_ANTENNAS; ants++)
	{
		delete [] corrPhiWave[ants];
		delete [] corrDT[ants];
		delete [] corrs[ants];
		delete [] adjacentAntPhiWave[ants];
		delete [] verticalAntPhiWave[ants];
		delete [] adjacentAntDeltaT[ants];
		delete [] verticalAntDeltaT[ants];
	}

	delete [] corrDT;
	delete [] corrPhiWave;
	delete [] corrs;
	delete [] adjacentAntPhiWave;
	delete [] verticalAntPhiWave;
	delete [] adjacentAntDeltaT;
	delete [] verticalAntDeltaT;

	//printf("meanA = %g\n", sumMeanAdj);
	//printf("rmsA = %g, rmsV = %g\n", sumRMSAdj, sumRMSVert);
	//printf("gradA = %g\n", sumGradAdj);

	double lambda = .01;
	//return (sumMeanAdj + sumRMSAdj + sumGradAdj + sumMeanVert + sumRMSVert + sumGradVert);
	//return (sumMeanAdj + sumGradAdj + sumMeanVert + sumGradVert);
	//return (sumMeanAdj + sumRMSAdj + sumMeanVert + sumRMSVert);
	return (sumMeanAdj);// + lambda*sumGradAdj + lambda*sumGradVert);
	//return sqrt(sumRMSAdj + sumRMSVert)/sigmaNorm;
	//return sqrt(sumRMSAdj + sumRMSVert)*sqrt(entry);
}

double findSlope(double **x, double **y, int *count)
{
	double sumX, sumY, sumXY, sumXX, sumRES, slope;
	sumX = 0; sumY = 0; sumXY = 0; sumXX = 0;
	double val = 0;

	Int_t wrappedPhi[16]={1,1,1,0,0,0,0,0,0,0,0,0,0,0,0};

	for (int ant = 0; ant < MAX_ANTENNAS; ant++)
	{
		sumX = 0; sumY = 0; sumXY = 0; sumXX = 0;
		for (int i = 0; i < count[ant]; i++)
		{
			double thisX = x[ant][i];
			
			if (wrappedPhi[ant%16])
			{
				if (thisX > TMath::Pi())
				{
					thisX -= TMath::TwoPi();
				}
			}
			sumX += thisX;
			sumY += y[ant][i];
			sumXY += thisX*y[ant][i];
			sumXX += thisX*thisX;
		}
		if (count[ant] == 0) continue;
		slope = (sumX*sumY - double(count[ant])*sumXY) / (sumX*sumX - double(count[ant])*sumXX);
		val += slope * slope;
	}
	return val;
}

double getMean(double **x, int *count, int **corrs)
{
	double mean = 0;
	double sumMean = 0;
	for (int ant = 0; ant < MAX_ANTENNAS; ant++)
	{
		for (int corr = 0; corr < 52; corr++)
		{
			mean = 0;
			for (int i = 0; i < count[ant]; i++)
			{
				//if(corrs[ant][i] == corr) mean += x[ant][i];
				if(corrs[ant][i] == corr) mean += x[ant][i]*x[ant][i];
			}
			if(count[ant] == 0) continue;
			//sumMean += mean*mean / (double(count[ant])*double(count[ant]));
			sumMean += mean / (double(count[ant]));
		}
	}
	return sumMean;
}

double getRMS(double **x, int *count)
{
	double rms = 0;
	double sumRMS = 0;
	for (int ant = 0; ant < MAX_ANTENNAS; ant++)
	{
		rms = 0;
		for (int i = 0; i < count[ant]; i++)
		{
			rms += x[ant][i]*x[ant][i];
		}
		if (count[ant] == 0) continue;
		sumRMS += rms / (count[ant]);
	}
	return sumRMS;
}

double getDeltaTExpected(int ant1, int ant2, double thetaWave, double phiWave, double *deltaR, double *deltaZ, double *deltaPhi)
{
	AnitaGeomTool *agt = AnitaGeomTool::Instance();
	agt->usePhotogrammetryNumbers(1);
	Double_t phi1 = agt->getAntPhiPositionRelToAftFore(ant1) + deltaPhi[ant1];
	Double_t r1 = agt->getAntR(ant1) + deltaR[ant1];
	Double_t z1 = agt->getAntZ(ant1) + deltaZ[ant1];

	Double_t phi2 = agt->getAntPhiPositionRelToAftFore(ant2) + deltaPhi[ant2];
	Double_t r2 = agt->getAntR(ant2) + deltaR[ant2];
	Double_t z2 = agt->getAntZ(ant2) + deltaZ[ant2];

	Double_t tanThetaW = TMath::Tan(thetaWave);
	Double_t part1 = z1 * tanThetaW - r1 * TMath::Cos(phiWave - phi1);
	Double_t part2 = z2 * tanThetaW - r2 * TMath::Cos(phiWave - phi2);

	return 1e9 * ((TMath::Cos(thetaWave) * (part1 - part2)) / C_LIGHT);
	//time in ns
}





