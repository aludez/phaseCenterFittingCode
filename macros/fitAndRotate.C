#include "AnitaConventions.h"
#include "PrettyAnitaEvent.h"
#include "UsefulAdu5Pat.h"
#include "AnitaGeomTool.h"
#include "TFile.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "RawAnitaHeader.h"
#include "CorrelationSummaryAnita3.h"
#include "TTreeIndex.h"
#include <iostream>
#include <fstream>

#define MAX_ANTENNAS 48

void fillArrays(double *eventNumberIndex, double *thetaWaveIndex, double *phiWaveIndex, int *antIndex1, int *antIndex2,  double *maxCorrTimeIndex, bool *adjacent, int *polIndex, int* centerIndex);

double fitObject(double *par, double *eventNumberIndex, double *thetaWaveIndex, double *phiWaveIndex, int *antIndex1, int *antIndex2, double *maxCorrTimeIndex, bool *adjacent, int *polIndex, int *centerIndex);

double findSlope(double **x, double **y, int *count);
double getMean(double **x, int *count);
double getRMS(double **x, int *count);
double getDeltaTExpected(int ant1, int ant2, double thetaWave, double phiWave, double *deltaR, double *deltaZ, double *deltaPhi);

const int arrayMax = 2000000;
double eventNumberIndex[arrayMax];
double thetaWaveIndex[arrayMax];
double phiWaveIndex[arrayMax];
int centerIndex[arrayMax];
int antIndex1[arrayMax];
int antIndex2[arrayMax];
double maxCorrTimeIndex[arrayMax];
bool adjacent[arrayMax];
int polIndex[arrayMax];

void fitFCN(Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t flag);

void fitAndRotate()
{
	double stepSize = 0.0001;

	gSystem->Load("libMathMore.so");
	gSystem->Load("libPhysics.so");
	gSystem->Load("libGeom.so");
	gSystem->Load("libMinuit.so");

	AnitaGeomTool *agt = AnitaGeomTool::Instance();
	agt->useKurtAnita3Numbers(1);

	fillArrays(eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent, polIndex, centerIndex);
	
	TMinuit * myMin[16];
	Double_t deltaR[MAX_ANTENNAS], deltaRErr[MAX_ANTENNAS];
	Double_t deltaZ[MAX_ANTENNAS], deltaZErr[MAX_ANTENNAS];
	Double_t deltaPhi[MAX_ANTENNAS], deltaPhiErr[MAX_ANTENNAS];

	for(int centAnt = 0; centAnt < 16; centAnt++)
	{
		myMin[centAnt] = new TMinuit(3*MAX_ANTENNAS + 1);
		myMin[centAnt]->SetPrintLevel(-1);
		//myMin[centAnt]->SetObjectFit(fitObject);
		myMin[centAnt]->SetFCN(fitFCN);
	
		for(int j = 0; j < MAX_ANTENNAS; j++)
		{
			deltaR[j] = 0.;
			deltaZ[j] = 0.;
			deltaPhi[j] = 0.;
	
			char name[30];
			sprintf(name, "r%d", j);
			myMin[centAnt]->DefineParameter(j, name, deltaR[j], stepSize, -.15, .15);
			sprintf(name, "z%d", j);
			myMin[centAnt]->DefineParameter(j+MAX_ANTENNAS,name, deltaZ[j], stepSize, -.15,.15);
			sprintf(name, "phi%d", j);
			myMin[centAnt]->DefineParameter(j+MAX_ANTENNAS*2, name, deltaPhi[j], stepSize, -.02, .02);
		}
		
		myMin[centAnt]->DefineParameter(3*MAX_ANTENNAS,"centerAnt",centAnt, stepSize, -1,1);
		myMin[centAnt]->FixParameter(3*MAX_ANTENNAS);
		myMin[centAnt]->Migrad();
	}
	
	std::time_t now = std::time(NULL);
	std::tm * ptm = std::localtime(&now);
	char timeBuffer[32];
	std::strftime(timeBuffer, 32, "%Y_%m_%d_time_%H_%M_%S", ptm);

	ofstream newfile(Form("TEST/phaseCenterNumbers_WAIS_%s.txt",timeBuffer));
	
	Double_t dR = 0;
	Double_t dRErr = 0;
	Double_t dZ = 0;
	Double_t dZErr = 0;
	Double_t dPhi = 0;
	Double_t dPhiErr = 0;
	
	for (int j = 0; j < MAX_ANTENNAS; j++)
	{
		deltaR[j] = 0.;
		deltaZ[j] = 0.;
		deltaPhi[j] = 0.;
	}
	
	for(int centAnt = 0; centAnt < 16; centAnt++)
	{
		std::cout << "Phi Sector " << centAnt << std::endl;
		for(int j = 0; j < MAX_ANTENNAS; j++)
		{
			
			myMin[centAnt]->GetParameter(j, dR, dRErr);
			deltaR[j] += dR/3.;
			std::cout << " deltaR[" << j << "] = " << dR << " +/- " << dRErr << std::endl; 

			myMin[centAnt]->GetParameter(j+MAX_ANTENNAS, dZ, dZErr);
			deltaZ[j] += dZ/3.;
			std::cout << " deltaZ[" << j << "] = " << dZ << " +/- " << dZErr << std::endl;

			myMin[centAnt]->GetParameter(j+MAX_ANTENNAS*2, dPhi, dPhiErr);
			deltaPhi[j] += dPhi/3.;
			std::cout << " deltaPhi[" << j << "] = " << dPhi << " +/- " << dPhiErr << std::endl; 

			newfile << j << "	" << dR << "	" << dZ << "	" << dPhi << std::endl;

		}
	}

	newfile << "Mean value of all fits" << std::endl;
	
	for (int j = 0; j < MAX_ANTENNAS; j++)
	{
		newfile << j << "	" << deltaR[j] << "	" << deltaZ[j] << "	" << deltaPhi[j] << std::endl;
		printf("R = %g, Z = %g, phi = %g\n", deltaR[j], deltaZ[j], deltaPhi[j]);
	}

}


void fitFCN(Int_t& npar, Double_t* gin, Double_t& f, Double_t* par, Int_t flag)
{
	double diffErr = fitObject(par, eventNumberIndex, thetaWaveIndex, phiWaveIndex, antIndex1, antIndex2, maxCorrTimeIndex, adjacent, polIndex, centerIndex);
	f = diffErr;
}


void fillArrays(double* eventNumberIndex, double* thetaWaveIndex, double* phiWaveIndex, int* antIndex1, int* antIndex2,  double* maxCorrTimeIndex, bool* adjacent, int* polIndex, int* centerIndex)
{
	char patName[FILENAME_MAX];
	char corrName[FILENAME_MAX];
	char headName[FILENAME_MAX];
	AnitaGeomTool * agt = AnitaGeomTool::Instance();
	agt->useKurtAnita3Numbers(1);

	Adu5Pat * pat = 0;
	CorrelationSummaryAnita3 * corr = 0;
	RawAnitaHeader * header = 0;
	int pol = 0;
	TChain * gpsChain = new TChain("adu5PatTree");
	TChain * corrChain = new TChain("corrTree");
	TChain * headChain = new TChain("headTree");
	
	for (int run = 120; run < 153; run++)
	{
		sprintf(patName,  "/project/kicp/avieregg/anitaIV/flight1617/root/run%d/gpsFile%d.root",run,run);
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
	AnitaPol::AnitaPol_t HPOL = AnitaPol::kHorizontal;
	for (int ant = 0; ant < MAX_ANTENNAS; ant++)
	{
		antPhi[ant] = agt->getAntPhiPositionRelToAftFore(ant, HPOL);
	}
	Double_t additionalPhi = 22.5*TMath::DegToRad();
	Double_t twoPi = TMath::TwoPi();
	
	for(Long64_t entry=0; entry<maxEntry; entry++)
	{
		corrChain->GetEntry(entry);
		if (pol == 0) continue; //V == 0 H == 1
		Long64_t gpsEntry = ind->GetEntryNumberWithIndex(corr->eventNumber, 0);
		if(gpsEntry < 0) continue;
		gpsChain->GetEntry(gpsEntry);

		UsefulAdu5Pat usefulPat(pat);
		usefulPat.getThetaAndPhiWaveWaisDivide(thetaWave, phiWave);

		for (int corrInd=0; corrInd < NUM_CORRELATIONS_ANITA3; corrInd++)
		{
			//this line ensures only close correlations are considered
			if (corrInd > 19 && corrInd!=37 && corrInd!=38 && corrInd!=39 && corrInd!=49 && corrInd!=50 && corrInd!=51) continue;

			maxCorrTime = corr->maxCorTimes[corrInd];
			ant1 = corr->firstAnt[corrInd];
			ant2 = corr->secondAnt[corrInd];

			deltaTExpected = usefulPat.getDeltaTExpected(ant1, ant2, AnitaLocations::LONGITUDE_WAIS, AnitaLocations::LATITUDE_WAIS, AnitaLocations::ALTITUDE_WAIS);

			lower = antPhi[ant1] - additionalPhi;
			upper = antPhi[ant2] + additionalPhi;
			if (lower<0) lower+=twoPi;
			if (upper>twoPi) upper-=twoPi;

			if (lower>upper)
			{
				if (phiWave<TMath::Pi()) lower -=twoPi;
				else upper +=twoPi;
			}

			if (phiWave > lower && phiWave < upper && (maxCorrTime - deltaTExpected) * (maxCorrTime - deltaTExpected) < 1)
			{
				//printf("delta t = %g\n", maxCorrTime - deltaTExpected);
				thetaWaveIndex[countIndex] = thetaWave;
				phiWaveIndex[countIndex] = phiWave;
				maxCorrTimeIndex[countIndex] = maxCorrTime;
				eventNumberIndex[countIndex] = corr->eventNumber;
				centerIndex[countIndex] = corr->centreAntenna;
				antIndex1[countIndex] = ant1;
				antIndex2[countIndex] = ant2;
				polIndex[countIndex] = pol;
				if ((corrInd > 5 && corrInd < 12) || corrInd>48) adjacent[countIndex]  = true;
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

	printf("%d events\n", countIndex);

	antIndex1[countIndex] = -999;
}

double fitObject(double *par, double *eventNumberIndex, double *thetaWaveIndex, double *phiWaveIndex, int *antIndex1, int *antIndex2, double *maxCorrTimeIndex, bool *adjacent, int *polIndex, int *centerIndex)
{
	int arrayNum = 100000;
	AnitaGeomTool *agt = AnitaGeomTool::Instance();
	agt->useKurtAnita3Numbers(1);
	
	Double_t deltaR[MAX_ANTENNAS] = {0};
	Double_t deltaZ[MAX_ANTENNAS] = {0};
	Double_t deltaPhi[MAX_ANTENNAS] = {0};

	for(int i = 0; i < MAX_ANTENNAS; i++)
	{
		deltaR[i] = par[i];
		deltaZ[i] = par[i + MAX_ANTENNAS];
		deltaPhi[i] = par[i + 2*MAX_ANTENNAS];
	}
	Double_t centerAnt = par[3*MAX_ANTENNAS];
	//printf("deltaphi = %g, deltaR = %g, deltaZ = %g\n", deltaPhi[0], deltaR[0], deltaZ[0]);

	double **adjacentAntDeltaT;
	double **adjacentAntPhiWave;
	double **verticalAntDeltaT;
	double **verticalAntPhiWave;

	adjacentAntDeltaT = new double*[48];
	adjacentAntPhiWave = new double*[48];
	verticalAntDeltaT = new double*[48];
	verticalAntPhiWave = new double*[48];
	for (int i = 0; i < 48; i++)
	{
		adjacentAntDeltaT[i] = new double[arrayNum];
		adjacentAntPhiWave[i] = new double[arrayNum];
		verticalAntDeltaT[i] = new double[arrayNum];
		verticalAntPhiWave[i] = new double[arrayNum];
	}

	Int_t countAdjacent[48] = {0};
	Int_t countVertical[48] = {0};
	AnitaPol::AnitaPol_t HPOL = AnitaPol::kHorizontal;

	Double_t additionalPhi = 22.5 * TMath::DegToRad();

	Double_t phiWave, thetaWave, deltaTExpected, maxCorrTime;

	Double_t lower, upper;
	Double_t twoPi = TMath::TwoPi();

	int ant1, ant2, vert;
	int arrayCount = 0;
	Int_t entry = 0;
	
	while(antIndex1[entry]!=-999)
	{
		if (centerAnt != centerIndex[entry])
		{
			entry++;
			continue;
		}
		thetaWave = thetaWaveIndex[entry];
		phiWave = phiWaveIndex[entry];

		maxCorrTime = maxCorrTimeIndex[entry];
		ant1 = antIndex1[entry];
		ant2 = antIndex2[entry];

		//printf("ant1 = %d, ant2 = %d, theta = %g, phi = %g, deltaR = %g, deltaZ = %g, deltaPhi = %g\n", ant1, ant2, thetaWave, phiWave, deltaR[ant1], deltaZ[ant1], deltaPhi[ant1]);
		deltaTExpected = getDeltaTExpected(ant1, ant2, thetaWave, phiWave, deltaR, deltaZ, deltaPhi);
		//printf("deltaT = %g\n", maxCorrTime - deltaTExpected);
		if(adjacent[entry])
		{
			adjacentAntDeltaT[ant1][countAdjacent[ant1]] = maxCorrTime - deltaTExpected;
			adjacentAntPhiWave[ant1][countAdjacent[ant1]] = phiWave;
			countAdjacent[ant1]++;
		}
		else
		{
			if (ant1 < 16)
			{
				if (ant2 <32) vert = ant1;
				else vert = ant2;
			}
			else vert = ant1;
			verticalAntDeltaT[vert][countVertical[vert]] = maxCorrTime - deltaTExpected;
			verticalAntPhiWave[vert][countVertical[vert]] = phiWave;
			countVertical[vert]++;
		}
		entry++;
	}

	Double_t sumMeanAdj = getMean(adjacentAntDeltaT, countAdjacent) * 10000;
	Double_t sumMeanVert = getMean(verticalAntDeltaT, countVertical) * 10000;
	Double_t sumRMSAdj = getRMS(adjacentAntDeltaT, countAdjacent) * 100;
	Double_t sumRMSVert = getRMS(verticalAntDeltaT, countAdjacent) * 100;
	Double_t sumGradAdj = findSlope(adjacentAntPhiWave, adjacentAntDeltaT, countAdjacent) * 10000;
	Double_t sumGradVert = findSlope(verticalAntPhiWave, verticalAntDeltaT, countVertical) * 10000;

	for (int ants = 0; ants < MAX_ANTENNAS; ants++)
	{
		delete [] adjacentAntPhiWave[ants];
		delete [] adjacentAntDeltaT[ants];
		delete [] verticalAntPhiWave[ants];
		delete [] verticalAntDeltaT[ants];
	}

	delete [] adjacentAntDeltaT;
	delete [] adjacentAntPhiWave;
	delete [] verticalAntDeltaT;
	delete [] verticalAntPhiWave;

	//printf("meanA = %g, meanV = %g\n", sumMeanAdj, sumMeanVert);
	//printf("rmsA = %g, rmsV = %g\n", sumRMSAdj, sumRMSVert);
	//printf("gradA = %g, gradV = %g\n", sumGradAdj, sumGradVert);

	//return (sumMeanAdj + sumRMSAdj + sumGradAdj + sumMeanVert + sumRMSVert + sumGradVert);
	//return (sumMeanAdj + sumRMSAdj + sumMeanVert + sumRMSVert);
	return (sumRMSAdj + sumRMSVert);
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
		slope = (sumX*sumY - count[ant]*sumXY) / (sumX*sumX - count[ant]*sumXX);
		val += slope * slope;
	}
	return val;
}

double getMean(double **x, int *count)
{
	double mean = 0;
	double sumMean = 0;
	for (int ant = 0; ant < MAX_ANTENNAS; ant++)
	{
		mean = 0;
		for (int i = 0; i < count[ant]; i++) mean += x[ant][i];
		if (count[ant] == 0) continue;
		sumMean += mean*mean / (count[ant]*count[ant]);
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
		for (int i = 0; i < count[ant]; i++) rms += x[ant][i]*x[ant][i];
		if (count[ant] == 0) continue;
		sumRMS += rms / count[ant];
	}
	return sumRMS;
}

double getDeltaTExpected(int ant1, int ant2, double thetaWave, double phiWave, double *deltaR, double *deltaZ, double *deltaPhi)
{
	AnitaGeomTool *agt = AnitaGeomTool::Instance();
	agt->useKurtAnita3Numbers(1);
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





