#include "FilteredAnitaEvent.h"
#include "AnitaEventSummary.h"
#include "AnitaConventions.h"
#include "Analyzer.h"
#include "AnitaDataset.h"
#include "BH13Filter.h"
#include "FilterStrategy.h"
#include "FFTtools.h"
#include "PrettyAnalysisWaveform.h"
#include "CorrelationSummaryAnita4.h"
#include "TChain.h"
#include "TTree.h"
#include "TString.h"
#include "Analyzer.h"
#include "AnitaGeomTool.h"
#include "AnitaEventCalibrator.h"
#include "TH2D.h"
#include "UCUtil.h"
#include "AnalysisConfig.h"
#include "AntennaPositions.h"

void doResolutions(const char * whichFile)
{
	TString inName = whichFile;
	FFTtools::loadWisdom("wisdom.dat");

	TString pointDir = "pointingGraphs/";
	TString psName = inName;
	TString blinding = "PolarityBlinded";
	TString fEnd = ".root";
	
	TString outDir = "summaryHistos/";
	TString outf = outDir + inName + fEnd;

	TChain* chain = new TChain("resolution");
	TChain* corrTree = new TChain("corrTree");

	for(int i = 123; i < 153; i++)
	//for(int i = 130; i < 133; i++)
	{
		if(i==129 || i==130 || i==131 || i==132) continue;
		//if(i==132) continue;
		TString name = pointDir + psName + TString::Itoa(i, 10) + fEnd;
		chain->Add(name);
		char corrName[FILENAME_MAX];
		sprintf(corrName, "corrTrees/run%dCorrTree.root", i);
		corrTree->Add(corrName);
	}

	AnitaEventSummary* sum = new AnitaEventSummary;
	chain->SetBranchAddress("summary", &sum);
	int isHoriz = 0;
	CorrelationSummaryAnita4* corr = 0;
	corrTree->SetBranchAddress("corr", &corr);
	corrTree->SetBranchAddress("pol", &isHoriz);

	Long_t nEntries = chain->GetEntries();
	chain->GetEntry(0);
	Double_t firstEvent = sum->eventNumber;
	chain->GetEntry(nEntries-1);
	Double_t lastEvent = sum->eventNumber;

	TH2D* HRES = new TH2D("HRES", "HRES", 100, -5, 5, 100, -2, 2);	
	TH2D* dThetavThetaH = new TH2D("Htheta", "Htheta", 200, 0, 40, 100, -2, 2);	
	TH2D* dPhivPhiH = new TH2D("Hphi", "Hphi", 360, 0, 360, 100, -5, 5);	
	TH2D* dThetavPhiH = new TH2D("HthetaPhi", "HthetaPhi", 200, 0, 360, 100, -2, 2);	
	TH2D* dPhivThetaH = new TH2D("HphiTheta", "HphiTheta", 200, 0, 40, 100, -5, 5);	
	TH2D* dThetavSNRH = new TH2D("HthetaSNR", "HthetaSNR", 200, 0, 60, 100, -2, 2);	
	TH2D* dPhivSNRH = new TH2D("HphiSNR", "HphiSNR", 200, 0, 60, 100, -5, 5);	
	TH2D* dThetavTimeH=new TH2D("HthetaTime", "HthetaTime", 200,firstEvent, lastEvent, 100, -2, 2);	
	TH2D* dPhivTimeH = new TH2D("HphiTime", "HphiTime", 200, firstEvent, lastEvent, 100, -5, 5);	
	TH2D* dThetavDistanceH=new TH2D("HthetaDist", "HthetaDist", 200,0, 700, 100, -2, 2);	
	TH2D* dPhivDistanceH = new TH2D("HphiDist", "HphiDist", 200, 0, 700, 100, -5, 5);	
	
	TH2D* VRES = new TH2D("VRES", "VRES", 100, -5, 5, 100, -2, 2);	
	TH2D* dThetavThetaV = new TH2D("Vtheta", "Vtheta", 200, 0, 40, 100, -2, 2);	
	TH2D* dPhivPhiV = new TH2D("Vphi", "Vphi", 360, 0, 360, 100, -5, 5);	
	TH2D* dThetavPhiV = new TH2D("VthetaPhi", "VthetaPhi", 200, 0, 360, 100, -2, 2);	
	TH2D* dPhivThetaV = new TH2D("VphiTheta", "VphiTheta", 200, 0, 40, 100, -5, 5);	
	TH2D* dThetavSNRV = new TH2D("VthetaSNR", "VthetaSNR", 200, 0, 60, 100, -2, 2);	
	TH2D* dPhivSNRV = new TH2D("VphiSNR", "VphiSNR", 200, 0, 60, 100, -5, 5);	
	TH2D* dThetavTimeV=new TH2D("VthetaTime", "VthetaTime", 200,firstEvent, lastEvent, 100, -2, 2);	
	TH2D* dPhivTimeV = new TH2D("VphiTime", "VphiTime", 200, firstEvent, lastEvent, 100, -5, 5);	
	TH2D* dThetavDistanceV=new TH2D("VthetaDist", "VthetaDist", 200,0, 700, 100, -2, 2);	
	TH2D* dPhivDistanceV = new TH2D("VphiDist", "VphiDist", 200, 0, 700, 100, -5, 5);	
	
	for (int j = 0; j < corrTree->GetEntries(); j ++)
	{
		int whichPeak = 0;
		corrTree->GetEntry(j);
		chain->GetEntry(j);
		AnitaEventSummary::PointingHypothesis point;

		if(isHoriz == 0) point = sum->peak[1][whichPeak]; //V = 1 H = 0 (1st index)
		if(isHoriz == 1) point = sum->peak[0][whichPeak]; //V = 1 H = 0 (1st index)
		whichPeak++;

		AnitaEventSummary::SourceHypothesis w = sum->wais;
		double waistheta = w.theta;
		double waisphi = w.phi;
		double theta1 = point.theta;
		double phi1 = point.phi;
		double dist = w.distance/1e3;
		double dtheta = theta1 - waistheta;
		double dphi = phi1 - waisphi;
		if(dphi < -180) dphi += 360;
		if(dphi >= 180) dphi -= 360;

		while(fabs(dphi) > 5 || fabs(dtheta) > 5) 
		{
			if(whichPeak == 3) break;
			if(isHoriz == 0) point = sum->peak[1][whichPeak];
			if(isHoriz == 1) point = sum->peak[0][whichPeak];
			whichPeak++;
			theta1 = point.theta;
			phi1 = point.phi;
			dphi = phi1 - waisphi;
			dtheta = theta1 - waistheta;
			if(dphi < -180) dphi += 360;
			if(dphi >= 180) dphi -= 360;
		}
		//if((fabs(dphi) > 2.5 && fabs(dphi) < 5) || (fabs(dtheta) > 1. && fabs(dtheta) < 2.5)) printf("ev = %u, dphi = %g, dtheta = %g\n", sum->eventNumber, dphi, waisphi);
		//if((fabs(dtheta) > 1. && fabs(dtheta) < 2.5)) printf("ev = %u, phi = %g, dtheta = %g\n", sum->eventNumber, waisphi, dtheta);
		//if(fabs(dphi) > 5 || fabs(dtheta) > 2.5) printf("ev = %u, dphi = %g, dtheta = %g\n", sum->eventNumber, dphi, dtheta);
		if(isHoriz == 1)
		{
			HRES->Fill(dphi, dtheta);
			dThetavThetaH->Fill(waistheta, dtheta);
			dPhivPhiH->Fill(waisphi, dphi);
			dPhivThetaH->Fill(waistheta, dphi); 
			dThetavPhiH->Fill(waisphi, dtheta);
			dThetavSNRH->Fill(corr->snr/2., dtheta);
			dPhivSNRH->Fill(corr->snr/2., dphi);
			dPhivTimeH->Fill(sum->eventNumber, dphi);
			dThetavTimeH->Fill(sum->eventNumber, dtheta);
			dPhivDistanceH->Fill(dist, dphi);
			dThetavDistanceH->Fill(dist, dtheta);
		}
		if(isHoriz == 0)
		{
			VRES->Fill(dphi, dtheta);
			dThetavThetaV->Fill(waistheta, dtheta);
			dPhivPhiV->Fill(waisphi, dphi);
			dPhivThetaV->Fill(waistheta, dphi); 
			dThetavPhiV->Fill(waisphi, dtheta);
			dThetavSNRV->Fill(corr->snr/2., dtheta);
			dPhivSNRV->Fill(corr->snr/2., dphi);
			dPhivTimeV->Fill(sum->eventNumber, dphi);
			dThetavTimeV->Fill(sum->eventNumber, dtheta);
			dPhivDistanceV->Fill(dist, dphi);
			dThetavDistanceV->Fill(dist, dtheta);
		}
	}

	TFile f2(outf.Data(), "RECREATE");
	f2.cd();
	HRES->Write();
	dThetavThetaH->Write();
	dPhivPhiH->Write();
	dThetavPhiH->Write(); 
	dPhivThetaH->Write();
	dThetavSNRH->Write();
	dPhivSNRH->Write();
	dPhivTimeH->Write();
	dThetavTimeH->Write();
	dPhivDistanceH->Write();
	dThetavDistanceH->Write();

	VRES->Write();
	dThetavThetaV->Write();
	dPhivPhiV->Write();
	dThetavPhiV->Write(); 
	dPhivThetaV->Write();
	dThetavSNRV->Write();
	dPhivSNRV->Write();
	dPhivTimeV->Write();
	dThetavTimeV->Write();
	dPhivDistanceV->Write();
	dThetavDistanceV->Write();

	f2.Close();
}



