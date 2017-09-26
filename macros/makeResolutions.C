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

void makeResolutions(const char * whichFile, int run)
{
	TString inName = whichFile;
	FFTtools::loadWisdom("wisdom.dat");

	TString pointDir = "pointingGraphs/";
	TString psDir = "phaseCenterFiles/";
	TString psName = inName;
	TString blinding = "PolarityBlinded";
	TString dat = ".dat";
	TString fEnd = ".root";

	TString phaseCenterFile = psDir + psName + dat;
	std::ifstream phaseCenterIn(phaseCenterFile.Data());
	TString outName = pointDir + psName + TString::Itoa(run,10) + fEnd;
		
	AnitaGeomTool* agt = AnitaGeomTool::Instance();
	agt->usePhotogrammetryNumbers(1);
	
	Double_t extraCableDelaysV[48];
	Double_t deltaRV[48];
	Double_t deltaPhiV[48];
	Double_t deltaZV[48];
	
	Double_t extraCableDelaysH[48];
	Double_t deltaRH[48];
	Double_t deltaPhiH[48];
	Double_t deltaZH[48];

	Int_t ant;
	Int_t pols;
	Double_t dr, dphi, dz, dt;

	int lineN = 1;

	dt = 0;
	string line;
	getline(phaseCenterIn, line);
	while(phaseCenterIn >> ant >> pols >> dr >> dphi >> dz >> dt)
	//while(phaseCenterIn >> ant >> pols >> dr >> dphi >> dz)
	{
		if(lineN%2)
		{
			extraCableDelaysH[ant] = dt;
			deltaRH[ant] = dr;
			deltaPhiH[ant] = dphi * TMath::Pi()/180.;
			deltaZH[ant] = dz;
		} 
		else 
		{
			extraCableDelaysV[ant] = dt;
			deltaRV[ant] = dr;
			deltaPhiV[ant] = dphi * TMath::Pi()/180.;
			deltaZV[ant] = dz;
		}
		lineN++;
	}

	AnitaPol::AnitaPol_t tempPol;
	AnitaEventCalibrator* cal = AnitaEventCalibrator::Instance();
	for(Int_t surf = 0; surf < NUM_SURF; surf++)
	{
		for(Int_t chan = 0; chan < NUM_CHAN-1; chan++)
		{
			agt->getAntPolFromSurfChan(surf, chan, ant, tempPol);
			if(ant!=-1){
				if(tempPol == 0) cal->relativePhaseCenterToAmpaDelays[surf][chan] = extraCableDelaysH[ant];
				if(tempPol == 1) cal->relativePhaseCenterToAmpaDelays[surf][chan] = extraCableDelaysV[ant];
			}
		}
	}

	for(Int_t iant = 0; iant < 48; iant++)
	{
		agt->deltaRPhaseCentre[iant][1] = deltaRV[iant];
		agt->deltaRPhaseCentre[iant][0] = deltaRH[iant];
		agt->deltaPhiPhaseCentre[iant][1] = deltaPhiV[iant];
		agt->deltaPhiPhaseCentre[iant][0] = deltaPhiH[iant];
		agt->deltaZPhaseCentre[iant][1] = deltaZV[iant];
		agt->deltaZPhaseCentre[iant][0] = deltaZH[iant];
	}

	agt->addPhaseCenters();
	agt->usePhotogrammetryNumbers(0);

	const UCorrelator::AntennaPositions* ap = UCorrelator::AntennaPositions::instance(4,agt);
	
	TTree* corrTree = 0;
	CorrelationSummaryAnita4* corr = 0;
	AnitaEventSummary* sum = new AnitaEventSummary;

	UCorrelator::AnalysisConfig cfg;
	//cfg.enable_group_delay = false;
	cfg.nmaxima = 3;
	UCorrelator::Analyzer analyzer(&cfg);
	uint64_t lshift(1ul<<45);
	analyzer.setDisallowedAntennas(0, lshift);
	
	AnitaPol::AnitaPol_t kpol;
	int isHoriz = 0;
	AnitaDataset dd(run);
	dd.firstEvent();
	double firstEvent = double(dd.header()->eventNumber);
	dd.lastEvent();
	double lastEvent = double(dd.header()->eventNumber);

	FilterStrategy strategy;
	strategy.addOperation(new UCorrelator::BH13Filter());
	
	char corrName[FILENAME_MAX];
	sprintf(corrName, "corrTrees/run%dCorrTree.root", run);
	TFile f(corrName);

	corrTree = (TTree*) f.Get("corrTree");
	corrTree->SetBranchAddress("corr",&corr);
	corrTree->SetBranchAddress("pol",&isHoriz);
		
	AnitaDataset d(run);
//	d.setStrategy(AnitaDataset::BlindingStrategy::kRandomizePolarity);
	TFile f2(outName.Data(), "RECREATE");
	TTree* tree = new TTree("resolution", "resolution");
	tree->Branch("summary", &sum);
	
	for (int j = 0; j < corrTree->GetEntries(); j ++)
	{
		corrTree->GetEntry(j);
		if (isHoriz == 1) kpol = AnitaPol::kHorizontal;
		if (isHoriz == 0) kpol = AnitaPol::kVertical; 
		int event = corr->eventNumber;
		d.getEvent(event);

		FilteredAnitaEvent ev(d.useful(), &strategy, d.gps(), d.header());
		analyzer.analyze(&ev, sum);
		
		tree->Fill();
	}
	f2.cd();
	tree->Write();
	f2.Close();
	FFTtools::saveWisdom("wisdom.dat");
}



