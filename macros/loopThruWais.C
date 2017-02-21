#include "AnitaConventions.h"
#include "AnitaEventCalibrator.h"
#include "CorrelationSummaryAnita3.h"
#include "RawAnitaEvent.h"
#include "UsefulAnitaEvent.h"
#include "FilteredAnitaEvent.h"
#include "FilterStrategy.h"
#include "TimedAnitaHeader.h"
#include "PrettyAnitaHk.h"
#include "UsefulAnitaEvent.h"
#include "PrettyAnitaEvent.h"
#include "AnitaGeomTool.h"
#include "Adu5Pat.h"
#include "UsefulAdu5Pat.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "AnitaDataset.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
#include "TTreeIndex.h"
#include "Util.h"
#include "AnitaEventSummary.h"
#include "Analyzer.h"
#include "AnalysisConfig.h"

void loopThruWais(int runStart, int runEnd)
{
	AnitaGeomTool * agt = AnitaGeomTool::Instance();
	agt->useKurtAnita3Numbers(1);
	AnitaEventCalibrator * cal = AnitaEventCalibrator::Instance();
	for(Int_t surf = 0; surf < NUM_SURF; surf++)
	{
		for(Int_t chan = 0; chan < NUM_CHAN; chan++)
		{
			cal->relativePhaseCenterToAmpaDelays[surf][chan] = 0;
		}
	}

	AnitaPol::AnitaPol_t pol;

	Double_t sourceLat, sourceLon, sourceAlt, timeOffset;
	Int_t cutTimeNs = 1.2e3;
	bool isLDB = false;

	Int_t delay = 0;
	Int_t constdelay = 500;
	Int_t ttnsDelta;

	RawAnitaEvent * event = 0;
	RawAnitaHeader * header = 0;
	Adu5Pat * pat = 0;

	TChain *evChain = new TChain("eventTree");
	TChain *adu5Chain = new TChain("adu5PatTree");
	TChain *headChain = new TChain("headTree");

	TString datadir = "/project/kicp/avieregg/anitaIV/flight1617/root/";
	for (int run = runStart; run < runEnd; run++)
	{
		TString runNum = TString::Itoa(run,10);
		TString evName = datadir + "run" + runNum + "/eventFile" + runNum + ".root";
		TString adu5Name = datadir + "run" + runNum + "/gpsEvent" + runNum + ".root";
		TString headName = datadir + "run" + runNum + "/headFile" + runNum + ".root";

		evChain->Add(evName.Data());
		adu5Chain->Add(adu5Name.Data());
		headChain->Add(headName.Data());
	}
	evChain->SetBranchAddress("event", &event);
	adu5Chain->SetBranchAddress("pat", &pat);
	headChain->SetBranchAddress("header", &header);
	
	CorrelationSummaryAnita3 * theCorr = 0;
	UsefulAdu5Pat * usefulPat = 0;
	Int_t isHoriz = 0;
	TTree * corrTree = new TTree("uPatTree", "Tree of Usefulpats");
	corrTree->Branch("usefulPat", "UsefulAdu5Pat", &usefulPat);

	Long64_t maxEntry = evChain->GetEntries();

	TGraph* g = new TGraph(maxEntry);
	TH1D* h = new TH1D("h", "h", 1000, -50000, 50000);

	Double_t thetaWave, phiWave;

	Double_t ttnsExpected;
	Double_t ttns;
	int ant;
	
	Double_t deltaT = 1./(2.6*40.);

	headChain->BuildIndex("eventNumber");
	TTreeIndex *ind = (TTreeIndex*) headChain->GetTreeIndex();

	printf("ev = %lld, adu5 = %lld, head = %lld\n", evChain->GetEntries(), adu5Chain->GetEntries(), headChain->GetEntries());



	for(Long64_t entry=0; entry < maxEntry; entry++)
	{
		evChain->GetEntry(entry);
		
		int entryIndex = ind->GetEntryNumberWithIndex(event->eventNumber, 0);
		if (entryIndex == -1) continue;
		headChain->GetEntry(entryIndex);
		adu5Chain->GetEntry(entryIndex);
		if ((header->trigType&1) == 0) continue;

		usefulPat = new UsefulAdu5Pat(pat);
		ttns = header->triggerTimeNs;
		ttnsExpected = usefulPat->getWaisDivideTriggerTimeNs();
		ttnsDelta = Int_t(ttns) - Int_t(ttnsExpected);
		double distToWais = usefulPat->getDistanceFromSource(AnitaLocations::LATITUDE_WAIS, AnitaLocations::LONGITUDE_WAIS, AnitaLocations::ALTITUDE_WAIS);
		double timeToWais = 1e9 * distToWais/299792458;
		double ttnsDeltaH = ttnsDelta + 10e3;
		g->SetPoint(entry, event->eventNumber, ttns - timeToWais);
		h->Fill(ttns-timeToWais);
		corrTree->Fill();
		delete theCorr;
		delete usefulPat;
	}
	TFile f(Form("outFile%dto%d.root", runStart, runEnd), "RECREATE");
	g->Write();
	h->Write();
	h->Draw();
	//g->Draw("ap");
	delete evChain;
	delete adu5Chain;
	delete headChain;
}

