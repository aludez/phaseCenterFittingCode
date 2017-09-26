#include "AnitaConventions.h"
#include "AnitaEventCalibrator.h"
#include "CorrelationSummaryAnita4.h"
#include "RawAnitaEvent.h"
#include "UsefulAnitaEvent.h"
#include "FilteredAnitaEvent.h"
#include "FilterStrategy.h"
#include "TimedAnitaHeader.h"
#include "PrettyAnitaHk.h"
#include "UsefulAnitaEvent.h"
#include "PrettyAnitaEvent.h"
#include "PrettyAnalysisWaveform.h"
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
#include "SpectrumAverage.h"
#include "BH13Filter.h"
#include "BasicFilters.h"
#include "UCFilters.h"

void makePulserTreeAnalysisWf(int run, int filterType, const char *outDir="corrTrees/");

void makePulserTreeAnalysisWf(int run, int filterType, const char *outDir)
{
	//AnalysisWaveform::enableDebug(1);
	AnitaGeomTool * agt = AnitaGeomTool::Instance();
	agt->usePhotogrammetryNumbers(1);
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
	TString runNum = TString::Itoa(run,10);
	TString evName = datadir + "run" + runNum + "/eventFile" + runNum + ".root";
	TString adu5Name = datadir + "run" + runNum + "/gpsEvent" + runNum + ".root";
	TString headName = datadir + "run" + runNum + "/headFile" + runNum + ".root";

	evChain->Add(evName.Data());
	adu5Chain->Add(adu5Name.Data());
	headChain->Add(headName.Data());

	evChain->SetBranchAddress("event", &event);
	adu5Chain->SetBranchAddress("pat", &pat);
	headChain->SetBranchAddress("header", &header);
	
	TString outName = outDir;
	TString theEnd = Form("CorrTree_%d.root", filterType);
	outName += "run" + TString::Itoa(run,10) + theEnd;
	TFile * outf = new TFile(outName.Data(), "RECREATE");

	CorrelationSummaryAnita4 * theCorr = 0;
	UsefulAdu5Pat * usefulPat = 0;
	Int_t isHoriz = 0;
	TTree * corrTree = new TTree("corrTree", "Tree of Correlation Summaries");
	corrTree->Branch("corr", "CorrelationSummaryAnita4", &theCorr);
	corrTree->Branch("pol", &isHoriz);

	Long64_t maxEntry = evChain->GetEntries();

	Double_t thetaWave, phiWave;

	Double_t ttnsExpected;
	Double_t ttns;
	int ant;
	Int_t adjust = -1000;

	//Double_t deltaT = 1./(2.6*40.);
	Int_t npadfreq = 39;

	headChain->BuildIndex("eventNumber");
	TTreeIndex *ind = (TTreeIndex*) headChain->GetTreeIndex();

	FilterStrategy strategy;
	strategy.addOperation(new UCorrelator::BH13Filter());
	if(filterType == 1) strategy.addOperation(new SimplePassBandFilter(.2, .4));
	if(filterType == 2) strategy.addOperation(new SimplePassBandFilter(.4, .65));
	if(filterType == 3) strategy.addOperation(new SimplePassBandFilter(.65, 1.2));
/*
	static UCorrelator::SpectrumAverageLoader avgldr;
	UCorrelator::SineSubtractFilter* ssf = new UCorrelator::SineSubtractFilter(.1,3);
	ssf->makeAdaptive(&avgldr, 2);
	strategy.addOperation(ssf);
	//strategy.addOperation(new UCorrelator::timePadFilter(1));
//	UCorrelator::AnalysisConfig cfg;
	//UCorrelator::Analyzer analyzer(&cfg);
	//FilterStrategy * strategy = new FilterStrategy(outf);
*/
	for(Long64_t entry=0; entry < maxEntry; entry++)
	{
		evChain->GetEntry(entry);
		
		int entryIndex = ind->GetEntryNumberWithIndex(event->eventNumber, 0);
		if (entryIndex == -1) continue;
		headChain->GetEntry(entryIndex);
		adu5Chain->GetEntry(entryIndex);

		if((header->trigType&1) == 0) continue;
		usefulPat = new UsefulAdu5Pat(pat);
		ttns = header->triggerTimeNs;
		//if (!isLDB)
		//{
		ttnsExpected = usefulPat->getWaisDivideTriggerTimeNs();
		ttnsDelta = Int_t(ttns) - Int_t(ttnsExpected) - adjust;
		double ttnsDeltaH = ttnsDelta + 10e3;
		//}
		//only doing WAIS for now
		/*
		else
		{
			ttnsExpected = usefulPat.getLDBTriggerTimeNs();
			ttnsDelta = Int_t(ttns) - Int_t(ttnsExpected);
		}
		*/

		if((TMath::Abs(ttnsDelta) > cutTimeNs) && (TMath::Abs(ttnsDeltaH) > cutTimeNs)) continue;
		
		if((TMath::Abs(ttnsDelta) < cutTimeNs)){
			isHoriz = 0; //should be V pulse
			pol = AnitaPol::kVertical;
		}
		if((TMath::Abs(ttnsDeltaH) < cutTimeNs)){
			isHoriz = 1; //should be H pulse
			pol = AnitaPol::kHorizontal;
		}
		
		usefulPat->getThetaAndPhiWave(AnitaLocations::LONGITUDE_WAIS_A4, AnitaLocations::LATITUDE_WAIS_A4, AnitaLocations::ALTITUDE_WAIS_A4, thetaWave, phiWave);
		UsefulAnitaEvent * uae = new UsefulAnitaEvent(event, WaveCalType::kDefault, header);
//want to add a pointing cut but havent gotten it to work
/*			
		AnitaEventSummary * aes = new AnitaEventSummary(header, usefulPat);
		int goodflag = 0;
		UsefulAnitaEvent * uae = new UsefulAnitaEvent(event, WaveCalType::kDefault, header);
		FilteredAnitaEvent * fae = new FilteredAnitaEvent(uae, strategy, pat, header);
		analyzer.analyze(fae, aes);	
		for (int peakNum = 0; peakNum<2; peakNum++)
		{
			AnitaEventSummary::PointingHypothesis ph = aes->peak[pol][peakNum];
			double phiDiff = ph.phi - phiWave*TMath::RadToDeg();
			double thetaDiff = ph.theta - thetaWave*TMath::RadToDeg();
			if(fabs(phiDiff) < 1 && fabs(thetaDiff) < 1) goodflag = 1;
			printf("phi diff = %g,theta diff = %g\n", fabs(phiDiff), fabs(thetaDiff));
		}
		if (goodflag == 0) continue;
*/
		PrettyAnalysisWaveform pae(uae, &strategy, pat, header);
		//else usefulPat.getThetaAndPhiWaveLDB(thetaWave, phiWave);
		ant = agt->getTopAntNearestPhiWave(phiWave, pol);
		theCorr = pae.getCorrelationSummaryAnita4(ant, pol, npadfreq, uae->eventNumber);
		corrTree->Fill();
		delete theCorr;
		delete usefulPat;
		delete uae;
		//delete aes;
	}
	outf->cd();
	corrTree->Write();
	outf->Close();
	delete evChain;
	delete adu5Chain;
	delete headChain;
}

void doAll()
{
	//for (int i = 125; i < 126; i++) 
	{
		makePulserTreeAnalysisWf(125, 0);
		//printf("%d done\n", i);
	}
}
