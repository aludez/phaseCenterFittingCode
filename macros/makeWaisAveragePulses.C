#include "FFTtools.h"
#include "AnitaConventions.h"
#include "AnitaDataset.h"
#include "PrettyAnitaEvent.h"
#include "UsefulAdu5Pat.h"
#include "AnitaGeomTool.h"
#include "TFile.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "RawAnitaHeader.h"
#include "CorrelationSummaryAnita3.h"
#include "TTreeIndex.h"
#include <iostream>
#include <fstream>

#define MAX_ANTENNAS 48

void makeAverageGraphs(TGraph* top, TGraph* mid, TGraph* bot, int center);

void makeWaisAveragePulses(int centAnt)
{
	TGraph* top = 0;
	TGraph* mid = 0;
	TGraph* bot = 0;
	makeAverageGraphs(top, mid, bot, centAnt);
}

void makeAverageGraphs(TGraph* top, TGraph* mid, TGraph* bot, int center)
{
	int i = 0;
	TGraph* tops[2] = {NULL};	
	TGraph* mids[2] = {NULL};	
	TGraph* bots[2] = {NULL};	
	TChain* corrChain = new TChain("corrTree");
	TChain* headChain = new TChain("headTree");
	TChain* evChain = new TChain("eventTree");
	WaveCalType::WaveCalType_t cal = WaveCalType::kDefault;
	TString datadir = "/project/kicp/avieregg/anitaIV/flight1617/root/";
	char corrName[FILENAME_MAX];
	RawAnitaEvent* event = 0;
	RawAnitaHeader* header = 0;
	TTreeIndex* ind = 0;
	CorrelationSummaryAnita3 * corr = 0;
	int pol = 0;

	for (int run = 120; run < 135; run++)
	{
		sprintf(corrName, "corrTrees/run%dCorrTree.root",run);
		corrChain->Add(corrName);
		TString runNum = TString::Itoa(run,10);
		TString headname = datadir + "run" + runNum + "/headFile" + runNum +".root";
		headChain->Add(headname.Data());
		TString evname = datadir + "run" + runNum + "/eventFile" + runNum +".root";
		evChain->Add(evname.Data());
	}

	corrChain->SetBranchAddress("corr", &corr);
	corrChain->SetBranchAddress("pol", &pol);
	evChain->SetBranchAddress("event", &event);
	headChain->SetBranchAddress("header", &header);
	
	headChain->BuildIndex("eventNumber");
	ind = (TTreeIndex*) headChain->GetTreeIndex();
	printf("starting one %lld\n", headChain->GetEntries());
	
	for(int j = 0; j < corrChain->GetEntries(); j++)
	{
		corrChain->GetEntry(j);
		if(corr->centreAntenna != center) continue;
		if(pol == 1) continue;//0 is V, 1 is H
		int entryIndex = ind->GetEntryNumberWithIndex(corr->eventNumber,0);
		headChain->GetEntry(entryIndex);
		evChain->GetEntry(entryIndex);
		UsefulAnitaEvent uae(event, cal, header);	
		
		TGraph* gtemptop = uae.getGraph(AnitaRing::kTopRing, center, AnitaPol::kVertical);
		TGraph* gtempmid = uae.getGraph(AnitaRing::kMiddleRing, center, AnitaPol::kVertical);
		TGraph* gtempbot = uae.getGraph(AnitaRing::kBottomRing, center, AnitaPol::kVertical);
		
		if(i == 0)
		{
			tops[0] = FFTtools::getInterpolatedGraph(gtemptop, 1/(2.6*16));
			mids[0] = FFTtools::getInterpolatedGraph(gtempmid, 1/(2.6*16));
			bots[0] = FFTtools::getInterpolatedGraph(gtempbot, 1/(2.6*16));
			
			i++;
			delete gtemptop; delete gtempmid; delete gtempbot;
			continue;
		}
		
		tops[1] = FFTtools::getInterpolatedGraph(gtemptop, 1/(2.6*16));
		mids[1] = FFTtools::getInterpolatedGraph(gtempmid, 1/(2.6*16));
		bots[1] = FFTtools::getInterpolatedGraph(gtempbot, 1/(2.6*16));
		
		gtemptop = FFTtools::correlateAndAverage(2, tops);
		gtempmid = FFTtools::correlateAndAverage(2, mids);
		gtempbot = FFTtools::correlateAndAverage(2, bots);
		delete tops[0]; delete mids[0]; delete bots[0]; 
		delete tops[1]; delete mids[1]; delete bots[1]; 
		
		tops[0] = new TGraph(gtemptop->GetN(), gtemptop->GetX(), gtemptop->GetY());
		mids[0] = new TGraph(gtempmid->GetN(), gtempmid->GetX(), gtempmid->GetY());
		bots[0] = new TGraph(gtempbot->GetN(), gtempbot->GetX(), gtempbot->GetY());
		
		for (int k = 0 ; k < tops[0]->GetN(); k++) tops[0]->GetY()[k] *= 2;
		for (int k = 0 ; k < mids[0]->GetN(); k++) mids[0]->GetY()[k] *= 2;
		for (int k = 0 ; k < bots[0]->GetN(); k++) bots[0]->GetY()[k] *= 2;
		
		i++;
		delete gtemptop; delete gtempmid; delete gtempbot;
	}
	
	top = new TGraph(tops[0]->GetN(), tops[0]->GetX(), tops[0]->GetY());
	mid = new TGraph(mids[0]->GetN(), mids[0]->GetX(), mids[0]->GetY());
	bot = new TGraph(bots[0]->GetN(), bots[0]->GetX(), bots[0]->GetY());


	for (int k = 0 ; k < top->GetN(); k++) top->GetY()[k] /= double(i);
	for (int k = 0 ; k < mid->GetN(); k++) mid->GetY()[k] /= double(i);
	for (int k = 0 ; k < bot->GetN(); k++) bot->GetY()[k] /= double(i);
	
	top->SetName(Form("top%d", center));
	mid->SetName(Form("mid%d", center));
	bot->SetName(Form("bot%d", center));
	
	TFile f(Form("averageGraphs/waisAverage%d_V.root", center), "RECREATE");
	top->Write();
	mid->Write();
	bot->Write();
	f.Close();
	
	delete tops[0]; 
	delete mids[0]; 
	delete bots[0]; 
	delete top; delete mid; delete bot;
	delete headChain; delete corrChain; delete evChain;
}
