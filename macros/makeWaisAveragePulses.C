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
#include "RawAnitaHeader.h"
#include "CorrelationSummaryAnita3.h"
#include "TTreeIndex.h"
#include <iostream>
#include <fstream>

#define MAX_ANTENNAS 48

void makeAverageGraphs(TGraph* top, TGraph* mid, TGraph* bot, int center);

void makeWaisAveragePulses(int centAnt)
{
	AnitaGeomTool *agt = AnitaGeomTool::Instance();
	agt->useKurtAnita3Numbers(1);

	TGraph* top = 0;
	TGraph* mid = 0;
	TGraph* bot = 0;
	makeAverageGraphs(top, mid, bot, centAnt);
	
	top->SetName(Form("top%d", centAnt));
	mid->SetName(Form("mid%d", centAnt));
	bot->SetName(Form("bot%d", centAnt));
	
	TFile f(Form("averageGraphs/waisAverage%d.root", centAnt), "RECREATE");
	top->Write();
	mid->Write();
	bot->Write();
	f.Close();

}

void makeAverageGraphs(TGraph* top, TGraph* mid, TGraph* bot, int center)
{
	int i = 0;
	TGraph* gtemptop = 0;
	TGraph* gtempmid = 0;
	TGraph* gtempbot = 0;
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
	UsefulAnitaEvent* uae = 0;

	for (int run = 120; run < 125; run++)
	{
		sprintf(corrName, "corrTrees/run%dCorrTree.root",run);
		corrChain->Add(corrName);
		TString runNum = TString::Itoa(run,10);
		TString headname = datadir + "run" + runNum + "/headFile" + runNum +".root";
		headChain->Add(headname.Data());
		TString evname = datadir + "run" + runNum + "/eventFile" + runNum +".root";
		evChain->Add(evname.Data());
		
		corrChain->SetBranchAddress("corr", &corr);
		corrChain->SetBranchAddress("pol", &pol);
		evChain->SetBranchAddress("event", &event);
		headChain->SetBranchAddress("header", &header);
		
		headChain->BuildIndex("eventNumber");
		ind = (TTreeIndex*) headChain->GetTreeIndex();
	}
	printf("starting one %lld\n", headChain->GetEntries());
	for(int j = 0; j < corrChain->GetEntries(); j++)
	{
		corrChain->GetEntry(j);
		if(corr->centreAntenna != center) continue;
		if(pol == 0) continue;
		int entryIndex = ind->GetEntryNumberWithIndex(corr->eventNumber,0);
		headChain->GetEntry(entryIndex);
		evChain->GetEntry(entryIndex);
		printf("stuck here %d\n", header->eventNumber);	
		uae = new UsefulAnitaEvent(event, cal, header);	
		printf("or here %d\n", j);	
		gtemptop = uae->getGraph(AnitaRing::kTopRing, center, AnitaPol::kHorizontal);
		gtempmid = uae->getGraph(AnitaRing::kMiddleRing, center, AnitaPol::kHorizontal);
		gtempbot = uae->getGraph(AnitaRing::kBottomRing, center, AnitaPol::kHorizontal);
		if(i == 0)
		{
			tops[0] = FFTtools::getInterpolatedGraph(gtemptop, 1/(2.6*16));
			mids[0] = FFTtools::getInterpolatedGraph(gtempmid, 1/(2.6*16));
			bots[0] = FFTtools::getInterpolatedGraph(gtempbot, 1/(2.6*16));
			delete uae;
			i++;
			continue;
		}
		tops[1] = FFTtools::getInterpolatedGraph(gtemptop, 1/(2.6*16));
		mids[1] = FFTtools::getInterpolatedGraph(gtempmid, 1/(2.6*16));
		bots[1] = FFTtools::getInterpolatedGraph(gtempbot, 1/(2.6*16));
		gtemptop = FFTtools::correlateAndAverage(2, tops);
		gtempmid = FFTtools::correlateAndAverage(2, mids);
		gtempbot = FFTtools::correlateAndAverage(2, bots);
		
		tops[0] = (TGraph*) gtemptop->Clone();
		mids[0] = (TGraph*) gtempmid->Clone();
		bots[0] = (TGraph*) gtempbot->Clone();
		for (int k = 0 ; k < tops[0]->GetN(); k++)
		{
			tops[0]->GetY()[k] *= 2;
			mids[0]->GetY()[k] *= 2;
			bots[0]->GetY()[k] *= 2;
		}
		delete uae;
		i++;
	}
	
	top = (TGraph*) tops[0]->Clone();
	mid = (TGraph*) mids[0]->Clone();
	bot = (TGraph*) bots[0]->Clone();

	for (int k = 0 ; k < top->GetN(); k++)
	{
		top->GetY()[k] /= i;
		mid->GetY()[k] /= i;
		bot->GetY()[k] /= i;
	}
	delete[] tops[0];
	delete[] tops[1];
	delete[] mids[0];
	delete[] mids[1];
	delete[] bots[0];
	delete[] bots[1];
	delete gtemptop;
	delete gtempmid;
	delete gtempbot;
}
