#include <stdio.h>
#include <string.h>
#include "TString.h"

void makePhaseCenterCalibfile(const char* name)
{
	TString namestr = name;
	ifstream vFile(Form("TEST/%s_V.txt", namestr.Data()));
	ifstream hFile(Form("TEST/%s_H.txt", namestr.Data()));
	ofstream newfile(Form("phaseCenterFiles/%s.dat", namestr.Data()));

	string line;
	getline(hFile, line);
	while(getline(vFile, line))
	{
		newfile << line <<"\n";
		getline(hFile, line);
		newfile << line <<"\n";
	}
}

void doAll()
{
	makePhaseCenterCalibfile("knownRIndepTLambda");
	makePhaseCenterCalibfile("knownRIndepTTenthLambda");
	makePhaseCenterCalibfile("knownRLambda");
	makePhaseCenterCalibfile("knownRTenthLambda");
}


