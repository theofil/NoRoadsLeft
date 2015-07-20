/*
compile this puppy with

g++ flash_remove_duplicates.C -o flash_remove_duplicates.exec `root-config --glibs --cflags`

### original author: Marco-Andrea Buchmann

*/
#include <iostream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <time.h>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <TTree.h>
#include <TFile.h>
#include <TSystem.h>

const int rMax=30;
const int jMax=30;
const int metMax=30;

using namespace std;

class femtoEventTree
{
// this is for a binary tree
// http://en.wikipedia.org/wiki/Binary_tree
  UInt_t RunNum;
  UInt_t lumi;
  ULong64_t EventNum;
  femtoEventTree *Left, *Right;
  
  int compare(UInt_t RunNum, ULong64_t EventNum, UInt_t lumi);//-1 if smaller, 0 if equal, 1 if larger
  
  public:
    bool insert(UInt_t RunNum, ULong64_t EventNum, UInt_t lumi);// 1: has been inserted, 0: was already there!
    femtoEventTree(UInt_t RunNum, ULong64_t EventNum, UInt_t lumi);
    ~femtoEventTree();
};

int femtoEventTree::compare(UInt_t RunNum, ULong64_t EventNum,UInt_t lumi)
{
  if (RunNum < this->RunNum) return -1;
  if (RunNum > this->RunNum) return 1;
//if we're dealing with two events in the same run
	if (lumi < this->lumi) return -1;
	if (lumi > this->lumi) return 1;
//if we're dealing with two events in the same lumi section
	if (EventNum < this->EventNum) return -1;
	if (EventNum > this->EventNum) return 1;
// we're dealing with a duplicate
  return 0;
}

bool femtoEventTree::insert(UInt_t RunNum, ULong64_t EventNum, UInt_t lumi)
{
  int compare_result = compare(RunNum, EventNum,lumi);
  if (compare_result==0) return false;
  femtoEventTree **Child;
  Child=(compare_result==-1)?&Left:&Right;
  
    if(*Child==NULL) // doesn't exist yet (we're at the bottom of the tree -- happy planting!)
    {
//      std::cout << this->EventNum << " " << this->RunNum << " " << EventNum << " " << RunNum << " " << compare_result << std::endl;
      *Child = new femtoEventTree(EventNum,RunNum,lumi);
      return true;
    }
    else
    {//does already exist
      return (*Child)->insert(EventNum,RunNum,lumi);
    }
}

femtoEventTree::femtoEventTree(UInt_t RunNum, ULong64_t EventNum, UInt_t lumi)
{
  this->RunNum=RunNum;
  this->EventNum=EventNum;
	this->lumi=lumi;
  Left=NULL;
  Right=NULL;
}

femtoEventTree::~femtoEventTree(){
if(Right != NULL) // not at the bottom
delete Right;
if(Left != NULL) // not at the bottom
delete Left;
}

struct femtoEvent
{
    ULong64_t eventNum;
    UInt_t runNum;
    UInt_t lumi;
} ;
femtoEvent fEvent;


void process_directory(TString directory, vector<string> &files)
{
	DIR *dp;
	struct dirent *ep;
	
	dp = opendir (directory);
	if (dp != NULL)
	{
		while (ep = readdir (dp))
		{
			string filename=(string)ep->d_name;
			if(filename.find(".root")!=-1) 
			{
				files.push_back(string(directory)+filename);
			}
		}
		(void) closedir (dp);
	}
	else
		perror ("Couldn't open the directory");
}


//________________________________________________________________________________________
// Print out usage
void usage( int status = 0 ) {
	std::cout << "Usage: removeDuplicates [-o filename] file1 file2 file3 ..." << std::endl;
	std::cout << "  where:" << std::endl;
	std::cout << "     filename is the merged output filename         " << std::endl;
	std::cout << "               default is 0 (quiet mode)            " << std::endl;
	std::cout << "     files to merge and clean        " << std::endl;
	std::cout << std::endl;
	exit(status);
}

int main(int argc, char * argv[])
{
    TString outputFileName;
	TString directory;
	bool directoryset=false;
	// Parse options
    char ch;
    while ((ch = getopt(argc, argv, "o:v:lh?d:")) != -1 ) {
		switch (ch) {
			case 'o': outputFileName = std::string(optarg); break;
			case '?':
			case 'h': usage(0); break;
			case 'd': directory=std::string(optarg); directoryset=true;break;
			default:
				std::cerr << "*** Error: unknown option " << optarg << std::endl;
				usage(-1);
		}
    }
	argc -= optind;
    argv += optind;
	vector<string> files;
	if(directoryset) process_directory(directory,files);
	
	for(int iarg=0;iarg<argc;iarg++)
		files.push_back(argv[iarg]);
	
    // Check arguments
    if( (files.size()<1) || outputFileName.Length()==0 ) { usage(-1); }
	
	
	/*
	 step 1 : merge all files into one
	 step 2 : open the one file and build a new one, while rejecting all duplicates.
	 */
	
	
	/// STEP 1
	stringstream mergecommand;
	mergecommand << "hadd -f allincluded.root";
	for(int ifile=0;ifile<files.size();ifile++)
	{
		mergecommand << " " << files[ifile];
	}
	if(files.size()>1) {
		cout << "Fusing all files using the command: " << mergecommand.str() << endl;
		gSystem->Exec(mergecommand.str().c_str());
	}
	else {
		stringstream alternativecommand;
		alternativecommand << "cp " << files[0] << " allincluded.root";
		gSystem->Exec(alternativecommand.str().c_str());
	}
		
	// STEP 2
	TFile* f1 = new TFile("allincluded.root");
        TFile* newFile = new TFile(outputFileName,"RECREATE");
        newFile->mkdir("demo")->cd();
	
	TTree* t1 = (TTree*)f1->Get("demo/events");
	TTree* newTree = t1->CloneTree(0);
	
	UInt_t runNum;
        UInt_t lumi;
	ULong64_t eventNum;
	t1->SetBranchAddress("eventNum",&eventNum);
	t1->SetBranchAddress("runNum",&runNum);
	t1->SetBranchAddress("lumi",&lumi);
	
	Int_t nentries = (Int_t)t1->GetEntries();
    Int_t duplicates = 0;
    int freq = nentries/100;
	
	femtoEventTree *binarytree=NULL;

	std::cout << "Processing ...   0% " << std::flush;
    for (Int_t i=0;i<nentries; i++) {
        if ( freq>0 && !(i%freq) ) { // Counter
			std::cout << "\b\b\b\b\b" << std::setprecision(0) << std::setw(3) << std::fixed 
			<< i/static_cast<double>(nentries)*100 << "% " << std::flush;
			std::cout << std::setprecision(4);
        }
        t1->GetEntry(i);
		if (i==0) {
			binarytree = new femtoEventTree(runNum,eventNum,lumi);
			newTree->Fill();
		}
		else {
			bool is_duplicate = !(binarytree->insert(runNum,eventNum,lumi));
			if (is_duplicate) {
				std::cout << "we've found a duplicate (our number " << duplicates << "): run number " << runNum << ",lumi " << lumi << " and event number " << eventNum << std::endl;
				++duplicates;
			}
			else {
				newTree->Fill();
			}
		}
		
    }
    std::cout << std::endl;
    newTree->AutoSave();
    cout << "We have found " << duplicates << " duplicates out of " << nentries << endl; 
    gSystem->Exec("rm allincluded.root");
     
    newFile->Close();

	return 0;
}

