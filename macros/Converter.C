class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;

#include <math.h>


void Converter(int nevts = 0, TString prefix = "muon_150k", TString out = "TREEOUT")
{
 	// *** some variables
	int i=0;
	gStyle->SetOptFit(1011);

	
	// *** the output file for FairRunAna
	TString OutFile="out_dummy.root";  
	
					
	// *** the files coming from the simulation
	TString inPidFile  = prefix+"_pid.root";    // this file contains the PndPidCandidates and McTruth
	TString inParFile  = prefix+"_par.root";
	
	// *** PID table with selection thresholds; can be modified by the user
	TString pidParFile = TString(gSystem->Getenv("VMCWORKDIR"))+"/macro/params/all.par";	
	
	// *** initialization
	FairLogger::GetLogger()->SetLogToFile(kFALSE);
	FairRunAna* fRun = new FairRunAna();
	fRun->AddFriend("sim_complete.root");
	FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
	fRun->SetSource(new FairFileSource("pid_complete.root"));
	
	// *** setup parameter database 	
	FairParRootFileIo* parIO = new FairParRootFileIo();
	parIO->open("simparams.root");
	FairParAsciiFileIo* parIOPid = new FairParAsciiFileIo();
	parIOPid->open(pidParFile.Data(),"in");
	
	rtdb->setFirstInput(parIO);
	rtdb->setSecondInput(parIOPid);
	rtdb->setOutput(parIO);  
	
	fRun->SetOutputFile(OutFile);
	fRun->Init(); 

	// *** the data reader object
	PndAnalysis* theAnalysis = new PndAnalysis();
	if (nevts==0) nevts= theAnalysis->GetEntries();
	
	// *** RhoCandLists for the analysis
	RhoCandList eminus, piminus, muminus, kminus, pminus, all;
	
//  ====================================================== Variables Begin ====================================================== //	
	// Basic Variables
     float momentumx, momentumy, momentumz, momentum, energy;
	 float positionx, positiony, positionz, position;
	 float fhitx, fhity, fhitz, lhitx, lhity,lhitz, fhit, lhit;
     int charge, MCindex, Trackindex;
	
	// Detector Variables
	
    // == MVD ==
     float MvdDEDX;
     int  MvdHits;
    
    // == STT ==
     float SttMeanDEDX;
    int    SttHits;
    
    // == GEM ==
      int   GemHits;
    
    // == TOF == 
     float TofStopTime;
     float TofM2;
     float TofTrackLength;
     float TofQuality;
        int TofIndex; 
      float TofBeta; 
    

    // == Barrel DIRC ==
    float  DrcThetaC;
     float DrcThetaCErr;
      float DrcQuality;
       int   DrcNumberOfPhotons;
      int   DrcIndex;    
    
    // == Disc DIRC ==
     float DiscThetaC;
     float DiscThetaCErr;
      float DiscQuality;
       int   DiscNumberOfPhotons;
     int    DiscIndex;
    
    // == RICH ==
     float RichThetaC;
    float RichThetaCErr;
     float RichQuality;
        int  RichNumberOfPhotons;
      int    RichIndex;
    
    // == EMC ==
     float EmcRawEnergy;
    float EmcCalEnergy;
     float EmcQuality;
      int   EmcNumberOfCrystals;
       int   EmcNumberOfBumps;
      int   EmcModule;
       int   EmcIndex;
    
       // EMC Cluster properties
        double  EmcZ20;
        double  EmcZ53;
         double EmcLat;
         double EmcE1;
         double EmcE9;
         double EmcE25;
    
    // == Muons ==
	 float MuoProbability;
	 float MuoQuality;
	 float MuoIron;
	 float MuoMomentumIn;
	  int   MuoNumberOfLayers;
	  int   MuoModule;
	 int    MuoHits;
	  int   MuoIndex;
	
	// == Tracking ==
	  int   DegreesOfFreedom;
	   int   FitStatus;
	 float ChiSquared;	
//  ====================================================== Variables end   ====================================================== //

//  ====================================================== Branches Begin   ====================================================== //
TFile f("background_tree.root","recreate");
TTree t1("t1","matched tree");	


t1.Branch("momentumx",&momentumx,"momentumx/F");
t1.Branch("momentumy",&momentumy,"momentumy/F");
t1.Branch("momentumz",&momentumz,"momentumz/F");
t1.Branch("momentum",&momentum,"momentum/F");
t1.Branch("energy",&energy,"energy/F");
t1.Branch("positionx",&positionx,"positionx/F");
t1.Branch("positiony",&positiony,"positiony/F");
t1.Branch("positionz",&positionz,"positionz/F");
t1.Branch("position",&position,"position/F");
t1.Branch("charge",&charge,"charge/I");
t1.Branch("MCindex",&MCindex,"MCindex/I");
t1.Branch("Trackindex",&Trackindex,"Trackindex/I");

t1.Branch("MvdDEDX",&MvdDEDX,"MvdDEDX/F");
t1.Branch("MvdHits",&MvdHits,"MvdHits/I");
t1.Branch("SttMeanDEDX",&SttMeanDEDX,"SttMeanDEDX/F");
t1.Branch("SttHits",&SttHits,"SttHits/I");
t1.Branch("GemHits",&GemHits,"GemHits/I");
t1.Branch("TofStopTime",&TofStopTime,"TofStopTime/F");
t1.Branch("TofM2",&TofM2,"TofM2/F");
t1.Branch("TofTrackLength",&TofTrackLength,"TofTrackLength/F");
t1.Branch("TofQuality",&TofQuality,"TofQuality/F");
t1.Branch("TofIndex",&TofIndex,"TofIndex/I");
t1.Branch("TofBeta",&TofBeta,"TofBeta/F");

t1.Branch("DrcThetaC",&DrcThetaC,"DrcThetaC/F");
t1.Branch("DrcThetaCErr",&DrcThetaCErr,"DrcThetaCErr/F");
t1.Branch("DrcQuality",&DrcQuality,"DrcQuality/F");
t1.Branch("DrcNumberOfPhotons",&DrcNumberOfPhotons,"DrcNumberOfPhotons/I");
t1.Branch("DrcIndex",&DrcIndex,"DrcIndex/I");

t1.Branch("DiscThetaC",&DiscThetaC,"DiscThetaC/F");
t1.Branch("DiscThetaCErr",&DiscThetaCErr,"DiscThetaCErr/F");
t1.Branch("DiscQuality",&DiscQuality,"DiscQuality/F");
t1.Branch("DiscNumberOfPhotons",&DiscNumberOfPhotons,"DiscNumberOfPhotons/I");
t1.Branch("DiscIndex",&DiscIndex,"DiscIndex/I");

t1.Branch("RichThetaC",&RichThetaC,"RichThetaC/F");
t1.Branch("RichThetaCErr",&RichThetaCErr,"RichThetaCErr/F");
t1.Branch("RichQuality",&RichQuality,"RichQuality/F");
t1.Branch("RichNumberOfPhotons",&RichNumberOfPhotons,"RichNumberOfPhotons/I");
t1.Branch("RichIndex",&RichIndex,"RichIndex/I");

t1.Branch("EmcRawEnergy",&EmcRawEnergy,"EmcRawEnergy/F");
t1.Branch("EmcCalEnergy",&EmcCalEnergy,"EmcCalEnergy/F");
t1.Branch("EmcQuality",&EmcQuality,"EmcQuality/F");
t1.Branch("EmcNumberOfCrystals",&EmcNumberOfCrystals,"EmcNumberOfCrystals/I");
t1.Branch("EmcNumberOfBumps",&EmcNumberOfBumps,"EmcNumberOfBumps/I");
t1.Branch("EmcModule",&EmcModule,"EmcModule/I");
t1.Branch("EmcIndex",&EmcIndex,"EmcIndex/I");
t1.Branch("EmcZ20",&EmcZ20,"EmcZ20/D");
t1.Branch("EmcZ53",&EmcZ53,"EmcZ53/D");
t1.Branch("EmcLat",&EmcLat,"EmcLat/D");
t1.Branch("EmcE1",&EmcE1,"EmcE1/D");
t1.Branch("EmcE9",&EmcE9,"EmcE9/D");
t1.Branch("EmcE25",&EmcE25,"EmcE25/D");

t1.Branch("MuoProbability",&MuoProbability,"MuoProbability/F");
t1.Branch("MuoQuality",&MuoQuality,"MuoQuality/F");
t1.Branch("MuoIron",&MuoIron,"MuoIron/F");
t1.Branch("MuoMomentumIn",&MuoMomentumIn,"MuoMomentumIn/F");
t1.Branch("MuoNumberOfLayers",&MuoNumberOfLayers,"MuoNumberOfLayers/I");
t1.Branch("MuoModule",&MuoModule,"MuoModule/I");
t1.Branch("MuoHits",&MuoHits,"MuoHits/I");
t1.Branch("MuoIndex",&MuoIndex,"MuoIndex/I");

t1.Branch("DegreesOfFreedom",&DegreesOfFreedom,"DegreesOfFreedom/I");
t1.Branch("FitStatus",&FitStatus,"FitStatus/I");
t1.Branch("ChiSquared",&ChiSquared,"ChiSquared/F");
//  ====================================================== Branches end   ====================================================== //
	
	// ***
	// the event loop
	// ***
	int I=0;
	//PndAnaPidSelector kaonSel("KaonSelector");
	//kaonSel.SetSelection("KaonLoose");
	while (theAnalysis->GetEvent() && i++<nevts)
	{
		if ((i%100)==0) cout<<"evt " << i << endl;
		
		// *** Select with no PID info ('All'); type and mass are set
		theAnalysis->FillList( eminus, "ElectronAllMinus");
		theAnalysis->FillList( piminus, "PionAllMinus");
		theAnalysis->FillList( muminus, "MuonAllMinus");
		theAnalysis->FillList( kminus, "KaonAllMinus");
		theAnalysis->FillList( pminus, "ProtonAllMinus");
		theAnalysis->FillList( all, "Charged");
		    
		for (int j=0; j<all.GetLength(); ++j)
		{
			//if (theAnalysis->McTruthMatch(all[j]))
			//{
				cout<<"PdgCode() is: "<<all[j]->PdgCode()<<endl;
				PndPidCandidate *myCand = (PndPidCandidate*) all[j]->GetRecoCandidate();                

                int index = myCand->GetMcIndex();
				if (index != 0) continue;
				
		        momentumx=myCand->GetMomentum().X();
		        momentumy=myCand->GetMomentum().Y();
		        momentumz=myCand->GetMomentum().Z();
		        
		        double mom = ((myCand->GetMomentum().X()*myCand->GetMomentum().X()) + (myCand->GetMomentum().Y()*myCand->GetMomentum().Y()) + (myCand->GetMomentum().Z()*myCand->GetMomentum().Z()));		        
		        momentum=sqrt(mom);
		        
		        energy=myCand->GetEnergy();
		        positionx=myCand->GetPosition().X();	
		        positiony=myCand->GetPosition().Y();
		        positionz=myCand->GetPosition().Z();	
		        
		        float pos = ((myCand->GetPosition().X()*myCand->GetPosition().X()) + (myCand->GetPosition().Y()*myCand->GetPosition().Y()) + (myCand->GetPosition().Z()*myCand->GetPosition().Z()));				        
		        position=sqrt(pos);
		        
		        charge=myCand->GetCharge();
		        MCindex=myCand->GetMcIndex();
		        Trackindex=myCand->GetTrackIndex();
		        
		        // == MVD ==
		        MvdDEDX=myCand->GetMvdDEDX();
		        MvdHits=myCand->GetMvdHits();		  
		        // == STT ==
		        SttMeanDEDX=myCand->GetSttMeanDEDX();
		        SttHits=myCand->GetSttHits();
		        // == GEM ==
		        GemHits=myCand->GetGemHits();   
		        // == TOF ==
		        TofStopTime=myCand->GetTofStopTime();
		        TofM2=myCand->GetTofM2();
		        TofTrackLength=myCand->GetTofTrackLength();
		        TofQuality=myCand->GetTofQuality();
		        TofIndex=myCand->GetTofIndex();
		        TofBeta=myCand->GetTofBeta();
		        // == Barrel DIRC ==
                DrcThetaC=myCand->GetDrcThetaC();
                DrcThetaCErr=myCand->GetDrcThetaCErr();
                DrcQuality=myCand->GetDrcQuality();
                DrcNumberOfPhotons=myCand->GetDrcNumberOfPhotons();
                DrcIndex=myCand->GetDrcIndex();    
                // == Disc DIRC ==
                DiscThetaC=myCand->GetDiscThetaC();
                DiscThetaCErr=myCand->GetDiscThetaCErr();
                DiscQuality=myCand->GetDiscQuality();
                DiscNumberOfPhotons=myCand->GetDiscNumberOfPhotons();
                DiscIndex=myCand->GetDiscIndex();
                // == RICH ==
                RichThetaC=myCand->GetRichThetaC();
                RichThetaCErr=myCand->GetRichThetaCErr();
                RichQuality=myCand->GetRichQuality();
                RichNumberOfPhotons=myCand->GetRichNumberOfPhotons();
                RichIndex=myCand->GetRichIndex();
                // == EMC ==
                EmcRawEnergy=myCand->GetEmcRawEnergy();
                EmcCalEnergy=myCand->GetEmcCalEnergy();
                EmcQuality=myCand->GetEmcQuality();
                EmcNumberOfCrystals=myCand->GetEmcNumberOfCrystals();
                EmcNumberOfBumps=myCand->GetEmcNumberOfBumps();
                EmcModule=myCand->GetEmcModule();
                EmcIndex=myCand->GetEmcIndex();
                  // EMC Cluster properties
                  EmcZ20=myCand->GetEmcClusterZ20();
                  EmcZ53=myCand->GetEmcClusterZ53();
                  EmcLat=myCand->GetEmcClusterLat();
                  EmcE1=myCand->GetEmcClusterE1();
                  EmcE9=myCand->GetEmcClusterE9();
                  EmcE25=myCand->GetEmcClusterE25();
               // == Muons ==
               MuoProbability=myCand->GetMuoProbability();
               MuoQuality=myCand->GetMuoQuality();
               MuoIron=myCand->GetMuoIron();
               MuoMomentumIn=myCand->GetMuoMomentumIn();
	           MuoNumberOfLayers=myCand->GetMuoNumberOfLayers();
	           MuoModule=myCand->GetMuoModule();
	           MuoHits=myCand->GetMuoHits();
	           MuoIndex=myCand->GetMuoIndex();
	           // == Tracking ==
               DegreesOfFreedom=myCand->GetDegreesOfFreedom();
               FitStatus=myCand->GetFitStatus();
               ChiSquared=myCand->GetChiSquared();	   
               t1.Fill();               
                
			//}
		}
		    
	}  // event loop
	t1.Write();
} // macro end

