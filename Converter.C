class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;

#include <math.h>
// *** routine to only keep PID matched candidates in list
int SelectTruePid(PndAnalysis *ana, RhoCandList &l)
{
        int removed = 0;

        for (int ii=l.GetLength()-1;ii>=0;--ii)
        {
                if ( !(ana->McTruthMatch(l[ii])) )
                {
                        l.Remove(l[ii]);
                        removed++;
                }
        }

        return removed;
}


void Converter(int nevts = 0, TString prefix = "", TString PreFix = "", TString x="", TString out = "TREEOUT")
{
     // *** some variables
    int i=0;
    gStyle->SetOptFit(1011);

   
    // *** the output file for FairRunAna
    TString OutFile="out_dummy.root";

    // *** the files coming from the simulation
    TString inPidFile  = prefix+"_pid.root";    // this file contains the PndPidCandidates and McTruth
    TString inParFile  = prefix+"_par.root";
    TString inSimFile  = prefix+"_sim.root";


        TFile *infile = new TFile(inPidFile, "read");
        TTree *t = (TTree*)infile->Get("pndsim");
       // Int_t entries = (Int_t)t->GetEntries();

    // *** PID table with selection thresholds; can be modified by the user
    TString pidParFile = TString(gSystem->Getenv("VMCWORKDIR"))+"/macro/params/all.par";   
   
    // *** initialization
    FairLogger::GetLogger()->SetLogToFile(kFALSE);
    FairRunAna* fRun = new FairRunAna();
    //fRun->AddFriend(inSimFile);
    FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
    fRun->SetSource(new FairFileSource(inPidFile));
   
    // *** setup parameter database    
    //FairParRootFileIo* parIO = new FairParRootFileIo();
    //parIO->open(inParFile);
    //FairParAsciiFileIo* parIOPid = new FairParAsciiFileIo();
    //parIOPid->open(pidParFile.Data(),"in");
   
    //rtdb->setFirstInput(parIO);
    //rtdb->setSecondInput(parIOPid);
    //rtdb->setOutput(parIO); 
   
    fRun->SetOutputFile(OutFile);
    fRun->Init();

    // *** the data reader object
    PndAnalysis* theAnalysis = new PndAnalysis();
    if (nevts==0) nevts= theAnalysis->GetEntries();
   
    // *** RhoCandLists for the analysis
    RhoCandList eminus, piminus, muminus, kminus, pminus, all;
   
//  ====================================================== Variables Begin ====================================================== //   
    // Basic Variables
    Float_t momentumx, momentumy, momentumz, momentum, energy, theta;
    float positionx, positiony, positionz, position;
    float fhitx, fhity, fhitz, lhitx, lhity,lhitz, fhit, lhit;
    int MCindex, Trackindex;
    float charge;
   
    // Detector Variables
   
    // == MVD ==
    Float_t MvdDEDX;
    Float_t MvdHits;
   
    // == STT ==
    Float_t SttMeanDEDX;
    Float_t SttHits;
   
    // == GEM ==
    Float_t GemHits;
   
    // == TOF ==
    Float_t TofStopTime;
    Float_t TofM2;
    Float_t TofTrackLength;
    Float_t TofQuality;
    Float_t TofIndex;
    Float_t TofBeta;
   

    // == Barrel DIRC ==
    Float_t DrcThetaC;
    Float_t DrcThetaCErr;
    Float_t DrcQuality;
    Float_t DrcNumberOfPhotons;
    Float_t DrcIndex;
   
    // == Disc DIRC ==
    Float_t DiscThetaC;
    Float_t DiscThetaCErr;
    Float_t DiscQuality;
    Float_t DiscNumberOfPhotons;
    Float_t DiscIndex;
   
    // == RICH ==
    Float_t RichThetaC;
    Float_t RichThetaCErr;
    Float_t RichQuality;
    Float_t RichNumberOfPhotons;
    Float_t RichIndex;
   
    // == EMC ==
    Float_t EmcRawEnergy;
    Float_t EmcCalEnergy;
    Float_t EmcQuality;
    Float_t EmcNumberOfCrystals;
    Float_t EmcNumberOfBumps;
    Float_t EmcModule;
    Float_t EmcIndex;
    Float_t EoverP;
   
     // EMC Cluster properties
     Float_t EmcZ20;
     Float_t EmcZ53;
     Float_t EmcLat;
     Float_t EmcE1;
     Float_t EmcE9;
     Float_t EmcE25;
   
   // == Muons ==
   Float_t MuoProbability;
   Float_t MuoQuality;
   Float_t MuoIron;
   Float_t MuoMomentumIn;
   Float_t MuoNumberOfLayers;
   Float_t MuoModule;
   Float_t MuoHits;
   Float_t MuoIndex;
   
   // == Tracking ==
   Float_t DegreesOfFreedom;
   Float_t FitStatus;
   Float_t ChiSquared;
   Float_t PidProb;
 
//  ====================================================== Variables end   ====================================================== //

//  ====================================================== Branches Begin   ====================================================== //
TFile f(PreFix+"_tree.root","recreate");
TTree t1("t1","matched tree");


t1.Branch("momentumx",&momentumx,"momentumx/F");
t1.Branch("momentumy",&momentumy,"momentumy/F");
t1.Branch("momentumz",&momentumz,"momentumz/F");
t1.Branch("momentum",&momentum,"momentum/F");
t1.Branch("energy",&energy,"energy/F");
t1.Branch("theta",&theta,"theta/F");
t1.Branch("positionx",&positionx,"positionx/F");
t1.Branch("positiony",&positiony,"positiony/F");
t1.Branch("positionz",&positionz,"positionz/F");
t1.Branch("position",&position,"position/F");
t1.Branch("charge",&charge,"charge/F");
//t1.Branch("MCindex",&MCindex,"MCindex/I");
//t1.Branch("Trackindex",&Trackindex,"Trackindex/I");

t1.Branch("MvdDEDX",&MvdDEDX,"MvdDEDX/F");
t1.Branch("MvdHits",&MvdHits,"MvdHits/F");
t1.Branch("SttMeanDEDX",&SttMeanDEDX,"SttMeanDEDX/F");
t1.Branch("SttHits",&SttHits,"SttHits/F");
t1.Branch("GemHits",&GemHits,"GemHits/F");
t1.Branch("TofStopTime",&TofStopTime,"TofStopTime/F");
t1.Branch("TofM2",&TofM2,"TofM2/F");
t1.Branch("TofTrackLength",&TofTrackLength,"TofTrackLength/F");
t1.Branch("TofQuality",&TofQuality,"TofQuality/F");
//t1.Branch("TofIndex",&TofIndex,"TofIndex/I");
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

t1.Branch("EoverP",&EoverP,"EoverP/F");
t1.Branch("EmcRawEnergy",&EmcRawEnergy,"EmcRawEnergy/F");
t1.Branch("EmcCalEnergy",&EmcCalEnergy,"EmcCalEnergy/F");
t1.Branch("EmcQuality",&EmcQuality,"EmcQuality/F");
t1.Branch("EmcNumberOfCrystals",&EmcNumberOfCrystals,"EmcNumberOfCrystals/F");
t1.Branch("EmcNumberOfBumps",&EmcNumberOfBumps,"EmcNumberOfBumps/F");
t1.Branch("EmcModule",&EmcModule,"EmcModule/I");
t1.Branch("EmcIndex",&EmcIndex,"EmcIndex/I");
t1.Branch("EmcZ20",&EmcZ20,"EmcZ20/F");
t1.Branch("EmcZ53",&EmcZ53,"EmcZ53/F");
t1.Branch("EmcLat",&EmcLat,"EmcLat/F");
t1.Branch("EmcE1",&EmcE1,"EmcE1/F");
t1.Branch("EmcE9",&EmcE9,"EmcE9/F");
t1.Branch("EmcE25",&EmcE25,"EmcE25/F");

t1.Branch("MuoProbability",&MuoProbability,"MuoProbability/F");
t1.Branch("MuoQuality",&MuoQuality,"MuoQuality/F");
t1.Branch("MuoIron",&MuoIron,"MuoIron/F");
t1.Branch("MuoMomentumIn",&MuoMomentumIn,"MuoMomentumIn/F");
t1.Branch("MuoNumberOfLayers",&MuoNumberOfLayers,"MuoNumberOfLayers/F");
t1.Branch("MuoModule",&MuoModule,"MuoModule/F");
t1.Branch("MuoHits",&MuoHits,"MuoHits/F");
t1.Branch("MuoIndex",&MuoIndex,"MuoIndex/I");

t1.Branch("DegreesOfFreedom",&DegreesOfFreedom,"DegreesOfFreedom/I");
t1.Branch("FitStatus",&FitStatus,"FitStatus/I");
t1.Branch("ChiSquared",&ChiSquared,"ChiSquared/F");
t1.Branch("PidProb",&PidProb,"PidProb/D");
//  ====================================================== Branches end   ====================================================== //
   
    // ***
    // the event loop
    // ***
    int I=0;
    double PI=3.14159;
    //PndAnaPidSelector kaonSel("KaonSelector");
    //kaonSel.SetSelection("KaonLoose");
    while (theAnalysis->GetEvent() && i++<nevts)
        //for (Int_t i=0; i<entries; i++)
    {
                //t->GetEntry(i);
        if ((i%100)==0) cout<<"evt " << i << endl;

        // *** Select with no PID info ('All'); type and mass are set
        //theAnalysis->FillList( eminus, "ElectronAllPlus");
        //theAnalysis->FillList( piminus, "PionAllMinus");
        //theAnalysis->FillList( muminus, "MuonAllMinus");
        //theAnalysis->FillList( kminus, "KaonAllMinus");
        //theAnalysis->FillList( pminus, "ProtonAllMinus");
                theAnalysis->FillList( all, x);
                //std::cout<<"->>>>>>>>>>>>>>>>>>>> " <<all.GetLength()<<std::endl;

                SelectTruePid(theAnalysis, all);

        for (int j=0; j<all.GetLength(); ++j)
        {
                   //PidProb = all[j]->GetPidInfo(0);
                   //cout<<"Pid " <<all[j]->GetPidInfo(0)<<endl;
           //if (theAnalysis->McTruthMatch(all[j]))
           //{

             PndPidCandidate *myCand = (PndPidCandidate*) all[j]->GetRecoCandidate();
                     //int index = myCand->GetMcIndex();
             //if (index != 0) continue;
             //cout<<"PdgCode(): "<<all[j]->PdgCode()<<endl;

                momentumx=myCand->GetMomentum().X();
                momentumy=myCand->GetMomentum().Y();
                momentumz=myCand->GetMomentum().Z();

                        momentum = TMath::Sqrt( (myCand->GetMomentum().X()*myCand->GetMomentum().X()) +
                                  (myCand->GetMomentum().Y()*myCand->GetMomentum().Y()) +
                                  (myCand->GetMomentum().Z()*myCand->GetMomentum().Z()) );

           
                theta = (((all[j]->P3().Theta())*180)/PI);
               
                energy=myCand->GetEnergy();
                positionx=myCand->GetPosition().X();   
                positiony=myCand->GetPosition().Y();
                positionz=myCand->GetPosition().Z();
               
                charge=myCand->GetCharge();
                //MCindex=myCand->GetMcIndex();
                //Trackindex=myCand->GetTrackIndex();
               
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
                //TofIndex=myCand->GetTofIndex();
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
            EoverP=(EmcCalEnergy/momentum);
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
            //MuoIndex=myCand->GetMuoIndex();

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
