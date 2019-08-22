class RhoCandList;
class RhoCandidate;
class PndAnaPidSelector;
class PndAnaPidCombiner;
class PndAnalysis;

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


void tut_ana(int nevts = 0, TString prefix = "signal")
{
 	// *** some variables
	int i=0,j=0, k=0, l=0;
	
    
	
	// *** the output file for FairRunAna
	TString OutFile="Out_dummy.root";  
					
	// *** the files coming from the simulation
	TString inPidFile  = prefix+"_pid.root";    // this file contains the PndPidCandidates and McTruth
	TString inParFile  = prefix+"_par.root";
	
	// *** PID table with selection thresholds; can be modified by the user
	TString pidParFile = TString(gSystem->Getenv("VMCWORKDIR"))+"/macro/params/all.par";	
	
	// *** initialization
	//FairLogger::GetLogger()->SetLogToFile(kFALSE);
	FairRunAna* fRun = new FairRunAna(); // fRun  FairRunAna()
	FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
    FairFileSource* source = new FairFileSource(inPidFile);
    source->AddFriend("signal_full_nn.root");
	fRun->SetSource(source);
	
	// *** setup parameter database 	
	FairParRootFileIo* parIO = new FairParRootFileIo();
	parIO->open(inParFile);
	FairParAsciiFileIo* parIOPid = new FairParAsciiFileIo();
	parIOPid->open(pidParFile.Data(),"in");
	
	rtdb->setFirstInput(parIO);
	rtdb->setSecondInput(parIOPid);
	rtdb->setOutput(parIO);  
	
	fRun->SetOutputFile(OutFile);
	fRun->Init(); 
	
	// *** create an output file for all histograms
	TFile *out = TFile::Open(prefix+"_ana_full.root","RECREATE");
	
	// *** create some histograms
	TH1F *hjpsim_all = new TH1F("hjpsim_all","J/#psi mass (all)",200,0,4.5);
	TH1F *hpsim_all  = new TH1F("hpsim_all","#psi(2S) mass (all)",200,0,5);
	
	TH1F *hjpsim_lpid = new TH1F("hjpsim_lpid","J/#psi mass (loose pid)",200,0,4.5);
	TH1F *hpsim_lpid  = new TH1F("hpsim_lpid","#psi(2S) mass (loose pid)",200,0,5);
	
	TH1F *hjpsim_tpid = new TH1F("hjpsim_tpid","J/#psi mass (tight pid)",200,0,4.5);
	TH1F *hpsim_tpid  = new TH1F("hpsim_tpid","#psi(2S) mass (tight pid)",200,0,5);
	
	TH1F *hjpsim_trpid = new TH1F("hjpsim_trpid","J/#psi mass (true pid)",200,0,4.5);
	TH1F *hpsim_trpid  = new TH1F("hpsim_trpid","#psi(2S) mass (true pid)",200,0,5);
	
	
	TH1F *hjpsim_ftm = new TH1F("hjpsim_ftm","J/#psi mass (full truth match)",200,0,4.5);
	TH1F *hpsim_ftm  = new TH1F("hpsim_ftm","#psi(2S) mass (full truth match)",200,0,5);
	
	TH1F *hjpsim_nm = new TH1F("hjpsim_nm","J/#psi mass (no truth match)",200,0,4.5);
	TH1F *hpsim_nm  = new TH1F("hpsim_nm","#psi(2S) mass (no truth match)",200,0,5);
	
	TH1F *hjpsim_diff = new TH1F("hjpsim_diff","J/#psi mass diff to truth",100,-2,2);
	TH1F *hpsim_diff  = new TH1F("hpsim_diff","#psi(2S) mass diff to truth",100,-2,2);
	
	
	TH1F *hjpsim_vf   = new TH1F("hjpsim_vf","J/#psi mass (vertex fit)",200,0,4.5);
	TH1F *hjpsim_4cf  = new TH1F("hjpsim_4cf","J/#psi mass (4C fit)",200,0,4.5);
	TH1F *hjpsim_mcf  = new TH1F("hjpsim_mcf","J/#psi mass (mass constraint fit)",200,0,4.5);
	
	TH1F *hjpsi_chi2_vf  = new TH1F("hjpsi_chi2_vf", "J/#psi: #chi^{2} vertex fit",100,0,10);
	TH1F *hpsi_chi2_4c   = new TH1F("hpsi_chi2_4c",  "#psi(2S): #chi^{2} 4C fit",100,0,250);
	TH1F *hjpsi_chi2_mf  = new TH1F("hjpsi_chi2_mf", "J/#psi: #chi^{2} mass fit",100,0,10);
	
	TH1F *hjpsi_prob_vf  = new TH1F("hjpsi_prob_vf", "J/#psi: Prob vertex fit",100,0,1);
	TH1F *hpsi_prob_4c   = new TH1F("hpsi_prob_4c",  "#psi(2S): Prob 4C fit",100,0,1);
	TH1F *hjpsi_prob_mf  = new TH1F("hjpsi_prob_mf", "J/#psi: Prob mass fit",100,0,1);
    
    TH1F *mu_prob        = new TH1F("mu_prob","Probability Distribution for Muons",100,0.5,1);
    TH1F *pi_prob        = new TH1F("pi_prob", "Probability Distribution for Pions",100,0.5,1);
    TH1F *el_prob        = new TH1F("el_prob","Probability Distribution for Electrons",100,0.5,1);
    TH1F *ka_prob        = new TH1F("ka_prob", "Probability Distribution for Kaons",100,0.5,1);
    TH1F *pr_prob        = new TH1F("pr_prob","Probability Distribution for Protons",100,0.5,1);
    
    TH1F *mu_prob_ml        = new TH1F("mu_prob_ml","Probability Distribution for Muons",100,0.5,1);
    TH1F *pi_prob_ml        = new TH1F("pi_prob_ml", "Probability Distribution for Pions",100,0.5,1);
    TH1F *el_prob_ml        = new TH1F("el_prob_ml","Probability Distribution for Electrons",100,0.5,1);
    TH1F *ka_prob_ml        = new TH1F("ka_prob_ml", "Probability Distribution for Kaons",100,0.5,1);
    TH1F *pr_prob_ml        = new TH1F("pr_prob_ml","Probability Distribution for Protons",100,0.5,1);
	
	TH2F *hvpos = new TH2F("hvpos","(x,y) projection of fitted decay vertex",100,-2,2,100,-2,2);
	
    
    //THStack *hs_muon = new THStack("hs_muon","Probability Distribution for Muons");
    //THStack *hs_pion = new THStack("hs_pion","Probability Distribution for Pions");
    //TCanvas *cs_muon = new TCanvas("cs_muon","cs_muon",10,10,600,800);
    //TCanvas *cs_pion = new TCanvas("cs_pion","cs_pion",10,10,600,800);
    
	//
	// Now the analysis stuff comes...
	//
	
	
	// *** the data reader object
	PndAnalysis* theAnalysis = new PndAnalysis();
	if (nevts==0) nevts= theAnalysis->GetEntries();
	
	// *** RhoCandLists for the analysis
	RhoCandList muplus, muminus, piplus, piminus, jpsi, psi2s, electron, kaon, proton;
	
	// *** Mass selector for the jpsi cands
	double m0_jpsi = TDatabasePDG::Instance()->GetParticle("J/psi")->Mass();   // Get nominal PDG mass of the J/psi
	RhoMassParticleSelector *jpsiMassSel=new RhoMassParticleSelector("jpsi",m0_jpsi,1.0);
	
	// *** the lorentz vector of the initial psi(2S)
	TLorentzVector ini(0, 0, 6.231552, 7.240065);
	// ***
	// the event loop
	// ***
	while (theAnalysis->GetEvent() && i++<nevts)
	{
		if ((i%100)==0) cout<<"evt " << i << endl;
				
		// *** Select with no PID info ('All'); type and mass are set 		
		theAnalysis->FillList(muplus,  "MuonAllPlus");
		theAnalysis->FillList(muminus, "MuonAllMinus");
		theAnalysis->FillList(piplus,  "PionAllPlus");
		theAnalysis->FillList(piminus, "PionAllMinus");
		
		// *** combinatorics for J/psi -> mu+ mu-
		jpsi.Combine(muplus, muminus);
		
		
		// ***
		// *** do the TRUTH MATCH for jpsi
		// ***
		jpsi.SetType(443);
				
		for (j=0;j<jpsi.GetLength();++j) 
		{
			hjpsim_all->Fill( jpsi[j]->M() );
			
			if (theAnalysis->McTruthMatch(jpsi[j]))
			{ 
				hjpsim_ftm->Fill( jpsi[j]->M() );
			 	hjpsim_diff->Fill( jpsi[j]->GetMcTruth()->M() - jpsi[j]->M() );
			}
			else 
				hjpsim_nm->Fill( jpsi[j]->M() );
		}
		
		// ***
		// *** do VERTEX FIT (J/psi)
		// ***
		for (j=0;j<jpsi.GetLength();++j) 
		{
			RhoKinVtxFitter vtxfitter(jpsi[j]);	// instantiate a vertex fitter
			vtxfitter.Fit();
			
			double chi2_vtx = vtxfitter.GetChi2();	// access chi2 of fit
			double prob_vtx = vtxfitter.GetProb();	// access probability of fit
			hjpsi_chi2_vf->Fill(chi2_vtx);
			hjpsi_prob_vf->Fill(prob_vtx);			
			
			if ( prob_vtx > 0.01 )				// when good enough, fill some histos
			{
				RhoCandidate *jfit = jpsi[j]->GetFit();	// access the fitted cand
				TVector3 jVtx=jfit->Pos();		// and the decay vertex position
				
				hjpsim_vf->Fill(jfit->M());            
				hvpos->Fill(jVtx.X(),jVtx.Y());
			}
		}
		
		// *** some rough mass selection
		jpsi.Select(jpsiMassSel);
		
		// *** combinatorics for psi(2S) -> J/psi pi+ pi-
		psi2s.Combine(jpsi, piplus, piminus);
		
		
		// ***
		// *** do the TRUTH MATCH for psi(2S)
		// ***
		psi2s.SetType(88888);

		for (j=0;j<psi2s.GetLength();++j) 
		{
			hpsim_all->Fill( psi2s[j]->M() );
			
			if (theAnalysis->McTruthMatch(psi2s[j])) 
			{
			 	hpsim_ftm->Fill( psi2s[j]->M() );
			 	hpsim_diff->Fill( psi2s[j]->GetMcTruth()->M() - psi2s[j]->M() );
			}
			else 
				hpsim_nm->Fill( psi2s[j]->M() );
		}			

		
		// ***
		// *** do 4C FIT (initial psi(2S) system)
		// ***
		for (j=0;j<psi2s.GetLength();++j) 
		{
			RhoKinFitter fitter(psi2s[j]);	// instantiate the kin fitter in psi(2S)
			fitter.Add4MomConstraint(ini);	// set 4 constraint
			fitter.Fit();		            // do fit
			
			double chi2_4c = fitter.GetChi2();	// get chi2 of fit
			double prob_4c = fitter.GetProb();	// access probability of fit
			hpsi_chi2_4c->Fill(chi2_4c);
			hpsi_prob_4c->Fill(prob_4c);			
			
			if ( prob_4c > 0.01 )			// when good enough, fill some histo
			{
				RhoCandidate *jfit = psi2s[j]->Daughter(0)->GetFit();	// get fitted J/psi
				
				hjpsim_4cf->Fill(jfit->M());
			}
		}		
		
		
		// ***
		// *** do MASS CONSTRAINT FIT (J/psi)
		// ***
		for (j=0;j<jpsi.GetLength();++j) 
		{
			RhoKinFitter mfitter(jpsi[j]);		// instantiate the RhoKinFitter in psi(2S)
			mfitter.AddMassConstraint(m0_jpsi);	// add the mass constraint
			mfitter.Fit();						// do fit
			
			double chi2_m = mfitter.GetChi2();	// get chi2 of fit
			double prob_m = mfitter.GetProb();	// access probability of fit
			hjpsi_chi2_mf->Fill(chi2_m);
			hjpsi_prob_mf->Fill(prob_m);			
			
			if ( prob_m > 0.01 )				// when good enough, fill some histo
			{
				RhoCandidate *jfit = jpsi[j]->GetFit();	// access the fitted cand
				hjpsim_mcf->Fill(jfit->M());
			}
		}		
		
		
		// ***
		// *** TRUE PID combinatorics
		// ***
		
		// *** do MC truth match for PID type
		SelectTruePid(theAnalysis, muplus);
		SelectTruePid(theAnalysis, muminus);
		SelectTruePid(theAnalysis, piplus);
		SelectTruePid(theAnalysis, piminus);
				
		// *** all combinatorics again with true PID
		jpsi.Combine(muplus, muminus);
		for (j=0;j<jpsi.GetLength();++j) hjpsim_trpid->Fill( jpsi[j]->M() );
		jpsi.Select(jpsiMassSel);
		
		psi2s.Combine(jpsi, piplus, piminus);
		for (j=0;j<psi2s.GetLength();++j) hpsim_trpid->Fill( psi2s[j]->M() );
		
		
		// ***
		// *** LOOSE PID combinatorics
		// ***
		
		// *** and again with PidAlgoMvd;PidAlgoStt;PidAlgoDrc and loose selection
		theAnalysis->FillList(muplus,  "MuonLoosePlus",  "PidAlgoMvd;PidAlgoStt;PidAlgoDrc;PidAlgoMdtHardCuts;PidAlgoDisc;PidAlgoEmcBayes;PidAlgoSciT;PidAlgoRich"); //PidAlgoMl
		theAnalysis->FillList(muminus, "MuonLooseMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc;PidAlgoMdtHardCuts;PidAlgoDisc;PidAlgoEmcBayes;PidAlgoSciT;PidAlgoRich");
		theAnalysis->FillList(piplus,  "PionLoosePlus",  "PidAlgoMvd;PidAlgoStt;PidAlgoDrc;PidAlgoMdtHardCuts;PidAlgoDisc;PidAlgoEmcBayes;PidAlgoSciT;PidAlgoRich");
		theAnalysis->FillList(piminus, "PionLooseMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc;PidAlgoMdtHardCuts;PidAlgoDisc;PidAlgoEmcBayes;PidAlgoSciT;PidAlgoRich");
		
		jpsi.Combine(muplus, muminus);
		for (j=0;j<jpsi.GetLength();++j) hjpsim_lpid->Fill( jpsi[j]->M() );
		jpsi.Select(jpsiMassSel);
		
		psi2s.Combine(jpsi, piplus, piminus);
		for (j=0;j<psi2s.GetLength();++j) hpsim_lpid->Fill( psi2s[j]->M() );
    
        
        
      
        
		
		
		// ***
		// *** TIGHT PID combinatorics
		// ***
		
		// *** and again with PidAlgoMvd;PidAlgoStt and tight selection
		theAnalysis->FillList(muplus,  "MuonTightPlus",  "PidAlgoMvd;PidAlgoStt;PidAlgoDrc;PidAlgoMdtHardCuts;PidAlgoDisc;PidAlgoEmcBayes;PidAlgoSciT;PidAlgoRich");
		theAnalysis->FillList(muminus, "MuonTightMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc;PidAlgoMdtHardCuts;PidAlgoDisc;PidAlgoEmcBayes;PidAlgoSciT;PidAlgoRich");
		theAnalysis->FillList(piplus,  "PionTightPlus",  "PidAlgoMvd;PidAlgoStt;PidAlgoDrc;PidAlgoMdtHardCuts;PidAlgoDisc;PidAlgoEmcBayes;PidAlgoSciT;PidAlgoRich");
		theAnalysis->FillList(piminus, "PionTightMinus", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc;PidAlgoMdtHardCuts;PidAlgoDisc;PidAlgoEmcBayes;PidAlgoSciT;PidAlgoRich");
		  
		jpsi.Combine(muplus, muminus);
		for (j=0;j<jpsi.GetLength();++j) hjpsim_tpid->Fill( jpsi[j]->M() );
		jpsi.Select(jpsiMassSel);
		
		psi2s.Combine(jpsi, piplus, piminus);
		for (j=0;j<psi2s.GetLength();++j) hpsim_tpid->Fill( psi2s[j]->M() );
        
        
                for (Int_t ii=0; ii<muplus.GetLength(); ii++){
            double prob_m = muplus[ii]->GetPidInfo(1);
            mu_prob->Fill(prob_m);
        }
        
        for (Int_t ii=0; ii<muminus.GetLength(); ii++){
            double prob_m = muminus[ii]->GetPidInfo(1);
            mu_prob->Fill(prob_m);
        }
        
        for (Int_t ii=0; ii<piplus.GetLength(); ii++){
            double prob_p = piplus[ii]->GetPidInfo(2);
            pi_prob->Fill(prob_p);
        }
        
        for (Int_t ii=0; ii<piminus.GetLength(); ii++){
            double prob_p = piminus[ii]->GetPidInfo(2);
            pi_prob->Fill(prob_p);
        }
        
        
        
        theAnalysis->FillList(electron,  "ElectronTight",  "PidAlgoMvd;PidAlgoStt;PidAlgoDrc;PidAlgoMdtHardCuts;PidAlgoDisc;PidAlgoEmcBayes;PidAlgoSciT;PidAlgoRich");
		theAnalysis->FillList(kaon, "KaonTight", "PidAlgoMvd;PidAlgoStt;PidAlgoDrc;PidAlgoMdtHardCuts;PidAlgoDisc;PidAlgoEmcBayes;PidAlgoSciT;PidAlgoRich");
		theAnalysis->FillList(proton,  "ProtonTight",  "PidAlgoMvd;PidAlgoStt;PidAlgoDrc;PidAlgoMdtHardCuts;PidAlgoDisc;PidAlgoEmcBayes;PidAlgoSciT;PidAlgoRich");
        
        
        
        for (Int_t ii=0; ii<electron.GetLength(); ii++){
            double prob_e = electron[ii]->GetPidInfo(0);
            el_prob->Fill(prob_e);
        }
        
        for (Int_t ii=0; ii<kaon.GetLength(); ii++){
            double prob_k = kaon[ii]->GetPidInfo(3);
            ka_prob->Fill(prob_k);
        }
        
        for (Int_t ii=0; ii<proton.GetLength(); ii++){
            double prob_p = proton[ii]->GetPidInfo(4);
            pr_prob->Fill(prob_p);
        }
        
        
        
        
        
        
        
        
        
        
        // ***
		// *** ML PID combinatorics
		// ***
		
		// *** 
		theAnalysis->FillList(muplus,  "MuonTightPlus",  "PidAlgoMl");
		theAnalysis->FillList(muminus, "MuonTightMinus", "PidAlgoMl");
		theAnalysis->FillList(piplus,  "PionTightPlus",  "PidAlgoMl");
		theAnalysis->FillList(piminus, "PionTightMinus", "PidAlgoMl");
		
		
        
        
        for (Int_t ii=0; ii<muplus.GetLength(); ii++){
            double prob_m = muplus[ii]->GetPidInfo(1);
            //std::cout << "ML probability for muonplus "<< prob_m << std::endl;
            //numMuonCandidates++;
            mu_prob_ml->Fill(prob_m);
        }
        
        for (Int_t ii=0; ii<muminus.GetLength(); ii++){
            double prob_m = muminus[ii]->GetPidInfo(1);
            //std::cout << "ML probability for muonminus "<< prob_m << std::endl;
            //numMuonCandidates++;
            mu_prob_ml->Fill(prob_m);
        }
        
        for (Int_t ii=0; ii<piplus.GetLength(); ii++){
            double prob_p = piplus[ii]->GetPidInfo(2);
            //std::cout << "ML probability for pionplus "<< prob_p << std::endl;
            pi_prob_ml->Fill(prob_p);
            //numPionCandidates++;
        }
        
        for (Int_t ii=0; ii<piminus.GetLength(); ii++){
            double prob_p = piminus[ii]->GetPidInfo(2);
            //std::cout << "ML probability for pionminus "<< prob_p << std::endl;
            pi_prob_ml->Fill(prob_p);
            //numPionCandidates++;
        }
        
        
        theAnalysis->FillList(electron,  "ElectronTight",  "PidAlgoMl");
		theAnalysis->FillList(kaon, "KaonTight", "PidAlgoMl");
		theAnalysis->FillList(proton,  "ProtonTight",  "PidAlgoMl");
        
        
        
        for (Int_t ii=0; ii<electron.GetLength(); ii++){
            double prob_e = electron[ii]->GetPidInfo(0);
            el_prob_ml->Fill(prob_e);
        }
        
        for (Int_t ii=0; ii<kaon.GetLength(); ii++){
            double prob_k = kaon[ii]->GetPidInfo(3);
            ka_prob_ml->Fill(prob_k);
        }
        
        for (Int_t ii=0; ii<proton.GetLength(); ii++){
            double prob_p = proton[ii]->GetPidInfo(4);
            pr_prob_ml->Fill(prob_p);
        }
        
        
        
        
        
		
	}
	// *** write out all the histos
	out->cd();
	

    hjpsim_all->SetLineWidth(2);
	hjpsim_all->Write();
    
    hpsim_all->SetLineWidth(2);
	hpsim_all->Write();
    
    hjpsim_lpid->SetLineWidth(2);
	hjpsim_lpid->Write();
    
    hpsim_lpid->SetLineWidth(2);
	hpsim_lpid->Write();
    
    hjpsim_tpid->SetLineWidth(2);
	hjpsim_tpid->Write();
    
    hpsim_tpid->SetLineWidth(2);
	hpsim_tpid->Write();
    
    hjpsim_trpid->SetLineWidth(2);
	hjpsim_trpid->Write();
    
    hpsim_trpid->SetLineWidth(2);
	hpsim_trpid->Write();
	
    hjpsim_ftm->SetLineWidth(2);
	hjpsim_ftm->Write();
    
    hpsim_ftm->SetLineWidth(2);
	hpsim_ftm->Write();
    
    hjpsim_nm->SetLineWidth(2);
	hjpsim_nm->Write();
    
    hpsim_nm->SetLineWidth(2);
	hpsim_nm->Write();
    
	hpsim_diff->SetLineWidth(2);
	hpsim_diff->Write();
    
    hjpsim_diff->SetLineWidth(2);
	hjpsim_diff->Write();
	
    hjpsim_vf->SetLineWidth(2);
	hjpsim_vf->Write();
    
    hjpsim_4cf->SetLineWidth(2);
	hjpsim_4cf->Write();
    
    hjpsim_mcf->SetLineWidth(2);
	hjpsim_mcf->Write();
	
    hjpsi_chi2_vf->SetLineWidth(2);
	hjpsi_chi2_vf->Write();
    
    hpsi_chi2_4c->SetLineWidth(2);
	hpsi_chi2_4c->Write();
    
    hjpsi_chi2_mf->SetLineWidth(2);
	hjpsi_chi2_mf->Write();
			
    hjpsi_prob_vf->SetLineWidth(2);
	hjpsi_prob_vf->Write();
    
    hpsi_prob_4c->SetLineWidth(2);
	hpsi_prob_4c->Write();
    
    hjpsi_prob_mf->SetLineWidth(2);
	hjpsi_prob_mf->Write();
    
    mu_prob->SetFillColor(kRed);
    mu_prob_ml->SetFillColor(kBlue);
    pi_prob->SetFillColor(kRed);
    pi_prob_ml->SetFillColor(kBlue);
    el_prob->SetFillColor(kRed);
    el_prob_ml->SetFillColor(kBlue);
    ka_prob->SetFillColor(kRed);
    ka_prob_ml->SetFillColor(kBlue);
    pr_prob->SetFillColor(kRed);
    pr_prob_ml->SetFillColor(kBlue);
    
    
    mu_prob->Write();
    pi_prob->Write();
    el_prob->Write();
    ka_prob->Write();
    pr_prob->Write();
    
    mu_prob_ml->Write();
    pi_prob_ml->Write();
    el_prob_ml->Write();
    ka_prob_ml->Write();
    pr_prob_ml->Write();
    /*hs_muon->Add(mu_prob);
    hs_muon->Add(mu_prob_ml);
    hs_muon->Write();
    hs_pion->Add(pi_prob);
    hs_pion->Add(pi_prob_ml);
    hs_pion->Write();*/
    
	hvpos->Write();
		
	out->Save();
	
}
