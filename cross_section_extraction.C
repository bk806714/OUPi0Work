#include "TMatrix.h"
#include "TH1.h"
#include "TH2.h"
#include "TEventList.h"
#include "TFile.h"
#include "TCutG.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLorentzVector.h"
#include "TNtuple.h"
#include "TLine.h"
#include <iostream>
#include <fstream>
#include "TStyle.h"
#include "TF1.h"
#include "TMath.h"
#include "TSystem.h"
//#include "binning.h"
#define me 0.000511
#define Mp 0.938272
#define mpi 0.13497
#define barn 1.0E-24
#define tonanobarn 1.0E-09
#define nvirt 1.070
#define pol 0.862
#define ntbins 5
#define nphibins 12
#define maxruns 209
//#define maxruns 234

Double_t Gamma_factor(Double_t Q2,Double_t xB,Double_t theta_e, double E0);
double luminosityFxn(int kin);
Double_t epsilon_factor(Double_t Q2,Double_t xB,Double_t theta_e, double E0);
void GetTandTmin(Float_t E0,Float_t kpx, Float_t kpy,Float_t kpz, Float_t xB,TLorentzVector pi0,Float_t *t,Float_t *tmin);
Double_t fitf(Double_t *x, Double_t *par);
bool flagFlipHelicity(float run);
void cross_section_extraction(Int_t kin)
{
  Int_t kin1=kin/10;
  Int_t kin2=kin%10;
  Double_t lum;
  ifstream runlist,iFilecharge;
  
   
  Double_t Corr,mm2cuthigh,mm2cutlow,enecut,mlow,mhigh,Ngen;
  Double_t targetoffset,rcut,rval,simscalef;
  Double_t w1,w2;
  Double_t edgecut=3.0;
  Double_t photon2ene;
  Double_t Minv=0.134977;
  int tot_sim;
  Double_t E0,Qsq,theta,eprime,trkeff;
  Double_t Clus3Corr;
  //  const Int_t maxruns=334;
  if(kin==601){
    
    // maxruns=234;
    Ngen=1000e6; // number of events generated in simulation
    w1 = 2.97;// w1 value for pre-shower
    w2 = 4.08; // w2 value for shower
    Corr=-13.7;// correlation factor for removing correlation between invariant and missing mass
    targetoffset=0.004; // target offset in cm
    rcut=0.004; // rcut
    mm2cuthigh=1.08;// missing mass higher end cut (after the correction)
    mm2cutlow=0.3;// missing mass lower end cut 
    mlow=0.11; // invariant mass lowe end cut
    mhigh=0.16;// invariant mass higher end cut
    tot_sim=1000; // just for calling the name 
    E0= 8.51849;// beam energy
    eprime=3.594; // scattered electron energy
    theta= 24.564;// scattering angle of electron
    photon2ene=1.1;// threshold for photon energy cut (both photons)
    lum=2.34922e+10; // luminosity per nb
    trkeff=0.9382; // tracking efficiency
    Clus3Corr=1.0398;// Three cluster correction
    }
 else if(kin==603){

   
    Ngen=1900e6;
    w1 = 2.64;
    w2 = 3.63;
    targetoffset=0.004;
    rcut=0.004;
    Corr=-17.64;
    mm2cuthigh=1.08;
    mm2cutlow=0.5;
    // enecut=1.10;
    mlow=0.10;
    mhigh=0.160;
    tot_sim=1900;
    E0=10.591;//;//10.6171;// energy
    eprime=3.154;
    theta= 29.003;
    photon2ene=1.1;
    lum= 7.66394e+10;
    trkeff=0.9122;
    Clus3Corr=1.055;
  }
  else{
    cout << "ERROR : No Pion Rejector corrections initiated! Run is out of range!" << endl;
  }
  enecut=photon2ene;
  int ThE=enecut*10;
  int update=mm2cuthigh*100;
  
  double th_rad=(theta/2.0)*TMath::DegToRad();
  Qsq=4*E0*eprime*pow(sin(th_rad),2);
  Double_t BR=0.98823;
  //switch(kin){
    //  case 601:
  double tedges[ntbins+1]={0.000, 0.167, 0.333, 0.51, 0.75, 1.02}; // for kin 601 final 
  //case 602: 
  //tedges[]={0.0,0.1500, 0.30,0.50,0.65,1.05};//final t bin for kin 603
  //default: cout<<"Locha hai vai kinematic mai"<<endl;
  // }
  double histo_scale=(1.0*Clus3Corr)/(trkeff*1.075*0.995*BR); // BR is branching ration
  //  1.075 is radiative correction
  //0.995 is CER*S2 efficiency
 
  
					      
  Double_t ndf=9.0;
  Int_t sizeV=ntbins*1;//No binning in \phi for vertex variables....Carlos thesis, Maxime communication
  Int_t sizeE=ntbins*nphibins;//Bin in both \phi and t for experiment data
  const Int_t NLambda=4;//number of parametrizing terms in cross section
			//  - unpolarized = 3;polarized=1;
			//  -lambda=1-- sigmaT+eps*sigmaL
			//  -lambda=2--sigmaTL
			//  -lamda=3--sigmaTT
			//  -lambda=4--sigmaTL'
			 
  TMatrixD *KernelUnp[NLambda];//(sizeE,sizeV)
  for(Int_t up=0;up<NLambda;up++) KernelUnp[up]=new TMatrixD(sizeE,sizeV);//KernelUnp[0]=sigma_T+eps*sigma_L
  Int_t tarraysize=sizeV*sizeV; //5*5
  TArrayD DATA(tarraysize);//120*120//TArrayD DATA(kernelsize);
  TMatrixD DATAExp(sizeE,1),DATAExp_acc(sizeE,1),DATAExp_error(sizeE,1);
  TMatrixD DATAExp_error_F(sizeE,1),DATAExp_F(sizeE,1);
  TMatrixD expcoincidence(sizeE,1),expacc1(sizeE,1),expacc2(sizeE,1),expacc3(sizeE,1);
  TMatrixD nplus(sizeE,1),nminus(sizeE,1),nplus_acc(sizeE,1),nminus_acc(sizeE,1);

  //// Bishnu is adding from next line
  TMatrixD npluscoin(sizeE,1),nplus_acc1(sizeE,1),nplus_acc2(sizeE,1),nplus_acc3(sizeE,1),nplus_err(sizeE,1);
   TMatrixD nminuscoin(sizeE,1),nminus_acc1(sizeE,1),nminus_acc2(sizeE,1),nminus_acc3(sizeE,1),nminus_err(sizeE,1);
  /////// end of my addition
  Double_t intialize_it_n[NLambda]={1./pow(Qsq,0),1./pow(Qsq,0),1./pow(Qsq,0)};
  //initialize matrices to zero
  for(Int_t up=0;up<NLambda;up++){
    for(Int_t myE=0;myE<sizeE;myE++){
      DATAExp(myE,0)=0.0;nplus(myE,0)=0.0;nminus(myE,0)=0.0;
      DATAExp_acc(myE,0)=0.0;DATAExp_error(myE,0)=0.0;
      expcoincidence(myE,0)=0.0;expacc1(myE,0)=0.0;
      expacc2(myE,0)=0.0;expacc3(myE,0)=0.0;nplus_acc(myE,0)=0.0;nminus_acc(myE,0)=0.0;

      npluscoin(myE,0)=0.0;nplus_acc1(myE,0)=0.0;nplus_acc2(myE,0)=0.0;
      nplus_acc3(myE,0)=0.0;nplus_err(myE,0)=0.0;
      nminuscoin(myE,0)=0.0;nminus_acc1(myE,0)=0.0;nminus_acc2(myE,0)=0.0;
      nminus_acc3(myE,0)=0.0;nminus_err(myE,0)=0.0;
      for(Int_t myV=0;myV<sizeV;myV++)
	{
	  (*KernelUnp[up])(myE,myV)=0.0;// myE=12*5    //myV=5
	}
    }
  }

  
 
  
  TCut helplus;TCut helminus;
 
  
  if(kin==601){
  helplus = Form("hel==+1");//
  helminus = Form("hel==-1");//
  }
  else if (kin==603)
    {
      helplus = Form("hel==-1");//
      helminus = Form("hel==+1");//
    }
   
  TCut finalGoodExclusive=Form("mm2-%f*(m-%f)<=%f&&mm2-%f*(m-%f)>=%f&&m>=%f&&m<=%f",Corr,Minv,mm2cuthigh,Corr,Minv,mm2cutlow,mlow,mhigh);
   TCut GoodVertexsimu=Form("v-%d>=-6.5&&v-%d<=6.5",0,0);//0.75,0.75);
  TCut Goodelectronsimu=Form("rval>%f",rcut);
   
  Double_t mint=-1.0;Double_t maxt=0.0;Double_t minphi=0.0;Double_t maxphi=360.0;
  // Double_t tedges[ntbins+1]={-0.0005,0.067,0.13,0.204,0.3,0.8};// MONGI'S KIN 361
  Double_t taverage[ntbins];
  for(Int_t tav=0;tav<ntbins;tav++)taverage[tav]=(tedges[tav]+tedges[tav+1])/2.;
  
  TH2D *hexp_main = new TH2D("hexp_main","t vs. #phi",nphibins,minphi,maxphi,ntbins,tedges);
  TH2D *hexp_acc1 = new TH2D("hexp_acc1","t vs. #phi",nphibins,minphi,maxphi,ntbins,tedges);
  TH2D *hexp_acc2 = new TH2D("hexp_acc2","t vs. #phi",nphibins,minphi,maxphi,ntbins,tedges);
  TH2D *hexp_acc3 = new TH2D("hexp_acc3","t vs. #phi",nphibins,minphi,maxphi,ntbins,tedges);

   // FOR HELICITY DEPENDENT STUDIES//////
 //   1D \phi bins for \pi^0 production asymmetries for positive helicity
  
  TH1D *hphi_helplus_main;
  TH1D *hphi_helplus_acc1;
  TH1D *hphi_helplus_acc2;
  TH1D *hphi_helplus_acc3;
  
  hphi_helplus_main = new TH1D("hphi_helplus_main","#phi hel^{+}",nphibins,minphi,maxphi);
  hphi_helplus_acc1 = new TH1D("hphi_helplus_acc1","#phi hel^{+}",nphibins,minphi,maxphi);
  hphi_helplus_acc2 = new TH1D("hphi_helplus_acc2","#phi hel^{+}",nphibins,minphi,maxphi);
  hphi_helplus_acc3 = new TH1D("hphi_helplus_acc3","#phi hel^{+}",nphibins,minphi,maxphi);
  //   1D \phi bins for \pi^0 production asymmetries for negative helicity 
  TH1D *hphi_helminus_main;
  hphi_helminus_main = new TH1D("hphi_helminus_main","#phi hel^{-}",nphibins,minphi,maxphi);
  TH1D *hphi_helminus_acc1;
  hphi_helminus_acc1 = new TH1D("hphi_helminus_acc1","#phi hel^{-}",nphibins,minphi,maxphi);
  TH1D *hphi_helminus_acc2;
  hphi_helminus_acc2 = new TH1D("hphi_helminus_acc2","#phi hel^{-}",nphibins,minphi,maxphi);
  TH1D *hphi_helminus_acc3;
  hphi_helminus_acc3 = new TH1D("hphi_helpminus_acc3","#phi hel^{-}",nphibins,minphi,maxphi);

  TH1D *hplus_main_run[maxruns],*hplus_acc1_run[maxruns],*hplus_acc2_run[maxruns],*hplus_acc3_run[maxruns];
  TH1D *hminus_main_run[maxruns],*hminus_acc1_run[maxruns],*hminus_acc2_run[maxruns],*hminus_acc3_run[maxruns];

  TH2D *hsimrec = new TH2D("hsimrec","Simulation t_{rec} vs. #phi_{rec}",nphibins,minphi,maxphi,ntbins,tedges);
  TH1D *hsimver = new TH1D("hsimver","Simulation vertex t_{ver}",ntbins,tedges);
  TH1F *htrecbin = new TH1F("htrecbin","t rec. ",ntbins,tedges);
  TH1F *htverbin = new TH1F("htverbin","t ver. ",ntbins,tedges);
  TH1F *hphirecbin = new TH1F("hphirecbin","#phi rec. ",nphibins,minphi,maxphi);
  TH1F *hphiverbin = new TH1F("hphiverbin","#phi ver. ",nphibins,minphi,maxphi);
  TH1F *hepsilon = new TH1F("hepsilon","#epsilon ",100,0.,1.);hepsilon->GetXaxis()->SetTitle("#epsilon");

  Float_t mm2,m,t,tmin,phi,cer,pr1,pr2;
  TLorentzVector gamma1,gamma2,pi0;
  Float_t q1x,q1y,q1z,q2x,q2y,q2z;
  Float_t ene1,ene2,kpx,kpy,kpz;
  Float_t Q2,xB,nu;
  

   TH1F *htbincoin = new TH1F("htbincoin","t ",ntbins,tedges);
  //  TH1F *htbincoin = new TH1F("htbincoin","t ",100,-10,10);
  TH1F *hphibincoin = new TH1F("hphibincoin","#phi  ",nphibins,minphi,maxphi);
  TH1F *hexpcount[ntbins],*hexpacc1[ntbins],*hexpacc2[ntbins],*hexpacc3[ntbins]; 
  for(Int_t md=0;md<ntbins;md++)
    {
    hexpcount[md]= new TH1F(Form("hexpcount_%d",(md+1)),Form("hexpcount tbin_%d",(md+1)),nphibins,minphi,maxphi);
    hexpacc1[md]= new TH1F(Form("hexpacc1_%d",(md+1)),Form("hexpacc1 tbin_%d",(md+1)),nphibins,minphi,maxphi);
    hexpacc2[md]= new TH1F(Form("hexpacc2_%d",(md+1)),Form("hexpacc2 tbin_%d",(md+1)),nphibins,minphi,maxphi);
    hexpacc3[md]= new TH1F(Form("hexpacc3_%d",(md+1)),Form("hexpacc3 tbin_%d",(md+1)),nphibins,minphi,maxphi);
    }

    TString file_path="./";
 
     
      TChain *ntu2_33 =  new TChain("ntu");//[-3,3] & [-3,3] - ccc
      TChain *ntu2_m115_p511 =  new TChain("ntu");//[-11,-5] & [5,11]- aaa/acc3
      TChain *ntu2_33_m115 =  new TChain("ntu");//[-3,3] & [-11,-5] - cac/acc2
      TChain *ntu2_m115_m115= new TChain("ntu");//[-11,-5] & [-11,-5]
     
      ntu2_33->Add(Form("%s/pi0_-3_3_%d_New.root",file_path.Data(),kin));//[-3,3]
      ntu2_33_m115->Add(Form("%s/pi0_acc2_%d_New.root",file_path.Data(),kin));// ACCIDENTAL 2
      ntu2_m115_p511->Add(Form("%s/pi0_acc3_%d_New.root",file_path.Data(),kin));//ACCIDENTAL 3
      ntu2_m115_m115->Add(Form("%s/pi0_-11_-5_%d_New.root",file_path.Data(),kin));
      
      // Loop on data files - fill main coinc. histo
      ntu2_33->SetBranchAddress("t",&t);
      ntu2_33->SetBranchAddress("tmin",&tmin);
     
      ntu2_33->SetBranchAddress("xB",&xB);
      ntu2_33->SetBranchAddress("phi",&phi);
      
      ntu2_33->Draw(">>eventlist", finalGoodExclusive);//ccc
      ntu2_33->Draw(">>eventlisthelminus", finalGoodExclusive&&helminus);//ccc
      ntu2_33->Draw(">>eventlisthelplus",finalGoodExclusive && helplus);//ccc
      TEventList *elist = (TEventList*)gDirectory->Get("eventlist");
      
      TEventList *elisthelminus = (TEventList*)gDirectory->Get("eventlisthelminus");
      TEventList *elisthelplus = (TEventList*)gDirectory->Get("eventlisthelplus");
      // cout<<"Check this   "<<xrun<<"    "<<elist->GetN()<<endl;
      //      cout<<"   where are you   "<<elist->GetN()<<endl;
      for(Int_t i=0;i<elist->GetN();i++)
	{
	  ntu2_33->GetEntry(elist->GetEntry(i));
	  
	  hexp_main->Fill(phi*TMath::RadToDeg(),t);
	  Int_t bint=htbincoin->Fill(tmin-t);
	  Int_t binphi=hphibincoin->Fill(phi*TMath::RadToDeg());
	  //if(binphi==12&&bint==5)cout<<"LAST BIN IN PHI: "<<binphi<<" BIN in t: "<<bint<<endl;
	  if((0<bint&&bint<(ntbins+1))&&(0<binphi&&binphi<(nphibins+1)))
	    {
	      hexpcount[bint-1]->Fill(phi*TMath::RadToDeg(),histo_scale);
	      Int_t finalbin=(bint-1)*nphibins+binphi-1;//(bintrec-1)*nphibins+(binphirec-1);
	      DATAExp[finalbin][0]+=1.0*histo_scale;//Unpolarized cross section
	      expcoincidence[finalbin][0]+=1.0*histo_scale;
	    }
	}
      
      // HELICITY DEPENDENT ANALYSIS
      
      for(Int_t i=0;i<elisthelminus->GetN();i++)
	{
	  ntu2_33->GetEntry(elisthelminus->GetEntry(i));
	 
	  hphi_helminus_main->Fill(phi*TMath::RadToDeg(),histo_scale);
	  // hminus_main_run->Fill(phi*TMath::RadToDeg(),histo_scale);
	  Int_t bint=htbincoin->Fill(tmin-t);
	  Int_t binphi=hphibincoin->Fill(phi*TMath::RadToDeg());
	  if((0<bint&&bint<ntbins+1)&&(0<binphi&&binphi<nphibins+1))
	    {
	      Int_t finalbin=(bint-1)*nphibins+binphi-1;//(bintrec-1)*nphibins+(binphirec-1);
	      nminus[finalbin][0]+=1.0*histo_scale/pol;
	      nminuscoin[finalbin][0]+=1.0*histo_scale/pol;
	      //cout<<"TMin-T BIN = "<<bint<<endl;
	    }
	  
	}
      
      for(Int_t i=0;i<elisthelplus->GetN();i++)
	{
	  ntu2_33->GetEntry(elisthelplus->GetEntry(i));
	 
	  hphi_helplus_main->Fill(phi*TMath::RadToDeg(),histo_scale);
	  // hplus_main_run->Fill(phi*TMath::RadToDeg(),histo_scale);
	  Int_t bint=htbincoin->Fill(tmin-t);
	  Int_t binphi=hphibincoin->Fill(phi*TMath::RadToDeg());
	  if((0<bint&&bint<ntbins+1)&&(0<binphi&&binphi<nphibins+1))
	    {
	      Int_t finalbin=(bint-1)*nphibins+binphi-1;//(bintrec-1)*nphibins+(binphirec-1);
	      nplus[finalbin][0]+=1.0*histo_scale/pol;
	      npluscoin[finalbin][0]+=1.0*histo_scale/pol;
	      //cout<<"TMin-T BIN = "<<bint<<endl;
	    }
	}
   
      delete ntu2_33;
      delete elist;
      delete elisthelminus;
      delete elisthelplus;
      
      //fill acc1 histo
       ntu2_m115_m115->SetBranchAddress("t",&t);
       ntu2_m115_m115->SetBranchAddress("tmin",&tmin);
      
      ntu2_m115_m115->SetBranchAddress("xB",&xB);
      ntu2_m115_m115->SetBranchAddress("phi",&phi);
      ntu2_m115_m115->Draw(">>eventlistacc1",finalGoodExclusive);//accidental 1 - unpolarized
      ntu2_m115_m115->Draw(">>eventlistacc1helplus",finalGoodExclusive&&helplus);
      ntu2_m115_m115->Draw(">>eventlistacc1helminus",finalGoodExclusive&&helminus);
       
      TEventList *elistacc1 = (TEventList*)gDirectory->Get("eventlistacc1");
      TEventList *elistacc1helplus = (TEventList*)gDirectory->Get("eventlistacc1helplus");
      TEventList *elistacc1helminus = (TEventList*)gDirectory->Get("eventlistacc1helminus");
       
      for(Int_t i=0;i<elistacc1->GetN();i++)
	{
	  ntu2_m115_m115->GetEntry(elistacc1->GetEntry(i));
	 
	  Int_t bint=htbincoin->Fill(tmin-t);
	  Int_t binphi=hphibincoin->Fill(phi*TMath::RadToDeg());
	  if((0<bint&&bint<ntbins+1)&&(0<binphi&&binphi<nphibins+1))
	    {
	      hexpacc1[bint-1]->Fill(phi*TMath::RadToDeg(),histo_scale);
	      Int_t finalbin=(bint-1)*nphibins+binphi-1;//(bintrec-1)*nphibins+(binphirec-1);
	      DATAExp_acc[finalbin][0]+=1.0*histo_scale;//Unpolarized cross section - acc1
	      expacc1[finalbin][0]+=1.0*histo_scale;
	      //cout<<"TMin-T BIN = "<<bint<<endl;
	    }
	  hexp_acc1->Fill(phi*TMath::RadToDeg(),t);
	}
      
      for(Int_t i=0;i<elistacc1helplus->GetN();i++)
	{
	  ntu2_m115_m115->GetEntry(elistacc1helplus->GetEntry(i));
	  hphi_helplus_acc1->Fill(phi*TMath::RadToDeg(),histo_scale);
	  // hplus_acc1_run->Fill(phi*TMath::RadToDeg(),histo_scale);
	  Int_t bint=htbincoin->Fill(tmin-t);
	  Int_t binphi=hphibincoin->Fill(phi*TMath::RadToDeg());
	  if((0<bint&&bint<ntbins+1)&&(0<binphi&&binphi<nphibins+1))
	    {
	      Int_t finalbin=(bint-1)*nphibins+binphi-1;//(bintrec-1)*nphibins+(binphirec-1);
	      nplus_acc[finalbin][0]+=1.0*histo_scale/pol;
	      // npluscoin[finalbin][0]+=1.0*histo_scale/pol;
	      nplus_acc1[finalbin][0]+=1.0*histo_scale/pol;// bishnu
	      //cout<<"TMin-T BIN = "<<bint<<endl;
	    }
	}
      for(Int_t i=0;i<elistacc1helminus->GetN();i++)
	{
	  ntu2_m115_m115->GetEntry(elistacc1helminus->GetEntry(i));
	  hphi_helminus_acc1->Fill(phi*TMath::RadToDeg(),histo_scale);
	  // hminus_acc1_run->Fill(phi*TMath::RadToDeg(),histo_scale);
	  Int_t bint=htbincoin->Fill(tmin-t);
	  Int_t binphi=hphibincoin->Fill(phi*TMath::RadToDeg());
	  if((0<bint&&bint<ntbins+1)&&(0<binphi&&binphi<nphibins+1))
	    {
	      Int_t finalbin=(bint-1)*nphibins+binphi-1;//(bintrec-1)*nphibins+(binphirec-1);
	      nminus_acc[finalbin][0]+=1.0*histo_scale/pol;
	       nminus_acc1[finalbin][0]+=1.0*histo_scale/pol;// bishnu correction
	    }
	}
      
      delete ntu2_m115_m115;
      delete elistacc1;
      delete elistacc1helplus;
      delete elistacc1helminus;
      //fill acc2 histo 
       ntu2_33_m115->SetBranchAddress("t",&t);
       ntu2_33_m115->SetBranchAddress("tmin",&tmin);
     
       ntu2_33_m115->SetBranchAddress("xB",&xB);
       ntu2_33_m115->SetBranchAddress("phi",&phi);
       ntu2_33_m115->Draw(">>eventlistacc2",finalGoodExclusive);//unpolarized
       ntu2_33_m115->Draw(">>eventlistacc2helplus",finalGoodExclusive&&helplus);
       ntu2_33_m115->Draw(">>eventlistacc2helminus",finalGoodExclusive&&helminus);

      TEventList *elistacc2 = (TEventList*)gDirectory->Get("eventlistacc2");
      TEventList *elistacc2helplus = (TEventList*)gDirectory->Get("eventlistacc2helplus");
      TEventList *elistacc2helminus = (TEventList*)gDirectory->Get("eventlistacc2helminus");
      
      for(Int_t i=0;i<elistacc2->GetN();i++)
	{
	  ntu2_33_m115->GetEntry(elistacc2->GetEntry(i));
	 
	  Int_t bint=htbincoin->Fill(tmin-t);
	  Int_t binphi=hphibincoin->Fill(phi*TMath::RadToDeg());
	  if((0<bint&&bint<ntbins+1)&&(0<binphi&&binphi<nphibins+1))
	    {
	      hexpacc2[bint-1]->Fill(phi*TMath::RadToDeg(),histo_scale);
	      Int_t finalbin=(bint-1)*nphibins+binphi-1;//(bintrec-1)*nphibins+(binphirec-1);
	      DATAExp_acc[finalbin][0]+=1.0*histo_scale;//Unpolarized cross section - acc2
	      expacc2[finalbin][0]+=1.0*histo_scale;
	      //cout<<"TMin-T BIN = "<<bint<<endl;
	    }
	  hexp_acc2->Fill(phi*TMath::RadToDeg(),t);
	}
      
      for(Int_t i=0;i<elistacc2helplus->GetN();i++)
	{
	  ntu2_33_m115->GetEntry(elistacc2helplus->GetEntry(i));
	  hphi_helplus_acc2->Fill(phi*TMath::RadToDeg(),histo_scale);
	  //  hplus_acc2_run->Fill(phi*TMath::RadToDeg(),histo_scale);
	  Int_t bint=htbincoin->Fill(tmin-t);
	  Int_t binphi=hphibincoin->Fill(phi*TMath::RadToDeg());
	  if((0<bint&&bint<ntbins+1)&&(0<binphi&&binphi<nphibins+1))
	    {
	      Int_t finalbin=(bint-1)*nphibins+binphi-1;//(bintrec-1)*nphibins+(binphirec-1);
	      nplus_acc[finalbin][0]+=1.0*histo_scale/pol;
	       nplus_acc2[finalbin][0]+=1.0*histo_scale/pol;
	      //cout<<"TMin-T BIN = "<<bint<<endl;
	    }
	}
      for(Int_t i=0;i<elistacc2helminus->GetN();i++)
	{
	  ntu2_33_m115->GetEntry(elistacc2helminus->GetEntry(i));
	  hphi_helminus_acc2->Fill(phi*TMath::RadToDeg(),histo_scale);
	  //	  hminus_acc2_run->Fill(phi*TMath::RadToDeg(),histo_scale);
	  Int_t bint=htbincoin->Fill(tmin-t);
	  Int_t binphi=hphibincoin->Fill(phi*TMath::RadToDeg());
	  if((0<bint&&bint<ntbins+1)&&(0<binphi&&binphi<nphibins+1))
	    {
	      Int_t finalbin=(bint-1)*nphibins+binphi-1;//(bintrec-1)*nphibins+(binphirec-1);
	      nminus_acc[finalbin][0]+=1.0*histo_scale/pol;
	      nminus_acc2[finalbin][0]+=1.0*histo_scale/pol;
	      //cout<<"TMin-T BIN = "<<bint<<endl;
	    }
	}
      
      delete ntu2_33_m115;
      delete elistacc2;
      delete elistacc2helplus;
      delete elistacc2helminus;
      //fill acc3 histo 
        ntu2_m115_p511->SetBranchAddress("t",&t);
       ntu2_m115_p511->SetBranchAddress("tmin",&tmin);
      
     
      ntu2_m115_p511->SetBranchAddress("xB",&xB);
      ntu2_m115_p511->SetBranchAddress("phi",&phi);
      ntu2_m115_p511->Draw(">>eventlistacc3",finalGoodExclusive);//unpolarized
      ntu2_m115_p511->Draw(">>eventlistacc3helplus",finalGoodExclusive&&helplus);
      ntu2_m115_p511->Draw(">>eventlistacc3helminus",finalGoodExclusive&&helminus);
      
      TEventList *elistacc3 = (TEventList*)gDirectory->Get("eventlistacc3");
      TEventList *elistacc3helplus = (TEventList*)gDirectory->Get("eventlistacc3helplus");
      TEventList *elistacc3helminus = (TEventList*)gDirectory->Get("eventlistacc3helminus");
      
      for(Int_t i=0;i<elistacc3->GetN();i++)
	{
	  ntu2_m115_p511->GetEntry(elistacc3->GetEntry(i));
	  // gamma1.SetPxPyPzE(q1x,q1y,q1z,ene1);
	  // gamma2.SetPxPyPzE(q2x,q2y,q2z,ene2);
	  // pi0=gamma1+gamma2;
	  // GetTandTmin( E0, kpx,  kpy, kpz, xB, pi0, &t, &tmin);
	  Int_t bint=htbincoin->Fill(tmin-t);
	  Int_t binphi=hphibincoin->Fill(phi*TMath::RadToDeg());
	  if((0<bint&&bint<ntbins+1)&&(0<binphi&&binphi<nphibins+1))
	    {
	      hexpacc3[bint-1]->Fill(phi*TMath::RadToDeg(),histo_scale);
	      Int_t finalbin=(bint-1)*nphibins+binphi-1;//(bintrec-1)*nphibins+(binphirec-1);
	      DATAExp_acc[finalbin][0]-=1.0*histo_scale;//Unpolarized cross section - acc3
	      expacc3[finalbin][0]+=1.0*histo_scale;
	      //cout<<"TMin-T BIN = "<<bint<<endl;
	    }
	  hexp_acc3->Fill(phi*TMath::RadToDeg(),t);
	}
      
      for(Int_t i=0;i<elistacc3helplus->GetN();i++)
	{
	  ntu2_m115_p511->GetEntry(elistacc3helplus->GetEntry(i));
	  hphi_helplus_acc3->Fill(phi*TMath::RadToDeg(),histo_scale);
	  //  hplus_acc3_run->Fill(phi*TMath::RadToDeg(),histo_scale);
	  Int_t bint=htbincoin->Fill(tmin-t);
	  Int_t binphi=hphibincoin->Fill(phi*TMath::RadToDeg());
	  if((0<bint&&bint<ntbins+1)&&(0<binphi&&binphi<nphibins+1))
	    {
	      Int_t finalbin=(bint-1)*nphibins+binphi-1;//(bintrec-1)*nphibins+(binphirec-1);
	      nplus_acc[finalbin][0]-=1.0*histo_scale/pol;
	      nplus_acc3[finalbin][0]+=1.0*histo_scale/pol;
	      //cout<<"TMin-T BIN = "<<bint<<endl;
	    }
	}
      for(Int_t i=0;i<elistacc3helminus->GetN();i++)
	{
	  ntu2_m115_p511->GetEntry(elistacc3helminus->GetEntry(i));
	  hphi_helminus_acc3->Fill(phi*TMath::RadToDeg(),histo_scale);
	  // hminus_acc3_run->Fill(phi*TMath::RadToDeg(),histo_scale);
	  Int_t bint=htbincoin->Fill(tmin-t);
	  Int_t binphi=hphibincoin->Fill(phi*TMath::RadToDeg());
	  if((0<bint&&bint<ntbins+1)&&(0<binphi&&binphi<nphibins+1))
	    {
	      Int_t finalbin=(bint-1)*nphibins+binphi-1;//(bintrec-1)*nphibins+(binphirec-1);
	      nminus_acc[finalbin][0]-=1.0*histo_scale/pol;
	      nminus_acc3[finalbin][0]+=1.0*histo_scale/pol;
	      //cout<<"TMin-T BIN = "<<bint<<endl;
	    }
	}
      
      delete ntu2_m115_p511;
      delete elistacc3;
      delete elistacc3helplus;
      delete elistacc3helminus;
      // }// if -3,3] is present
       //}// for loop for the run
  //Subtract accidentals from exp data
       double tot=0.0;
   for(Int_t k=0;k<sizeE;k++)
     {
       DATAExp[k][0]-=DATAExp_acc[k][0];
       tot+=DATAExp[k][0];
       nminus[k][0]-=nminus_acc[k][0];
       nplus[k][0]-=nplus_acc[k][0];
       DATAExp_error[k][0]=expcoincidence[k][0]+expacc1[k][0]+expacc2[k][0]+expacc3[k][0];//sig2;else sig=Sqrt();
       nplus_err[k][0]=npluscoin[k][0]+nplus_acc1[k][0]+nplus_acc2[k][0]+nplus_acc3[k][0];
       nminus_err[k][0]=nminuscoin[k][0]+nminus_acc1[k][0]+nminus_acc2[k][0]+nminus_acc3[k][0];
     }
   // DATAExp.Print();
   //  DATAExp_error.Print();
   // nminus.Print();
   // nplus.Print();
   //  htbincoin->Draw();
   

///////////////////////////////////// SIMULATION /////////////////////////////
/////////////////////////////////////////////////////////////////////////////
   TChain *simu = new TChain("tr1");
   // Int_t ThE=enecut*10;    
   //  simu->Add(Form("/home/bishnu/pi0/simFiles/trialFiles/dvcs_sim_new_smeared_%d_%dM_corr_Plustheta_newDB.root",kin,1000));
   //  simu->Add(Form("/home/bishnu/pi0/simFiles/trialFiles/dvcs_sim_new_smeared_%d_%dM_corr_Plustheta_newDB_1.root",kin,1000));
   
   simu->Add(Form("%s/dvcs_smeared_simu_%d_%dM_%dThresh.root",file_path.Data(),kin,tot_sim,ThE));
  
   //TString file_path="/media/bishnu/0172A59E02A74CC1/From_PC/pi0/simFiles/trialFiles";     

   

  
   Double_t phisim_rec,phisim_ver,psf,kz_v,ky_v,kx_v,kpx_v;
   Float_t tsim_rec,tsim_min_rec,tsim_ver,tsim_min_ver;
   Double_t kpy_v,kpz_v,xB_v,Q2_v;
   /////////////////////////////////////
   Double_t kpy_m,kpz_m,kpx_m,xB_m,Q2m;
   double q1x_v,q1y_v, q1z_v,q2x_v,q2y_v, q2z_v;
   double q1x_m,q1y_m, q1z_m,q2x_m,q2y_m, q2z_m;
   double ene1_m,ene2_m;
   TLorentzVector gamma1_m,gamma2_m,pi0_m,gamma1_v,gamma2_v,pi0_v;
   kx_v=0.;ky_v=0.;
   Double_t gammasimu, eps,E0s,mm2BS,minvBS,mm3; 
   Bool_t finalGoodExclusiveB;
   simu->SetBranchAddress("phi",&phisim_rec);
   simu->SetBranchAddress("phi_v",&phisim_ver);
   simu->SetBranchAddress("psf",&psf);
   ///////////////////////////////////
   simu->SetBranchAddress("q2_v",&Q2_v);
   simu->SetBranchAddress("xB_v",&xB_v);
   simu->SetBranchAddress("xB",&xB_m);
   simu->SetBranchAddress("mm2",&mm2BS);
   simu->SetBranchAddress("m",&minvBS);
   simu->SetBranchAddress("tRec",&tsim_rec);
   simu->SetBranchAddress("tmin_Rec",&tsim_min_rec);
   simu->SetBranchAddress("t_v",&tsim_ver);
   simu->SetBranchAddress("t_min",&tsim_min_ver);
   simu->SetBranchAddress("t_v",&tsim_ver);
   simu->SetBranchAddress("gammaFac",&gammasimu);
   simu->SetBranchAddress("eE_v",&E0s);
   simu->SetBranchAddress("epsilon",&eps);
   // simu->SetBranchAddress("mm3",&mm3);



  
   Int_t nsimevent=simu->GetEntries();
   //  cout<<"sim events     "<<nsimevent<<"  lum   "<<lum<<endl;
   TH1F *hvertbins= new TH1F("hvertbins","hvertex t bins",5,0,5);
   TH1F *hthetae=new TH1F("hthetae","theta_e",20,10,30);
   TH1F *htv_tm=new TH1F("htv_tm","tv_tm",200,-0.2,0.2);
   TH1F *hsimestimate[ntbins]; for(Int_t md=0;md<ntbins;md++)hsimestimate[md]= new TH1F(Form("hsimestimate_%d",(md+1)),Form("simestimate tbin_%d",(md+1)),nphibins,minphi,maxphi);
 
   // Double_t lum=luminosityFxn(kin);
   Double_t lumhelicity_p=lum*0.985/2.; // the factor 0.985 has to do the uncertainty on timing of helicity flip 
   Double_t lumhelicity_n=lum*0.985/2.;
   // cout<<"My Luminosity = "<<lum<<" /cm2"<<endl;
   

 for(Int_t i=0;i<simu->GetEntries();i++)
    {
      simu->GetEntry(i);
      mm3=(mm2BS-Corr*(minvBS-Minv));
       finalGoodExclusiveB= (mm3<=mm2cuthigh &&mm3>=mm2cutlow && minvBS>=mlow && minvBS<=mhigh);
      if(finalGoodExclusiveB){

       Double_t fac = psf/Ngen;
       Double_t epsL=eps;
       hepsilon->Fill(eps);
       hsimrec->Fill(phisim_rec*TMath::RadToDeg(),tsim_rec);
       hsimver->Fill(tsim_ver);
       Int_t bintrec=htrecbin->Fill(tsim_min_rec-tsim_rec);
       Int_t bintver=htverbin->Fill(tsim_min_ver-tsim_ver);
       Int_t binphirec=hphirecbin->Fill(phisim_rec*TMath::RadToDeg());
       Int_t binphiver=hphiverbin->Fill(phisim_ver*TMath::RadToDeg());
      
       if((0<bintrec&&bintrec<ntbins+1)&&(0<bintver&&bintver<ntbins+1)&&(0<binphirec&&binphirec<nphibins+1))
	 {
	   Int_t recbin=(bintrec-1)*nphibins+binphirec-1;
	   Int_t verbin=bintver-1;
	   htv_tm->Fill(tsim_min_ver-tsim_min_rec);
	   hsimestimate[bintrec-1]->Fill(phisim_ver*TMath::RadToDeg(),fac*gammasimu/(2.*TMath::Pi()));
	   //cout<<"REC BIN : "<<binphirec<<"  VER BIN  : "<<verbin<<endl;
	   hvertbins->Fill(verbin);
	   ((*KernelUnp[0])(recbin,verbin))+=fac*gammasimu/(2.*TMath::Pi());//sigma_T+eps*sigma_L
	   ((*KernelUnp[1])(recbin,verbin))+=fac*gammasimu*TMath::Sqrt(2.*epsL*(1.+eps))*TMath::Cos(phisim_ver)/(2.*TMath::Pi());//sqrt(2*eps*(1+eps))*sigma_LT*cos(phi)
	   ((*KernelUnp[2])(recbin,verbin))+=fac*gammasimu*eps*TMath::Cos(2.*phisim_ver)/(2.*TMath::Pi());//eps*sigma_TT*cos(2*phi)
	   ((*KernelUnp[3])(recbin,verbin))+=fac*gammasimu*TMath::Sqrt(2.*epsL*(1.-eps))*TMath::Sin(phisim_ver)/(2.*TMath::Pi());//eps*sigma_TT*cos(2*phi)....
	 }
      
     }
    }
   delete simu;
   // delete elistsimu;
   //  htverbin->Draw();
			 
   cout<<" done simu  "<<"Done Here"<<endl;
   // Minimization here
   TMatrixD *alpha=new TMatrixD((NLambda-1)*sizeV,(NLambda-1)*sizeV);
   TMatrixD *X    =new TMatrixD((NLambda-1)*sizeV,1);
   TMatrixD *beta =new TMatrixD((NLambda-1)*sizeV,1);
   // TMatrixD *Err= new TMatrixD((NLambda-1)*sizeV,1);// =new TMatrixD(NLambda*sizeV,1);
     TMatrixD Err((NLambda-1)*sizeV,1);// =new TMatrixD(NLambda*sizeV,1);
  
   TMatrixD *alphaP=new TMatrixD(sizeV,sizeV);
   TMatrixD *XP    =new TMatrixD(sizeV,1);
   TMatrixD *betaP =new TMatrixD(sizeV,1);
   // TMatrixD *ErrP= new TMatrixD(sizeV,1);// =new TMatrixD(NLambda*sizeV,1);
     TMatrixD ErrP(sizeV,1);// =new TMatrixD(NLambda*sizeV,1);
 
 for(Int_t lambda=0;lambda<NLambda-1;lambda++)
    {
    for(Int_t jv=0;jv<sizeV;jv++){
      (*beta)(jv*(NLambda-1)+lambda,0)=0.;
      for(Int_t ie=0;ie<sizeE;ie++)
	{
	  // mm2high1>>DATAExp[ie][0];
	  // DATAExp_error[ie][0]=DATAExp[ie][0];
	  if(DATAExp_error[ie][0]==0)DATAExp_error[ie][0]=1;
	  (*beta)(jv*(NLambda-1)+lambda,0)+=lum*DATAExp[ie][0]*((*KernelUnp[lambda])(ie,jv))/(DATAExp_error[ie][0]);
	}
      //  cout<<"(jv*(NLambda-1)+lambda,0)    "<<jv*(NLambda-1)+lambda<<endl;
    
         }
}

   
  // cout<<"// There you go sir////////////////////"<<endl;
  // beta->Print();
for(Int_t lambda=0;lambda<NLambda-1;lambda++){
  for(Int_t jv=0;jv<sizeV;jv++){
      for(Int_t lambdap=0;lambdap<NLambda-1;lambdap++){
       for(Int_t jvp=0;jvp<sizeV;jvp++){
	(*alpha)(jv*(NLambda-1)+lambda,jvp*(NLambda-1)+lambdap)=0.;
	for(Int_t ie=0;ie<sizeE;ie++){
	  // cout<<"THIS IS WHAT I READ "<<DATAExp[ie][0]<<endl;
	   if(DATAExp_error[ie][0]==0)DATAExp_error[ie][0]=1;
	   (*alpha)(jv*(NLambda-1)+lambda,jvp*(NLambda-1)+lambdap)+=lum*lum*((*KernelUnp[lambda])(ie,jv))*((*KernelUnp[lambdap])(ie,jvp))/(DATAExp_error[ie][0]); 
	}
	//	cout<<"alpha term "<<jv*(NLambda-1)+lambda<<" , "<<jvp*(NLambda-1)+lambdap<<endl;
      }
     }
  }
}
// alpha->Print();
     
/////////////////////////// For POLARIZED COMPONENT ///////////////////////////////////////////////
    for(Int_t jv=0;jv<sizeV;jv++){
      //  (*betaP)(jv*1,0)=0.;
      for(Int_t ie=0;ie<sizeE;ie++)
	{
	   DATAExp_error_F[ie][0]=0;
	   DATAExp_F[ie][0]=0;
	   DATAExp_F[ie][0]=0.5*(nplus[ie][0]-nminus[ie][0]);
	   DATAExp_error_F[ie][0]=0.5*( nplus_err[ie][0]+ nminus_err[ie][0]);
	    
	    if(DATAExp_error_F[ie][0]==0)DATAExp_error_F[ie][0]=1;
	  	(*betaP)(jv*1,0)+=lum*DATAExp_F[ie][0]*((*KernelUnp[3])(ie,jv))/(DATAExp_error_F[ie][0]);
	}
    }
     //  DATAExp_F.Print();
          
  for(Int_t jv=0;jv<sizeV;jv++){
         for(Int_t jvp=0;jvp<sizeV;jvp++){
	(*alphaP)(jv*1,jvp*1)=0.;
	for(Int_t ie=0;ie<sizeE;ie++){
	  //  DATAExp_error_F[ie][0]=0;
	   // DATAExp_error_F[ie][0]=0.5*( nplus_err[ie][0]+ nminus_err[ie][0]);
	    	   
	  if(DATAExp_error_F[ie][0]==0)DATAExp_error_F[ie][0]=1;
	  (*alphaP)(jv*1,jvp*1)+=lum*lum*((*KernelUnp[3])(ie,jv))*((*KernelUnp[3])(ie,jvp))/(DATAExp_error_F[ie][0]);///////
	 	   
	}
	
    }
  }
  // alphaP->Print();
  //  betaP->Print();
  
//DATAExp.Print();
 alpha->Invert();
 // alpha->Print();			       
 alphaP->Invert();
 
 for(Int_t Ic=0;Ic<(NLambda-1)*sizeV;Ic++)
   {
     Err[Ic][0]=0.0; 
     Err[Ic][0]=(*alpha)(Ic,Ic);
     Err[Ic][0]=sqrt((Err)[Ic][0]);
   }
 for(Int_t Icc=0;Icc<sizeV;Icc++)
   {
     ErrP[Icc][0]=0.0;
     ErrP[Icc][0]=(*alphaP)(Icc,Icc);
     ErrP[Icc][0]=sqrt((ErrP)[Icc][0]);
   }
 
 //Err.Print();

 Double_t sigma_params[NLambda][ntbins];
 Double_t params[NLambda][ntbins];
 //Display results
 X->Mult((*alpha),(*beta));
 XP->Mult((*alphaP),(*betaP));

 //X is 20x1 matrix (0,0),(1,0),2,0) so on (0,0) to 3,0) on t bin
 // (0,0) , (4,0), (8,0),(12,0), (16,0) 5 different t bin with lambda=0
 for(Int_t lambda=0;lambda<(NLambda-1);lambda++){
   for(Int_t jv=0;jv<sizeV;jv++){
     params[lambda][jv]=(*X)(jv*(NLambda-1)+lambda,0)/(1.0);///(2.*TMath::Pi//barn*tonanobarn
     sigma_params[lambda][jv]=(Err)(jv*(NLambda-1)+lambda,0)/(1.0);///(2.*TMath::Pi());
      if(lambda==0) cout<< params[lambda][jv]<<"   "<<   sigma_params[lambda][jv]<<endl;
   }
 }
 TFile *fi;
 fi=new TFile(Form("./Test%d.root",kin),"recreate");
 TTree *tr1=new TTree("tr1","tr1");
 // tr1->Branch("Params",*params);
 // tr1->Branch("SigmaParams",*sigma_params);
 // tr1->Branch("experiment",*experimentcount);
 // tr1->Branch("Nsim_pertbin",*Nsim_pertbin);
 for(Int_t lambda=3;lambda<NLambda;lambda++){
   for(Int_t jv=0;jv<sizeV;jv++){
     //cout<<Form("Coef %d for bin %d",lambda,jv)<<endl;
     //cout<<(*X)(jv*NLambda+lambda,0)<<endl;
     params[lambda][jv]=(*XP)(jv*1,0)/(1.0);///(2.*TMath::Pi());//barn*tonanobarn
     //   cout<<"  "<<lambda<<"   "<<jv<<"    "<<params[lambda][jv]<<endl;
     sigma_params[lambda][jv]=(ErrP)(jv*1,0)/(1.0);///(2.*TMath::Pi());barn*tonanobarn
     // tr1->Fill();
   }
 }
 

 //Calculate the number of simulation events per t[\phi] bin
 TMatrixD *Nsim[NLambda];
 TMatrixD *Nsimtot[NLambda];// =  new TMatrixD(sizeE,1);
 for(Int_t par=0;par<NLambda;par++)
   {
     Nsim[par] = new TMatrixD(sizeV,1);
     Nsimtot[par] = new TMatrixD(sizeE,1);
   }
 for(Int_t lambda=0;lambda<NLambda;lambda++)
   {
     for(Int_t par=0;par<sizeV;par++)
       {
	 (*Nsim[lambda])(par,0)=params[lambda][par];
       }
   }

 for(Int_t lambda=0;lambda<NLambda;lambda++)
   {
     Nsimtot[lambda]->Mult((*KernelUnp[lambda]),(*Nsim[lambda]));
   }
 Double_t Nsimfinal[ntbins*nphibins];//sum up Nsim for all t-bins+phibins
 for(Int_t i=0;i<sizeE;i++)
   {
     Nsimfinal[i] = ((*Nsimtot[0])(i,0)+(*Nsimtot[1])(i,0)+(*Nsimtot[2])(i,0))*1.0*lum;//barn*tonanobarn
   }
 
 Double_t Nsim_pertbin[ntbins][nphibins];
 Double_t Nsim_data_residual_pertbin[ntbins][nphibins];
 // totsig is total sigma for each of tbins
 // totsigerror is toal error in totsig
 // params[] are each parameter in xsec
 //sigma_params is error in each xsec terms
 if(kin==601)
 eps =0.6658;
 else if (kin==603)
   eps=0.50;
 else cout<<" epsilon value is not defined for this kinematic"<<endl;
 Double_t epsL=eps;
 Double_t totsig[ntbins][nphibins],plotphi[nphibins],totsigerror[ntbins][nphibins];
 Double_t totsigpoints[ntbins][nphibins];
 for(Int_t tb=0;tb<ntbins;tb++)
   {
     for(Int_t cr=0;cr<(nphibins);cr++)
       {
	 Int_t phi2 = hphirecbin->GetXaxis()->GetBinCenter(cr+1);
	 plotphi[cr] = phi2;
	 //cout<<"phi in total cross sec = "<<phi2<<endl;
	 totsig[tb][cr] = params[0][tb]+params[1][tb]*TMath::Cos(phi2*TMath::DegToRad())*TMath::Sqrt(2.*epsL*(1.+eps))+params[2][tb]*TMath::Cos(2.*phi2*TMath::DegToRad())*eps;//To be plotted as a function (line) calculated using the xsec params
	 totsigerror[tb][cr] = TMath::Sqrt(pow(sigma_params[0][tb],2)+pow(sigma_params[1][tb]*TMath::Cos(phi2*TMath::DegToRad())*TMath::Sqrt(2.*epsL*(1.+eps)),2)+pow(sigma_params[2][tb]*TMath::Cos(2.*phi2*TMath::DegToRad())*eps,2));
       }
   }

 Double_t simulation_estimate[ntbins][nphibins],experimentcount[ntbins][nphibins],experror[ntbins][nphibins],totsigmapoints_error[ntbins][nphibins];

  


 
 for(Int_t jv=0;jv<sizeV;jv++)
   {
     for(Int_t phib=0;phib<nphibins;phib++)
       {
	 simulation_estimate[jv][phib]=totsig[jv][phib]*(hsimestimate[jv]->GetBinContent(phib+1))*1.0*lum;//barn*tonanobarn
	 experimentcount[jv][phib]=hexpcount[jv]->GetBinContent(phib+1)-hexpacc1[jv]->GetBinContent(phib+1)-hexpacc2[jv]->GetBinContent(phib+1)+hexpacc3[jv]->GetBinContent(phib+1);
	 experror[jv][phib]=TMath::Sqrt(hexpcount[jv]->GetBinContent(phib+1)+hexpacc1[jv]->GetBinContent(phib+1)+hexpacc2[jv]->GetBinContent(phib+1)+hexpacc3[jv]->GetBinContent(phib+1));
	
       }
   }

 for(Int_t jv=0;jv<sizeV;jv++)
   {
     for(Int_t phib=0;phib<nphibins;phib++)
       {
	 Int_t choosetbin;
	 choosetbin=jv*nphibins+phib;// for each t bin 12 phi bins
	 Nsim_pertbin[jv][phib]= Nsimfinal[choosetbin];
	 Nsim_data_residual_pertbin[jv][phib]= experimentcount[jv][phib] - Nsim_pertbin[jv][phib];
	 totsigpoints[jv][phib] = (experimentcount[jv][phib]/Nsim_pertbin[jv][phib])*totsig[jv][phib];
	 totsigmapoints_error[jv][phib] = (experror[jv][phib]/Nsim_pertbin[jv][phib])*totsig[jv][phib];
	 tr1->Fill();
       }
   }
 Double_t chi2[ntbins]={0};
 for(Int_t jv=0;jv<ntbins;jv++)
   {
     for(Int_t phib=0;phib<nphibins;phib++)
       {
	 if(experror[jv][phib]==0)experror[jv][phib]=1;
	 chi2[jv]+=pow((Nsim_pertbin[jv][phib]-experimentcount[jv][phib])/(experror[jv][phib]),2);
       }
   }






 ///////////// WRITING THE CROSS-SECTION TERMS TO OUTPUT FILE
 /*
 // ofstream outxsecparams;
 //  outxsecparams.open(Form("/media/bishnu/0172A59E02A74CC1/From_PC/pi0/output%d/xsecparams_ETh%d.txt",kin,ThE)); //output file
 for(Int_t ip=0;ip<ntbins;ip++)
   {
    cout<<taverage[ip]<<"\t"<<params[0][ip]<<"\t"<<sigma_params[0][ip]<<"\t"<<params[1][ip]<<"\t"<<sigma_params[1][ip]<<"\t"<<params[2][ip]<<"\t"<<sigma_params[2][ip]<<"\t"<<params[3][ip]<<"\t"<<sigma_params[3][ip]<<"\t"<<endl;
   }
 */
 TMultiGraph *mg[ntbins];
 for(Int_t j=0;j<ntbins;j++)
   {
     mg[j]=new TMultiGraph();
     mg[j]->SetTitle(Form(" t_{min}-t = %.4f GeV^{2},#chi^{2}/ndf = %.2f;#phi (deg.);Number of events",taverage[j],chi2[j]/ndf));
     mg[j]->SetMinimum(0);
     mg[j]->SetMaximum(2200);
   }
 TGraphErrors *grsim[ntbins],*grexp[ntbins],*grexpsim_residue[ntbins];
 TCanvas *csimu_estimate[ntbins],*cexpsim_residue[ntbins];
 for(Int_t ip=0;ip<ntbins;ip++)
   {
     grsim[ip]= new TGraphErrors(nphibins,plotphi,Nsim_pertbin[ip],0,0);
     grexp[ip]= new TGraphErrors(nphibins,plotphi,experimentcount[ip],0,experror[ip]);
     grexpsim_residue[ip]= new TGraphErrors(nphibins,plotphi,Nsim_data_residual_pertbin[ip],0,experror[ip]);
     grexp[ip]->SetMarkerStyle(20);
     grexp[ip]->SetMarkerColor(kBlack);
     grexp[ip]->SetMarkerSize(1.1);
     csimu_estimate[ip] = new TCanvas();  
     csimu_estimate[ip]->cd();
     grsim[ip]->SetFillColor(2);
     mg[ip]->Add(grsim[ip],"AB");// draw option "AB"
     mg[ip]->Add(grexp[ip],"AP");
     mg[ip]->Draw("A");
     mg[ip]->GetYaxis()->SetTitleOffset(1.1);
     mg[ip]->GetYaxis()->SetLabelSize(0.032);
     mg[ip]->GetYaxis()->SetLabelFont(32);
     mg[ip]->GetXaxis()->SetTitleOffset(1.1);
     mg[ip]->GetXaxis()->SetLabelSize(0.034);
     mg[ip]->GetXaxis()->SetLabelFont(32);
     mg[ip]->GetXaxis()->SetTitleFont(32);
     mg[ip]->GetYaxis()->SetTitleFont(32);
     csimu_estimate[ip]->Update();
     // csimu_estimate[ip]->Print(Form("./outputfiles%d/compare_counts_t_av_%d.pdf",kin,ip));
     // csimu_estimate[ip]->Print(Form("./outputfiles%d/compare_counts_t_av_%d.png",kin,ip));
     cexpsim_residue[ip] = new TCanvas();
     cexpsim_residue[ip]->cd();
     grexpsim_residue[ip]->SetMarkerStyle(33);
     grexpsim_residue[ip]->SetMarkerSize(1.5);
     grexpsim_residue[ip]->SetTitle(Form("E_{B} = %f GeV, Q^{2} = %f GeV^{2}, x_{B}=0.60, t_{min}-t = %f GeV^{2};#phi (deg.);N^{exp.} - N^{sim.};",E0,Qsq,taverage[ip]));
     grexpsim_residue[ip]->Draw("AP");
     //  cexpsim_residue[ip]->Print(Form("./outputfiles%d/compare_counts_residue_t_av_%d.pdf",kin,ip));
     //  cexpsim_residue[ip]->Print(Form("./outputfiles%d/compare_counts_residue_t_av_%d.png",kin,ip));
   }

 TGraphErrors *grparameters[NLambda],*grtotxsec[ntbins],*grtotxsecpoints[ntbins];
 TString xsecparams[NLambda]={"#frac{d#sigma_{T}+#epsilon#sigma_{L}}{dt}","#frac{d#sigma_{TL}}{dt}","#frac{d#sigma_{TT}}{dt}","#frac{d#sigma_{TL'}}{dt}"};
 TString xsecparamsnamefile[NLambda]={"sigma_T_epsilon_sigma_L","sigma_TL","sigma_TT","sigma_TLp"};
 TString xsecparams2[NLambda]={"sigmaT+#epsilon#sigma_{L}","#sigma_{TL}","#sigma_{TT}","#sigma_{TL'}"};
 TCanvas *cparameters[NLambda],*ctotxsec[ntbins];
 gStyle->SetOptFit(1111);
 //TVirtualFitter::SetMaxIterations(300000);
 TF1 *fitfunc = new TF1("fitfunc",fitf,15.*TMath::DegToRad(),345.*TMath::DegToRad(),3);// 15 t0 345
 // fitfunc->SetParName(0,"#sigma_{T}+#epsilon#sigma_{L}");
 // fitfunc->SetParName(1,"#sigma_{TL}");
 // fitfunc->SetParName(2,"#sigma_{TT}");
 // fitfunc->SetParameter(0,40,0.25,0,60);
 //fitfunc->SetParameter(1,"#sigma_{TL}",0,0.2,-10,30);
 // fitfunc->SetParameter(2,"#sigma_{TT}",-5,0.2,-10,-20);
 // fitfunc->SetParameters(0.,40.,0.,-5.);

 
 TMultiGraph *mgtotxsec[ntbins];
 for(Int_t j=0;j<ntbins;j++)
   {
     mgtotxsec[j]=new TMultiGraph();
     mgtotxsec[j]->SetTitle(Form("E_{B} = %.2f GeV, Q^{2} = %.2f GeV^{2}, t_{min} - t = %f GeV^{2};#phi (deg.);#frac{d#sigma}{dt} (nb/GeV^{2})",E0,Qsq,taverage[j]));
     mgtotxsec[j]->SetMinimum(0);
   }

 for(Int_t ip=0;ip<NLambda;ip++)
   {
     cparameters[ip] = new TCanvas();
     cparameters[ip]->cd();
     grparameters[ip] = new TGraphErrors(ntbins,taverage,params[ip],0,sigma_params[ip]);
     grparameters[ip]->SetMarkerStyle(33);
     grparameters[ip]->SetMarkerColor(kRed);
     grparameters[ip]->SetMarkerSize(1.5);
     grparameters[ip]->SetTitle(Form("E_{B} = %.2f GeV, Q^{2} =%.2f GeV^{2};t_{min}-t (GeV^{2});%s(nbarn/GeV^{2});",E0,Qsq,xsecparams[ip].Data()));
     grparameters[ip]->Draw("AP"); 
     grparameters[ip]->GetYaxis()->SetTitleOffset(1.1);
     grparameters[ip]->GetYaxis()->SetLabelSize(0.032);
     grparameters[ip]->GetYaxis()->SetLabelFont(32);
     grparameters[ip]->GetXaxis()->SetTitleOffset(1.0);
     grparameters[ip]->GetXaxis()->SetLabelSize(0.035);
     grparameters[ip]->GetXaxis()->SetLabelFont(32);
     grparameters[ip]->GetXaxis()->SetTitleFont(32);
     grparameters[ip]->GetXaxis()->SetTitleSize(0.045);
     grparameters[ip]->GetYaxis()->SetTitleFont(32);
     
     cparameters[ip]->Update();
     //   cparameters[ip]->Print(Form("./kin361output/%s_t_av.pdf",xsecparamsnamefile[ip].Data()));
     //  cparameters[ip]->Print(Form("./kin361output/%s_t_av.png",xsecparamsnamefile[ip].Data()));
   }
 
 for(Int_t ip=0;ip<ntbins;ip++)
   {
     ctotxsec[ip] = new TCanvas();
     ctotxsec[ip]->cd();
     grtotxsec[ip] = new TGraphErrors(nphibins,plotphi,totsig[ip],0,0);
     grtotxsecpoints[ip] = new TGraphErrors(nphibins,plotphi,totsigpoints[ip],0,totsigmapoints_error[ip]);
     grtotxsecpoints[ip]->SetMarkerStyle(33);
     grtotxsecpoints[ip]->SetMarkerColor(kBlue);
     grtotxsecpoints[ip]->SetMarkerSize(1.5);
     grtotxsec[ip]->SetLineColor(2);
     grtotxsec[ip]->SetLineWidth(4);
     mgtotxsec[ip]->Add(grtotxsecpoints[ip],"AP");
     mgtotxsec[ip]->Add(grtotxsec[ip],"AC");
     mgtotxsec[ip]->Draw("A"); 
     mgtotxsec[ip]->GetYaxis()->SetTitleOffset(1.1);
     mgtotxsec[ip]->GetYaxis()->SetLabelSize(0.032);
     mgtotxsec[ip]->GetYaxis()->SetLabelFont(32);
     mgtotxsec[ip]->GetYaxis()->SetTitleFont(32);
     mgtotxsec[ip]->GetXaxis()->SetTitleFont(32);
     mgtotxsec[ip]->GetXaxis()->SetTitleOffset(1.1);
     mgtotxsec[ip]->GetXaxis()->SetLabelSize(0.034);
     mgtotxsec[ip]->GetXaxis()->SetLabelFont(32);
     mgtotxsec[ip]->GetYaxis()->SetRangeUser(0,120);
     ctotxsec[ip]->Update();
     //  ctotxsec[ip]->Print(Form("./outputfiles%d/xsectot_t_av_%d.pdf",kin,ip));
     // ctotxsec[ip]->Print(Form("./outputfiles%d/xsectot_t_av_%d.png",kin,ip));
        
     
   }
 for(Int_t lambda=0;lambda<NLambda;lambda++)
   {
     Nsimtot[lambda]->Write(Form("Nsimtot%d",lambda));
   }
 // Nsimtot->Write("Nsimtot");
 X->Write("X");
 XP->Write("XP");  
 alpha->Write("alpha");
 alphaP->Write("alphaP");
 for(Int_t hcnt=0;hcnt<ntbins;hcnt++)
   {
 hexpcount[hcnt]->Write(Form("hexpcount%d",hcnt));
 hexpacc1[hcnt]->Write(Form("hexpacc1%d",hcnt));
  hexpacc2[hcnt]->Write(Form("hexpacc2%d",hcnt));
   hexpacc3[hcnt]->Write(Form("hexpacc3%d",hcnt));
   }
 // //params->Write("params");
 // tr1->Write();
 fi->Close();
  
}
   

double luminosityFxn(int kin)
{
  //  const int maxruns;
  ifstream iFilecharge;
  //if(kin==601)maxruns=234;
  //if(kin==603)maxruns=334;
  iFilecharge.open(Form("/home/bishnu/pi0/chargeFiles/charge_LT_%d.dat",kin));
  // string titlelinecharge;
  //  std::getline(iFilecharge,titlelinecharge);
  //  cout<<titlelinecharge<<endl;

  Double_t runcharge[maxruns],runnumber[maxruns],trkeff[maxruns],livetime[maxruns],luminosity[maxruns];
  Double_t tot_luminosityF=0.;
  Double_t histo_scale[maxruns];
   Double_t chargeQ=0.0;
   for(Int_t i = 0;i<maxruns;i++)
    {
      iFilecharge>>runnumber[i]>>runcharge[i]>>livetime[i]>>luminosity[i]>>trkeff[i];
      //  chargeQ+=runcharge[i]/(livetime[i]);
      tot_luminosityF+=luminosity[i]/livetime[i];
      // histo_scale[i]=1.0/(trkeff[i]*livetime[i]);

    }
   return tot_luminosityF;
}

Double_t Gamma_factor(Double_t Q2,Double_t xB,Double_t theta_e,double E0)
{
  Double_t nu = Q2/(2*Mp*xB);
  Double_t q = sqrt(pow(nu,2)+Q2);
  Double_t s = pow(Mp,2)+Q2*((1.-xB)/xB);
  Double_t kgamma = (s-pow(Mp,2))/(2*Mp);
  Double_t alphaQED = 0.0072982;
  Double_t eps=epsilon_factor(Q2,xB,theta_e,E0);
  Double_t Gamma = (alphaQED*Q2*(1.-xB))/(8.*TMath::Pi()*pow(Mp*E0,2)*pow(xB,3)*(1.-eps));//From Maxime, makes dimensional sense
  return Gamma;
}

Double_t epsilon_factor(Double_t Q2,Double_t xB,Double_t theta_e, double E0)
{
  Double_t nu = Q2/(2*Mp*xB);
  Double_t q = sqrt(pow(nu,2)+Q2);
  Double_t y=nu/E0;
  //Double_t eps = 1/(1+2*pow(q*tan(theta_e/2),2)/Q2);
  Double_t eps = (1.-y-(Q2/(4*pow(E0,2))))/(1.-y+pow(y,2)/2+(Q2/(4*pow(E0,2))));//From Maxime
  return eps;
}

Double_t fitf(Double_t *x, Double_t *par)
{
  Double_t fitval =(par[0]+par[1]*cos(x[0]*TMath::DegToRad())+par[2]*cos(2.*x[0]*TMath::DegToRad()));///(1.0+par[1]*cos(x[0]*TMath::DegToRad()));//+par[2]*cos(2.*x[0]*TMath::DegToRad()));
  return fitval;
}
void GetTandTmin(Float_t E0,Float_t kpx, Float_t kpy,Float_t kpz, Float_t xB,TLorentzVector pi0,Float_t *t,Float_t *tmin)
{
  Float_t eprime,nu;
  eprime=sqrt(kpx*kpx+ kpy*kpy+kpz*kpz);
  nu=E0-eprime;
  TLorentzVector ki,kf,q2;
  Float_t WB,nuCM,qCM,QpCM;
  float q_exp,s,W,F2;
  float  q2Mag,ph,Eh,cosqqp;
  float Q2;
  float FacA,FacB;
  ki.SetPxPyPzE(0.,0.,E0,E0);
  kf.SetPxPyPzE(kpx, kpy, kpz,eprime);
  q2=ki-kf;
  Q2=-q2.Mag2();
  cosqqp=(q2.Vect().Dot((pi0).Vect()))/(q2.Vect().Mag()*((pi0).Vect()).Mag());
  s = Mp*Mp-q2.Mag2()*((1.-xB)/xB);
  F2=(TMath::Power(s-mpi*mpi-Mp*Mp,2)/4)-mpi*mpi*Mp*Mp;
  ph=(s+mpi*mpi-Mp*Mp)*q2.P()*cosqqp+2*(nu+Mp)*TMath::Sqrt(F2-mpi*mpi*q2.P()*q2.P()*(1-cosqqp*cosqqp));
  ph=0.5*ph/(s+q2.P()*q2.P()*(1-cosqqp*cosqqp));
  Eh=TMath::Sqrt(ph*ph+mpi*mpi);
  *t=(-Q2-2*nu*Eh+2*q2.P()*ph*cosqqp+mpi*mpi);
  
  WB=sqrt(s);//pow(Mp,2.)+Q2*((1.-xB)/xB));
  QpCM=(s-Mp*Mp+mpi*mpi)/(2*WB);// Epicm
  nuCM=(s-Mp*Mp-Q2)/(2*WB);
   
 FacA=  sqrt(nuCM*nuCM+ Q2);//qCM0
 FacB= sqrt(QpCM*QpCM - pow(mpi,2));// qpiCM
  // *tmin=(-Q2-2*QpCM*(nuCM-qCM)+pow(mpi,2)); 
  *tmin= (-Q2-2*QpCM*nuCM+ 2*(FacA)*(FacB)+ pow(mpi,2));
  
}


bool flagFlipHelicity(float run){

  bool flip_flag = false;

  bool hw_flag = false;

  bool pr_flag = false;

  

  //--- half wave plate IN

  if((run>=12649 && run<=12661) || (run>=12931 && run<=13418) || (run>=14227 && run<=14516) || (run>=14522 && run<=15117)) hw_flag = true;

  //--- spin flip (precession in accelerator)

  if((run>=14150 && run<=14423) || (run>=12508 && run<=12661)) pr_flag = true;


  flip_flag = hw_flag ^ pr_flag;

  return flip_flag;

}

