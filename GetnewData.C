/////////////////////////////
// CHANGE PHOTON THRESHOLD///
// FILES GOES TO PION DIRECT//


#include "/home/bishnu/pi0/accidental_subtraction_function.h"
#define me 0.000511
#define Mp 0.938272
#define mpi 0.13497

void GetTandTmin(Float_t E0,Float_t kpx, Float_t kpy,Float_t kpz, Float_t xB,TLorentzVector pi0,Float_t *t,Float_t *tmin);

void GetnewData(Int_t kin,Int_t mode)
{
  TChain *ntu2 = new TChain("ntu2");
  TChain *ntuacp = new TChain("ntu2"); 
  TChain *ntuacc =  new TChain("ntu2");
  TChain *ntuaaa = new TChain("ntu2");
  Double_t rcut,E0,Corr,targetoffset;
  Double_t edgecut=3.0;
  

  int Thresh;
  //////////////CHECK HERE /////////////////////////
   /// CHANGE THE  PHOTON THRESHOLD
  Float_t photon2ene=1.1;
  Float_t enecut=photon2ene;
  // int opt=0; ///// for 0 = photonE = clustering
   int opt=1; ///// for 1 = photonE = 1.5
   /* 
    if(opt==0){
       Thresh=enecut*10;
    }
    else if (opt==1)
      {
	Thresh=10;  // from which folder the files to be read
      }
 
   */
 
    
 

 ////////////////////////////////////////////////////

 
  Double_t Minv=0.134977;
  Float_t phi,hel,mm3;
  int totEvent;
   TString filename;
   if(kin==601){
    w1 = 2.97;
    w2 = 4.08;
    rcut=0.004;
    E0= 8.51849;// energy
    totEvent=1000;
    Corr=-13.7;
    targetoffset=0.0004;
    Thresh=8;
    }
  if(kin==603){
    w1 = 2.64;
    w2 = 3.63;
    Corr=-17.64;
    rcut=0.004;
    E0=10.591;//;//10.6171;// energy
    totEvent=1900;
    targetoffset=0.0004;
    Thresh=10;
  }
  // TString file_path=Form("/media/bishnu/0172A59E02A74CC1/WF/Kin%d/threshE%d",kin,Thresh);// where the files is
   TString file_path=Form("/media/bishnu/0172A59E02A74CC1/WF/Kin%d/3Clus",kin);// where the files is
   cout<<" THRESH   "<<Thresh<<"    path      "<<file_path<<endl;
 TLorentzVector gamma1,gamma2,pi0;
 vector<int> Vecruns;
  //const int maxrun=334; //234 for601
  int run;
  bool flip_flag;
  int i=0;
   Int_t kin1=kin/10;
  int kin2=kin%10;
 ifstream runlist;
 runlist.open(Form("/home/bishnu/pi0/QualityRunlist/list_run_kin%d_%d.txt",kin1,kin2));
  while (!runlist.fail() && !runlist.eof())
{
      runlist >> run;
      Vecruns.push_back(run);
    
}
  int count=0;
  
  for(int it=0;it<Vecruns.size()-1;it++)
   {
     if (!gSystem->AccessPathName(Form("%s/ntu_-3_3_%d.root",file_path.Data(),Vecruns[it])))
       {
	 count++;
	 if(mode==0){
	 	 ntu2->Add(Form("%s/ntu_-3_3_%d.root",file_path.Data(),Vecruns[it]));
		 filename="-3_3";
	 }
	 else if(mode==2){
	 	 ntu2->Add(Form("%s/ntu_acc2_%d.root",file_path.Data(),Vecruns[it]));
	 filename="acc2";	 
	 }
	 else if(mode==3){
	   	 ntu2->Add(Form("%s/ntu_acc3_%d.root",file_path.Data(),Vecruns[it]));
	 filename="acc3";	 
	 }
	 else if(mode==1){
	 	 ntu2->Add(Form("%s/ntu_-11_-5_%d.root",file_path.Data(),Vecruns[it]));
	 filename="-11_-5";	 
       }
	 else cout<<"Wrong mode "<<endl;
	}
     else cout<<"Run  "<<Vecruns[it]<< "  is not present"<<endl;
}
  cout<<"Main coincidence window event" << kin<<" Total runs   "<<count<<endl;

 ntu2->SetBranchAddress("m",&m);
 ntu2->SetBranchAddress("mm2",&mm2);
 ntu2->SetBranchAddress("q2",&qs);
 ntu2->SetBranchAddress("rvalAlexa",&rval);
 ntu2->SetBranchAddress("q1x",&q1x_f);
 ntu2->SetBranchAddress("q1y",&q1y_f); 
 ntu2->SetBranchAddress("q1z",&q1z_f);  
 ntu2->SetBranchAddress("q2x",&q2x_f);
 ntu2->SetBranchAddress("q2y",&q2y_f); 
 ntu2->SetBranchAddress("q2z",&q2z_f);   
 ntu2->SetBranchAddress("ene1",&ene1); 
 ntu2->SetBranchAddress("ene2",&ene2);
 ntu2->SetBranchAddress("ntr",&ntr);
 ntu2->SetBranchAddress("xc1",&xc1);
 ntu2->SetBranchAddress("xc2",&xc2);
 ntu2->SetBranchAddress("yc1",&yc1);
 ntu2->SetBranchAddress("yc2",&yc2); 
 ntu2->SetBranchAddress("vz",&v_f);// vertex
 ntu2->SetBranchAddress("t1moy",&t1moy);
 ntu2->SetBranchAddress("t2moy",&t2moy);
 ntu2->SetBranchAddress("nclustu1",&u1clst);
 ntu2->SetBranchAddress("nclustu2",&u2clst);
 ntu2->SetBranchAddress("nclustv1",&v1clst);
 ntu2->SetBranchAddress("nclustv2",&v2clst);
 ntu2->SetBranchAddress("cer",&cer);
 ntu2->SetBranchAddress("pr1",&pr1);
 ntu2->SetBranchAddress("pr2",&pr2);
 ntu2-> SetBranchAddress("xB",&xB);
 ntu2->SetBranchAddress("kpx",&kpx);
 ntu2->SetBranchAddress("kpy",&kpy);
 ntu2->SetBranchAddress("kpz",&kpz);
 ntu2->SetBranchAddress("phi",&phi);
 ntu2->SetBranchAddress("hel",&hel);

 TString out_filename;
 if(opt==1)
   //  out_filename=(Form("/home/bishnu/pi0/kin%dClus/threshE%d/pi0_%s_%d_%dThresh_E15.root",kin,Thresh,filename.Data(),kin,Thresh));
  out_filename=(Form("/home/bishnu/pi0/simFiles/trialFiles/pi0_%s_%d_New.root",filename.Data(),kin));
 else if (opt==0)
    out_filename=(Form("/home/bishnu/pi0/kin%dClus/threshE%d/pi0_%s_%d_%dThresh.root",kin,Thresh,filename.Data(),kin,Thresh));
 else cout<<" K BHAIL GAVA    "<<endl;
  TFile *froot = new TFile(out_filename,"recreate");
  TTree *ntu = new TTree("ntu","ntu");
  TBranch *rg1=ntu->Branch("mm2",&mm2);	 
  TBranch *rg2=ntu->Branch("m",&m);
  TBranch *rg3=ntu->Branch("t",&t);	 
  TBranch *rg4=ntu->Branch("tmin",&tmin);
  TBranch *rg5=ntu->Branch("phi",&phi);	 
  TBranch *rg6=ntu->Branch("xB",&xB);
  TBranch *rg7=ntu->Branch("hel",&hel);
  TBranch *rg8=ntu->Branch("q2",&qs);	 
  TBranch *rg9=ntu->Branch("xc1",&xc1);
  TBranch *rg10=ntu->Branch("xc2",&xc2);	 
  TBranch *rg11=ntu->Branch("yc1",&yc1);
  TBranch *rg12=ntu->Branch("yc2",&yc2);
  TBranch *rg13=ntu->Branch("rvalAlexa",&rval);
  TBranch *rg14=ntu->Branch("mm3",&mm3);
  TBranch *rg15=ntu->Branch("flip",&flip_flag);
 
   Bool_t GoodVertex;
  cout<<"Entries from true = "<<ntu2->GetEntries()<<endl;
  
  for (Int_t i=0;i<ntu2->GetEntries();i++){
    ntu2->GetEntry(i);
    VDC0M4SClustTrk = (u1clst==1 && u2clst==1 && v1clst==1 && v2clst==1 );
    VDC1M3SClustTrk = (u1clst>1 && u2clst==1 && v1clst==1 && v2clst==1) || (u1clst==1 && u2clst>1 && v1clst==1 && v2clst==1) || (u1clst==1 && u2clst==1 && v1clst>1 && v2clst==1) || (u1clst==1&& u2clst==1   && v1clst==1 && v2clst>1);
    
    GoodSingleTrk = (VDC0M4SClustTrk || VDC1M3SClustTrk);
    
    PID = cer>=150 && (pr1/w1+pr2/w2)>600 && pr1/w1>200;
    
    acceptance =rval>rcut;
    
    Goode=PID  && GoodSingleTrk && acceptance;
    goodphotonE=(ene1>=enecut && ene2>=photon2ene)/*||(ene1>photon2ene&&ene2>enecut)*/  ;
    photonPos=xc1<=(14.-edgecut)&&xc1>=(-23.4+edgecut)&&yc1<=(24.-edgecut)&&yc1>=(-24.+edgecut)&&xc2<=(14.-edgecut)&&xc2>=(-23.4+edgecut)&&yc2<=(24.-edgecut)&&yc2>=(-24.+edgecut);
    goodphoton=goodphotonE && photonPos;
    mm3=(mm2-Corr*(m-Minv));
    GoodVertex=(v_f-targetoffset)<=0.065 && (v_f-targetoffset)>=-0.065;
    if(Goode && goodphoton &&GoodVertex )
 {
 
       gamma1.SetPxPyPzE(q1x_f,q1y_f,q1z_f,ene1);
       gamma2.SetPxPyPzE(q2x_f,q2y_f,q2z_f,ene2);
       pi0=gamma1+gamma2;
       GetTandTmin( E0, kpx,  kpy, kpz, xB, pi0, &t, &tmin);
       ntu->Fill();
 }
 
}
ntu->Write();
froot->Close();
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
  
  ////////This is the CHARLES's defination of TMIN////////////////
  WB=sqrt(s);
  QpCM=(s-Mp*Mp+mpi*mpi)/(2*WB);
  nuCM=(s-Mp*Mp-Q2)/(2*WB);
  FacA=  sqrt(nuCM*nuCM+ Q2);//qCM0
  FacB= sqrt(QpCM*QpCM - pow(mpi,2));// qpiCM
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
