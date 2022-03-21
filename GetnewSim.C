/////////////////////////
//  CHANGE THE THRESHOLD //
// AND FACTOR *100 OR 10 

#define me 0.000511
#define Mp 0.938272
#define mpi 0.13497
void GetTandTmin(Float_t E0,Float_t kpx, Float_t kpy,Float_t kpz, Float_t xB,TLorentzVector pi0,Float_t *t,Float_t *tmin);
Double_t Gamma_factor(Double_t Q2,Double_t xB, double E0);
Double_t epsilon_factor(Double_t Q2,Double_t xB, double E0);
void GetnewSim(Int_t kin)
{
  Double_t rcut,E0,Corr;
  Double_t edgecut=3.0;
 
 
  //////////////////////////////////////////
  Double_t Minv=0.134977;
  int totEvent;
   if(kin==601){
   
    rcut=0.004;
    E0= 8.51849;// energy
    totEvent=1000;
   
    }
  if(kin==603){
   
    //  Corr=-20.64;
    rcut=0.004;
    E0=10.591;//;//10.6171;// energy
    totEvent=1900;
    
  }
 //////////////CHECK THIS PLEASE ///////
  TH1F *histNew= new TH1F ("histNew","",100,-0.5,2.0);
  TH1F *hist= new TH1F ("hist","",100,-0.5,2.0);
   
  Float_t photon2ene=1.1;///// PHOTON THRESHOLD
  Float_t enecut=photon2ene;
  Int_t Thresh;
  Thresh=enecut*10;
  // TString filename=(Form("/home/bishnu/pi0/kin%dClus/threshE%d/dvcs_smeared_simu_%d_%dM_%dThresh.root",kin,Thresh,kin,totEvent,Thresh));
  // TString filename=(Form("/home/bishnu/pi0/kin%dClus/dvcs_smeared_simu_%d_%dM_%dThresh.root",kin,kin,totEvent,Thresh));
   TString filename=(Form("/home/bishnu/pi0/simFiles/trialFiles/dvcs_smeared_simu_%d_%dM_%dThresh.root",kin,totEvent,Thresh));
   TChain *simu = new TChain("tr1");
   if(kin==603){
   simu->Add("dvcs_sim_new_smeared_603_1900M_corr_Plustheta_newDB_tweak098.root");
   simu->Add("dvcs_sim_new_smeared_603_1900M_corr_Plustheta_newDB_tweak098_1.root");
   }
   else if(kin==601)
     {
 simu->Add("dvcs_sim_new_smeared_601_1000M_corr_Plustheta_newDB_tweak098.root");
 simu->Add("dvcs_sim_new_smeared_601_1000M_corr_Plustheta_newDB_tweak098_1.root");
     }
  
     


   cout<<"this is threshold "<<Thresh<<endl;
   

   
   TFile *froot = new TFile(filename,"recreate");
   TTree *tr1 = new TTree("tr1","tr1");
   Double_t phisim_rec,phisim_ver,psf,kz_v,ky_v,kx_v,kpx_v;
   Float_t tsim_rec,tsim_min_rec,tsim_ver,tsim_min_ver;
   Float_t tsim_min_recNew,tsim_min_verNew;
   Double_t kpy_v,kpz_v,xB_v,Q2_v,Q2;
   /////////////////////////////////////
   
   /////////////////////////////////////
   Double_t kpy_m,kpz_m,kpx_m,xB_m,Q2m;
   double q1x_v,q1y_v, q1z_v,q2x_v,q2y_v, q2z_v;
   double q1x_m,q1y_m, q1z_m,q2x_m,q2y_m, q2z_m,mm2BS,minvBS;
   double ene1_m,ene2_m,rval,xc1,xc2,yc1,yc2,v_v;
   TLorentzVector gamma1_m,gamma2_m,pi0_m,gamma1_v,gamma2_v,pi0_v;
  kx_v=0.;ky_v=0.;
  // TCut GoodVertexsimu=Form("v_v-%d>=-6.5&&v_v-%d<=6.5",0,0);//0.75,0.75);
  //  TCut Goodelectronsimu=Form("rval>%f",rcut);
  // TCut goodphoton=Form("((ene1>=%f&&ene2>=%f)||(ene1>=%f&&ene2>=%f))&&xc1<=(14.-%f)&&xc1>=(-23.4+%f)&&yc1<=(24.-%f)&&yc1>=(-24.+%f)&&xc2<=(14.-%f)&&xc2>=(-23.4+%f)&&yc2<=(24.-%f)&&yc2>=(-24.+%f)",enecut,photon2ene,photon2ene,enecut,edgecut,edgecut,edgecut,edgecut,edgecut,edgecut,edgecut,edgecut);
  
  
  simu->SetBranchAddress("phi",&phisim_rec);
  simu->SetBranchAddress("phi_v",&phisim_ver);
  simu->SetBranchAddress("psf",&psf);
  
  simu->SetBranchAddress("kz_v",&kz_v);
  simu->SetBranchAddress("kpx_v",&kpx_v);
  simu->SetBranchAddress("kpy_v",&kpy_v);
  simu->SetBranchAddress("kpz_v",&kpz_v);
  simu->SetBranchAddress("q2_v",&Q2_v);
  simu->SetBranchAddress("xB_v",&xB_v);
  simu->SetBranchAddress("xB",&xB_m);
  //////////////////////////////////////
  simu->SetBranchAddress("kpx",&kpx_m);
  simu->SetBranchAddress("kpy",&kpy_m);
  simu->SetBranchAddress("kpz",&kpz_m);
  simu->SetBranchAddress("ene1",&ene1_m);
  simu->SetBranchAddress("ene2",&ene2_m);

  simu->SetBranchAddress("q1x",&q1x_m);
  simu->SetBranchAddress("q1y",&q1y_m);
  simu->SetBranchAddress("q1z",&q1z_m);
  simu->SetBranchAddress("q1x_v",&q1x_v);
  simu->SetBranchAddress("q1y_v",&q1y_v);
  simu->SetBranchAddress("q1z_v",&q1z_v);

  simu->SetBranchAddress("q2x",&q2x_m);
  simu->SetBranchAddress("q2y",&q2y_m);
  simu->SetBranchAddress("q2z",&q2z_m);
  simu->SetBranchAddress("q2x_v",&q2x_v);
  simu->SetBranchAddress("q2y_v",&q2y_v);
  simu->SetBranchAddress("q2z_v",&q2z_v);
  simu->SetBranchAddress("mm2",&mm2BS);
  simu->SetBranchAddress("m",&minvBS);
  simu->SetBranchAddress("rval",&rval);
  simu->SetBranchAddress("xc1",&xc1);
  simu->SetBranchAddress("xc2",&xc2);
  simu->SetBranchAddress("yc1",&yc1);
  simu->SetBranchAddress("yc2",&yc2);
  simu->SetBranchAddress("ene1",&ene1_m);
  simu->SetBranchAddress("ene2",&ene2_m);
   simu->SetBranchAddress("v",&v_v);
 
  


  Double_t gammasimu,eps,E0s;
  TBranch *rg1=tr1->Branch("mm2",&mm2BS);	 
  TBranch *rg2=tr1->Branch("m",&minvBS);
  TBranch *rg3=tr1->Branch("tRec",&tsim_rec);	 
  TBranch *rg4=tr1->Branch("tmin_Rec",&tsim_min_rec);
  TBranch *rg5=tr1->Branch("t_v",&tsim_ver);	 
  TBranch *rg6=tr1->Branch("t_min",&tsim_min_ver);
  
  TBranch *rg7=tr1->Branch("phi",&phisim_rec);	 
  TBranch *rg8=tr1->Branch("phi_v",&phisim_ver);
  TBranch *rg9=tr1->Branch("psf",&psf);	 
  TBranch *rg10=tr1->Branch("xB",&xB_m);
  
  TBranch *rg11=tr1->Branch("q2_v",&Q2_v);	 
  TBranch *rg12=tr1->Branch("xB_v",&xB_v);
  TBranch *rg13=tr1->Branch("q2",&Q2);	 
  TBranch *rg14=tr1->Branch("xB",&xB_m);
  TBranch *rg15=tr1->Branch("epsilon",&eps);	 
  TBranch *rg16=tr1->Branch("gammaFac",&gammasimu);
  TBranch *rg17=tr1->Branch("eE_v",&E0s);

  TBranch *rg18=tr1->Branch("xc1",&xc1);
  TBranch *rg19=tr1->Branch("xc2",&xc2);	 
  TBranch *rg20=tr1->Branch("yc1",&yc1);
  TBranch *rg21=tr1->Branch("yc2",&yc2);
  TBranch *rg22=tr1->Branch("rval",&rval);
  TBranch *rg23=tr1->Branch("ene2",&ene2_m);
  TBranch *rg24=tr1->Branch("ene1",&ene1_m);

  TBranch *rg25=tr1->Branch("t_minNew",&tsim_min_verNew);
  TBranch *rg26=tr1->Branch("tmin_RecNew",&tsim_min_recNew);
  
  
 
 // simu->SetBranchAddress("tRec",&tsim_rec);
  // simu->SetBranchAddress("tmin_Rec",&tsim_min_rec);
  // simu->SetBranchAddress("t_v",&tsim_ver);
  // simu->SetBranchAddress("tmin_v",&tsim_min_ver);
  Bool_t GoodVertexsimu,Goodelectronsimu,goodphoton,goodEne,goodPos;
   for(Int_t i=0;i<simu->GetEntries();i++)
    {
      simu->GetEntry(i);
      GoodVertexsimu=v_v>=-6.5&&v_v<=6.5;//0.75,0.75);
      Goodelectronsimu=rval>rcut;
      goodEne=((ene1_m>=enecut&&ene2_m>=photon2ene)||(ene1_m>=photon2ene&&ene2_m>=enecut));
      goodPos=xc1<=(14.-edgecut)&&xc1>=(-23.4+edgecut)&&yc1<=(24.-edgecut)&&yc1>=(-24.+edgecut)&&xc2<=(14.-edgecut)&&xc2>=(-23.4+edgecut)&&yc2<=(24.-edgecut)&&yc2>=(-24.+edgecut);
      goodphoton=goodPos && goodEne;
      if(goodphoton && GoodVertexsimu && Goodelectronsimu)
	{
      TVector3 eini(kx_v,ky_v,kz_v);
      TVector3 eprime(kpx_v,kpy_v,kpz_v);
       E0s=sqrt(kx_v*kx_v + ky_v*ky_v + kz_v*kz_v);
      //double E0_m=sqrt(kx_m*kpx_m + kpy_m*kpy_m + kpz_m*kpz_m);// energy measured
      float ene1_v=sqrt(q1x_v*q1x_v+q1y_v*q1y_v+q1z_v*q1z_v);
      float ene2_v=sqrt(q2x_v*q2x_v+q2y_v*q2y_v+q2z_v*q2z_v);
      // htrecbinB->Fill(ene1_v-ene1_m);
      
      gamma1_v.SetPxPyPzE(q1x_v,q1y_v,q1z_v,ene1_v);
      gamma2_v.SetPxPyPzE(q2x_v,q2y_v,q2z_v,ene2_v);

      gamma1_m.SetPxPyPzE(q1x_m,q1y_m,q1z_m,ene1_m);
      gamma2_m.SetPxPyPzE(q2x_m,q2y_m,q2z_m,ene2_m);

      pi0_m=gamma1_m+gamma2_m;
      pi0_v=gamma1_v+gamma2_v;
      float kpx_vU=(float)kpx_v; float kpy_vU=(float)kpy_v; float kpz_vU=(float)kpz_v;
      float kpx_mU=(float)kpx_m; float kpy_mU=(float)kpy_m; float kpz_mU=(float)kpz_m;
      float xB_vU=(float)xB_v; float xB_mU=(float)xB_m;
      float E0sVU=(float)E0s;
      GetTandTmin(E0sVU, kpx_vU,  kpy_vU, kpz_vU, xB_vU, pi0_v, &tsim_ver, &tsim_min_ver);//,&tsim_min_verNew );// vertex
	GetTandTmin(E0,  kpx_mU,  kpy_mU, kpz_mU, xB_mU, pi0_m, &tsim_rec, &tsim_min_rec);//, &tsim_min_recNew );
     gammasimu = Gamma_factor(Q2_v,xB_v,E0s);
     eps = epsilon_factor(Q2_v,xB_v,E0s);
     histNew->Fill(tsim_min_recNew-tsim_rec);
     hist->Fill(tsim_min_rec-tsim_rec);
      tr1->Fill();
	}
}
 tr1->Write();
 froot->Close();

TCanvas *c1= new TCanvas();
hist->Draw();
hist->SetLineColor(4);
histNew->Draw("same");
histNew->SetLineColor(2);
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
  QpCM=(s-Mp*Mp+mpi*mpi)/(2*WB); // SAME AS EPI'CM
  nuCM=(s-Mp*Mp-Q2)/(2*WB);
  FacA=  sqrt(nuCM*nuCM+ Q2);//qCM0
  FacB= sqrt(QpCM*QpCM - pow(mpi,2));// qpiCM
 *tmin= (-Q2-2*QpCM*nuCM+ 2*(FacA)*(FacB)+ pow(mpi,2));

 // qCM=sqrt(Q2+(pow((s-Mp*Mp-Q2),2))/(4*s));
  //  *tminold=(-Q2-2*QpCM*(nuCM-qCM)+pow(mpi,2));


 
}
Double_t epsilon_factor(Double_t Q2,Double_t xB, double E0)
{
  Double_t nu = Q2/(2*Mp*xB);
  Double_t q = sqrt(pow(nu,2)+Q2);
  Double_t y=nu/E0;
  //Double_t eps = 1/(1+2*pow(q*tan(theta_e/2),2)/Q2);
  Double_t eps = (1.-y-(Q2/(4*pow(E0,2))))/(1.-y+pow(y,2)/2+(Q2/(4*pow(E0,2))));//From Maxime
  return eps;
}
Double_t Gamma_factor(Double_t Q2,Double_t xB,double E0)
{
  Double_t nu = Q2/(2*Mp*xB);
  Double_t q = sqrt(pow(nu,2)+Q2);
  Double_t s = pow(Mp,2)+Q2*((1.-xB)/xB);
  Double_t kgamma = (s-pow(Mp,2))/(2*Mp);
  Double_t alphaQED = 0.0072982;
  Double_t eps=epsilon_factor(Q2,xB,E0);
  Double_t Gamma = (alphaQED*Q2*(1.-xB))/(8.*TMath::Pi()*pow(Mp*E0,2)*pow(xB,3)*(1.-eps));//From Maxime, makes dimensional sense
  // double gammaTest=((alphaQED*(1./xB-Q2))/(16.*pow(TMath::Pi(),2)))*pow(Mp*E0,2)*Q2*(1.-eps);
  // cout<<" Gamma     "<<Gamma<<"    there here and everywhere"<<gammaTest<<endl;
  return Gamma;
}
/////////////////////////////NEW TMIN FROM CHARLES/////////////
