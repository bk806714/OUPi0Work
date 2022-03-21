#include "TH1.h"
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

void xsec_modelcompare(int kin)
{
    
  const Int_t tbins=5;
  const Int_t phibins=12;
  const Int_t params =4;
  const Int_t GKpoints=12;
  Double_t sigmaLTsyst[tbins]={0,0,0,0,0};
  Double_t sigmaTTHsyst[tbins]={0,0,0,0,0};
  Double_t sigmaTTLsyst[tbins]={0,0,0,0,0};
  Double_t sigmaTTAvsyst[tbins];
  Double_t sigmaUsyst[tbins]={4,3.,0.,4.,2};// systematic bin to bin for kin 601
  Double_t sigmaLTpsyst[tbins]={9.3,3.0,9.1,13.5,7.0};
  Double_t ymin[tbins],ymax[tbins],yminU[tbins],ymaxU[tbins],yminLTp[tbins];
  Double_t yminTT[tbins],ymaxTT[tbins],ymaxLTp[tbins];
  Double_t GKepsilon,GKepsilonE;
  Double_t Q2;
  Double_t x_B;
  TString filepath="./";
  TString filepathDrp="/home/bishnu/Dropbox/pi0Scripts/GKPredction";
  Double_t ClusCorr;
  Double_t minT,minTT,minLT,minLTP;
  Double_t maxT,maxTT,maxLT,maxLTP;
  if(kin==601)
    {
      GKepsilon = 0.663;
      GKepsilonE = 0.663;
      Q2=5.54;
      x_B=0.60;
      ClusCorr=1.039;
      minT=0;
      maxT=120;
      
      minTT=-30;
      maxTT=5;
      
      minLT=-25;
      maxLT=5;
      
      minLTP=-3;
      maxLTP=5;
    }
  else if(kin==603)
    {
      GKepsilon = 0.495;
      GKepsilonE = 0.501;
      Q2=8.4;
      x_B=0.60;
      ClusCorr=1.036;
      minT=0;
      maxT=35;

      minLT=-8;
      maxLT=1;

       minTT=-10;
      maxTT=2;

      minLTP=-2;
      maxLTP=3;
      
    }
  ifstream experiment(Form("./xsecparams_%d.txt",kin));
  ifstream GKmodel(Form("./GK_kin%d.dat",kin));
  ifstream GKmodelE(Form("./GK_kin%d.dat",kin));

  for(Int_t i=0;i<tbins;i++)
    {
      sigmaTTAvsyst[i]=3.6;
      sigmaLTsyst[i] = 3.6;
      sigmaUsyst[i] = 3.6;/*TMath::Sqrt(pow(sigmaUsyst[i],2)+pow(4.15,2)+pow(sigmaUsyst[i],2) )*/;
      sigmaLTpsyst[i] = 3.75;
      sigmaTTAvsyst[i] = TMath::Sqrt(pow(sigmaUsyst[i],2)+pow( sigmaTTAvsyst[i],2));
    }

  Double_t xsec1[tbins],err1[tbins],xsec2[tbins],err2[tbins],xsec3[tbins],err3[tbins],xsec4[tbins],err4[tbins],dummyexp[tbins];//experiment
  //1==sigma_T+epsilon*sigma_L,2=sigma_TT,3=sigma_TL
  Double_t GKsigL[GKpoints],GKsigT[GKpoints],GKsigLTp[GKpoints],GKsigTT[GKpoints],GKsigLT[GKpoints],GKsigunseparated[GKpoints];
   Double_t GKsigLE[GKpoints],GKsigTE[GKpoints],GKsigLTpE[GKpoints],GKsigTTE[GKpoints],GKsigLTE[GKpoints],GKsigunseparatedE[GKpoints];
  Double_t t[tbins],tGKmodel[GKpoints],epsGKsigL[GKpoints],dummy;
 Double_t tGKmodelE[GKpoints],epsGKsigLE[GKpoints],dummyE;
 string titlelineData;
 std::getline(experiment,titlelineData);
 cout<<titlelineData<<endl;
  for(Int_t tbin=0;tbin<(tbins);tbin++)
    {
      experiment>>t[tbin]>>xsec1[tbin]>>err1[tbin]>>xsec2[tbin]>>err2[tbin]>>xsec3[tbin]>>err3[tbin]>>xsec4[tbin]>>err4[tbin];
      xsec1[tbin]=xsec1[tbin]*ClusCorr;
      xsec2[tbin]=xsec2[tbin]*ClusCorr;
      xsec3[tbin]=xsec3[tbin]*ClusCorr;
      xsec4[tbin]=xsec4[tbin]*ClusCorr;
      cout<<t[tbin]<<" & "<<setprecision(2)<<fixed<<xsec1[tbin]<<" $\\pm$ "<<err1[tbin]<<" $\\pm$ "<<0.01*sigmaTTAvsyst[tbin]*xsec1[tbin]<<" & "<<xsec2[tbin]<<" $\\pm$ "<<err2[tbin]<<" & "<<xsec3[tbin]<<" $\\pm$ "<<err3[tbin]<<" & "<<xsec4[tbin]<<" $\\pm$ "<<err4[tbin]<<"\t"<<endl;
      cout<<"\\hline"<<endl;
    }

	string titlelineGK1,titlelineGK2,titlelineGK3,titlelineGK4,titlelineGK5,titlelineGK6;
  std::getline(GKmodel,titlelineGK1);
  cout<<titlelineGK1<<endl;
  // std::getline(GKmodel,titlelineGK2);
  // cout<<titlelineGK2<<endl;
 
  for(Int_t i=0;i<GKpoints;i++)
    {
      GKmodel>>tGKmodel[i]>>dummy>>GKsigL[i]>>GKsigT[i]>>GKsigTT[i]>>GKsigLT[i]>>GKsigLTp[i];
      
      GKsigunseparated[i]=(GKsigT[i]+GKepsilon*GKsigL[i])*1.0;
      epsGKsigL[i] =GKepsilon*GKsigL[i];
       GKmodelE>>tGKmodelE[i]>>dummyE>>GKsigLE[i]>>GKsigTE[i]>>GKsigTTE[i]>>GKsigLTE[i]>>GKsigLTpE[i];
        GKsigunseparatedE[i]=(GKsigTE[i]+GKepsilonE*GKsigLE[i])*1.0;
      //  cout<<"============================================="<<endl;
        cout<<" GK sigma_T+epsilon*sigma_L = "<<GKsigunseparatedE[i]<<endl;
    }

  TGraph *grshadeLT = new TGraph(2*tbins);
  TGraph *grshadeLTp = new TGraph(2*tbins);
  TGraph *grshadeU = new TGraph(2*tbins);
  TGraph *grshadeTT = new TGraph(2*tbins);
  
  for (Int_t i=0;i<tbins;i++) 
    {
      ymin[i]=xsec2[i]-(xsec2[i]*sigmaLTsyst[i])/100.;
      yminU[i]=xsec1[i]-(xsec1[i]*sigmaUsyst[i])/100.;//
      yminTT[i]=xsec3[i]-(xsec3[i]*sigmaTTAvsyst[i])/100.;
      yminLTp[i]=xsec4[i]-(xsec4[i]*sigmaLTpsyst[i])/100.;
      ymax[i]=xsec2[i]+(xsec2[i]*sigmaLTsyst[i])/100.;//
      ymaxU[i]=xsec1[i]+(xsec1[i]*sigmaUsyst[i])/100.;
      ymaxTT[i]=xsec3[i]+(xsec3[i]*sigmaTTAvsyst[i])/100.;
      ymaxLTp[i]=xsec4[i]+(xsec4[i]*sigmaLTpsyst[i])/100.;
      //  cout<<ymin[i]<<"  "<<ymax[i]<<endl;
    }

  for (Int_t i=0;i<tbins;i++) 
    {
      grshadeLT->SetPoint(i,t[i],ymax[i]);
      grshadeLT->SetPoint(tbins+i,t[tbins-i-1],ymin[tbins-i-1]);
      grshadeLTp->SetPoint(i,t[i],ymaxLTp[i]);
      grshadeLTp->SetPoint(tbins+i,t[tbins-i-1],yminLTp[tbins-i-1]);
      grshadeTT->SetPoint(i,t[i],ymaxTT[i]);
      grshadeTT->SetPoint(tbins+i,t[tbins-i-1],yminTT[tbins-i-1]);
       grshadeU->SetPoint(i,t[i],ymaxU[i]);
       grshadeU->SetPoint(tbins+i,t[tbins-i-1],yminU[tbins-i-1]);
    }
  grshadeLT->SetFillStyle(3014);
  grshadeLT->SetFillColor(7);
  grshadeLTp->SetFillStyle(3025);
  grshadeLTp->SetFillColor(7);
  grshadeTT->SetFillStyle(3014);
  grshadeTT->SetFillColor(7);
  grshadeU->SetFillStyle(3001);
  grshadeU->SetFillColor(13);
 
  TString xsecparams[params]={"#sigma_{T}+#epsilon#sigma_{L}","#sigma_{LT}","#sigma_{TT}","#sigma_{LT'}"};
  TString xsecparamsnamefile[params]={"sigma_T_epsilon_sigma_L","sigma_LT","sigma_TT","sigma_LTp"};
  TMultiGraph *mgxsec1,*mgxsec2,*mgxsec3,*mgxsec4;

  mgxsec1=new TMultiGraph();
   mgxsec2=new TMultiGraph();
   mgxsec3=new TMultiGraph();
  mgxsec4=new TMultiGraph();


mgxsec1->SetTitle(Form(" ; t_{min} - t (GeV^{2});#frac{d#sigma_{T}+#epsilon#sigma_{L}}{dt} (nbarn.GeV^{-2})"));
 
  mgxsec2->SetTitle(Form("  ; t_{min} - t (GeV^{2});#frac{d#sigma_{TL}}{dt} (nbarn.GeV^{-2})"));
  
  mgxsec3->SetTitle(Form(" ; t_{min} - t (GeV^{2});#frac{d#sigma_{TT}}{dt} (nbarn.GeV^{-2})"));	

  mgxsec4->SetTitle(Form(" ; t_{min} - t (GeV^{2});#frac{d#sigma_{TL'}}{dt} (nbarn.GeV^{-2})"));
	       /*
  mgxsec1->SetTitle(Form(" Q^{2} =%.2f GeV^{2}, x_{B} =%.2f; t_{min} - t (GeV^{2});#frac{d#sigma_{T}+#epsilon#sigma_{L}}{dt} (nbarn.GeV^{-2})",Q2,x_B));
 
  mgxsec2->SetTitle(Form(" Q^{2} =%.2f GeV^{2}, x_{B} =%.2f ; t_{min} - t (GeV^{2});#frac{d#sigma_{TL}}{dt} (nbarn.GeV^{-2})",Q2,x_B));
  
  mgxsec3->SetTitle(Form(" Q^{2} =%.2f GeV^{2}, x_{B} =%.2f; t_{min} - t (GeV^{2});#frac{d#sigma_{TT}}{dt} (nbarn.GeV^{-2})",Q2,x_B));	

  mgxsec4->SetTitle(Form(" Q^{2} =%.2f GeV^{2}, x_{B} =%.2f; t_{min} - t (GeV^{2});#frac{d#sigma_{TL'}}{dt} (nbarn.GeV^{-2})",Q2,x_B));
	       */
  TGraphErrors *gr0exp,*gr0GK,*gr1exp,*gr1GK,*gr2GK,*gr2exp,*gr3GK,*gr3exp,*grGKsigT,*grGKepsigL,*grGKsigTE,*gr3GKE;;
  TGraphErrors *gr0GKE,*gr1GKE,*gr2GKE;
  
  TCanvas *c00,*c01,*c02,*c03;
  c00=new TCanvas();
  c00->cd();
  c00->SetGrid();
  gr0exp= new TGraphErrors(tbins,t,xsec1,0,err1);
  gr0exp->SetMarkerStyle(20);
  //  gr0exp->SetMarkerSize(2);
  gr0exp->SetMarkerColor(2);
  gr0exp->SetDrawOption("p");
  
  // gr0GK= new TGraphErrors(GKpoints,tGKmodel,GKsigunseparated,0,0);
  // gr0GK->SetDrawOption("L");
  // gr0GK->SetMarkerStyle(33);
  // gr0GK->SetMarkerColor(1);
  // gr0GK->SetLineWidth(3);
  // gr0GK->SetLineStyle(2);

  gr0GKE= new TGraphErrors(GKpoints,tGKmodelE,GKsigunseparatedE,0,0);
  gr0GKE->SetDrawOption("L");
  gr0GKE->SetMarkerStyle(33);
  gr0GKE->SetMarkerColor(4);
  gr0GKE->SetLineWidth(2);
  gr0GKE->SetLineStyle(2);
  gr0GKE->SetLineColor(2);
  
  grGKsigT= new TGraphErrors(GKpoints,tGKmodelE,GKsigTE,0,0);
  grGKsigT->SetDrawOption("p");
  grGKsigT->SetMarkerStyle(7);
  grGKsigT->SetMarkerSize(2);
  grGKsigT->SetMarkerColor(6);
  grGKsigT->SetLineWidth(3);
  grGKsigT->SetLineStyle(2);
  
  grGKepsigL= new TGraphErrors(GKpoints,tGKmodelE,epsGKsigLE,0,0);
  grGKepsigL->SetDrawOption("P");
  grGKepsigL->SetMarkerStyle(29);
  grGKepsigL->SetMarkerSize(2);
  grGKepsigL->SetMarkerColor(3);
  grGKepsigL->SetLineWidth(3);
  grGKepsigL->SetLineStyle(2);

  mgxsec1->SetMaximum(maxT);
  mgxsec1->SetMinimum(minT);
  mgxsec1->Add(grGKepsigL,"p");
  mgxsec1->Add(grGKsigT,"p");
  mgxsec1->Add(gr0exp,"p");
  //  mgxsec1->Add(gr0GK,"L");
  mgxsec1->Add(gr0GKE,"L");
  //  mgxsec1->Add(grshadeU,"f");
  mgxsec1->Draw("A");
  //  grshadeU->Draw("f");
  TLegend *leg1 = new TLegend(0.5,0.65,0.9,0.9);                             
  leg1->AddEntry(gr0exp,"d#sigma_{T}+#epsilon#sigma_{L} - This work","p");//
  leg1->AddEntry(grGKsigT,"d#sigma_{T} - GK model","p");
  leg1->AddEntry(grGKepsigL,"#epsilond#sigma_{L} - GK model","p");
  leg1->AddEntry(gr0GKE,"d#sigma_{T}+#epsilon#sigma_{L} - GK model","L");
  leg1->SetFillColor(kWhite);
  leg1->Draw();
  c00->Update();
  //  c00->Print(Form("./kin361output/%s_modelcompare.pdf",xsecparamsnamefile[0].Data()));
  //  c00->Print(Form("./kin361output/%s_modelcompare.png",xsecparamsnamefile[0].Data()));
  
  c01=new TCanvas();
  c01->cd();
  c01->SetGrid();
  gr1exp= new TGraphErrors(tbins,t,xsec2,0,err2);
  gr1exp->SetMarkerStyle(20);
  //  gr1exp->SetMarkerSize(2);
  gr1exp->SetMarkerColor(2);
  gr1exp->SetDrawOption("P");
  
  // gr1GK= new TGraphErrors(GKpoints,tGKmodel,GKsigLT,0,0);
  // gr1GK->SetDrawOption("P");
  // gr1GK->SetMarkerStyle(33);
  // gr1GK->SetMarkerColor(1);
  // gr1GK->SetLineWidth(3);
  // gr1GK->SetLineStyle(2);


  gr1GKE= new TGraphErrors(GKpoints,tGKmodelE,GKsigLTE,0,0);
  gr1GKE->SetDrawOption("P");
  gr1GKE->SetMarkerStyle(20);
  gr1GKE->SetMarkerColor(4);
  gr1GKE->SetLineColor(4);
  gr1GKE->SetLineWidth(2);
  gr1GKE->SetLineStyle(2);
  
  mgxsec2->SetMaximum(maxLT);
  mgxsec2->SetMinimum(minLT);
  mgxsec2->Add(gr1exp,"p");
  // mgxsec2->Add(gr1GK,"L");
  mgxsec2->Add(gr1GKE,"L"); 

   
  
  // mgxsec2->Add(grshadeLT,"f");
  mgxsec2->Draw("A");
  TLegend *leg2 = new TLegend(0.7,0.55,0.9,0.7);                             
  leg2->AddEntry(gr1exp,"This work","p");
  // leg2->AddEntry(gr1GK,"GK","L");
  leg2->AddEntry(gr1GKE,"GK model","L");
  leg2->SetFillColor(kWhite);
  //leg2->SetTextSize(0.04);
  leg2->Draw();
  c01->Update();
  //  c01->Print(Form("./kin361output/%s_modelcompare.pdf",xsecparamsnamefile[1].Data()));
  //  c01->Print(Form("./kin361output/%s_modelcompare.png",xsecparamsnamefile[1].Data()));
  
  c02=new TCanvas();
  c02->cd();
   c02->SetGrid();
  gr2exp= new TGraphErrors(tbins,t,xsec3,0,err3);
  gr2exp->SetMarkerStyle(20);
  // gr2exp->SetMarkerSize(2);
  gr2exp->SetMarkerColor(2);
  gr2exp->SetDrawOption("P");
  
  // gr2GK= new TGraphErrors(GKpoints,tGKmodel,GKsigTT,0,0);
  // gr2GK->SetDrawOption("P");
  // gr2GK->SetMarkerStyle(33);
  // gr2GK->SetMarkerColor(1);
  // gr2GK->SetLineWidth(3);
  // gr2GK->SetLineStyle(2);

  gr2GKE= new TGraphErrors(GKpoints,tGKmodelE,GKsigTTE,0,0);
  gr2GKE->SetDrawOption("P");
  gr2GKE->SetMarkerStyle(20);
  gr2GKE->SetMarkerColor(2);
  gr2GKE->SetLineColor(4);
  gr2GKE->SetLineWidth(2);
  gr2GKE->SetLineStyle(2);

  mgxsec3->SetMaximum(maxTT);
  mgxsec3->SetMinimum(minTT);
  mgxsec3->Add(gr2exp,"p");
  // mgxsec3->Add(gr2GK,"L");
  mgxsec3->Add(gr2GKE,"L");
  //  mgxsec3->Add(grshadeTT,"f");
  mgxsec3->Draw("A");
  TLegend *leg3 = new TLegend(0.7,0.75,0.9,0.9);                             
  leg3->AddEntry(gr2exp,"This work","p");
  // leg3->AddEntry(gr2GK,"GK","L");
  leg3->AddEntry(gr2GKE,"GK model","L");
  leg3->SetFillColor(kWhite);
  //leg3->SetTextSize(0.04);
  leg3->Draw();
  c02->Update();
  // c02->Print(Form("./kin361output/%s_modelcompare.pdf",xsecparamsnamefile[2].Data()));
  //  c02->Print(Form("./kin361output/%s_modelcompare.png",xsecparamsnamefile[2].Data()));

  c03=new TCanvas();
  c03->cd();
  c03->SetGrid();
  gr3exp= new TGraphErrors(tbins,t,xsec4,0,err4);
  gr3exp->SetMarkerStyle(20);
  // gr3exp->SetMarkerSize(2);
  gr3exp->SetMarkerColor(2);
  gr3exp->SetDrawOption("P");
  
  // gr3GK= new TGraphErrors(GKpoints,tGKmodel,GKsigLTp,0,0);
  // gr3GK->SetDrawOption("P");
  // gr3GK->SetMarkerStyle(33);
  // gr3GK->SetMarkerColor(1);
  // gr3GK->SetLineWidth(3);
  // gr3GK->SetLineStyle(2);

  gr3GKE= new TGraphErrors(GKpoints,tGKmodelE,GKsigLTpE,0,0);
  gr3GKE->SetDrawOption("P");
  gr3GKE->SetMarkerStyle(20);
  gr3GKE->SetMarkerColor(2);
  gr3GKE->SetLineColor(4);
  gr3GKE->SetLineWidth(2);
  gr3GKE->SetLineStyle(2);

  mgxsec4->SetMaximum(maxLTP);
  mgxsec4->SetMinimum(minLTP);
  mgxsec4->Add(gr3exp,"p");
  // mgxsec4->Add(gr3GK,"L");
  mgxsec4->Add(gr3GKE,"L");
  // mgxsec4->Add(grshadeLTp,"f");
  mgxsec4->Draw("A");

  
  TLegend *leg4 = new TLegend(0.7,0.75,0.9,0.9);                             
  leg4->AddEntry(gr3exp,"This work","p");
  // leg4->AddEntry(gr3GK,"GK","L");
  leg4->AddEntry(gr3GKE,"GK model","L");
  leg4->SetFillColor(kWhite);
  // leg4->SetTextSize(0.04);
  leg4->Draw();
  c03->Update();
  //  c03->Print(Form("./kin361output/%s_modelcompare.pdf",xsecparamsnamefile[3].Data()));
  // c03->Print(Form("./kin361output/%s_modelcompare.png",xsecparamsnamefile[3].Data()));

  /////////////////////////////////////////////////////////////////
  TMultiGraph *mgT;//,*mgxsec2,*mgxsec3,*mgxsec4;

  mgT=new TMultiGraph();
  mgT->SetTitle(Form(";t_{min}-t (GeV^{2});#frac{d#sigma}{dt} (nb GeV^{-2})"));
  TCanvas *cT0;//,*c01,*c02,*c03;
  cT0=new TCanvas();
  // mgT->Add(grGKepsigL,"p");
  // mgT->Add(grGKsigT,"p");
  mgT->Add(gr0exp,"p");
  mgT->Add(gr0GKE,"L");
  mgT->Add(grshadeU,"f");


   mgT->Add(gr1exp,"p");
   mgT->Add(gr1GKE,"L");


   mgT->Add(gr2exp,"p");
   mgT->Add(gr2GKE,"L");
   // mgT->Add(grshadeTT,"f");


     mgT->Add(gr3exp,"p");
     mgT->Add(gr3GKE,"L");
   // mgT->Add(grshadeLTp,"f");
 mgT->SetMaximum(120);
  mgT->SetMinimum(-40);
   
   mgT->Draw("A");
   mgT->GetXaxis()->SetRangeUser(0,0.8);
   TLegend *le1 = new TLegend(0.6,0.65,0.8,0.9);                             
   le1->AddEntry(gr0exp,"d#sigma_{T}+#epsilon#sigma_{L}","p");//
   le1->AddEntry(gr1exp,"d#sigma_{TL}","p");
   le1->AddEntry(gr2exp,"d#sigma_{TT}","p");
    le1->AddEntry(gr3exp,"d#sigma_{TL'}","p");
   
  le1->SetFillColor(kWhite);
  //le1->SetTextSize(0.06);
  le1->Draw();
  // c00->Update();
  TLatex * la= new TLatex(0.6,100,Form("#splitline{Q^{2}=%.2f GeV^{2}}{x_{B}=%.2f}",Q2,x_B));
  la->SetTextSize(0.05);
  la->Draw();
  
  
}


void DrawCol()
{
   Int_t i,n;
   Double_t x,y;
   TLatex *l;
   Double_t colors[8]={6,9,5,13,1,46,7,2};

   TGraph *g = (TGraph*)gPad->GetListOfPrimitives()->FindObject("Graph");
   n = g->GetN();
   TMarker *m;
   for (i=0; i<n; i++) {
      g->GetPoint(i,x,y);
      m = new TMarker(x,y,20);
      m->SetMarkerColor(colors[i]);;
      m->Paint();
   }
}
