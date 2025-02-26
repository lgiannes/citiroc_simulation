



double HGol(double *x,double *p){
    double s = x[0];
    double wc1 = p[0];
    double Aol1 = p[1];

    double Fs_HGol = Aol1 * wc1/(s+wc1); // HG preamp
    return  abs(Fs_HGol);

}

double HGcl(double *x,double *p){
    double s = x[0];
    double wc1 = p[0];
    double Aol1 = p[1];
    double beta = p[2];

    double Fs_HG = -Aol1 * beta/( 1+beta+Aol1+s/wc1*(1+beta) ); // HG preamp
    return  abs(Fs_HG);

}

double FSH(double *x,double *p){
    double s = x[0];
    double wc2 = p[0];
    double Aol2 = p[1];

    double tau = 2.5e-9; // 400 kHz feedback circuit
    double R_F = 25e3;
    double R_in = 500;

    double Fs_FSH = -( Aol2 * R_F/R_in * tau * s)/
        ( Aol2 * pow((tau*s+1),2) + ( (R_F/R_in+2) * tau * s +tau*tau*s*s + 1 )*(1+s/wc2) ) ;

    return  abs(Fs_FSH);

}
double F(double *x,double *p){
    // HG preamp
    double wc1 = p[0];
    double Aol1 = p[1];
    // FSH
    double wc2 = p[2];
    double Aol2 = p[3];

    // HG loopback
    double beta = p[4];

    double p_HG[4] = {wc1,Aol1,beta};
    double p_FSH[2] = {wc2,Aol2};

    return HGcl(x,p_HG) * FSH(x,p_FSH); 
}

void transfer_function_fsh(){
    

    double Aol1 = pow(10,3.25); // 65 dB 
    double Aol2 = 1000; // 60 dB

    double wc1 = 450e3*2*ROOT::Math::Pi(); // 450 kHz cutoff frequency for the HG opa in OL
    double wc2 = 1.2e6*2*ROOT::Math::Pi(); // 1.2 MHz cutoff frequency for the FSH opa in OL

    double HG_Cin = 15e-12; // 15 pF
    int HG_DAC = 44;
    double HG_Cfb = 25e-15*(63-HG_DAC); // 1.5 pF
    double beta = HG_Cin/HG_Cfb; // feedback factor for the HG opa in CL

    cout<<"HG Gain: "<<beta<<endl;

    TF1 *f = new TF1("f",F,1e1,1e12,5);
    f->SetParameters(wc1,Aol1,wc2,Aol2,beta);
    f->SetLineColor(kRed);
    f->SetLineWidth(3);
    f->SetNpx(1000);

    TCanvas *c = new TCanvas("c","c",800,600);
    f->Draw();
    c->SetLogx();
    c->SetLogy();
    f->GetXaxis()->SetTitle("f [Hz]");
    f->GetYaxis()->SetTitle("Gain");

    TF1 *f_HGol = new TF1("f_HGoL",HGol,1e1,1e12,2);
    f_HGol->SetParameters(wc1,Aol1);
    f_HGol->SetLineColor(kBlue);
    f_HGol->Draw("same");

    TF1 *f_HGcl = new TF1("f_HGcl",HGcl,1e1,1e12,3);
    f_HGcl->SetParameters(wc1,Aol1,beta);
    f_HGcl->SetLineColor(kGreen);
    f_HGcl->Draw("same");

    TF1 *f_FSH = new TF1("f_FSH",FSH,1e1,1e12,2);
    f_FSH->SetParameters(wc2,Aol2);
    f_FSH->SetLineColor(kBlack);
    f_FSH->Draw("same");

    TLegend *leg = new TLegend(0.1,0.7,0.3,0.9);
    leg->AddEntry(f,"HG (CL) x FSH (CL)","l");
    leg->AddEntry(f_HGol,"HG OL","l");
    leg->AddEntry(f_HGcl,"HG CL","l");
    leg->AddEntry(f_FSH,"FSH CL","l");
    leg->Draw();





    


}