

void plot_pairpt() {
    TFile *fin = new TFile( "GENV2_TEST.root" );

    TH1* hcoh = ((TH1*)fin->Get( "coh_pt" ));
    hcoh->Draw();
    hcoh->GetXaxis()->SetRangeUser(0, 0.25);
    ((TH1*)fin->Get( "inc_pt" ))->Draw("same");
    ((TH1*)fin->Get( "pimu_pt" ))->Draw("same");
    ((TH1*)fin->Get( "pipi_pt" ))->Draw("same");
    gPad->SetLogy(1);
}