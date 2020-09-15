

TProfile * makeProfile( TH2*h2D, string profileName, int nRebin = 1 ){
    auto profileH = h2D->ProfileX(profileName.c_str());

    // this does the right thing, checked (maintains average)
    profileH->RebinX( nRebin );
    profileH->Scale(2.0);
    return  profileH;
}

void proj( float ptmin = 0.7, float ptmax = 0.9 ){

    TFile *fin = new TFile( "GENV2_TEST.root" );

    TH3 * h3 = (TH3*)fin->Get( "c2ptmass" );

    h3->GetXaxis()->SetRangeUser(ptmin, ptmax);
    TH2 * h2 = (TH2*)h3->Project3D( "zy" );
    h2->Draw("colz");

    TH1* h = makeProfile( h2, "TEMP", 10 );
    h->Fit( "pol0", "R", "", 0.15, 0.25 );
}