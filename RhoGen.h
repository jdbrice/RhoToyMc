#ifndef RHO_GEN_H
#define RHO_GEN_H

#include "HistoAnalyzer.h"

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TProfile.h"

#include "loguru.h"

TF1 * f1FFQ = nullptr;
TF1 * fspin1 = nullptr;

template <class T>
class MinMaxPair {
public:
    MinMaxPair( ) {
        max = std::numeric_limits<T>::max();
        min = std::numeric_limits<T>::min();
    }
    MinMaxPair( T _min, T _max ) : max(_max), min(_min) {
    }
    ~MinMaxPair() {}
    
    T min;
    T max;

    void set( T _min, T _max ) {
        this->min = _min;
        this->max = _max;
    }

    bool fully_in( T _v ){
        if ( _v >= this->max ) return false;
        if ( _v <= this->min ) return false;
        return true;
    }

    bool in( T _v ){
        if ( _v > this->max ) return false;
        if ( _v < this->min ) return false;
        return true;
    }

};

// Cannot get partial template specialization to work :/
template <>
MinMaxPair<float> XmlConfig::get<MinMaxPair<float>>( string path ) const {
    MinMaxPair<float> mmp;
    mmp.min = get<float>( path + ":min", std::numeric_limits<float>::min() );
    mmp.max = get<float>( path + ":max", std::numeric_limits<float>::max() );
    return mmp;
}
template <>
MinMaxPair<float> XmlConfig::get<MinMaxPair<float>>( string path, MinMaxPair<float> def ) const {
    MinMaxPair<float> mmp;
    mmp.min = get<float>( path + ":min", def.min );
    mmp.max = get<float>( path + ":max", def.max );
    return mmp;
}

template <>
MinMaxPair<int> XmlConfig::get<MinMaxPair<int>>( string path ) const {
    MinMaxPair<int> mmp;
    mmp.min = get<int>( path + ":min", std::numeric_limits<int>::min() );
    mmp.max = get<int>( path + ":max", std::numeric_limits<int>::max() );
    return mmp;
}
template <>
MinMaxPair<int> XmlConfig::get<MinMaxPair<int>>( string path, MinMaxPair<int> def ) const {
    MinMaxPair<int> mmp;
    mmp.min = get<int>( path + ":min", def.min );
    mmp.max = get<int>( path + ":max", def.max );
    return mmp;
}

double ff( double *x, double *par ){
    return f1FFQ->Eval(x[0]) * pow( x[0], 1.2 ) ; 
}

class RhoGen : public HistoAnalyzer {
public:

    RhoGen() {}
    ~RhoGen() {}

    TH1 * hMass = nullptr;
    TH1 * hIncohPt = nullptr;
    TF1 * ff1 = nullptr;
    TF1 * fphoton = nullptr;
    TH1 * hvm = nullptr;
    TH1 * htoffEff = nullptr;
    double momres;
    double shift;
    bool apply_tof = true;
    
    MinMaxPair<float> pi_pt;
    MinMaxPair<float> pi_eta;
    
    float mu_max_p;

    double pt_range0_min, pt_range1_min;
    double pt_range0_max, pt_range1_max;

    MinMaxPair<float> rho_pt, rho_eta;


    virtual void initialize(){

        gRandom = new TRandom3();
        gRandom->SetSeed( 1331 );
        // f1FFQ = get<TF1>( "fFFQ", "FF" );
        // fphoton = get<TF1>( "fPhotonkt", "FF" );
        // ff1 = new TF1( "ff1", &ff, 0, 1, 0 );
        hvm = get<TH1>( "hvm", "FF" );
        hMass = get<TH1F>( "mc_mass", "mass" );
        hIncohPt = get<TH1F>( "mPt", "incoh" );

        htoffEff = get<TH1>( "tof_eff_pt", "TOFEFF" );

        LOG_F( INFO, "f1FFQ = %p", f1FFQ );
        LOG_F( INFO, "hMass = %p", hMass );

        book->cd();
        book->makeAll( config, nodePath + ".histograms" );

        momres = config.get<double>( "p.momres", 0.02 );

        apply_tof = config.get<bool>( "p.applyTof", true );

        
        mu_max_p = config.get<float>( "p.mu_max_p", 0.0 );
        LOG_F( INFO, "Reject muons with p below p=%0.3f", mu_max_p );

        pt_range0_min = config.get<float>( "p.pt_range0:min", 0 );
        pt_range0_max = config.get<float>( "p.pt_range0:max", 0.04 );

        pt_range1_min = config.get<float>( "p.pt_range1:min", 0.11 );
        pt_range1_max = config.get<float>( "p.pt_range1:max", 0.15 );

        rho_pt.set( 0, 100 );
        rho_pt = config.get<MinMaxPair<float>>( "p.rho.pt", rho_pt );
        rho_eta.set( -1, 1 );
        rho_eta = config.get<MinMaxPair<float>>( "p.rho.eta", rho_eta );

        // force the PDF to zero outside the range
        for ( int i = 1; i <= hvm->GetXaxis()->GetNbins(); i++ ){
            if ( rho_pt.in( hvm->GetBinCenter(i) ) == false )
                hvm->SetBinContent( i, 0 );
        }
        for ( int i = 1; i <= hIncohPt->GetXaxis()->GetNbins(); i++ ){
            if ( rho_pt.in( hIncohPt->GetBinCenter(i) ) == false )
                hIncohPt->SetBinContent( i, 0 );
        }

        pi_pt.set( 0.2, 100 );
        pi_pt = config.get<MinMaxPair<float>>( "p.pi.pt", pi_pt );
        pi_eta.set( -1, 1 );
        pi_eta = config.get<MinMaxPair<float>>( "p.pi.eta", pi_eta );

    }

    float tof_weight( float pt ){
        assert( htoffEff != nullptr && "TOF Efficiency is NULL" );
        float pts = pt;

        if ( pts > 0.7 ) pts = 0.7;
        float ev = htoffEff->Interpolate( pt );
        if ( ev <= 0.001 ) return 0;
        return ev;
    }

    void applyBoost( TLorentzVector &_parent_lv, TLorentzVector &_d_lv ){

        float betaX = _parent_lv.Px() / _parent_lv.E();
        float betaY = _parent_lv.Py() / _parent_lv.E();
        float betaZ = _parent_lv.Pz() / _parent_lv.E();

        _d_lv.Boost(betaX,betaY,betaZ);
    }

    /* Two Body Decay in the rest frame of the parent
     *
     * ONLY Valid for m1 == m2! The daughter mass must be equal for this simplified form.
     * The decay is computed and then the daughters are boosted into the frame of the parent
     */
    vector<TLorentzVector> twoBodyDecay( TLorentzVector _parent_lv, double M_d1, double M_d2 ){

        double dm2 = pow(M_d1 - M_d2,2);
        double sm2 = pow(M_d1 + M_d2,2);

        double M_p2 = pow(_parent_lv.M(), 2);

        double p = sqrt( (M_p2 - sm2) * (M_p2 - dm2) ) / ( 2. * _parent_lv.M() );

        if (fspin1 == nullptr){
            fspin1 = new TF1( "fspin1", "1 + pow( cos(x),2 )" );
            fspin1->SetRange( 0, TMath::Pi() );
            fspin1->SetNpx(1000);
        }


        double costheta = gRandom->Uniform(-1.,1.);
        // costheta = cos( fspin1->GetRandom() );
        double phi = gRandom->Uniform(0,TMath::Pi()*2);

        // make sure that the magnitude of the mom vector is correct!
        // May allow these distributions to be input?
        double pz = p*costheta;
        double px = p*sqrt(1.-costheta*costheta)*cos(phi);
        double py = p*sqrt(1.-costheta*costheta)*sin(phi);

        TLorentzVector daughter1( px, py, pz, sqrt( p*p + M_d1*M_d1 ));
        TLorentzVector daughter2( -px, -py, -pz, sqrt( p*p + M_d2*M_d2 ) );

        applyBoost(_parent_lv,daughter1);
        applyBoost(_parent_lv,daughter2);

        // LOG_F( INFO, "daughter1.M() %f == M %f", daughter1.M(), M_d1 );
        // LOG_F( INFO, "daughter2.M() %f == M %f", daughter2.M(), M_d2 );

        vector<TLorentzVector> dlv;
        dlv.push_back( daughter1 );
        dlv.push_back( daughter2 );
        return dlv;
    }

    // vector<TLorentzVector> twoBodyDecayInAcc( TLorentzVector _parent_lv, double M_d1, double M_d2 ) {
    //     vector<TLorentzVector> ds = twoBodyDecay( _parent_lv, M_d1, M_d2 );
    //     while( ds[0].Pt() < min_pt || ds[1].Pt() < min_pt ){
    //         ds = twoBodyDecay( _parent_lv, M_d1, M_d2 );
    //     }
    //     return ds;
    // }

    TLorentzVector genRho( bool coherent = true ){
        double pt = 0;

        if ( coherent ){
            pt = hvm->GetRandom();
        } else {
            pt = hIncohPt->GetRandom();
        }


        double mass = hMass->GetRandom();
        TLorentzVector lv;
        lv.SetPtEtaPhiM( 
            pt,
            gRandom->Uniform( rho_eta.min, rho_eta.max ),
            gRandom->Uniform(0,TMath::Pi()*2),
            mass
         );
        return lv;
    }

    bool in_star( float pt, float eta ){
        return ( pi_pt.in( pt ) && pi_eta.in( eta ) );
    }
    bool in_star( TLorentzVector &lv ){
        return in_star( lv.Pt(), lv.Eta() );
    }

    TProfile * makeProfile( string name, string profileName, int nRebin = 1 ){
        TH2 * h2D = book->get2D(name);
        // h2D->Sumw2();
        auto profileH = h2D->ProfileX(profileName.c_str());

        // this does the right thing, checked (maintains average)
        profileH->RebinX( nRebin );
        profileH->Scale(2.0);
        return  profileH;
    }

    TLorentzVector momSmear( TLorentzVector &lv ){
        double nPt = lv.Pt() * (1 + gRandom->Gaus( 0, momres ) );
        lv.SetPerp( nPt );
        return lv;
    }

    void fill_pol( TLorentzVector &lv, TLorentzVector &lvn,
        string prefix = "", float weight = 1.0 ){

        double pol = lvn.DeltaPhi( lv );
        if ( prefix == "DATA" ) {
            book->fill( "polptmass", lv.M(), lv.Pt(), pol, weight );
            book->fill( "c1ptmass", lv.M(), lv.Pt(), cos( 1 * pol), weight );
            book->fill( "c2ptmass", lv.M(), lv.Pt(), cos( 2 * pol), weight );
            book->fill( "c3ptmass", lv.M(), lv.Pt(), cos( 3 * pol), weight );
            book->fill( "c4ptmass", lv.M(), lv.Pt(), cos( 4 * pol), weight );
            book->fill( "c6ptmass", lv.M(), lv.Pt(), cos( 6 * pol), weight );
            return;
        }

        if ( lv.M() > 0.6 && lv.M() < 1.0 ){
            book->fill( prefix + "c1ppt", lv.Pt(), cos( pol ), weight );
            book->fill( prefix + "c1ppt", lv.Pt(), cos( pol ), weight );
            book->fill( prefix + "c2ppt", lv.Pt(), cos( 2*pol ), weight );
            book->fill( prefix + "c4ppt", lv.Pt(), cos( 4*pol ), weight );
        }

        if ( lv.Pt() > pt_range0_min && lv.Pt() < pt_range0_max  ){
            book->fill( prefix + "c1pmass0", lv.M(), cos( pol ), weight );
            book->fill( prefix + "c2pmass0", lv.M(), cos( 2*pol ), weight );
            book->fill( prefix + "c4pmass0", lv.M(), cos( 4*pol ), weight );
        }
        if ( lv.Pt() > pt_range1_min && lv.Pt() < pt_range1_max  ){
            book->fill( prefix + "c1pmass1", lv.M(), cos( pol ), weight );
            book->fill( prefix + "c2pmass1", lv.M(), cos( 2*pol ), weight );
            book->fill( prefix + "c4pmass1", lv.M(), cos( 4*pol ), weight );
        }
    }



    virtual void make(){
        LOG_SCOPE_FUNCTION( INFO );

        TLorentzVector lv, lvpipi, lvpimu, lvnpipi, lvnpimu;
        TLorentzVector lvi, lvipipi, lvinpipi, lvipimu, lvinpimu;

        vector<TLorentzVector> ds, dsi;
        vector<TLorentzVector> dsmu, dsimu;
        size_t nevents = config.get<size_t>( "p.nevents", 10000 );
        double mufrac = config.get<double>( "p.mufraction", 0.025 );
        double winc = config.get<double>( "p.winc", 1.0 );

        double mass_min = config.get<double>( "p.mass:min", 0.0 );
        double mass_max = config.get<double>( "p.mass:max", 2.0 );

        double w0 = 1.0, w0i = 1.0;
        double wmu = 1.0;

        LOG_F( INFO, "Starting event loop" );
        for ( size_t ievent = 0; ievent < nevents; ievent++ ){
            if ( ievent % 10000 == 0 ) {
                cout << "." << std::flush;
            }

            w0 = 1.0;
            wmu = 1.0;
            w0i = 1.0;

            // GENERATE A COHERENT RHO
            lv = genRho();
            book->fill( "mass", lv.M() );
            book->fill( "pt", lv.Pt() );

            // GENERATE AN INCOHERENT RHO
            lvi = genRho( false );
            book->fill( "mass", lvi.M(), winc );
            book->fill( "pt", lvi.Pt(), winc );


            // clear the vectors
            ds.clear();
            dsi.clear();
            dsmu.clear();

            // generate two body decay of the RHO and also of one pion to muon
            ds = twoBodyDecay( lv, 0.139, 0.139 );
            dsmu = twoBodyDecay( ds[0], 0.105, 1e-14 );

            float tofw = tof_weight( ds[0].Pt() ) * tof_weight( ds[1].Pt() );
            if (false == apply_tof)
                tofw=1.0;

            book->fill( "tof_weight_pt", ds[0].Pt(), tof_weight( ds[0].Pt() ) );
            

            book->fill( "pi_pt", ds[0].Pt() );
            book->fill( "pi_eta", ds[0].Eta() );
            
            book->fill( "mu_pt", dsmu[0].Pt(), mufrac );
            book->fill( "mu_eta", dsmu[0].Eta(), mufrac );

            if ( in_star( ds[0] ) ){
                book->fill( "star_pi_pt", ds[0].Pt() );
                book->fill( "star_pi_eta", ds[0].Eta() );
            }
            if ( in_star( dsmu[0] ) ){
                book->fill( "star_mu_pt", dsmu[0].Pt() );
                book->fill( "star_mu_eta", dsmu[0].Eta() );
            }

            lvpipi = momSmear(ds[0]) + momSmear(ds[1]);
            lvnpipi = ds[0] - ds[1];

            if ( !in_star( ds[0] ) || !in_star( ds[1] ) )
                w0 = 0.0;
            w0 *= tofw;

            book->fill( "tof_pair_weight_pt", lvpipi.Pt(), tofw );

            book->fill( "pipi_mass", lvpipi.M(), w0 );
            book->fill( "pipi_pt", lvpipi.Pt(), w0 );
            book->fill( "coh_mass", lvpipi.M(), w0 );
            book->fill( "coh_pt", lvpipi.Pt(), w0 );

            fill_pol( lvpipi, lvnpipi, "", w0 );
            fill_pol( lvpipi, lvnpipi, "pipi_", w0 );
            fill_pol( lvpipi, lvnpipi, "DATA", w0 );

            if ( lvpipi.Pt() > 0.150 && lvpipi.Pt() < 0.200 ){
                book->fill( "phi", lvpipi.Phi(), w0 );
                book->fill( "phi1", ds[0].Phi(), w0 );
                book->fill( "phi2", ds[1].Phi(), w0 );
            }


            book->fill( "pi_inp_pt", ds[0].Pt(), w0 );
            book->fill( "pi_inp_eta", ds[0].Eta(), w0 );

            book->fill( "pi_inp_pt", ds[1].Pt(), w0 );
            book->fill( "pi_inp_eta", ds[1].Eta(), w0 );
            
            

            // pi mu
            lvpimu = momSmear(ds[1]) + momSmear(dsmu[0]);
            lvnpimu = ds[1] - dsmu[0];
            // if ( gRandom->Rndm() > 0.5 ){
            //     lvnpimu = dsmu[0] - ds[1];
            // }
            if ( !in_star( dsmu[0] ) || !in_star( ds[1] ) )
                wmu = 0.0;
            if ( dsmu[0].P() < mu_max_p ) 
                wmu = 0.0;
            float tofmuw = pow(0.72256, 2);
            if ( false == apply_tof )
                tofmuw = 1.0;
            wmu *= tofmuw;

            book->fill( "wpi_pt", ds[0].Pt(), w0 );
            book->fill( "wmu_pt", dsmu[0].Pt(), mufrac * wmu );

            book->fill( "mu_inp_pt", dsmu[0].Pt(), mufrac * wmu );
            book->fill( "mu_inp_eta", dsmu[0].Eta(), mufrac * wmu );

            book->fill( "pi_inp_pt", ds[1].Pt(), mufrac * wmu );
            book->fill( "pi_inp_eta", ds[1].Eta(), mufrac * wmu );

            book->fill( "pimu_mass", lvpimu.M(), mufrac * wmu );
            book->fill( "pimu_pt", lvpimu.Pt(), mufrac * wmu );

            fill_pol( lvpimu, lvnpimu, "", mufrac * wmu );
            fill_pol( lvpimu, lvnpimu, "pimu_", mufrac * wmu );
            fill_pol( lvpimu, lvnpimu, "DATA", mufrac * wmu );
            

            dsi = twoBodyDecay( lvi, 0.139, 0.139 );
            dsimu = twoBodyDecay( ds[0], 0.105, 1e-14 );

            lvipipi = momSmear(dsi[0]) + momSmear(dsi[1]);
            lvinpipi = dsi[0] - dsi[1];

            // if (gRandom->Rndm() > 0.5){
            //     lvinpipi = dsi[1] - dsi[0];
            // }

            float tofiw = tof_weight( dsi[0].Pt() ) * tof_weight( dsi[1].Pt() );
            if ( false == apply_tof )
                tofiw = 1.0;

            if ( !in_star( dsi[0] ) || !in_star( dsi[1] ) )
                w0i = 0.0;
            w0i *= tofiw;

            book->fill( "pipi_mass", lvipipi.M(), winc * w0i );
            book->fill( "pipi_pt", lvipipi.Pt(), winc * w0i );

            fill_pol( lvipipi, lvinpipi, "", winc * w0i );
            fill_pol( lvipipi, lvinpipi, "inc_", winc * w0i );
            fill_pol( lvipipi, lvinpipi, "DATA", winc * w0i );

            book->fill( "inc_mass", lvipipi.M(), winc* w0i );
            book->fill( "inc_pt", lvipipi.Pt(), winc* w0i );


            lvipimu = momSmear(dsi[1]) + momSmear(dsimu[0]);
            lvinpimu = dsi[1] - dsimu[0];

            if ( !in_star( dsimu[0] ) || !in_star( dsi[1] ) )
                wmu = 0.0;
            if ( dsimu[0].P() < mu_max_p ) 
                wmu = 0.0;
            tofmuw = pow(0.72256, 2);
            tofmuw = 1.0;
            wmu *= tofmuw;


            fill_pol( lvipipi, lvinpipi, "", winc * mufrac * wmu );
            fill_pol( lvipipi, lvinpipi, "inc_", winc * mufrac * wmu );
            fill_pol( lvipipi, lvinpipi, "DATA", winc * mufrac * wmu );

            book->fill( "pimu_mass", lvpimu.M(), winc * mufrac * wmu );
            book->fill( "pimu_pt", lvpimu.Pt(), winc * mufrac * wmu );


        }  



        book->get( "pipi_mass" )->SetLineColor(kRed);
        book->get( "pipi_pt" )->SetLineColor(kRed);
        
        book->get( "inc_mass" )->SetLineColor(kBlack);
        book->get( "inc_pt" )->SetLineColor(kBlack);

        book->get( "pimu_mass" )->SetLineColor(kBlue);
        book->get( "pimu_pt" )->SetLineColor(kBlue);



        makeProfile( "c1ppt", "c1p" );
        makeProfile( "c2ppt", "c2p" );
        makeProfile( "c4ppt", "c4p" );

        makeProfile( "c1pmass0", "c1m0" );
        makeProfile( "c2pmass0", "c2m0" );
        makeProfile( "c4pmass0", "c4m0" );

        makeProfile( "c1pmass1", "c1m1" );
        makeProfile( "c2pmass1", "c2m1" );
        makeProfile( "c4pmass1", "c4m1" );

        makeProfile( "pipi_c1ppt", "pic1p" );
        makeProfile( "pipi_c2ppt", "pic2p" );
        makeProfile( "pipi_c4ppt", "pic4p" );

        makeProfile( "pimu_c1ppt", "muc1p" );
        makeProfile( "pimu_c2ppt", "muc2p" );
        makeProfile( "pimu_c4ppt", "muc4p" );

        makeProfile( "inc_c1ppt", "incc1p" );
        makeProfile( "inc_c2ppt", "incc2p" );
        makeProfile( "inc_c4ppt", "incc4p" );


        makeProfile( "pipi_c1pmass0", "pic1m0" );
        makeProfile( "pipi_c2pmass0", "pic2m0" );
        makeProfile( "pipi_c4pmass0", "pic4m0" );

        makeProfile( "pimu_c1pmass0", "muc1m0" );
        makeProfile( "pimu_c2pmass0", "muc2m0" );
        makeProfile( "pimu_c4pmass0", "muc4m0" );

        makeProfile( "inc_c1pmass0", "incc1m0" );
        makeProfile( "inc_c2pmass0", "incc2m0" );
        makeProfile( "inc_c4pmass0", "incc4m0" );

        makeProfile( "pipi_c1pmass1", "pic1m1" );
        makeProfile( "pipi_c2pmass1", "pic2m1" );
        makeProfile( "pipi_c4pmass1", "pic4m1" );

        makeProfile( "pimu_c1pmass1", "muc1m1" );
        makeProfile( "pimu_c2pmass1", "muc2m1" );
        makeProfile( "pimu_c4pmass1", "muc4m1" );

        makeProfile( "inc_c1pmass1", "incc1m1" );
        makeProfile( "inc_c2pmass1", "incc2m1" );
        makeProfile( "inc_c4pmass1", "incc4m1" );


        // write the config to the ROOT file
        TNamed n( "config", config.toXml() );
        n.Write();

    }
};


#endif