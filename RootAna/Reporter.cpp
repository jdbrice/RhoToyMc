
#include "Reporter.h"
#include "Utils.h"
#include "LoggerConfig.h"
#include "TStyle.h"

namespace jdb{

	int Reporter::instances = 0;

    

	Reporter::Reporter( string filename, int w, int h ) : mt(0.1f), mr(0.1f), mb(0.13f), ml(0.13f) {
		DEBUG( classname(), "(" <<  filename <<  ", width=" <<  w <<  ", height=" <<  h << ")" )
	
        this->canvas = nullptr;
		this->filename = filename;
        this->w = w;
        this->h = h;

		setup();

		DEBUG( classname(), " Instance #" << instances );
	}

	Reporter::Reporter( string filename, TCanvas * _canvas ) : mt(0.1f), mr(0.1f), mb(0.13f), ml(0.13f) {
		DEBUG( classname(), "(" <<  filename <<  ", canvas=" << _canvas << ")" )
	
        this->canvas = nullptr;
		this->filename = filename;
        this->w = -1;
        this->h = -1;
		this->canvas = _canvas;
		
        setup();

		DEBUG( classname(), " Instance #" << instances );
	}

	Reporter::Reporter( XmlConfig &config, string np, string prefix ) : mt(0.1f), mr(0.1f), mb(0.13f), ml(0.13f) {
		DEBUG( classname(), "( config, np=" << np << ", prefix=" << prefix << ")" )

        this->canvas = nullptr;
		this->config = config;
		this->nodePath = config.basePath( np );

		this->filename = prefix + config.get<string>( nodePath + ".output:url",
            prefix + config.get<string>( nodePath + ":url", "rp.pdf" )
         );
        INFOC( "filename : " << this->filename << endl );

		this->w = config.getInt( nodePath + ".output:width", config.get<int>( nodePath + ":w", 400 ) );
		this->h = config.getInt( nodePath + ".output:height", config.get<int>( nodePath + ":h", 400 ) );

        INFOC( "w=" << w << ", h= " << h);

        vector<float> margins = config.getVector<float>( nodePath + ".Reporter:margins" );
        if ( margins.size() >= 4 ){
            this->mt = margins[0];
            this->mr = margins[1];
            this->mb = margins[2];
            this->ml = margins[3];
            gStyle->SetStatX( 0.99 );
            gStyle->SetStatW( this->mr/2.0 );
            gStyle->SetStatY( 1 - this->mt );
        }

		setup();

		DEBUG( classname(), " Instance #" << instances )
	}

    void Reporter::setup(){
        // unified setup entry point

        if ( nullptr == canvas ){
            canvas = new TCanvas( ("Reporter"+ts(instances)).c_str() , "canvas", this->w, this->h);
            instances++;
        }

        canvas->Print( ( filename + "[" ).c_str() );
        INFO( classname(), " Opening report " << filename );
        instances++;

        isOpen = true;

        margins();

    }

	Reporter::~Reporter() {
		// properly close the report file
		DEBUG( classname(), "" )
		
		if ( isOpen )
			canvas->Print( ( filename + "]" ).c_str() );
		delete canvas;

		DEBUG( classname(), filename << " Closed " );
	}

	void Reporter::close(){
		canvas->Print( ( filename + "]" ).c_str() );
		isOpen = false;
	}

	void Reporter::newPage( int dx, int dy, float marginX, float marginY ){
		DEBUG( classname(), " New Page ( " << dx << ", " << dy << ", " << marginX << ", " << marginY << " ) " );

		this->dx = dx;
		this->dy = dy;

		canvas->Clear();
		canvas->Divide( dx, dy, marginX, marginY );
		cd( 1 );

		// make the global margins set
		margins( );
	}

	void Reporter::cd( int pad ){
		DEBUG( classname(), " Current Pad " << pad  );
		canvas->cd( pad );
		currentPad = pad;
	}

	void Reporter::cd( int x, int y){
		cd( ( (y - 1 ) * dx) + x );
	}

	void Reporter::next(){
		DEBUG( classname(), " Pushing to next pad " );
		currentPad++;
		if ( currentPad > dx*dy){
			savePage();
			newPage( dx, dy);
		} else {
			cd( currentPad );
		}
	}

	void Reporter::savePage( string name ) {
		DEBUG( classname(), "Save : " << name  );
		if ( "" == name ){
			canvas->Print( ( filename ).c_str() );
			DEBUG( classname(), " Saving Page to " << filename );
		}
		else {
			canvas->Print( name.c_str() );
			DEBUG( classname(), " Saving Page to " << name );
		}
	}

	void Reporter::savePage( TCanvas *_canvas ){
		if ( nullptr == _canvas ) return;
		_canvas->Print( ( filename ).c_str() );
	}

	void Reporter::saveImage( string name ){
		DEBUG( classname(), " Saving Image to " << name );
		canvas->SaveAs( name.c_str() );
	}
	
}
