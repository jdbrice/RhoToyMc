//
// Created by James Brandenburg on 6/4/15.
//

#ifndef ROOBARB_XMLCANVAS_H
#define ROOBARB_XMLCANVAS_H

// STL
#include <map>
using namespace std;

#include "Logger.h"
#include "XmlConfig.h"
#include "XmlPad.h"

#include "TCanvas.h"


namespace jdb{
	class XmlCanvas : public IObject{

	public:
		string name, title;

		int pxWidth, pxHeight;
		int nCol, nRow;

		TCanvas * rootCanvas = nullptr;
		map<string, shared_ptr<XmlPad> > pads;

	
	virtual const char * classname() const { return "XmlCanvas"; }
		string getName() { return name; }

		XmlCanvas( XmlConfig &cfg, string _nodePath ) {
			if ( cfg.exists( _nodePath ) ) {
				// string preNode = cfg.cn(_nodePath);

				pxWidth 	= cfg.getInt( _nodePath + ":width", cfg.getInt( _nodePath + ":w", 800 ) );
				pxHeight	= cfg.getInt( _nodePath + ":height", cfg.getInt( _nodePath + ":h", 800 ) );

				name 		= cfg.getString( _nodePath + ":name", cfg.getString( _nodePath + ":n", "c" ) );
				title 		= cfg.getString( _nodePath + ":title", cfg.getString( _nodePath + ":t", "c" ) );

				nCol 		= cfg.getInt( _nodePath + ":columns", cfg.getInt( _nodePath + ":cols", cfg.getInt( _nodePath + ":nCols", 12 ) ) );
				nRow 		= cfg.getInt( _nodePath + ":rows", 12 );

				DEBUGC( "name=" << name << ", title=" << title <<", width=" << pxWidth << ", height=" << pxHeight );
				rootCanvas = new TCanvas( name.c_str(), title.c_str(), pxWidth, pxHeight );


				// createPads( cfg, _nodePath );

				// cfg.cn(preNode);
			} else {
				ERRORC("No XmlConfig Given : Creating Default Canvas" );
				rootCanvas = new TCanvas( "XmlCanvas", "XmlCanvas", 800, 1200 );
			}

		}

		void cd( int iPad = 0 ) {
			if ( !rootCanvas )
				return;
			rootCanvas->cd(iPad);
		}

		void saveImage( string iname ) {
			if ( !rootCanvas )
				return;

			rootCanvas->Print( iname.c_str() );
		}

		XmlPad* activatePad( string padName ){

			// cout << "ACTIVATE PAD " << padName << endl;
			if ( pads.find( padName ) != pads.end() && pads[ padName ] != nullptr ){
				// cout << "CD PAD " << padName << endl;
				pads[ padName ]->cd();
				// cout << "RETURN PAD " << endl;
				return pads[ padName ]->getPad();
			} else {
				ERRORC( "Could not find pad!" );
				rootCanvas->cd();
			}

			ERRORC( "Return nullptr" );
			return nullptr;
		}

		TCanvas * getCanvas(){
			return rootCanvas;
		}

		void createPad( XmlConfig &cfg, string _path ){
			this->cd();
			
			string name = cfg.getString( _path + ":name", cfg.getString( _path + ":n", "" ) );
			if ( "" != name ){
				cout << "MAKING THE XPAD (name=" << name << " )" << endl; 
				pads[ name ] = shared_ptr<XmlPad>( new XmlPad( cfg, _path, nCol, nRow ) );
			}
			else {
				cout << "CANNOT MAKE THE XPAD (name=" << name << " )" << endl; 
			}
		}

		

	protected:
		// void createPads( XmlConfig &cfg, string _nodePath ){
			
		// 	vector<string> children = cfg.childrenOf( _nodePath, "Pad" );
		// 	DEBUGC( "Found " << children.size() );
		// 	for ( string path : children ){
		// 		DEBUGC( "Creating Pad From " << path );
		// 		this->cd();
		// 		// cfg.cn( path );
		// 		string name = cfg.getString( path + ":name", cfg.getString( path + ":n", "" ) );
		// 		DEBUGC( "Creating Pad named " << name );
		// 		if ( "" != name ){
		// 			cout << "MAKING THE XPAD (name=" << name << ", path=" << path << " )" << endl; 
		// 			pads[ name ] = shared_ptr<XmlPad>( new XmlPad( cfg, path, nCol, nRow ) );
		// 		} else {
		// 			cout << "CANNOT MAKE THE XPAD (name=" << name << " )" << endl; 
		// 			ERRORC( "Cannot Create XPad" );
		// 		}
		// 	}

		// }

	};
}// namespace

#endif //ROOBARB_XMLCANVAS_H
