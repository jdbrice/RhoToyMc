

// RooBarb
#include "XmlConfig.h"
#include "TaskEngine.h"
using namespace jdb;

// STL
#include <iostream>
#include <exception>

#include "RhoGen.h"


#define LOGURU_IMPLEMENTATION 1
#include "loguru.h"

int main( int argc, char* argv[] ) {
	loguru::add_file("everything.log", loguru::Truncate, loguru::Verbosity_MAX);
	Logger::setGlobalLogLevel( "none" );

    TaskFactory::registerTaskRunner<RhoGen>( "RhoGen" );
    TaskEngine engine( argc, argv, "PairDstMaker" );
	
	return 0;
}
