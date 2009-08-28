/*
 *  mwMaskStimulusPlugin.h
 *  MWPlugin
 *
 *  Created by Brett Graham on 8/27/09.
 *  Copyright 2009 Harvard University. All rights reserved.
 *
 */
#ifndef mwMaskStimulus_PLUGIN_H_
#define mwMaskStimulus_PLUGIN_H_

#include <MonkeyWorksCore/Plugin.h>
using namespace mw;

extern "C"{
    Plugin *getPlugin();
}

class mwMaskStimulusPlugin : public Plugin {
    
	virtual void registerComponents(shared_ptr<mw::ComponentRegistry> registry);	
};


#endif