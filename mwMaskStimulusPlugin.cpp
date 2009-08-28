/*
 *  mwMaskStimulusPlugin.cpp
 *  MWPlugin
 *
 *  Created by Brett Graham on 8/27/09.
 *  Copyright 2009 Harvard University. All rights reserved.
 *
 */

#include "mwMaskStimulusPlugin.h"
#include "mwMaskStimulusFactory.h"
#include "MonkeyWorksCore/ComponentFactory.h"
using namespace mw;

Plugin *getPlugin(){
    return new mwMaskStimulusPlugin();
}


void mwMaskStimulusPlugin::registerComponents(shared_ptr<mw::ComponentRegistry> registry) {
	
    // TODO: you need to customize the "signature" of the object your plugin will create
    //       The signature is of the form component/type –(e.g. stimulus/circle, or iodevice/NIDAQ)
    registry->registerFactory(std::string("stimulus/image_mask"),
							  (ComponentFactory *)(new mwMaskStimulusFactory()));
}

