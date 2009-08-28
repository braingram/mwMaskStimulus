/*
 *  mwMaskStimulus.h
 * mwMaskStimulus
 *  MWPlugin
 *
 *  Created by Brett Graham on 8/27/09.
 *  Copyright 2009 Harvard University. All rights reserved.
 *
 */

#ifndef mwMaskStimulus_H_
#define mwMaskStimulus_H_

#include <MonkeyWorksCore/Plugin.h>
#include <MonkeyWorksCore/StandardStimuli.h>

using namespace mw;

class mwMaskStimulus : public ImageStimulus{

protected:
    shared_ptr<Variable> random_seed;
    float *image_data; // of size 4 * height * width (for RGBA Float image format)
public:
	mwMaskStimulus(std::string _tag, std::string filename,
                                        shared_ptr<Variable> _xoffset,
                                        shared_ptr<Variable> _yoffset,
                                        shared_ptr<Variable> _xscale,
                                        shared_ptr<Variable> _yscale,
                                        shared_ptr<Variable> _rot,
                                        shared_ptr<Variable> _alpha,
                                        shared_ptr<Variable> _random_seed);
	mwMaskStimulus(const mwMaskStimulus &tocopy);
	~mwMaskStimulus();
    shared_ptr<Variable> getRandomSeed();
    
    virtual void makeMask(StimulusDisplay *display);
    virtual void load(StimulusDisplay *display);
//    virtual void drawInUnitSquare(StimulusDisplay *display);
    virtual Data getCurrentAnnounceDrawData();
};

#endif 
