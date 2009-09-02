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
#include <complex>
#include <fftw3.h>
#include <boost/random.hpp>

//#define STIM_SEED "seed"
//#define STIM_PER_CHAN "per_chan"

using namespace mw;

class mwMaskStimulus : public ImageStimulus{

protected:
    //shared_ptr<Variable> random_seed;
    shared_ptr<Variable> random_seed;
    //bool random_phase_per_channel;
    shared_ptr<Variable> random_phase_per_channel;
    //uint32_t random_seed;
    float *image_data; // of size 4 * height * width (for RGBA Float image format)
    //float *mask_data; // of size 4 * height * width (for RGBA Float image format)
    float *(channel_modulus[4]); // of size height * width * 3 ordered RGB
    //fftwf_complex (*channel_fft)[3];
    //fftwf_complex *fft_phase;
    //fftwf_complex *fft_mask;
    //Lockable* lock;
    boost::mt19937 rng;
    boost::uniform_real<> phase_distribution;
    boost::variate_generator<boost::mt19937, boost::uniform_real<> > random_phase_gen;
    //fftwf_plan fft_mask_plan[3];
    //fftwf_plan ifft_mask_plan[3];
    bool imageLoaded;
public:
	mwMaskStimulus(std::string _tag, std::string filename,
                                        shared_ptr<Variable> _xoffset,
                                        shared_ptr<Variable> _yoffset,
                                        shared_ptr<Variable> _xscale,
                                        shared_ptr<Variable> _yscale,
                                        shared_ptr<Variable> _rot,
                                        shared_ptr<Variable> _alpha,
                                        shared_ptr<Variable> _random_seed,
                                        shared_ptr<Variable> _random_phase_per_channel);
	mwMaskStimulus(const mwMaskStimulus &tocopy);
	~mwMaskStimulus();
    //shared_ptr<Variable> getRandomSeed();
    
    virtual void makeMask(StimulusDisplay *display);
    virtual void load(StimulusDisplay *display);
    //virtual bool isLoaded();
    virtual Data getCurrentAnnounceDrawData();
};

#endif 
