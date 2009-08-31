/*
 *  mwMaskStimulus.cpp
 *  MWPlugin
 *
 *  Created by Brett Graham on 8/27/09.
 *  Copyright 2009 Harvard University. All rights reserved.
 *
 */

#include "mwMaskStimulus.h"


mwMaskStimulus::mwMaskStimulus(std::string _tag, std::string _filename,
                               shared_ptr<Variable> _xoffset,
                               shared_ptr<Variable> _yoffset,
                               shared_ptr<Variable> _xscale,
                               shared_ptr<Variable> _yscale,
                               shared_ptr<Variable> _rot,
                               shared_ptr<Variable> _alpha,
                               shared_ptr<Variable> _random_seed)
                                : ImageStimulus (_tag, _filename, _xoffset, _yoffset, _xscale, _yscale, _rot, _alpha) {
    random_seed = _random_seed;
}

mwMaskStimulus::mwMaskStimulus(const mwMaskStimulus &tocopy)
    :ImageStimulus((ImageStimulus&) tocopy) {
        random_seed = tocopy.random_seed;
}

mwMaskStimulus::~mwMaskStimulus(){
    if(loaded){
        free(image_data);
        for (int i = 0; i < 4; i++) {
            free(channel_modulus[i]);
        }
        // !!! TODO !!! free fftwf data
    }
}

shared_ptr<Variable> mwMaskStimulus::getRandomSeed() {
    return random_seed;
}

//void printStats(float* things, int N) {
//    float max = things[0];
//    float min = things[0];
//    float mean = things[0];
//    for (int i = 1; i < N; i++) {
//        max = (things[i] > max)? things[i]: max;
//        min = (things[i] < min)? things[i]: min;
//        mean += things[i];
//    }
//    mean /= N;
//    mprintf("\tMax: %f\tMin: %f\tMean: %f",max,min,mean);
//}
//
//void printStats(fftwf_complex* things, int N) {
//    float rMax = things[0][0]; float iMax = things[0][1];
//    float rMin = things[0][0]; float iMin = things[0][1];
//    float rMean = things[0][0]; float iMean = things[0][1];
//    for (int i = 1; i < N; i++) {
//        rMax = (things[i][0] > rMax)? things[i][0]: rMax;
//        iMax = (things[i][1] > iMax)? things[i][1]: iMax;
//        rMin = (things[i][0] < rMin)? things[i][0]: rMin;
//        iMin = (things[i][1] < iMin)? things[i][1]: iMin;
//        rMean += things[i][0];
//        iMean += things[i][1];
//    }
//    rMean /= N; iMean /= N;
//    mprintf("\tReal:");
//    mprintf("\t\tMax: %f\tMin: %f\tMean: %f",rMax,rMin,rMean);
//    mprintf("\tImag:");
//    mprintf("\t\tMax: %f\tMin: %f\tMean: %f",iMax,iMin,iMean);
//}

void mwMaskStimulus::makeMask(StimulusDisplay *display) {
    if(!loaded){
        mprintf("Attempted to make mask when image was not loaded, no mask made");
        return;
    }
    
    // generate random phase
    boost::mt19937 rng;
    // !!! TODO !!! should I use an existing pi value?
    boost::uniform_real<> phase_distribution(-3.14,3.14);
    //boost::uniform_real<> phase_distribution(0.0,6.2831853071795862);
    // make generator, which when called, generates a random number between 0 and 1.0
    boost::variate_generator<boost::mt19937, boost::uniform_real<> > random_phase_gen(rng,phase_distribution);
    // allocate space for random phase !!! possibly move this to load !!!
    fftwf_complex *random_phase = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * height * width);
    // use the generator to fill the phase array
    for (int i = 0; i < (height * width); i++) {
        complex<float> temp_phase(0.0,random_phase_gen());
        temp_phase = exp(temp_phase);
        random_phase[i][0] = real(temp_phase);
        random_phase[i][1] = imag(temp_phase);
        //random_phase[i] = reinterpret_cast<fftwf_complex>(exp(temp_phase));
    }
    // make space for mask
    float *mask_data = (float*)calloc(height*width*4,sizeof(float));
    // combine with modulus and do inverse
    for (int i = 0; i < 3; i++) {
        
        // make space for fftw input !!! possibly move this out of the loop? !!!
        //fftwf_complex *in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * height * width);
        fftwf_complex *in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * ((height * width)/ 2 + 1));
        // fill input
        //for (int j = 0; j < (height * width); j++) {
        for (int j = 0; j < ((height * width)/2 + 1); j++) {
            complex<float> temp_modulus(channel_modulus[i][j],0.0);
            complex<float> temp_phase(random_phase[j][0],random_phase[j][1]);
            complex<float> temp_result;
            temp_result = temp_modulus * temp_phase;
            //in[j] = random_phase[j] * temp_modulus;
            in[j][0] = real(temp_result);
            in[j][1] = imag(temp_result);
        }
        // make space for fftw output !!! possibly move this out of the loop? !!!
        float *out = (float*) calloc(height*width,sizeof(float));
        // make fftw plan !!! possibly move this out of the loop? !!!
        // using FFTW_MEASURE overwrites whatever is in the input array
        fftwf_plan fft_mask_plan = fftwf_plan_dft_c2r_2d(height, width, in, out, FFTW_ESTIMATE);
        // exectue fft
        fftwf_execute(fft_mask_plan);
        // copy out to mask_data
        for (int j = 0; j < (height * width); j++) {
            //mask_data[(j*4)+i] = out[j];
            mask_data[(j*4)+i] = out[j]/(height * width);
        }
        
        fftwf_destroy_plan(fft_mask_plan);
        fftwf_free(in);
        free(out);
        //lock->unlock();
    }
    
    // normalize
    for (int i = 0; i < 3; i++) {
        float maskMax = mask_data[i];
        float maskMin = mask_data[i];
        float maskMean = 0.0;
        float imageMax = image_data[i];
        float imageMin = image_data[i];
        float imageMean = 0.0;
        for (int j = 0; j < (height * width); j++) {
            maskMean += mask_data[(j*4)+i];
            maskMax = (mask_data[(j*4)+i] > maskMax) ? mask_data[(j*4)+i] : maskMax;
            maskMin = (mask_data[(j*4)+i] < maskMin) ? mask_data[(j*4)+i] : maskMin;
            imageMean += image_data[(j*4)+i];
            imageMax = (image_data[(j*4)+i] > imageMax) ? image_data[(j*4)+i] : imageMax;
            imageMin = (image_data[(j*4)+i] < imageMin) ? image_data[(j*4)+i] : imageMin;
        }
        maskMean /= (height * width);
        imageMean /= (height * width);
        // normalize
        float scale = imageMean/(maskMean-maskMin);
        scale = (scale > ((1-imageMean)/(maskMax-maskMean)))? ((1-imageMean)/(maskMax-maskMean)) : scale;
        //mScale = min(mean(image)/(mean(mask)-mask.min()), (1-mean(image))/(mask.max()-mean(mask)))
        for (int j = 0; j < (height * width); j++) {
            mask_data[(j*4)+i] = scale * mask_data[(j*4)+i] + (imageMean - maskMean * scale);
        }
    }
    
    // copy over alpha
    for (int i = 0; i < (height * width); i++) {
        // -- copy alpha --
        mask_data[(i*4)+3] = channel_modulus[3][i];
    }
    
    // delete old textures
    // move 'masks' (original images right now) to gpu
    for(int i = 0; i < display->getNContexts(); i++){
		display->setCurrent(i);
        GLuint texture_map;
        texture_map = *(texture_maps->getElement(i));
        //texture_maps.getElement(i,&texture_map);
        // delete old texture
        glDeleteTextures(1,&texture_map);
        // replace with new texture
        glGenTextures(1,&texture_map);
        glBindTexture(GL_TEXTURE_2D, texture_map);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_FLOAT, mask_data);
        gluBuild2DMipmaps(GL_TEXTURE_2D,GL_RGBA,width,height,GL_RGBA,GL_FLOAT,mask_data);
		if(texture_map){
			mprintf("Image Mask loaded into texture_map %d", texture_map);
		}
        // !!! do I need to readd the texture to the list !!!
        texture_maps->addElement(i, texture_map);
        glBindTexture(GL_TEXTURE_2D, 0);
	}
    
    // replace with new textures
    fftwf_free(random_phase);
    free(mask_data);
}

void mwMaskStimulus::load(StimulusDisplay *display) {
    if(loaded){
        // if stimulus is already loaded, just generate a new mask
        makeMask(display);
		return;
	}
	
    // check if file exists
    if(filename == ""){
		throw SimpleException("Cannot load image (NULL filename).");
		return;
	}
    FILE *test = fopen(filename.c_str(),"r");
	if(!test){
		throw SimpleException("Image file does not exist", filename);
	}
	fclose(test);
    
    // start up Devil
    ilInit();
    ilutInit(); // !!! necessary?
    ilutRenderer(ILUT_OPENGL); // !!! necessary?
    ilutEnable(ILUT_OPENGL_CONV); // !!! necessary?
    
    //load image from FILE
    GLuint texture_map;
    ILuint il_image_name;
	ILenum il_error;
    ilGenImages(1,&il_image_name);
    ilBindImage(il_image_name);
    ilLoadImage(filename.c_str());
    ilConvertImage(IL_RGBA, IL_FLOAT);
    if((il_error = ilGetError()) != IL_NO_ERROR) {
        // TODO HANDLE ERROR
        merror(M_DISPLAY_MESSAGE_DOMAIN,
               "IL Image Library Error: %x", 
               il_error);
		
		throw SimpleException("Cannot load image", filename);
		//lock->unlock();
        return;
    }
    
    width = (int)ilGetInteger(IL_IMAGE_WIDTH);
    height = (int)ilGetInteger(IL_IMAGE_HEIGHT);
    image_data = (float*) calloc(4 * height * width, sizeof(float));
    
    // copy data to image_data
    // BJG: the image data gets loaded upside down (this is probably mac specific)
    ilOriginFunc(IL_ORIGIN_LOWER_LEFT);
    ilEnable(IL_ORIGIN_SET);
    ilCopyPixels(0, 0, 0, (ILuint)width, (ILuint)height, 1, IL_RGBA, IL_FLOAT, (ILvoid*) image_data);
    if ((il_error = ilGetError()) != IL_NO_ERROR) {
        throw SimpleException(
                              (boost::format("Cannot copy image data from ilImage") % il_error).str());
        //image_loaded = false;
        free(image_data);
        return;
    }
    
    // break up image data into R G B A
    for (int i = 0; i < 4; i++) {
        channel_modulus[i] = (float*) calloc(height * width, sizeof(float));
    }
    for (int i = 0; i < (height * width); i++) {
        channel_modulus[0][i] = image_data[(i*4)]; // R
        channel_modulus[1][i] = image_data[(i*4)+1]; // G
        channel_modulus[2][i] = image_data[(i*4)+2]; // B
        channel_modulus[3][i] = image_data[(i*4)+3];// A
    }
    
    mprintf("=!= Modulus =!= %s",tag.c_str());
    // store magnitude (not phase)
    for (int i = 0; i < 3; i++) {
        mprintf("  %i::",i);
        //printStats(channel_modulus[i],(height*width));
        // make space for fftw result
        fftwf_complex *channel_fft = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * ((height * width)/2 +1));
        //channel_fft[i] = (fftwf_complex *)fftwf_malloc(sizeof(fftwf_complex) * height * width);
        // plan fft for this channel
        // using FFTW_MEASURE overwrites whatever is in the input array
        fftwf_plan fft_mask_plan = fftwf_plan_dft_r2c_2d(height, width, channel_modulus[i], channel_fft, FFTW_ESTIMATE);
        // need to do some paddings
        
        // perform fft storing output in channel_fft[i]
        fftwf_execute(fft_mask_plan);
        mprintf("\t-post-fft-");
        //printStats(channel_fft,((height*width)/2 + 1));
        // calculate modulus of fft: sqrt(real ** 2 + imag ** 2)
        for (int j = 0; j < ((height * width)/2+1); j++) {
            channel_modulus[i][j] = sqrt(channel_fft[j][0] * channel_fft[j][0] + channel_fft[j][1] * channel_fft[j][1]);
        }
        // !!! TODO !!! do I need scaling? 1/sqrt(rows * cols)
        // clean up
        fftwf_destroy_plan(fft_mask_plan);
        fftwf_free(channel_fft);
        mprintf("\t-post-modulus-");
        //printStats(channel_modulus[i],((height*width)/2+1));
    }
    
    // move 'masks' (original images right now) to gpu
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    for(int i = 0; i < display->getNContexts(); i++){
		display->setCurrent(i);
        
        glGenTextures(1,&texture_map);
        glBindTexture(GL_TEXTURE_2D, texture_map);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_FLOAT, mask_data);
        gluBuild2DMipmaps(GL_TEXTURE_2D,GL_RGBA,width,height,GL_RGBA,GL_FLOAT,image_data);
		if(texture_map){
			mprintf("Image Mask loaded into texture_map %d", texture_map);
		}
        // !!! do I need to readd the texture to the list !!!
        texture_maps->addElement(i, texture_map);
        glBindTexture(GL_TEXTURE_2D, 0);
	}
    
    if((il_error = ilGetError()) != IL_NO_ERROR) {
        throw SimpleException(
                              (boost::format("Cannot bind image texture. IL Image Library Error") % il_error).str());
        free(image_data);
        return;
    }
	
	ilDeleteImage(il_image_name);
	if((il_error = ilGetError()) != IL_NO_ERROR) {
        throw SimpleException(M_DISPLAY_MESSAGE_DOMAIN,
                              (boost::format("Cannot delete image. IL Image Library Error: %s") % il_error).str());
        free(image_data);
        return;
    }
    
    loaded = true;
    
    // make mask
    makeMask(display);
}

Data mwMaskStimulus::getCurrentAnnounceDrawData() {
    
    //mprintf("getting announce DRAW data for image stimulus %s",tag );
    
    Data announceData(M_DICTIONARY, 9);
    announceData.addElement(STIM_NAME,tag);        // char
    announceData.addElement(STIM_ACTION,STIM_ACTION_DRAW);
    announceData.addElement(STIM_TYPE,STIM_TYPE_IMAGE);
    announceData.addElement(STIM_FILENAME,filename);  
    announceData.addElement(STIM_POSX,last_posx);  
    announceData.addElement(STIM_POSY,last_posy);  
    announceData.addElement(STIM_SIZEX,last_sizex);  
    announceData.addElement(STIM_SIZEY,last_sizey);  
    announceData.addElement(STIM_ROT,last_rot);
    //TODO    announceData.addElement(STIM_ALPHA,last_alpha);  
    
    return (announceData);
}