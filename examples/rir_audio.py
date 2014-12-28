""" Computes the response from one point to another in a room. Load a .wav
file, convolves it with the response and saves the result.

The default settings of the script corresponds to a large reverberant space
like a church or similar. Try playing around with the parameters!

Beware that you have to supply a .wav file yourself! Also, this script may
take a minute or more to run! """

import numpy as np
import matplotlib.pyplot as plt
import mimpy
from scipy.io import wavfile
from scipy.signal import fftconvolve


#-----------------------------------------------------------------------------#
#                               Parameters                                    #
#-----------------------------------------------------------------------------#
# Audio file (not supplied!!)
inputname = 'audio.wav'
outputname = 'audio_reverberated.wav'

# Impulse response length in seconds
implen = 3.0

# The dimensions of the room in meters
dims = np.array([30., 40., 12.])

# Reflection coefficients for the six walls
# Making the coefficients negative provides a more realistic response although
# this is not physically motivated
beta = np.array([-.75, -.75, -.75, -.75, -.65, -.65])

# The two points (source and receiver)
a = np.array([1.0, 1.0, 2.0])
b = np.array([20.0, 25.0, 2.0])

#-----------------------------------------------------------------------------#
#                               Load audio                                    #
#-----------------------------------------------------------------------------#
fs, x = wavfile.read(inputname)

#-----------------------------------------------------------------------------#
#                               Computations                                  #
#-----------------------------------------------------------------------------#
# Create room object
room = mimpy.room(dims, beta, fs=fs, impsamples=int(implen*fs))

# Get impulse response from a to b
rir = room.get_rir(a, b)

#-----------------------------------------------------------------------------#
#                               Reverberate audio                             #
#-----------------------------------------------------------------------------#
# Carry out convolution separately for each channel
y = np.zeros((x.shape[0] + rir.shape[0] - 1, x.shape[1]))
for k in xrange(x.shape[1]):
    y[:,k] = fftconvolve(x[:,k], rir)

# Renormalize output
y = y * np.max(x) / np.max(y)

# Save to .wav file
wavfile.write(outputname, fs, y.astype(np.int16))
