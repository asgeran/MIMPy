""" Simply generates and plots the impulse response from one point to another
in an arbitrary space. The plot is saved as .pdf in the current folder. May
take several seconds to run. """

import numpy as np
import matplotlib.pyplot as plt
import mimpy


#-----------------------------------------------------------------------------#
#                               Parameters                                    #
#-----------------------------------------------------------------------------#
# Sampling frequency
fs = 20000

# Impulse response length in samples
impsamples = 8000

# The dimensions of the room in meters
dims = np.array([5., 6., 3.])

# Reflection coefficients for the six walls
# Making the coefficients negative provides a more realistic response although
# this is not physically motivated
beta = np.array([-.7, -.7, -.7, -.7, -.5, -.5])

# The two points
a = np.array([1.0, 1.0, 1.0])
b = np.array([4.0, 5.0, 2.0])

#-----------------------------------------------------------------------------#
#                               Computations                                  #
#-----------------------------------------------------------------------------#
# Create room object
room = mimpy.room(dims, beta, fs, impsamples)

# Get impulse response from a to b
rir = room.get_rir(a, b)

#-----------------------------------------------------------------------------#
#                               Plotting                                      #
#-----------------------------------------------------------------------------#
t = 1. / fs * np.arange(impsamples)
plt.plot(t, rir)
plt.grid()
plt.xlabel('Time [s]')
plt.ylabel('Value [-]')
plt.title('MIMPy room impulse response')
plt.savefig('rir_plot.pdf')
plt.show()
