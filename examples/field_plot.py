""" We consider a rectangular room with a single emitting source. We can
compute the response of any location in the room to this source. We compute a
single sample of the frequency response for a regular grid of points on a
horizontal plane spanning the entire width/length of the room (in some
specified distance from the floor). When plotted, this greatly visualizes the
response of the room to a source emitting a sinusoidal signal. This is plotted
for three different sets of reflection coefficients:

    1) No reflections (anechoic).
    2) Only one walls reflects.
    3) All walls reflect.

The amplitude and phase is plotted separately in each condition. The resulting
plot is saved to pdf. Note that this script can take several minutes to run!
I went a little crazy with matplotlib settings on this one. Sorry."""

from __future__ import print_function
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import mimpy


#-----------------------------------------------------------------------------#
#                               Parameters                                    #
#-----------------------------------------------------------------------------#
# Sampling frequency in Hz
fs = 20000

# Length of computed impulse responses in samples
impsamples = 5000

# Frequency at which the response is plotted in Hz
plot_freq = 1000

# Height of the plane which is plotted in meters
height = 2.0

# The dimensions of the room in meters
dims = np.array([4., 4., 4.])

# Reflection coefficients for the six walls
# Making the coefficients negative provides a more realistic response although
# this is not physically motivated
beta_1 = np.array([0., 0., 0., 0., 0., 0.])
beta_2 = np.array([-.7, 0., 0., 0., 0., 0.])
beta_3 = np.array([-.7, -.7, -.7, -.7, -.5, -.5])

# The position of the source
source = np.array([0.8, 2.0, 2.1])

#-----------------------------------------------------------------------------#
#                               Computations                                  #
#-----------------------------------------------------------------------------#
# Create room object
room = mimpy.room(dims, beta_1, fs, impsamples)

# Get impulse response from a to b
print('Computing responses 1/3 ...')
grid_1 = room.get_plane_resp_at_freq(source, plot_freq, height, 0.0, 100, 100)
print('Computing responses 2/3 ...')
room.beta = beta_2
grid_2 = room.get_plane_resp_at_freq(source, plot_freq, height, 0.0, 100, 100)
print('Computing responses 3/3 ...')
room.beta = beta_3
grid_3 = room.get_plane_resp_at_freq(source, plot_freq, height, 0.0, 100, 100)

#-----------------------------------------------------------------------------#
#                               Plotting                                      #
#-----------------------------------------------------------------------------#
# Matplotlib setup
try:
    mpl.rc('text', usetex=False)
    mpl.rc('font', family='sans-serif')
    mpl.rc('font', size=8)
except AttributeError:
    None

# Color for text
textcolor = (239./255, 239./255, 239./255)

# Generate colormaps for amplitude and phase
cmap_ampl = plt.get_cmap('YlOrRd')
cmap = plt.get_cmap('RdYlBu')
new_cmap_vals = ([cmap(i) for i in range(cmap.N)]
               + [cmap(i) for i in reversed(range(cmap.N))])
cmap_phase = mcolors.ListedColormap(new_cmap_vals)

# Make figure
fig = plt.figure(figsize=(20 / 2.54, 24 / 2.54), dpi=300)

# Plot condition 1 amplitude
ax1_a = fig.add_subplot(3, 2, 1)
plt.imshow(20*np.log10(np.abs(grid_1)), cmap=cmap_ampl, vmin=-40, vmax=0)
cbar1_a = plt.colorbar()
cbar1_a.set_label('Amplitude [dB]', color=textcolor)
plt.xlabel('x [m]', color=textcolor)
plt.ylabel('y [m]', color=textcolor)

# Plot condition 1 phase
ax1_p = fig.add_subplot(3, 2, 2)
plt.imshow(np.angle(grid_1), cmap=cmap_phase)
cbar1_p = plt.colorbar()
cbar1_p.set_label('Phase [rad]', color=textcolor)
plt.xlabel('x [m]', color=textcolor)
plt.ylabel('y [m]', color=textcolor)

# Plot condition 2 amplitude
ax2_a = fig.add_subplot(3, 2, 3)
plt.imshow(20*np.log10(np.abs(grid_2)), cmap=cmap_ampl, vmin=-40, vmax=0)
cbar2_a = plt.colorbar()
cbar2_a.set_label('Amplitude [dB]', color=textcolor)
plt.xlabel('x [m]', color=textcolor)
plt.ylabel('y [m]', color=textcolor)

# Plot condition 2 phase
ax2_p = fig.add_subplot(3, 2, 4)
plt.imshow(np.angle(grid_2), cmap=cmap_phase)
cbar2_p = plt.colorbar()
cbar2_p.set_label('Phase [rad]', color=textcolor)
plt.xlabel('x [m]', color=textcolor)
plt.ylabel('y [m]', color=textcolor)

# Plot condition 3 amplitude
ax3_a = fig.add_subplot(3, 2, 5)
plt.imshow(20*np.log10(np.abs(grid_3)), cmap=cmap_ampl, vmin=-40, vmax=0)
cbar3_a = plt.colorbar()
cbar3_a.set_label('Amplitude [dB]', color=textcolor)
plt.xlabel('x [m]', color=textcolor)
plt.ylabel('y [m]', color=textcolor)

# Plot condition 3 phase
ax3_p = fig.add_subplot(3, 2, 6)
plt.imshow(np.angle(grid_3), cmap=cmap_phase)
cbar3_p = plt.colorbar()
cbar3_p.set_label('Phase [rad]', color=textcolor)
plt.xlabel('x [m]', color=textcolor)
plt.ylabel('y [m]', color=textcolor)

# Change the colors of ... everything
for ax in [ax1_a, ax2_a, ax3_a, ax1_p, ax2_p, ax3_p]:
    ax.spines['bottom'].set_color(textcolor)
    ax.spines['top'].set_color(textcolor)
    ax.spines['left'].set_color(textcolor)
    ax.spines['right'].set_color(textcolor)
    ax.tick_params(axis='x', colors=textcolor)
    ax.tick_params(axis='y', colors=textcolor)

for cbar in [cbar1_a, cbar2_a, cbar3_a, cbar1_p, cbar2_p, cbar3_p]:
    cbar.ax.tick_params(axis='y', colors=textcolor)
    cbar.ax.spines['left'].set_color(textcolor)

plt.savefig('field_plot.pdf', facecolor=(43./255,58./255,66./255))
