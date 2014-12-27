import numpy as np
import time
import sys
from multiprocessing import Pool
import room_impulse as rimp
import parameters as params


#-----------------------------------------------------------------------------#
#                               Callback functions                            #
#-----------------------------------------------------------------------------#
def _f1(x):
    imp, _ = (
        rimp.get_imp_response_approx(x['fs'], x['impsamples'],
                                     x['c'],
                                     x['sinchalfwidth'],
                                     x['beta'], x['lspks'],
                                     x['mics'], x['dims']))

    return np.fft.fft(imp)


def _f2(x):
    g, _ = rimp.get_grid(x['f'], x['c'], x['height'], x['maxdelay'], x['resx'],
                         x['resy'], x['beta'], x['lspks'], x['dims'],
                         x['lims'])
    return g

#-----------------------------------------------------------------------------#
#                               Room Class                                    #
#-----------------------------------------------------------------------------#
class room:
    """
    The room object contains important acoustical features of a rectangular room
    and is capable of carrying out various simulations based on mirror image
    modelling.
    """
    # Default parameters
    c = 344.0           # Speed of sound [m/s].
    #fs = 8000           # Sampling frequency [Hz].
    #impsamples = 4000   # Number of samples for impulse respones [-].
    sinchalfwidth = 30  # Parameter for simulation precision [samples]. Higher
                        # is better but more expensive.


    def __init__(self, dims, beta, fs=8000, impsamples=4000):
        """
        Initialise room object.

        Parameters
        ----------
        dims : ndarray(ndim=1)
            Room dimension [m].
        beta : ndarray
            Reflection coefficients of walls as
            ndarray([left, right, front, back, floor, ceiling]).
            Should be contained in interval [-1, 1].
            Negative if reflections are to cause phase inversions.
        fs : integer
            Sampling frequency [Hz].
        impsamples : integer
            Length of computed impulse responses [samples].

        Returns
        -------
        -
        """
        self.dims = dims
        self.beta = beta
        self.fs = 8000
        self.impsamples = 4000


    def get_rir(self, source, receiver):
        """
        Get the transfer function from one point in space to another. Result is
        returned in frequency domain with the same grid as the room impulse
        responses are measured with.

        Parameters
        ----------
        source : ndarray(ndim=1)
            Point from which the rir should be simulated [m].
        receiver : ndarray(ndim=1)
            Point to which the rir should be simulated [m].

        Returns
        -------
        rir : ndarray(ndim=1)
            Impulse response at point.
        """

        rir, _ = rimp.get_imp_response_approx(self.fs, self.impsamples, self.c,
                                              self.sinchalfwidth, self.beta, 
                                              source, receiver, self.dims)

        return rir


    def get_plane_resp_at_freq(self, source, f, height, maxdelay=0.0, resx=50,
                               resy=50, lims=np.zeros(4)):
        """
        Get one frequency sample of the response from one source to a grid of
        points on a horizontal plane. This can be used to visualize the sound
        field at a frequency.

        Parameters
        ----------
        source : ndarray(ndim=1)
            Point from which the sound should emanate [m].
        f : float
            The frequency at which the simulation should be caried out [Hz].
        height : float
            The height of the plane in which the simulation is caried out [m].
        maxdelay : float
            The maximum time delay of reflections account for [s].
        resx : int
            Resolution of grid, x-axis [-].
        resy : int
            Resolution of grid, y-axis [-].
        lims : array
            Array with the grid limits [xmin, xmax, ymin, ymax] [m].

        Returns
        -------
        grid : ndarray(ndim=2)
            Complex frequency response at chosen grid / frequency [-].
        """
        # Default limits
        if np.linalg.norm(lims) < 0.0001:
            lims = np.array([0.0, self.dims[0], 0.0, self.dims[1]])

        grid, _ = rimp.get_grid(f, self.c, height, maxdelay, resx, resy,
                                self.beta, source, self.dims, lims)

        return grid
