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
    c = 340.0           # Speed of sound [m/s].
    fs = 8000           # Sampling frequency [Hz].
    impsamples = 4000   # Number of samples for impulse respones [-].
    sinchalfwidth = 30  # Parameter for simulation precision [samples].
    responses = []      # Responses


    def __init__(self, dims, lspks, mics, beta):
        """
        Initialise room object.

        Parameters
        ----------
        dims : ndarray(ndim=1)
            Room dimension [m].
        lspks : ndarray(ndim=2)
            Loudspeaker positions [m].
        mics : ndarray(ndim=2)
            Microphone positions [m].
        beta : ndarray
            Reflection coefficients of walls as
            ndarray([left, right, front, back, floor, ceiling]).
            Should be contained in interval [-1, 1].
            Negative if reflections are to cause phase inversions.

        Returns
        -------
        -
        """
        self.dims = dims
        self.lspks = lspks
        self.mics = mics
        self.beta = beta


    def compute_responses(self, verbose=False):
        """
        Update responses from all loudspeakers to all microphones contained in
        the object.

        Parameters
        ----------
        verbose : bool
            Chooses wheter information is printed out along the way [-].

        Returns
        -------
        -
        """
        self.responses = np.empty((len(self.lspks), len(self.mics),
                                  self.impsamples), dtype=np.complex128)
        p = Pool(params.threads)
        X = []
        for k in xrange(len(self.lspks)):
            for l in xrange(len(self.mics)):
                X = X + [{'fs':self.fs, 'impsamples':self.impsamples,
                            'c':self.c, 'sinchalfwidth':self.sinchalfwidth,
                            'beta':self.beta, 'lspks':self.lspks[k],
                            'mics':self.mics[l], 'dims':self.dims}]

        mp = p.map_async(_f1, X)
        p.close()
        if verbose:
            while not mp.ready():
                rem = mp._number_left
                print "\r%d tasks left..." % (rem),
                sys.stdout.flush()
                time.sleep(0.5)
            print "\r",
        p.join()
        self.responses = np.array(mp.get()).reshape((len(self.lspks),
                                                          len(self.mics),
                                                          self.impsamples))


    def set_source_filters(self, filters):
        """
        Set filters in fron of loudspeakers.

        Parameters
        ----------
        filters : ndarray(ndim=2)
            Array of len(lspks) FIR filters. The array sould simply contain the
            impulse responses  [-].

        Returns
        -------
        -
        """
        self.filters = filters
        self.filtsamples = filters.shape[1]


    def get_resp_at_point(self, point):
        """
        Get the transfer function from all loudspeakers combined (with filters), 
        same filter input, to a point in space. Result is returned in frequency
        domain with the same grid as the room impulse responses are measured with.
        Simulation is caried out either in the frequency or time domain depending
        on the value of room.simtype.

        Parameters
        ----------
        point : ndarray(ndim=1)
            Point at which reponse should be simulated [m].

        Returns
        -------
        resp : ndarray(ndim=1)
            Complex frequency response at point.
        """
        resp = np.zeros((self.impsamples + self.filtsamples - 1,),
                        dtype=np.complex128)

        for k in xrange(len(self.lspks)):
            imp, _ = rimp.get_imp_response_approx(self.fs, self.impsamples,
                                                  self.c,
                                                  self.sinchalfwidth,
                                                  self.beta, self.lspks[k],
                                                  point, self.dims)

            filt = np.zeros((self.impsamples + self.filtsamples - 1,),
                            dtype=np.complex128)
            filt[0:self.filtsamples] = self.filters[k]
            resp = (resp +
                    np.fft.fft(imp, n=self.impsamples + self.filtsamples - 1)
                    * np.fft.fft(filt))

        return resp


    def get_resp_at_mics(self, n, length=0):
        """
        Get the transfer function at a selection of the contained microphones. This
        function relies on compute_responses() having been called after the room
        object has been initialised.

        Parameters
        ----------
        n : ndarray(ndim=1)
            Indices for microphones for which responses should be returned [-].
        length : positive integer
            Size of return frequency grid. Defaults to <impsamples>.
        Returns
        -------
        resp : ndarray(ndim=2)
            Complex frequency response at chosen michrophones [-].
        """
        if length == 0:
            length = self.impsamples

        resp = np.zeros((len(n), int(length)), dtype=np.complex128)
        for k in xrange(len(self.lspks)):
            resp = (resp
                 + np.fft.fft(np.fft.ifft(self.responses[k,n]), n=int(length))
                 * np.fft.fft(self.filters[k], n=int(length)))

        return resp

    def get_plane_resp_at_freq(self, f, height, maxdelay=0.0, resx=50, resy=50,
                               lims=np.zeros(4), verbose=False):

        """
        Get the response of a horizontal plane in the room to all the
        loudspeakers/filters contained at a single frequency. The function relies
        on the control filters having been set.

        Parameters
        ----------
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
        verbose : bool
            Chooses wheter information is printed out along the way [-].

        Returns
        -------
        grid : ndarray(ndim=2)
            Complex frequency response at chosen grid / frequency [-].
        """
        # Default limits
        if np.linalg.norm(lims) < 0.0001:
            lims = np.array([0.0, self.dims[0], 0.0, self.dims[1]])

        grid = np.zeros((resx, resy), dtype=np.complex128)

        p = Pool(params.threads)
        X = []
        for k in xrange(len(self.lspks)):
            X = X + [{'f':f, 'c':self.c, 'beta':self.beta,
                      'lspks':self.lspks[k], 'dims':self.dims, 'resx':resx,
                      'resy':resy, 'height':height, 'maxdelay':maxdelay,
                      'lims':lims}]

        mp = p.map_async(_f2, X)
        p.close()
        if verbose:
            while not mp.ready():
                rem = mp._number_left
                print "\r%d tasks left..." % (rem),
                sys.stdout.flush()
                time.sleep(0.5)
            print "\r",
        p.join()
        res = mp.get()

        for k in xrange(len(self.lspks)):
            exp = np.exp(-1J * 2 * np.pi * f / self.fs 
                        * np.arange(len(self.filters[k])))
            resp = np.inner(self.filters[k], exp)
            grid = grid + resp * res[k]

        return grid


    def convolve_at_point(self, sig, point):
        """
        Convolve a signal with the impulse response at a ceratin point in the
        room. Both filters and room response for all loudspeakers are taken into
        account.

        Parameters
        ----------
        sig : ndarray(ndim=1)
            Array with the samples of the signal [-].
        point : ndarray(ndim=1)
            Array with the position of the the point in the room [m].

        Returns
        -------
        csig : ndarray(ndim=1)
            Signal which has been convolved with room impulse at point [-].
        """

        imp = np.fft.ifft(self.get_resp_at_point(point)).real
        return np.convolve(sig, imp)
