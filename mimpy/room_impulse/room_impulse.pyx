import numpy as np
cimport numpy as np

#-----------------------------------------------------------------------------#
#                           Import of C Functions                             #
#-----------------------------------------------------------------------------#
cdef extern from "_room_impulse.h": 
    int _get_imp_response_approx(double fs, int length, double c, int sincwidth,
                                 double* beta, double* lspk, double* mic,
                                 double* room, double* res, int* nrefl)

cdef extern from "_room_impulse.h": 
    int _get_grid_at_freq(double f, double c, double maxdelay, double* beta,
                          double* lspk, double height, double* room,
                          double * lims, int resx, int resy, double* gridx,
                          double* gridy, int* nrefl)

#-----------------------------------------------------------------------------#
#                       Python Interface to Functions                         #
#-----------------------------------------------------------------------------#
def get_imp_response_approx(fs, length, c, sinchalfwidth,
                            np.ndarray[np.double_t,ndim=1] beta,
                            np.ndarray[np.double_t,ndim=1] lspk,
                            np.ndarray[np.double_t,ndim=1] mic,
                            np.ndarray[np.double_t,ndim=1] room):
    """
    Get the simulated approximated discrete impulse response of an
    omnidirectional michrophone excited by an omnidirectional loudspeaker in a
    room. A sinc function multiplied by a hanning window is used to approximate
    the response of impulses which fall in between samples.

    Parameters
    ----------
    fs : float
        Sampling frequency [Hz].
    length : int
        Length of response [samples].
    c : float
        The speed of sound [m/s].
    sinchalfwidth : int
        Half the width of the sinc used to implement impulses between samples.
    beta : ndarray
        Reflection coefficients of walls as
        ndarray([left, right, front, back, floor, ceiling]).
        Should be contained in interval [-1, 1].
        Negative if reflections are to cause phase inversions.
    lspk: ndarray
        Loudspeaker position [x, y, z] [m].
    mic: ndarray
        Microphone position [x, y, z] [m].
    room: ndarray
        Room dimensions [x, y, z] [m].

    Returns
    -------
    impulse : ndarray
        Impulse response.
    nrefl : int
        Number of reflections included.

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import room_impulse
    >>> beta = np.array([-0.8, -0.8, -0.8, -0.8, -0.5, -0.7])
    >>> lspk = np.array([1.0, 1.5, 1.0])
    >>> mic = np.array([3.0, 4.0, 1.0])
    >>> room = np.array([5.0, 5.0, 3.0])
    >>> imp, n = room_impulse.get_imp_response_approx(44100.0, 10000, 340.0,
                                                      beta, lspk, mic, room)
    >>> print n     # Number of reflections
    >>> plt.plot(imp)
    >>> plt.show()
    """
    # Ensure that ndarrays are contiguous
    beta = np.ascontiguousarray(beta)
    lspk = np.ascontiguousarray(lspk)
    mic = np.ascontiguousarray(mic)
    room = np.ascontiguousarray(room)

    # Allocate space for results
    cdef np.ndarray[np.double_t, ndim=1, mode="c"] res = np.zeros(length)
    cdef int nrefl = 0

    # Compute response
    err = _get_imp_response_approx(fs, length, c, sinchalfwidth,
                                   <double*> beta.data, <double*> lspk.data,
                                   <double*> mic.data, <double*> room.data,
                                   <double*> res.data, &nrefl)
    return res, nrefl


def get_grid_at_freq(f, c, height, maxdelay, resx, resy,
                     np.ndarray[np.double_t,ndim=1] beta,
                     np.ndarray[np.double_t,ndim=1] lspk,
                     np.ndarray[np.double_t,ndim=1] room, 
                     np.ndarray[np.double_t,ndim=1] lims):
    """
    Get the simulated response to a monopole across a horizontal plane. 

    Parameters
    ----------
    f : float
        Frequency [Hz].
    c : float
        The speed of sound [m/s].
    height : float
        The height of the plane at which the simulation is performed.
    maxdelay : float
        Maximum time delay of reflections.
    resx : int
        Resolution, x axis.
    resy : int
        Resolution, y axis.
    beta : ndarray
        Reflection coefficients of walls as
        ndarray([left, right, front, back, floor, ceiling]).
        Should be contained in interval [-1, 1].
        Negative if reflections are to cause phase inversions.
    lspk: ndarray
        Loudspeaker position [x, y, z] [m].
    room: ndarray
        Room dimensions [x, y, z] [m].
    room: ndarray
        Grid limits [xmin, xmax, ymin, ymax] [m].

    Returns
    -------
    impulse : ndarray
        Impulse response.
    nrefl : int
        Number of reflections included.
    """
    # Ensure that ndarrays are contiguous
    beta = np.ascontiguousarray(beta)
    lspk = np.ascontiguousarray(lspk)
    room = np.ascontiguousarray(room)

    # Allocate space for results
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] gridx = np.zeros((resx, resy))
    cdef np.ndarray[np.double_t, ndim=2, mode="c"] gridy = np.zeros((resx, resy))
    cdef int nrefl = 0

    err = _get_grid_at_freq(f, c, maxdelay, <double*> beta.data,
                            <double*> lspk.data, height, <double*> room.data,
                            <double*> lims.data, resx, resy,
                            <double*> gridx.data, <double*> gridy.data, &nrefl)

    return np.copy(gridx + 1J * gridy), nrefl
