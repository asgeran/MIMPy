MIMPy
=====

A framework for simple acoustical simulations in rectangular rooms.

### Functionality
The package computes room impulse responses from ideal loudspeakers to ideal microphones in rectangular rooms. MIMPy can also help graphically visualize the sound field in a rectangular room. Have a look at the examples which are supplied with the source code.

It is the hope that MIMPy will eventually also be able to compute binaural room impulse responses and aid in producing binauralized audio   .

![Example of sound field visualization](/images/field_plot.png)

### Algorithms
The computations are based on classical mirror image modelling (MIM) as proposed in [1] and with the modification proposed in [2]. See also [3] if you are interested in the specifics of how the algorithm is implemented.

### TODO
* Add some system for connecting MIMPy to the CIPIC HRTF database. E.g. env-variable.
* Add room.py functionality to load CIPIC HRTFs.
* Add room.py functionality to compute average HRTFs.
* Add C function to get HRTFs.
* Add C function to simulate binaural RIRs.
* Add C function to simulate RIRs for an interval of time delays (from x ms to y ms).
* Add Cython function to combine the two above to compute approximate BRIRs at reasonable computational expense.
* Add room.py interface to the above.
* Add room.py functionality to binauralize audio with a moving source.
* Add example script(s) where binaural reverb is added to audio.
* Make everything PEP8.
* Provide some propper means of building and installing.
* Get the project on the PyPI repository.

### References
[1] J. B. Allen, D. A. Berkley, “Image Method for Efficiently Simulating Small-Room Acoustics,” The Journal of the Acoustical Society of America, vol. 65, no. 4, pp. 943–950, Apr. 1979.  
[2] P. M. Peterson, “Simulating the Response of Multiple Microphones to a Single Acoustic Source in a Reverberant Room,” The Journal of the Acoustical Society of America, vol. 80, no. 5, pp. 1527–1529, Nov. 1986.  
[3] E. A. Lehmann, A. M. Johansson, “Prediction of Energy Decay in Room Impulse Responses Simulated with an Image-Source Model,” The Journal of the Acoustical Society of America, vol. 124, no. 1, pp. 269–277, Jul. 2008.  
[4] V. R. Algazi, R. O. Duda, D. M. Thompson and C. Avendano, “The CIPIC HRTF Database,” Proc. 2001 IEEE Workshop on Applications of Signal Processing to Audio and Electroacoustics, pp. 99-102, Mohonk Mountain House, New Paltz, NY, Oct. 21-24, 2001.
