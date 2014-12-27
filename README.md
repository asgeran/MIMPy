MIMPy
=====

A framework for simple acoustical simulations in rectangular rooms.

### Functionality
The package computes room impulse responses from ideal loudspeakers to ideal microphones in rectangular rooms.

### Algorithms
The computations are based on classical mirror image modelling (MIM) as proposed in [1] and with the modification proposed in [2]. See also [3] if you are interested in the specifics of how the algorithm is implemented.

### TODO
* The API is to be restructured before the package is really of any use.
* A number of demonstrations should be added to the package.
* I hope to include functionality for simulating binaural RIRs through the use of a HRTF database.

### References
[1] J. B. Allen and D. A. Berkley, “Image Method for Efficiently Simulating Small-Room Acoustics,” The Journal of the Acoustical Society of America, vol. 65, no. 4, pp. 943–950, Apr. 1979.  
[2] P. M. Peterson, “Simulating the Response of Multiple Microphones to a Single Acoustic Source in a Reverberant Room,” The Journal of the Acoustical Society of America, vol. 80, no. 5, pp. 1527–1529, Nov. 1986.  
[3] E. A. Lehmann and A. M. Johansson, “Prediction of Energy Decay in Room Impulse Responses Simulated with an Image-Source Model,” The Journal of the Acoustical Society of America, vol. 124, no. 1, pp. 269–277, Jul. 2008.  
