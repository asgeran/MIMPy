#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "_room_impulse.h"


int _get_imp_response_approx(double fs, int length, double c, int sinchalfwidth,
                             double* beta, double* lspk, double* mic,
                             double* room, double* res, int* nrefl)
{
    // Variables
    int a=0, b=0, d=0, l=0, m=0, n=0, i=0, lr=0, fb=0, tb=0;  // For iteration
    int nits=0, lits=0, totalrefl=0, ii=0, midsample=0; 
    double dist=0, maxdist=0, dx=0, dy=0, dz=0, amplitude=0, tau=0;

    maxdist = length * 1.0 / fs * c;

    // Compute reponse
    // Outer loops (a,b,d) run over the 2^3 room inversions
    // Middle loops (lr, fb, tb) run over mirror directions
    // Inner loops (m,n,l) run over the individual mirror images
    for(a = 0; a <= 1; a++)for(b = 0; b <= 1; b++)for(d = 0; d <= 1; d++)
    {
        m = 0;
        for(lr = -1; lr <= 1; lr=lr+2)
        {
        n = 0;
        for(fb = -1; fb <= 1; fb=fb+2)
        {
        l = 0;
        for(tb = -1; tb <= 1; tb=tb+2)
        {
            if(lr==-1) m=0; else m=1;
            while(1)
            {
                nits = 0;
                if(fb==-1) n=0; else n=1;
                while(1)
                {
                    lits = 0;
                    if(tb==-1) l=0; else l=1;
                    while(1)
                    {
                        dx = (2 * a - 1) * lspk[0] + mic[0] - 2 * m * room[0];
                        dy = (2 * b - 1) * lspk[1] + mic[1] - 2 * n * room[1];
                        dz = (2 * d - 1) * lspk[2] + mic[2] - 2 * l * room[2];
                        dist = sqrt(dx * dx + dy * dy + dz * dz);
                        if(dist > maxdist) break;
                        amplitude = pow(beta[0], abs(m-a))
                                  * pow(beta[1], abs(m))
                                  * pow(beta[2], abs(n-b))
                                  * pow(beta[3], abs(n))
                                  * pow(beta[4], abs(l-d))
                                  * pow(beta[5], abs(l)) / (4 * M_PI * dist);

                        tau = dist / c;
                        midsample = tau * fs;

                        // Compute endpoints for impulse/sinc convolution
                        if(midsample - sinchalfwidth < 0)
                            i = 0;
                        else i = midsample - sinchalfwidth;
                        if(midsample + sinchalfwidth >= length) 
                            ii = length;
                        else ii = midsample + sinchalfwidth;

                        // Carry out impulse/sinc convolution and add to reult
                        for(; i<ii; i++)
                        {
                            if(fabs(tau*fs-i) > 0.0001)
                            {
                                res[i] += amplitude * sin(M_PI*(tau*fs-i))
                                   / (M_PI*(tau*fs-i)) * 0.5
                                   * (1.0 + cos(M_PI*(tau*fs-i) / sinchalfwidth));
                            }
                            else
                            {
                                res[i] += amplitude * 0.5 
                                * (1.0 + cos(M_PI*(tau*fs-i) / sinchalfwidth));
                            }
                        }

                        totalrefl++;
                        lits++;
                        l = l + tb;
                    }

                    l = 0;
                    if(lits == 0) break;
                    nits++;
                    n = n + fb;
                }

                if(nits == 0) break;
                n = 0;
                m = m + lr;
            }
        }
        }
        }
    }
    *nrefl = totalrefl;
    return 0;
}


int _get_grid_at_freq(double f, double c, double maxdelay, double* beta,
                      double* lspk, double height, double* room, double* lims,
                      int resx, int resy, double* gridx, double* gridy,
                      int* nrefl)
{
    // Variables
    int a=0, b=0, d=0, l=0, m=0, n=0, i=0, lr=0, fb=0, tb=0; // For iteration
    int  x=0, y=0; // For iteration
    int nits=0, lits=0, totalrefl=0, distexceed=0; 
    double dist=0, maxdist=0, dx=0, dy=0, dz=0, amplitude=0, tau=0;

    complex double* grid =
        (complex double*)malloc(resx * resy * sizeof(complex double));
    if(!grid) return -1;
    for(i=0; i<resx*resy; i++)
        grid[i] = 0;
    maxdist = maxdelay * c;

    // Compute grid points
    double* xpts = (double*)malloc(resx * sizeof(double));
    if(!xpts) return -1;
    double* ypts = (double*)malloc(resy * sizeof(double));
    if(!ypts) return -1;
    for(x=0; x<resx; x++)
    {
        xpts[x] = lims[0] + (lims[1] - lims[0]) * (x + 0.5) / resx;
    }
    for(y=0; y<resy; y++)
    {
        ypts[y] = lims[2] + (lims[3] - lims[2]) * (y + 0.5) / resy;
    }

    // Compute reponse
    // Outer loops (a,b,d) run over the 2^3 room inversions
    // Middle loops (lr, fb, tb) run over mirror directions
    // Inner loops (m,n,l) run over the individual mirror images
    // Innermost loops (x, y) run over dimensions of grid
    for(a = 0; a <= 1; a++)for(b = 0; b <= 1; b++)for(d = 0; d <= 1; d++)
    {
        m = 0;
        for(lr = -1; lr <= 1; lr=lr+2)
        {
        n = 0;
        for(fb = -1; fb <= 1; fb=fb+2)
        {
        l = 0;
        for(tb = -1; tb <= 1; tb=tb+2)
        {
            if(lr==-1) m=0; else m=1;
            while(1)
            {
                nits = 0;
                if(fb==-1) n=0; else n=1;
                while(1)
                {
                    lits = 0;
                    if(tb==-1) l=0; else l=1;
                    while(1)
                    {
                        for(y=0; y < resy; y++)for(x=0; x < resx; x++)
                        {
                            // Compute distance from loadspeaker to microphone
                            // through reflection
                            dx = (2 * a - 1) * lspk[0] + xpts[x]
                               - 2 * m * room[0];
                            dy = (2 * b - 1) * lspk[1] + ypts[y]
                               - 2 * n * room[1];
                            dz = (2 * d - 1) * lspk[2] + height - 2 * l * room[2];
                            dist = sqrt(dx * dx + dy * dy + dz * dz);

                            // Compute amplitude of received signal
                            if(dist > maxdist) distexceed = 1;
                            amplitude = pow(beta[0], abs(m-a))
                                        * pow(beta[1], abs(m))
                                        * pow(beta[2], abs(n-b))
                                        * pow(beta[3], abs(n))
                                        * pow(beta[4], abs(l-d))
                                        * pow(beta[5], abs(l)) / (4 * M_PI * dist);

                            // Compute time delay
                            tau = dist / c;

                            // Compute grid point
                            grid[resy * x + y] += amplitude
                                            * cexp(-I * 2 * M_PI * f * tau);
                        }
                        totalrefl++;
                        if(distexceed)
                        {
                            distexceed = 0;
                            break;
                        }
                        lits++;
                        l = l + tb;
                    }

                    l = 0;
                    if(lits == 0) break;
                    nits++;
                    n = n + fb;
                }

                if(nits == 0) break;
                n = 0;
                m = m + lr;
            }
        }
        }
        }
    }

    for(y=0; y < resy; y++)for(x=0; x < resx; x++)
    {
        gridx[x + resx * y] = creal(grid[x + resx * y]);
        gridy[x + resx * y] = cimag(grid[x + resx * y]);
    }
    free(grid);
    *nrefl = totalrefl;
    return 0;
}
