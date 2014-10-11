int _get_imp_response_approx(double fs, int length, double c, int sincwidth,
                             double* beta, double* lspk, double* mic,
                             double* room, double* res, int* nrefl);

int _get_grid_at_freq(double f, double c, double maxdelay, double* beta,
                      double* lspk, double height, double* room, double* lims,
                      int resx, int resy, double* gridx, double* gridy,
                      int* nrefl);
