def get_csigma_cqr(cmatrix):
    # C11 = cmatrix[0, 0]
    C12 = cmatrix[0, 1] if cmatrix.ndim == 2 else cmatrix[:, 0, 1]
    C13 = cmatrix[0, 2] if cmatrix.ndim == 2 else cmatrix[:, 0, 2]
    C22 = cmatrix[1, 1] if cmatrix.ndim == 2 else cmatrix[:, 1, 1]
    C23 = cmatrix[1, 2] if cmatrix.ndim == 2 else cmatrix[:, 1, 2]
    C33 = cmatrix[2, 2] if cmatrix.ndim == 2 else cmatrix[:, 2, 2]

    # Formulas adapted from https://qudev.phys.ethz.ch/static/content/science/Documents/semester/Burkhard_Simon_SemesterThesis_130211.pdf
    c_sigma = ((C33 + C13)*(C22 + C12))/(C33 + C22 + C13 + C12) + C23
    beta = (C33*C12 - C22*C13)/((C33+C13)*(C22 + C12) + (C33 + C22 + C13 + C12)*C23)
    # Note that c_qr can also be written (perhaps more intuitively) as:
    # (C12*C33 - C13*C22) / (C22 + C12 + C33 + C13)
    c_qr = beta*c_sigma
    return c_sigma, c_qr