from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares


class ResFit():
    """
    Resonator fitter with cubic complex background and improved numerics & init.

    Improvements over previous version (addressing your poor-fit screenshot):
    - polynomial background is expressed in normalized frequency units x=(f-f0)/f0
      to avoid extremely large coefficients and numerical instability;
    - a simple automatic initial guess for Ql is computed from the measured
      linewidth (half-depth method) so curve_fit starts in a reasonable region;
    - bounds on Q's are enforced to avoid pathological values during optimization;
    - the cubic background initial estimate is computed on normalized x as well.

    Usage is the same as before. Call `fit_circle()` for IQ geometry, then
    `fit_resonance_complex()` to run the complex fit. Use `plot_iq()` and
    `plot_resonance()` to visualise.
    """

    def __init__(self, freqs, s21, id_min, id_max, ideal=False, fit_type="S21", poly_noise=True):
        self.raw_freqs = np.asarray(freqs)[id_min:id_max]
        self.raw_s21 = np.asarray(s21)[id_min:id_max]
        self.fit_type = fit_type
        self.ideal = ideal
        self.poly_noise = poly_noise

        # Detect whether s21 is complex (preferred) or real magnitude
        self.is_complex = np.iscomplexobj(self.raw_s21)

        # Store convenience versions
        if self.is_complex:
            self.s21_complex = self.raw_s21
            self.s21_mag = np.abs(self.s21_complex)
            self.s21_db = 20 * np.log10(self.s21_mag)
        else:
            # treated as linear magnitude (fallback)
            self.s21_complex = None
            self.s21_mag = np.asarray(self.raw_s21)
            self.s21_db = 20 * np.log10(np.maximum(self.s21_mag, 1e-20))

        # estimate fmin as frequency of minimum magnitude (typical resonance dip)
        self.fmin = self.raw_freqs[np.argmin(self.s21_mag)]
        # center frequency axis around fmin for diagnostics
        self.freqs_centered = self.raw_freqs - self.fmin

        # placeholders for fit results
        self.circle_center = None
        self.circle_radius = None
        self.bg_coeffs = None
        self.popt = None
        self.pcov = None

    # ----------------------- Circle fit utilities -----------------------
    def _algebraic_circle_fit(self, x, y):
        A = np.column_stack([x, y, np.ones_like(x)])
        b = x ** 2 + y ** 2
        params, *_ = np.linalg.lstsq(A, b, rcond=None)
        D, E, F = params
        xc = D / 2.0
        yc = E / 2.0
        r = np.sqrt(max(xc * xc + yc * yc + F, 0.0))
        return xc, yc, r

    def fit_circle(self, use_subset=0.6):
        if not self.is_complex:
            raise ValueError("Circle fit requires complex S21 data. Provide complex S21.")

        f = self.raw_freqs
        S = self.s21_complex

        mags = np.abs(S)
        mid = np.argmin(mags)
        n = len(f)
        half = max(3, int(n * use_subset / 2))
        i0 = max(0, mid - half)
        i1 = min(n, mid + half)

        x = np.real(S[i0:i1])
        y = np.imag(S[i0:i1])

        xc, yc, r = self._algebraic_circle_fit(x, y)
        self.circle_center = xc + 1j * yc
        self.circle_radius = r
        return self.circle_center, self.circle_radius

    # --------------------- Background estimation (cubic, normalized x) ------------------------
    def _fit_complex_cubic_bg_norm(self, f, S, f0_guess):
        """
        Fit a complex cubic polynomial in the normalized variable
        x = (f - f0_guess)/f0_guess. Returns complex coefficients a0..a3.
        """
        x = (f - f0_guess) / f0_guess
        X = np.vstack([np.ones_like(x), x, x ** 2, x ** 3]).T
        coeffs_real, *_ = np.linalg.lstsq(X, S.real, rcond=None)
        coeffs_imag, *_ = np.linalg.lstsq(X, S.imag, rcond=None)
        a0 = coeffs_real[0] + 1j * coeffs_imag[0]
        a1 = coeffs_real[1] + 1j * coeffs_imag[1]
        a2 = coeffs_real[2] + 1j * coeffs_imag[2]
        a3 = coeffs_real[3] + 1j * coeffs_imag[3]
        return a0, a1, a2, a3

    # ------------------- Simple initial Q estimate from linewidth --------------------
    def _estimate_Q_from_width(self, f, mags, f0_guess):
        """
        Estimate loaded Q (Ql) via a simple half-depth width measurement.
        This is robust for moderately clean data and gives a reasonable starting
        point for the optimizer.
        """
        n = len(f)
        idx_min = np.argmin(mags)
        # estimate baseline from edges (10% on either side)
        edge = max(1, int(0.1 * n))
        baseline = np.mean(np.concatenate([mags[:edge], mags[-edge:]]))
        min_val = mags[idx_min]
        depth = baseline - min_val
        if depth <= 0:
            # fallback
            return 1e4
        half_level = min_val + depth / 2.0
        # find left crossing
        left_idx = np.argmax(mags[:idx_min] > half_level) if idx_min > 0 else 0
        # find right crossing
        right_candidates = np.where(mags[idx_min + 1:] > half_level)[0]
        if right_candidates.size:
            right_idx = idx_min + 1 + right_candidates[0]
        else:
            right_idx = n - 1
        f_left = f[left_idx]
        f_right = f[right_idx]
        width = max(1e-12, f_right - f_left)
        Ql_est = abs(f0_guess / width) if width > 0 else 1e4
        # clamp to reasonable range
        Ql_est = float(np.clip(Ql_est, 10.0, 1e8))
        return Ql_est

    # ------------------- Complex resonator model fit (with cubic bg, normalized x) --------------------
    @staticmethod
    def _resonator_model_complex_cubic_norm(f, f0, Ql, Qc, phi,
                                           a0_r, a0_i, a1_r, a1_i, a2_r, a2_i, a3_r, a3_i):
        a0 = a0_r + 1j * a0_i
        a1 = a1_r + 1j * a1_i
        a2 = a2_r + 1j * a2_i
        a3 = a3_r + 1j * a3_i
        x = (f - f0) / f0
        P = a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3
        denom = 1.0 + 2j * Ql * (x)
        S_res = 1.0 - (Ql / Qc) * np.exp(1j * phi) / denom
        S = P * S_res
        return S

    def fit_resonance_complex(self, initial_params=None, bounds=None, maxfev=int(1e6)):
        if not self.is_complex:
            raise ValueError("Complex resonance fit requires complex S21 data. Provide complex S21.")

        f = self.raw_freqs
        S = self.s21_complex

        # initial f0 guess from dip
        f0_0 = self.fmin

        # estimate Ql from linewidth
        Ql_0 = self._estimate_Q_from_width(f, np.abs(S), f0_0)
        # set Qc guess somewhat larger than Ql (weakly coupled initial guess)
        Qc_0 = Ql_0 * 2.0

        # initial cubic background estimate in normalized units
        a0, a1, a2, a3 = self._fit_complex_cubic_bg_norm(f, S, f0_0)

        if initial_params is None:
            phi_0 = 0.0
            p0 = [f0_0, Ql_0, Qc_0, phi_0,
                  a0.real, a0.imag, a1.real, a1.imag, a2.real, a2.imag, a3.real, a3.imag]
        else:
            p0 = initial_params

        # sensible default bounds (protect against runaway coefficients)
        if bounds is None:
            lower = [min(f) * 0.9, 1.0, 1.0, -np.pi] + [-np.inf] * 8
            upper = [max(f) * 1.1, 1e9, 1e9, np.pi] + [np.inf] * 8
            bounds = (lower, upper)

        def fitfun_concat(f_in, f0, Ql, Qc, phi,
                          a0_r, a0_i, a1_r, a1_i, a2_r, a2_i, a3_r, a3_i):
            S_fit = ResFit._resonator_model_complex_cubic_norm(f_in, f0, Ql, Qc, phi,
                                                              a0_r, a0_i, a1_r, a1_i, a2_r, a2_i, a3_r, a3_i)
            return np.concatenate([np.real(S_fit), np.imag(S_fit)])

        ydata = np.concatenate([S.real, S.imag])

        popt, pcov = curve_fit(fitfun_concat, f, ydata, p0=p0, bounds=bounds, maxfev=maxfev)

        self.popt = popt
        self.pcov = pcov

        f0 = popt[0]
        Ql = popt[1]
        Qc = popt[2]
        phi = popt[3]
        try:
            Qi = 1.0 / (1.0 / Ql - 1.0 / Qc)
        except Exception:
            Qi = np.nan

        a0_fit = popt[4] + 1j * popt[5]
        a1_fit = popt[6] + 1j * popt[7]
        a2_fit = popt[8] + 1j * popt[9]
        a3_fit = popt[10] + 1j * popt[11]

        self.fit_results = {
            'f0': f0,
            'Ql': Ql,
            'Qc': Qc,
            'Qi': Qi,
            'phi': phi,
            'a0': a0_fit,
            'a1': a1_fit,
            'a2': a2_fit,
            'a3': a3_fit
        }

        print(f"Resonant frequency f0: {f0:.6e} Hz")
        print(f"Loaded Q (Ql): {Ql:.3e}")
        print(f"Coupling Q (Qc): {Qc:.3e}")
        print(f"Internal Q (Qi, derived): {Qi:.3e}")
        print(f"Complex phase phi: {phi:.4f} rad")

        return popt, pcov

    # ---------------------------- Plotting -----------------------------
    def plot_iq(self, show_fit=True, nsamples=300):
        if not self.is_complex:
            raise ValueError("IQ plot requires complex S21 data.")

        S = self.s21_complex
        plt.figure()
        plt.scatter(S.real, S.imag, s=8, label='Measured IQ')

        if self.circle_center is not None and self.circle_radius is not None:
            cc = self.circle_center
            r = self.circle_radius
            theta = np.linspace(0, 2 * np.pi, 400)
            circle = cc + r * (np.cos(theta) + 1j * np.sin(theta))
            plt.plot(circle.real, circle.imag, '-', label='Algebraic circle fit')

        if hasattr(self, 'fit_results'):
            fgrid = np.linspace(self.raw_freqs.min(), self.raw_freqs.max(), nsamples)
            p = self.popt
            Sfit = ResFit._resonator_model_complex_cubic_norm(fgrid, p[0], p[1], p[2], p[3],
                                                              p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11])
            plt.plot(Sfit.real, Sfit.imag, '-', label='Model fit')

        plt.xlabel('Re(S21)')
        plt.ylabel('Im(S21)')
        plt.axis('equal')
        plt.legend()
        plt.title('IQ plane')
        plt.grid(True)

    def plot_resonance(self, params=None, show_dB=True):
        if not self.is_complex and self.s21_db is None:
            raise ValueError("No magnitude data available to plot.")

        f = self.raw_freqs
        plt.figure()
        plt.plot(f * 1e-9, self.s21_db, label='Measured')

        if params is None and hasattr(self, 'popt'):
            params = self.popt
        if params is not None:
            p = params
            Sfit = ResFit._resonator_model_complex_cubic_norm(f, p[0], p[1], p[2], p[3],
                                                              p[4], p[5], p[6], p[7], p[8], p[9], p[10], p[11])
            mag_db = 20 * np.log10(np.abs(Sfit))
            plt.plot(f * 1e-9, mag_db, label='Fit')

        plt.xlabel('Frequency [GHz]')
        plt.ylabel('S21 [dB]')
        plt.legend()
        plt.grid(True)

    from scipy.optimize import least_squares

    def robust_fit_resonance(self,
                            n_starts=6,
                            loss='soft_l1',
                            f_tol=1e-9,
                            max_nfev=20000,
                            outlier_sigma=4.0,
                            do_outlier_rejection=True):
        """
        Robust fit using least_squares with multiple starts and log-QL/QC parameterization.

        Returns: best_result dict with keys 'popt', 'cost', 'success', 'message', 'pcov_approx' (approx),
                and fills self.popt/self.pcov/self.fit_results similarly to existing API.
        """
        # ensure complex data
        if not self.is_complex:
            raise ValueError("robust_fit_resonance requires complex S21 arrays")

        f = np.asarray(self.raw_freqs).ravel()
        S = np.asarray(self.s21_complex).ravel()

        # initial f0 guess from dip
        f0_guess = self.fmin

        # heuristic Q estimate from linewidth
        def _estimate_Ql():
            # use half-depth width estimate (same as earlier helper)
            mags = np.abs(S)
            n = len(f)
            edge = max(1, int(0.08 * n))
            baseline = np.mean(np.concatenate([mags[:edge], mags[-edge:]]))
            min_val = np.min(mags)
            depth = baseline - min_val
            if depth <= 0:
                return 1e4
            half_level = min_val + depth / 2.0
            mid = np.argmin(mags)
            # find crossings by simple index search
            left_idx = np.max(np.where(mags[:mid] > half_level)[0]) if mid>0 and np.any(mags[:mid] > half_level) else 0
            right_candidates = np.where(mags[mid+1:] > half_level)[0]
            right_idx = mid + 1 + right_candidates[0] if right_candidates.size else n-1
            width = max(1e-12, f[right_idx] - f[left_idx])
            return float(np.clip(abs(f0_guess / width), 10.0, 1e8))

        Ql_guess0 = _estimate_Ql()
        Qc_guess0 = Ql_guess0 * 2.0
        phi_guess0 = 0.0

        # estimate cubic polynomial in normalized x around f0_guess
        x_norm = (f - f0_guess) / f0_guess
        Xmat = np.vstack([np.ones_like(x_norm), x_norm, x_norm**2, x_norm**3]).T
        coef_real, *_ = np.linalg.lstsq(Xmat, S.real, rcond=None)
        coef_imag, *_ = np.linalg.lstsq(Xmat, S.imag, rcond=None)
        # a0..a3 initial
        a_init = np.array([coef_real[0], coef_imag[0],
                        coef_real[1], coef_imag[1],
                        coef_real[2], coef_imag[2],
                        coef_real[3], coef_imag[3]])

        # model builder: params = [f0, log_Ql, log_Qc, phi, a0_r, a0_i, a1_r, a1_i, a2_r, a2_i, a3_r, a3_i]
        def model_complex(params, freqs):
            f0 = params[0]
            Ql = np.exp(params[1])
            Qc = np.exp(params[2])
            phi = params[3]
            a0 = params[4] + 1j * params[5]
            a1 = params[6] + 1j * params[7]
            a2 = params[8] + 1j * params[9]
            a3 = params[10] + 1j * params[11]
            x = (freqs - f0) / f0
            P = a0 + a1 * x + a2 * x**2 + a3 * x**3
            denom = 1.0 + 2j * Ql * x
            S_res = 1.0 - (Ql / Qc) * np.exp(1j * phi) / denom
            return P * S_res

        # residual function (real+imag concatenated). We can weight by magnitude (optional)
        def residuals(params, freqs, data, weights=None):
            Sfit = model_complex(params, freqs)
            res = np.concatenate([ (Sfit.real - data.real), (Sfit.imag - data.imag) ])
            if weights is not None:
                return res * np.concatenate([weights, weights])
            return res

        # bounds: keep f0 inside measured window with a small margin, logs of Q bounded
        f_min, f_max = f.min(), f.max()
        f_margin = 0.02 * (f_max - f_min)
        lower = [f_min + f_margin, np.log(1.0), np.log(1.0), -np.pi] + [-np.inf]*8
        upper = [f_max - f_margin, np.log(1e9), np.log(1e9), np.pi] + [np.inf]*8

        # create multiple starting guesses: scale Ql by factors
        q_factors = np.geomspace(0.25, 4.0, n_starts)
        starts = []
        for qf in q_factors:
            p0 = [f0_guess,
                np.log(Ql_guess0 * qf),
                np.log(Qc_guess0 * qf),
                phi_guess0] + list(a_init)
            starts.append(np.array(p0, dtype=float))

        best = None
        best_cost = np.inf

        # optionally compute simple per-point weights: emphasize dip region
        mags = np.abs(S)
        # weight by inverse distance from dip index (more weight near dip)
        mid_idx = np.argmin(mags)
        idx = np.arange(len(self.raw_freqs))
        dist = np.abs(idx - mid_idx).astype(float) + 1.0
        weights = 1.0 / np.sqrt(dist)   # heuristic
        weights /= np.mean(weights)

        for p0 in starts:
            try:
                res = least_squares(residuals, p0, args=(f, S, weights), loss=loss,
                                    bounds=(lower, upper), max_nfev=max_nfev, xtol=1e-12, ftol=1e-12)
            except Exception as e:
                # try without weights as fallback
                res = least_squares(residuals, p0, args=(f, S, None), loss=loss,
                                    bounds=(lower, upper), max_nfev=max_nfev, xtol=1e-12, ftol=1e-12)
            if res.cost < best_cost and res.success:
                best = res
                best_cost = res.cost

        if best is None:
            # final fallback: run single start without weights with looser tolerances
            p0 = starts[0]
            best = least_squares(residuals, p0, args=(f, S, None), loss=loss, bounds=(lower, upper),
                                max_nfev=max_nfev)

        # optional outlier rejection & re-fit
        if do_outlier_rejection:
            rvec = residuals(best.x, f, S)
            # compute per-sample complex residual norm (sqrt(re^2+im^2))
            re = rvec[:len(f)]
            im = rvec[len(f):]
            res_norm = np.sqrt(re**2 + im**2)
            sigma = np.std(res_norm)
            mask = res_norm <= outlier_sigma * sigma
            if np.sum(mask) < len(mask):
                # re-fit using only masked points (weights set to None for simplicity)
                try:
                    f_mask = f[mask]
                    S_mask = S[mask]
                    # re-init with best.x
                    best2 = least_squares(residuals, best.x, args=(f_mask, S_mask, None),
                                        loss=loss, bounds=(lower, upper), max_nfev=max_nfev)
                    # accept if improved
                    if best2.cost < best.cost:
                        best = best2
                except Exception:
                    pass

        # store results in same format as original fit_resonance_complex
        popt = best.x.copy()
        # convert back Ql/Qc
        popt_full = popt.copy()
        popt_full[1] = np.exp(popt[1])
        popt_full[2] = np.exp(popt[2])

        # fill self.popt similar to earlier code (but note different ordering)
        # We'll put popt_full with same order used elsewhere:
        self.popt = popt_full
        # approximate covariance from Jacobian: cov = inv(J^T J) * cost / dof
        try:
            J = best.jac
            JTJ = J.T.dot(J)
            pcov = np.linalg.inv(JTJ)
        except Exception:
            pcov = None
        self.pcov = pcov

        # fill convenience dict
        f0 = popt_full[0]
        Ql = popt_full[1]
        Qc = popt_full[2]
        phi = popt_full[3]
        a0 = popt_full[4] + 1j * popt_full[5]
        a1 = popt_full[6] + 1j * popt_full[7]
        a2 = popt_full[8] + 1j * popt_full[9]
        a3 = popt_full[10] + 1j * popt_full[11]
        try:
            Qi = 1.0 / (1.0 / Ql - 1.0 / Qc)
        except Exception:
            Qi = np.nan
        self.fit_results = dict(f0=f0, Ql=Ql, Qc=Qc, Qi=Qi, phi=phi, a0=a0, a1=a1, a2=a2, a3=a3)

        # return a summary
        return {'popt': popt_full, 'cost': best.cost, 'success': best.success,
                'message': best.message, 'pcov_approx': pcov}


# Replace / add the following methods in your ResFit class

    # ------------------- Phase unwrap and delay estimation --------------------
    def _estimate_delay_from_phase(self, f, S):
        """
        Estimate an initial time delay tau by unwrapping measured phase and
        fitting a line: phase(f) ~ slope * f + intercept. Then tau = -slope/(2*pi).
        Returns tau_est (float, seconds). If estimation fails returns 0.0.
        """
        try:
            ph = np.angle(S)
            ph_un = np.unwrap(ph)
            # linear fit ph_un = m * f + b
            A = np.vstack([f, np.ones_like(f)]).T
            m, b = np.linalg.lstsq(A, ph_un, rcond=None)[0]
            tau_est = -m / (2.0 * np.pi)
            # clamp to reasonable delay magnitude (e.g. +-1e-6 s) to avoid extremes
            tau_est = float(np.clip(tau_est, -1e-6, 1e-6))
            return tau_est
        except Exception:
            return 0.0

    # ------------------- Complex resonator model fit (with cubic bg, normalized x) --------------------
    @staticmethod
    def _resonator_model_complex_cubic_norm(f, f0, Ql, Qc, phi, tau,
                                           a0_r, a0_i, a1_r, a1_i, a2_r, a2_i, a3_r, a3_i):
        """
        Complex resonator model extended with time delay 'tau' (seconds).
        Model: S(f) = P(x) * S_res(f) * exp(-1j*2*pi*f*tau)
        where P(x) is cubic complex background in normalized x=(f-f0)/f0.
        """
        a0 = a0_r + 1j * a0_i
        a1 = a1_r + 1j * a1_i
        a2 = a2_r + 1j * a2_i
        a3 = a3_r + 1j * a3_i
        x = (f - f0) / f0
        P = a0 + a1 * x + a2 * x ** 2 + a3 * x ** 3
        denom = 1.0 + 2j * Ql * (x)
        S_res = 1.0 - (Ql / Qc) * np.exp(1j * phi) / denom
        delay_factor = np.exp(-1j * 2.0 * np.pi * f * tau)
        S = P * S_res * delay_factor
        return S

    def fit_resonance_complex(self, initial_params=None, bounds=None, maxfev=int(1e6)):
        if not self.is_complex:
            raise ValueError("Complex resonance fit requires complex S21 data. Provide complex S21.")

        f = self.raw_freqs
        S = self.s21_complex

        # initial f0 guess from dip
        f0_0 = self.fmin

        # estimate Ql from linewidth
        Ql_0 = self._estimate_Q_from_width(f, np.abs(S), f0_0)
        Qc_0 = Ql_0 * 2.0

        # initial cubic background estimate in normalized units
        a0, a1, a2, a3 = self._fit_complex_cubic_bg_norm(f, S, f0_0)

        # estimate initial tau from unwrapped phase
        tau_0 = self._estimate_delay_from_phase(f, S)

        if initial_params is None:
            phi_0 = 0.0
            # New param order: f0, Ql, Qc, phi, tau, a0_r, a0_i, a1_r, a1_i, a2_r, a2_i, a3_r, a3_i
            p0 = [f0_0, Ql_0, Qc_0, phi_0, tau_0,
                  a0.real, a0.imag, a1.real, a1.imag, a2.real, a2.imag, a3.real, a3.imag]
        else:
            p0 = initial_params

        # sensible default bounds (protect against runaway coefficients)
        if bounds is None:
            # lower/upper for tau added (e.g. +/- 1 microsecond default)
            lower = [min(f) * 0.9, 1.0, 1.0, -np.pi, -1e-6] + [-np.inf] * 8
            upper = [max(f) * 1.1, 1e9, 1e9, np.pi, 1e-6] + [np.inf] * 8
            bounds = (lower, upper)

        def fitfun_concat(f_in, f0, Ql, Qc, phi, tau,
                          a0_r, a0_i, a1_r, a1_i, a2_r, a2_i, a3_r, a3_i):
            S_fit = ResFit._resonator_model_complex_cubic_norm(f_in, f0, Ql, Qc, phi, tau,
                                                              a0_r, a0_i, a1_r, a1_i, a2_r, a2_i, a3_r, a3_i)
            return np.concatenate([np.real(S_fit), np.imag(S_fit)])

        ydata = np.concatenate([S.real, S.imag])

        popt, pcov = curve_fit(fitfun_concat, f, ydata, p0=p0, bounds=bounds, maxfev=maxfev)

        self.popt = popt
        self.pcov = pcov

        # Extract results (note new ordering)
        f0 = popt[0]
        Ql = popt[1]
        Qc = popt[2]
        phi = popt[3]
        tau = popt[4]
        try:
            Qi = 1.0 / (1.0 / Ql - 1.0 / Qc)
        except Exception:
            Qi = np.nan

        a0_fit = popt[5] + 1j * popt[6]
        a1_fit = popt[7] + 1j * popt[8]
        a2_fit = popt[9] + 1j * popt[10]
        a3_fit = popt[11] + 1j * popt[12]

        self.fit_results = {
            'f0': f0,
            'Ql': Ql,
            'Qc': Qc,
            'Qi': Qi,
            'phi': phi,
            'tau': tau,
            'a0': a0_fit,
            'a1': a1_fit,
            'a2': a2_fit,
            'a3': a3_fit
        }

        print(f"Resonant frequency f0: {f0:.6e} Hz")
        print(f"Loaded Q (Ql): {Ql:.3e}")
        print(f"Coupling Q (Qc): {Qc:.3e}")
        print(f"Internal Q (Qi, derived): {Qi:.3e}")
        print(f"Complex phase phi: {phi:.4f} rad")
        print(f"Estimated delay tau: {tau:.4e} s")

        return popt, pcov

    # ---------------------------- robust_fit_resonance (least_squares) ------------------------------
    def robust_fit_resonance(self,
                            n_starts=6,
                            loss='soft_l1',
                            f_tol=1e-9,
                            max_nfev=20000,
                            outlier_sigma=4.0,
                            do_outlier_rejection=True):
        """
        Robust fit using least_squares with multiple starts and log-QL/QC parameterization.
        Extended to include time delay 'tau' with auto initial guess via phase unwrapping.
        """
        if not self.is_complex:
            raise ValueError("robust_fit_resonance requires complex S21 arrays")

        f = np.asarray(self.raw_freqs).ravel()
        S = np.asarray(self.s21_complex).ravel()

        # initial f0 guess from dip
        f0_guess = self.fmin

        # heuristic Q estimate
        def _estimate_Ql():
            mags = np.abs(S)
            n = len(f)
            edge = max(1, int(0.08 * n))
            baseline = np.mean(np.concatenate([mags[:edge], mags[-edge:]]))
            min_val = np.min(mags)
            depth = baseline - min_val
            if depth <= 0:
                return 1e4
            half_level = min_val + depth / 2.0
            mid = np.argmin(mags)
            left_idx = np.max(np.where(mags[:mid] > half_level)[0]) if mid>0 and np.any(mags[:mid] > half_level) else 0
            right_candidates = np.where(mags[mid+1:] > half_level)[0]
            right_idx = mid + 1 + right_candidates[0] if right_candidates.size else n-1
            width = max(1e-12, f[right_idx] - f[left_idx])
            return float(np.clip(abs(f0_guess / width), 10.0, 1e8))

        Ql_guess0 = _estimate_Ql()
        Qc_guess0 = Ql_guess0 * 2.0
        phi_guess0 = 0.0
        tau_guess0 = self._estimate_delay_from_phase(f, S)

        # estimate cubic polynomial in normalized x around f0_guess
        x_norm = (f - f0_guess) / f0_guess
        Xmat = np.vstack([np.ones_like(x_norm), x_norm, x_norm**2, x_norm**3]).T
        coef_real, *_ = np.linalg.lstsq(Xmat, S.real, rcond=None)
        coef_imag, *_ = np.linalg.lstsq(Xmat, S.imag, rcond=None)
        a_init = np.array([coef_real[0], coef_imag[0],
                           coef_real[1], coef_imag[1],
                           coef_real[2], coef_imag[2],
                           coef_real[3], coef_imag[3]])

        # model builder: params = [f0, log_Ql, log_Qc, phi, tau, a0_r, a0_i, a1_r, a1_i, a2_r, a2_i, a3_r, a3_i]
        def model_complex(params, freqs):
            f0 = params[0]
            Ql = np.exp(params[1])
            Qc = np.exp(params[2])
            phi = params[3]
            tau = params[4]
            a0 = params[5] + 1j * params[6]
            a1 = params[7] + 1j * params[8]
            a2 = params[9] + 1j * params[10]
            a3 = params[11] + 1j * params[12]
            x = (freqs - f0) / f0
            P = a0 + a1 * x + a2 * x**2 + a3 * x**3
            denom = 1.0 + 2j * Ql * x
            S_res = 1.0 - (Ql / Qc) * np.exp(1j * phi) / denom
            delay_factor = np.exp(-1j * 2.0 * np.pi * freqs * tau)
            return P * S_res * delay_factor

        def residuals(params, freqs, data, weights=None):
            Sfit = model_complex(params, freqs)
            res = np.concatenate([ (Sfit.real - data.real), (Sfit.imag - data.imag) ])
            if weights is not None:
                return res * np.concatenate([weights, weights])
            return res

        # bounds and starts
        f_min, f_max = f.min(), f.max()
        f_margin = 0.02 * (f_max - f_min)
        lower = [f_min + f_margin, np.log(1.0), np.log(1.0), -np.pi, -1e-6] + [-np.inf]*8
        upper = [f_max - f_margin, np.log(1e9), np.log(1e9), np.pi, 1e-6] + [np.inf]*8

        q_factors = np.geomspace(0.25, 4.0, n_starts)
        starts = []
        for qf in q_factors:
            p0 = [f0_guess,
                  np.log(Ql_guess0 * qf),
                  np.log(Qc_guess0 * qf),
                  phi_guess0,
                  tau_guess0] + list(a_init)
            starts.append(np.array(p0, dtype=float))

        best = None
        best_cost = np.inf

        # heuristic weights
        mags = np.abs(S)
        mid_idx = np.argmin(mags)
        idx = np.arange(len(self.raw_freqs))
        dist = np.abs(idx - mid_idx).astype(float) + 1.0
        weights = 1.0 / np.sqrt(dist)
        weights /= np.mean(weights)

        for p0 in starts:
            try:
                res = least_squares(residuals, p0, args=(f, S, weights), loss=loss,
                                    bounds=(lower, upper), max_nfev=max_nfev, xtol=1e-12, ftol=1e-12)
            except Exception:
                res = least_squares(residuals, p0, args=(f, S, None), loss=loss,
                                    bounds=(lower, upper), max_nfev=max_nfev, xtol=1e-12, ftol=1e-12)
            if res.cost < best_cost and res.success:
                best = res
                best_cost = res.cost

        if best is None:
            p0 = starts[0]
            best = least_squares(residuals, p0, args=(f, S, None), loss=loss, bounds=(lower, upper),
                                max_nfev=max_nfev)

        # outlier rejection & re-fit (unchanged, but model has tau now)
        if do_outlier_rejection:
            rvec = residuals(best.x, f, S)
            re = rvec[:len(f)]
            im = rvec[len(f):]
            res_norm = np.sqrt(re**2 + im**2)
            sigma = np.std(res_norm)
            mask = res_norm <= outlier_sigma * sigma
            if np.sum(mask) < len(mask):
                try:
                    f_mask = f[mask]
                    S_mask = S[mask]
                    best2 = least_squares(residuals, best.x, args=(f_mask, S_mask, None),
                                        loss=loss, bounds=(lower, upper), max_nfev=max_nfev)
                    if best2.cost < best.cost:
                        best = best2
                except Exception:
                    pass

        # final results: convert logs back for Q's and store
        popt = best.x.copy()
        popt_full = popt.copy()
        popt_full[1] = np.exp(popt[1])
        popt_full[2] = np.exp(popt[2])

        # store
        self.popt = popt_full
        try:
            J = best.jac
            JTJ = J.T.dot(J)
            pcov = np.linalg.inv(JTJ)
        except Exception:
            pcov = None
        self.pcov = pcov

        # convenience dict using new ordering (f0, Ql, Qc, phi, tau, a0,...)
        f0 = popt_full[0]
        Ql = popt_full[1]
        Qc = popt_full[2]
        phi = popt_full[3]
        tau = popt_full[4]
        a0 = popt_full[5] + 1j * popt_full[6]
        a1 = popt_full[7] + 1j * popt_full[8]
        a2 = popt_full[9] + 1j * popt_full[10]
        a3 = popt_full[11] + 1j * popt_full[12]
        try:
            Qi = 1.0 / (1.0 / Ql - 1.0 / Qc)
        except Exception:
            Qi = np.nan
        self.fit_results = dict(f0=f0, Ql=Ql, Qc=Qc, Qi=Qi, phi=phi, tau=tau, a0=a0, a1=a1, a2=a2, a3=a3)

        return {'popt': popt_full, 'cost': best.cost, 'success': best.success,
                'message': best.message, 'pcov_approx': pcov}
