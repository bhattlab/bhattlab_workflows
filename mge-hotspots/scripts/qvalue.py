# -*- python -*-
# -*- coding utf-8 -*-

#  This file is part of GDSCTools software
#
#  Copyright (c) 2015 - Wellcome Trust Sanger Institute
#  All rights reserved
#
#  File author(s): Thomas Cokelaer <cokelaer@gmail.comWE HERE>
#
#  Distributed under the BSD 3-Clause License.
#  See accompanying file LICENSE.txt distributed with this software
#
#  website: http://github.com/CancerRxGene/gdsctools
#
##############################################################################
"""Implementation of qvalue estimate

Author: Thomas Cokelaer
"""
# inspiration
# source: https://github.com/nfusi/qvalue/blob/master/qvalue/qvalue.py
# and qvalue (R package)
import numpy as np
import scipy.interpolate


class QValue(object):
    """Compute Q-value for a given set of P-values



    """
    def __init__(self, pv, lambdas=None,
        pi0=None, df=3, method='smoother', 
        smooth_log_pi0=False, verbose=True):
        """.. rubric:: Constructor

        The q-value of a test measures the proportion of false 
        positives incurred (called the false discovery rate or FDR) 
        when that particular test is called significant.

        :param pv: A vector of p-values (only necessary input)
        :param lambdas: The value of the tuning parameter to estimate pi_0.
            Must be in [0,1). Can be a single value or a range of values.
            If none, the default is a range from 0 to 0.9 with a step size of
            0.05 (inluding 0 and 0.9)
        :param method: Either "smoother" or "bootstrap"; the method for
            automatically choosing tuning parameter in the estimation of
            pi_0, the proportion of true null hypotheses. Only smoother 
            implemented for now.
        :param df: Number of degrees-of-freedom to use when
            estimating pi_0 with a smoother (default to 3 i.e., cubic 
            interpolation.)
        :param float pi0: if None, it's estimated as suggested in Storey and
            Tibshirani, 2003. May be provided, which is convenient for testing.
        :param smooth_log_pi0: If True and 'pi0_method' = "smoother",
            pi_0 will be estimated by applying a smoother to a
            scatterplot of log pi_0 rather than just pi_0

        .. note:: Estimation of pi0 differs slightly from the one given in R
            (about 0.3%) due to smoothing.spline function differences between
            R and SciPy.

        If no options are selected, then the method used to estimate pi_0
        is the smoother method described in Storey and Tibshirani (2003).
        The bootstrap method is described in Storey, Taylor & Siegmund
        (2004) but not implemented yet.


        .. seealso:: :class:`gdsctools.stats.MultipleTesting`
        """
        try:
            self.pv = np.array(pv)
        except:
            self.pv = pv.copy()
        assert(self.pv.min() >= 0 and self.pv.max() <= 1), \
            "p-values should be between 0 and 1"

        if lambdas is None:
            epsilon = 1e-8
            lambdas = scipy.arange(0,0.9+1e-8,0.05)

        if len(lambdas)>1 and len(lambdas)<4:
            raise ValueError("""if length of lambda greater than 1, you need at least 4 values""")

        if len(lambdas) >= 1 and (min(lambdas)<0 or max(lambdas)>=1):
            raise ValueError("lambdas must be in the range[0, 1[")
        self.m = float(len(self.pv))

        self.df = df 
        self.lambdas = lambdas
        self.method = method
        self.verbose = verbose
        self.smooth_log_pi0 = smooth_log_pi0
        self.pi0 = self.estimate_pi0(pi0)

    def estimate_pi0(self, pi0):
        """Estimate pi0 based on the pvalues"""
        pv = self.pv.ravel() # flatten array

        if pi0 is not None:
            pass
        elif len(self.lambdas) == 1:
            pi0 = np.mean(pv >= self.lambdas[0])/(1-self.lambdas[0])
            pi0 = min(pi0, 1)
        else:
            # evaluate pi0 for different lambdas
            pi0 = [np.mean(pv>=this)/(1-this) for this in self.lambdas]
            # in R
            # lambda = seq(0,0.09, 0.1)
            # pi0 = c(1.0000000, 0.9759067, 0.9674164, 0.9622673, 0.9573241,
            #         0.9573241 0.9558824, 0.9573241, 0.9544406, 0.9457901)
            # spi0 = smooth.spline(lambda, pi0, df=3, all.knots=F, spar=0)
            # predict(spi0, x=max(lambda))$y  --> 0.9457946
            # spi0 = smooth.spline(lambda, pi0, df=3, all.knots=F)
            # predict(spi0, x=max(lambda))$y  --> 0.9485383
            # In this function, using pi0 and lambdas, we get 0.9457946
            # this is not too bad, the difference on the v17 data set
            # is about 0.3 %
            if self.method == 'smoother':
                if (self.smooth_log_pi0):
                    pi0 = np.log(pi0)
                # In R, the interpolation is done with smooth.spline
                # within qvalue. However this is done with default
                # parameters, and this is different from the Python
                # code. Note, however, that smooth.spline has a parameter
                # called spar. If set to 0, then we would get the same
                # as in scipy. It looks like scipy has no equivalent of
                # the smooth.spline function in R if spar is not 0
                tck = scipy.interpolate.splrep(self.lambdas, pi0, 
                        k = self.df)
                pi0 = scipy.interpolate.splev(self.lambdas[-1], tck)
                if (self.smooth_log_pi0):
                    pi0 = np.exp(pi0)
                pi0 = min(pi0, 1.)
            elif self.method == 'bootstrap':
                raise NotImplementedError
                """minpi0 = min(pi0)
                mse = rep(0, len(lambdas))
                pi0.boot = rep(0, len(lambdas))
                for i in range(1,100):
                    p.boot = sample(p, size = m, replace = TRUE)
                    for i in range(0,len(lambdas)):
                        pi0.boot[i] <- mean(p.boot > lambdas[i])/(1 - lambdas[i])
                    mse = mse + (pi0.boot - minpi0)^2

                pi0 = min(pi0[mse == min(mse)])
                pi0 = min(pi0, 1)"""

            if pi0 > 1:
                if self.verbose:
                    print("got pi0 > 1 (%.3f) while estimating qvalues, setting it to 1" % pi0)

                pi0 = 1.0
        assert(pi0 >= 0 and pi0 <= 1), "pi0 is not between 0 and 1: %f" % pi0
        return pi0

    def qvalue(self):
        """Return the qvalues using pvalues stored in :attr:`pv` attribute"""
        pv = self.pv.ravel()
        p_ordered = scipy.argsort(pv)
        pv = pv[p_ordered]
        qv = self.pi0 * self.m/len(pv) * pv
        qv[-1] = min(qv[-1],1.0)

        for i in range(len(pv)-2, -1, -1):
            qv[i] = min(self.pi0*self.m*pv[i]/(i+1.0), qv[i+1])
        # reorder qvalues
        qv_temp = qv.copy()
        qv = scipy.zeros_like(qv)
        qv[p_ordered] = qv_temp

        # reshape qvalues
        original_shape = self.pv.shape
        qv = qv.reshape(original_shape)
        return qv

