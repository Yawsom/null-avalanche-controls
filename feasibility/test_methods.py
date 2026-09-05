"""Independent numerical and behavioral checks; python -m unittest feasibility.test_methods."""
import unittest
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp
from bp_analyze import fit_mle
from .methods import (absolute_fit, block_fit, fit_counts, finite_refractory, lag_product,
                      measure, scalar_family, surrogate)


class MethodsTests(unittest.TestCase):
    def test_reject_bad_input(self):
        for times,dt in [(np.array([]),1),(np.array([1]),0),(np.array([1]),np.nan),
                         (np.array([-1]),1),(np.array([2,1]),1)]:
            with self.assertRaises(ValueError): measure(times,dt)

    def test_partition_and_float_bins(self):
        f=measure(np.array([0.,0.,1.,2.,5.,7.]),1.)
        self.assertEqual(list(f['size']),[4,1,1])
        self.assertEqual(int(measure(np.array([0.,1.,2.,4.]),1.5)['size'].sum()),4)

    def test_mle_independent_optimizer(self):
        counts=np.arange(30,0,-1)
        x=np.arange(1,31)
        theta,_=scalar_family(counts,np.log(x))
        oracle=minimize_scalar(lambda t: counts.sum()*logsumexp(-t*np.log(x))+t*(counts@np.log(x)),
                               bounds=(-5,12),method='bounded')
        self.assertAlmostEqual(theta[0],oracle.x,places=5)
        fit,_=fit_counts(counts)
        legacy=fit_mle(np.repeat(x,counts),1,30)
        self.assertAlmostEqual(fit['tau_mle'],legacy[0],places=6)
        self.assertAlmostEqual(fit['vuong_z'],legacy[2],places=5)

    def test_exact_powerlaw_and_alternative(self):
        rng=np.random.default_rng(1);x=np.arange(1,31)
        p=x**-1.5;p/=p.sum()
        fit=absolute_fit(rng.multinomial(10000,p),1,rng,99)
        self.assertLess(abs(fit['tau_mle']-1.5),.06)
        self.assertGreater(fit['ks_iid_p'],.01)
        self.assertLess(fit['aic_powerlaw_minus_lognormal'],8)
        p=np.exp(-.2*x);p/=p.sum()
        fit=absolute_fit(rng.multinomial(10000,p),1,rng,99)
        self.assertLess(fit['vuong_z'],-2)
        self.assertLessEqual(fit['ks_iid_p'],.02)

    def test_surrogate_preserves_window_and_count(self):
        times=np.array([0,1,19,20,39,99])
        moved=surrogate(times,np.random.default_rng(0),100,'interval',20)
        np.testing.assert_array_equal(times//20,moved//20)
        self.assertEqual(len(moved),len(times))

    def test_lag_statistic_against_dense(self):
        times=np.array([0,0,4,4,4,12,16,16])
        dense=np.bincount(times//4)
        self.assertEqual(lag_product(times),int(dense[:-1]@dense[1:]))

    def test_real_refractory_and_finite_span(self):
        ids,times,meta=finite_refractory(np.random.default_rng(3),.1,1.2)
        self.assertGreater(times.size,1000)
        self.assertLess(times.max(),360000)
        for electrode in np.unique(ids):
            self.assertTrue(np.all(np.diff(times[ids==electrode])>=20))
        self.assertEqual(meta['generation_limit_hits'],0)
        self.assertLess(meta['actual_offspring'],meta['attempted_offspring'])

    def test_block_resampling_identical_blocks(self):
        found=pd.DataFrame({'size':[1,2,3]*3,
                            'start_bin':[0,10,20,1000,1010,1020,2000,2010,2020]})
        fit=block_fit(found,1,3000,np.random.default_rng(3),high=3,reps=49)
        self.assertEqual(fit['block_count'],3)
        self.assertAlmostEqual(fit['block_tau_low'],0.,places=10)
        self.assertAlmostEqual(fit['block_tau_high'],0.,places=10)


if __name__=='__main__': unittest.main()
