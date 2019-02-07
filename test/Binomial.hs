
module Binomial (testBinomial) where

import Statistics.Matrix (fromList)
import Test.Tasty
import Test.Tasty.HUnit
import qualified Data.Vector.Unboxed as VU

import Statistics.Models.Linear.Generalized.Fit
import Statistics.Models.Linear.Generalized.Types

veq :: Double -> VU.Vector Double -> VU.Vector Double -> Bool
veq tol v0 v1 = VU.length v0 == VU.length v1
  && VU.and (VU.zipWith (\vi0 vi1 -> abs (vi0 - vi1) < tol) v0 v1)

random1 :: Assertion
random1 = assertBool "" $ veq 1e-4 ans b
  where
    ans = VU.fromList [-35.56101,  35.64225,  43.60247]
    b = glmFit fam x y
    dat = GLMData {
        glmY = y
      , glmX = x
      , glmWt = VU.replicate 4 1
      , glmFamily = fam
      }
    fam = binomialFamily linkLogit
    x = fromList 4 3 [
        1, -0.13085361, 0.2403199
      , 1,  1.25961852, 0.3239020
      , 1, -0.38175267, 1.6668202
      , 1, -0.01973108, 0.3110317
      ]
    y = VU.fromList [0, 1, 1, 0]

testBinomial :: TestTree
testBinomial = testGroup "binomial" [
    testCase "random 1" random1
  ]
