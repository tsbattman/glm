
module Binomial (testBinomial) where

import Statistics.Matrix (fromList, fromRowLists, generate)
import Test.Tasty
import Test.Tasty.HUnit
import qualified Data.Vector.Unboxed as VU

import Statistics.Models.Linear.Generalized.Fit
import Statistics.Models.Linear.Generalized.Types

veq :: Double -> VU.Vector Double -> VU.Vector Double -> Bool
veq tol v0 v1 = VU.length v0 == VU.length v1
  && VU.and (VU.zipWith (\vi0 vi1 -> abs (vi0 - vi1) < tol) v0 v1)

random1 :: Assertion
random1 = assertBool "" $ maybe False (veq 1e-4 ans) b
  where
    ans = VU.fromList [-35.56101,  35.64225,  43.60247]
    b = glmFit fam x y
    fam = familyBinomial linkLogit
    x = fromList 4 3 [
        1, -0.13085361, 0.2403199
      , 1,  1.25961852, 0.3239020
      , 1, -0.38175267, 1.6668202
      , 1, -0.01973108, 0.3110317
      ]
    y = VU.fromList [0, 1, 1, 0]

example1 :: Assertion
example1 = assertBool "" $ maybe False (veq 1e-4 ans) b
  where
    ldose = VU.map fromIntegral $ VU.generate 12 (`mod` 6)
    male = VU.generate 12 $ \i -> if i < 6 then 1.0 else 0.0
    x = generate 12 4 $ \i j -> case j of
      0 -> 1.0
      1 -> male VU.! i
      2 -> ldose VU.! i
      3 -> male VU.! i * ldose VU.! i
      _ -> error "bad dimension"
    numdead = VU.fromList [1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16]
    y = VU.map (/ 20.0) numdead
    b = glmFit (familyBinomial linkLogit) x y
    ans = VU.fromList [-2.9935418, 0.1749868, 0.9060364, 0.3529130]

degen :: Assertion
degen = assertBool "" $ veq 1e-4 ans x
  where
    a = fromRowLists [[1, 2], [1, 2], [1, 2]]
    b = VU.fromList [1, 0, 1]
    ans = VU.fromList [0.6931472, 0]
    x = glm dfltGLMControl GLMData {
        glmY = b
      , glmX = a
      , glmWt = VU.replicate (VU.length b) 1
      , glmFamily = familyBinomial linkLogit
      }

testBinomial :: TestTree
testBinomial = testGroup "binomial" [
    testCase "random 1" random1
  , testCase "example 1" example1
  , testCase "degenerate" degen
  ]
