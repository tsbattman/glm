{-# LANGUAGE ViewPatterns #-}

module Binomial (testBinomial) where

import Test.Tasty
import Test.Tasty.HUnit
import Numeric.LinearAlgebra ((#>), Vector, build, cmap, konst, matrix, maxElement, size, vector)
import qualified Data.Vector.Storable as VS

import Statistics.Models.Linear.Generalized.Fit
import Statistics.Models.Linear.Generalized.Types

veq :: Double -> Vector Double -> Vector Double -> Bool
veq tol v0 v1 = size v0 == size v1 && maxElement (abs (v0 - v1)) < tol

random1 :: Assertion
random1 = assertBool "" $ maybe False (veq 1e-4 ans) b
  where
    ans = vector [-35.56101,  35.64225,  43.60247]
    b = glmFit fam x y
    fam = familyBinomial linkLogit
    x = matrix 3 [
        1, -0.13085361, 0.2403199
      , 1,  1.25961852, 0.3239020
      , 1, -0.38175267, 1.6668202
      , 1, -0.01973108, 0.3110317
      ]
    y = vector [0, 1, 1, 0]

example1 :: Assertion
example1 = assertBool "" $ maybe False (veq 1e-4 ans) b
  where
    ldose = build 12 $ \i -> fromIntegral $ round i `mod` (6 :: Int)
    male = build 12 $ \i -> if i < 6 then 1.0 else 0.0
    x = build (12, 4) $ \(round -> i) (round -> j) -> case (j :: Int) of
      0 -> 1.0
      1 -> male VS.! i
      2 -> ldose VS.! i
      3 -> male VS.! i * ldose VS.! i
      _ -> error "bad dimension"
    numdead = vector [1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16]
    y = cmap (/ 20.0) numdead
    b = glmFit (familyBinomial linkLogit) x y
    ans = vector [-2.9935418, 0.1749868, 0.9060364, 0.3529130]

-- check against R solution, but check eta instead of coef
-- directly because we solve using a different method, svd vs qr
degen :: Assertion
degen = do
  print x
  assertBool "" $ veq 1e-4 (a #> ans) (a #> x)
  where
    a = matrix 2 [1, 2, 1, 2, 1, 2]
    b = vector [1, 0, 1]
    ans = vector [0.6931472, 0]
    x = glm dfltGLMControl GLMData {
        glmY = b
      , glmX = a
      , glmWt = konst 1 (size b)
      , glmFamily = familyBinomial linkLogit
      }

testBinomial :: TestTree
testBinomial = testGroup "binomial" [
    testCase "random 1" random1
  , testCase "example 1" example1
  , testCase "degenerate" degen
  ]
