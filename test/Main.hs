
module Main (main) where

import Test.Tasty

import Binomial

main :: IO ()
main = defaultMain $ testGroup "all tests" [
    testBinomial
  ]
