
module Statistics.Models.Linear.Generalized.Fit (
    GLMData(..)
  , equalWeightData
  , GLMIter(..)
  , glmIterate
  , glmFit
  , glm
  ) where

import Numeric.Sum (kbn, sumVector)
import Statistics.Matrix (Matrix, Vector, dimension, generate, multiplyV, unsafeIndex)
import Statistics.Regression (ols)
import qualified Data.Vector.Unboxed as VU

import Statistics.Models.Linear.Generalized.Types

data GLMData = GLMData {
    glmY :: Vector
  , glmX :: Matrix
  , glmWt :: Vector
  , glmFamily :: Family
  }

equalWeightData :: Family -> Matrix -> Vector -> Maybe GLMData
equalWeightData fam x y
  | n == nr = Just $ GLMData y x (VU.replicate n 1) fam
  | otherwise = Nothing
  where
    (nr, _) = dimension x
    n = VU.length y

data GLMControl = GLMControl {
    glmControlEpsilon :: {-# UNPACK #-}!Double
  , glmControlMaxIt :: {-# UNPACK #-}!Int
  } deriving (Eq, Show, Read)

dfltGLMControl :: GLMControl
dfltGLMControl = GLMControl 1e-8 25

data GLMIter = GLMIter {
    glmIterMu :: Vector
  , glmIterEta :: Vector
  , glmIterCoef :: Vector
  , glmIterDev :: {-# UNPACK #-} !Double
  } deriving (Eq, Show, Read)

glmIterBreak :: GLMControl -> GLMIter -> GLMIter -> Bool
glmIterBreak control g0 g1 = abs (dev1 - dev0) / (0.1 + abs dev1) < eps
  where
    dev1 = glmIterDev g1
    dev0 = glmIterDev g0
    eps = glmControlEpsilon control

rowScale :: Matrix -> Vector -> Matrix
rowScale m0 v = generate nr nc $ \i j -> unsafeIndex m0 i j * v VU.! i
  where (nr, nc) = dimension m0

devResiduals :: GLMData -> Vector -> Double
devResiduals (GLMData y _x wt fam) mu = sumVector kbn $ VU.zipWith3 devresid y mu wt
  where devresid = familyDevResid fam

glmStepValidate :: GLMData -> Vector -> Vector -> GLMIter
glmStepValidate dat@(GLMData _y x _wt fam) coef0 start
  | not (isFinite dev) = glmStepValidate dat coef0 midp
  | not (VU.all valideta eta1 && VU.all validmu mu1) = glmStepValidate dat coef0 midp
  | otherwise = GLMIter {
      glmIterEta = eta1
    , glmIterMu = mu1
    , glmIterCoef = start
    , glmIterDev = dev
    }
  where
    lnk = familyLink fam
    validmu = familyValidMu fam
    valideta = linkValidEta lnk
    linvf = linkInv lnk
    eta1 = x `multiplyV` start
    mu1 = VU.map linvf eta1
    dev = devResiduals dat mu1
    midp = VU.zipWith (\si ci -> (si + ci) / 2) start coef0

glmIter :: GLMData -> GLMIter -> GLMIter
glmIter dat@(GLMData y x wt fam) (GLMIter mu0 eta0 coef0 _dev0) = glmStepValidate dat coef0 start
  where
    varf = familyVariance fam
    lnk = familyLink fam
    dmudeta = linkDmuDeta lnk
    -- varmu = VU.map varf mu0
    mueta = VU.map dmudeta eta0
    z = VU.zipWith4 (\yi mui etai muetai -> etai + (yi - mui) / muetai) y mu0 eta0 mueta
    w = VU.zipWith3 (\wi muetai mui -> sqrt $ (wi * sq muetai) / varf mui) wt mueta mu0
    start = ols (rowScale x w) (VU.zipWith (*) z w)

runLoop :: GLMControl -> [GLMIter] -> Vector
runLoop _ [] = VU.empty
runLoop _ (it0:[]) = glmIterCoef it0
runLoop contr (it0:it1:ls) = if glmIterBreak contr it0 it1
  then glmIterCoef it1
  else runLoop contr (it1:ls)

glmIterate :: GLMData -> [GLMIter]
glmIterate dat = iterate (glmIter dat) (GLMIter mu eta (VU.replicate n 0.0) dev)
  where
    n = VU.length y
    y = glmY dat
    wt = glmWt dat
    fam = glmFamily dat
    lnkf = linkFun $ familyLink fam
    linvf = linkInv $ familyLink fam
    mustart = VU.zipWith (familyInitialize fam) y wt
    eta = VU.map lnkf mustart
    mu = VU.map linvf eta
    dev = devResiduals dat mu

glmFit :: Family -> Matrix -> Vector -> Maybe Vector
glmFit fam x y = glm dfltGLMControl <$> equalWeightData fam x y

glm :: GLMControl -> GLMData -> Vector
glm contr = runLoop contr . take (glmControlMaxIt contr) . glmIterate
