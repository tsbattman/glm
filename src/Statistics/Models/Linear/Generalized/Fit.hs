
module Statistics.Models.Linear.Generalized.Fit (
    GLMData(..)
  , equalWeightData
  , GLMControl(..)
  , dfltGLMControl
  , GLMIter(..)
  , glmIterate
  , glmFit
  , glm
  , glmEvalLink
  , glmEvalResponse
  ) where

import Numeric.Sum (kbn, sumVector)
import Numeric.LinearAlgebra ((#>), Matrix, Vector, asColumn, cmap, cols, flatten, konst, linearSolveSVD, rows, size)
import qualified Data.Vector.Storable as VS

import Statistics.Models.Linear.Generalized.Types

data GLMData = GLMData {
    glmY :: Vector Double
  , glmX :: Matrix Double
  , glmWt :: Vector Double
  , glmFamily :: Family
  }

equalWeightData :: Family -> Matrix Double -> Vector Double -> Maybe GLMData
equalWeightData fam x y
  | n == nr = Just $ GLMData y x (konst 1 n) fam
  | otherwise = Nothing
  where
    nr = rows x
    n = size y

data GLMControl = GLMControl {
    glmControlEpsilon :: {-# UNPACK #-}!Double
  , glmControlMaxIt :: {-# UNPACK #-}!Int
  } deriving (Eq, Show, Read)

dfltGLMControl :: GLMControl
dfltGLMControl = GLMControl 1e-8 25

data GLMIter = GLMIter {
    glmIterMu :: Vector Double
  , glmIterEta :: Vector Double
  , glmIterCoef :: Vector Double
  , glmIterDev :: {-# UNPACK #-} !Double
  } deriving (Eq, Show, Read)

glmIterBreak :: GLMControl -> GLMIter -> GLMIter -> Bool
glmIterBreak control g0 g1 = abs (dev1 - dev0) / (0.1 + abs dev1) < eps
  || VS.any (not . isFinite) (glmIterCoef g1)
  where
    dev1 = glmIterDev g1
    dev0 = glmIterDev g0
    eps = glmControlEpsilon control

devResiduals :: GLMData -> Vector Double -> Double
devResiduals (GLMData y _x wt fam) mu = sumVector kbn $ VS.zipWith3 devresid y mu wt
  where devresid = familyDevResid fam

glmStepValidate :: GLMData -> Vector Double -> Vector Double -> GLMIter
glmStepValidate dat@(GLMData _y x _wt fam) coef0 start
  | not (isFinite dev) = glmStepValidate dat coef0 midp
  | not (VS.all valideta eta1 && VS.all validmu mu1) = glmStepValidate dat coef0 midp
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
    eta1 = x #> start
    mu1 = cmap linvf eta1
    dev = devResiduals dat mu1
    midp = VS.zipWith (\si ci -> (si + ci) / 2) start coef0

glmIter :: GLMData -> GLMIter -> GLMIter
glmIter dat@(GLMData y x wt fam) (GLMIter mu0 eta0 coef0 _dev0) = glmStepValidate dat coef0
  . flatten
  $ linearSolveSVD (x * asColumn w) (asColumn (z * w))
  where
    varf = familyVariance fam
    lnk = familyLink fam
    dmudeta = linkDmuDeta lnk
    -- varmu = VU.map varf mu0
    mueta = cmap dmudeta eta0
    z = VS.zipWith4 (\yi mui etai muetai -> etai + (yi - mui) / muetai) y mu0 eta0 mueta
    w = VS.zipWith3 (\wi muetai mui -> sqrt $ (wi * sq muetai) / varf mui) wt mueta mu0

runLoop :: GLMControl -> [GLMIter] -> Vector Double
runLoop _ [] = VS.empty
runLoop _ [it0] = glmIterCoef it0
runLoop contr (it0:it1:ls) = if glmIterBreak contr it0 it1
  then glmIterCoef it1
  else runLoop contr (it1:ls)

glmIterate :: GLMData -> [GLMIter]
glmIterate dat = iterate (glmIter dat) (GLMIter mu eta (konst 0.0 m) dev)
  where
    m = cols (glmX dat)
    y = glmY dat
    wt = glmWt dat
    fam = glmFamily dat
    lnkf = linkFun $ familyLink fam
    linvf = linkInv $ familyLink fam
    mustart = VS.zipWith (familyInitialize fam) y wt
    eta = cmap lnkf mustart
    mu = cmap linvf eta
    dev = devResiduals dat mu

glmFit :: Family -> Matrix Double -> Vector Double -> Maybe (Vector Double)
glmFit fam x y = glm dfltGLMControl <$> equalWeightData fam x y

glm :: GLMControl -> GLMData -> Vector Double
glm contr = runLoop contr . take (glmControlMaxIt contr) . glmIterate

glmEvalLink :: Vector Double -> Vector Double -> Double
glmEvalLink a = sumVector kbn . VS.zipWith(*) a

glmEvalResponse :: Link -> Vector Double -> Vector Double -> Double
glmEvalResponse lnk x = linkInv lnk . glmEvalLink x
