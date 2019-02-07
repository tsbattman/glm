
module Statistics.Models.Linear.Generalized.Types (
    Family(..)
  , binomialFamily
  , Link(..)
  , linkLogit, linkProbit, linkCauchit, linkCloglog
  , linkIdentity, linkLog, linkSqrt, linkInvSq, linkInverse
  , isFinite, sq
  ) where

import Numeric (expm1)

import Statistics.Distribution (ContDistr(..), Distribution(..))
import Statistics.Distribution.CauchyLorentz (standardCauchy)
import Statistics.Distribution.Normal (standard)

data Family = Family {
    familyName :: String
  , familyLink :: Link
  , familyVariance :: Double -> Double
  , familyDevResid :: Double -> Double -> Double -> Double
  , familyValidMu :: Double -> Bool
  , familyInitialize :: Double -> Double -> Double
  }

ylogy :: Double -> Double -> Double
ylogy y mu = if y /= 0 then y * log (y / mu) else 0

binomialFamily :: Link -> Family
binomialFamily lnk = Family {
    familyName = "binomial"
  , familyLink = lnk
  , familyVariance = \mu -> mu * (1 - mu)
  , familyDevResid = \y mu wt -> 2 * wt * (ylogy y mu + ylogy (1 - y) (1 - mu))
  , familyValidMu = \mu -> isFinite mu && 0 < mu && mu < 1
  , familyInitialize = \y wt -> (y * wt + 0.5) / (wt + 1)
  }

data Link = Link {
    linkName :: String
  , linkFun :: Double -> Double
  , linkInv :: Double -> Double
  , linkDmuDeta :: Double -> Double
  , linkValidEta :: Double -> Bool
  }

linkLogit, linkProbit, linkCauchit, linkCloglog, linkIdentity, linkLog
  , linkSqrt, linkInvSq, linkInverse :: Link
linkLogit = Link {
    linkName = "logit"
  , linkFun = \mu -> log $ mu / (1 - mu)
  , linkInv = \eta -> recip $ 1 + exp (negate eta)
  , linkDmuDeta = \eta -> let e = exp eta in e / sq (1 + e)
  , linkValidEta = const True
  }
-- thresh <- -qnorm(.Machine$double.eps)
-- eta <- pmin(pmax(eta, -thresh), thresh)
linkProbit = Link {
    linkName = "probit"
  , linkFun = quantile standard
  , linkInv = cumulative standard
  , linkDmuDeta = density standard
  , linkValidEta = const True
  }
-- thresh <- -qcauchy(.Machine$double.eps)
-- eta <- pmin(pmax(eta, -thresh), thresh)
linkCauchit = Link {
    linkName = "cauchit"
  , linkFun = quantile standardCauchy
  , linkInv = cumulative standardCauchy
  , linkDmuDeta = density standardCauchy
  , linkValidEta = const True
  }
-- DmuDeta -> pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
linkCloglog = Link {
    linkName = "cloglog"
  , linkFun = \mu -> log . negate . log $ 1 - mu
  , linkInv = negate . expm1 . negate . exp
  , linkDmuDeta = \eta -> let e = min eta 700 in exp e * exp (negate (exp e))
  , linkValidEta = const True
  }
linkIdentity = Link {
    linkName = "identity"
  , linkFun = id
  , linkInv = id
  , linkDmuDeta = const 1
  , linkValidEta = const True
  }
linkLog = Link {
    linkName = "log"
  , linkFun = log
  , linkInv = exp
  , linkDmuDeta = exp
  , linkValidEta = const True
  }
linkSqrt = Link {
    linkName = "sqrt"
  , linkFun = sqrt
  , linkInv = \eta -> eta * eta
  , linkDmuDeta = (2 *)
  , linkValidEta = \eta -> isFinite eta && eta > 0
  }
linkInvSq = Link {
    linkName = "1/mu^2"
  , linkFun = \mu -> recip $ mu * mu
  , linkInv = \eta -> recip $ sqrt eta
  , linkDmuDeta = \eta -> negate . recip $ 2 * eta * sqrt eta
  , linkValidEta = \eta -> isFinite eta && eta > 0
  }
linkInverse = Link {
    linkName = "inverse"
  , linkFun = recip
  , linkInv = recip
  , linkDmuDeta = \eta -> negate . recip $ eta * eta
  , linkValidEta = \eta -> isFinite eta && eta /= 0
  }

isFinite :: Double -> Bool
isFinite a = not $ isInfinite a || isNaN a

sq :: Double -> Double
sq x = x * x
