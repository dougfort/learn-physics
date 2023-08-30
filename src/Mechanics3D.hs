{-# OPTIONS -Wall #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Mechanics3D where

import qualified Graphics.Gloss as G
import Graphics.Gnuplot.Simple
  ( Aspect (..),
    Attribute (..),
    plotFunc,
    plotPaths,
  )
import Mechanics1D
  ( Diff (..),
    NumericalMethod,
    RealVectorSpace (..),
    Time,
    TimeStep,
    rungeKutta4,
    solver,
  )
import SimpleVec
  ( PosVec,
    R,
    Vec,
    iHat,
    jHat,
    kHat,
    magnitude,
    sumV,
    vec,
    xComp,
    yComp,
    zComp,
    zeroV,
    (*^),
    (<.>),
    (><),
    (^*),
    (^+^),
    (^-^),
    (^/),
  )
import SpatialMath
  ( Euler (..),
    V3 (..),
  )

-- import qualified Vis as V

data ParticleState = ParticleState
  { mass :: R,
    charge :: R,
    time :: R,
    posVec :: Vec,
    velocity :: Vec
  }
  deriving (Show)

defaultParticleState :: ParticleState
defaultParticleState =
  ParticleState
    { mass = 1,
      charge = 1,
      time = 0,
      posVec = zeroV,
      velocity = zeroV
    }

type OneBodyForce = ParticleState -> Vec

data DParticleState = DParticleState
  { dmdt :: R,
    dqdt :: R,
    dtdt :: R,
    drdt :: Vec,
    dvdt :: Vec
  }
  deriving (Show)

newtonSecondPS :: [OneBodyForce] -> ParticleState -> DParticleState -- difeq
newtonSecondPS fs st =
  let fNet = sumV [f st | f <- fs]
      m = mass st
      v = velocity st
      acc = fNet ^/ m
   in DParticleState
        { dmdt = 0,
          dqdt = 0,
          dtdt = 1,
          drdt = v,
          dvdt = acc
        }

-- z direction is toward the sky
-- assumes SI units
earthSurfaceGravity :: OneBodyForce
earthSurfaceGravity st =
  let g = 9.80665 -- m/s^2
   in (-mass st * g) *^ kHat

-- origin at center of sun
-- assumes SI units
sunGravity :: OneBodyForce
sunGravity (ParticleState m _q _t r _v) =
  let bigG = 6.67408e-11 -- N m^2/kg^2
      sunMass = 1.98848e30 -- kg
   in (-bigG * sunMass * m) *^ r ^/ magnitude r ** 3

airResistance ::
  R -> -- drag coefficient
  R -> -- air density
  R -> -- cross sectional area of object
  OneBodyForce
airResistance drag rho area (ParticleState _m _q _t _r v) =
  (-0.5 * drag * rho * area * magnitude v) *^ v
