{-# OPTIONS -Wall #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module MultipleObjects where

import Mechanics1D
  ( Diff (..),
    Mass,
    NumericalMethod,
    RealVectorSpace (..),
    TimeStep,
    euler,
  )
import Mechanics3D
  ( DParticleState (..),
    HasTime (..),
    OneBodyForce,
    ParticleState (..),
    defaultParticleState,
    newtonSecondPS,
  )
import SimpleVec
  ( R,
    Vec,
    magnitude,
    zeroV,
    (*^),
    (<.>),
    (><),
    (^*),
    (^+^),
    (^-^),
    (^/),
  )
