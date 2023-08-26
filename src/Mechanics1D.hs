{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# OPTIONS_GHC -Wall #-}

module Mechanics1D where

import Newton2 (fAir)
import Types (R)

type Time = R

type TimeStep = R

type Mass = R

type Position = R

type Velocity = R

type Force = R

type State1D = (Time, Position, Velocity)

newtonSecond1D ::
  Mass ->
  [State1D -> Force] -> -- force funcs
  State1D -> -- current state
  (R, R, R) -- deriv of state
newtonSecond1D m fs (t, x0, v0) =
  let fNet = sum [f (t, x0, v0) | f <- fs]
      acc = fNet / m
   in (1, v0, acc)

euler1D ::
  R -> -- time step dt
  (State1D -> (R, R, R)) -> -- differential equation
  (State1D -> State1D) --state update function
euler1D dt deriv (t0, x0, v0) =
  let (_, _, dvdt) = deriv (t0, x0, v0)
      t1 = t0 * dt
      x1 = x0 + v0 * dt
      v1 = v0 + dvdt * dt
   in (t1, x1, v1)

updateTXV ::
  R -> -- time interval dt
  Mass ->
  [State1D -> Force] -> -- list of force funcs
  (State1D -> State1D) -- state update func
updateTXV dt m fs = euler1D dt (newtonSecond1D m fs)

statesTXV ::
  R -> -- time interval dt
  Mass ->
  State1D -> -- initial state
  [State1D -> Force] -> -- list of force funcs
  [State1D] -- infinite list of states
statesTXV dt m txv0 fs = iterate (updateTXV dt m fs) txv0

-- assume that dt is the same between adjacent pairs
velocity1D ::
  [State1D] -> -- infinite list
  (Time -> Velocity) -- velocity function
velocity1D sts t =
  let (t0, _, _) = sts !! 0
      (t1, _, _) = sts !! 1
      dt = t1 - t0
      numSteps = abs $ round (t / dt)
      (_, _, v0) = sts !! numSteps
   in v0

velocityFtxv ::
  R -> -- time step
  Mass ->
  State1D -> -- initial state
  [State1D -> Force] -> -- list of force funcs
  (Time -> Velocity) -- velocity function
velocityFtxv dt m txv0 fs = velocity1D (statesTXV dt m txv0 fs)

-- assume that dt is the same between adjacent pairs
position1D ::
  [State1D] -> -- infinite list
  (Time -> Position) -- position function
position1D sts t =
  let (t0, _, _) = sts !! 0
      (t1, _, _) = sts !! 1
      dt = t1 - t0
      numSteps = abs $ round (t / dt)
      (_, x0, _) = sts !! numSteps
   in x0

positionFtxv ::
  R -> -- time step
  Mass ->
  State1D -> -- initial state
  [State1D -> Force] -> -- list of force funcs
  (Time -> Position) -- position function
positionFtxv dt m txv0 fs = position1D (statesTXV dt m txv0 fs)

springForce :: R -> State1D -> Force
springForce k (_, x0, _) = -k * x0
