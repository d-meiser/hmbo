module Main where

import qualified Data.Vector.Unboxed as VU
import Data.Complex

import HMbo

detuning = 4.0
g = 6.0 
Just hamiltonian = ((0.5 * detuning) `scale` sigmaZ) `add` ((0.5 * g) `scale` sigmaX)

basisState :: Int -> Int -> Ket
basisState d i = VU.fromList [delta i j | j <- [0..(d - 1)]]
  where
    delta i j | i == j = 1.0
              | otherwise = 0.0

matrixElement :: Ket -> LinearOp -> Ket -> Amplitude
matrixElement psi a phi = VU.foldl1 (+) $ VU.zipWith (*) psi' aPhi 
  where
    Just aPhi = a `apply` phi
    psi' = VU.map conjugate psi

main :: IO ()
main = do
  putStr "Detuning: "
  putStrLn $ show $ realPart detuning
  putStr "Coupling constant: "
  putStrLn $ show $ realPart g
  let basis = [basisState 2 i | i <- [0..(2-1)]]
  let matrix = [[matrixElement psi hamiltonian phi |
                 phi <- basis] | psi <- basis]
  putStrLn "Matrix:"
  putStrLn $ unlines $ map (unwords . map (show . magnitude) ) matrix
