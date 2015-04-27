\section{Ket}
\label{sec:Ket}

\begin{code}
module HMbo.Ket (
     Ket
    ,transpose) where

import HMbo.Amplitude
import qualified Data.Vector.Unboxed as VU

type Ket = VU.Vector Amplitude
transpose :: Int -> Ket -> Ket
transpose n v = VU.fromList
                  [(VU.!) v (i * n + j) | j <- [0..(n - 1)], i <- [0..(m - 1)]]
  where
    m = VU.length v `div` n
\end{code}
