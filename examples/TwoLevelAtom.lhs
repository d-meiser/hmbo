\section{Hamiltonian of a two level atom}
\label{sec:TwoLevelAtom}

This basic example illustrates the usage of the hmbo library for the
representation of the Hamiltonian of a simple two level system.  The
operator we're trying to express is
\begin{equation}
\hat H = \frac{\delta}{2} \hat \sigma_z + \frac{g}{2} \hat \sigma_x\;.
\label{eqn:Hamiltonian}
\end{equation}
The $2\times 2$ matrices $\hat \sigma_z$ and $\hat \sigma_x$ are Pauli
matrices, $\delta$ is the energy difference between the excited and
ground state, and $g$ is the coupling strength of the ground state to
the excited state.

First, we need to include some modules:
\begin{code}
module Main where

import qualified Data.Vector.Unboxed as VU
import Data.Complex
import HMbo
\end{code}
{\tt Ket}'s are represented by unboxed vectors from the module {\tt
Data.Vector.Unboxed}.  The {\tt Data.Complex} module is needed to
construct quantum mechanical amplitudes.  Although the parameters
$\delta$ and $g$ in Eq. (\ref{eqn:Hamiltonian}) are real in our example,
they are still represented by complex numbers internally because that is
the general case.

We construct the Hamiltonian out of elementary operators by means of
the {\tt scale} and {\tt add}:
\begin{code}
buildHamiltonian :: Amplitude -> Amplitude -> LinearOp
buildHamiltonian delta g = fromJust $
  ((0.5 * delta) `scale` sigmaZ) `add` ((0.5 * g) `scale` sigmaX)
\end{code}
The operators {\tt sigmaZ} and {\tt sigmaX} are provided by the hmbo
library.  Note that the {\tt add} function does not directly return a
{\tt LinearOp}.  This is because this function fails if the two
operators to be added to one another have different dimensions.  In that
case we cannot construct a meaningful sum of the two operators.
Therefore the type of the {\tt add} function is
\begin{spec}
add :: LinearOp -> LinearOp -> Maybe LinearOp
\end{spec}
We return {\tt Nothing} if the addition of the two linear operators
fails.  For the purposes of this example we get rid of the {\tt Just} in
the result by means of the function
\begin{code}
fromJust :: Maybe a -> a
fromJust (Just a) = a
fromJust Nothing = error "Value was nothing."
\end{code}
This means that our program crashes with an error message if the
addition in {\tt buildHamiltonian} were to fail.  In more complicated
situations, especially when writing a library, it may be necessary to
handle failures more safely.

Our program simply builds the {\tt LinearOp} corresponding to the two
level atom Hamiltonian, computes its matrix in the canonical basis, and
prints the matrix to the screen:
\begin{code}
main :: IO ()
main = do
  let detuning = 4.0
  putStrLn ""
  putStr "Detuning: "
  print (realPart detuning)
  let g = 6.0
  putStrLn ""
  putStr "Coupling constant: "
  print (realPart g)
  let hamiltonian = buildHamiltonian detuning g
  let basis = [basisState 2 i | i <- [0..(2-1)]]
  let matrix = computeMatrix hamiltonian basis
  putStrLn ""
  putStrLn "Matrix:"
  printMatrix matrix
\end{code}

The basis states are straight forward to construct:
\begin{code}
basisState :: Int -> Int -> Ket
basisState d i = VU.fromList [kroneckerDelta i j | j <- [0..(d - 1)]]
  where
    kroneckerDelta m n | m == n = 1.0
                       | otherwise = 0.0
\end{code}

We obtain the matrix element $\langle \psi|\hat H|\phi\rangle$ by
applying the Hamiltonian to the state $|\phi\rangle$ and taking the
inner product of the complex conjugate of $|\psi\rangle$ with the result:
\begin{code}
matrixElement :: Ket -> LinearOp -> Ket -> Amplitude
matrixElement psi a phi = VU.foldl1 (+) $ VU.zipWith (*) psi' aPhi
  where
    aPhi = fromJust $ a `apply` phi
    psi' = VU.map conjugate psi
\end{code}
Note that we encounter a similar issue with {\tt apply} as we did with
{\tt add} earlier: this function can fail if the dimension of the linear
operator doesn't agree with the length of the vector.  In that case {\tt
apply} returns {\tt Nothing}.

From the individual matrix elements we build up the matrix by taking
matrix element between all pairs of basis states:
\begin{code}
computeMatrix :: LinearOp -> [Ket] -> Matrix
computeMatrix op b = [[matrixElement psi op phi | phi <- b] | psi <- b]
\end{code}
We naively represent the {\tt Matrix} by a list of lists of matrix
elements:
\begin{code}
type Matrix = [[Amplitude]]
\end{code}
We print the (magnitude) of the matrix to the screen by interpreting
each inner list as a row in the matrix:
\begin{code}
printMatrix :: Matrix -> IO ()
printMatrix m = putStrLn $ unlines $
  map (unwords . map (show .  magnitude) ) m
\end{code}

The program generates the following output:
\begin{spec}
*Main> :l TwoLevelAtom.lhs
[1 of 1] Compiling Main             ( TwoLevelAtom.lhs, interpreted )
Ok, modules loaded: Main.
*Main> main

Detuning: 4.0

Coupling constant: 6.0

Matrix:
2.0 3.0
3.0 2.0
\end{spec}
