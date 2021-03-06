\section{Hamiltonian of a Heisenberg spin chain}
\label{sec:SpinChain}

In this example we conside the Heisenberg model for a one dimensional
magnetic material.  The Heisenberg model consists of $N$ spin $1/2$
systems.  The Hamiltonian for the Heisenberg spin chain is given by
\begin{equation}
\hat H = -\frac{1}{2}\sum_{j=1}^N\left(
    J_x \hat \sigma_x^{(j)} \hat \sigma_x^{(j+1)} +
    J_y \hat \sigma_y^{(j)} \hat \sigma_y^{(j+1)} +
    J_z \hat \sigma_z^{(j)} \hat \sigma_z^{(j+1)} -
    h\hat \sigma_z^{(j)}\right)\;.
\label{eqn:HeisenbergSpinChain}
\end{equation}
The $2\times 2$ operators $\hat \sigma_x^{(j)}$, $\hat \sigma_y^{(j)}$,
and $\hat \sigma_z^{(j)}$ are the usual Pauli matices acting on the
degrees of freedom of the $j$-th spin.  We assume periodic boundary
conditions so that the last spin with index $j=N$ couples to the first
spin.The coupling constants $J_x$, $J_y$, and $J_z$ describe the
coupling of neighbouring spins.  For a ferromagnetic system the coupling
constants are positive.  For an antiferromagnetic system they are
negative.  The last term corresponds to an external magnetic field along
the $z$ axis.

First, we include some modules:
\begin{code}
module Main where

import qualified Data.Vector.Unboxed as VU
import Data.Complex (conjugate)
import Data.Maybe (fromJust)
import Data.List (foldl')
import Data.Monoid (mconcat)
import HMbo
\end{code}
We use {\tt Data.Vector.Unboxed} to deal with hmbo's {\tt Ket}s and
Data.Complex for quantum mechanical amplitudes.  We use the standard
{\tt fromJust} function to deal with the possibility of failure in some
of the HMbo operations.  The {\tt foldlM} function is useful for
expressing the sums over terms in Eq.~(\ref{eqn:HeisenbergSpinChain}).
Finally, we obviously need the HMbo module.

Each spin is 2 dimensional quantum system.  We represent the Hilbert
space of the total system by a list of the dimensions of each subsystem.
To build the Hilbert space for our spin chain we simply generate a list
with $N$ copies of the single spin Hilbert space:
\begin{code}
dim :: Dim
dim = fromJust $ toDim 2

type HilbertSpace = [Dim]

spinSpace :: Int -> HilbertSpace
spinSpace numSpins = replicate numSpins dim
\end{code}
Note that in our particular case a simpler representation of the total
system is possible because the subsystems are identical.  We do not take
advantage of this for now.

A central operation in our problem and in most many body quantum
mechanical calculations is the embedding of single particle operators
into the many body Hilbert space.  For example, the single particle
Pauli matrix $\sigma_x$ becomes an operator acting on spin $j$ by the
construction
\begin{equation}
\hat \sigma_x^{(j)} =
\underbrace{\mathds{1}_2 \otimes\cdots\otimes\mathds{1}_2}_{j - 1}
\otimes\sigma_x\otimes
\underbrace{\mathds{1}_2\otimes\cdots\otimes\mathds{1}_2}_{N - j}\;,
\end{equation}
where $\mathds{1}_2$ is the $2\times 2$ identity matrix.

This construction can be expressed directly using the hmbo library:
\begin{code}
embed :: ManyBodyOperator -> Int -> HilbertSpace -> ManyBodyOperator
embed op j dims = mconcat $ zipWith factor [0..] dims
  where
    factor i d | i == j = op
               | otherwise = eye d
\end{code}
We wish to embed the single particle operator {\tt op} into the $j$-th
``slot'' in the product space.  We recurse through the dimensions in the
product space.  When we reach the index $j$ we insert the single
particle operator.  Otherwise we insert an identity of the appropriate
dimension.  The operators are multiplied together using the {\tt kron}
function.  The recursion is terminated with an identity of dimension
one, which is the identity element for {\tt kron}.

Now we have all the ingredients in place to construct the Heisenberg
spin chain Hamiltonian.  We start with the external magnetic field
piece,
\begin{equation}
\hat H_z = \frac{1}{2}h\sum_{j=1}^N\sigma_z^{(j)}\;.
\end{equation}
\begin{code}
hz :: Int -> Amplitude -> ManyBodyOperator
hz n h = scale (0.5 * h) $ foldl' add (zero totalDim)
  [embed sigmaZ j space | j <- [0..(n - 1)]]
  where
    space = spinSpace n
    totalDim = fromJust $ toDim (2^n)
\end{code}
We insert {\tt sigmaZ} into each slot $j=0,1,\ldots,N-1$, add the
results up using {\tt foldlM add (zero totalDim)}, and scale the result
by $h/2$.

The interaction pieces have a rather similar form for the $x$, $y$, and
$z$ components suggesting that we can deal with them using the same
function.
\begin{code}
buildInteractionPiece :: ManyBodyOperator -> Int -> Amplitude ->
                         ManyBodyOperator
buildInteractionPiece op n coupling = scale ((-0.5) * coupling) $
  foldl' add (zero totalDim)
    [embed op j space `mul` embed op ((j + 1) `mod` n) space
    | j <- [0..(n - 1)]]
  where
    space = spinSpace n
    totalDim = fromJust $ toDim (2^n)
\end{code}
For each $j$ we construct two embeddings, $\hat op^{(j)}$ and $\hat
op^{(j+1)}$, and multiply the two together.  All these terms are then
added up exactly as in the case of $\hat H_z$ above.  We can now build
the complete Hamiltonian:
\begin{code}
buildH :: Int -> Amplitude -> Amplitude -> Amplitude -> Amplitude ->
          ManyBodyOperator
buildH n h jx jy jz = foldl' add (zero totalDim)
  [hz n h
  ,buildInteractionPiece sigmaX n jx
  ,buildInteractionPiece sigmaY n jy
  ,buildInteractionPiece sigmaZ n jz
  ]
  where
    totalDim = fromJust $ toDim (2^n)
\end{code}

To illustrate how the operator might be used we compute the expectation
value of the Hamiltonian in the state with all spins in the down state:
\begin{code}
main :: IO ()
main = print
  [matrixElement
     (basisState (2^n) 0)
     (buildH n 1.0 2.0 3.0 4.0)
     (basisState (2^n) 0) | n <- [1..10]]
\end{code}

The functions {\tt matrixElement} and {\tt basisState} were introduced
already in {\it TwoLevelAtom.lhs}.  For completeness, we reprint their
code in the following:
\begin{code}
basisState :: Int -> Int -> Ket
basisState d i = VU.fromList [kroneckerDelta i j | j <- [0..(d - 1)]]
  where
    kroneckerDelta m n | m == n = 1.0
                       | otherwise = 0.0

matrixElement :: Ket -> ManyBodyOperator -> Ket -> Amplitude
matrixElement psi a phi = VU.foldl1 (+) $ VU.zipWith (*) psi' aPhi
  where
    aPhi = a `apply` phi
    psi' = VU.map conjugate psi
\end{code}
