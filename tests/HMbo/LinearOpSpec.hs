module HMbo.LinearOpSpec (spec) where

import Test.Hspec
import HMbo
import Data.Maybe
import qualified Data.Vector.Unboxed as VU

spec :: Spec
spec = do
  describe "Dim" $ do
    it "Can't be built from negative integeter." $ do
      toDim (-1) `shouldBe` Nothing
    it "Can't be built from 0." $ do
      toDim 0 `shouldBe` Nothing
    it "Can be built from positive integer." $ do
      toDim 2 `shouldSatisfy` (Nothing /=)
    it "Has the right dimension." $ do
      toDim 2 `shouldSatisfy` ((Just 2 ==) . (fmap fromDim))

  describe "Null Matrix" $ do
    it "Can be built from a dimension." $ do
      getDim (zero (fromJust $ toDim 2)) `shouldBe` (fromJust $ toDim 2)
    it "Does not change dimension when scaled." $ do
      getDim (scale 3.7 $ zero (fromJust $ toDim 2)) 
        `shouldBe` (fromJust $ toDim 2)

  describe "Identity Matrix" $ do
    it "Can be built from a dimension." $ do
      getDim (identity (fromJust $ toDim 2)) `shouldBe` (fromJust $ toDim 2)
    it "Does not change dimension when scaled." $ do
      getDim (scale 3.7 $ identity (fromJust $ toDim 2))
        `shouldBe` (fromJust $ toDim 2)

  describe "kron" $ do
    it "Multiplies dimensions." $ do
      getDim (
          kron (identity (fromJust $ toDim 3))
          (identity (fromJust $ toDim 3))) `shouldBe`
            (fromJust $ toDim 9)
    it "Results in ScaledId when both of the operators are identites." $ do
      kron
        (identity (fromJust $ toDim 3))
        (identity (fromJust $ toDim 3))
          `shouldBe`
            identity (fromJust $ toDim 9)

  describe "add" $ do
    it "Results in Nothing if dimensions don't match." $ do
      add
        (identity (fromJust $ toDim 2))
        (identity (fromJust $ toDim 3))
          `shouldBe`
            Nothing

  describe "apply" $ do
    it "Returns nothing when dimensions don't match." $ do
      let v = (VU.fromList [1.3, 3.4]) :: Ket
      let d = fromJust $ toDim 3
      apply (identityMatrix d) v `shouldBe` Nothing
    it "Returns something when dimensions match." $ do
      let v = (VU.fromList [1.3, 3.4]) :: Ket
      let d = fromJust $ toDim 2
      apply (identityMatrix d) v `shouldSatisfy` not . isNothing

    describe "Identity" $ do
      it "Returns vectors unchanged." $ do
        let v = (VU.fromList [1.3, 3.4]) :: Ket
        let d = fromJust $ toDim 2
        apply (identityMatrix d) v `shouldBe` Just v

    describe "Zero operator" $ do
      it "Returns 0." $ do
        let v = (VU.fromList [1.3, 3.4]) :: Ket
        let w = (VU.replicate 2 0) :: Ket
        let d = fromJust $ toDim 2
        (fromJust $ apply (zero d) v) `shouldBe` w
