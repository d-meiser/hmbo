module HMbo.LinearOpSpec (spec) where

import Test.Hspec
import HMbo
import Data.Maybe

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
    it "Can be build from a dimension." $ do
      getDim (zero (fromJust $ toDim 2)) `shouldBe` (fromJust $ toDim 2)

