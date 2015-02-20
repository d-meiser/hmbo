module HMbo.DimSpec (spec) where

import Test.Hspec
import HMbo

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
