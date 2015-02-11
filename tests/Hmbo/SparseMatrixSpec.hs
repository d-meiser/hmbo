module Hmbo.SparseMatrixSpec (spec) where

import Test.Hspec
import Hmbo.SparseMatrix

spec :: Spec
spec = do
  describe "SparseMatrix" $ do
    it "Zero Matrix has no entries" $ do
      numEntries zeroMatrix `shouldBe` 0
