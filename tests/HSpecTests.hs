module Main where
 
import HMbo
import Test.Hspec
 
main :: IO ()
main = hspec $ do
  describe "Import HMbo module" $ do
    it "Dummy test" $ do
      (1 :: Int) `shouldBe` (1 ::Int)
