import fire,sys

class Solution:
  def isValid(self,string):
    pair={"(":")","[":"]","{":"}"}
    lst=list()
    if len(string)%2: return False
    for i in string:
      if i in pair.values():
        if not lst: return False
        if pair[lst[-1]]!=i: return False
        lst.pop()
      else:
        lst.append(i)
    return True
  
  def maxLeng(self,string):
    target=set("aeiouAEIOU")
    ans,anchor=0,-1
    for i,j in enumerate(string):
      if j in target:
        ans=max(ans,i-anchor)
      else:
        anchor=i
    return ans

  def skipHouseII(self,string,step):
    string.split(",")


if __name__=="__main__":
  fire.Fire(Solution)


