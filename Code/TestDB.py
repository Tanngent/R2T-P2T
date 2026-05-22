testArray = []
with open("Test.txt") as f:
  for x in f:
    v = float(x.split(" ")[0])
    testArray.append(v)

testArray2 = []
with open("../Information/TPCH/Q18_1.txt") as f:
  for x in f:
    v = float(x.split(" ")[0])
    testArray2.append(v)

print(sorted(testArray) == sorted(testArray2))
