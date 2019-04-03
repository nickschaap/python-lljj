import math
sum = 0

for number in range(981,1001,1) :
    sum = sum + number

print(sum + 1)
print(math.ceil(math.pow(2,20)/sum))