import numpy as np

a = np.arange(1, 10)
print(a[1:])
for i in enumerate(a[1:]):
    print(i)
print(a)