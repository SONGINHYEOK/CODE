
from random import random, randrange

import time
import numpy as np
import pandas as pd
import random
start = time.time()

c = np.random.randint(0, 1000000000000000000, size=(4000000,24))

df = pd.DataFrame(c)
df[23]=1
df.head()
df  = df.drop(0, axis=1)
df.head()

df.to_csv("./sample_400.csv", index=False, header=False)
end = time.time() - start
print('Running time;', end)
print("Job done")    




