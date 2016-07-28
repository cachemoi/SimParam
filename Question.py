import pandas as pd
import numpy


df = pd.DataFrame(numpy.random.randn(2, 2), columns=['A', 'B'])

df.columns = pd.MultiIndex.from_tuples(zip(['AA', 'BB'], df.columns))

print(df)