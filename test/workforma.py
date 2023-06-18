from cmath import nan
import pandas as pd
import os
pah = './ma'
by = '出资比例'
excels = [pd.read_excel(pah+'/'+fname)
          for fname in os.listdir(pah) if '.xls' in fname]
df = pd.concat(excels)
df[by] = pd.to_numeric(df[by].str.replace('-', '').str.replace('%', 'e-2'))
df = df.sort_values(by=by, ascending=False, na_position='last')
df.to_excel('汇总.xlsx', index=False)
