import argparse, pandas as pd
import matplotlib.pyplot as plt
p = argparse.ArgumentParser()
p.add_argument('--in', dest='inp', required=True)
p.add_argument('--out', dest='out', required=True)
a = p.parse_args()
df = pd.read_csv(a.inp, sep='\t')
ax = df.plot(x=df.columns[0], y=df.columns[1], kind='bar')
ax.set_title('Example read counts')
plt.tight_layout()
plt.savefig(a.out)
print('Saved', a.out)
