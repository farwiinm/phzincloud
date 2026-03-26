import pandas as pd
df = pd.read_csv('results/results_cloudrun_pkad.csv')
tier1 = df[df['pka_tier'] == 1]
print(f'Tier 1 hits: {len(tier1)}')
if len(tier1) > 0:
    print(tier1[['pdb_id','residue_name','residue_seq','pka_value','pka_source']].to_string(index=False))
print(f'Tier 2: {len(df[df["pka_tier"]==2])} ({len(df[df["pka_tier"]==2])/len(df)*100:.1f}%)')
print(f'Tier 3: {len(df[df["pka_tier"]==3])} ({len(df[df["pka_tier"]==3])/len(df)*100:.1f}%)')
