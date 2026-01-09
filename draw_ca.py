from rdkit import Chem
from rdkit.Chem import Draw
import pandas as pd
import matplotlib.pyplot as plt

# === Step 1. 定义 SMARTS 列表 ===
smarts_list = [
    "[#8][#8]",
    "[#6;+]",
    "[#16][#16]",
    "[#7;!n][S;!$(S(=O)=O)]",
    "[#7;!n][#7;!n]",
    "[C,c]=,#[C,c]",
    "[NX3&H0]",
    "[C,O,S]C(=[O,S])[O,S]",
    "[#7;!n]C(=[O,S])[S]",
    "[#7;!n]C(=S)[O]",
    "S(=O)(=O)O",
    "[#7;!n][C;!$(C(=[O,N])[N,O])][#16;!s]",
    "[#7;!n][C;!$(C(=[O,N])[N,O])][#7;!n]",
    "[#7;!n][C;!$(C(=[O,N])[N,O])][#8;!o]",
    "[#8;!o][C;!$(C(=[O,N])[N,O])][#16;!s]",
    "[#8;!o][C;!$(C(=[O,N])[N,O])][#8;!o]",
    "[#16;!s][C;!$(C(=[O,N])[N,O])][#16;!s]"
]

# === Step 2. 生成分子对象 ===
mols = []
for s in smarts_list:
    try:
        mol = Chem.MolFromSmarts(s)
        mols.append(mol)
    except:
        mols.append(None)

# === Step 3. 生成分子图 ===
img = Draw.MolsToGridImage(
    mols,
    molsPerRow=4,
    subImgSize=(250, 250),
    legends=[f"{i+1}. {s}" for i, s in enumerate(smarts_list)]
)

# === Step 4. 保存为图片 ===
img.save("custom_alerts.png")
print("✅ 图片已保存为 custom_alerts.png")

# === Step 5. 可选：生成表格（带编号和 SMARTS） ===
df = pd.DataFrame({"Index": range(1, len(smarts_list)+1), "SMARTS": smarts_list})
df.to_csv("custom_alerts.csv", index=False)
print("✅ SMARTS 列表已保存为 custom_alerts.csv")
