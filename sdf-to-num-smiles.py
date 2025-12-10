import sys

from rdkit import Chem  
from rdkit import Chem
from rdkit.Chem import rdmolops 
from rdkit.Chem import AllChem 

input_sdf = sys.argv[1] #'DRGXXXX_dH.sdf'  
output_smi = input_sdf[:-4]+'.smi'
output_sdf_recover = input_sdf[:-4]+'_recovered.sdf'

# 读取SDF文件  
suppl = Chem.SDMolSupplier(input_sdf)  
mol = suppl[0]  # 获取第一个分子  
  
# 为每个原子设置索引  
for atom in mol.GetAtoms():  
    atom.SetAtomMapNum(atom.GetIdx() + 1)  # 索引从1开始  
  
# 生成带原子索引的SMILES  
smiles_with_indices = Chem.MolToSmiles(mol)  
print(smiles_with_indices)
f1 = open(output_smi,'w')
f1.write(smiles_with_indices)
f1.close()



# validation 

# 从 SMILES 创建分子对象  
mol = Chem.MolFromSmiles(smiles_with_indices)  
  
params = AllChem.ETKDGv3()  
params.randomSeed = 42  # 确保可重现 
params.trackFailures = True  # 启用失败跟踪  
params.useRandomCoords = False

# 直接生成3D坐标（不添加H）  
cid = AllChem.EmbedMolecule(mol, params)  
  
print('cid:', cid)  
if cid < 0:  
    failures = params.GetFailureCounts()  
    print("失败统计:")  
    for i, count in enumerate(failures):  
        if count > 0:  
            print(f"  失败类型 {i}: {count} 次")

    params.useRandomCoords = True  
    AllChem.EmbedMolecule(mol, params)  

# 可选：使用力场优化结构  
cid = AllChem.UFFOptimizeMolecule(mol) 

print('cid:', cid)  
if cid < 0:  
    failures = params.GetFailureCounts()  
    print("失败统计:")  
    for i, count in enumerate(failures):  
        if count > 0:  
            print(f"  失败类型 {i}: {count} 次")

# 获取所有原子及其映射号  
atoms_with_map = []  
for atom in mol.GetAtoms():  
    map_num = atom.GetAtomMapNum()  
    atoms_with_map.append((map_num, atom.GetIdx()))  
  
# 按映射号排序  
atoms_with_map.sort(key=lambda x: x[0])  
  
# 创建新的原子顺序（基于映射号）  
new_order = [idx for (_, idx) in atoms_with_map]  
  
# 重新编号原子  
mol_renumbered = rdmolops.RenumberAtoms(mol, new_order)  
  
# 写入SDF  
with Chem.SDWriter(output_sdf_recover) as w:  
    w.write(mol_renumbered)
