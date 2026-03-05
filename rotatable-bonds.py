from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Lipinski import RotatableBondSmarts




# find all rotatable bonds 
l0 = mol.GetSubstructMatches(RotatableBondSmarts)

# find all amide peptide bonds 
amide_smarts = '[H]-[N;H1]-[CX3]=[OX1]'
amide_bonds = Chem.MolFromSmarts(amide_smarts)
l1 = mol.GetSubstructMatches(amide_bonds)

# remove amide peptide bonds from rotatable bonds (the middle two atoms in amide bonds should be removed from rotatable bonds)
# note the random ordering of two bonded atoms 
l2 = list(l0)
for b in l1:
    if (b[1],b[2]) in l0:
        l2.remove((b[1],b[2]))
    if (b[2],b[1]) in l0:
        l2.remove((b[2],b[1]))

print('verify number of rotatable bonds:', len(l0),'=',len(l1),'+',len(l2))

