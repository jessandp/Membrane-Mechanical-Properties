import os
protein_file = os.path.join("Users", "jessicapallarez", "Desktop", "Lab Research", "Jessica", "P-pdbs", "CG_DOPC_1Ch25.quarter.pdb")
print(protein_file)

with open(protein_file, "r") as outfile:
#Places file contents into a list where each element is a line of the file
    data = outfile.readlines()

#Finding number of atoms
for line in data:
        if 'PROTEIN ATOMS' in line:
            PROTEIN_ATOM_line = line
            words = PROTEIN_ATOM_line.split(':')
            print(PROTEIN_ATOM_line)
