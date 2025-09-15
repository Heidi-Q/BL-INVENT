from rdkit import Chem
import os
import argparse
import re


def get_mol_label_idx(mol):
    for idx, atom in enumerate(mol.GetAtoms()):
        result_idx = 0
        if atom.GetIsotope() != 0:
            return idx

def add_atom_label(file, warhead_ori_smi):
    f = open(file)
    lines = f.readlines()
    f.close()

    warhead_ori_mol = Chem.MolFromSmarts(warhead_ori_smi)
    # E3_ligand_ori_mol = Chem.MolFromSmarts(E3_ligand_ori_smi)

    warhead_isotope_idx = get_mol_label_idx(warhead_ori_mol)

    # E3_ligand_isotope_idx = get_mol_label_idx(E3_ligand_ori_mol)


    warhead_true_smi = re.sub('\[[0-9]+(.*?)\]', r"\1", warhead_ori_smi)
    # E3_ligand_true_smi = re.sub('\[[0-9]+(.*?)\]', r"\1", E3_ligand_ori_smi)

    warhead_true_patt = Chem.MolFromSmarts(warhead_true_smi)
    # E3_ligand_true_patt = Chem.MolFromSmarts(E3_ligand_true_smi)


    result_list = []

    for idx, line in enumerate(lines):
        print(idx)
        protac_mol = Chem.MolFromSmiles(line)

        warhead_hit_list = protac_mol.GetSubstructMatches(warhead_true_patt)

        # warhead_hits = list(protac_mol.GetSubstructMatches(warhead_true_patt))
        # E3_ligand_hits = list(protac_mol.GetSubstructMatch(E3_ligand_true_patt))

        protac_warhead_label_idx = warhead_hits[warhead_isotope_idx]
        # protac_E3_label_idx = E3_ligand_hits[E3_ligand_isotope_idx]

        warhead_atom = protac_mol.GetAtomWithIdx(protac_warhead_label_idx)
        # E3_atom = protac_mol.GetAtomWithIdx(protac_E3_label_idx)

        warhead_atom.SetIsotope(100)
        # E3_atom.SetIsotope(101)

        for neigh_atom in warhead_atom.GetNeighbors():
            # print(neigh_atom.GetIdx())
            if neigh_atom.GetIdx() not in warhead_hits:
                neigh_atom.SetIsotope(100)

        # for neigh_atom in E3_atom.GetNeighbors():
        #     if neigh_atom.GetIdx() not in E3_ligand_hits:
        #         neigh_atom.SetIsotope(101)
        new_protac_smi = Chem.MolToSmiles(protac_mol)

        result_list.append(new_protac_smi+'\n')

    file_dir, file_name = os.path.split(file)
    result_file = os.path.join(file_dir, file_name.split('.')[0]+'_add_label.csv')

    w = open(result_file, 'w')
    w.writelines(result_list)
    w.close()


if __name__ == '__main__':
    # base_dir = '/data/baiqing/PycharmProjects/Reinvent-master-3.2/data/protac/BAF'
    #
    # file_dir = os.path.join(base_dir, 'linkinvent_7_20_scaffold_memory_top100.csv')
    #
    # result_file = os.path.join(base_dir, 'linkvent_7_20_scaffold_memory_top100_add_label.csv')


    # file_dir = os.path.join(base_dir, 'linkinvent_TL_piperidine_scaffold_memory_top100.csv')
    # result_file = os.path.join(base_dir, 'linkinvent_TL_piperidine_scaffold_memory_top100_add_label.csv')

    # warhead = 'Nc1nnc(cc1N1CC[15N](*)CC1)-c1ccccc1O'
    # E3_ligand = 'Cc1ncsc1-c1ccc(CNC(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@@H](NC(=O)C2(F)CC2)C(C)(C)C)c([18O](*))c1'

    # warhead = 'Nc1nnc(cc1N1CC[15N]CC1)-c1ccccc1O'
    # E3_ligand = 'Cc1ncsc1-c1ccc(CNC(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@@H](NC(=O)C2(F)CC2)C(C)(C)C)c([18O])c1'



    # add_atom_label(file_dir, warhead, E3_ligand, warhead_true, E3_ligand_true, result_file)

    parser = argparse.ArgumentParser(description='Generate scaffold memory add label')
    parser.add_argument('-i', '--input', type=str, default=None, help='Specify scaffold_memory csv file')

    parser.add_argument('-w', '--warhead', type=str, default=None, help='Specify warhead smiles')
    # parser.add_argument('-e3', '--E3_ligand', type=str, default=None, help='Specify E3-ligand smiles')

    # parser.add_argument('-o', '--output', type=str, default=None, help='Specify scaffold_memory_add_label csv file')
    # parser.add_argument('-f', '--file', type=str, help='Specify the filename for the output_file')

    # args = parser.parse_args(['-i', '/data/baiqing/PycharmProjects/Reinvent-master-3.2/data/protac/BAF/linkinvent_7_20_scaffold_memory_top100.csv',
    #                           '-w', 'Nc1nnc(cc1N1CC[15N]CC1)-c1ccccc1O',
    #                           '-e3', 'Cc1ncsc1-c1ccc(CNC(=O)[C@@H]2C[C@@H](O)CN2C(=O)[C@@H](NC(=O)C2(F)CC2)C(C)(C)C)c([18O])c1',
    #                           ])

    # args = parser.parse_args(['-i', '/data/baiqing/PycharmProjects/Reinvent-master-3.2/data/protac/BTK/linkinvent_7_15_3_scaffold_memory_top100.csv',
    #                           '-w', 'NC(=O)c1c(N)n(nc1-c1ccc(Oc2ccc(F)cc2F)cc1)[C@@H]1CCCN(C1)[100C]=O',
    #                           '-e3', 'CN[C@@H](C)C(=O)N[C@H](C(=O)N1Cc2cc([101O])ccc2C[C@H]1C(=O)N[C@@H]1CCCc2ccccc12)C(C)(C)C',
    #                           ])

    # args = parser.parse_args(['-i', '/data/baiqing/PycharmProjects/Reinvent-master-3.2/data/protac/BTK/linkinvent_no_constraint_3_scaffold_memory_top100.csv',
    #                           '-w', 'NC(=O)c1c(N)n(nc1-c1ccc(Oc2ccc(F)cc2F)cc1)[C@@H]1CCCN(C1)[100C]=O',
    #                           '-e3', 'CN[C@@H](C)C(=O)N[C@H](C(=O)N1Cc2cc([101O])ccc2C[C@H]1C(=O)N[C@@H]1CCCc2ccccc12)C(C)(C)C',
    #                           ])

    # args = parser.parse_args(['-i', '/data/baiqing/PycharmProjects/Reinvent-master-3.2/data/protac/BTK/linkinvent_sample.csv',
    #                           '-w', 'NC(=O)c1c(N)n(nc1-c1ccc(Oc2ccc(F)cc2F)cc1)[C@@H]1CCCN(C1)[100C]=O',
    #                           '-e3', 'CN[C@@H](C)C(=O)N[C@H](C(=O)N1Cc2cc([101O])ccc2C[C@H]1C(=O)N[C@@H]1CCCc2ccccc12)C(C)(C)C',
    #                           ])

    # args = parser.parse_args(['-i', '/data/baiqing/PycharmProjects/Protac-invent/data/protac/BTK/linkinvent_7_15_3_scaffold_memory_101_200.csv',
    #                           '-w', 'NC(=O)c1c(N)n(nc1-c1ccc(Oc2ccc(F)cc2F)cc1)[C@@H]1CCCN(C1)[100C]=O',
    #                           '-e3', 'CN[C@@H](C)C(=O)N[C@H](C(=O)N1Cc2cc([101O])ccc2C[C@H]1C(=O)N[C@@H]1CCCc2ccccc12)C(C)(C)C',
    #                           ])

    # args = parser.parse_args(['-i', '/data/baiqing/PycharmProjects/Protac-invent/data/protac/BAF/linkinvent_weight1_scaffold_memory_top200.csv',
    #                           '-w', 'C1C[100N]CCN1c2cc(nnc2N)-c3c(O)cccc3',
    #                           '-e3', 'C1CC1(F)C(=O)N[C@@H](C(C)(C)C)C(=O)N(C[C@@H](C2)O)[C@@H]2C(=O)NCc3ccc(cc3[101O])-c(c4C)scn4',
    #                           ])

    # args = parser.parse_args(['-i', '/data/baiqing/PycharmProjects/Protac-invent/data/protac/5T35/linkinvent_weight1_scaffold_memory_top200.csv',
    #                           '-w', 'c1cc(Cl)ccc1C(=N[C@H]2C[100C]=O)c(c(C)c(s3)C)c3n(c24)c(nn4)C',
    #                           '-e3', 'N1CSC(=C1C)c2ccc(cc2)CNC(=O)[C@@H]3C[C@@H](O)CN3C(=O)[C@H](C(C)(C)C)N[101C]=O',
    #                           ])
    args = parser.parse_args(['-i', '/data/baiqing/PycharmProjects/Protac-invent/data/MsbA/top200.csv',
                              '-w', 'O=C([O-])c1c[101c](*)ccc1NS(=O)(=O)c2ccc(cc2)CCCC'
                              ])



    add_atom_label(args.input, args.warhead)





