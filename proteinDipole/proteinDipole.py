#!/usr/bin/env python

import re
import argparse
import numpy as np
import sys

from datetime import datetime
startTime = datetime.now()



def distance(x0, x1):
    '''
    calculates distance (no PBC).
    '''
    x0 = np.array(x0)
    x1 = np.array(x1)
    delta = np.abs(x1 - x0)
    #delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
    dist = np.sqrt((delta ** 2).sum(axis=-1))
    return dist


def com(coord_list, mass_list):
    '''
    calculates centre of mass from a list of coords
    and a list of accompanying masses.
    '''
    x_cm = 0
    y_cm = 0
    z_cm = 0
    tot_mass = 0
    for coord, mass in zip(coord_list, mass_list):
        x_cm += mass * coord[0]
        y_cm += mass * coord[1]
        z_cm += mass * coord[2]
        tot_mass += mass
    x_cm = float(x_cm) / float(tot_mass)
    y_cm = float(y_cm) / float(tot_mass)
    z_cm = float(z_cm) / float(tot_mass)

    cm = (x_cm, y_cm, z_cm)
    #print ('com: ', cm, tot_mass)

    return cm



def atomicCharges(pdbList, pH):
    '''
    Get fraction of charge on ionisable atoms based on input pH,
    used in calcMultiBinding function.

    pdList format: [atom_name, resname, resid, x, y, z, 
                        occupancy, beta]
    '''

    residueTypesDict = {
            'aliphatics': ['ALA', 'GLY', 'ILE', 'LEU', 'PRO', 'VAL'],
            'aromatics': ['PHE', 'TRP', 'TYR'],
            'acidics': [['ASP', 'CG', 4], ['GLU', 'CD', 4.4]],
            'basics': [['ARG', 'CZ', 12], ['LYS', 'NZ', 10.4], 
                    ['HIS', 'CE1', 6.3]],
            'positives': [['ARG', 'CZ', 12], ['LYS', 'NZ', 10.4]],
            'hydroxylics': ['SER', 'THR'],
            'sulfites': ['CYS', 'MET'],
            'amidics': ['ASN', 'GLN'],
            'terminals': [['NT', 7.5], ['CT', 3.8]]
    }


    sequence = []

    for atom in pdbList:
        atomid = atom[0]
        atom_name = atom[1]
        resname = atom[2]
        resid = atom[3]
        coords = [atom[4], atom[5], atom[6]]
        occupancy = atom[7]
        beta = atom[8]
        for key, value in residueTypesDict.items():
            if key == 'acidics':
                for v in value:
                    pka = v[2]
                    if (resname, atom_name) == (v[0], v[1]):
                        exp_here = pH - pka
                        ratio_neg = 10 ** exp_here
                        frac_neg = -(float(1) / float(1 + (float(1) / 
                                float(ratio_neg))))
                        sequence.append((atomid, atom_name, resname, resid, 
                                coords, occupancy, beta, frac_neg))
                    else:
                        continue
            if key == 'basics':
                for v in value:
                    pka = v[2]
                    if (resname, atom_name) == (v[0], v[1]):
                        exp_here = pka - pH
                        ratio_pos = 10 ** exp_here
                        frac_pos = float(1) / float(1 + (float(1) / 
                                float(ratio_pos)))
                        sequence.append((atomid, atom_name, resname, resid, 
                                coords, occupancy, beta, frac_pos))
                    else:
                        continue
            if key == 'terminals':
                for v in value:
                    if (atom_name) == (v[0]):
                        pka = v[1]
                        if v[0] == 'NT':
                            exp_here = pka - pH
                            ratio_pos = 10 ** exp_here
                            frac_pos = float(1) / float(1 + (float(1) / 
                                    float(ratio_pos)))
                            sequence.append((atomid, atom_name, resname, resid, 
                                    coords, occupancy, beta, frac_pos))
                        if v[0] == 'CT':
                            exp_here = pH - pka
                            ratio_neg = 10 ** exp_here
                            frac_neg = -(float(1) / float(1 + (float(1) / 
                                    float(ratio_neg))))
                        sequence.append((atomid, atom_name, resname, resid, 
                                coords, occupancy, beta, frac_neg))
                    else:
                        continue

    return sequence




def proteinDipole(pdbList, pH, pdb, pqr, fileName):
    '''
    Get protein dipole
    '''

    atom_mass_dict = {
            'H': 1.0079,
            'C': 12.0107,
            'N': 14.0067,
            'O': 15.9994,
            'S': 32.065,

            }

    print ('Output/s for protein dipole at pH %s' % (pH))
    sequence = atomicCharges(pdbList, pH) 

    coord_list = []
    mass_list = []
    for atom in pdbList:
        #print (atom)
        atom_name = atom[1]
        try:
            mass = atom_mass_dict[atom_name[0]]
        except IndexError:
            print ('Atom type \'%s\' is not in \'atom_mass_dict\', '\
                    'setting mass to 1' % (atom[1][0]))
            mass = 1

        coords = [atom[4], atom[5], atom[6]]

        #if atom_name == 'CA':
        coord_list.append(coords)
        mass_list.append(mass)
        #else:
            #continue

    centre = com(coord_list, mass_list) #include mass info


    #####Atom charge dipole

    output = open('%s_dipole.pdb' % (fileName), 'w')

    output.write('\n'.join(['REMARK ATOMIC CHARGE DIPOLES']) + '\n')

    output.write('\n'.join(['{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}'\
            '   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2d}'\
            .format('ATOM', 0, 'CM', ' ', 'COM', 
            'A', 0, ' ', centre[0], centre[1], centre[2], 
            0.00, 0.00, ' ', 0)]) + '\n')

    if pdb == True:
        dipole_vector = [0, 0, 0]
        for atom in sequence:
            atomid = atom[0]
            atom_name = atom[1]
            resname = atom[2]
            resid = atom[3]
            coords = atom[4]
            occupancy = atom[5]
            beta = atom[6]
            charge = atom[7]

            output.write('\n'.join(['{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}'\
                    '{:1s}'\
                    '   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          '\
                    '{:>2s}{:2d}'\
                    .format('ATOM', atomid, atom_name, ' ', resname, 
                    'A', resid, ' ', coords[0], coords[1], coords[2], 
                    round(charge, 2), 0, ' ', 0)]) + '\n')


            dist_diff = np.array(coords) - np.array(centre)
            d = float(charge) * dist_diff
            dipole_vector += d


        dipole_vector_debye = np.array(dipole_vector) * 4.803 #convert to Debyes
        dipole_magnitude = distance(centre, dipole_vector)
        dipole_magnitude_debye = distance(centre, dipole_vector_debye)
        dipole_vector_norm = np.divide(dipole_vector, dipole_magnitude)

        dipole_vector = np.array(centre) + np.array(dipole_vector_norm) * \
                dipole_magnitude
        dipole_vector_norm = np.array(centre) + np.array(dipole_vector_norm)


        #print (centre, dipole_vector, dipole_vector_norm, 
                #dipole_magnitude_debye)

        output.write('\n'.join(['REMARK STANDARD pKa DIPOLE at pH %s' 
                % (pH)]) + '\n')

        output.write('\n'.join(['{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}'\
                '   {:8.0f}{:8.0f}{:8.0f}{:6.2f}{:6.2f}          {:>2s}{:2d}'\
                .format('ATOM', 0, 'DV', ' ', 'STD', 
                'A', 0, ' ', round(dipole_vector[0], 0), 
                round(dipole_vector[1], 0), round(dipole_vector[2], 0), 
                0.00, 0.00, ' ', 0)]) + '\n')
        output.write('\n'.join(['{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}'\
                '   {:8.0f}{:8.0f}{:8.0f}{:6.2f}{:6.2f}          {:>2s}{:2d}'\
                .format('ATOM', 0, 'DN', ' ', 'STD', 
                'A', 0, ' ', round(dipole_vector_norm[0], 0), 
                round(dipole_vector_norm[1], 0), 
                round(dipole_vector_norm[2], 0), 
                0.00, 0.00, ' ', 0)]) + '\n')
        output.write('\n'.join(['{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}'\
                '   {:8.0f}{:8.0f}{:8.0f}{:6.2f}{:6.2f}          {:>2s}{:2d}'\
                .format('REMARK', 0, 'DM', ' ', 'STD', 
                'A', 0, ' ', round(dipole_magnitude_debye, 0), 
                round(dipole_magnitude_debye, 0), 
                round(dipole_magnitude_debye, 0), 
                0.00, 0.00, ' ', 0)]) + '\n')

        print ('Protein dipole moment from standard pKa\'s is %s Debye' % 
                (int(dipole_magnitude_debye)))


    if pqr == True:
        dipole_vector = [0, 0, 0]
        for atomC in pdbList:
            charge = atomC[7]
            coords = [atomC[4], atomC[5], atomC[6]]
            dist_diff = np.array(coords) - np.array(centre)
            d = float(charge) * dist_diff
            dipole_vector += d



        dipole_vector_debye = np.array(dipole_vector) * 4.803 #convert to Debyes
        dipole_magnitude = distance(centre, dipole_vector)
        dipole_magnitude_debye = distance(centre, dipole_vector_debye)
        dipole_vector_norm = np.divide(dipole_vector, dipole_magnitude)

        dipole_vector = np.array(centre) + np.array(dipole_vector_norm) * \
                dipole_magnitude
        dipole_vector_norm = np.array(centre) + np.array(dipole_vector_norm)


        #print (centre, dipole_vector, dipole_vector_norm, 
                #dipole_magnitude_debye)

        output.write('\n'.join(['REMARK DIPOLE FROM PDB2PQR FILE']) + '\n')

        output.write('\n'.join(['{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}'\
                '   {:8.0f}{:8.0f}{:8.0f}{:6.2f}{:6.2f}          {:>2s}{:2d}'\
                .format('ATOM', 0, 'DV', ' ', 'PQR', 
                'A', 0, ' ', round(dipole_vector[0], 0), 
                round(dipole_vector[1], 0), round(dipole_vector[2], 0), 
                0.00, 0.00, ' ', 0)]) + '\n')
        output.write('\n'.join(['{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}'\
                '   {:8.0f}{:8.0f}{:8.0f}{:6.2f}{:6.2f}          {:>2s}{:2d}'\
                .format('ATOM', 0, 'DN', ' ', 'PQR', 
                'A', 0, ' ', round(dipole_vector_norm[0], 0), 
                round(dipole_vector_norm[1], 0), 
                round(dipole_vector_norm[2], 0), 
                0.00, 0.00, ' ', 0)]) + '\n')
        output.write('\n'.join(['{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}'\
                '   {:8.0f}{:8.0f}{:8.0f}{:6.2f}{:6.2f}          {:>2s}{:2d}'\
                .format('REMARK', 0, 'DM', ' ', 'PQR', 
                'A', 0, ' ', round(dipole_magnitude_debye, 0), 
                round(dipole_magnitude_debye, 0), 
                round(dipole_magnitude_debye, 0), 
                0.00, 0.00, ' ', 0)]) + '\n')

        print ('Protein dipole moment from PDB2PQR calculated charges '\
                'is %s Debye' % 
                (int(dipole_magnitude_debye)))


    return output




def readPDB(fileName, *args, **kwargs):
    '''
    read pdb file and extract data. 
    '''

    pH = kwargs.get('pH', 7)

    pdbList = []
    resNumList = []
    residueNumDict = {}
    num_pos_neg_points = [[], []]


    with open(fileName) as data:
        count = 0
        patch_point_total = 0
        patchCoordsList = []
        patch_num = None

        for line in data:
            if line[0:6] == 'ATOM  ':
                atomid = int(float(re.sub('[^-0-9.]', '', line[6:11])))
                atom_name = line[12:16].replace(' ','')
                resname = line[17:20].replace(' ','')
                chain = line[21:22].replace(' ','')
                resid = int(float(line[22:26].replace(' ','')))
                x = float(re.sub('[^-0-9.]', '', line[30:38]))
                y = float(re.sub('[^-0-9.]', '', line[38:46]))
                z = float(re.sub('[^-0-9.]', '', line[46:54]))

                try:
                    occupancy = float(line[54:60].replace(' ',''))
                    beta = float(line[60:66].replace(' ',''))
                except ValueError:
                    occupancy = None
                    beta = None

                pdbList.append([atomid, atom_name, resname, resid, x, y, z, 
                        occupancy, beta])


                if resname == 'ARG' and 'CZ' in atom_name and pH < 12:
                    num_pos_neg_points[0].append([resname, resid, 
                        atom_name])
                if resname == 'ASP' and 'CG' in atom_name and pH > 4:
                    num_pos_neg_points[1].append([resname, resid, 
                        atom_name])
                if resname == 'GLU' and 'CD' in atom_name and pH > 4.4:
                    num_pos_neg_points[1].append([resname, resid, 
                        atom_name])
                if resname == 'LYS' and 'NZ' in atom_name and pH < 10.4:
                    num_pos_neg_points[0].append([resname, resid, 
                        atom_name])

                if resname == 'HIS' and 'CE1' in atom_name and pH < 6:
                    num_pos_neg_points[0].append([resname, resid, 
                        atom_name])


                if resname not in residueNumDict and \
                        [resid, chain] not in resNumList:
                    residueNumDict[resname] = 1
                    resNumList.append([resid, chain])
                elif resname in residueNumDict and \
                        [resid, chain] not in resNumList:
                    residueNumDict[resname] += 1
                    resNumList.append([resid, chain])
                else:
                    continue


    data.close()



    if len(pdbList) != 0:
        return pdbList, residueNumDict, num_pos_neg_points




def outputData(fileName, *args, **kwargs):
    '''
    Use input data to run sequence binder and output protein
    residue info and binding info.

    '''

    pH = kwargs.get('pH', 7)

    fileName1 = fileName.replace('..', '')
    fileName1 = fileName1.replace('/', '')
    fileName_split = fileName1.split('.')
    #print (fileName_split)

    pdb, pqr = True, False

    if fileName_split[1] == 'pqr':
        pqr = True

    pdbList, residueNumDict, num_pos_neg_points = readPDB(fileName, pH=pH) 
    proteinDipole(pdbList, pH, pdb, pqr, fileName_split[0])



def main():

    try:
        usage = 'proteinDipole.py [-h]'
        parser = argparse.ArgumentParser(description='Tool for calculating '
                'the dipole moment of proteins at a given pH. '\
                '\nExample:'\
                '\npython proteinDipole.py --fileName file.pdb --pH 7', 
                usage=usage, 
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        group = parser.add_argument_group('Options')
        group = parser.add_argument('-f', '--fileName', action='store', 
                default='file.pdb', help='input pdb list file name')
        group = parser.add_argument('-pH', '--pH', action='store', 
                type=float, default=7, help='input pH')

        op = parser.parse_args()
    except argparse.ArgumentError:
        print ('ERROR - command line arguments are ill-defined, '\
                'please check the arguments')
        raise
        sys.exit(1)


    print (startTime)
       
    outputData(fileName=op.fileName, pH=op.pH)

    print(datetime.now() - startTime)


if __name__ == '__main__':
    main()
   

