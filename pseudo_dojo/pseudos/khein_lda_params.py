# Parameteres used to generate the periodic table (A. Khein and D.C. Allan).
# See http://www.abinit.org/downloads/psp-links/psp-links/lda_tm
from __future__ import print_function

_PARAMS = {\

"H":
{'Z': 1,
 'Z_val': 1,
 'l_local': 0,
 'l_max': 0,
 'nlcc': False,
 'projectors': [{'l': 0, 'n': 1, 'rcut': 1.5855604, 'scheme': 'tm'}],
 'reference_conf': '1s1'}
,
"He":
{'Z': 2,
 'Z_val': 2,
 'l_local': 0,
 'l_max': 0,
 'nlcc': False,
 'projectors': [{'l': 0, 'n': 1, 'rcut': 1.976299, 'scheme': 'tm'}],
 'reference_conf': '1s2'}
,
"Li":
{'Z': 3,
 'Z_val': 1,
 'l_local': 1,
 'l_max': 1,
 'nlcc': {'rcore': 2.44451305764117},
 'projectors': [{'l': 0, 'n': 2, 'rcut': 2.4315963, 'scheme': 'tm'},
                {'l': 1, 'n': 2, 'rcut': 2.4315963, 'scheme': 'tm'}],
 'reference_conf': '[He] 2s1'}
,
"Be":
{'Z': 4,
 'Z_val': 2,
 'l_local': 1,
 'l_max': 1,
 'nlcc': {'rcore': 1.53889850660196},
 'projectors': [{'l': 0, 'n': 2, 'rcut': 2.2840139, 'scheme': 'tm'},
                {'l': 1, 'n': 2, 'rcut': 2.2840139, 'scheme': 'tm'}],
 'reference_conf': '[He] 2s2'}
,
"B":
{'Z': 5,
 'Z_val': 3,
 'l_local': 1,
 'l_max': 1,
 'nlcc': {'rcore': 1.10004537463277},
 'projectors': [{'l': 0, 'n': 2, 'rcut': 1.5924135, 'scheme': 'tm'},
                {'l': 1, 'n': 2, 'rcut': 1.5924135, 'scheme': 'tm'}],
 'reference_conf': '[He] 2s2 2p1'}
,
"C":
{'Z': 6,
 'Z_val': 4,
 'l_local': 1,
 'l_max': 1,
 'nlcc': {'rcore': 0.83985002509544},
 'projectors': [{'l': 0, 'n': 2, 'rcut': 1.4850707, 'scheme': 'tm'},
                {'l': 1, 'n': 2, 'rcut': 1.4850707, 'scheme': 'tm'}],
 'reference_conf': '[He] 2s2 2p2'}
,
"N":
{'Z': 7,
 'Z_val': 5,
 'l_local': 1,
 'l_max': 1,
 'nlcc': {'rcore': 0.67622446232424},
 'projectors': [{'l': 0, 'n': 2, 'rcut': 1.4975844, 'scheme': 'tm'},
                {'l': 1, 'n': 2, 'rcut': 1.4975844, 'scheme': 'tm'}],
 'reference_conf': '[He] 2s2 2p3'}
,
"O":
{'Z': 8,
 'Z_val': 6,
 'l_local': 1,
 'l_max': 1,
 'nlcc': {'rcore': 0.56990156784787},
 'projectors': [{'l': 0, 'n': 2, 'rcut': 1.4482335, 'scheme': 'tm'},
                {'l': 1, 'n': 2, 'rcut': 1.4482335, 'scheme': 'tm'}],
 'reference_conf': '[He] 2s2 2p4'}
,
"F":
{'Z': 9,
 'Z_val': 7,
 'l_local': 1,
 'l_max': 1,
 'nlcc': {'rcore': 0.4940614870118},
 'projectors': [{'l': 0, 'n': 2, 'rcut': 1.3876018, 'scheme': 'tm'},
                {'l': 1, 'n': 2, 'rcut': 1.3876018, 'scheme': 'tm'}],
 'reference_conf': '[He] 2s2 2p5'}
,
"Ne":
{'Z': 10,
 'Z_val': 8,
 'l_local': 1,
 'l_max': 1,
 'nlcc': {'rcore': 0.49764572558398},
 'projectors': [{'l': 0, 'n': 2, 'rcut': 2.4838633, 'scheme': 'tm'},
                {'l': 1, 'n': 2, 'rcut': 2.4838633, 'scheme': 'tm'}],
 'reference_conf': '[He] 2s2 2p6'}
,
"Na":
{'Z': 11,
 'Z_val': 1,
 'l_local': 2,
 'l_max': 2,
 'nlcc': {'rcore': 2.0948818708049},
 'projectors': [{'l': 0, 'n': 3, 'rcut': 2.9359409, 'scheme': 'tm'},
                {'l': 1, 'n': 3, 'rcut': 3.1646217, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 2.9359409, 'scheme': 'tm'}],
 'reference_conf': '[Ne] 3s1'}
,
"Mg":
{'Z': 12,
 'Z_val': 2,
 'l_local': 2,
 'l_max': 2,
 'nlcc': {'rcore': 2.54196289048337},
 'projectors': [{'l': 0, 'n': 3, 'rcut': 2.5922174, 'scheme': 'tm'},
                {'l': 1, 'n': 3, 'rcut': 2.5922174, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 2.5922174, 'scheme': 'tm'}],
 'reference_conf': '[Ne] 3s2'}
,
"Al":
{'Z': 13,
 'Z_val': 3,
 'l_local': 2,
 'l_max': 2,
 'nlcc': {'rcore': 2.09673076353074},
 'projectors': [{'l': 0, 'n': 3, 'rcut': 2.2761078, 'scheme': 'tm'},
                {'l': 1, 'n': 3, 'rcut': 2.2761078, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 2.2761078, 'scheme': 'tm'}],
 'reference_conf': '[Ne] 3s2 3p1'}
,
"Si":
{'Z': 14,
 'Z_val': 4,
 'l_local': 2,
 'l_max': 2,
 'nlcc': {'rcore': 1.80626423934776},
 'projectors': [{'l': 0, 'n': 3, 'rcut': 2.0872718, 'scheme': 'tm'},
                {'l': 1, 'n': 3, 'rcut': 2.0872718, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 2.0872718, 'scheme': 'tm'}],
 'reference_conf': '[Ne] 3s2 3p2'}
,
"P":
{'Z': 15,
 'Z_val': 5,
 'l_local': 2,
 'l_max': 2,
 'nlcc': {'rcore': 1.56401531320115},
 'projectors': [{'l': 0, 'n': 3, 'rcut': 1.8300797, 'scheme': 'tm'},
                {'l': 1, 'n': 3, 'rcut': 1.8300797, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 1.8300797, 'scheme': 'tm'}],
 'reference_conf': '[Ne] 3s2 3p3'}
,
"S":
{'Z': 16,
 'Z_val': 6,
 'l_local': 2,
 'l_max': 2,
 'nlcc': {'rcore': 1.3774138103074},
 'projectors': [{'l': 0, 'n': 3, 'rcut': 1.7591368, 'scheme': 'tm'},
                {'l': 1, 'n': 3, 'rcut': 1.7591368, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 1.7591368, 'scheme': 'tm'}],
 'reference_conf': '[Ne] 3s2 3p4'}
,
"Cl":
{'Z': 17,
 'Z_val': 7,
 'l_local': 2,
 'l_max': 2,
 'nlcc': {'rcore': 1.23315314129373},
 'projectors': [{'l': 0, 'n': 3, 'rcut': 1.6350894, 'scheme': 'tm'},
                {'l': 1, 'n': 3, 'rcut': 1.6350894, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 1.6350894, 'scheme': 'tm'}],
 'reference_conf': '[Ne] 3s2 3p5'}
,
"Ar":
{'Z': 18,
 'Z_val': 8,
 'l_local': 2,
 'l_max': 2,
 'nlcc': {'rcore': 1.12177160690209},
 'projectors': [{'l': 0, 'n': 3, 'rcut': 3.0330928, 'scheme': 'tm'},
                {'l': 1, 'n': 3, 'rcut': 3.0330928, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 3.0330928, 'scheme': 'tm'}],
 'reference_conf': '[Ne] 3s2 3p6'}
,
"K":
{'Z': 19,
 'Z_val': 7,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 1.21940820299569},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 3.1756741, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 2.024852, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 2.2944752, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 4s1'}
,
"Ca":
{'Z': 20,
 'Z_val': 2,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 3.75160324990621},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 2.8697489, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 3.4616035, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 3.0168904, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 4s2'}
,
"Sc":
{'Z': 21,
 'Z_val': 3,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 3.23292695886301},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 2.6324967, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 3.2153649, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 2.1285098, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 3d1 4s2'}
,
"Ti":
{'Z': 22,
 'Z_val': 4,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.82741504495369},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 2.4203472, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 2.9934301, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 1.9326638, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 3d2 4s2'}
,
"V":
{'Z': 23,
 'Z_val': 5,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.54061776212796},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 2.1748425, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 2.689797, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 1.9926227, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 3d3 4s2'}
,
"Cr":
{'Z': 24,
 'Z_val': 6,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.37464045745395},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 2.4828437, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 3.0707215, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 2.2465601, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 3d5 4s1'}
,
"Mn":
{'Z': 25,
 'Z_val': 7,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 1.96209638733297},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 2.4746132, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 2.9478926, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 2.0773155, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 3d5 4s2'}
,
"Fe":
{'Z': 26,
 'Z_val': 8,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 1.56404770202776},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 2.2918558, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 2.8345121, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 2.2918558, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 3d6 4s2'}
,
"Co":
{'Z': 27,
 'Z_val': 9,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 1.39728307798009},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 2.3788675, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 2.9421245, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 2.0474975, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 3d7 4s2'}
,
"Ni":
{'Z': 28,
 'Z_val': 10,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 1.28166113122492},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 2.9454613, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 2.9454613, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 2.293908, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 3d8 4s2'}
,
"Cu":
{'Z': 29,
 'Z_val': 11,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 1.13377888705965},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 2.0547672, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 2.2994432, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 2.0547672, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 3d10 4s1'}
,
"Zn":
{'Z': 30,
 'Z_val': 12,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 1.1521850289427},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 2.5825328, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 2.5825328, 'scheme': 'tm'},
                {'l': 2, 'n': 3, 'rcut': 2.395922, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 3d10 4s2'}
,
"Ga":
{'Z': 31,
 'Z_val': 3,
 'l_local': 0,
 'l_max': 1,
 'nlcc': {'rcore': 2.5764986203335},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 1.9956558, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 1.9956558, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 3d10 4s2 4p1'}
,
"Ge":
{'Z': 32,
 'Z_val': 4,
 'l_local': 1,
 'l_max': 1,
 'nlcc': {'rcore': 2.25844779088698},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 1.982235, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 1.982235, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 3d10 4s2 4p2'}
,
"As":
{'Z': 33,
 'Z_val': 5,
 'l_local': 1,
 'l_max': 1,
 'nlcc': {'rcore': 2.0573171556401},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 2.530616, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 2.530616, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 3d10 4s2 4p3'}
,
"Se":
{'Z': 34,
 'Z_val': 6,
 'l_local': 1,
 'l_max': 1,
 'nlcc': {'rcore': 1.85251755370876},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 1.8424568, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 1.8424568, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 3d10 4s2 4p4'}
,
"Br":
{'Z': 35,
 'Z_val': 7,
 'l_local': 1,
 'l_max': 1,
 'nlcc': {'rcore': 1.71181633387051},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 1.8351264, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 1.8351264, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 3d10 4s2 4p5'}
,
"Kr":
{'Z': 36,
 'Z_val': 8,
 'l_local': 1,
 'l_max': 1,
 'nlcc': {'rcore': 1.58309363830263},
 'projectors': [{'l': 0, 'n': 4, 'rcut': 2.9786136, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 3.2106057, 'scheme': 'tm'}],
 'reference_conf': '[Ar] 3d10 4s2 4p6'}
,
"Rb":
{'Z': 37,
 'Z_val': 7,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 1.76736606338703},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 3.539775, 'scheme': 'tm'},
                {'l': 1, 'n': 4, 'rcut': 2.4944175, 'scheme': 'tm'},
                {'l': 2, 'n': 4, 'rcut': 2.5897367, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 5s1'}
,
"Sr":
{'Z': 38,
 'Z_val': 8,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 1.59650627399569},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 3.489977, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 2.2816189, 'scheme': 'tm'},
                {'l': 2, 'n': 4, 'rcut': 2.5854217, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 5s2'}
,
"Y":
{'Z': 39,
 'Z_val': 3,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 3.77869318742782},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 2.890462, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 3.5304324, 'scheme': 'tm'},
                {'l': 2, 'n': 4, 'rcut': 2.890462, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 4d1 5s2'}
,
"Zr":
{'Z': 40,
 'Z_val': 4,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 3.37554947056294},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 2.8182005, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 3.3994115, 'scheme': 'tm'},
                {'l': 2, 'n': 4, 'rcut': 2.5183299, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 4d2 5s2'}
,
"Nb":
{'Z': 41,
 'Z_val': 5,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 3.05525591962919},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 2.890435, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 3.316499, 'scheme': 'tm'},
                {'l': 2, 'n': 4, 'rcut': 2.2793753, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 4d4 5s1'}
,
"Mo":
{'Z': 42,
 'Z_val': 6,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.63204696240783},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 2.8930462, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 3.1973167, 'scheme': 'tm'},
                {'l': 2, 'n': 4, 'rcut': 2.197463, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 4d5 5s1'}
,
"Tc":
{'Z': 43,
 'Z_val': 7,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.41507221543774},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 2.687949, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 2.9706491, 'scheme': 'tm'},
                {'l': 2, 'n': 4, 'rcut': 2.1733578, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 4d5 5s2'}
,
"Ru":
{'Z': 44,
 'Z_val': 8,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.13557494563996},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 2.6933599, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 2.976629, 'scheme': 'tm'},
                {'l': 2, 'n': 4, 'rcut': 2.1239633, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 4d7 5s1'}
,
"Rh":
{'Z': 45,
 'Z_val': 9,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 1.91316594310899},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 2.5684846, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 2.9841622, 'scheme': 'tm'},
                {'l': 2, 'n': 4, 'rcut': 2.183245, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 4d8 5s1'}
,
"Pd":
{'Z': 46,
 'Z_val': 10,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 1.6934637418257},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 2.4201665, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 3.1862405, 'scheme': 'tm'},
                {'l': 2, 'n': 4, 'rcut': 2.2452898, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 4d10'}
,
"Ag":
{'Z': 47,
 'Z_val': 11,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 1.57659479851128},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 2.3392488, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 2.9663571, 'scheme': 'tm'},
                {'l': 2, 'n': 4, 'rcut': 2.3392488, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 4d10 5s1'}
,
"Cd":
{'Z': 48,
 'Z_val': 12,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 1.48692763900022},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 2.3780416, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 2.8684763, 'scheme': 'tm'},
                {'l': 2, 'n': 4, 'rcut': 2.3780416, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 4d10 5s2'}
,
"In":
{'Z': 49,
 'Z_val': 3,
 'l_local': 0,
 'l_max': 1,
 'nlcc': {'rcore': 3.24177624483359},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 2.1883691, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 2.6396875, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 4d10 5s2 5p1'}
,
"Sn":
{'Z': 50,
 'Z_val': 4,
 'l_local': 0,
 'l_max': 1,
 'nlcc': {'rcore': 2.9107668115427},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 2.2829199, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 2.2829199, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 4d10 5s2 5p2'}
,
"Sb":
{'Z': 51,
 'Z_val': 5,
 'l_local': 0,
 'l_max': 1,
 'nlcc': {'rcore': 2.64748978083368},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 2.0764318, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 2.3236831, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 4d10 5s2 5p3'}
,
"Te":
{'Z': 52,
 'Z_val': 6,
 'l_local': 0,
 'l_max': 1,
 'nlcc': {'rcore': 2.43925356543835},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 2.0880559, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 2.2789968, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 4d10 5s2 5p4'}
,
"I":
{'Z': 53,
 'Z_val': 7,
 'l_local': 0,
 'l_max': 1,
 'nlcc': {'rcore': 2.24822719245523},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 1.9980759, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 2.2359969, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 4d10 5s2 5p5'}
,
"Xe":
{'Z': 54,
 'Z_val': 8,
 'l_local': 0,
 'l_max': 1,
 'nlcc': {'rcore': 2.12537586803331},
 'projectors': [{'l': 0, 'n': 5, 'rcut': 1.9857424, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 1.9857424, 'scheme': 'tm'}],
 'reference_conf': '[Kr] 4d10 5s2 5p6'}
,
"Cs":
{'Z': 55,
 'Z_val': 7,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.51709401117676},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 3.5971947, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 2.4722992, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.7666869, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 6s1'}
,
"Ba":
{'Z': 56,
 'Z_val': 8,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.29351246205421},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 3.3606525, 'scheme': 'tm'},
                {'l': 1, 'n': 5, 'rcut': 2.4896211, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.7860712, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 6s2'}
,
"La":
{'Z': 57,
 'Z_val': 3,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 4.37059395158683},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 3.4278596, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.9331384, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.7371927, 'scheme': 'tm'},
                {'l': 3, 'n': 4, 'rcut': 3.2201735, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 5d1 6s2'}
,
"Ce":
{'Z': 58,
 'Z_val': 4,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 4.29523888345602},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 3.3687586, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.8173093, 'scheme': 'tm'},
                {'l': 2, 'n': 4, 'rcut': 2.6899998, 'scheme': 'tm'},
                {'l': 3, 'n': 5, 'rcut': 2.972914, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f1 5d1 6s2'}
,
"Pr":
{'Z': 59,
 'Z_val': 5,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 4.49476615739633},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 3.2298947, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.8960054, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.7799903, 'scheme': 'tm'},
                {'l': 3, 'n': 4, 'rcut': 2.9225256, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f3 6s2'}
,
"Nd":
{'Z': 60,
 'Z_val': 6,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 4.41985338810639},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 3.1760632, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.7834811, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.7680428, 'scheme': 'tm'},
                {'l': 3, 'n': 4, 'rcut': 2.9465689, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f4 6s2'}
,
"Pm":
{'Z': 61,
 'Z_val': 7,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 4.29339178794363},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 3.1239966, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.7214568, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.7226651, 'scheme': 'tm'},
                {'l': 3, 'n': 4, 'rcut': 2.8982645, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f5 6s2'}
,
"Sm":
{'Z': 62,
 'Z_val': 8,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 4.22414353329938},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 3.0736095, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.6614333, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.7465651, 'scheme': 'tm'},
                {'l': 3, 'n': 4, 'rcut': 2.8873865, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f6 6s2'}
,
"Eu":
{'Z': 63,
 'Z_val': 9,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 4.15709363594542},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.9872466, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.6945347, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.6693915, 'scheme': 'tm'},
                {'l': 3, 'n': 4, 'rcut': 2.841555, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f7 6s2'}
,
"Gd":
{'Z': 64,
 'Z_val': 10,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.74929015376274},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.9405709, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.4594365, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.6942033, 'scheme': 'tm'},
                {'l': 3, 'n': 5, 'rcut': 2.83234, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f7 5d1 6s2'}
,
"Tb":
{'Z': 65,
 'Z_val': 11,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.97913103610795},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.9317506, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.5363741, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.686122, 'scheme': 'tm'},
                {'l': 3, 'n': 4, 'rcut': 2.8238444, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f9 6s2'}
,
"Dy":
{'Z': 66,
 'Z_val': 12,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.91884117192449},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.9236488, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.4827927, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.678699, 'scheme': 'tm'},
                {'l': 3, 'n': 4, 'rcut': 2.8160408, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f10 6s2'}
,
"Ho":
{'Z': 67,
 'Z_val': 13,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.86035100517935},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.9162388, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.4739655, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.6719099, 'scheme': 'tm'},
                {'l': 3, 'n': 4, 'rcut': 2.7740103, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f11 6s2'}
,
"Er":
{'Z': 68,
 'Z_val': 14,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.75633161450348},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.7675961, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.4659327, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.5999138, 'scheme': 'tm'},
                {'l': 3, 'n': 4, 'rcut': 2.7675961, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f12 6s2'}
,
"Tm":
{'Z': 69,
 'Z_val': 15,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.70189202588749},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.7274861, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.4157018, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.5944632, 'scheme': 'tm'},
                {'l': 3, 'n': 4, 'rcut': 2.7965336, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f13 6s2'}
,
"Yb":
{'Z': 70,
 'Z_val': 16,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.6490078540891},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.688522, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.409257, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.4942562, 'scheme': 'tm'},
                {'l': 3, 'n': 4, 'rcut': 2.7912571, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f14 6s2'}
,
"Lu":
{'Z': 71,
 'Z_val': 17,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.3376585931622},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.6177281, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.2375256, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.3984089, 'scheme': 'tm'},
                {'l': 3, 'n': 4, 'rcut': 2.683997, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f14 5d1 6s2'}
,
"Hf":
{'Z': 72,
 'Z_val': 4,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 3.1701620283436},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.6467193, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.2733809, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.4863606, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f14 5d2 6s2'}
,
"Ta":
{'Z': 73,
 'Z_val': 5,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.90080448978014},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.4831474, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.2285401, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.3326991, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f14 5d3 6s2'}
,
"W":
{'Z': 74,
 'Z_val': 6,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.65483123902653},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.3891099, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.1849111, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.2443591, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f14 5d4 6s2'}
,
"Re":
{'Z': 75,
 'Z_val': 7,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.43015877533117},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.2990535, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.0648576, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.1869257, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f14 5d5 6s2'}
,
"Os":
{'Z': 76,
 'Z_val': 8,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.25288147499031},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.2406189, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 2.9869588, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.0787172, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f14 5d6 6s2'}
,
"Ir":
{'Z': 77,
 'Z_val': 9,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.08889782480262},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.1301224, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 2.8753758, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 2.0262336, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f14 5d7 6s2'}
,
"Pt":
{'Z': 78,
 'Z_val': 10,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 1.84269728317135},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.1292637, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 2.8742165, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 1.9754083, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f14 5d9 6s1'}
,
"Au":
{'Z': 79,
 'Z_val': 11,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 1.84225747565994},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.0761953, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 3.0973349, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 1.9261744, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f14 5d10 6s1'}
,
"Hg":
{'Z': 80,
 'Z_val': 12,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 1.75226975450466},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 2.2377245, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 2.9460428, 'scheme': 'tm'},
                {'l': 2, 'n': 5, 'rcut': 1.8784686, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f14 5d10 6s2'}
,
"Tl":
{'Z': 81,
 'Z_val': 3,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 4.20386791783521},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 1.8786145, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 2.4732633, 'scheme': 'tm'},
                {'l': 2, 'n': 6, 'rcut': 4.1290338, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f14 5d10 6s2 6p1'}
,
"Pb":
{'Z': 82,
 'Z_val': 4,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 3.2747188384661},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 1.9266159, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 2.3827805, 'scheme': 'tm'},
                {'l': 2, 'n': 6, 'rcut': 3.9285604, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f14 5d10 6s2 6p2'}
,
"Bi":
{'Z': 83,
 'Z_val': 5,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.96420654001913},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 1.8333467, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 2.2959493, 'scheme': 'tm'},
                {'l': 2, 'n': 6, 'rcut': 3.5560522, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f14 5d10 6s2 6p3'}
,
"Po":
{'Z': 84,
 'Z_val': 6,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.75146149201665},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 1.744846, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 2.1851177, 'scheme': 'tm'},
                {'l': 2, 'n': 6, 'rcut': 3.3423507, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f14 5d10 6s2 6p4'}
,
"At":
{'Z': 85,
 'Z_val': 7,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.58647757315343},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 1.746008, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 2.1325854, 'scheme': 'tm'},
                {'l': 2, 'n': 6, 'rcut': 3.0643615, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f14 5d10 6s2 6p5'}
,
"Rn":
{'Z': 86,
 'Z_val': 8,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.43172293492207},
 'projectors': [{'l': 0, 'n': 6, 'rcut': 1.6830971, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 2.0816042, 'scheme': 'tm'},
                {'l': 2, 'n': 6, 'rcut': 2.9172542, 'scheme': 'tm'}],
 'reference_conf': '[Xe] 4f14 5d10 6s2 6p6'}
,
"Fr":
{'Z': 87,
 'Z_val': 7,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.82792118185733},
 'projectors': [{'l': 0, 'n': 7, 'rcut': 3.8926282, 'scheme': 'tm'},
                {'l': 1, 'n': 7, 'rcut': 2.7775844, 'scheme': 'tm'},
                {'l': 2, 'n': 7, 'rcut': 3.4784392, 'scheme': 'tm'}],
 'reference_conf': '[Rn] 7s1'}
,
"Ra":
{'Z': 88,
 'Z_val': 8,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.69288363766496},
 'projectors': [{'l': 0, 'n': 7, 'rcut': 3.4389115, 'scheme': 'tm'},
                {'l': 1, 'n': 6, 'rcut': 2.547601, 'scheme': 'tm'},
                {'l': 2, 'n': 6, 'rcut': 2.9971259, 'scheme': 'tm'}],
 'reference_conf': '[Rn] 7s2'}
,
"Ac":
{'Z': 89,
 'Z_val': 3,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 4.22837965059553},
 'projectors': [{'l': 0, 'n': 7, 'rcut': 3.6195725, 'scheme': 'tm'},
                {'l': 1, 'n': 7, 'rcut': 4.4765713, 'scheme': 'tm'},
                {'l': 2, 'n': 6, 'rcut': 3.0766907, 'scheme': 'tm'}],
 'reference_conf': '[Rn] 6d1 7s2'}
,
"Th":
{'Z': 90,
 'Z_val': 4,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 4.2872512858933},
 'projectors': [{'l': 0, 'n': 7, 'rcut': 3.0807756, 'scheme': 'tm'},
                {'l': 1, 'n': 7, 'rcut': 4.0559443, 'scheme': 'tm'},
                {'l': 2, 'n': 6, 'rcut': 2.9305231, 'scheme': 'tm'},
                {'l': 3, 'n': 5, 'rcut': 3.0807756, 'scheme': 'tm'}],
 'reference_conf': '[Rn] 6d2 7s2'}
,
"Pa":
{'Z': 91,
 'Z_val': 5,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.98323914580088},
 'projectors': [{'l': 0, 'n': 7, 'rcut': 2.9716915, 'scheme': 'tm'},
                {'l': 1, 'n': 7, 'rcut': 3.7215237, 'scheme': 'tm'},
                {'l': 2, 'n': 6, 'rcut': 2.8983195, 'scheme': 'tm'},
                {'l': 3, 'n': 5, 'rcut': 3.0852468, 'scheme': 'tm'}],
 'reference_conf': '[Rn] 5f2 6d1 7s2'}
,
"U":
{'Z': 92,
 'Z_val': 6,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.70123153635492},
 'projectors': [{'l': 0, 'n': 7, 'rcut': 2.8312036, 'scheme': 'tm'},
                {'l': 1, 'n': 7, 'rcut': 3.6810723, 'scheme': 'tm'},
                {'l': 2, 'n': 6, 'rcut': 2.5939999, 'scheme': 'tm'},
                {'l': 3, 'n': 5, 'rcut': 2.8312036, 'scheme': 'tm'}],
 'reference_conf': '[Rn] 5f3 6d1 7s2'}
,
"93":
{'Z': 93,
 'Z_val': 7,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.61594979562936},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.8716628, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.463892, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.5342304, 'scheme': 'tm'},
                {'l': 3, 'n': None, 'rcut': 2.6310694, 'scheme': 'tm'}],
 'reference_conf': None}
,
"94":
{'Z': 94,
 'Z_val': 8,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.53304155508527},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.7365434, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.3844704, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.4761243, 'scheme': 'tm'},
                {'l': 3, 'n': None, 'rcut': 2.6358223, 'scheme': 'tm'}],
 'reference_conf': None}
,
"95":
{'Z': 95,
 'Z_val': 9,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.45242499123138},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.6408827, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.390968, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.4196244, 'scheme': 'tm'},
                {'l': 3, 'n': None, 'rcut': 2.6408827, 'scheme': 'tm'}],
 'reference_conf': None}
,
"96":
{'Z': 96,
 'Z_val': 10,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.1695977824795},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.4859165, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.1919874, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.3353009, 'scheme': 'tm'},
                {'l': 3, 'n': None, 'rcut': 2.6462461, 'scheme': 'tm'}],
 'reference_conf': None}
,
"97":
{'Z': 97,
 'Z_val': 11,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.13692151668075},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.4912355, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.1198372, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.2825148, 'scheme': 'tm'},
                {'l': 3, 'n': None, 'rcut': 2.5864315, 'scheme': 'tm'}],
 'reference_conf': None}
,
"98":
{'Z': 98,
 'Z_val': 12,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.10491211344931},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.4968312, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.0880022, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.231159, 'scheme': 'tm'},
                {'l': 3, 'n': None, 'rcut': 2.5282377, 'scheme': 'tm'}],
 'reference_conf': None}
,
"99":
{'Z': 99,
 'Z_val': 13,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.03536865488811},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.4716106, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.0952605, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.1811858, 'scheme': 'tm'},
                {'l': 3, 'n': None, 'rcut': 2.4716106, 'scheme': 'tm'}],
 'reference_conf': None}
,
"100":
{'Z': 100,
 'Z_val': 14,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.6247476226788},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.8786493, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.4291849, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.1865358, 'scheme': 'tm'},
                {'l': 3, 'n': None, 'rcut': 2.477673, 'scheme': 'tm'}],
 'reference_conf': None}
,
"101":
{'Z': 101,
 'Z_val': 15,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.54427705215569},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.8147424, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.3952325, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.1921182, 'scheme': 'tm'},
                {'l': 3, 'n': None, 'rcut': 2.4839986, 'scheme': 'tm'}],
 'reference_conf': None}
,
"102":
{'Z': 102,
 'Z_val': 16,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.50952923791887},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.7183314, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.3201829, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.1979303, 'scheme': 'tm'},
                {'l': 3, 'n': None, 'rcut': 2.4905845, 'scheme': 'tm'}],
 'reference_conf': None}
,
"103":
{'Z': 103,
 'Z_val': 17,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 3.22432918619433},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.6919399, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.2879481, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.1765911, 'scheme': 'tm'},
                {'l': 3, 'n': None, 'rcut': 2.4974279, 'scheme': 'tm'}],
 'reference_conf': None}
,
"104":
{'Z': 104,
 'Z_val': 18,
 'l_local': 0,
 'l_max': 3,
 'nlcc': {'rcore': 2.92578256844262},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.6329373, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.2158821, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.1827776, 'scheme': 'tm'},
                {'l': 3, 'n': None, 'rcut': 2.4734142, 'scheme': 'tm'}],
 'reference_conf': None}
,
"105":
{'Z': 105,
 'Z_val': 5,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 3.16291339217158},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.4806736, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.5202533, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.189184, 'scheme': 'tm'}],
 'reference_conf': None}
,
"106":
{'Z': 106,
 'Z_val': 6,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.90668700536708},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.48818, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.4437263, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.1958084, 'scheme': 'tm'}],
 'reference_conf': None}
,
"107":
{'Z': 107,
 'Z_val': 7,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.70505820018047},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.4959312, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.2859773, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.1752868, 'scheme': 'tm'}],
 'reference_conf': None}
,
"108":
{'Z': 108,
 'Z_val': 8,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.5176355989302},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.4728207, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.1751712, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.1822539, 'scheme': 'tm'}],
 'reference_conf': None}
,
"109":
{'Z': 109,
 'Z_val': 9,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.37287630593593},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.4809535, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.0683646, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.1894311, 'scheme': 'tm'}],
 'reference_conf': None}
,
"110":
{'Z': 110,
 'Z_val': 10,
 'l_local': 0,
 'l_max': 2,
 'nlcc': {'rcore': 2.23662857110534},
 'projectors': [{'l': 0, 'n': None, 'rcut': 2.4893225, 'scheme': 'tm'},
                {'l': 1, 'n': None, 'rcut': 3.078715, 'scheme': 'tm'},
                {'l': 2, 'n': None, 'rcut': 2.1968167, 'scheme': 'tm'}],
 'reference_conf': None}
,
}

def get_params(symbol):
    try:
        return _PARAMS[symbol]
    except KeyError:
        return None

def ape_pseudize(symbol, wave_equation, xc_functional):
    from pseudo_dojo.ppcodes.ape import ApeAtomicConfiguration, ApeAeSolver
    params = get_params(symbol)
    ae_conf = ApeAtomicConfiguration.from_string(params["Z"], params["reference_conf"])
    print(ae_conf.to_input())
    #ae_solver = ApeAeSolver.from_aeconf("hello_ape", ae_conf, verbose=1)
    #ae_solver.solve()
    #inpgen = ApeInputGenerator(template=None, newvars=None)
    #pp_components = ApePPComponents.from_string()
    #setup = ApePPSetup(pp_components, core_correction=0, llocal=-1)

if __name__ == "__main__":
    ape_pseudize("Si", "None", "None")
    #from  pseudo_dojo.ppcodes.nist import nist_database
    #from pprint import pprint
    #for k, d in _PARAMS.items():
    #    try:
    #        symbol = nist_database.symbol_from_Z(d["Z"])
    #    except KeyError:
    #        symbol = None
    #    d["symbol"] = symbol
    #pprint(_PARAMS)

