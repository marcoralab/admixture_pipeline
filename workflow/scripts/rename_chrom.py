from pyfaidx import Fasta
import json

fasta = snakemake.input['fasta']
mapping_out = snakemake.output['mapping']
json_out = snakemake.output['json']

print('Starting')

seqs = Fasta(fasta).keys()
if any(x in seqs for x in ['1', 'chrUn_gl000247', 'GL000247.1']):
    mapping = {'23': 'X', '25': 'X', 'XY': 'X', '24': 'Y',
               '26': 'MT', 'M': 'MT', '0M': 'MT', 'chrM': 'MT',
               'chrMT': 'MT'}
    mapping_chr = {'chr{}'.format(x): str(x) for x in range(1, 23)}
    json_ = {'build': 37}
else:
    mapping = {'23': 'chrX', '25': 'chrX', 'XY': 'chrX', '24': 'chrY',
               'X': 'chrX', 'Y': 'chrY', '26': 'chrM', 'M': 'chrM',
               'MT': 'chrM', '0M': 'chrM', 'chrMT': 'chrM'}
    mapping_chr = {str(x): 'chr{}'.format(x) for x in range(1, 23)}
    json_ = {'build': 38,
             'map': {**mapping_chr, 'Y': 'chrY', 'X': 'chrX',
                     'MT': 'chrM', 'M': 'chrM'}}

print(f'Build {json_["build"]} detected')

mapping.update(mapping_chr)

with open(json_out, "w") as f:
    json.dump(json_, f)

with open(mapping_out, 'w') as f_out:
    for k, v in mapping.items():
        print(f'{k} {v}', file=f_out)