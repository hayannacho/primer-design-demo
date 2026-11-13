"""Primer3 design around a one-based genomic position."""
import argparse,json,random
from pathlib import Path
import primer3

def read_fasta(path, chromosome):
    records={};name=None
    for line in Path(path).read_text().splitlines():
        line=line.strip()
        if not line:continue
        if line.startswith('>'):
            name=line[1:].split()[0]
            if name in records:raise ValueError('duplicate FASTA identifier')
            records[name]=''
        elif name is None:raise ValueError('FASTA sequence has no header')
        else:records[name]+=line
    if chromosome not in records:raise ValueError('chromosome not found')
    return records[chromosome]

def design(sequence,position):
    sequence=sequence.upper()
    if not sequence or set(sequence)-set('ACGT'):raise ValueError('sequence must contain only ACGT')
    if type(position) is not int or not 1<=position<=len(sequence):raise ValueError('position outside sequence')
    start=max(1,position-500);end=min(len(sequence),position+500)
    region=sequence[start-1:end];center=position-start
    target_start=max(0,center-100);target_length=min(200,len(region)-target_start)

    result={'position':position,'window_start':start,'window_end':end,'target_start':start+target_start,'target_length':target_length,'pairs_found':0}
    if len(region)<600:
        result['explanation']={'window':'Window shorter than minimum 600 bp product'}
        return result
    raw=primer3.bindings.design_primers({'SEQUENCE_ID':'target','SEQUENCE_TEMPLATE':region,'SEQUENCE_TARGET':[target_start,target_length]}, {'PRIMER_PRODUCT_SIZE_RANGE':[[600,800]],'PRIMER_NUM_RETURN':1})
    result['pairs_found']=raw['PRIMER_PAIR_NUM_RETURNED']
    if result['pairs_found']:
        left,ll=raw['PRIMER_LEFT_0'];right,rl=raw['PRIMER_RIGHT_0']
        result.update(left_primer=raw['PRIMER_LEFT_0_SEQUENCE'],right_primer=raw['PRIMER_RIGHT_0_SEQUENCE'],left_start=start+left,right_end=start+right,product_size=raw['PRIMER_PAIR_0_PRODUCT_SIZE'])
    else:result['explanation']={k:v for k,v in raw.items() if k.endswith('_EXPLAIN')}
    return result

def synthetic_sequence():
    rng=random.Random(42)
    return ''.join(rng.choice('ACGT') for _ in range(2000))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fasta');parser.add_argument('--chromosome');parser.add_argument('--position',type=int,default=1000)
    args=parser.parse_args()
    if bool(args.fasta)!=bool(args.chromosome):parser.error('--fasta and --chromosome must be used together')
    try:result=design(read_fasta(args.fasta,args.chromosome) if args.fasta else synthetic_sequence(),args.position)
    except ValueError as e:parser.error(str(e))
    print(json.dumps(result,indent=2))
