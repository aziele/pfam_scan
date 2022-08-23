#!/usr/bin/env python3
"""A script to scan one or more protein sequences against the Pfam HMMs library.

Copyright 2022 Andrzej Zielezinski (a.zielezinski@gmail.com)
https://github.com/aziele/pfam_scan
"""

from __future__ import annotations
import argparse
import collections
import csv
import json
import multiprocessing
import pathlib
import subprocess
import sys
import tempfile
import typing
import uuid

__version__ = '1.0'


DOMAIN_INFO = [
    'aln_start',
    'aln_end',
    'env_start',
    'env_end',
    'hmm_acc', 
    'hmm_name',
    'type',
    'hmm_start',
    'hmm_end',
    'hmm_length',
    'score',
    'evalue',
    'significance',
    'clan',
]
Domain = collections.namedtuple('Domain', DOMAIN_INFO)


def get_parser() -> argparse.ArgumentParser:
    """Returns an argument parser."""
    parser = argparse.ArgumentParser(description=f'pfam_scan.py v.{__version__}'
    ': search a FASTA file against a library of Pfam HMMs', add_help=False)
    # Required arguments
    p = parser.add_argument_group('Positional arguments (required)')
    p.add_argument('fasta_file',
                   help='Input FASTA file name')
    p.add_argument('pfam_dir', 
                   help='Directory containing Pfam data files')
    # Optional arguments
    p = parser.add_argument_group('Optional arguments')
    p.add_argument('-out', type=argparse.FileType('w'), default=sys.stdout,
                   help='Output file name, otherwise send to STDOUT')
    p.add_argument('-outfmt', choices=['csv', 'json'], default='csv',
                   help='Output format [default: %(default)s]')
    p.add_argument('-evalue', type=float,
                   help='E-value threshold for a predicted domain. By default, '
                   'this option is disabled and the tool uses the bit score '
                   'gathering (GA) threshold recommended by Pfam')
    p.add_argument('-cpu', type=int, 
                   default=min(multiprocessing.cpu_count(), 16),
                   help='Number of parallel CPU workers to use for multithreads'
                   ' (hmmscan) [default: %(default)s]')
    # Other arguments
    p = parser.add_argument_group('Other arguments')
    p.add_argument('-h', '--help', action='help',
                   default=argparse.SUPPRESS,
                   help='Print this help message and exit')
    p.add_argument('-v', '--version', action='version',
                   version=f'{__version__}',
                   help="Print version information and exit")
    # Display help if the script is run without arguments.
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        #parser.print_usage()
        parser.exit()
    return parser


def validate_args(parser: argparse.ArgumentParser) -> argparse.Namespace:
    """Validates arguments provided by the user.

    Returns:
        Arguments provided by the users.

    Raises:
        argparse.ArgumentParser.error if arguments are invalid.
    """
    args = parser.parse_args()
    # Check if the input FASTA file exists.
    fasta_path = pathlib.Path(args.fasta_file)
    if not fasta_path.exists():
        parser.error(f'The input file does not exist: {fasta_path}')

    # Check if the Pfam database directory exists.
    dir_path = pathlib.Path(args.pfam_dir)
    if not dir_path.exists():
        parser.error(f'The Pfam directory does not exist: {dir_path}')

    # Find the Pfam HMM database in the Pfam directory.
    dir_files = list(dir_path.iterdir())
    db_paths = [f for f in dir_files if f.suffix == '.hmm']
    if not db_paths:
        parser.error(f'Cannot find the HMM library file in: {dir_path}')
    if (n := len(db_paths)) > 1:
        parser.error(f'Found {n} HMM files in {dir_path}, expect one.')
    db_path = db_paths[0]

    # Make sure we have all binaries for the HMM database.
    for suffix in ['.h3f', '.h3i', '.h3m', '.h3p', '.dat']:
        path = dir_path / (db_path.name + suffix)
        if suffix == '.dat': args.dat = path
        if not path.exists():
            parser.error(f'Cannot find: {path}')
    args.fasta = fasta_path
    args.db = db_path
    # Create the name of the temporary file for storing the hmmscan output.
    args.temp_file = pathlib.Path(tempfile.gettempdir(), str(uuid.uuid4().hex))
    return args


def read_pfam_data(
        filename: Union[str, pathlib.Path]
        ) -> dict[str, collections.namedtuple]:
    """Reads the Pfam data file to dictionary.

    Args:
        filename: Name/Path of the Pfam data file (Pfam-A.hmm.dat).

    Returns:
        A dict mapping HMM profile name to the corresponding information.
        For example:

        {'1-cysPrx_C': Data(type='Domain', clan=None, ga_seq=21.1, ga_dom=21.1),
         'RRM': Data(type='Domain', clan=None, ga_seq=21.0, ga_dom=21.0),
         'SOXp': Data(type='Family', clan=None, ga_seq=22.1, ga_dom=22.1)}
    """
    data = {}
    Data = collections.namedtuple('Data', ['type', 'clan', 'ga_seq', 'ga_dom'])
    with open(filename) as fh:
        clan = None   # Not all domains have clan assigned.
        for line in fh:
            if line.startswith('#=GF ID'):
                hmm_name = line[10:-1]
            elif line.startswith('#=GF TP'):
                typ = line[10:-1]
            elif line.startswith('#=GF CL'):
                clan = line[10:-1]
            elif line.startswith('#=GF GA'):
                scores = line[10:-1].strip().rstrip(';').split(';')
                ga_seq = float(scores[0])
                ga_dom = float(scores[1])
            elif line.startswith('//'):
                data[hmm_name] = Data(typ, clan, ga_seq, ga_dom)
                clan = None
    return data


def run_hmmscan(args) -> subprocess.CompletedProcess:
    """Runs hmmscan search using the supplied arguments."""
    filtering = ['--cut_ga']
    if args.evalue:
        #filtering = ['-E', f'{args.evalue}', '--domE', f'{args.evalue}']
        filtering = ['--domE', f'{args.evalue}']
    cmd = [
        'hmmscan',
        '--notextw',
        '--cpu',
        f'{args.cpu}',
        *filtering,
        '--domtblout',
        f'{args.temp_file}',
        f'{args.db}',
        f'{args.fasta}',
    ]
    return subprocess.run(cmd, stdout=subprocess.DEVNULL, 
                          stderr=subprocess.PIPE, text=True)

def parse_hmmscan_output(
        filename: Union[str, pathlib.Path],
        pfam_data: dict[str, collections.namedtuple]
        ) -> dict[str, list[collections.namedtuple]]:
    """Parses hmmscan output and returns a dictionary with the domain results.

    Args:
        filename: Name/Path of the query FASTA file.
        pfam_data: A dict containing information on HMM profiles.

    Returns:
        A dict mapping protein ids to the corresponding list of domains.
        For example:

        {'E0SP36': [
            Domain(
                aln_start=65, aln_end=204, env_start=65, env_end=215,
                hmm_acc='PF00004.32', hmm_name='AAA', type='Domain',
                hmm_start=1, hmm_end=120, hmm_length=132, score=33.3,
                evalue=5.9e-08, significance=1, clan='CL0023'
            ),
            Domain(
                aln_start=312, aln_end=391, env_start=312, env_end=392,
                hmm_acc='PF09079.14', hmm_name='Cdc6_C', type='Domain',
                hmm_start=1, hmm_end=83, hmm_length=84, score=83.1,
                evalue=1.1e-23, significance=1, clan='CL0123'
            )
        ], 
        'G0ECS7': [
            Domain(
                aln_start=164, aln_end=337, env_start=164, env_end=339,
                hmm_acc='PF01170.21', hmm_name='UPF0020', type='Domain',
                hmm_start=1, hmm_end=195, hmm_length=197, score=96.8,
                evalue=1.4e-27, significance=1, clan='CL0063'
            )
        ]}
    """
    results = {}
    with open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.split()
            hmm_name = cols[0]
            seq_id = cols[3]
            score_dom = float(cols[13])
            score_seq = float(cols[7])
            # Determine which domain hits are significant. The significance 
            # value is 1 if the bit scores for a domain and a sequence are 
            # greater than or equal to the curated gathering thresholds for 
            # the matching domain, 0 otherwise. 
            significance = 0
            if (score_dom >= pfam_data[hmm_name].ga_dom and 
                score_seq >= pfam_data[hmm_name].ga_seq):
                # Since both conditions are true, the domain hit is significant.
                significance = 1
            dom = Domain(
                int(cols[17]),             # aln_start
                int(cols[18]),             # aln_end
                int(cols[19]),             # env_start
                int(cols[20]),             # env_end
                cols[1],                   # hmm_acc
                hmm_name,                  # hmm_name
                pfam_data[hmm_name].type,  # type
                int(cols[15]),             # hmm_start
                int(cols[16]),             # hmm_end
                int(cols[2]),              # hmm_length
                score_dom,                 # score
                float(cols[12]),           # evalue
                significance,              # significance
                pfam_data[hmm_name].clan,  # clan
            )
            if seq_id not in results:
                results[seq_id] = []
            results[seq_id].append(dom)
    # For each protein, sort domains by their start position.
    for seq_id in results:
        results[seq_id].sort(key=lambda x: x.aln_start)
    return results


def resolve_overlapping_domains(results):
    """Resolves overlapping domains belonging to the same clan.

    When a protein sequence region has overlapping matches to more than one
    domains within the same clan, only the best scoring domain match is shown.

    Args:
        results: a dict with the domain results
    """
    # TODO: Consider resolving overlapping domains belonging to different clans.
    for seq_id, domains in results.items():
        is_overlap = True
        while is_overlap:
            is_overlap = False
            indexes = set()
            for i in range(len(domains)-1):
                j = i + 1
                di = domains[i]  # Domain i (previous domain)
                dj = domains[j]  # Domain j (next domain)
                if di.clan:
                    # Overlapping domains within the same clan
                    if di.clan == dj.clan and dj.aln_start <= di.aln_end:
                        idx = j if di.evalue < dj.evalue else i
                        indexes.add(idx)
                        is_overlap = True
            # Filter out overlapping domains that have weak scores.
            domains = [d for i, d in enumerate(domains) if i not in indexes]
        results[seq_id] = domains
    return results


def main():
    parser = get_parser()
    args = validate_args(parser)
    pfam_data = read_pfam_data(args.dat)
    
    # Run hmmscan
    process = run_hmmscan(args)
    if process.returncode:
        print(f'Error running hmmscan: {process.stderr.strip()}')
        sys.exit(1)
    
    # Get results
    results = parse_hmmscan_output(args.temp_file, pfam_data)
    results = resolve_overlapping_domains(results)
    
    # Write/show results in csv or json format.
    if args.outfmt == 'csv':
        csv_out = csv.writer(args.out)
        csv_out.writerow(['seq_id', *DOMAIN_INFO])
        for seq_id in results:
            for domain in results[seq_id]:
                if not domain.clan: domain = domain._replace(clan='No_Clan')
                domain_list = list(domain)
                domain_list.insert(0, seq_id)
                csv_out.writerow(domain_list)
    elif args.outfmt == 'json':
        for seq_id, domains in results.items():
            results[seq_id] = [domain._asdict() for domain in domains]
        json.dump(results, args.out, indent=2)
    
    # Remove a temporary file containing hmmscan output.
    args.temp_file.unlink()


if __name__ == '__main__':
    main()