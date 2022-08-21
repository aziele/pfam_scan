# pfam_scan

A Python script to identify protein domains in one or more protein sequences. The script uses the HMMER software to scan query protein sequences against the Pfam library of hidden Markov models (HMMs). Searching a protein sequence against the Pfam library of HMMs allows you to find out the domain architecture of the protein and thus can provide insight into protein's function.

## Requirements

* Python >= 3.8
* [HMMER software](http://hmmer.wustl.edu/) >= 3.3
   > If you can run `hmmscan -h` without getting an error, you should be good to go.

### Prepare library of Pfam HMMs

You will need to have a local copy of the Pfam's HMMs library. If you are using bash, you can follow these steps:

1. Download two files from the [Pfam FTP site](ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/):

	```
	wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
	wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
	```

2. Unpack the downloaded files to the same directory.

	```
	mkdir pfamdb
	gunzip -c Pfam-A.hmm.dat.gz > pfamdb/Pfam-A.hmm.dat
	gunzip -c Pfam-A.hmm.gz > pfamdb/Pfam-A.hmm
	rm Pfam-A.hmm.gz Pfam-A.hmm.dat.gz
	```

3. Prepare Pfam database for HMMER by creating binary files.
  
	```
	hmmpress pfamdb/Pfam-A.hmm
	```


## Installation
Since `pfam_scan.py`is just a single script, no installation is required. You can simply clone it and run it:

```
git clone https://github.com/aziele/pfam_scan
cd pfam_scan
./pfam_scan.py --help
```

## Quick usage
The script is easy to run. Just give it an input FASTA-format file containing your protein sequences and a directory location of Pfam files.:

```
./pfam_scan.py test/test.fasta pfamdb/
```

## Full usage

```
usage: pfam_scan.py [-out OUT] [-outfmt {csv,json}] [-evalue EVALUE] [-cpu CPU] [-h] [-v]
                    fasta_file pfam_dir

pfam_scan.py v.1.0: search a FASTA file against a library of Pfam HMMs

Positional arguments (required):
  fasta_file          Input FASTA file name
  pfam_dir            Directory containing Pfam data files

Optional arguments:
  -out OUT            Output file name, otherwise send to STDOUT
  -outfmt {csv,json}  Output format [default: csv]
  -evalue EVALUE      E-value threshold for a predicted domain. By default, this option is
                      disabled and the tool uses the bit score gathering (GA) threshold
                      recommended by Pfam
  -cpu CPU            Number of parallel CPU workers to use for multithreads (hmmscan)
                      [default: 16]

Other arguments:
  -h, --help          Print this help message and exit
  -v, --version       Print version information and exit
```

Positional arguments:

- `fasta_file`: input FASTA file containing protein sequences
- `pfam_dir`: directory containing Pfam data files with the following extensions: *.hmm, .h3f, .h3i, .h3m, .h3p*, and *.dat*. The file with the *.hmm* extension is a library of Pfam HMMs and it is usually named *Pfam-A.hmm*. The *.h3\** files are binaries of that library and are created by the *hmmpress* command. The data file *.dat* contains information about each Pfam HMM.

Optional arguments:

- `-out`: output file name. If not provided, the results of pfam_scan will be printed to the screen.
- `-outfmt`: output format (csv or json)
- `-evalue`: E-value domain cutoff for Pfam searches (default: Pfam defined). The E-value is the expected number of false positives (non-homologous sequence regions) that scored this well or better. The E-value is a measure of statistical significance of the domain match. The lower the E-value, the more significant the domain match. Typically, consider domains with E-values < 10<sup>−3</sup> to be significant domain matches. However, E-values are dependent on the size of the database searched, so Pfam uses its own system for maintaining Pfam models, based on a bit score, which is independent of the size of the database searched. For each Pfam family, Pfam curators set a bit score gathering (GA) threshold by hand, such that all sequences scoring at or above this threshold appear in the full alignment. The threshold is usually conservative, so that no known false positives appear in the family. Setting your own E-value threshold is possible, and may allow the presence of more distantly related domains to be found on your sequence; however, interpreting such hits requires caution.
- `-cpu`: the number of parallel worker threads used by hmmscan. By default, this is the number of CPU cores in your machine.

## Output

The output format is CSV or JSON. The output includes the following information, which should be familiar to anyone who used Pfam. Most of these values are derived from the output of hmmscan (see [HMMER3 documentation](http://eddylab.org/software/hmmer/Userguide.pdf) for details).


| # | Column | Column description |
| --- | --- | --- |
| 1 | `seq_id` | Sequence identifier of a query protein sequence |
| 2 | `aln_start` | Start position of the reported region in the protein sequence that aligns to the HMM profile |
| 3 | `aln_end` | End position |
| 4 | `env_start` | Most probable start position of the domain on the protein sequence. The envelope coordinates delineate the region on the sequence where the domain match has been probabilistically determined to lie. The envelope is almost always a little wider than the alignment coordinates. |
| 5 | `env_end` | Most probable end position of the domain on the protein sequence. |
| 6 | `hmm_acc` | Domain accession number |
| 7 | `hmm_name` | Domain name |
| 8 | `type` | Type of the predicted region (Domain / Family / Repeat / Motif) |
| 9 | `hmm_start` | Start position of the reported local alignment with respect to the HMM profile |
| 10 | `hmm_end` | End position of the reported local alignment with respect to the HMM profile |
| 11 | `hmm_length` | Length of the HMM profile |
| 12 | `score` | Score of a domain match |
| 13 | `evalue` | Evalue of a domain match |
| 14 | `significance` | Significance of a domain match (0 or 1). The significance is 1 if the score of the domain is greater than or equal to the [Pfam's gathering threshold](https://pfam-docs.readthedocs.io/en/latest/glossary.html#gathering-threshold-ga) for the matching domain, 0 otherwise. |
| 15 | `clan` | A collection of related Pfam domains. The relationship is defined by Pfam based on the similarity of sequence, tertiary structure or HMM profile. |
> Of note, when a protein sequence region has overlapping matches to more than one domain within the same clan, only the best scoring domain match is shown.


### CSV output

```csv
aln_start,aln_end,env_start,env_end,hmm_acc,hmm_name,type,hmm_start,hmm_end,hmm_length,score,evalue,significance,clan
O43347,23,90,22,91,PF00076.25,RRM_1,Domain,2,69,70,59.4,2.5e-16,1,CL0221
O43347,111,177,111,180,PF00076.25,RRM_1,Domain,1,67,70,55.1,5.6e-15,1,CL0221
E0SP36,65,204,65,215,PF00004.32,AAA,Domain,1,120,132,33.3,5.9e-08,1,CL0023
E0SP36,312,391,312,392,PF09079.14,Cdc6_C,Domain,1,83,84,83.1,1.1e-23,1,CL0123
G0ECS7,164,337,164,339,PF01170.21,UPF0020,Domain,1,195,197,96.8,1.4e-27,1,CL0063
```

### JSON output

```json
{
  "O43347": [
    {
      "aln_start": 23,
      "aln_end": 90,
      "env_start": 22,
      "env_end": 91,
      "hmm_acc": "PF00076.25",
      "hmm_name": "RRM_1",
      "type": "Domain",
      "hmm_start": 2,
      "hmm_end": 69,
      "hmm_length": 70,
      "score": 59.4,
      "evalue": 2.5e-16,
      "significance": 1,
      "clan": "CL0221"
    },
    {
      "aln_start": 111,
      "aln_end": 177,
      "env_start": 111,
      "env_end": 180,
      "hmm_acc": "PF00076.25",
      "hmm_name": "RRM_1",
      "type": "Domain",
      "hmm_start": 1,
      "hmm_end": 67,
      "hmm_length": 70,
      "score": 55.1,
      "evalue": 5.6e-15,
      "significance": 1,
      "clan": "CL0221"
    }
  ],
  "E0SP36": [
    {
      "aln_start": 65,
      "aln_end": 204,
      "env_start": 65,
      "env_end": 215,
      "hmm_acc": "PF00004.32",
      "hmm_name": "AAA",
      "type": "Domain",
      "hmm_start": 1,
      "hmm_end": 120,
      "hmm_length": 132,
      "score": 33.3,
      "evalue": 5.9e-08,
      "significance": 1,
      "clan": "CL0023"
    },
    {
      "aln_start": 312,
      "aln_end": 391,
      "env_start": 312,
      "env_end": 392,
      "hmm_acc": "PF09079.14",
      "hmm_name": "Cdc6_C",
      "type": "Domain",
      "hmm_start": 1,
      "hmm_end": 83,
      "hmm_length": 84,
      "score": 83.1,
      "evalue": 1.1e-23,
      "significance": 1,
      "clan": "CL0123"
    }
  ],
  "G0ECS7": [
    {
      "aln_start": 164,
      "aln_end": 337,
      "env_start": 164,
      "env_end": 339,
      "hmm_acc": "PF01170.21",
      "hmm_name": "UPF0020",
      "type": "Domain",
      "hmm_start": 1,
      "hmm_end": 195,
      "hmm_length": 197,
      "score": 96.8,
      "evalue": 1.4e-27,
      "significance": 1,
      "clan": "CL0063"
    }
  ]
}
```

## Test
You can run tests to ensure that the module works as expected.

```
./test.py
```

## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
