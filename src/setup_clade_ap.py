""" setup_clade_ap.py """
import os
import sys
from datetime import datetime
import argparse as ap
from clint.textui import colored
from conf import DI
from conf import py
import emoticons
import tree_reader

def generate_argparser ():
    """ args parser """
    parser = ap.ArgumentParser(prog="setup_clade_ap.py",
                               formatter_class=ap.ArgumentDefaultsHelpFormatter
                               )
    parser.add_argument("-t", "--taxon",
                        type=str, nargs=1, required=True, metavar=("ID/NAME"),
                        help=("The id or name of the taxon to be processes.")
                        )
    parser.add_argument("-b", "--database",
                        type=str, nargs=1, required=True,
                        help=("Location of database.")
                        )
    parser.add_argument("-s", "--seqgzfolder",
                        type=str, nargs=1, required=True,
                        help=("Location of the gzseqs directory")
                        )
    parser.add_argument("-o", "--outdir",
                        type=str, nargs=1, required=True,
                       help=("Location of the output directory (must already exist).")
                        )
    parser.add_argument("-l", "--logfile",
                        type=str, nargs=1, required=True,
                        help=("Where to write the logfile.")
                        )
    parser.add_argument("-f", "--tlistf",
                        type=str, nargs=1, required=False,
                        help=("Taxon list file.")
                        )
    return parser

if __name__ == "__main__":
    p = generate_argparser()
    args = p.parse_args(sys.argv[1:])

    print(colored.blue("STARTING PYPHLAWD "+emoticons.get_ran_emot("excited")))
    start = datetime.now()

    # Output directory (must exists!)
    dirl = args.outdir[0]
    if not os.path.isdir(dirl):
        print(f"Error: '{dirl}' does not exist")
        sys.exit(1)

    if dirl[-1] == "/":
        dirl = dirl[:-1]

    # Taxon string
    taxon = args.taxon[0]

    # Location of database file
    db = args.database[0]
    if not os.path.isfile(db):
        print(f"Error: '{db}' does not exist")
        sys.exit(1)

    # This will be used to limit the taxa
    TAXALISTF = None
    if args.tlistf is not None:
        TAXALISTF = args.tlistf[0]
        print(colored.yellow("LIMITING TO TAXA IN"), TAXALISTF)

    # Log file
    logfile = args.logfile[0]
    if logfile[-len(".md.gz"):] != ".md.gz":
        logfile += ".md.gz"

    # Folder with ?????
    gzfiles = args.seqgzfolder[0]
    if not os.path.isdir(gzfiles):
        print(f"Error: '{gzfiles}' does not exist")
        sys.exit(1)

    if gzfiles[-1] != "/":
        gzfiles += "/"

    # run get_ncbi_tax_tree_no_species.py
    print(colored.yellow("MAKING TREE"),taxon,colored.yellow(emoticons.get_ran_emot("excited")))
    tname = dirl+"/"+taxon+".tre"

    if TAXALISTF is not None:
        cmd = py+" "+DI+"get_ncbi_tax_tree_no_species.py "+taxon+" "+db+" "+TAXALISTF+" > "+tname
    else:
        cmd = py+" "+DI+"get_ncbi_tax_tree_no_species.py "+taxon+" "+db+" > "+tname

    os.system(cmd)

    trn = tree_reader.read_tree_file_iter(tname).__next__().label

    # run make_dirs.py
    print(colored.yellow("MAKING DIRS IN"), dirl, colored.yellow(emoticons.get_ran_emot("excited")))
    cmd = py+" "+DI+"make_dirs.py "+tname+" "+dirl
    os.system(cmd)

    # run populate_dirs_first.py
    print(colored.yellow("POPULATING DIRS"), dirl,colored.yellow(emoticons.get_ran_emot("excited")))
    cmd = py+" "+DI+"populate_dirs_first.py "+tname+" "+dirl+" "+db+" "+gzfiles

    if TAXALISTF is not None:
        cmd += " "+TAXALISTF

    os.system(cmd)

    if os.path.isfile("log.md.gz"):
        os.remove("log.md.gz")

    # run cluster_tree.py
    cmd = py+" "+DI+"cluster_tree.py "+dirl+"/"+trn+"/ "+logfile
    os.system(cmd)

    print(colored.blue("PYPHLAWD DONE "+emoticons.get_ran_emot("excited")))
    end = datetime.now()
    print(colored.blue("Total time (H:M:S): "+str(end-start)+" "+emoticons.get_ran_emot("excited")))
    from utils import bcolors
    print(bcolors.HEADER, end=' ')
    emoticons.animate(emoticons.glasses_animated)
    print(bcolors.ENDC)
