""" setup_clade_ap.py """
import os
import sys
from datetime import datetime
import argparse as ap
from clint.textui import colored
from utils import bcolors
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
    return parser.parse_args()

if __name__ == "__main__":
    args = generate_argparser()

    print(colored.blue("STARTING PYPHLAWD "+emoticons.get_ran_emot("excited")))
    start = datetime.now()

    # Output directory (must exists!)
    dirl = args.outdir
    dirl.rstrip('/')
    if not os.path.isdir(dirl):
        print(f"Error: '{dirl}' does not exist")
        sys.exit(1)

    # Taxon string
    taxon = args.taxon

    # Location of database file
    db = args.database
    if not os.path.isfile(db):
        print(f"Error: '{db}' does not exist")
        sys.exit(1)

    # This will be used to limit the taxa
    TAXALISTF = None
    if args.tlistf is not None:
        TAXALISTF = args.tlistf
        print(colored.yellow("LIMITING TO TAXA IN"), TAXALISTF)

    # Log file
    logfile = args.logfile
    if not logfile.endswith(".md.gz"):
        logfile += ".md.gz"

    # Folder with ?????
    gzfiles = args.seqgzfolder
    gzfiles.rstrip('/')
    if not os.path.isdir(gzfiles):
        print(f"Error: '{gzfiles}' does not exist")
        sys.exit(1)

    # run get_ncbi_tax_tree_no_species.py
    print(colored.yellow("MAKING TREE"), taxon, colored.yellow(emoticons.get_ran_emot("excited")))
    tname = dirl+"/"+taxon+".tre"

    if TAXALISTF is not None:
        cmd = py+" "+DI+"get_ncbi_tax_tree_no_species.py "+taxon+" "+db+" "+TAXALISTF+" > "+tname
    else:
        cmd = py+" "+DI+"get_ncbi_tax_tree_no_species.py "+taxon+" "+db+" > "+tname

    os.system(cmd)

    trn = tree_reader.read_tree_file_iter(tname).__next__().label
    ## DEBUG
    print("DEBUG done with get_ncbi_tax_tree_no_species.py")
    print(f"DEBUG trn: {trn}")
    input("DEBUG Press Enter to continue...")

    # run make_dirs.py
    print(colored.yellow("MAKING DIRS IN"), dirl, colored.yellow(emoticons.get_ran_emot("excited")))
    cmd = py+" "+DI+"make_dirs.py "+tname+" "+dirl
    os.system(cmd)

    # run populate_dirs_first.py
    print(colored.yellow("POPULATING DIRS"), dirl, colored.yellow(emoticons.get_ran_emot("excited")))
    cmd = py+" "+DI+"populate_dirs_first.py "+tname+" "+dirl+" "+db+" "+gzfiles

    if TAXALISTF is not None:
        cmd += " "+TAXALISTF

    os.system(cmd)
    ## DEBUG
    print("DEBUG done with populate_dirs_first.py")
    input("DEBUG Press Enter to continue...")

    if os.path.isfile("log.md.gz"):
        os.remove("log.md.gz")

    # run cluster_tree.py
    cmd = py+" "+DI+"cluster_tree.py "+dirl+"/"+trn+"/ "+logfile
    os.system(cmd)
    ## DEBUG
    print("DEBUG done with cluster_tree.py")
    input("DEBUG Press Enter to continue...")

    print(colored.blue("PYPHLAWD DONE "+emoticons.get_ran_emot("excited")))
    end = datetime.now()
    print(colored.blue("Total time (H:M:S): "+str(end-start)+" "+emoticons.get_ran_emot("excited")))
    print(bcolors.HEADER, end=' ')
    emoticons.animate(emoticons.glasses_animated)
    print(bcolors.ENDC)
