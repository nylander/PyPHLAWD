"""
get_subset_genbank

Usage: python get_subset_genbank.py tid db outfile

Positional arguments

    tid: taxon id
    db: database file
    outfile: outfile base name

This version of the file is updated to take into account the change from gi to acc

JN: rewrite Nov 2025

"""

import sys
import os
import sqlite3
import gzip
from bad_seqs import gids
from bad_taxa import taxonids
from exclude_patterns import patterns
from exclude_desc_patterns import desc_patterns
from conf import smallest_size
from conf import filternamemismatch

def clean_name(nm):
    """ clean_name """
    nm = nm.replace(",", " ")
    return nm

def get_seqs_from_gz(gzdir, filename, idstoget):
    """ JN: New get_seqs_from_gz """
    try:
        with gzip.open(gzdir+"/"+filename, 'rb') as fl:
            #fl = gzip.open(gzdir+"/"+filename,"r")
            idtoseq = {tid: "" for tid in idstoget}
            current_id = None
            in_sequence_block = False
            for line in fl:
                line = line.decode()
                if line.startswith("LOCUS"):
                    # Reset for a new record
                    current_id = None
                    in_sequence_block = False
                if line.startswith("ACCESSION"):
                    # Extract the accession ID (e.g., OZ185191)
                    parts = line.strip().split()
                    if len(parts) > 1:
                        current_id = parts[1]
                # Check if we are in the sequence block for an ID we want
                if current_id in idstoget:
                    if line.strip().startswith("ORIGIN"):
                        in_sequence_block = True
                        continue
                    if in_sequence_block:
                        # Sequence lines start with a number and then the sequence data
                        # Filter out line numbers and spaces
                        line_parts = line.strip().split()
                        if line_parts and line_parts[0].isdigit():
                            sequence_data = "".join(line_parts[1:])
                            idtoseq[current_id] += sequence_data
                        elif line.strip().startswith("//"):
                            # End of record
                            in_sequence_block = False
                            current_id = None
            fl.close()
        # Final check and cleanup
        for tid in list(idtoseq.keys()):
            if idtoseq[tid] == "":
                # Sequence not found, set to None as per original script's logic for "too big"
                idtoseq[tid] = None
        return idtoseq
    except (OSError, ValueError):
        print(f"ERROR: {OSError} {ValueError}")
        sys.exit(1)

def make_files_with_id(taxonid, DB, outfilen, outfile_tbln, gzfileloc,
        remove_genomes=False, limitlist=None, excludetax=None):
    """
    make_files_with_id
    If outfilen and outfile_tbln are None, the results will be returned
    TODO: Q: gzfileloc: where is it provided?
          A: when running make_files_with_id from inside populate_dirs_first.py via setup_clade_ap.py!
    """
    if outfilen is not None and outfile_tbln is not None:
        outfile = open(outfilen, "w")
        outfileg = None
        if remove_genomes:
            outfileg = open(outfilen+".genomes", "w")
        outfile_tbl = open(outfile_tbln, "w")
    retseqs = [] # return if filenames aren't given
    rettbs = [] # returning if filenames aren't given
    conn = sqlite3.connect(DB)
    try:
        c = conn.cursor()
    except Exception as ex:
        print(f"Error: could not connect to {DB}: {ex}")
        sys.exit(1)
    species = []
    stack = []
    stack.append(str(taxonid))
    # DEBUG JN
    #print(f"DEBUG {__file__}: stack: {stack}", file=sys.stderr)
    files_ids = {} # key is the file, value is a list of ids
    ids_props = {} # key is id, value is list of properties
    while len(stack) > 0:
        id = stack.pop()
        # DEBUG JN
        #print(f"DEBUG {__file__}: id: {id}", file=sys.stderr)
        if id in species:
            continue
        species.append(id)
        if str(id) in taxonids:
            continue
        if excludetax is not None:
            if str(id) in excludetax:
                continue
        c.execute("select name from taxonomy where ncbi_id = ? and name_class = 'scientific name'", (id, ))
        l = c.fetchall()
        for j in l:
            tname = str(j[0])
            # DEBUG JN
            #print(f"DEBUG {__file__}: tname: {tname}", file=sys.stderr)
        badpattern = False
        for i in patterns:
            if i in tname:
                badpattern = True
                break
        if badpattern:
            continue
        c.execute("select * from sequence where ncbi_id = ?", (id, )) # this will give the filename in the folder
        l = c.fetchall()
        for j in l:
            # if the title sequence name is not the same as the id name (first part)
            # then we skip it. sorry sequence! you are outta here
            if filternamemismatch:
                try:
                    if tname.split(" ")[0]+tname.split(" ")[1] != str(j[4]).split(" ")[0]+str(j[4]).split(" ")[1]:
                        continue
                except:
                    continue
            # catch bad seqs
            if str(j[3]) in gids or str(j[2]) in gids:
                continue
            # bad description
            bad_desc = False
            for k in desc_patterns:
                if k in str(j[4]):
                    bad_desc = True
                    break
            if bad_desc:
                continue
            # str(j[7]) is the seq but now it is the seq file
            # now str(j[5]) is the file
            tfilen = "seqs."+str(j[5])
            if tfilen not in files_ids:
                files_ids[tfilen] = []
            files_ids[tfilen].append(str(j[2]))
            ids_props[str(j[2])] = [str(j[0]), str(j[1]), str(j[2]), str(j[3]), str(clean_name(tname)), str(j[5])]
        c.execute("select ncbi_id from taxonomy where parent_ncbi_id = ?", (id, ))
        childs = []
        l = c.fetchall()
        for j in l:
            childs.append(str(j[0]))
            stack.append(str(j[0]))
    # get all the seqs from the file at once
    for fn in files_ids:
        # DEBUG JN
        #print(f"DEBUG {__file__}: fn: {fn}", file=sys.stderr)
        idstoseq = get_seqs_from_gz(gzfileloc, fn, files_ids[fn])
        for tid in idstoseq:
            # DEBUG JN
            #print(f"DEBUG {__file__}: tid: {tid}", file=sys.stderr)
            seqstr = idstoseq[tid]
            if seqstr is None: # too big
                continue
            # DEBUG JN
            print(f"DEBUG {__file__}: for tid {tid}, len(seqstr): {len(seqstr)}", file=sys.stderr)
            if len(seqstr) < smallest_size:
                continue
            # exclude bad taxa
            if ids_props[tid][1] in taxonids:
                continue
            badpattern = False
            # MAKE SURE THIS IS OK. TODO: THIS SHOULD BE DONE ABOVE
            #for i in patterns:
            #    if i in tname:
            #        badpattern = True
            #        break
            #if badpattern:
            #    continue
            if limitlist is not None and ids_props[tid][1] not in limitlist:
                continue
            # we are writing
            seqst = ">"+str(ids_props[tid][3]+"\n"+seqstr)
            # DEBUG JN
            #print(f"DEBUG {__file__}: ids_props[tid][3]: {ids_props[tid][3]}", file=sys.stderr)
            print(f"DEBUG {__file__}: seqstr: {seqstr}", file=sys.stderr)
            tblst = "\t".join(ids_props[tid])
            if outfilen is not None and outfile_tbln is not None:
                if remove_genomes:
                    if len(seqstr) > 10000:
                        outfileg.write(seqst+"\n")
                    else:
                        outfile.write(seqst+"\n")
                else:
                    outfile.write(seqst+"\n")
                outfile_tbl.write(tblst+"\n")
            # we are returning
            else:
                if len(seqstr) < 10000:
                    retseqs.append(seqst)
                    rettbs.append(tblst)
    # we are writing
    if outfilen is not None and outfile_tbln is not None:
        outfile.close()
        if remove_genomes:
            outfileg.close()
        outfile_tbl.close()
    # we are returning
    else:
        return retseqs,rettbs

def make_files_with_id_internal(taxonid, DB, outfilen, outfile_tbln, gzfileloc,
        remove_genomes=False, limitlist=None):
    """
    make_files_with_id_internal
    If outfilen and outfile_tbln are None, the results will be returned
    """
    if outfilen is not None and outfile_tbln is not None:
        outfile = open(outfilen, "w")
        outfileg = None
        if remove_genomes:
            outfileg = open(outfilen+".genomes", "w")
        outfile_tbl = open(outfile_tbln, "w")
    retseqs = [] # return if filenames aren't given
    rettbs = [] # returning if filenames aren't given
    files_ids = {} # key is the file, value is a list of ids
    ids_props = {} # key is id, value is list of properties
    conn = sqlite3.connect(DB)
    try:
        c = conn.cursor()
    except Exception as ex:
        print(f"Error: could not connect to {DB}: {ex}")
        sys.exit(1)
    c.execute("select name from taxonomy where ncbi_id = ? and name_class = 'scientific name'", (str(taxonid), ))
    l = c.fetchall()
    for j in l:
        tname = str(j[0])
    c.execute("select * from sequence where ncbi_id = ?", (str(taxonid), ))
    l = c.fetchall()
    for j in l:
        # if the title sequence name is not the same as the id name (first part)
        # then we skip it. sorry sequence! you are outta here
        if filternamemismatch:
            try:
                if tname.split(" ")[0]+tname.split(" ")[1] != str(j[4]).split(" ")[0]+str(j[4]).split(" ")[1]:
                    continue
            except:
                continue
        # catch bad seqs
        if str(j[3]) in gids or str(j[2]) in gids:
            continue
        # bad description
        bad_desc = False
        for k in desc_patterns:
            if k in str(j[4]):
                bad_desc = True
                break
        if bad_desc:
            continue
        # now str(j[5]) is the file
        tfilen = "seqs."+str(j[5])
        if tfilen not in files_ids:
            files_ids[tfilen] = []
        files_ids[tfilen].append(str(j[2]))
        ids_props[str(j[2])] = [str(j[0]), str(j[1]), str(j[2]), str(j[3]), str(clean_name(tname)), str(j[5])]
    # get the children of the taxon that have no children (and so the sequences would go here)
    keepers = []
    c.execute("select ncbi_id from taxonomy where parent_ncbi_id = ?", (str(taxonid), ))
    l = c.fetchall()
    for j in l:
        nt = str(j[0])
        c.execute("select ncbi_id from taxonomy where parent_ncbi_id = ?", (str(nt), ))
        m = c.fetchall()
        if len(m) == 0:
            keepers.append(nt)
    # get everything else for the table
    species = []
    stack = []
    stack.append(str(taxonid))
    # DEBUG JN
    print(f"DEBUG {__file__}: stack: {stack}", file=sys.stderr)
    while len(stack) > 0:
        id = stack.pop()
        # DEBUG JN
        #print(f"DEBUG {__file__}: id: {id}", file=sys.stderr)
        if id in species:
            continue
        species.append(id)
        if str(id) in taxonids:
            continue
        c.execute("select name from taxonomy where ncbi_id = ? and name_class = 'scientific name'", (id, ))
        l = c.fetchall()
        for j in l:
            tname = str(j[0])
            # DEBUG JN
            #print(f"DEBUG {__file__}: in make_files_with_id_internal, tname: {tname}", file=sys.stderr)
        badpattern = False
        for i in patterns:
            if i in tname:
                badpattern = True
                break
        if badpattern:
            continue
        c.execute("select * from sequence where ncbi_id = ?", (id, ))
        l = c.fetchall()
        # only record everything for the table
        for j in l:
            # if the title sequence name is not the same as the id name (first part)
            # then we skip it. sorry sequence! you are outta here
            if filternamemismatch:
                try:
                    if tname.split(" ")[0]+tname.split(" ")[1] != str(j[4]).split(" ")[0]+str(j[4]).split(" ")[1]:
                        continue
                except:
                    continue
            # catch bad seqs
            if str(j[3]) in gids or str(j[2]) in gids:
                continue
            if limitlist is not None and str(j[1]) not in limitlist:
                continue
            tfilen = "seqs."+str(j[5])
            if tfilen not in files_ids:
                files_ids[tfilen] = []
            if str(j[1]) in keepers:
                files_ids[tfilen].append(str(j[2]))
            ids_props[str(j[2])] = [str(j[0]), str(j[1]), str(j[2]), str(j[3]), str(clean_name(tname)), str(j[5])]
            tblst = "\t".join(ids_props[str(j[2])])
            # DEBUG JN
            print(f"DEBUG {__file__}: in make_files_with_id_internal, tblst: {tblst}", file=sys.stderr)
            outfile_tbl.write(tblst+"\n")
        c.execute("select ncbi_id from taxonomy where parent_ncbi_id = ?", (id, ))
        childs = []
        l = c.fetchall()
        for j in l:
            childs.append(str(j[0]))
            stack.append(str(j[0]))
    for fn in files_ids:
        idstoseq = get_seqs_from_gz(gzfileloc, fn, files_ids[fn])
        for tid in idstoseq:
            # DEBUG JN
            #print(f"DEBUG {__file__}: in make_files_with_id_internal, tid: {tid}", file=sys.stderr)
            seqstr = idstoseq[tid]
            if seqstr is None: # too big
                continue
            if len(seqstr) < smallest_size:
                continue
            # exclude bad taxa
            if ids_props[tid][1] in taxonids:
                continue
            # MAKE SURE THIS IS OK TODO
            #badpattern = False
            #for i in patterns:
            #    if i in tname:
            #        badpattern = True
            #        break
            #if badpattern:
            #    continue
            if limitlist is not None and ids_props[tid][1] not in limitlist:
                continue
            # we are writing
            seqst = ">"+str(ids_props[tid][3]+"\n"+seqstr)
            # DEBUG JN
            #print(f"DEBUG {__file__}: in make_files_with_id_internal, ids_props[tid][3] : {ids_props[tid][3]}", file=sys.stderr)
            print(f"DEBUG {__file__}: in make_files_with_id_internal; tid: {tid} seqstr : {seqstr}", file=sys.stderr)
            tblst = "\t".join(ids_props[tid])
            if outfilen is not None and outfile_tbln is not None:
            # if ids_props[tid][1] in keepers:
                if remove_genomes:
                    if len(seqstr) > 10000:
                        outfileg.write(seqst+"\n")
                    else:
                        outfile.write(seqst+"\n")
                else:
                    outfile.write(seqst+"\n")
                #outfile_tbl.write(tblst+"\n")
            # we are returning
            else:
                if len(seqstr) < 10000:
                    retseqs.append(seqst)
                    rettbs.append(tblst)
    # we are writing
    if outfilen is not None and outfile_tbln is not None:
        outfile.close()
        if remove_genomes:
            outfileg.close()
        outfile_tbl.close()
    # we are returning
    else:
        return retseqs,rettbs

def make_files_with_id_justtable(taxonid, DB, outfile_tbln):
    """
    make_files_with_id_justtable
    If you send outfile_tbln as None, it will return the results
    """
    if outfile_tbln is not None:
        outfile_tbl = open(outfile_tbln, "w")
    conn = sqlite3.connect(DB)
    try:
        c = conn.cursor()
    except Exception as ex:
        print(f"Error: could not connect to {DB}: {ex}")
        sys.exit(1)
    species = []
    stack = []
    tbl = []
    stack.append(str(taxonid))
    # DEBUG JN
    #print(f"DEBUG {__file__}: stack : {stack}", file=sys.stderr)
    while len(stack) > 0:
        id = stack.pop()
        # DEBUG JN
        #print(f"DEBUG {__file__}: id : {id}", file=sys.stderr)
        if id in species:
            continue
        species.append(id)
        c.execute("select name from taxonomy where ncbi_id = ? and name_class = 'scientific name'", (id, ))
        l = c.fetchall()
        for j in l:
            tname = str(j[0])
        c.execute("select * from sequence where ncbi_id = ?", (id, ))
        l = c.fetchall()
        for j in l:
            tbls = str(j[0])+"\t"+str(j[1])+"\t"+str(j[2])+"\t"+str(j[3])+"\t"+str(clean_name(tname))+"\t"+str(j[5])+"\t"+str(j[6])
            # DEBUG JN
            #print(f"DEBUG {__file__}: in make_files_with_id_justtable, tbls : {tbls}", file=sys.stderr)
            if outfile_tbln is not None:
                outfile_tbl.write(tbls+"\n")
            else:
                tbl.append(tbls)
        c.execute("select ncbi_id from taxonomy where parent_ncbi_id = ?", (id, ))
        childs = []
        l = c.fetchall()
        for j in l:
            childs.append(str(j[0]))
            stack.append(str(j[0]))
    if outfile_tbln is not None:
        outfile_tbl.close()
    else:
        return tbl

def make_files(taxon, DB, outfilen, outfile_tbln):
    """ make_files """
    outfile = open(outfilen, "w")
    outfile_tbl = open(outfile_tbln, "w")
    conn = sqlite3.connect(DB)
    try:
        c = conn.cursor()
    except Exception as ex:
        print(f"Error: could not connect to {DB}: {ex}")
        sys.exit(1)
    species = []
    stack = []
    c.execute("select ncbi_id from taxonomy where name = ?", (taxon, ))
    for j in c:
        stack.append(str(j[0]))
    # DEBUG JN
    #print(f"DEBUG {__file__}: stack : {stack}", file=sys.stderr)
    while len(stack) > 0:
        id = stack.pop()
        # DEBUG JN
        #print(f"DEBUG {__file__}: in make_files, id : {id}", file=sys.stderr)
        if id in species:
            continue
        species.append(id)
        c.execute("select * from sequence where ncbi_id = ?", (id, ))
        l = c.fetchall()
        for j in l:
            outfile.write(">"+str(j[3])+"\n")
            outfile.write(str(j[7])+"\n")
            # TODO: Undefined variable 'tname' below!
            outfile_tbl.write(str(j[0])+"\t"+str(j[1])+"\t"+str(j[2])+"\t"+str(j[3])+"\t"+str(clean_name(tname))+"\t"+str(j[4])+"\n")
        c.execute("select ncbi_id from taxonomy where parent_ncbi_id = ?", (id, ))
        childs = []
        l = c.fetchall()
        for j in l:
            childs.append(str(j[0]))
            stack.append(str(j[0]))
    outfile.close()
    outfile_tbl.close()

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage: python "+sys.argv[0]+" tid db outfile")
        sys.exit(0)
    tid = sys.argv[1]
    DB = sys.argv[2]
    if not os.path.isfile(DB):
        print(f"Error: '{DB}' does not exist")
        sys.exit(1)
    outfilen = sys.argv[3]
    outfile_tbln = sys.argv[3]+".table"
    make_files_with_id(tid, DB, outfilen, outfile_tbln) #  No value for argument 'gzfileloc' in function call (no-value-for-parameter)!
