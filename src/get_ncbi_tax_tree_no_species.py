""" get_ncbi_tax_tree_no_species """
import sys
import os
import sqlite3
import conf
if conf.usecython:
    import cnode as node
else:
    import node

exclude_un_en = True # unused?

def clean_name(name):
    """ clean name """
    return name.replace(" ", "_").replace("&", "_").replace(":", "_").replace("(", "_").replace(")", "_")

def get_all_included(taxalist, c):
    """ get_all_included """
    inc = set()
    for i in taxalist:
        cur = str(i)
        going = True
        while going:
            c.execute("select ncbi_id, parent_ncbi_id from taxonomy where ncbi_id = ?", (cur, ))
            count = 0
            for j in c:
                count += 1
                inc.add(str(j[0]))
                par = str(j[1])
                if par != "1" and par not in inc:
                    cur = par
                else:
                    going = False
                    break
            if count == 0:
                print("delete "+i)
    return inc

def construct_tree(taxon, db, taxalist=None):
    """ construct tree """
    conn = sqlite3.connect(db)
    c = conn.cursor()
    includelist = None
    if taxalist is not None:
        tl = set()
        tlf = open(taxalist,"r")
        for i in tlf:
            tl.add(i)
        tlf.close()
        includelist = get_all_included(tl,c)
    species = []
    stack = []
    done = set()
    rt = None
    # node_ids = {}  # id is key, value is parent id
    nodes = {}  # id is key, value is node
    if taxon.isdigit():
        c.execute("select ncbi_id from taxonomy where ncbi_id = ?", (taxon, ))
    else:
        c.execute("select ncbi_id from taxonomy where name = ? and node_rank != 'species'", (taxon, ))
    for j in c:
        stack.append(str(j[0]))
        rt = node.Node()
        rt.label = taxon+"_"+str(j[0])
        rt.data["id"] = str(j[0])
        nodes[str(j[0])] = rt
    while len(stack) > 0:
        id = stack.pop()
        if id in done:
            continue
        # added, if some error, remove but add on line 81
        if id not in nodes:
            continue
        done.add(id)
        c.execute("select ncbi_id,name,name_class,edited_name from taxonomy where parent_ncbi_id = ? and node_rank != 'species'", (id, ))
        childs = []
        for j in c:
            tid = str(j[0])
            if includelist is not None and tid not in includelist:
                continue
            childs.append(tid)
            stack.append(tid)
            if "unclassified" in str(j[1]) or "environmental" in str(j[1]):
                continue
            if str(j[2]) == "scientific name":
                #name = str(j[1])
                edname = str(j[3])
                nn = node.Node()
                nn.label = clean_name(edname)+"_"+str(tid)
                nn.data["id"] = tid
                nodes[tid] = nn
                # node_ids[tid] = id
                nn.parent = nodes[id]
                nodes[id].add_child(nn)
        if len(childs) == 0 and id not in species:
            species.append(id)
    return rt


if __name__ == "__main__":
    if len(sys.argv) != 3 and len(sys.argv) != 4:
        print("usage: python "+sys.argv[0]+" taxon db [taxalist]")
        sys.exit(0)
    taxon = sys.argv[1]
    if not taxon:
        print(f"Error: taxon '{taxon}' is empty")
        sys.exit(1)
    DB = sys.argv[2]
    if not os.path.isfile(DB):
        print(f"Error: DB '{DB}' does not exist")
        sys.exit(1)
    TAXALIST = None
    if len(sys.argv) == 4:
        TAXALIST = sys.argv[3]
    tree = construct_tree(taxon, DB, TAXALIST)
    ## DEBUG
    print(f"DEBUG {__file__}: done with construct_tree")
    print(f"DEBUG {__file__}: tree: {tree}")
    input(f"DEBUG {__file__}: Press Enter to continue...")
    print(tree.get_newick_repr(False)+";")
