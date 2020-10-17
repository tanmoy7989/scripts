"""
Create a typically requested multiscale rendering of an integrative model by
combining detailed structures from the PDB files of rigid bodies used
in the model construction as well as flexible beads extracted from a
representative RMF file (usually the centroid of the top cluster after
structural clustering).

RMF is a specialized file format for storing
multiscale, multi-representation coordinate data from integrative models
generated with the Salilab software IMP (https://integrativemodeling.org/). Further details about the RMF format can be found at: 
https://github.com/salilab/rmf

This script will create:
a) A PDB file constructed by stitching together the different rigid
body PDB files after aligning them with their corresponding particles in
the given RMF file.

b) If requested, will write out a chimerax script that displays the
multiscale (PDB + unstructured parts from RMF) nature of the model.

c) Additionally, if crosslink (XL) data is provided, then it will also
map the crosslinks on to the multiscale structure.

Requirements:
IMP (https://integrativemodeling.org/)
Biopython (https://biopython.org/)
pandas (https://pandas.pydata.org/)

Author: Tanmoy Sanyal, Sali lab, UCSF
"""

import os
import numpy as np
import string
import pandas as pd
import argparse

from Bio.PDB import PDBParser, PDBIO
from Bio.PDB import Residue, Chain, Model
from Bio.SVDSuperimposer import SVDSuperimposer

import RMF
import IMP
import IMP.rmf
import IMP.core
from IMP.pmi.topology import TopologyReader

# constants
CHIMERAX_PDB_MODEL_NUM = 1
CHIMERAX_RMF_MODEL_NUM = 2
STRUCT_COLOR = "tan"
UNSTRUCT_COLOR = "tan"
SAT_XL_COLOR = "blue"
VIOL_XL_COLOR = "red"
XL_RADIUS = 1.0


def assembl e_multiscale_visualization(topology_fn, rmf_fn, pdb_dir,
                                      outprefix=None, chimerax=True,
                                      xl_fn=None):
    """
    Render multiscale versions of rigid bodies from PDB files + flexible
    beads from RMF files w/o mapped crosslinks.
    
    Args: 
    topology_fn (str): Topolgy file in pipe-separated-value (PSV) format
    as required in integrative modeling using IMP. For details on how
    to write a topology file, see:
    https://integrativemodeling.org/2.13.0/doc/ref/classIMP_1_1pmi_1_1topology_1_1TopologyReader.html
        
    rmf_fn (str): Name of the RMF file.
    
    pdb_dir (str): Directory containing all the PDB files for the rigid
    bodies used in modeling.
    
    outprefix (str, optional): Prefix for output files. Defaults to None.
    
    chimerax (bool, optional): If true, a Chimerax script will be written (extension ".cxc"). Defaults to True.
    
    xl_fn (str, optional): A file containing a XL dataset. Defaults to None.
    If this dataset is supplied, then it will be mapped on to the overall 
    structure with satisfied XLs drawn in blue and violated XLs drawn in red.
    A XL dataset should be supplied in a comma-separated-value (CSV) format
    containing at least the following fields
    
    protein1, residue1, protein2, residue2, sat
    
    where the last field <sat> is a boolean 1 or 0 depending on whether
    the particular XL is satisfied (in the ensemble sense) as a result of the
    integrative modeling exercise.
    """
    
    # -------------------------------------------
    # read the RMF file and extract all particles
    # -------------------------------------------
    of = RMF.open_rmf_file_read_only(rmf_fn)
    rmf_model = IMP.Model()
    hier = IMP.rmf.create_hierarchies(of, rmf_model)[0]
    IMP.rmf.load_frame(of, 0)
    particles = IMP.core.get_leaves(hier)
    rmf_ps = {}
    for p in particles:
        molname = p.get_parent().get_parent().get_parent().get_name().strip()
        name = p.get_name().strip()
        coord = IMP.core.XYZ(p).get_coordinates()
        rmf_ps[(molname, name)] = coord
        
    # --------------------------------------------------------------
    # map pdb residues to rmf particles for each rigid body pdb file
    # --------------------------------------------------------------
    # read the topology file
    t = TopologyReader(topology_fn, pdb_dir=pdb_dir)
    components = t.get_components()

    map_pdb2rmf = {}
    rigid_body_models = {}
    rigid_body_residues = {}
    chain_ids = {} # these are matched to the chimerax rmf plugin
    chain_id_count = 0
    for c in components:
        # ignore unstructured residues
        if c.pdb_file == "BEADS": continue
        mol = c.molname
        pdb_prefix = os.path.basename(c.pdb_file).split(".pdb")[0]
        chain_id = c.chain
        resrange = c.residue_range
        offset = c.pdb_offset
        
        if mol not in chain_ids:
            chain_ids[mol] = string.ascii_uppercase[chain_id_count]
            chain_id_count += 1
        
        if pdb_prefix not in map_pdb2rmf:
            map_pdb2rmf[pdb_prefix] = {}
            this_rigid_body_model = PDBParser().get_structure("x", c.pdb_file)[0]
            this_rigid_body_residues = {(r.full_id[2], r.id[1]): r for r in this_rigid_body_model.get_residues()}
            rigid_body_models[pdb_prefix] = this_rigid_body_model
            rigid_body_residues[pdb_prefix] = this_rigid_body_residues
            
        r0 = resrange[0] + offset
        r1 = resrange[1] + 1 + offset
        for r in range(r0, r1):
            key = (chain_id, r)
            val = (mol, r)
            if key in rigid_body_residues[pdb_prefix]:
                map_pdb2rmf[pdb_prefix][key] = val
    
    # --------------------------------
    # align all pdb files with the rmf
    # --------------------------------
    print("\nAligning all rigid body structures...")
    align = SVDSuperimposer()
    for pdb_prefix, mapper in map_pdb2rmf.items():
        pdb_coords = []
        pdb_atoms = []
        rmf_coords = []
        
        residues = rigid_body_residues[pdb_prefix]
        for (chain, pdbres), r in residues.items():
            pdb_coords.append(r["CA"].coord)
            pdb_atoms.extend([a for a in r.get_atoms()])
                
            mol, rmf_res = mapper[(chain, pdbres)]
            rmf_coords.append(rmf_ps[(mol, str(rmf_res))])        
                 
        pdb_coords = np.array(pdb_coords)
        rmf_coords = np.array(rmf_coords)
        align.set(rmf_coords, pdb_coords)
        align.run()
        rotmat, vec = align.get_rotran()
        [a.transform(rotmat, vec) for a in pdb_atoms]
  
    # --------------------------
    # assemble the composite pdb
    # --------------------------
    mols = set(sorted([c.molname for c in components]))
    print("\nChain IDs by molecule:")
    for k, v in chain_ids.items():
        print("molecule %s, chain ID %s" % (k, v))
    
    reslists = {mol: [] for mol in mols}
    for pdb_prefix, mapper in map_pdb2rmf.items():
        residues = rigid_body_residues[pdb_prefix]
        for (chain, pdbres), r in residues.items():
            mol, resid = mapper[(chain, pdbres)]
            new_id = (r.id[0], resid, r.id[2])
            new_resname = r.resname
            new_segid = r.segid
            new_atoms = r.get_atoms()
            new_residue = Residue.Residue(id=new_id, resname=new_resname, segid=new_segid)
            [new_residue.add(a) for a in new_atoms]
            reslists[mol].append(new_residue)
    
    composite_model = Model.Model(0)
    for mol, chain_id in chain_ids.items():
        this_residues = sorted(reslists[mol], key=lambda r: r.id[1])
        this_chain = Chain.Chain(chain_id)
        [this_chain.add(r) for r in this_residues]
        composite_model.add(this_chain)
    
    # save the composite pdb to file
    io = PDBIO()
    io.set_structure(composite_model)
    if outprefix is None:
        outprefix = "centroid_model"
    io.save(outprefix + ".pdb")

    # -------------------------------------------------------------------
    # chimerax rendering (hide most of the rmf except unstructured beads)
    # -------------------------------------------------------------------
    if not chimerax: exit()
    print("\nWriting UCSF Chimerax script...")
    s = ""
    s += "open %s\n" % (outprefix + ".pdb")
    s += "open %s\n" % rmf_fn
    s += "color #%d %s\n" % (CHIMERAX_PDB_MODEL_NUM, STRUCT_COLOR)
    s += "color #%d %s\n" % (CHIMERAX_RMF_MODEL_NUM, UNSTRUCT_COLOR)
    s += "hide #%d\n" % CHIMERAX_RMF_MODEL_NUM
    
    struct_residues = []
    for key, val in map_pdb2rmf.items():
        struct_residues.extend(list(val.values()))
    
    unstruct_atomspec = {}
    for p in rmf_ps:
        molname, particle_name = p
        rmf_chain_id = chain_ids[molname]
        if "bead" in particle_name:
            r0, r1 = particle_name.split("_")[0].split("-")
            r0 = int(r0) ; r1 = int(r1)
            this_atomspec = "#%d/%s:%d-%d" % \
                            (CHIMERAX_RMF_MODEL_NUM, rmf_chain_id, r0, r1)
            for r in range(r0, r1+1):
                unstruct_atomspec[(molname, r)] = this_atomspec
        else:
            if (molname, int(particle_name)) not in struct_residues:
                r = int(particle_name)
                this_atomspec = "#%d/%s:%d" % \
                (CHIMERAX_RMF_MODEL_NUM, rmf_chain_id, r)
                unstruct_atomspec[(molname, r)] = this_atomspec
                
    s += "show %s\n" % (" ".join(set(unstruct_atomspec.values())))

    # ----------------------------------------------------------
    # if crosslink data is supplied, write out a pseudobond file
    # ----------------------------------------------------------
    if xl_fn is not None:
        # parse XL data
        df = pd.read_csv(os.path.abspath(xl_fn))
        xls = []
        for i in range(len(df)):
            this_df = df.iloc[i]
            p1 = this_df["protein1"] ; r1 = this_df["residue1"]
            p2 = this_df["protein2"] ; r2 = this_df["residue2"]
            #sat = this_df["sat"]
            sat = 1
            xls.append((p1, r1, p2, r2, sat))
        
        # get lists of struct atomspecs
        atomspec = {}
        for (mol, particle_name) in rmf_ps:
            if "bead" in particle_name: continue
            if (mol, int(particle_name)) in unstruct_atomspec: continue
            chain_id = chain_ids[mol]
            resid = int(particle_name)
            atomspec[(mol, resid)] = "#%d/%s:%d@CA" % \
                                     (CHIMERAX_PDB_MODEL_NUM, chain_id, resid)
        
        # now add in all the unstruct atomspecs
        atomspec.update(unstruct_atomspec)

        # write pseudobond script
        s_pb = ""
        s_pb += "; radius = %2.2f\n" % XL_RADIUS
        s_pb += "; dashes = 0\n"
        for xl in xls:
            p1, r1, p2, r2, sat = xl
            atomspec_1 = atomspec[(p1, r1)]
            atomspec_2 = atomspec[(p2, r2)]
            if atomspec_1 == atomspec_2:
                continue
            color = SAT_XL_COLOR if sat else VIOL_XL_COLOR
            s_pb += "%s %s %s\n" % (atomspec_1, atomspec_2, color)
        s_pb += "\n"
        pb_fn = outprefix + "_XLs.pb"
        with open(pb_fn, "w") as of:
            of.write(s_pb)        
        s += "open %s\n" % pb_fn
            
    s += "preset 'overall look' publication\n"
    chimerax_out_fn = outprefix + ".cxc"
    with open(chimerax_out_fn, "w") as of:
        of.write(s)


### SCRIPT MODE ###
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument("-t", "--topology_fn", 
                        help="IMP (pmi) style topology file.")
    
    parser.add_argument("-r", "--rmf_fn", help="RMF file.")
    
    parser.add_argument("-p", "--pdb_dir", 
                        help="Directory containing all PDB files.")
    
    parser.add_argument("-cx", "--chimerax", action="store_true",
                        help="True to write out UCSF Chimerax script.")
    
    parser.add_argument("-x", "--xl_fn", help="Crosslink dataset CSV file.")
    
    parser.add_argument("-o", "--outprefix", help="Prefix of all output files.")
    
    args = parser.parse_args()
    xl_fn = None if args.xl_fn is None else os.path.abspath(args.xl_fn)
    
    assemble_multiscale_visualization(
        topology_fn=os.path.abspath(args.topology_fn),
        rmf_fn=os.path.abspath(args.rmf_fn),
        pdb_dir=os.path.abspath(args.pdb_dir),
        outprefix=args.outprefix,
        chimerax=args.chimerax,
        xl_fn=xl_fn)
