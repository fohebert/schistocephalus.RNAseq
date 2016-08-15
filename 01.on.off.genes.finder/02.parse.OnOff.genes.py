#!/usr/bin/python

""" 
v.1.0                 User's commands              v.1.0

\033[1mDESCRIPTION\033[0m
    Takes a read-count matrix file in which the raw
    read counts have been transformed into their CPM
    value (Counts Per Millions). These CPM values are
    usually obtained with the cpm() function from the
    limma/edgeR package and the transformation is performed
    on the DGE object. The command is used that way:
    > cpm(dge).
    
    This program will thus look at all the genes in all
    the samples to see if they are above or below the 
    pre-determined low expression cutoff. This cutoff
    is obtained by looking at the CPM value that corres-
    ponds to the 5th percentile of the distribution of
    all the observed CPM values. Below this CPM threshold
    the gene is considered not expressed at all. The goal
    here is to find the genes that have on/off patterns
    of expression, i.e. above the CPM threshold in one
    specific experimental group and below the CPM cutoff
    in all other groups. E.g.: a gene is above the CPM
    cutoff in at least 50% of the samples in the non-infective
    group AND below the CPM cutoff in ALL other samples
    from the other groups (infective and adult). This would
    mean that the gene is non-infective specific, i.e.
    it follows an on/off pattern and is only expressed at
    the non-infective stage. This program will also output
    the genes that follow an on/off pattern in a host-specific
    manner, i.e. systematically expressed in a given host
    and completely non-expressed in the other host(s). For
    this, we use the exact same CPM threshold.

\033[1mUSAGE\033[0m
    %program <in.CPM.file> <out.on.Inf> <out.on.Ninf> <on.Adult>
        <on.Fish>

\033[1mCREDITS\033[0m
    Doc Pants 2016 \m/
"""

import sys

try:
    in_file = sys.argv[1]
    out_on_inf = sys.argv[2]
    out_on_ninf = sys.argv[3]
    out_on_adult = sys.argv[4]
    out_on_fish = sys.argv[5]
except:
    print __doc__
    sys.exit(1)

sp_col = {} # Dictionnary with info on which sample = which column
on_ninf = {} # Dict() contains CPM values for genes = ON in non-infective
on_inf = {} # Dict() contains CPM values for genes = ON in infective
on_adult = {} # Dict() contains CPM values for genes = ON in adults
on_fish = {} # Dict() contains CPM values for genes = ON in fish only
with open(in_file, "rU") as i_f:
    for line in i_f:
        line = line.strip()
        
        # IF the line is actually the first one
        # keeps the info of which column is which
        # sample.
        if line.startswith("NI"):
        
            # Total number of samples
            num_sp = len(line.split("\t"))
            # We start with column number 1
            curr_col = 1
            count = 0
            # For each column, adds the column number
            # and the sample ID into the "sp.col" dict.
            while count <= num_sp-1:
                sp_col[curr_col] = line.split("\t")[count]
                curr_col += 1
                count += 1
        
        # If the line = read counts for a given transcript
        elif line.startswith("NI") == False:
            
            # Keeps the name of the transcript in an object
            trans = line.split("\t")[0] # Transcript ID = column 0
            
            # Objects that will contain the number of samples
            # in which we detect expression levels above the
            # pre-determined threshold.
            ninf = 0
            inf = 0
            adult = 0
            
            # Total number of samples
            num_sp = len(line.split("\t"))-1
            # We start with column number 1 (column 1 = sample 1)
            curr_col = 1
            # For each column, adds the column number
            # and the sample ID into the "sp.col" dict.
            while curr_col <= num_sp:
                # If the expression level for the current transcript,
                # in this particular individual that we are looking
                # at with the current column number, is above the pre-
                # defined CPM threshold (here it's 0.55).
                if float(line.split("\t")[curr_col]) > 0.55:
                    # A count is added in the appropriate category
                    # (i.e. ninf, inf or adult)
                    if sp_col[curr_col].split(".")[0] == "NI":
                        ninf += 1
                    elif sp_col[curr_col].split(".")[0] == "I":
                        inf += 1
                    elif sp_col[curr_col].split(".")[0] == "A":
                        adult += 1
                curr_col += 1
                
            # If the transcript is expressed above the threshold in
            # the non-infective samples only, it is added in the
            # corresponding dict().
            if ninf >= 3 and inf == 0 and adult == 0:
                # Only expressed in at least 50% of non-infective samples
                on_ninf[trans] = line.split("\t")[1:]
            elif inf >= 2 and ninf == 0 and adult == 0:
                # Only expressed in at least 50% of the infective samples
                on_inf[trans] = line.split("\t")[1:]
            elif adult >= 2 and ninf == 0 and inf == 0:
                # Only expressed in at least 66% of the adult samples
                on_adult[trans] = line.split("\t")[1:]
            elif ninf >= 3 and inf >= 2 and adult == 0:
                # Only expressed in the fish and not in the bird
                on_fish[trans] = line.split("\t")[1:]

# Creating the output file for the transcripts "ON" in infective
with open(out_on_inf, "w") as o_on_inf:
    
    # Outputs the header, i.e. sample names in the right order
    # For that, we use the sp_col dict() previously constructed.
    for col_num in sp_col:
        o_on_inf.write("\t" + sp_col[col_num])
    
    # Outputs the transcripts that are "ON" in infective right
    # after the header. We use the on_inf dict().
    for trans in sorted(on_inf):
        o_on_inf.write("\n" + trans + "\t" + "\t".join(on_inf[trans]))
    
# Creating the output file for the transcripts "ON" in non-infective
with open(out_on_ninf, "w") as o_on_ninf:
    
    # Outputs the header, i.e. sample names in the right order
    # For that, we use the sp_col dict() previously constructed.
    for col_num in sp_col:
        o_on_ninf.write("\t" + sp_col[col_num])
    
    # Outputs the transcripts that are "ON" in non-infective right
    # after the header. We use the on_ninf dict().
    for trans in sorted(on_ninf):
        o_on_ninf.write("\n" + trans + "\t" + "\t".join(on_ninf[trans]))

# Creating the output file for the transcripts "ON" in adult
with open(out_on_adult, "w") as o_on_adult:
    
    # Outputs the header, i.e. sample names in the right order
    # For that, we use the sp_col dict() previously constructed.
    for col_num in sp_col:
        o_on_adult.write("\t" + sp_col[col_num])
    
    # Outputs the transcripts that are "ON" in adult right
    # after the header. We use the on_adult dict().
    for trans in sorted(on_adult):
        o_on_adult.write("\n" + trans + "\t" + "\t".join(on_adult[trans]))

# Creating the output file for the transcripts "ON" in fish
with open(out_on_fish, "w") as o_on_fish:
    
    # Outputs the header, i.e. sample names in the right order
    # For that, we use the sp_col dict() previously constructed.
    for col_num in sp_col:
        o_on_fish.write("\t" + sp_col[col_num])
    
    # Outputs the transcripts that are "ON" in fish only right
    # after the header. We use the on_adult dict().
    for trans in sorted(on_fish):
        o_on_fish.write("\n" + trans + "\t" + "\t".join(on_fish[trans]))

print "\n\033[1mJOB COMPLETED\033[0m\n"
