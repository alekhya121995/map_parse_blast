###########################################
# Program: parse_blast_de.py
# Name: Alekhya Akkunuri
# Date: 11-8-2017
# Description: This program parses the blast and DE file
# to write a report of proteins and their corresponding
# stress condition data samples.
##############################################

# open the blast file
blast_file = open("blastp.outfmt6")
# read the blast file line by line such that each line is its own element
blast_line = blast_file.readlines()

# open output file
blast_output = open("map_parse_blast_output.txt", "w")

# function that parses BLAST line
def parse_blastline(line):
    # split at tabs
    line_split = line.rstrip("\n").split("\t")
    # store each required element into a variable
    transcriptId_isoform = line_split[0]
    swissprotId = line_split[1]
    # split the variables further
    subsplit1 = transcriptId_isoform.split("|")
    subsplit2 = swissprotId.split("|")
    transcript = subsplit1[0]
    swissprot = subsplit2[3]
    return transcript, swissprot
       
# store values in dictionary
transcript_swissprot = {}
for line in blast_line:
    (transcript, swissprot) = parse_blastline(line)
    transcript_swissprot[transcript] = swissprot   

# unpacking contents of differential expression file
def parse_matrix_contents(matrix):
    matrix_split = matrix.rstrip("\n").split("\t")
    (matrix_transcript, sp_ds, sp_hs, sp_log, sp_plat) = matrix_split
    
    # check for BLAST matches and print to output file
    if matrix_transcript in transcript_swissprot:
       protein = transcript_swissprot.get(matrix_transcript)
       output = (protein, sp_ds, sp_hs, sp_log, sp_plat)
       return(tuple(output))
    
    # use transcript as protein ID if no BLAST match is found and print to output file
    else: 
        output2 = (matrix_transcript, sp_ds, sp_hs, sp_log, sp_plat)
        return(tuple(output2))

def formatted_tuple(list_tuple):
     return ("\t".join(list_tuple))


# open the differentially expressed transcript file
de_file = open("diffExpr.P1e-3_C2.matrix")

# read the file line by line such that each line is its own element
de_contents = de_file.readlines()

map_tuples = map(parse_matrix_contents, de_contents)
formatted_map_tuples = map(formatted_tuple, map_tuples)
blast_output.write("\n".join(formatted_map_tuples))
