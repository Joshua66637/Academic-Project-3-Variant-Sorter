#BCPY Mutation Categorizing Assignment for GUID: 3047039
#Basic Modules
import os
import math
import sys
import argparse
import logging
import sqlite3

#Biological Computing Modules
import vcf
import gffutils as gf
import seaborn.objects as so
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError #Used to check if a sequence is a cds

#User defined function
def isExistingInRightFormat(file, format, logobj):
    """
    Checks that a file is present and of the expected format and logs the result
    :param file: Name of file to be checked
    :param format: Expected format
    :param logobj: Logger object to log the result to
    :return: True if the file format is as expected and the file exists, else False
    """
    if(str(file).endswith(format)): #Check if file is of the expected format
        if(os.path.isfile(file)): #Check if the file exists
            logobj.info(format[1:] + " file " + str(file) + " accepted")
            return True
        else:
            logobj.error(format[1:] + " file " + str(file) + " does not exist. Please try again using a different file.")
            return False  
    else:
        logobj.error(format[1:] + " file " + str(file) + " is not of the correct format. Please try again using a different file.")
        return False

#Setting up argparse command line interface
parser = argparse.ArgumentParser(description = "This script categorises variants from a VCF file as Synonymous, Non-synonymous or Non-coding using information from a GFF file and a Fasta File\n", 
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#Defining an argument for each input file type
parser.add_argument('--vcfFile', required = True, help = "Bgzipped Variant Call Format File(.vcf.gz)")
parser.add_argument('--gffFile', required = True, help = "General Feature Format file(.gff)")
parser.add_argument('--fastaFile', required = True, help = "Genome fasta file(.fasta)")
files = parser.parse_args()

#Setting up log file
logger = logging.getLogger("Variant Sorting Logger")
logger.setLevel(logging.INFO)
fh = logging.FileHandler("3047039_Variant_Sorting_logfile.log") #Logfile name
fh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(fh)

#Accepting files and confirming their formats and existence
if(isExistingInRightFormat(files.vcfFile, ".vcf.gz", logger) and isExistingInRightFormat(files.gffFile, ".gff", logger)
    and isExistingInRightFormat(files.fastaFile, ".fasta", logger)):
    
    vcread = vcf.Reader(filename= files.vcfFile) #Initialising VCF file
    tbi = files.vcfFile.replace(".vcf.gz", ".vcf.gz.tbi")
    if os.path.isfile(tbi):
        logger.info("Detected Tabix Index for zipped vcf file")
    else:
        logger.error("Tabix Index for zipped vcf file is missing")
        sys.exit() #Exit in case of missing tabix file
    
    database = files.gffFile.replace(".gff", ".db") #Initialising GFF file
    if not os.path.isfile(database):
        logger.info("Creating database " + database + "\n") #Create database only if it does not already exist
        try:
            db = gf.create_db(files.gffFile, dbfn = database, force = False, 
                              keep_order = True, merge_strategy = "merge", sort_attribute_values = True)
        except sqlite3.OperationalError:
            logger.error("Failed to create database " + database + "\n")
            sys.exit() #Exit due to lack of database
    else:
        logger.info("Connecting to existing database " + database + "\n")
        try:
            db =  gf.FeatureDB(database, keep_order=True)
        except ValueError:
            logger.error("Failed to read existing database " + database + "\n")
            sys.exit() #Exit due to lack of database
        except sqlite3.OperationalError:
            logger.error("Database " + database + " is corrupted\n")
            sys.exit() #Exit due to malformatted database

    fasta = files.fastaFile #Initialising Fasta File

else:
    sys.exit() #Exit in case of improper file formats

#Initialising Counters
Hcount = 0 #Counter for High Quality VCF records
Lcount = 0 #Counter for Low Quality VCF records
Muttcount = {"Synonymous": 0, "Non-synonymous": 0, "Non-coding":0} #Counters for each category of mutation

#Defining outputs for mutation table
Outputtable = "3047039_Variant_Table.tsv"  #Output file
output = {"CHROM": "NA", "POS": "NA", "REF": "NA", "ALT": "NA", "Type": "NA", "Transcript": "NA", 
              "Protein Location": "NA", "Ref AA": "NA", "Alt AA": "NA"} #Default Values for each row
header = list(output.keys())

#Opening output file
with(open(Outputtable, 'w') as out):
    out.write("\t".join(header))
    #Reading VCF Records
    for record in vcread:
        #Only work on a record if its high quality
        if record.QUAL > 20 and record.call_rate == 1.0 and record.is_snp:
            Hcount += 1
            #Initialising default row values
            output = {"Chrom": str(record.CHROM), "Pos": str(record.POS), "Ref": str(record.REF), "Alt": str(record.ALT[0]), 
                    "Type": "Non-coding", "Transcript": "NA", "Protein Location": "NA", "Ref AA": "NA", "Alt AA": "NA"} 
            
            for CDS in db.region(seqid = record.CHROM, start = record.POS, end = record.POS, featuretype = "CDS"):
                #If the snp lies in a coding region then get its parent transcript
                for mRNA in db.parents(CDS.id, featuretype = "mRNA"):
                    mRNAPOS = 0 #Position tracker for the variant on the transcript
                    seq = "" #To append the transcript so that the frame of the variant can be ascertained
                    
                    output["Transcript"] = mRNA.id
                    if mRNA.strand == '+':
                        #For positive strand get all CDS children of the mRNA and order by start position
                        for child in db.children(mRNA.id, featuretype = "CDS", order_by = "start"):
                            #Build transcript to track snp postion both on transcript and translated protein
                            seq = seq + child.sequence(fasta, use_strand=True)
                            if CDS == child:
                                #For the CDS containing the snp record its position on the transcript
                                mRNAPOS = mRNAPOS + (int(record.POS) - child.start + 1)
                            elif child.start < record.POS:
                                #For CDS prior to snp sum up the lengths
                                mRNAPOS = mRNAPOS + (child.end - child.start + 1)
                            
                    elif mRNA.strand == '-':
                        #For negative strand get all the CDS children of the mRNA and reverse order by start position
                        for childrev in db.children(mRNA.id, featuretype = "CDS", order_by = "start", reverse=True):
                            #Build transcript to track snp postion both on transcript and translated protein
                            seq = seq + childrev.sequence(fasta, use_strand=True)
                            if CDS == childrev:
                                #For the CDS containing the mutation record its position on the transcript from 3' end
                                mRNAPOS = mRNAPOS + (childrev.end - int(record.POS) + 1)
                            elif childrev.end > record.POS:
                                #For CDS prior to snp sum up the lengths
                                mRNAPOS = mRNAPOS + (childrev.end - childrev.start + 1)
                    else:
                        logger.warning("GFF file has an mRNA with a malformatted strand (ID: " + mRNA.id + " ). Skipping row")
                        continue
                    
                    #Translating Reference and alternate sequences
                    transFlag = True #Flag to signal potential errors in translating reference sequence
                    try:
                        trans = Seq(seq).translate(cds = True) #Translate the reference sequence and check if it is a complete CDS
                    except TranslationError as e:
                        trans = Seq(seq).translate() #Translate the sequence and warn user that its not a complete CDS
                        logger.warning(str(e) + " in reference sequence for " + output["Transcript"] + "\n")
                        transFlag = False
                    
                    #Calculate position and name of reference Amino Acid
                    protPOS= math.ceil(mRNAPOS/3) 
                    output["Protein Location"] = str(protPOS)
                    output["Ref AA"] = trans[protPOS-1] 
                    
                    #Building the alternate sequence
                    alt = record.ALT
                    if mRNA.strand == '+': #For + Strand directly replace the reference base with the alternate base
                        altseq = seq[:mRNAPOS-1] + str(alt[0]) + seq[mRNAPOS:]      
                    elif mRNA.strand == '-': #For - Strand replace the reference base with reverse compliment of alternate base 
                        altseq = seq[:mRNAPOS-1] + Seq(str(alt[0])).reverse_complement() + seq[mRNAPOS:] 
                    try:
                        alttrans = Seq(altseq).translate(cds = True) #Translate the alternate sequence and check if it is a complete CDS
                    except TranslationError as e:
                        alttrans = Seq(altseq).translate() #Translate the sequence and warn user that its not a complete CDS
                        if(transFlag): #Dont give a warning if one has already been given for the reference sequence
                            logger.warning(str(e) + " in " + output["Transcript"] + " after adding alt base\n")
                    
                    output["Alt AA"] = alttrans[protPOS-1] #Find the alternate animo acid
                    
                    #Decide the type of mutation based on their definitions
                    if(output["Ref AA"] == output["Alt AA"]):
                        output["Type"] = "Synonymous"
                    else: 
                        output["Type"] = "Non-synonymous"
            
            Muttcount[output["Type"]] += 1 #Count each type of mutation
            
            #Write each row to the file
            out.write('\n')
            out.write("\t".join(list(output.values())))
            
        else:
            Lcount += 1 #Counting low quality variants

#Logging the mutation statistics
logger.info(str(Hcount) + " high quality VCF records (QUAL>20)")
logger.info(str(Lcount) + " low quality VCF records (QUAL<=20)")
logger.info(str(Muttcount["Non-coding"]) + " mutations are non-coding")
logger.info(str(Muttcount["Non-synonymous"]) + " mutations are non-synomymous")
logger.info(str(Muttcount["Synonymous"]) + " mutations are synonymous\n")

#Generating the Bargraph
OutputGraph = "3047039_muttGraph.png"
#Getting Values and mutation categories from the mutation counter dictionary
barkeys = list(Muttcount.keys()) 
datapoints = [float(Muttcount[i]) for i in barkeys]

#Plotting and saving the graph
p = so.Plot(x = barkeys, y = datapoints)
p = p.add(so.Bar(color = "green")).label(title = "Types of Mutations", y = "Number of Mutations")
p.save(OutputGraph)

#Logging the locations of the output files
logger.info("Tab Separated Variants Table is present in current working directory. Filename: " + Outputtable)
logger.info("Bargraph of mutations is present in current working directory. Filename: " + OutputGraph + "\n")