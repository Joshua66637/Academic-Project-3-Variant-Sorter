#python3 Script/Variant Sorter.py --vcfFile Data/assessmentData.vcf.gz --gffFile Data/PlasmoDB-54_Pfalciparum3D7.gff --fastaFile Data/PlasmoDB-54_Pfalciparum3D7_Genome.fasta

import os
import math
import sys
import argparse
import logging

from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError
import vcfpy
import gffutils as gf

import seaborn.objects as so

#User defined functions
def isRightFormat(file, format, logobj):
    """
    Checks that a file is of the right format and logs the result
    :param file: Name of file to be checked
    :param format: Expected format
    :logobj: Logger object to log the result to
    """
    if(str(file).endswith(format)):    
        logobj.info(format[1:] + " file " + str(file) + " accepted")
        return True
    else:
        logobj.error(format[1:] + " file " + str(file) + " is not of the correct format. Please try again using a different file")
        return False

#Setting up command line interface
parser = argparse.ArgumentParser(description='Categorises variants as Synonymous, Non-synonymous or Non-coding\n', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--vcfFile', required=True, help='Bgzipped Variant Call Format File(.vcf.gz)')
parser.add_argument('--gffFile', required=True, help='General Feature Format file(.gff)')
parser.add_argument('--fastaFile', required=True, help='Genome fasta file(.fasta)')
files = parser.parse_args()

#setting up logger
logger = logging.getLogger("Variant Sorting Logger")
logger.setLevel(logging.INFO)
fh = logging.FileHandler("3047039_Variant_Sorting_logfile.log")
fh.setFormatter(logging.Formatter('%(levelname)s - %(asctime)s - %(message)s'))
logger.addHandler(fh)

#Accepting files and checking their formats
if(isRightFormat(files.vcfFile, ".vcf.gz", logger) and isRightFormat(files.gffFile, ".gff", logger) and isRightFormat(files.fastaFile, ".fasta", logger)):
    #vcread = vcfpy.Reader(filename= files.vcfFile) #Initialising VCF file
    vcread = vcfpy.Reader.from_path(files.vcfFile)
    tbi = files.vcfFile.replace('.vcf.gz', '.vcf.gz.tbi')
    database = files.gffFile.replace('.gff', '.db') #Initialising GFF file
    if os.path.isfile(tbi):
        logger.info("Detected Tabix Index for zipped vcf file")
    else:
        logger.error("Tabix Index for zipped vcf file is missing")
        sys.exit #Exit in case of missing file
    if not os.path.isfile(database):
        logger.info("Creating database" + database + "\n")
        db = gf.create_db(files.gffFile, dbfn=database, 
                      force=False, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    else:
        logger.info("Connecting to existing database" + database + "\n")
        db =  gf.FeatureDB(database, keep_order=True)
    fasta = files.fastaFile #Initialising Fasta File
else:
    sys.exit #Exit in case of improper file formats

Hcount = 0 #Counter for High Quality VCF records
Lcount = 0 #Counter for Low Quality VCF records

barchart = {"Synonymous": 0, "Non-synonymous": 0, "Non-coding":0}

Outputtable = "3047039_Variant_Table.tsv" #Defining output table
outline = []
outputdict = {"CHROM": "NA", "POS": "NA", "REF": "NA", "ALT": "NA", "Type": "NA", "Transcript": "NA", 
              "Protein Location": "NA", "Ref AA": "NA", "Alt AA": "NA"}
header = list(outputdict.keys())
outline.append("\t".join(header))

for record in vcread:
    output = {"Chrom": "NA", "Pos": "NA", "Ref": "NA", "Alt": "NA", "Type": "NA", "Transcript": "NA", 
              "Protein Location": "NA", "Ref AA": "NA", "Alt AA": "NA"}
    call_rate = sum(1 for call in record.calls if call.called) / len(record.calls)
    if record.QUAL > 20 and call_rate == 1.0 and record.is_snv:
        Hcount+=1
        output["Chrom"] = str(record.CHROM)
        output["Pos"] = str(record.POS)
        output["Ref"] = str(record.REF)
        output["Alt"] = str(record.ALT[0])
        output["Type"] = "Non-coding"
        for CDS in db.region(seqid=record.CHROM, start=record.POS, end= record.POS, featuretype= "CDS"):
            #If the snp lies in a coding region then get its parent
            for mRNA in db.parents(CDS.id, featuretype='mRNA'):
                mRNAPOS = 0 #variable to track the position of the variant on the transcript
                output["Transcript"] = mRNA.id
                seq = '' #to append the transcript

                if mRNA.strand == '+':
                    #for positive strand get all the CDS children of the mRNA and order by start position
                    for child in db.children(mRNA.id, featuretype='CDS', order_by='start'):
                        #build transcript to track snp postion both on transcript and translated protein
                        seq = seq + child.sequence(fasta, use_strand=True)
                        if CDS == child:
                            #for the CDS containing the mutation record its position on the transcript
                            mRNAPOS = mRNAPOS + (int(record.POS) - child.start + 1)
                        elif child.start < record.POS:
                            #for CDS prior to snp sum up the lengths
                            mRNAPOS = mRNAPOS + (child.end - child.start + 1)
                        
                elif mRNA.strand == '-':
                    #for negative strand get all the CDS children of the mRNA and reverse order by start position
                    for childrev in db.children(mRNA.id, featuretype='CDS', order_by='start', reverse=True):
                        #build transcript till snp to track snp postion both on transcript and translated protein
                        seq = seq + childrev.sequence(fasta, use_strand=True)
                        if CDS == childrev:
                            #for the CDS containing the mutation record its position on the transcript from 3' end
                            mRNAPOS = mRNAPOS + (childrev.end - int(record.POS) + 1)
                        elif childrev.end > record.POS:
                            #for CDS prior to snp sum up the lengths
                            mRNAPOS = mRNAPOS + (childrev.end - childrev.start + 1)
                transFlag = True
                try:
                    trans = Seq(seq).translate(cds=True)
                except TranslationError as e:
                    trans = Seq(seq).translate()
                    logger.warning(str(e) + " in reference sequence for " + output["Transcript"] + "\n")
                    transFlag = False
                protPOS= math.ceil(mRNAPOS/3)
                output["Protein Location"] = str(protPOS)
                output["Ref AA"] = str(trans[protPOS-1]) 
                alt = str(record.ALT[0]) if record.ALT else ""
                alt = "".join(filter(lambda x: x in "ACGTNacgtn", alt))  # clean sequence

                if mRNA.strand == '+':
                    altseq = seq[:mRNAPOS-1] + alt + seq[mRNAPOS:]
            
                elif mRNA.strand == '-':
                    altseq = seq[:mRNAPOS-1] + Seq(str(alt[0])).reverse_complement() + seq[mRNAPOS:] 
                try:
                    alttrans = Seq(altseq).translate(cds=True)
                except TranslationError as e:
                    alttrans = Seq(altseq).translate()
                    if(transFlag):
                        logger.warning(str(e) + " in " + output["Transcript"] + " after adding alt base\n")
                output["Alt AA"] = alttrans[protPOS-1]
                if(trans[protPOS-1] == alttrans[protPOS-1]):
                    output["Type"] = "Synonymous"
                else: 
                    output["Type"] = "Non-synonymous"
        barchart[output["Type"]] +=1
        outline.append("\t".join(list(output.values())))
    else:
        Lcount+=1

logger.info(str(Hcount) + " high quality VCF records (QUAL>20)")
logger.info(str(Lcount) + " low quality VCF records (QUAL<=20)")
logger.info(str(barchart["Non-coding"]) + " mutations are non-coding")
logger.info(str(barchart["Non-synonymous"]) + " mutations are non-synomymous")
logger.info(str(barchart["Synonymous"]) + " mutations are synonymous\n")

with(open(Outputtable, 'w') as out): 
    out.write('\n'.join(outline))

OutputGraph = "3047039_mutts.png"
barkeys = list(barchart.keys())
datapoints = [float(barchart[i]) for i in barkeys]
p = so.Plot(x= barkeys, y= datapoints)
p = p.add(so.Bar(color= "green")).label(title= "Types of Mutations", y = "Number of Mutations")
p.save(OutputGraph)

logger.info("Tab Separated Variants Table present in current working directory. Filename: " + Outputtable)
logger.info("Bargraph of mutations present in current working directory. Filename: " + OutputGraph)
