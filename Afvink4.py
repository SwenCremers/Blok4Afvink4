from flask import Flask, render_template, request
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Blast import NCBIWWW, NCBIXML
import re

app = Flask(__name__)

@app.route('/', methods=["GET"])
def sequence():
    ''' Passes trough correct data to html file according to what type
    of sequence is provided. Makes use of the BLAST function if that
    sequence is a protein sequence
    :return: renders html template
    '''

    # retrieves sequence from page
    seq = request.args.get("sequentie", "").strip()

    # Standard translation table
    standard_trans_table = CodonTable.ambiguous_dna_by_name["Standard"]

    # empty strings
    seqtype = ""
    protein = ""
    rna = ""
    gene = ""

    # if sequence exists
    if seq:

        # regexes for different seq types
        dna_re = '[^ACTG]'
        rna_re = '[^ACUG]'
        protein_re = '[^ARNDCEQHILKMFPSTWYCXUVG]'

        # DNA sequence
        if not re.search(dna_re, seq):
            seqtype = "DNA"
            bio_dna = Seq(seq, IUPAC.ambiguous_dna)
            protein = bio_dna.translate(table=standard_trans_table)
            rna = bio_dna.transcribe()

        # RNA sequence
        elif not re.search(rna_re, seq):
            seqtype = "RNA"

        # Protein sequence
        elif not re.search(protein_re, seq):
            seqtype = "Protein"
            gene = str(BLAST(seq))

        # No valid sequence
        else:
            seqtype = "Geen geldige sequentie"
    return render_template("Afvink4.html", sequentie=seq, type=seqtype, dna=seq
                           , rna=rna, protein=protein, gene=gene)

def BLAST(seq):
    ''' BLASTS a sequence and returns the description of the first hit
    :param seq: sequence string
    :return: description string
    '''
    try:

        # BLAST
        result = NCBIWWW.qblast('blastp', 'nr', seq, hitlist_size=1)

        # writes result to xml
        with open("BLAST.xml", 'w') as out_result:
            out_result.write(result.read())

        # opens the xml
        result_handle = open("BLAST.xml")

        # parses the xml
        blast_record = NCBIXML.parse(result_handle)

        # loops over the records of the result
        # and returns the description
        for record in blast_record:
            gene_desc = record.alignments[0].title
        gene_desc = gene_desc.split(" >")[0]
        return gene_desc
    except:
        return "Geen resultaten"




def main():
    app.run(debug=True)

main()