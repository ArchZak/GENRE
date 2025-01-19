from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def translate_fasta(input_file, output_file):
    protein_records = []
    
    for record in SeqIO.parse(input_file, "fasta"):
        protein_seq = record.seq.translate()
        
        protein_record = SeqRecord(
            protein_seq,
            id=record.id,
            description=f"translated from {record.description}"
        )
        protein_records.append(protein_record)
    
    SeqIO.write(protein_records, output_file, "fasta")

if __name__ == "__main__":
    translate_fasta("app/Jollymon_gene_2.fasta", "output_protein.fasta")