from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML
import ssl

#really slow so it may require the use of setting up a local database

def run_blast_nucleotide(fasta_file="app/Jollymon_gene_2.fasta", database="nr"):
    ssl._create_default_https_context = ssl._create_unverified_context
    
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    
    if not sequences:
        raise ValueError(f"no sequences found in {fasta_file}")
    
    print(f"found {len(sequences)} sequences in {fasta_file}")
    
    for i, sequence in enumerate(sequences, 1):
        print(f"processing sequence {i} out of {len(sequences)}: {sequence.id}")
        
        try:
            protein_sequence = sequence.seq.translate(to_stop=True)
            
            result_handle = NCBIWWW.qblast(
                program="blastp",
                database=database,
                sequence=protein_sequence,
                hitlist_size=10
            )
            
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)
            
            print(f"\nResults for {sequence.id}:")
            print("-" * 100)
            print(f"{'Hit Description':<60} {'E-value':<15} {'Score':<10} {'Identity %':<10}")
            print("-" * 100)
            if len(blast_record.alignments) == 0: 
                print("No hits found")
                continue
                
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    identity_pct = (hsp.identities / hsp.align_length) * 100
                    desc = alignment.title[:57] + "..." if len(alignment.title) > 60 else alignment.title
                    print(f"{desc:<60} {hsp.expect:<15.2e} {hsp.score:<10.1f} {identity_pct:<10.1f}")
            
        except Exception as e:
            print(f"error processing sequence {sequence.id}: {str(e)}")
            continue

if __name__ == "__main__":
    run_blast_nucleotide()
