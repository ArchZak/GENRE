from bioblend.galaxy import GalaxyInstance
from dotenv import load_dotenv
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
import concurrent.futures
import subprocess
import time
import os
import ssl

load_dotenv()

gi = GalaxyInstance('https://usegalaxy.org', key=os.getenv("GALAXY_API_KEY"))
history_id = os.getenv("HISTORY_ID")

def run_blastp(fasta_file=None, database="nr"):
    if not os.path.exists(fasta_file): # SHOULD LEAVE ERROR HANDLING TO RUN FUNCTION AS MUCH AS POSSIBLE
        raise FileNotFoundError(f"given fasta file not found: {fasta_file}")
    
    ssl._create_default_https_context = ssl._create_unverified_context # ssl verification required to use blastp over the internet
    
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    
    if not sequences:
        raise ValueError(f"no sequences found in {fasta_file}") #check for fasta file
    
    for sequence in sequences:
        try:
            protein_sequence = sequence.seq.translate(to_stop=True) #convert to protein seq so blastp can run
            
            result_handle = NCBIWWW.qblast( 
                program="blastp", 
                database=database,
                sequence=protein_sequence,
                hitlist_size=10 
            )
            
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)
            
            with open("blastp_downloads/blastp_results.txt", "w") as file: #creates blastp file with top 10 hits
                file.write(f"\nResults for {sequence.id}:\n")
                file.write("-" * 100 + "\n")
                file.write(f"{'Hit Description':<60} {'E-value':<15} {'Score':<10} {'Identity %':<10}\n")
                file.write("-" * 100 + "\n")
                
                if len(blast_record.alignments) == 0: 
                    file.write("No hits found\n")
                else:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            identity_pct = (hsp.identities / hsp.align_length) * 100
                            desc = alignment.title
                            file.write(f"{desc:<60} {hsp.expect:<15.2e} {hsp.score:<10.1f} {identity_pct:<10.1f}\n")

        except Exception as e:
            print(f"error processing sequence: {str(e)}")
            continue

def run_hhsearch(file=None):
    tool_id = "toolshed.g2.bx.psu.edu/repos/guerler/hhsearch/hhsearch/3.2.0+galaxy0"

    protein_file = translate_fasta(file)

    upload_response = gi.tools.upload_file(protein_file, history_id) 
    dataset_id = upload_response.get('outputs')[0].get('id') 
    upload_job_id = upload_response.get('jobs', [])[0]['id']

    upload_job_state = ""
    while upload_job_state != 'ok':  # checks job status before moving on # UNDERSTAND WHAT DIFFERENT JOB STATES EXIST TO HAVE BETTER ERROR HANDLING
        upload_job = gi.jobs.show_job(upload_job_id, True)
        upload_job_state = upload_job.get('state')
        time.sleep(10)

    inputs = {
        'i': {'src': 'hda', 'id': dataset_id},  
        'method': 'hhsearch',  
        'db_source': {
            'db_source_selector': 'indexed',  
            'ffindex': 'pdb70_hmm_2021-03' 
        },
        'e': 0.0000001  # 1*10^-7 cut off
    }

    hhsearch_run = gi.tools.run_tool(history_id,tool_id,inputs)
    hhsearch_job_id = hhsearch_run.get('jobs', [])[0].get('id')

    hhsearch_job_state = ""
    while hhsearch_job_state != 'ok': # checks job status before moving on # UNDERSTAND WHAT DIFFERENT JOB STATES EXIST TO HAVE BETTER ERROR HANDLING
        hhsearch_job = gi.jobs.show_job(hhsearch_job_id, True)
        hhsearch_job_state = hhsearch_job.get('state')
        time.sleep(10)

    hhsearch_results = gi.jobs.show_job(hhsearch_job_id, True)
    output_dataset_id = hhsearch_results.get('outputs', {}).get('output', {}).get('id')

    gi.datasets.download_dataset(output_dataset_id, "galaxy_downloads") # downloads results

def translate_fasta(input_file):
    protein_records = []
    
    for record in SeqIO.parse(input_file, "fasta"):
        protein_seq = record.seq.translate()
        
        protein_record = SeqRecord(
            protein_seq,
            id=record.id,
            description=f"translated from {record.description}"
        )
        protein_records.append(protein_record)
    
    SeqIO.write(protein_records, "output_protein.fasta", "fasta")
    return "output_protein.fasta"

def run_augustus(file=None):
    tool_id = "toolshed.g2.bx.psu.edu/repos/bgruening/augustus/augustus/3.4.0+galaxy3"

    upload_response = gi.tools.upload_file(file, history_id)
    dataset_id = upload_response.get('outputs')[0].get('id')
    upload_job_id = upload_response.get('jobs', [])[0]['id']

    upload_job_state = ""
    while upload_job_state != 'ok':   # checks job status before moving on # UNDERSTAND WHAT DIFFERENT JOB STATES EXIST TO HAVE BETTER ERROR HANDLING
        upload_job = gi.jobs.show_job(upload_job_id, True)
        upload_job_state = upload_job.get('state')
        time.sleep(10)

    inputs = {
        'input_genome': {'values': [{'id': dataset_id, 'src': 'hda'}]},
        'noInFrameStop': False,
        'singlestrand': False,
        'utr': False,
        'model': {
            'augustus_mode': 'builtin',
            '__current_case__': 1,
            'organism': 'human'  # ASK WHAT RECOMMENDED DEFAULT VALUE SHOULD BE
        },
        'softmasking': True,
        'strand': 'both',
        'genemodel': 'complete',
        'hints': {
            'usehints': 'F',
            '__current_case__': 1
        },
        'range': {
            'userange': 'F',
            '__current_case__': 1
        },
        'gff': False,
        'outputs': ['protein', 'codingseq', 'start', 'stop', 'cds']  
    }

    augustus_run = gi.tools.run_tool(history_id, tool_id, inputs)
    augustus_job_id = augustus_run.get('jobs', [])[0].get('id')

    augustus_job_state = ""
    while augustus_job_state != 'ok':  # checks job status before moving on # UNDERSTAND WHAT DIFFERENT JOB STATES EXIST TO HAVE BETTER ERROR HANDLING
        augustus_job = gi.jobs.show_job(augustus_job_id, True)
        augustus_job_state = augustus_job.get('state')
        time.sleep(10)

    augustus_results = gi.jobs.show_job(augustus_job_id, True)
    output_dataset_id = augustus_results.get('outputs', {}).get('output', {}).get('id')

    gi.datasets.download_dataset(output_dataset_id, "galaxy_downloads")

def run_genemark(fasta_file: str):
    if not os.path.exists(fasta_file): #checks to see if fasta_file exists, should return page error later
        raise FileNotFoundError(f"given fasta file not found: {fasta_file}")
    
    command = ["perl", "../gms2_macos/gms2.pl"] # run subprocess command to access tool
    command.extend([
        "--seq", fasta_file,
        "--genome-type", "auto"
    ])
    
    try:
        subprocess.run(
            command,
            text=True,
            capture_output=True,
            check=True
        )
            
    except subprocess.CalledProcessError as e: # subprocess error handling
        print(f"Error running GeneMarkS: {e}")
        if e.stderr: 
            print("\nError output:") 
            print(e.stderr)
        if e.stdout:
            print("\nStandard output:")
            print(e.stdout)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def run_all_tools():
    with concurrent.futures.ThreadPoolExecutor() as executor: #python library to run all functions simultaneously
        executor.submit(run_blastp, "Jollymon_gene_2.fasta", "nr")
        executor.submit(run_hhsearch, "Jollymon_gene_2.fasta")
        executor.submit(run_augustus, "Jollymon.fasta")
        executor.submit(run_genemark, "Jollymon.fasta")

run_all_tools()