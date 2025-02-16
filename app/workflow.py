from bioblend.galaxy import GalaxyInstance
from dotenv import load_dotenv
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
import concurrent.futures
import linecache
import subprocess
import time
import os
import ssl

load_dotenv()

gi = GalaxyInstance('https://usegalaxy.org', key=os.getenv("GALAXY_API_KEY"))
history_id = os.getenv("HISTORY_ID")

def run_prokka(file=None):
    tool_id = "toolshed.g2.bx.psu.edu/repos/crs4/prokka/prokka/1.14.6+galaxy1"

    upload_response = gi.tools.upload_file(file, history_id)
    dataset_id = upload_response.get('outputs')[0].get('id')
    upload_job_id = upload_response.get('jobs', [])[0].get('id')

    upload_job_state = ""
    while upload_job_state != 'ok':   # checks job status before moving on # UNDERSTAND WHAT DIFFERENT JOB STATES EXIST TO HAVE BETTER ERROR HANDLING
        upload_job = gi.jobs.show_job(upload_job_id, True)
        upload_job_state = upload_job.get('state')
        time.sleep(1)

    inputs = {
        "input": {"values": [{"id": dataset_id, "src": "hda"}]},  
        "locustag": "GENOME",
        "increment": 1,
        "gffver": "3",
        "compliant_select": "no",
        "addgenes": False,
        "mincontig": 200,
        "centre": "XYZ",
        "genus": "Escherichia",
        "species": "coli",
        "strain": "K12",
        "plasmid": "",
        "kingdom_select": "Bacteria",
        "gcode": 11,
    }

    prokka_run = gi.tools.run_tool(history_id, tool_id, inputs)
    prokka_job_id = prokka_run.get('jobs', [])[0].get('id')

    prokka_job_state = ""
    while prokka_job_state != 'ok':
        prokka_job = gi.jobs.show_job(prokka_job_id, True)
        prokka_job_state = prokka_job.get('state')
        time.sleep(1)

    prokka_results = gi.jobs.show_job(prokka_job_id, True)
    output_dataset_id = prokka_results.get('outputs', {}).get('out_gff', {}).get('id')
    prokka_file = gi.datasets.download_dataset(output_dataset_id, "galaxy_downloads")
    return prokka_file

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

def count_genemark_genes(genemark_file):
    parsed_data = {}
    
    with open(genemark_file, 'r') as file:
        for line in file:
            if not line.strip():
                continue
            
            parts = line.split()
            if not parts[0].isnumeric():
                continue # make sure it's the line we want
            if len(parts) == 6:
                gene_number = int(parts[0])  # The gene number
                info = {
                    'strand_direction': parts[1],    
                    'start': (parts[2]) if parts[1] == '+' else (parts[3]),     
                    'stop': (parts[3]) if parts[1] == '+' else (parts[2]),     
                    'gene_length': int(parts[4]),
                }
                
                parsed_data[gene_number] = info
            if len(parts) < 9:
                continue # HANDLE THE < AND > LINES AHHHHH

            gene_number = int(parts[0])  # The gene number
            info = {
                'strand_direction': parts[1],    
                'start': int(parts[2]) if parts[1] == '+' else int(parts[3]),     
                'stop': int(parts[3]) if parts[1] == '+' else int(parts[2]),     
                'gene_length': int(parts[4]),
            }
            
            parsed_data[gene_number] = info
    
    return parsed_data

def count_prokka_genes(prokka_file):
    parsed_data = {}
    i=1
    with open(prokka_file, 'r') as file:
        for line in file:
            if not line.strip():
                continue
            
            parts = line.split()
            if len(parts) < 9:
                continue
            if parts[2] == 'CDS':
                info = {
                    'strand_direction': parts[6],
                    'start': int(parts[3]) if parts[6] == '+' else int(parts[4]),
                    'stop': int(parts[4]) if parts[6] == '+' else int(parts[3]), #basically checks for strand and assigns appropriately
                    'gene_length': int(parts[4])-int(parts[3])+1,
                }

                parsed_data[i] = info
                i+=1
    
    return parsed_data

def run_blastp(fasta_file=None):
    if not os.path.exists(fasta_file): #checks to see if fasta_file exists, should return page error later
        raise FileNotFoundError(f"given fasta file not found: {fasta_file}")
    
    ssl._create_default_https_context = ssl._create_unverified_context #wont run without ssl verification
    
    sequences = list(SeqIO.parse(fasta_file, "fasta")) #parse fasta file
    
    if not sequences:
        raise ValueError(f"no sequences found in {fasta_file}") #check for fasta file
    
    for sequence in sequences:
        
        try:
            protein_sequence = sequence.seq.translate(to_stop=True) #convert to protein seq
            
            result_handle = NCBIWWW.qblast( 
                program="blastp", #run blastp to find gene function
                database="nr",
                sequence=protein_sequence,
                hitlist_size=10
            )
            
            blast_records = NCBIXML.parse(result_handle)
            blast_record = next(blast_records)
            
            with open("blastp_downloads/blastp_results.txt", "w") as file:
                file.write("Hit Description,E-value,Score,Identity %\n")  # CSV-style header

                if len(blast_record.alignments) == 0:
                    raise ValueError("No hits found")
                else:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            identity_pct = (hsp.identities / hsp.align_length) * 100
                            desc = alignment.title
                            file.write(f"{desc},{hsp.expect:.2e},{hsp.score:.1f},{identity_pct:.1f}\n")

            
        except Exception as e:
            print(f"error processing sequence : {str(e)}")
            continue

def parse_top_10_blast_hits():
    top_10_hits = {}
    i=1
    file = "app/blastp_downloads/blastp_results.txt"

    for line_number in range(2,12):
        line = linecache.getline(file, line_number)
        if line:
            split_line = line.split(',')
            info = {
                'description': split_line[0],
                'e-value': split_line[1],
                'identity': float(split_line[3])
            }

            top_10_hits[i] = info
            i+=1

    linecache.clearcache()
    
    return top_10_hits
            
def parse_top_10_hh_hits(file_path):
    top_10_hits = {}
    i=1

    with open(file_path, "r") as file:
        for line in file:
            match = re.match(r"\s*\d+\s+(.+?)\s+([\d.]+)\s+", line)
            if match:
                description = match.group(1).strip()
                probability = float(match.group(2))
                
                info = {
                    'description': description,
                    'probability': probability
                }

                top_10_hits[i] = info
                i+=1
            
            if i > 10:
                return top_10_hits
    

def run_hhsearch(fasta_file=None):
    tool_id = "toolshed.g2.bx.psu.edu/repos/guerler/hhsearch/hhsearch/3.2.0+galaxy0"

    upload_response = gi.tools.upload_file(fasta_file, history_id)
    dataset_id = upload_response.get('outputs')[0].get('id')
    upload_job_id = upload_response.get('jobs', [])[0].get('id')


    upload_job_state = ""
    while upload_job_state != 'ok':
        upload_job = gi.jobs.show_job(upload_job_id, True)
        upload_job_state = upload_job.get('state')
        time.sleep(1)

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
    while hhsearch_job_state != 'ok':  
        hhsearch_job = gi.jobs.show_job(hhsearch_job_id, True)
        hhsearch_job_state = hhsearch_job.get('state')
        time.sleep(1)

    hhsearch_results = gi.jobs.show_job(hhsearch_job_id, True)
    output_dataset_id = hhsearch_results.get('outputs', {}).get('output', {}).get('id')

    hhsearch_file = gi.datasets.download_dataset(output_dataset_id, "galaxy_downloads")
    return hhsearch_file

def temp(prokka_g,genemark_g): #temporary structure to show an example of user workflow
    sentinel = True

    print(f"{len(genemark_g)} genes were found by GeneMarkS2 in the nucleotide sequence, which ones would you like to submit for LLM feedback?")
    print(f"{len(prokka_g)} genes were found by Prokka in the nucleotide sequence, which ones would you like to submit for LLM feedback?")
    while sentinel:
        gene_number = int(input("Enter gene number to view start and stop of: "))

        print(f'{genemark_g[gene_number]}')
        print(f'{prokka_g[gene_number]}')
        option = input('Would you like to see the potential function of the gene? (y or n): ')

        if option == 'y':
            sentinel = False
            return gene_number

def main():
    # user will be asked to input their fasta nuc of choice here
    print("Please wait for Prokka and GeneMarkS2 to finish running, this may take a moment")
    with concurrent.futures.ThreadPoolExecutor() as executor:
        prokka_future = executor.submit(run_prokka, "Jollymon.fasta")
        executor.submit(run_genemark, "Jollymon.fasta")
    genemark_genes = count_genemark_genes("gms2.lst")
    prokka_file = prokka_future.result()
    prokka_genes = count_prokka_genes(prokka_file)
    gene_number = temp(prokka_genes,genemark_genes) #add something to grab gene 2 seq
    with concurrent.futures.ThreadPoolExecutor() as executor2:
        executor2.submit(run_blastp, "Jollymon_gene_2.fasta") #running gene 2 for example sake right now, 
        future_hhsearch = executor2.submit(run_hhsearch,"output_protein.fasta")
    top_10_blast_hits = parse_top_10_blast_hits()
    top_10_hh_hits = parse_top_10_hh_hits(future_hhsearch)
    



main()