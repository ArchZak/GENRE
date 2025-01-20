from bioblend.galaxy import GalaxyInstance
from dotenv import load_dotenv
import time
import os

load_dotenv()

gi = GalaxyInstance('https://usegalaxy.org', key=os.getenv("GALAXY_API_KEY"))
history_id = os.getenv("HISTORY_ID")
tool_id = "toolshed.g2.bx.psu.edu/repos/bgruening/augustus/augustus/3.4.0+galaxy3"
file = "app/Jollymon.fasta"

upload_response = gi.tools.upload_file(file, history_id)
dataset_id = upload_response.get('outputs')[0].get('id')
upload_job_id = upload_response.get('jobs', [])[0]['id']

i=0
upload_job_state = ""
while upload_job_state != 'ok':  #understand what diff job states exist
    upload_job = gi.jobs.show_job(upload_job_id, True)
    upload_job_state = upload_job.get('state')
    i+=1
    time.sleep(10)
    print(f'you have been waiting for {i*10} seconds')

inputs = {
    'input_genome': {'values': [{'id': dataset_id, 'src': 'hda'}]},
    'noInFrameStop': False,
    'singlestrand': False,
    'utr': False,
    'model': {
        'augustus_mode': 'builtin',
        '__current_case__': 1,
        'organism': 'human'  #what should i set it as?
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

i=0
augustus_job_state = ""
while augustus_job_state != 'ok':  #understand what diff job states exist
    augustus_job = gi.jobs.show_job(augustus_job_id, True)
    augustus_job_state = augustus_job.get('state')
    i+=1
    time.sleep(10)
    print(f'you have been waiting for {i*10} seconds')

augustus_results = gi.jobs.show_job(augustus_job_id, True)
output_dataset_id = augustus_results['outputs']['output']['id'] #change to .get

output_dataset = gi.datasets.download_dataset(output_dataset_id, "galaxy_downloads")