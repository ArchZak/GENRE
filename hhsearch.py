from bioblend.galaxy import GalaxyInstance
from dotenv import load_dotenv
import time
import os

load_dotenv()

gi = GalaxyInstance('https://usegalaxy.org', key=os.getenv("GALAXY_API_KEY"))
history_id = os.getenv("HISTORY_ID")
tool_id = "toolshed.g2.bx.psu.edu/repos/guerler/hhsearch/hhsearch/3.2.0+galaxy0"
file = "output_protein.fasta"

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
    'i': {'src': 'hda', 'id': dataset_id},  
    'method': 'hhsearch',  
    'db_source': {
        'db_source_selector': 'indexed',  
        'ffindex': 'pdb70_hmm_2021-03' 
    },
    'e': 0.0000001  
}

hhsearch_run = gi.tools.run_tool(history_id,tool_id,inputs)
hhsearch_job_id = hhsearch_run.get('jobs', [])[0].get('id')

i=0
hhsearch_job_state = ""
while hhsearch_job_state != 'ok':  #understand what diff job states exist
    hhsearch_job = gi.jobs.show_job(hhsearch_job_id, True)
    hhsearch_job_state = hhsearch_job.get('state')
    i+=1
    time.sleep(10)
    print(f'you have been waiting for {i*10} seconds')

hhsearch_results = gi.jobs.show_job(hhsearch_job_id, True)
output_dataset_id = hhsearch_results['outputs']['output']['id'] #change to .get

output_dataset = gi.datasets.download_dataset(output_dataset_id, "galaxy_downloads", )