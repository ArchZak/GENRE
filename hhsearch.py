from bioblend.galaxy import GalaxyInstance
from dotenv import load_dotenv
import os

load_dotenv()

gi = GalaxyInstance('https://usegalaxy.org', key=os.getenv("GALAXY_API_KEY"))
history_id = os.getenv("HISTORY_ID")
tool_id = "toolshed.g2.bx.psu.edu/repos/guerler/hhsearch/hhsearch/3.2.0+galaxy0"
file = "output_protein.fasta"

upload_response = gi.tools.upload_file(file, history_id) #add stopper to check file status
dataset_id = upload_response.get('outputs')[0].get('id') #40 seconds (slower than usual?)
input('ready?')

inputs = {
    'i': {'src': 'hda', 'id': dataset_id},  
    'method': 'hhsearch',  
    'db_source': {
        'db_source_selector': 'indexed',  
        'ffindex': 'pdb70_hmm_2021-03' 
    },
    'e': 0.001  
}

results = gi.tools.run_tool(history_id,tool_id,inputs) #3.5 minutes (slower than usual?)
print(results) #stopper again for results...


# temp = gi.tools.build(tool_id,None,None,history_id)
# print(temp)