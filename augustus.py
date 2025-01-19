from bioblend.galaxy import GalaxyInstance
from dotenv import load_dotenv
import os

load_dotenv()

gi = GalaxyInstance('https://usegalaxy.org', key=os.getenv("GALAXY_API_KEY"))
history_id = os.getenv("HISTORY_ID")
tool_id = "toolshed.g2.bx.psu.edu/repos/bgruening/augustus/augustus/3.4.0+galaxy3"
file = "app/Jollymon.fasta"

upload_response = gi.tools.upload_file(file, history_id)
dataset_id = upload_response.get('outputs')[0].get('id')
input('ready?')

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
    'outputs': ['protein', 'codingseq', 'start', 'stop', 'cds']  #takes like 10 secs
}

results = gi.tools.run_tool(history_id, tool_id, inputs)
print(results)


# temp = gi.tools.build(tool_id,None,None,history_id)
# print(temp)