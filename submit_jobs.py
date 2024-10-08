import os
import json
import subprocess

# Paths to template and output files
TEMPLATE_FILE = 'job_template.pbs'
PARAMS_FILE = 'job_parameters.json'
OUTPUT_DIR = 'generated_jobs'

# GLOBAL CONSTANT
NCPUS = 5
MEMORY = 8
WALLTIME = '24:00:00'

# Create the output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Read the job parameters from the JSON file
with open(PARAMS_FILE, 'r') as file:
    job_configs = json.load(file)

# Read the PBS template
with open(TEMPLATE_FILE, 'r') as file:
    template_content = file.read()

# Loop over each job configuration and create the job script
for config in job_configs:
    job_script = template_content
    job_name = config['job_name']
    
    # Set default values for ncpus and mem if not provided
    ncpus = config.get('ncpus', NCPUS)  # Default to 4 CPUs
    mem = config.get('mem', MEMORY)   # Default to 8 GB of memory
    walltime = config.get('walltime', WALLTIME)  # Default to 12 hours

    # Replace placeholders with actual values from the configuration
    job_script = job_script.replace('{{job_name}}', job_name)
    job_script = job_script.replace('{{ncpus}}', str(ncpus))
    job_script = job_script.replace('{{mem}}', mem)
    job_script = job_script.replace('{{walltime}}', walltime)
    job_script = job_script.replace('{{script_name}}', config['script_name'])

    # Write the job script to a file
    job_file_path = os.path.join(OUTPUT_DIR, f'{job_name}.pbs')
    with open(job_file_path, 'w') as job_file:
        job_file.write(job_script)

    # Submit the job to the HPC scheduler
    print(f'Submitting job: {job_name}')
    subprocess.run(['qsub', job_file_path])

