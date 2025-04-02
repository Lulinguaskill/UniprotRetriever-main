import subprocess
from config_manager import get_config


def run_blast(query_file, db_path, output_file, threshold):
    """ run the blast query in the cmd"""
    cpu = get_config('config.yml')['cpu']
    command = [
        "blastp",  # Assuming blastp for protein sequences, change if needed
        "-query", query_file,
        "-db", db_path,
        "-out", output_file,
        "-outfmt", "5",
        "-evalue", str(threshold),
        "-num_threads", str(cpu)]

    # Run the shell command
    try:
        subprocess.run(command, check=True)
        print("BLAST search completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running BLAST: {e}")









