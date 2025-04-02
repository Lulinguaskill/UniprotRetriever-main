import os
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import time
import threading
from queue import Queue


def blast(query_sequence, result_queue, counter):
    result = NCBIWWW.qblast("blastp", "nr", query_sequence)  # BLAST search from NCBI request
    result_queue.put(result)
    counter.put(1)


def display_runtime(thread, counter, total_sequences):
    start_time = time.time()
    loading = ['   ', '.  ', '.. ', '...']
    loading_index = 0
    while True:
        elapsed_time = time.time() - start_time
        hours = int(elapsed_time // 3600)
        minutes = int((elapsed_time % 3600) // 60)
        seconds = int(elapsed_time % 60)
        print(f"Loading BLAST {counter}/{total_sequences}{loading[loading_index]}  Elapsed time : {hours:02d}:{minutes:02d}:{seconds:02d}", end='\r')
        loading_index = (loading_index + 1) % len(loading)
        time.sleep(1)
        if not thread.is_alive():
            break

                            
def blast_fasta(query_file):   # query_file is the path to your FASTA file .txt
    
    total_sequences = sum(1 for _ in SeqIO.parse(query_file, "fasta"))
    counter = Queue()
    counter.put(0)

    # Step 1: Perform BLAST search for each sequence

    for i, record in enumerate(SeqIO.parse(query_file, "fasta"), start=1):

        query_sequence = record.seq

        result_queue = Queue()
        thread_blast = threading.Thread(target=blast, args=(query_sequence, result_queue, counter,))
        thread_time = threading.Thread(target=display_runtime, args=(thread_blast, i, total_sequences))

        thread_blast.start()
        thread_time.start()

        thread_blast.join()
        thread_time.join()

        result_handle = result_queue.get()
        print(f"\nBLAST {i}/{total_sequences} completed : {record.id}")

        # Step 2: Parse the BLAST result
        blast_record = NCBIXML.read(result_handle)

        # Step 3: Retrieve similar protein sequences
        similar_sequences = []

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                sequence_data = {
                    "title": alignment.title,
                    "length": alignment.length,
                    "e_value": hsp.expect,
                    "sequence": alignment.hsps[0].sbjct
                }
                similar_sequences.append(sequence_data)

        # Step 4: Remove duplicates based on sequence
        unique_sequences = []
        seen_sequences = set()

        for sequence_data in similar_sequences:
            sequence = sequence_data["sequence"]
            if sequence not in seen_sequences:
                unique_sequences.append(sequence_data)
                seen_sequences.add(sequence)

        # Step 5: Save unique sequences to a file
        output_file = f"{record.id}result.txt"
        output_file = output_file.replace("|", "_")

        with open(output_file, "w") as output:
            for seq_data in unique_sequences:
                output.write(f"Title: {seq_data['title']}\n")
                output.write(f"Length: {seq_data['length']}\n")
                output.write(f"E-value: {seq_data['e_value']}\n")
                output.write(f"Sequence: {seq_data['sequence']}\n\n")

        # Step 6: Save unique sequences into a FASTA file to align
        fasta_output_file = f"{record.id}fasta.txt"
        fasta_output_file = fasta_output_file.replace("|", "_")

        with open(fasta_output_file, "w") as output:
            for seq_data in unique_sequences:
                output.write(f">{seq_data['title']}\n")
                output.write(f"{seq_data['sequence']}\n\n")
