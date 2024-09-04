import logging
import configparser
from pathlib import Path
from dataclasses import dataclass
import csv
import re
import os
import warnings
import time
import shutil
import subprocess
import argparse
import statistics
import requests
from collections import defaultdict
from typing import Optional, List, Set, Tuple, Dict
import warnings

from openmm import *
import openmm as mm
from openmm.app import *
from openmm.unit import *
from openmm.app import PDBFile
import openmm.app as app
from openmm import Platform, LangevinIntegrator, unit
from openmm.app import PDBFile, StateDataReporter
from mdtraj.reporters import XTCReporter, DCDReporter
#import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.topology import guessers
from MDAnalysis.transformations import unwrap

import numpy as np
from Bio import PDB, SeqIO, pairwise2
from Bio.PDB import PDBParser, NeighborSearch, Selection, PDBIO, Select
from Bio.PDB.PDBParser import PDBConstructionWarning
from Bio.Seq import Seq
from io import StringIO
from pdbfixer import PDBFixer
import freesasa
import pandas as pd

# Setup logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger(__name__)
timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
#print(f"[{timestamp}] {message}")

# Load configuration
config = configparser.ConfigParser()
config.read('config.ini')

@dataclass
class Protein:
    def __init__(self, name: str, species: str, fasta_path: Path, pdb_path: Path, files_path: Path, residues: List[List[str]] = None):
        self.name = name
        self.species = species
        self.fasta_path = fasta_path
        self.pdb_path = pdb_path
        self.files_path = files_path
        self.residues = residues if residues is not None else []
class DownloadError(Exception):
    pass

class CSVSummarizer:
    def __init__(self, input_dir: str, crosslinker: str, file:str):
        self.input_dir = input_dir
        self.output_dir = input_dir + "_summarized"
        self.crosslinker = crosslinker
        self.file = file
        self.summary = defaultdict(dict)
        self.summarized_file_header = ["IDx", "Protein1", "Peptide1", "AA1", "Protein2", "Peptide2", "AA2", "Copies" ,"Crosslinker"]
        self.fieldnames, self.fieldnames2 = None, None
        self._parse_files()

    def _parse_files(self):
        os.makedirs(self.output_dir, exist_ok=True)
        for dirName, subdirList, fileList in os.walk(self.input_dir):
            for fname in fileList:
                if self.file is not None and fname != self.file:
                    continue

                elif fname.endswith('.csv'):
                    csv_file = os.path.join(dirName, fname)
                    self.fieldnames, self.fieldnames2 = self._read_first_two_lines(csv_file)

                    if self.fieldnames == self.summarized_file_header or fname.endswith('_summarized.csv'):
                        shutil.copyfile(csv_file, os.path.join(self.output_dir, fname))
                    
                    else:
                        with open(csv_file, 'r') as f:
                            reader = csv.reader(f)
                            logger.info(f"Processing {fname}")
                            self._handle_file(csv_file, reader, self.output_dir, fname, self.crosslinker)

    def _read_first_two_lines(self, csv_file):
        with open(csv_file, 'r') as f:
            # Read the first two lines
            first_line = f.readline().strip()
            second_line = f.readline().strip()
            
            # Use csv.reader to parse these lines
            fieldnames = next(csv.reader([first_line]))
            fieldnames2 = next(csv.reader([second_line]))
            
        return fieldnames, fieldnames2

    def _find_column(self, headers, possible_names):
        for i, name in enumerate(headers):
            if name.lower() in possible_names:
                return i
        return None

    def _find_multiple_columns(self, fieldnames, patterns):
        matching_indices = []
        for i, fieldname in enumerate(fieldnames):
            for pattern in patterns:
                if re.match(pattern, fieldname, re.IGNORECASE):
                    matching_indices.append(i)
        return matching_indices

    def _extract_protein_pairs(self, cell_content, link_type):
        entries = cell_content.split('/')
        protein_pairs = []
        proteins = []

        if link_type == None:
            link_type = "crosslink"     # Default to crosslink if link type is not provided
        
        if "loop" in link_type or "intra" in link_type:
            for entry in entries:
                proteins += re.findall(r'(?:sp|gi|tr)\|[A-Za-z][A-Za-z0-9_.-]+\|([A-Za-z0-9_:-]+)', entry)
                proteins += re.findall(r'([A-Za-z0-9-]+)\(\d+\)', entry)
                proteins += re.findall(r'XP_\d+\.\d+-(\w+)\s*\(\d+\)', entry)
                proteins += re.findall(r'(\w+)-XP_\d+\.\d+\s*\(\d+\)', entry)
            for protein in proteins:
                if (protein.strip("=-"), protein.strip("=-")) not in protein_pairs:
                    protein_pairs.append((protein.strip("=-"), protein.strip("=-")))  
                
        elif "cross" in link_type or "inter" in link_type:
            for entry in entries:
                proteins += re.findall(r'(?:sp|gi|tr)\|[A-Za-z][A-Za-z0-9_.-]+\|([A-Za-z0-9_:-]+)', entry)
                proteins += re.findall(r'([A-Za-z0-9-]+)\(\d+\)', entry)
                proteins += re.findall(r'XP_\d+\.\d+-(\w+)\s*\(\d+\)', entry)
                proteins += re.findall(r'(\w+)-XP_\d+\.\d+\s*\(\d+\)', entry)
                for i in range(0, len(proteins), 2):
                    if i + 1 < len(proteins):
                        prot1 = proteins[i].strip("=-")
                        prot2 = proteins[i + 1].strip("=-")
                        if (prot1, prot2) not in protein_pairs:
                            protein_pairs.append((prot1, prot2))
        
        elif link_type == "mono":
            #logger.info("Can't parse mono links")
            return None

        return protein_pairs

    def _extract_peptides(self, peptide_cell, link_type):
        peptide_pairs = []
        seq1, pos1, seq2, pos2 = None, None, None, None
        peptides = re.findall(r'[a-zA-Z]{2,}', peptide_cell)
        positions =  re.findall(r'\d+', peptide_cell)

        if link_type == None:
            link_type = "crosslink"     # Default to crosslink if link type is not provided

        if "mono" in link_type:
            #logger.info("Can't parse mono links")
            return None

        if len(peptides) == 2:
            seq1 = peptides[0]
            seq2 = peptides[1]

        elif len(peptides) == 1:
            seq1 = peptides[0]
            seq2 = peptides[0]

        pos1, pos2 = positions
        peptide_pairs += [seq1, pos1, seq2, pos2]
        
        return peptide_pairs

    def _separate_search(self, row, fieldnames):
        protein1_column = self._find_column(fieldnames, ['protein1', 'protein_a', 'alpha peptide protein'])
        protein2_column = self._find_column(fieldnames, ['protein2', 'protein_b', 'beta peptide protein'])
        peptide1_column = self._find_column(fieldnames, ['peptide1', 'pepseq1', 'peptide_a'])
        peptide2_column = self._find_column(fieldnames, ['peptide2', 'pepseq2', 'peptide_b'])
        peptide_pos1 = self._find_column(fieldnames, ['linkpos1', 'xl_a', 'pep1_link_pos'])
        peptide_pos2 = self._find_column(fieldnames, ['linkpos2', 'xl_b', 'pep2_link_pos'])

        protein_pairs = []
        peptide_pairs = []

        if peptide1_column is None or peptide2_column is None or peptide_pos1 is None or peptide_pos2 is None:
            peptide_column = row[self._find_column(fieldnames, ['peptide', 'id'])]
            
            peptides = re.findall(r'[a-zA-Z]{2,}', peptide_column)
            positions =  re.findall(r'\d+', peptide_column)

            if len(peptides) == 2:
                pep1 = peptides[0]
                pep2 = peptides[1]
            elif len(peptides) == 1:
                pep1 = peptides[0]
                pep2 = peptides[0]

            peptide_pos1, peptide_pos2 = positions

        else:
            pep1 = row[peptide1_column]
            pep2 = row[peptide2_column]
            peptide_pos1 = row[peptide_pos1]
            peptide_pos2 = row[peptide_pos2]
        
        if row[protein2_column] == "" or row[protein2_column] == "-":
            protein2_column = protein1_column

        proteins = []
        for protein in [row[protein1_column], row[protein2_column]]:
            proteins += re.findall(r'(?:sp|gi|tr)\|[A-Za-z][A-Za-z0-9_.-]+\|([A-Za-z0-9_:-]+)', protein)
            proteins += re.findall(r'([A-Za-z0-9-]+)\(\d+\)', protein)
            proteins += re.findall(r'XP_\d+\.\d+-(\w+)\s*\(\d+\)', protein)
            proteins += re.findall(r'(\w+)-XP_\d+\.\d+\s*\(\d+\)', protein)
        
        if proteins == []:
            for protein in [row[protein1_column], row[protein2_column]]:
                proteins += re.findall(r'([A-Za-z0-9]+)', protein)

        for i in range(0, len(proteins), 2):
            if i + 1 < len(proteins):  
                protein_pairs.append((proteins[i].strip("=-"), proteins[i + 1].strip("=-")))

        peptide_pairs = [pep1, peptide_pos1, pep2, peptide_pos2]
        
        return protein_pairs, peptide_pairs

    def _handle_file(self, csv_file, reader, out_dir, fname, crosslinker_def):
        with open(f"{out_dir}/{fname[:-4]}_summarized.csv", 'w', newline='') as summarized_file:
            writer = csv.writer(summarized_file)
            writer.writerow(self.summarized_file_header)
            fieldnames, fieldnames2 = self._read_first_two_lines(csv_file)

            next(reader)    #skip header
            if fieldnames2[0] == "":
                next(reader)    #skip header
            summarized_row = []
            new_row = []
            protein_pairs = []
            peptide_pairs = []
            dual_header = False

            for row in reader:
                #logger.info(f"Row: {row}")
                if row[0].isdigit():
                    protein_pairs = []
                    peptide_pairs = []

                protein_column = self._find_column(fieldnames, ['protein', 'proteins', 'prot'])
                peptide_column = self._find_column(fieldnames, ['peptide', 'pep', 'seq'])
                if peptide_column is None:
                    peptide_column = self._find_column(fieldnames2, ['peptide'])
                    if peptide_column is not None:
                        dual_header = True
                
                if crosslinker_def is None:
                    crosslinker_column = self._find_column(fieldnames, ['linker', 'crosslinker'])
                    try:
                        crosslinker = row[crosslinker_column]
                    except TypeError:
                        crosslinker = "Unknown"
                else:
                    crosslinker = crosslinker_def    # Overwrite the crosslinker if provided

                try:
                    link_type = row[self._find_column(fieldnames, ['peptide_type', 'protein_type', 'type'])].lower()
                    if "mono" in link_type.lower() or "common" in link_type.lower():
                        #logger.info("Skipping mono links")
                        continue
                except TypeError:
                    link_type = None

                decoy_column = self._find_multiple_columns(fieldnames, [r'^decoy\d*$', r'^Decoy\d*$', r'^isDecoy\d*$', r'^is_decoy\d*$', "target_decoy"])
                if len(decoy_column) != 0:    
                    for index in decoy_column:
                        if row[index].lower() == "true" or row[index].lower() == "yes" or row[index].lower() == "decoy":
                            #logger.info("Decoy found. Skipping...")
                            continue
                
                if dual_header:
                    if row[0] != "":
                        returned_protein = self._extract_protein_pairs(row[protein_column], link_type)
                        if returned_protein not in protein_pairs:
                            protein_pairs += returned_protein
                        continue
                    
                    else:
                        peptide_pairs = self._extract_peptides(row[peptide_column], link_type)

                elif protein_column is None or peptide_column is None:
                    protein_pairs, peptide_pairs = self._separate_search(row, fieldnames)

                else:    
                    if row[peptide_column].isdigit():
                        value = int(row[peptide_column])
                        summarized_row[-1][7] = value
                        protein_pairs = []
                        peptide_pairs = []
                        continue
                    else:
                        protein_pairs += self._extract_protein_pairs(row[protein_column], link_type)
                        peptide_pairs = self._extract_peptides(row[peptide_column], link_type)

                for entry in protein_pairs:
                    protein1, protein2 = entry
                    pep1, aa1, pep2, aa2 = peptide_pairs

                    if "decoy" in protein1.lower() or "decoy" in protein2.lower():
                        #logger.info("Decoy found. Skipping...")
                        continue

                    new_row = [
                        len(summarized_row) + 1,
                        protein1,
                        pep1,
                        aa1,
                        protein2,
                        pep2,
                        aa2,
                        1,  # Number of repeats
                        f"{crosslinker}" 
                    ]

                    index = next((i for i, existing_row in enumerate(summarized_row) if existing_row[1:7] == new_row[1:7]), -1)
                    if index == -1:
                        # If new_row does not exist in summarized_row, append it
                        summarized_row.append(new_row)
                        new_row = []

                    else:
                        # If new_row exists in summarized_row, increment the number of repeats and update the existing row
                        pop_counter = summarized_row[index][7] + 1
                        summarized_row[index][7] = pop_counter
                        new_row = []

            writer.writerows(summarized_row)

class ProteinDownloader:
    def __init__(self, prot1: Protein, prot2: Protein):
        self.prot1 = prot1
        self.prot2 = prot2

    def download_fasta(self):
        for prot in [self.prot1, self.prot2]:
            if not prot.fasta_path.exists():
                protein_name = prot.name
                if ":" in protein_name:
                    protein_id_type, protein_id = protein_name.split(":", 1)
                    if protein_id_type == "TREMBL":
                        old_name = f"{prot.name}_{prot.species}"
                        self._download_fasta_from_uniprot_by_ID(protein_id, prot)
                    elif protein_id_type == "ENSEMBL":
                        logger.info(f"ENSEMBL based ID is not supported. ENSEMBL:{protein_id}.")
                        #self._download_fasta_from_ensembl(protein_id, prot.fasta_path, prot)
                    elif protein_id_type == "SWISS-PROT":
                        self._download_fasta_from_uniprot_by_ID(protein_id, prot)
                    else:
                        logger.info(f"Database search for '{protein_name}' is not supported.")
                
                else:
                    self._download_fasta_from_uniprot_by_name(protein_name, prot.fasta_path, prot.species)
        
        current_directory = os.getcwd()
        current_folder_name = os.path.basename(current_directory)
        desired_folder_name = f"{self.prot1.name}_{self.prot1.species}_{self.prot2.name}_{self.prot2.species}"

        if current_folder_name != desired_folder_name:
            os.chdir("..")
            
            new_folder_name = desired_folder_name
            counter = 1
            while os.path.exists(new_folder_name):
                new_folder_name = f"{desired_folder_name}({counter})"
                counter += 1
                
            os.rename(current_folder_name, new_folder_name)
            os.chdir(new_folder_name)


    def _download_fasta_from_uniprot_by_ID(self, protein_id, prot):
        url = f"https://www.uniprot.org/uniprot/{protein_id}.fasta"
        response = requests.get(url)
        if response.status_code == 200 and response.text.strip():
            lines = response.text.split('\n')
            header = lines.pop(0)
            # Extract the name and species from the header
            parts = header.split()
            protein_identifier = parts[0].split('|')[2]
            protein_name, species = protein_identifier.split('_')
            
            sequence = ''.join(lines)
            output_file = Path(f"pdb_fa_files/{protein_name}_{species}.fa")
            with open(output_file, "w", encoding='utf-8') as file:
                file.write(f">{header}\n")
                file.write(sequence + '\n')
            logger.info(f"FASTA file for {protein_id} downloaded successfully from Swiss-Prot.")
            prot.name = protein_name
            prot.species = species 
            prot.fasta_path = output_file
            prot.pdb_path = Path(f"pdb_fa_files/{protein_name}_{species}.pdb")

            logger.info(f"ID {protein_id} updated to {protein_name} and species {species}.")
        else:
            error_message = f"Failed to download FASTA file by ID for {protein_id} from Swiss-Prot. "
            if response.status_code != 200:
                error_message += f"Status code: {response.status_code}"
            else:
                error_message += "File could not be downloaded (blank page)."
            logger.error(error_message)


    def _download_fasta_from_uniprot_by_name(self, protein_name, output_file, species):
        url = f"https://www.uniprot.org/uniprot/{protein_name}_{species}.fasta"
        response = requests.get(url)
        if response.status_code == 200 and response.text.strip():
            lines = response.text.split('\n')
            header = lines.pop(0)
            sequence = ''.join(lines)
            with open(output_file, "w", encoding='utf-8') as file:
                file.write(f">{header}\n")
                file.write(sequence + '\n')
            logger.info(f"FASTA file for {protein_name} downloaded successfully.")
        else:
            error_message = f"Failed to download FASTA file by name for {protein_name}_{species}. "
            if response.status_code != 200:
                error_message += f"Status code: {response.status_code}"
            else:
                error_message += "File could not be downloaded (blank page)."
            logger.error(error_message)
            # Niche case: I3LJZ9_PIG doesnt have a uniprot entry, though it can be found
            # here: https://rest.uniprot.org/unisave/I3LJZ9?format=txt&versions=1
            # This is a case where the protein is no longer being annotated by uniprot
            # This specific case has a PDB entry in alphafold, but it cannot be relied upon
            # Currently, we just throw an error, as this is likely contamination from MS-MS
            # and not important for our purposes

    def download_pdb(self, download_from: str):
        for prot in [self.prot1, self.prot2]:
            if not prot.pdb_path.exists():
                if download_from.upper() == "AF":
                    self._download_pdb_from_alphafold(prot, None)
                elif download_from.upper() == "SWISS":
                    self._download_pdb_from_swissmodel(prot.name, prot.species, prot, None)

    def _download_pdb_from_alphafold(self, prot, coming_from):
        attempt = "AF"
        get_accession_url = f"https://rest.uniprot.org/uniprotkb/search?query={prot.name}_{prot.species}&format=json"
        
        response = requests.get(get_accession_url)
        if response.status_code != 200:
            logger.error(f"Failed to query UniProt for {prot.name}_{prot.species}. Status code: {response.status_code}")
            logger.info("Attempting to download PDB file from SwissModel...")
            if coming_from == "SM":
                self._download_pdb_file_from_uniprot(prot.name, prot.species, prot)
            else:
                self._download_pdb_from_swissmodel(prot.name, prot.species, prot, attempt)
            return

        data = response.json()
        results = data.get('results', [])
        if not results:
            logger.error(f"No results found for {prot.name}_{prot.species}.")
            logger.info("Attempting to download PDB file from SwissModel...")
            if coming_from == "SM":
                self._download_pdb_file_from_uniprot(prot.name, prot.species, prot)
            else:
                self._download_pdb_from_swissmodel(prot.name, prot.species, prot, attempt)
            return

        accession = results[0].get('primaryAccession')
        if not accession:
            logger.error(f"No accession number found for {prot.name}_{prot.species}.")
            logger.info("Attempting to download PDB file from SwissModel...")
            if coming_from == "SM":
                self._download_pdb_file_from_uniprot(prot.name, prot.species, prot)
            else:
                self._download_pdb_from_swissmodel(prot.name, prot.species, prot, attempt)
            return
    
        # Construct the URL for AlphaFold API
        url = f"https://alphafold.ebi.ac.uk/files/AF-{accession}-F1-model_v4.pdb"
        response = requests.get(url)

        if response.status_code != 200:
            logger.error(f"Failed to download AlphaFold model for {accession}. Status code: {response.status_code}")
            if coming_from == "SM":
                self._download_pdb_file_from_uniprot(prot.name, prot.species, prot)
            else:
                self._download_pdb_from_swissmodel(prot.name, prot.species, prot, attempt)
            return

        pdb_filename = f"pdb_fa_files/AF-{accession}-F1-model_v4.pdb"
        with open(pdb_filename, 'wb') as file:
            file.write(response.content)
        
        logger.info(f"AlphaFold model for {accession} downloaded successfully as {pdb_filename}")

        if not prot.fasta_path.exists():
            logger.info(f"FASTA file for {prot.name}_{prot.species} not found. Re-attempting to download...")
            self._download_fasta_from_uniprot_by_ID(accession, prot)
        
        self._extract_matching_chain(pdb_filename, f"pdb_fa_files/{prot.name}_{prot.species}.pdb", prot.name, prot.species, prot)

    def _download_pdb_from_swissmodel(self, protein_name, species, prot, coming_from):
        attempt = "SM"
        # Query SwissModel
        url = f"https://swissmodel.expasy.org/repository/uniprot/{protein_name}_{species}.json"
        response = requests.get(url)

        if response.status_code != 200:
            logger.error(f"Failed to query SwissModel for {protein_name}. Status code: {response.status_code}")
            if coming_from == "AF":
                self._download_pdb_file_from_uniprot(protein_name, species, prot)
            else:
                self._download_pdb_from_alphafold(prot, attempt)
            return

        data = response.json()
        models = data.get('result', {}).get('structures', [])

        if not models:
            logger.error(f"No models found for {protein_name}_{species}.")
            if coming_from == "AF":
                self._download_pdb_file_from_uniprot(protein_name, species, prot)
            else:
                self._download_pdb_from_alphafold(prot, attempt)
            return

        # Sort models by sequence identity percentage
        sorted_models = sorted(models, key=lambda x: x.get('identity', 0), reverse=True)

        for model in sorted_models:
            prot_pdb_url = model.get('coordinates')
            pdb_template_name = model.get('template')

            if prot_pdb_url:
                # Download the PDB file
                pdb_response = requests.get(prot_pdb_url)

                if pdb_response.status_code == 200:
                    with open(f"pdb_fa_files/{pdb_template_name}.pdb", 'wb') as f:
                        f.write(pdb_response.content)
                    self._extract_matching_chain(f"pdb_fa_files/{pdb_template_name}.pdb", f"pdb_fa_files/{protein_name}_{species}.pdb", protein_name, species, prot)
                    return
                else:
                    logger.error(f"Failed to download PDB file from {prot_pdb_url}. Status code: {pdb_response.status_code}")
                    if coming_from == "AF":
                        self._download_pdb_file_from_uniprot(protein_name, species, prot)
                    else:
                        self._download_pdb_from_alphafold(prot, attempt)
                    return
            else:
                logger.error(f"No PDB URL found for the model with template {pdb_template_name}.")
                if coming_from == "AF":
                    self._download_pdb_file_from_uniprot(protein_name, species, prot)
                else:
                    self._download_pdb_from_alphafold(prot, attempt)
                return

        logger.error(f"All models failed to provide a valid PDB URL for {protein_name}. Please provide your own PDB file under the folder and name <pdb_fa_files/{protein_name}_{species}.pdb>.")

    def _download_pdb_file_from_uniprot(self, protein_name, species, prot):
        uniprot_id = ""
        try:
            with open(prot.fasta_path, 'r') as file:
                for line in file:
                    if line.startswith('>'):
                        # Assuming the UniProt ID is the first word in the header line after '>'
                        uniprot_id = line.split('|')[1]
        except FileNotFoundError:
            return
        
        pdb_url = f"https://files.rcsb.org/download/{uniprot_id}.pdb"
        response = requests.get(pdb_url)
        if response.status_code == 200:
            with open(f"pdb_fa_files/{uniprot_id}.pdb", 'wb') as file:
                file.write(response.content)
            logger.info(f"PDB file {uniprot_id}.pdb downloaded successfully.")
            self._extract_matching_chain(f"pdb_fa_files/{uniprot_id}.pdb", f"pdb_fa_files/{protein_name}_{species}.pdb", protein_name, species, prot)
        
        else:
            logger.info(f"All other attempts failed to fetch PDB file for {protein_name}_{species}, UniProt ID: {uniprot_id}. Please provide your own PDB file under the name <{protein_name}_{species}.pdb>.")

    def _extract_matching_chain(self, input_file, output_file, protein_name, species, prot):
        # Fetch the sequence of the target protein
        target_sequence = self._get_protein_sequence(protein_name, species)
        
        if not target_sequence:
            logger.error(f"Failed to fetch sequence for {protein_name}_{species}")
            return

        # Parse the PDB file
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", input_file)

        # Find the chain that matches the beginning of the target sequence
        matching_chain = None

        for chain in structure.get_chains():
            chain_sequence = self._get_chain_sequence(chain)
            
            # Check if the first 30 residues of the chain sequence are in the target sequence
            if len(chain_sequence) >= 30 and chain_sequence[:30] in target_sequence:
                matching_chain = chain.id
                break

        if matching_chain is None:
            logger.error(f"No matching chain found for {protein_name}_{species}")
            return

        # Define a selector for the matching chain
        class ChainSelect(Select):
            def __init__(self, chain_id):
                self.chain_id = chain_id

            def accept_chain(self, chain):
                return chain.id == self.chain_id

        # Save only the selected chain
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_file, ChainSelect(matching_chain))

        logger.info(f"Extracted matching chain {matching_chain} from {input_file} to {output_file}")
        prot.pdb_path = Path(f"{output_file}")


    def _get_chain_sequence(self, chain):
        """Extract the amino acid sequence from a PDB chain."""
        return ''.join(
            self._convert_to_1_letter_code(residue.resname)
            for residue in chain
            if residue.id[0] == ' '  # Check if it's a standard amino acid
        )

    def _get_protein_sequence(self, protein_name, species):
        # Fetch the protein sequence from UniProt
        uniprot_id = f"{protein_name}_{species}"
        url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
        response = requests.get(url)
        
        if response.status_code != 200:
            logger.error(f"Failed to fetch sequence for {uniprot_id}")
            return None

        # Parse the FASTA file
        fasta_sequences = SeqIO.parse(StringIO(response.text), "fasta")
        for fasta in fasta_sequences:
            return str(fasta.seq)

        return None
    
    def _parse_fasta(self, fasta_path):
        with open(fasta_path, "r") as fasta_file:
            records = list(SeqIO.parse(fasta_file, "fasta"))
            if len(records) != 1:
                raise ValueError("FASTA file should contain exactly one sequence.")
            return records[0].seq

    def _extract_pdb_sequence(self, pdb_path):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_path)
        seq = []
        res_ids = []

        for model in structure:
            for chain in model:
                for residue in chain:
                    residue_id = residue.id
                    if residue_id[0] == ' ':
                        one_letter_code = self._convert_to_1_letter_code(residue.resname)
                        seq.append(one_letter_code)
                        res_ids.append(residue_id)

        return Seq("".join(seq)), res_ids, structure
    
    def _renumber_fixed_pdb(self, fixed_pdb_path, fasta_seq, prot):
        fixed_seq, fixed_res_ids, structure = self._extract_pdb_sequence(fixed_pdb_path)

        # Get the first 30 residues from the fixed PDB sequence
        pdb_start_seq = str(fixed_seq[:30])

        # Find the start position of the PDB sequence in the FASTA sequence
        start_pos = str(fasta_seq).find(pdb_start_seq)  # Add 1 to convert from 0-based to 1-based index

        if start_pos == -1:
            raise ValueError("Could not find the PDB sequence in the FASTA file.")

        logger.info(f"Found PDB sequence at position {start_pos} in the FASTA sequence.")
        
        # Calculate the end position based on the sequence length
        end_pos = start_pos + len(fixed_seq)

        # Renumber residues starting from end_pos (converting to 1-based index)
        new_residue_number = end_pos

        for model in structure:
            for chain in model:
                # Get the residues in reverse order
                residues = list(chain.get_residues())[::-1] # Reverse to start from the end to avoid ID errors
                for residue in residues:
                    # Update residue ID tuple (chain identifier, residue number, insertion code)
                    new_id = (residue.id[0], new_residue_number, residue.id[2])
                    residue.id = new_id
                    new_residue_number -= 1  # Decrement for the next residue

        # Write the modified structure to a new PDB file
        io = PDBIO()
        io.set_structure(structure)
        io.save(fixed_pdb_path)

    def process_pdb_with_residue_validation(self, fix_partial_prot):
        logger.info("Populating files with missing residues and hydrogens. This may take a while...")
        for prot in [self.prot1, self.prot2]:
            fasta_path = prot.fasta_path
            fixed_pdb_path = f"{str(prot.pdb_path)[:-4]}_fixed.pdb"

            if os.path.exists(fixed_pdb_path):
                logger.info(f"Fixed PDB file already exists for {prot.name}. Skipping fixing process.")
                self._replace_cd1_with_cd(fixed_pdb_path)
                logger.info(f"Replaced CD1 with CD for ILE residues in {prot.name}")
                continue

            self._add_fixer_tag(prot)

            fixer = PDBFixer(filename=f"{str(prot.pdb_path)}")
            fixer.findMissingResidues()
            #logger.info(f"Missing residues: {fixer.missingResidues}")
            if fix_partial_prot:
                chains = list(fixer.topology.chains())
                keys_to_delete = []
                for key in fixer.missingResidues.keys():
                    chain = chains[key[0]]
                    if key[1] == 0 or key[1] == len(list(chain.residues())):
                        keys_to_delete.append(key)

                for key in keys_to_delete:
                    del fixer.missingResidues[key]

            fixer.findNonstandardResidues()
            fixer.replaceNonstandardResidues()
            fixer.removeHeterogens(True)
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(7.0)

            PDBFile.writeFile(fixer.topology, fixer.positions, open(fixed_pdb_path, 'w'))
            logger.info(f"Fixed PDB file saved for {prot.name}")
            fasta_seq = self._parse_fasta(fasta_path)
            self._renumber_fixed_pdb(fixed_pdb_path, fasta_seq, prot)
            logger.info(f"Corrected residue numbering for {prot.name}")

            self._replace_cd1_with_cd(fixed_pdb_path)
            logger.info(f"Replaced CD1 with CD for ILE residues in {prot.name}")

            if not fix_partial_prot:
                # Attempt a short simulation to fix clashes
                self._fix_clashes(fixed_pdb_path)

            prot.pdb_path = fixed_pdb_path

    def _fix_clashes(self, fixed_pdb_path):
        pdb = app.PDBFile(fixed_pdb_path)
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        system = forcefield.createSystem(pdb.topology, nonbondedCutoff=1.0* unit.nanometer)
        integrator = LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 0.002*unit.picoseconds)
        simulation = app.Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)

        # Minimize the energy
        logger.info('Attempting to fix clashes through energy minimization...')
        simulation.minimizeEnergy(maxIterations=50000, tolerance=0.1)
        energies = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        logger.info(f"System minimized at {energies}")

        # Save the minimized structure
        positions = simulation.context.getState(getPositions=True).getPositions()
        with open(fixed_pdb_path, 'w') as f:
            app.PDBFile.writeFile(simulation.topology, positions, f)
        logger.info(f'Saved simulation endpoint to {fixed_pdb_path}')

          
    def _add_fixer_tag(self, prot):
        with open(prot.pdb_path, 'r') as f:
            pdb_lines = f.readlines()

        # Find the index where the SEQRES section starts
        seqres_index = next((i for i, line in enumerate(pdb_lines) if line.startswith('SEQRES')), None)

        # If SEQRES section is not found, set the index to start adding missing residues
        if seqres_index is None:
            seqres_index = 0
        
        # Create lines for SEQRES section
        seqres_lines = []
        fasta_residues = self._get_fasta_residues(prot.fasta_path)

        # Sort the residues by residue number before adding to SEQRES
        sorted_fasta_residues = sorted(fasta_residues, key=lambda x: x[1])
        
        # Initialize a counter for SEQRES lines
        seqres_counter = 1
        for i, (residue_name, residue_number) in enumerate(sorted_fasta_residues):
            # Check if a new SEQRES line needs to be started
            if i % 13 == 0:
                seqres_line = f'SEQRES{seqres_counter:>4} A {len(sorted_fasta_residues):>4}  '
                seqres_counter += 1
            seqres_line += f'{residue_name:<3} '

            # If this is the last residue for the current line or the last residue overall, append the line to seqres_lines
            if (i + 1) % 13 == 0 or i == len(sorted_fasta_residues) - 1:
                seqres_line = seqres_line.rstrip() + '\n'
                seqres_lines.append(seqres_line)

        # Change character at index 22 to 'N' for each line
        pdb_lines = [line[:21] + 'A' + line[22:] for line in pdb_lines if line.startswith("ATOM")]

        #print(f"Seqres lines: {seqres_lines}")
        # Write the modified PDB file
        with open(f"{str(prot.pdb_path)[:-4]}_fixed.pdb", 'w') as f:
            # Write SEQRES section
            f.writelines(seqres_lines)
            # Write original PDB lines
            f.writelines(pdb_lines)

    def _get_fasta_residues(self, filename):
        fasta_residues = set()
        with open(filename, 'r') as fasta_file:
            sequence = ""
            for line in fasta_file:
                line = line.strip()
                if line.startswith('>'):
                    continue  # Ignore header lines
                sequence += line  # Concatenate sequence lines
            # Create tuples of residue and its position
            for i, residue in enumerate(sequence, start=1):
                residue = self._convert_to_3_letter_code(residue)  # Convert to 3-letter code
                fasta_residues.add((residue, i))
        # Sort the set based on position
        fasta_AAs = sorted(fasta_residues, key=lambda x: x[1])
        return set(fasta_AAs)
    
    def _convert_to_1_letter_code(self, residue_name):
        # Dictionary mapping 3-letter amino acid names to 1-letter codes
        residue_codes = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }
        # Convert residue name to uppercase and use the dictionary to get the 1-letter code
        return residue_codes.get(residue_name.upper(), 'X')
    
    def _convert_to_3_letter_code(self, residue_name):
        # Dictionary mapping 3-letter amino acid names to 1-letter codes
        residue_codes = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
        'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
        'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
        'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
        }
        # Convert residue name to uppercase and use the dictionary to get the 1-letter code
        return residue_codes.get(residue_name.upper(), 'X')  # Return 'X' if the code is not found
    
    def _replace_cd1_with_cd(self, pdb_file):
        with open(pdb_file, 'r') as f:
            lines = f.readlines()

        modified_lines = []
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                residue_name = line[17:20]
                if atom_name == "CD1" and residue_name == "ILE":
                    modified_line = line[:12] + " CD  ILE" + line[20:]
                    modified_lines.append(modified_line)
                else:
                    modified_lines.append(line)
            else:
                modified_lines.append(line)

        with open(pdb_file, "w") as f:
            f.writelines(modified_lines)

class OpenMMSimulation:
    def __init__(self, prot1: Protein, prot2: Protein):
        self.prot1 = prot1
        self.prot2 = prot2
        self.full_data = []
    
    def run_simulation(self, sim_type: str):
        prot1 = self.prot1.name
        prot1_residues = self.prot1.residues
        prot2 = self.prot2.name
        prot2_residues = self.prot2.residues

        os.makedirs("OpenMM_sim", exist_ok=True)

        if sim_type == "Single":
            work_dir = f"{prot1}_Mono_sim"
            file_name = f"{prot1}"
            shutil.copy(self.prot1.pdb_path, f"OpenMM_sim/{file_name}.pdb")
        elif sim_type == "Dimer":
            work_dir = f"{prot1}_{prot2}_{'Homo' if prot1 == prot2 else 'Hetero'}dimer_sim"
            file_name = f"{prot1}_{prot2}"
            shutil.copy("Haddock/model_1.pdb", f"OpenMM_sim/{file_name}.pdb")
        
        os.chdir("OpenMM_sim")
        os.makedirs(work_dir, exist_ok=True)
        os.makedirs(f"{work_dir}/Data", exist_ok=True)

        # Load the PDB file
        pdb = PDBFile(f"{file_name}.pdb")

        # Set up force field
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        platform = Platform.getPlatformByName("OpenCL")
        properties = {'Precision': 'mixed'}

        # Create system
        system = forcefield.createSystem(
            pdb.topology,
            nonbondedCutoff=1* unit.nanometer,
            constraints=HBonds
        )

        # Set up simulation parameters
        temperature = 300*unit.kelvin
        #pressure = 1*unit.bar
        friction_coeff = 1/unit.picosecond
        timestep = 0.5*unit.femtoseconds

        # Create integrator
        integrator = LangevinMiddleIntegrator(temperature, friction_coeff, timestep)

        if not os.path.exists(f"{work_dir}/equi.state") or not os.path.exists(f"{work_dir}/equi.chk"):
            # Set up simulation
            simulation = Simulation(pdb.topology, system, integrator, platform, properties)
            simulation.context.setPositions(pdb.positions)

            # Minimize energy
            logger.info("Minimizing energy...")
            simulation.minimizeEnergy()
            energies = simulation.context.getState(getEnergy=True).getPotentialEnergy()
            logger.info(f"System minimized at {energies}")

            # Equilibrate
            logger.info("Equilibrating...")
            simulation.context.setVelocitiesToTemperature(temperature)
            simulation.step(1000)  # 20 ps

            # Production run
            logger.info("Starting production run...")
            simulation.reporters.append(StateDataReporter(f"{work_dir}/prod.log", 100,
                                                        step=True,
                                                        time=True,
                                                        potentialEnergy=True,
                                                        temperature=True,
                                                        progress=True,
                                                        remainingTime=True,
                                                        speed=True,
                                                        totalSteps=50000,
                                                        separator='\t'))

            simulation.reporters.append(DCDReporter(f'{work_dir}/trajectory.dcd', 100))
            xtc_reporter = XTCReporter(f'{work_dir}/prod.xtc', 100)
            simulation.reporters.append(xtc_reporter)

            simulation.saveState(f"{work_dir}/equi.state")
            simulation.saveCheckpoint(f"{work_dir}/equi.chk")

        else:
            logger.info("Previous minimized state found. Loading state...")
            simulation = Simulation(pdb.topology, system, integrator)
            simulation.loadState(f"{work_dir}/equi.state")

        # Run production simulation
        sim_accepted_residues = []
        sim_rejected_residues = []
        step_limit = 5000

        residue_pairs = zip(prot1_residues, prot2_residues)
        
        for prot1_residue, prot2_residue in residue_pairs:
            simulation.loadCheckpoint(f"{work_dir}/equi.chk")
            idx, aa1, peptide1, copies, crosslinker, max_dist, original_file_path = prot1_residue
            idx, aa2, peptide2, copies, crosslinker, max_dist, original_file_path = prot2_residue

            if aa1 == aa2 and sim_type == "Single":
                logger.info(f"Same residue pair: {aa1} and {aa2}. Logging as rejected. Retrying as Homodimer case.")
                sim_rejected_residues.append([idx, prot1, peptide1, aa1, prot2, peptide2, aa2, copies, crosslinker, "0", max_dist, "Same residue: Likely a homodimer case", "N/A", original_file_path])
                continue

            logger.info(f"Working on residues: {aa1} and {aa2}")

            atom1, atom2 = self._extract_highest_index_atoms_from_pdb(pdb.topology, aa1, aa2, sim_type)
            distance = self._compute_distance(simulation, atom1, atom2, sim_type)
            state = simulation.context.getState(getEnergy=True).getPotentialEnergy() / unit.kilocalories_per_mole
            logger.info(f"Starting distance: {distance:.3f} nm")
            logger.info(f"Target distance set to: {max_dist/10} nm")

            # Harmonic force parameters
            force_constant = 1500.0 #* copy_modifier  # kJ/mol/nm^2
            equilibrium_distance = max_dist/10  # nm
            harmonic_force = mm.CustomBondForce("k*(r-r0)^2")
            harmonic_force.addGlobalParameter("k", force_constant)
            harmonic_force.addGlobalParameter("r0", equilibrium_distance)
            new_force_id = simulation.system.addForce(harmonic_force)

            if distance > (max_dist / 10):
                steps_taken = 0
                # Add harmonic restraint
                harmonic_force.addBond(atom1, atom2, [force_constant, equilibrium_distance])
                
                while distance > (max_dist / 10) and steps_taken < step_limit and state < 0:
                    # Run simulation
                    simulation.step(100)  # 20 ps
                    steps_taken += 100
                    distance = self._compute_distance(simulation, atom1, atom2, sim_type)

                    # Checking the state of the simulation
                    state = simulation.context.getState(getEnergy=True).getPotentialEnergy() / unit.kilocalories_per_mole
                    if steps_taken % 1000 == 0:
                        logger.info(f"Step: {steps_taken:.4f}, Distance: {distance:.3f} nm, Potential Energy: {state:.2f}")

                # Simulation ended - Remove restraint
                simulation.system.removeForce(new_force_id)


                if steps_taken <= step_limit and distance <= (max_dist / 10):
                    sim_accepted_residues.append([idx, prot1, peptide1, aa1, prot2, peptide2, aa2, copies, crosslinker, distance*10, max_dist, f"Accepted after {steps_taken} steps", state, original_file_path])
                    
                    with open(f"{work_dir}/Data/idx_{idx}_{self.prot1.name}{'_' + self.prot2.name if sim_type == 'Dimer' else ''}_residues_{aa1}_{aa2}.pdb", 'w') as file:
                        PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), file)
                    
                    logger.info(f"Residues {aa1} and {aa2} reached acceptable distances. PDB saved at {work_dir}/Data.")
                    continue
                
                else:
                    sim_rejected_residues.append([idx, prot1, peptide1, aa1, prot2, peptide2, aa2, copies, crosslinker, distance*10, max_dist, "Reached max steps of 5000", state, original_file_path])
                    logger.info(f"Step limit reached for residues {aa1} and {aa2}.")
                    continue
            
            else:
                logger.info(f"Initial distance is already under max distance.")
                sim_accepted_residues.append([idx, prot1, peptide1, aa1, prot2, peptide2, aa2, copies, crosslinker, distance*10, max_dist, "Pre-sim accepted", state, original_file_path])

        # Save final structure and results
        self._unwrap_trajectory(f"{file_name}.pdb", f"{work_dir}/trajectory.dcd", f"{work_dir}/trajectory_unwrapped.dcd", sim_type)
        self._save_final_structure(simulation, work_dir)
        self._save_results(sim_accepted_residues, sim_rejected_residues, work_dir, sim_type)

    def _extract_highest_index_atoms_from_pdb(self, topology, aminoacid_number1, aminoacid_number2, sim_type):
        highest_atom1 = None
        highest_atom2 = None
        highest_index1 = -1
        highest_index2 = -1

        for atom in topology.atoms():
            if sim_type == "Single":
                if atom.residue.index + 1 == aminoacid_number1:
                    if atom.index > highest_index1:
                        highest_index1 = atom.index
                        highest_atom1 = atom.index
                elif atom.residue.index + 1 == aminoacid_number2:
                    if atom.index > highest_index2:
                        highest_index2 = atom.index
                        highest_atom2 = atom.index
            elif sim_type == "Dimer":
                if atom.residue.chain.index == 0 and atom.residue.index + 1 == aminoacid_number1:
                    if atom.index > highest_index1:
                        highest_index1 = atom.index
                        highest_atom1 = atom.index
                elif atom.residue.chain.index == 1 and atom.residue.index + 1 == aminoacid_number2:
                    if atom.index > highest_index2:
                        highest_index2 = atom.index
                        highest_atom2 = atom.index

        if highest_atom1 is None or highest_atom2 is None:
            raise ValueError("One or both of the specified amino acids were not found in the topology.")

        return highest_atom1, highest_atom2

    def _unwrap_trajectory(self, pdb_file, input_dcd, output_dcd, sim_type):
        if sim_type == "Single" or not os.path.exists(input_dcd):
            return
        u = mda.Universe(pdb_file, input_dcd)
        # Guess bonds if they're not present
        if not hasattr(u.atoms, 'bonds') or len(u.atoms.bonds) == 0:
            guessers.guess_bonds(u.atoms, u.atoms.positions)
        all_atoms = u.select_atoms("all")

        # Function to fix jumps
        def fix_jumps(coords, box):
            for i in range(1, len(coords)):
                diff = coords[i] - coords[i-1]
                for dim in range(3):
                    if abs(diff[dim]) > box[dim]/2:
                        coords[i, dim] -= np.sign(diff[dim]) * box[dim]
            return coords

        # Create unwrap transformation
        unwrapper = unwrap(all_atoms)

        with mda.Writer(output_dcd, all_atoms.n_atoms) as W:
            # Store coordinates of the first frame
            prev_coords = all_atoms.positions.copy()
            
            for ts in u.trajectory:
                # Apply unwrapping
                unwrapper(ts)
                
                # Fix any remaining jumps
                all_atoms.positions = fix_jumps(np.vstack((prev_coords, all_atoms.positions)), ts.dimensions[:3])[-all_atoms.n_atoms:]
                
                # Update previous coordinates
                prev_coords = all_atoms.positions.copy()
                
                # Write the frame
                W.write(all_atoms)

        logger.info("Unwrapping completed. Check the 'unwrapped_fixed.dcd' file.")

    def _compute_distance(self, simulation, atom_index1, atom_index2, sim_type: str):
        positions = simulation.context.getState(getPositions=True).getPositions()
        topology = simulation.topology

        if sim_type == "Single":
            pos1 = positions[atom_index1]
            pos2 = positions[atom_index2]
        elif sim_type == "Dimer":
            chain_a_atoms = [atom.index for atom in topology.atoms() if atom.residue.chain.id == "1"]
            chain_b_atoms = [atom.index for atom in topology.atoms() if atom.residue.chain.id != "1"]      
            if atom_index1 in chain_a_atoms and atom_index2 in chain_b_atoms:
                pos1 = positions[atom_index1]
                pos2 = positions[atom_index2]

            else:
                raise ValueError(f"For Dimer simulations, atoms must be in different chains. Atom1: {atom_index1}, Atom2: {atom_index2}")
        else:
            raise ValueError(f"Invalid sim_type: {sim_type}. Must be 'Single' or 'Dimer'.")

        distance_qt = np.linalg.norm(np.array(pos1) - np.array(pos2))
        distance = distance_qt.value_in_unit(unit.nanometer)
        return distance

    def _save_final_structure(self, simulation, work_dir):
        output_pdb_filename = f'{work_dir}/{self.prot1.name}_{self.prot2.name}_post_sim.pdb'            
        new_sim_state = f'{work_dir}/post_sim_sim.state'
        new_sim_checkpoint = f'{work_dir}/post_sim_sim.chk'
        
        with open(output_pdb_filename, 'w') as file:
            PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), file)
            logger.info(f"Final structure saved to {output_pdb_filename}")
        simulation.saveState(new_sim_state)
        simulation.saveCheckpoint(new_sim_checkpoint)

    def _save_results(self, sim_accepted_residues, sim_rejected_residues, work_dir, sim_type):
        prot2_name = "_" + self.prot2.name if sim_type == "Dimer" else "_Single"
        
        # Write accepted residues to CSV
        with open(f"{work_dir}/CG_Simulation_accepted_residues_{self.prot1.name}{prot2_name}.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["IDx", "Protein1", "Peptide1", "AA1", "Protein2", "Peptide2", "AA2", "Copies", "Crosslinker", "Final Distance", "Max Dist", "Notes", "Potential Energy", "Original File Path"])
            writer.writerows(sim_accepted_residues)
        
        # Write rejected residues to CSV
        with open(f"{work_dir}/CG_Simulation_rejected_residues_{self.prot1.name}{prot2_name}.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["IDx", "Protein1", "Peptide1", "AA1", "Protein2", "Peptide2", "AA2", "Copies", "Crosslinker", "Final Distance", "Max Dist", "Notes", "Potential Energy", "Original File Path"])
            writer.writerows(sim_rejected_residues)

        os.chdir("..")
            
        buried_residues_file = f"pdb_fa_files/{self.prot1.name}_{self.prot2.name}_buried_residues.csv"
        if os.path.exists(buried_residues_file):
            shutil.copy(buried_residues_file, f"CG_OpenMM/{work_dir}/{self.prot1.name}_{prot2_name}_buried_residues.csv")
        else:
            logger.info(f"No buried residues file found for {self.prot1.name}{prot2_name}.")
        
        logger.info(f"Results saved to {work_dir}")

class HaddockSimulation:
    def __init__(self, prot1: Protein, prot2: Protein):
        self.prot1 = prot1
        self.prot2 = prot2
        self.full_data = []

    def generate_haddock_files(self):
        logger.info(f"{self.prot1.name} and {self.prot2.name} HADDOCK simulation started.")
        logger.info(f"PDB files: {self.prot1.pdb_path} and {self.prot2.pdb_path}")
        logger.info(f"Fasta files: {self.prot1.fasta_path} and {self.prot2.fasta_path}")
        # Generate active and passive residues based on crosslink data
        active_residues = self._get_active_residues()
        self._get_passive_residues(active_residues)

        self._generate_param_file(self.prot1.name, self.prot2.name, self.prot1.species, self.prot2.species)

        self._transfer_files()

    def _get_active_residues(self):
        restraint_lines = []
        copies_all_data = []
        filtered_lines = []
        
        for residues1, residues2 in zip(self.prot1.residues, self.prot2.residues):
                _, aa1, _, copies, _, max_dist, _ = residues1
                _, aa2, _, copies, _, max_dist, _ = residues2
                restraint_lines.append([aa1, aa2, max_dist, int(copies)])
                copies_all_data.append(int(copies))

        # Filter the restraints based on the average copies score
        sum_copies = sum(copies_all_data)
        copies_score = [entry/sum_copies for entry in copies_all_data]
        average_copies_score = statistics.mean(copies_score)
        for score, entry in zip(copies_score, restraint_lines):
            if score >= average_copies_score:
                filtered_lines.append(entry)

        with open(f"{self.prot1.name}_{self.prot2.name}_active_restraints.tbl", 'w') as file:
            for line in filtered_lines:
                file.write(f"assign (resid {line[0]} and name CA and segid A) (resid {line[1]} and name CA and segid B) {line[2]} {line[2]} 0\n")

        logger.info(f"Active restraints written to file.")
        return filtered_lines

    def _get_passive_residues(self, active_residues):
        prot1 = self.prot1
        prot2 = self.prot2
        distance_threshold = 8

        with open(f"{prot1.name}_{prot2.name}_passive_restraints.tbl", 'w') as file:

            for residue1, residue2, max_dist, _ in active_residues:
                pass_res1 = self._find_nearby_residues(residue1, distance_threshold, prot1.pdb_path)
                pass_res2 =self._find_nearby_residues(residue2, distance_threshold, prot2.pdb_path)

                file.write("! HADDOCK AIR restraints for 1st partner\n")
                file.write("!\n")

                # Generate the ambiguous restraints
                file.write(f"assign (resid {residue1} and segid A)\n")
                file.write("(\n")

                for res2 in pass_res2:
                    file.write(f"    (resid {res2} and segid B)\n")
                    # If it's not the last residue, write "or"
                    if res2 != pass_res2[-1]:
                        file.write("    or\n")

                # Write the final line
                file.write(f") {max_dist} {max_dist} 0.0\n\n")
                        
                #Write second header
                file.write("! HADDOCK AIR restraints for 2nd partner\n")
                file.write("!\n")
                    
                # Generate the ambiguous restraints
                file.write(f"assign (resid {residue2} and segid B)\n")
                file.write("(\n")
                for res1 in pass_res1:
                    file.write(f"    (resid {res1} and segid A)\n")
                    # If it's not the last residue, write "or"
                    if res1 != pass_res1[-1]:
                        file.write("    or\n")
                    
                # Write the final line
                file.write(f") {max_dist} {max_dist} 0.0\n\n")
        
        logger.info(f"Passive restraints written to file.")

    def _find_nearby_residues(self, amino_acid_number, distance_threshold, pdb_file):
        parser = PDBParser(QUIET=True) # QUIET=True avoids comments on errors in the pdb.
        
        structures = parser.get_structure('protein', pdb_file)
        structure = structures[0] # 'structures' may contain several proteins in this case only one.
        target_atom = structure["A"][int(amino_acid_number)]['CA']
        
        atoms  = Selection.unfold_entities(structure, "A")
        ns = NeighborSearch(atoms)
        
        residue_search = ns.search(target_atom.coord, distance_threshold, level="R")
        close_residues = [residue.get_id()[1] for residue in residue_search]
        #logger.info(f"Residues close to {amino_acid_number}: {close_residues}")
        return close_residues


    def _generate_param_file(self, prot1, prot2, species1, species2):
        if os.path.exists(f"{prot1}_{prot2}_haddock.txt"):
            logger.info(f"HADDOCK param file already exists for {prot1}_{prot2}. Skipping...")
            return
        with open(f"{prot1}_{prot2}_haddock.txt", 'w') as file:
            file.writelines("# directory in which the scoring will be done\n")
            file.writelines(f"run_dir = \"{prot1}_{prot2}_run/\"\n\n")
            
            file.writelines("# execution mode\n")
            file.writelines("mode = \"local\"\n")
            file.writelines("ncores = 40\n\n")
            
            file.writelines("# molecules to be docked\n")
            file.writelines("molecules = [\n")
            file.writelines(f"    \"{prot1}_{species1}_fixed_A.pdb\",\n")
            file.writelines(f"    \"{prot2}_{species2}_fixed_B.pdb\"\n")
            file.writelines("]\n\n")
            
            file.writelines("# Parameters for each stage are defined below, prefer full paths\n")
            file.writelines("# ====================================================================\n")
            file.writelines("[topoaa]\n\n")
            
            file.writelines("[rigidbody]\n")
            file.writelines("tolerance = 20\n")
            file.writelines(f"unambig_fname = \"{prot1}_{prot2}_active_restraints.tbl\"\n")
            file.writelines(f"ambig_fname = \"{prot1}_{prot2}_passive_restraints.tbl\"\n")
            file.writelines("sampling = 20\n")
            file.writelines("amb_scale = 10\n")
            file.writelines("unamb_scale = 100\n\n")
            
            file.writelines("[caprieval]\n\n")
            
            file.writelines("# ====================================================================\n")

        logger.info("HADDOCK restraints file generated")


    def _transfer_files(self):
        shutil.move(f"{self.prot1.name}_{self.prot2.name}_active_restraints.tbl", 
                    f"Haddock3/{self.prot1.name}_{self.prot2.name}_active_restraints.tbl")
        shutil.move(f"{self.prot1.name}_{self.prot2.name}_passive_restraints.tbl", 
                    f"Haddock3/{self.prot1.name}_{self.prot2.name}_passive_restraints.tbl")
        shutil.move(f"{self.prot1.name}_{self.prot2.name}_haddock.txt", 
                    f"Haddock3/{self.prot1.name}_{self.prot2.name}_haddock.txt")
        shutil.copy(f"{self.prot1.pdb_path}",
                    f"Haddock3/{self.prot1.name}_{self.prot1.species}_fixed_A.pdb")
        shutil.copy(f"{self.prot2.pdb_path}",
                    f"Haddock3/{self.prot2.name}_{self.prot2.species}_fixed_B.pdb")
        self._change_pdb_chainA_to_B(f"Haddock3/{self.prot2.name}_{self.prot2.species}_fixed_B.pdb")
        logger.info("HADDOCK restraint files transferred to Haddock3 directory")
        
    def _change_pdb_chainA_to_B(self, input_file):
        with open(input_file, 'r') as f:
            lines = f.readlines()

        with open(input_file, 'w') as f:
            for line in lines:
                if line.startswith('ATOM'):
                    if line[21] == 'A':
                        line = line[:21] + 'B' + line[22:]
                f.write(line)

    def analyze_crosslinks(self):
        parser = PDBParser(QUIET=True)
        structure1 = parser.get_structure("prot1", str(self.prot1.pdb_path))
        structure2 = parser.get_structure("prot2", str(self.prot2.pdb_path))

        for idx, prot1, _, aa1, prot2, _, aa2, _, _ in self.full_data:
            res1 = structure1[0][prot1][(' ', aa1, ' ')]
            res2 = structure2[0][prot2][(' ', aa2, ' ')]

            distance = res1['CA'] - res2['CA']
            return distance
        

    def run_haddock(self):
        os.chdir(f"Haddock3")
        os.makedirs(f"{self.prot1.name}_{self.prot2.name}_run/", exist_ok=True)
        logger.info("Running Haddock3...")
        command = f"haddock3 {self.prot1.name}_{self.prot2.name}_haddock.txt"
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        # Check if there were any errors
        if process.returncode != 0:
            logger.info("Error occurred:")
            logger.info(stderr.decode("utf-8"))
        else:
            logger.info("Haddock3 completed successfully.")
        
        os.chdir("..")

    def fetch_haddock_results(self, out_path: Path):
        os.chdir(f"Haddock3")
        command = f"tar zxvf {self.prot1.name}_{self.prot2.name}_run/analysis/2_caprieval_analysis/summary.tgz"
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        
        # Check if there were any errors
        if process.returncode != 0:
            logger.info("Error occurred:")
            logger.info(stderr.decode("utf-8"))
        else:
            logger.info("File Unzipped successfully.")

            

        os.chdir("..")
        os.makedirs("CG_OpenMM", exist_ok=True)
        shutil.copy("Haddock3/model_1.pdb", f"CG_OpenMM/{self.prot1.name}_{self.prot2.name}_docked.pdb")

class CrossLinkAnalyzer:
    def __init__(self, csv_path: Path, species: str = None, prot1: Optional[Protein] = None, prot2: Optional[Protein] = None):
        self.prot1 = prot1
        self.prot2 = prot2
        self.species = species
        self.csv_path = csv_path
        self.prot1_values: Set[str] = set()
        self.prot2_values: Set[str] = set()
        self.full_data: List[Dict[str, str]] = []

    def parse_csv(self, species: str):
        for file_name in os.listdir(self.csv_path):
            # Check if the file ends with .csv
            if file_name.endswith(".csv"):
                file_path = os.path.join(self.csv_path, file_name)
                
                # Open and process the CSV file
                with open(file_path, 'r') as file:
                    reader = csv.DictReader(file)
                    for idx, row in enumerate(reader, start=1):
                        protein1 = row['Protein1']
                        peptide1 = row['Peptide1']
                        aa1 = int(row['AA1'])
                        protein2 = row['Protein2']
                        peptide2 = row['Peptide2']
                        aa2 = int(row['AA2'])
                        copies = int(row['Copies'])
                        crosslinker = row['Crosslinker']

                        self.prot1_values.add(protein1)
                        self.prot2_values.add(protein2)

                        self.full_data.append([idx, protein1, peptide1, aa1, protein2, peptide2, aa2, copies, crosslinker, file_path])
                
        logger.info("The following proteins pairs were found in the CSV directory:")
        for protein1, protein2 in zip(self.prot1_values, self.prot2_values):
            logger.info(f"Protein 1: {protein1} -> Protein 2:{protein2}")

    def validate_aa_positions(self):
        prot1 = self.prot1
        prot2 = self.prot2
        for protein in [prot1, prot2]:
            logger.info(f"Validating protein: {protein.name}")
            for entry in protein.residues:
                _, aa, peptide, _, _, _, _ = entry
                
                protein_fasta = protein.fasta_path

                # Validate and correct AA position
                corrected_pos = self._correct_aa_position(protein_fasta, aa, peptide)
                if corrected_pos is not None:
                    corrected_aa = corrected_pos + aa - 1  # Adjust aa based on corrected position                    
                    entry[1] = corrected_aa  # Update the entry with the corrected AA position
                else:
                    logger.info(f"Peptide '{peptide}' not found in {protein_fasta}")
        
        self._validate_presence_in_pdb(prot1, prot2)
    
    def _validate_presence_in_pdb(self, prot1, prot2):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", prot1.pdb_path)

        res_ids_prot1 = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    residue_id = residue.id
                    if residue_id[0] == ' ':
                        res_ids_prot1.append(residue_id[1])
        
        res_ids_prot2 = []
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", prot2.pdb_path)
        for model in structure:
            for chain in model:
                for residue in chain:
                    residue_id = residue.id
                    if residue_id[0] == ' ':
                        res_ids_prot2.append(residue_id[1])

        invalid_indices = []
        invalid_entry  = []
        # First loop to identify invalid entries
        for i, (entry1, entry2) in enumerate(zip(prot1.residues, prot2.residues)):
            set1, set2 = True, True
            if entry1[1] not in res_ids_prot1:
                set1 = False
            if entry2[1] not in res_ids_prot2:
                set2 = False
            if not set1 or not set2:
                invalid_indices.append(i)
                invalid_entry.append([entry1[0], 
                                        prot1.name, entry1[2], entry1[1], 
                                        prot2.name, entry2[2], entry2[1], 
                                        entry1[3], entry2[6], f"{entry1[1]}: {set1}, {entry2[1]}: {set2}"])

        # Second loop to remove invalid entries
        for index in sorted(invalid_indices, reverse=True):
            del prot1.residues[index]
            del prot2.residues[index]

        if invalid_entry != []:
            with open(f"{prot1.name}_{prot2.name}_invalid_residues.csv", 'w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(["IDx", "Protein1", "Peptide1", "AA1", "Protein2", "Peptide2", "AA2", "Copies", "File Path", "Residue Present?"])
                writer.writerows(invalid_entry)
            logger.info(f"Invalid residues removed from consideration. Residues saved to {prot1.name}_{prot2.name}_invalid_residues.csv")
 
    def _correct_aa_position(self, fasta_path, aa_pos, peptide):
        with open(fasta_path, "r") as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequence = str(record.seq)
                # Find the correct position of the peptide in the sequence
                corrected_pos = sequence.find(peptide)
                if corrected_pos != -1:
                    return corrected_pos + 1  # Convert to 1-based index
        return None
    
    def check_crosslinker_distance(self, crosslinker_name):
        # Dictionary mapping common crosslinker names to their approximate sizes in angstroms (Å)
        crosslinker_sizes = {
                "BS3": 11.4,                "BS(PEG)5": 21.7,
                "BS(PEG)9": 35.8,           "BS2G": 7.7,
                "DSG": 7.7,                 "DSS": 11.4,
                "DSSO": 10.3,               "DSP": 12.0,
                "EGS": 16.1,                "SDA": 3.9,
                "SIA": 3.9,                 "SMPB": 7.7,
                "SMCC": 8.3,                "SPDP": 6.8,
                "Sulfo-EGS": 16.1,          "Sulfo-SDA": 3.9,
                "Sulfo-SIA": 3.9,           "Sulfo-SMPB": 7.7,
                "Sulfo-SMCC": 8.3,          "Sulfo-LC-SPDP": 15.7,
                "EGS-d0": 16.1,             "EGS-d36": 16.1,
                "DTSSP": 12.0,              "DMP": 9.2,
                "BM(PEG)2": 14.7,           "BM(PEG)3": 17.8,
                "Sulfo-SBED": 12.5,         "EDC": 0,  # Zero-length crosslinker
                "GMBS": 7.3,                "Sulfo-GMBS": 7.3,
                "SPDP": 6.8,                "LC-SPDP": 15.7,
                "Sulfo-LC-SPDP": 15.7,      "BMPS": 5.9,
                "EMCS": 9.4,                "Sulfo-EMCS": 9.4,
                "MBS": 9.9,                 "Sulfo-MBS": 9.9,
                "SIAB": 10.6,               "Sulfo-SIAB": 10.6,
                "SBAP": 16.8,               "EDC-DE": 0,                
                "PhoX": 4.8,                "DMTMM": 25.0,    
                "ADH": 21.0,                
                "Unknown": 11.4     # Default value for unknown crosslinkers
            }
        
        if crosslinker_name in crosslinker_sizes:
            return crosslinker_sizes[crosslinker_name]
        else:
            logger.info("Crosslinker size not found in the database. Please provide the size in Å: ")
            raise logger.info(f"Crosslinker '{crosslinker_name}' not found. ")
        
    def process_buried_residues(self, files_path: Path, remove_buried: bool): 
        if not os.path.exists(f"{files_path}/sasa_files"):
            os.makedirs(f"{files_path}/sasa_files")
        
        out_path = f"{files_path}/sasa_files"
        for protein in [self.prot1, self.prot2]:
            if not os.path.exists(f"{out_path}/{protein.name}_buried_residues.txt"):
                input_file = str(protein.pdb_path)
                output_file = f"{out_path}/{protein.name}_sasa.pdb"
                self._write_sasa_files(input_file, output_file)
                self._extract_residues_with_sasa_over_threshold(output_file, f"{out_path}/{protein.name}_nonburied_residues.txt", threshold=40)
                logger.info(f"Non Buried residues file saved for {protein.name}")
            
        buried_residues = []
        for residue1, residue2 in zip(self.prot1.residues, self.prot2.residues):
            idx, aa1, peptide1, copies, _, _, _ = residue1
            idx, aa2, peptide2, copies, _, _, _ = residue2
            is_1_buried, is_2_buried = "Not Buried", "Not Buried"
            
            with open(f"{out_path}/{self.prot1.name}_nonburied_residues.txt", 'r') as file:
                nonburied_residues1 = file.read().split()
            with open(f"{out_path}/{self.prot2.name}_nonburied_residues.txt", 'r') as file:
                nonburied_residues2 = file.read().split()
            if aa1 in nonburied_residues1 and aa2 in nonburied_residues2:
                #logger.info(f"Residues {aa1} and {aa2} are not buried. Continuing with the simulation.")
                continue
            if aa1 not in nonburied_residues1:
                is_1_buried = "is Buried"
            if aa2 not in nonburied_residues2:
                is_2_buried = "is Buried"

            else:
                buried_residues.append([idx, self.prot1.name, aa1, peptide1, self.prot2.name, aa2, peptide2, copies, f"Residue {aa1} {is_1_buried} ", f"Residue {aa2} {is_2_buried}"])
                if remove_buried:
                    self.prot1.residues.remove(residue1)
                    self.prot2.residues.remove(residue2)
                    logger.info(f"Residues {aa1} and {aa2} removed from simulations.")
                continue
        
        with open(f"{files_path}/{self.prot1.name}_{self.prot2.name}_buried_residues.csv", 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["IDx", "Protein1", "Peptide1", "AA1", "Protein2", "Peptide2", "AA2", "Copies", "Notes"])
            writer.writerows(buried_residues)
        logger.info(f"Buried residues file saved for {self.prot1.name} and {self.prot2.name}")

    def _write_sasa_files(self, input_file, output_file):
        warnings.filterwarnings("ignore", category=Warning, module="freesasa")
        if os.path.exists(output_file):
            logger.info(f"SASA file already exists for {input_file}. Skipping...")
            return
        
        #options  = {'hydrogen' : True}
        try:
            structure = freesasa.Structure(input_file)
            result = freesasa.calc(structure)

            result.write_pdb(output_file)
        except Exception as e:
            pass

    def _extract_residues_with_sasa_over_threshold(self, file_path, output_file, threshold=40):
        residues_list = []
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith("ATOM"):
                    residue_name = line[17:20].strip()
                    residue_number = int(line[22:26].strip())
                    sasa = float(line[60:66].strip())  # SASA value position
                    if sasa > threshold:
                        residues_list.append((residue_name, residue_number))
        with open(output_file, 'w') as file:
            for residue in sorted(residues_list):
                file.write(f"{residue[1]} ")

class FinalResultVisualizer:
    def __init__(self, base_directory):
        self.base_directory = base_directory
        self.results = {}

    def analyze_folders(self):
        for folder_name in os.listdir(self.base_directory):
            folder_path = os.path.join(self.base_directory, folder_name)
            if os.path.isdir(folder_path):
                self.analyze_folder(folder_path)

    def analyze_folder(self, folder_path):
        folder_name = os.path.basename(folder_path)
        
        invalid_residues_count = 0
        all_crosslinks_count = 0

        for file_name in os.listdir(folder_path):
            if file_name.endswith("_invalid_residues.csv"):
                invalid_residues_count = self.count_lines(os.path.join(folder_path, file_name))
            elif file_name.endswith("_all_crosslinks.csv"):
                all_crosslinks_count = self.count_lines(os.path.join(folder_path, file_name))

        self.results[folder_name] = {
            "invalid_residues": invalid_residues_count,
            "all_crosslinks": all_crosslinks_count
        }

    def count_lines(self, file_path):
        with open(file_path, 'r') as csvfile:
            reader = csv.reader(csvfile)
            next(reader)  # Skip the header
            return sum(1 for row in reader)

    def get_results(self):
        return self.results

    def print_results(self):
        for folder, counts in self.results.items():
            print(f"Folder: {folder}")
            print(f"  Invalid Residues: {counts['invalid_residues']}")
            print(f"  All Crosslinks: {counts['all_crosslinks']}")
            print()
        

def main():
    parser = argparse.ArgumentParser(description="Cross-link analysis and protein processing")

    parser.add_argument('-i', '--input_dir', required=True, 
                        help="Path to the input directory containing CSV files. If this is the second time running the script, avoid repeated summarization step by providing the folder with summarized files ending in <_summarized>.")
    parser.add_argument("-o", "--out_dir", type=str, help="Path to the intermediate results.", required=True)

    parser.add_argument('-f', '--file', required=False, default = None, help="Skips to the specified file")
    parser.add_argument('-x', '--crosslinker', required=False, default = None,  help="Used to overwrite the identified cross-linker, if any.")
    parser.add_argument("-m", "--database", default="SWISS", help="Choose from which database the PDB file will be taken from: <SWISS> for SwissModel, or <AF> for AlphaFold. Default is <SWISS>.", type=str)
    
    parser.add_argument("--remove_buried", default= True, help="Default True. If False, will consider buried residues as available for crosslinking.", type=bool)
    parser.add_argument("--reject_interspecies", default= "False", help="Default False. If True, will reject crosslinks between proteins of different species. Give a species in UniProt format (i.e. HUMAN) to only analyze crosslinks of that species.", type=str)
    parser.add_argument("--protein", default= None, help="Will limit the analysis protein pairs that contain the desired protein", type=str)

    parser.add_argument("--species", default=None, help="Species for protein analysis. Will overwrite species name present in csv file.", type=str)
    parser.add_argument("--fix_partial_prot", default= True, help="If False, PDBFixer will attempt to fix the full protein + short simulation to fix clashes. Suggest trying to download from AlphaFold first.", type=bool)
    parser.add_argument("--skip_simulation1", default= False, help="Skip pre-docking simulation. Will assume all crosslinks are inter-protein.", type=bool)
    parser.add_argument("--skip_haddock", default= False, help="Skip generating HADDOCK files. Assumes docked structure will be provided under Haddock3/<protein1>_<protein2>_docked.pdb, or post-docking simulation will be skipped.", type=bool)
    parser.add_argument("--skip_simulation2", default= False, help="Skip post-docking simulation. Will not calculate inter-protein crosslinks.", type=bool)
    parser.add_argument("--crosslinker_size", default=None , help="Size of the crosslinker in Å. Will overwrite crosslinker size in csv file, if present.", type=int)
    parser.add_argument("--continue_from", default= False, help="Will skip proteins until it reaches the desired pair. Format <Protein1_Protein2>", type=str)
    
    args = parser.parse_args()

    original_path = os.getcwd()
    if not args.input_dir.endswith("_summarized"):
        CSVSummarizer(args.input_dir, args.crosslinker, args.file)
        input_dir = f"{args.input_dir}_summarized"
    else:
        input_dir = args.input_dir
    
    analyzer = CrossLinkAnalyzer(input_dir)
    analyzer.parse_csv(args.species)

    protein_pairs = defaultdict(list)
    # Iterate over each entry in full_data.
    # This loop should handle cases where the same pair of proteins is flipped.
    for data in analyzer.full_data:
        idx, protein1, peptide1, aa1, protein2, peptide2, aa2, copies, crosslinker, file_path = data
        
        if args.crosslinker_size == None:
            max_dist = analyzer.check_crosslinker_distance(crosslinker) + 3 # Add 3 Å to the crosslinker size for flexibility
        else:
            max_dist = args.crosslinker_size + 3

        # Use sorted tuple to ensure the same key for the same pair of proteins
        protein_pair = tuple(sorted([protein1, protein2]))
        
        # Determine the order of amino acid positions based on the sorted protein order
        if protein_pair == (protein1, protein2):
            aa_pair = (aa1, peptide1, aa2, peptide2)
        else:
            aa_pair = (aa2, peptide2, aa1, peptide1)
        
        # Append the data
        protein_pairs[protein_pair].append((idx, *aa_pair, copies, crosslinker, max_dist, file_path)) # List of lists will be [idx, aa1, peptide1, aa2, peptide2, copies, max_dist]
    
    os.makedirs(Path(args.out_dir), exist_ok=True)
    os.chdir(Path(args.out_dir))   
    
    for pair, data_list in protein_pairs.items():
        protein1, protein2 = pair
        protein1 = protein1.split("_")[0]
        protein2 = protein2.split("_")[0]
        if args.continue_from != False:
            prot1, prot2 = args.continue_from.split("_")
            if protein1 != prot1 or protein2 != prot2:
                continue
            else:
                args.continue_from = False
        
        if args.protein != None:
            if args.protein != protein1 and args.protein != protein2:
                continue

        logger.info(f"Working on proteins {protein1} and {protein2}")
        prot1_residues = []
        prot2_residues = []
        species1, species2 = "HUMAN", "HUMAN"
        if args.species == None:
            try:
                species1 = protein1.split("_")[1]
                species2 = protein2.split("_")[1]
            except IndexError:
                logger.info("Species not found within provided CSV file.")
                logger.info("Trying to fetch through UniProt API. Otherwise, it will be defaulted to HUMAN.")
        else:
            species1, species2 = args.species, args.species

        os.makedirs(f"{protein1}_{species1}_{protein2}_{species2}", exist_ok=True)
        files_path = f"{protein1}_{species1}_{protein2}_{species2}/pdb_fa_files"
        os.makedirs(files_path, exist_ok=True)
        os.chdir(f"{protein1}_{species1}_{protein2}_{species2}")

        with open(f"{protein1}_{species1}_{protein2}_{species2}_all_crosslinks.csv", 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["IDx", "Protein1", "Peptide1", "AA1", "Protein2", "Peptide2", "AA2", "Copies", "Crosslinker", "Crosslinker arm length", "Original File"])
            for data in data_list:
                idx, aa1, peptide1, aa2, peptide2, copies, crosslinker, max_dist, file_path = data
                writer.writerow([idx, f"{protein1}_{species1}", peptide1, aa1, f"{protein2}_{species2}", peptide2, aa2, copies, crosslinker, max_dist - 3, file_path])


        for data in data_list:
            idx, aa1, peptide1, aa2, peptide2, copies, crosslinker, max_dist, file_path = data
            
            # Append data to prot1_residues
            prot1_residues.append([idx, aa1, peptide1, copies, crosslinker, max_dist, file_path]) 
            
            # Append data to prot2_residues
            prot2_residues.append([idx, aa2, peptide2, copies, crosslinker, max_dist, file_path])

        # Extract data for prot1_obj
        prot1_obj = Protein(name=protein1.split("_")[0], species=species1,
                            files_path = files_path,
                            fasta_path=Path(f"pdb_fa_files/{protein1}_{species1}.fa"),
                            pdb_path=Path(f"pdb_fa_files/{protein1}_{species1}_fixed.pdb"),
                            residues=prot1_residues)
        
        # Extract data for prot2_obj
        prot2_obj = Protein(name=protein2.split("_")[0], species=species2,
                            files_path = files_path,
                            fasta_path=Path(f"pdb_fa_files/{protein2}_{species2}.fa"),
                            pdb_path=Path(f"pdb_fa_files/{protein2}_{species2}_fixed.pdb"),
                            residues=prot2_residues)   
    
        # Download FASTA and PDB files for both proteins
        downloader = ProteinDownloader(prot1_obj, prot2_obj)
        downloader.download_fasta()
        downloader.download_pdb(args.database)

        if args.reject_interspecies == "True":
            if prot1_obj.species != prot2_obj.species:
                os.chdir("..")
                continue
        elif args.reject_interspecies != "False":
            if prot1_obj.species.lower() != args.reject_interspecies.lower() or prot2_obj.species.lower() != args.reject_interspecies.lower():
                os.chdir("..")
                continue

        try:
            if not prot1_obj.fasta_path.is_file() or not prot1_obj.pdb_path.is_file():
                raise DownloadError(f"Missing FASTA or PDB file for {protein1}")
            if not prot2_obj.fasta_path.is_file() or not prot2_obj.pdb_path.is_file():
                raise DownloadError(f"Missing FASTA or PDB file for {protein2}")
        except DownloadError as de:
            logger.error(de)
            os.chdir("..")
            continue  

        downloader.process_pdb_with_residue_validation(args.fix_partial_prot)

        analyzer = CrossLinkAnalyzer(input_dir, prot1_obj.species, prot1_obj, prot2_obj)
        analyzer.validate_aa_positions()
        analyzer.process_buried_residues("pdb_fa_files", args.remove_buried)

        if prot1_obj.residues == [] or prot2_obj.residues == []:
            logger.info(f"No residues left for {protein1} and/or {protein2}. Skipping...")
            os.chdir("..")
            continue

        OpenMM_Sim = OpenMMSimulation(prot1_obj, prot2_obj)
        
        if protein1 == protein2 and not args.skip_simulation1:
            logger.info("Running simulation for single protein")
            OpenMM_Sim.run_simulation(sim_type="Single")


        Haddock_Sim = HaddockSimulation(prot1_obj, prot2_obj)
        if not args.skip_haddock:
            logger.info("Running HADDOCK docking")
            os.makedirs("Haddock3", exist_ok=True)
                         
            try:
                Haddock_Sim.generate_haddock_files()
            except statistics.StatisticsError as e:
                logger.info(f"Error {e} indicates no residues are left to guide docking. Skipping to next protein pair.")
                os.chdir("..")
                continue

            Haddock_Sim.run_haddock()
            Haddock_Sim.fetch_haddock_results(args.out_dir)    # Path with Haddock3 folder, OpenMM folder, and fasta/pdb files folder

        if not args.skip_simulation2:
            logger.info("Running simulation for docked protein")
            OpenMM_Sim.run_simulation(sim_type="Dimer")

        os.chdir("..")

    #Fetching of final results goes here
    FinalResultVisualizer(args.out_dir)

    os.chdir(original_path)
    logger.info("Done :D")

if __name__ == "__main__":
    main()
