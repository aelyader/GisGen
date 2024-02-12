import hashlib
import sys
from Bio import SeqIO
from PySide6.QtGui import QFont, QPixmap
from PySide6.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout,
                               QPushButton, QLabel, QLineEdit, QCheckBox, QTextEdit, QFileDialog, QMessageBox, QStyle,
                               QProgressBar)
from PySide6.QtCore import QThread, Signal, QTime, QTimer, Qt
import requests
from io import StringIO
import cProfile
from concurrent.futures import as_completed
import tarfile
import pandas as pd
import os
from concurrent.futures import ThreadPoolExecutor
from rapidfuzz import fuzz, process
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from memory_profiler import profile
import gc

warnings.filterwarnings("ignore", category=FutureWarning, message="When grouping with a length-1 list-like")
warnings.filterwarnings("ignore", category=FutureWarning, message="SeriesGroupBy.grouper is deprecated")


class WorkerThread(QThread):
    update_message = Signal(str)
    pause_execution = Signal()
    completed = Signal()  # Signal to indicate completion
    error_occurred = Signal(str)  # Signal to indicate an error with a message
    # Signal to indicate resuming of work
    resume_signal = Signal()
    progress_update = Signal(int)

    def __init__(self):
        super().__init__()
        self.gisaid_data_file = ""
        self.genbank_data_file = ""
        self.filter_by_omicron = False
        self.output_folder = ""
        self.should_install_fasta = False
        self.gisaid_fasta_folder_path = ""
        self.GIS_GEN_final = None
        self.sample_ids_file_path = ""
        self.sample_ids = []
        self.perform_downstream_analysis = False

    @staticmethod
    def split_gisaid_fasta(gisaid_fasta_file, output_folder, sample_id_list):
        # Check if the output folder exists; if not, create it
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Iterate through each record in the GISAID FASTA file
        for record in SeqIO.parse(gisaid_fasta_file, "fasta"):
            # Extract the sample ID from the record ID
            sample_id = record.id.split('|')[1]

            # Check if the sample ID is in the provided list
            if sample_id in sample_id_list:
                # Construct the output file path for the individual sample
                output_file_path = os.path.join(output_folder, f"{sample_id}.fasta")

                # Open the output file in write mode
                with open(output_file_path, "w") as f:
                    # Write the individual sample record to the output file
                    SeqIO.write(record, f, "fasta")

    @staticmethod
    def standardize_text(text):
        return text.strip().upper() if isinstance(text, str) else text

    def compare_columns(self, df, col1, col2):
        mismatch_count = 0
        mismatch_details = Counter()
        total_count = 0

        for value1, value2 in zip(df[col1], df[col2]):
            value1 = self.standardize_text(value1) if pd.notna(value1) else ''
            value2 = self.standardize_text(value2) if pd.notna(value2) else ''

            if value1 != value2:
                mismatch_count += 1
                mismatch_details[(value1, value2)] += 1
            total_count += 1

        mismatch_percentage = (mismatch_count / total_count) * 100 if total_count > 0 else 0
        # Get the top 3 mismatches
        top_mismatches = mismatch_details.most_common(3)
        return mismatch_percentage, top_mismatches

    def custom_print(self, msg):
        formatted_msg = msg + "\n"
        self.update_message.emit(formatted_msg)

    @profile
    def run(self):
        # Step 1: Reading the data
        self.custom_print("Step 1: Reading the data")
        tar_xz_directory = os.path.dirname(self.gisaid_data_file)
        # Check if the file is a .tar.xz file
        GISAID = pd.DataFrame()
        if self.gisaid_data_file.endswith('.tar.xz'):
            self.custom_print("The .tar.xz file exceeds 700MB, indicating that the TSV table size surpasses 12GB. "
                              "( The extraction and TSV upload process will require approximately 3 to 4 minutes. )")
            # Extract the .tar.xz file
            self.custom_print("-> Running extraction on the .tar.xz file")
            with tarfile.open(self.gisaid_data_file, 'r:xz') as tar:
                # Extract all the contents into the output folder
                self.progress_update.emit(5)
                tar.extractall(tar_xz_directory)
                self.progress_update.emit(10)
                # Find the .tsv file in the extracted files
                tsv_file = next((f for f in tar.getnames() if f.endswith('.tsv')), None)
                if tsv_file:
                    extracted_file_path = os.path.join(tar_xz_directory, tsv_file)
                    self.custom_print("-> Reading extracted GISAID TSV table")
                    GISAID = pd.read_csv(extracted_file_path, sep='\t', low_memory=False)
                    self.progress_update.emit(20)
                else:
                    self.custom_print("-> No .tsv file found in the extracted contents.")
                    return  # Exit the function if no .tsv file is found
        else:
            # Read the file directly if it's not a .tar.xz file
            self.custom_print("Please ensure you select the appropriate file for the (GISAID Meta Table) "
                              "The currently selected file is not in the .tar.xz format.")

        GenBank = pd.read_csv(self.genbank_data_file, engine='python')
        self.progress_update.emit(25)

        # Step 2: Pattern matching and filtering
        if self.filter_by_omicron:
            self.custom_print("Step 2: Creating Key and filtering")
            patterns = ["^B\\.1\\.1\\.529", "^BA\\.1", "^BA\\.2", "^BA\\.4", "^BA\\.5", "^XBB", "^XBB\\.1",
                        "^XBB\\.1\\.5", "^XBB\\.1\\.16"]
            regex_pattern = "|".join(patterns)
            GISAID = GISAID[GISAID['Pango lineage'].str.match(regex_pattern, na=False)]
            gc.collect()
        else:
            self.custom_print("Skipping Step 2: No filtering will be applied.")
        self.progress_update.emit(30)
        # Step 3: Renaming columns
        self.custom_print("Step 3: Renaming columns and Organizing Data")
        GISAID = GISAID.rename(columns={
            'Collection date': 'Collection_Date_GisGen',
            'Sequence length': 'Length_GisGen',
            'Pango lineage': 'Pangolin',
            'Pango version': 'PangoVersions',
            'Type': 'Genus',
            'Submission date': 'Submission_data',
            'Accession ID': 'AccessionID'
        })

        GenBank = GenBank.rename(columns={
            'Isolate': 'Isolate_Match_GisGen',
            'Release_Date': 'Submission_data',
            'Accession': 'AccessionID',
            'Length': 'Length_GisGen',
            'Collection_Date': 'Collection_Date_GisGen',
            'BioProject': 'Accession_BioPro',
            'BioSample': 'Accession_BioSam',
            'SRA_Accession': 'Accession_SRA'
        })
        self.progress_update.emit(35)
        # Step 4: Splitting columns
        self.custom_print("Step 4: Splitting columns and Normalizing")
        GISAID['Isolate_Match_GisGen'] = GISAID['Virus name'].str.split('/').str[2]
        self.custom_print("-> Finished Splitting")
        self.progress_update.emit(40)
        # Assign the split elements to new columns, with checks to avoid IndexError
        self.custom_print("-> Assigning the split elements to new columns, with checks to avoid IndexError")
        split_locations = GISAID['Location'].str.split('/')
        GISAID['Continent'] = split_locations.str.get(0)
        GISAID['Country'] = split_locations.str.get(1)
        GISAID['County_State'] = split_locations.str.get(2)
        self.progress_update.emit(45)
        # For 'County', check if the split list has at least 4 elements
        GISAID['County'] = split_locations.apply(lambda x: x[3] if len(x) > 3 else None)

        GISAID['Isolate'] = GISAID['Isolate_Match_GisGen']
        GenBank['Isolate'] = GenBank['Isolate_Match_GisGen']
        self.progress_update.emit(50)
        # Extracting country, state, and city information from the 'Geo_Location' column
        self.custom_print("-> Extracting country, state, and city from the 'Geo_Location' column")
        geo_location_pattern = r'^(?P<Country>[^:,]+)(?::\s*(?P<State>[^,]+))?(?:,\s*(?P<City>.*))?$'
        GenBank[['Country', 'County_State', 'County']] = GenBank['Geo_Location'].str.extract(geo_location_pattern)
        self.progress_update.emit(55)
        # Remove the 'Country' and 'USA' columns from the GenBank DataFrame
        if 'Location' in GISAID.columns:
            GISAID.drop('Location', axis=1, inplace=True)
        if 'USA' in GenBank.columns:
            GenBank.drop('USA', axis=1, inplace=True)
        if 'Geo_Location' in GenBank.columns:
            GenBank.drop('Geo_Location', axis=1, inplace=True)
        if 'Org_location' in GenBank.columns:
            GenBank.drop('Org_location', axis=1, inplace=True)
        gc.collect()
        gisaid_excluded_columns = ['Length_GisGen', 'Collection_Date_GisGen', 'Isolate_Match_GisGen']
        GISAID.rename(
            columns={col: col + '_GIS' if col not in gisaid_excluded_columns else col for col in GISAID.columns},
            inplace=True)
        GenBank.rename(
            columns={col: col + '_GEN' if col not in gisaid_excluded_columns else col for col in GenBank.columns},
            inplace=True)
        self.progress_update.emit(60)
        # Step 5: Preparing for fuzzy matching
        self.custom_print("Step 5: Preparing for fuzzy matching")
        common_columns = ['Length_GisGen', 'Collection_Date_GisGen',
                          'Isolate_Match_GisGen']  # Add other common columns here
        GIS_GEN_initial_merge = GISAID.merge(GenBank, on=common_columns)
        del GISAID, GenBank
        gc.collect()
        # Extract isolates for fuzzy matching
        gisaid_isolates = GIS_GEN_initial_merge['Isolate_GIS']  # Assuming 'Isolate_x' is from GISAID
        genbank_isolates = GIS_GEN_initial_merge['Isolate_GEN']  # Assuming 'Isolate_y' is from GenBank
        self.progress_update.emit(70)
        # Step 6: Starting parallelized fuzzy matching
        self.custom_print("Step 6: Starting parallelized fuzzy matching")
        with ThreadPoolExecutor() as executor:
            fuzzy_matches = list(executor.map(
                lambda x: self.fuzzy_match(x[0], x[1]),
                zip(gisaid_isolates, genbank_isolates)
            ))

        # Filter the merged DataFrame based on fuzzy match results
        self.custom_print("-> Filter the merged DataFrame based on fuzzy match results")
        GIS_GEN_initial_merge['Fuzzy_Match'] = fuzzy_matches
        self.GIS_GEN_final = GIS_GEN_initial_merge.dropna(subset=['Fuzzy_Match'])

        # Sort the columns alphabetically, excluding the ones that should remain at the beginning
        sorted_columns = sorted([col for col in self.GIS_GEN_final.columns if col not in common_columns])
        # Place the common columns at the beginning
        final_column_order = common_columns + sorted_columns
        self.GIS_GEN_final = self.GIS_GEN_final[final_column_order]

        if 'Fuzzy_Match' in self.GIS_GEN_final.columns:
            self.GIS_GEN_final.drop('Fuzzy_Match', axis=1, inplace=True)
        self.progress_update.emit(80)
        # List of column pairs to compare
        column_pairs = [
            ('Country_GIS', 'Country_GEN'),
            ('County_GIS', 'County_GEN'),
            ('County_State_GIS', 'County_State_GEN'),
            ('Genus_GIS', 'Genus_GEN'),
            ('Isolate_GIS', 'Isolate_GEN'),
            ('PangoVersions_GIS', 'PangoVersions_GEN'),
            ('Pangolin_GIS', 'Pangolin_GEN'),
            ('Submission_data_GIS', 'Submission_data_GEN')
        ]

        self.custom_print("Step 7: Generating Value Mismatch report between GISAID and GenBank")
        mismatch_results = {}
        for col1, col2 in column_pairs:
            mismatch_percentage, top_mismatches = self.compare_columns(self.GIS_GEN_final, col1, col2)
            mismatch_results[f'{col1} vs {col2}'] = {
                'Mismatch_Percentage': mismatch_percentage,
                'Top_Mismatches': top_mismatches
            }

        if self.output_folder:
            mismatch_report_path = os.path.join(self.output_folder, 'mismatch_report.txt')
            with open(mismatch_report_path, 'w') as report_file:
                for pair, result in mismatch_results.items():
                    report_file.write(f'Mismatch Percentage for {pair}: {result["Mismatch_Percentage"]:.2f}%\n')
                    for i, (mismatch_pair, count) in enumerate(result['Top_Mismatches'], start=1):
                        report_file.write(f' - Top Mismatch {i}: {mismatch_pair}, Count: {count}\n')
            self.custom_print(f"Mismatch report saved to {mismatch_report_path}")
        self.progress_update.emit(90)
        if self.output_folder:
            self.custom_print(
                f"Step 8: Saving merged tables as CSV files in the selected directory: {self.output_folder}")
            self.GIS_GEN_final.to_csv(os.path.join(self.output_folder, 'GIS_GEN_final.csv'), index=False)
            self.custom_print("-> Files successfully saved.")
        else:
            self.custom_print("-> No output folder selected. Please select an output folder.")

            # Extract GISAID IDs as a list
        self.sample_ids = self.GIS_GEN_final['AccessionID_GIS'].tolist()

        # Section to save the sample_ids list to a text file
        self.sample_ids_file_path = os.path.join(self.output_folder, "sample_ids.txt")
        with open(self.sample_ids_file_path, 'w') as file:
            for sample_id in self.sample_ids:
                file.write(sample_id + '\n')
        self.progress_update.emit(100)
        # Notify the user through the custom_print method
        self.custom_print(f"-> GISAID Sample IDs .txt file saved to  {self.sample_ids_file_path}")
        # Emit the pause_execution signal
        self.pause_execution.emit()

    @staticmethod
    def download_data(url):
        response = requests.get(url)
        if response.status_code == 200:
            data = StringIO(response.text)
            df = pd.read_csv(data)
            return df
        else:
            raise Exception(f"-> Failed to download data, status code: {response.status_code}")

    @profile
    def resume_work(self):
        # Emit signal instead of directly performing tasks
        self.resume_signal.emit()

    @profile
    def perform_heavy_tasks(self):
        self.custom_print("-> Let's Continue to the MD5 validation steps")

        # Determine the original GISAID FASTA file
        fasta_files = [f for f in os.listdir(self.gisaid_fasta_folder_path) if f.endswith('.fasta')]
        if fasta_files:
            original_gisaid_fasta_file = fasta_files[0]  # Assuming the first file is the original
        else:
            self.custom_print("-> Error: No FASTA files found in the GISAID FASTA folder.")
            return
        self.custom_print("Step 9: Splitting multi-sample GISAID FASTA")
        self.progress_update.emit(10)
        # Ensure the GISAID FASTA folder path is provided
        if not self.gisaid_fasta_folder_path:
            self.custom_print("-> Error: GISAID FASTA folder path is not provided.")
            return
        self.progress_update.emit(20)
        # Find the FASTA file in the selected folder
        fasta_files = [f for f in os.listdir(self.gisaid_fasta_folder_path) if f.endswith('.fasta')]

        if len(fasta_files) == 1:
            gisaid_fasta_file_path = os.path.join(self.gisaid_fasta_folder_path, fasta_files[0])
            self.split_gisaid_fasta(gisaid_fasta_file_path, self.gisaid_fasta_folder_path, self.sample_ids)
        else:
            self.custom_print(
                f"-> Error: The folder {self.gisaid_fasta_folder_path} should contain exactly one FASTA file.")
            return
        self.progress_update.emit(30)
        # original_gisaid_fasta_file = os.path.basename(self.gisaid_fasta_folder_path)  # Corrected attribute name
        self.custom_print("Step 10: Convert FASTA files to raw format")

        for fasta_file in os.listdir(self.gisaid_fasta_folder_path):
            self.reformat_fasta(os.path.join(self.gisaid_fasta_folder_path, fasta_file), original_gisaid_fasta_file)

        # Insert the MD5 validation code here

        self.progress_update.emit(45)
        fasta_files = [os.path.join(self.gisaid_fasta_folder_path, f)
                       for f in os.listdir(self.gisaid_fasta_folder_path)
                       if f.endswith('.fasta')]
        md5_results = {}
        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = {executor.submit(self.generate_md5_parallel, fasta_file): fasta_file for fasta_file in
                       fasta_files}
            for future in as_completed(futures):
                accession_id, gisaid_md5 = future.result()
                md5_results[accession_id] = gisaid_md5
        self.progress_update.emit(55)
        # Download GenBank MD5 table
        self.custom_print("Step 11: Download GenBank FASTA MD5 table")
        genbank_md5_table_url = "https://zenodo.org/record/10578584/files/updated_genbank_sequences.csv"
        genbank_md5_table = self.download_data(genbank_md5_table_url)
        gc.collect()
        self.progress_update.emit(75)
        # Compare MD5 checksums
        md5_comparison_results = pd.DataFrame(columns=['GISAID_ID', 'GenBank_ID', 'GISAID_MD5', 'GenBank_MD5', 'Match'])
        self.custom_print("Step 12: Compare MD5 checksums between GISAID & GenBank")

        for index, row in self.GIS_GEN_final.iterrows():
            gisaid_id = row['AccessionID_GIS']
            gisaid_md5 = md5_results.get(gisaid_id, None)

            # Find the corresponding GenBank MD5 from the loaded table
            genbank_md5_row = genbank_md5_table[genbank_md5_table['Accession'] == row['AccessionID_GEN']]
            genbank_md5 = genbank_md5_row['GenBank_MD5'].iloc[0] if not genbank_md5_row.empty else None
            # Compare MD5 checksums
            is_match = 'Yes' if gisaid_md5 == genbank_md5 else 'No'

            # Append result to the DataFrame
            new_row = {'GISAID_ID': gisaid_id, 'GenBank_ID': row['AccessionID_GEN'],
                       'GISAID_MD5': gisaid_md5, 'GenBank_MD5': genbank_md5, 'Match': is_match}
            new_row_df = pd.DataFrame([new_row])
            md5_comparison_results = pd.concat([md5_comparison_results, new_row_df], ignore_index=True)
            # Logging

        self.progress_update.emit(85)
        # Save the md5_results DataFrame to CSV
        csv_file_path = os.path.join(self.output_folder, 'md5_results.csv')
        md5_comparison_results.to_csv(csv_file_path, index=False)

        self.custom_print("MD5 Validation Done!")

        if self.perform_downstream_analysis:
            self.custom_print("Step 13: Resuming work. Performing downstream analysis...")
            self.perform_downstream_analysis_function()

        self.progress_update.emit(100)

    def perform_downstream_analysis_function(self):
        # self.custom_print("Performing downstream analysis...")
        # Load your DataFrame
        df = self.GIS_GEN_final
        # Perform analyses
        self.temporal_analysis_of_variants(df)
        self.geographical_distribution_of_variants(df)
        self.sequence_length_variation_analysis(df)
        self.amino_acid_substitutions_analysis(df)
        self.submission_date_lag_analysis(df)
        self.custom_print("Downstream analysis completed.")

    def temporal_analysis_of_variants(self, df):
        variant_counts = df.groupby(['Collection_Date_GisGen', 'Variant_GIS']).size().unstack(fill_value=0)
        # Create a larger figure to give more space for labels
        plt.figure(figsize=(15, 8))  # Adjust the figure size as needed
        sns.lineplot(data=variant_counts)
        # Rotate x-axis labels to avoid overlap
        plt.xticks(rotation=90, ha='center')  # Rotate the labels vertically
        # Optional: Set a custom font size for x-axis labels if needed
        plt.tick_params(axis='x', which='major', labelsize=8)
        # Adjust the frequency of labels shown on the x-axis
        plt.gca().xaxis.set_major_locator(plt.MaxNLocator(10))  # Show only 10 labels
        plt.ylabel('Count')
        plt.title('Temporal Analysis of Virus Variants')
        plt.tight_layout()  # This will make sure the labels fit into the figure area
        plt.savefig(os.path.join(self.output_folder, 'temporal_analysis_variants.png'))
        variant_counts.to_csv(os.path.join(self.output_folder, 'temporal_analysis_variants.csv'))

    def geographical_distribution_of_variants(self, df):
        # Assuming 'Country' and 'State' columns are available
        variant_distribution = df.groupby(['Country_GIS', 'Variant_GIS']).size().unstack(fill_value=0)
        # Saving as CSV
        variant_distribution.to_csv(os.path.join(self.output_folder, 'geographical_distribution_variants.csv'))

    def sequence_length_variation_analysis(self, df):
        # Determine the most common Pangolin lineages
        top_lineages = df['Pangolin_GEN'].value_counts().nlargest(10).index.tolist()

        # Create a new column for the top lineages and 'Other'
        df['Pangolin_GEN_Top'] = df['Pangolin_GEN'].apply(lambda x: x if x in top_lineages else 'Other')

        # Define a color palette for the top lineages, with a default color for 'Other'
        palette = sns.color_palette('husl', n_colors=len(top_lineages) + 1)
        lineage_colors = {lineage: color for lineage, color in zip(top_lineages + ['Other'], palette)}

        # Box Plot for Sequence Length Variation
        plt.figure(figsize=(14, 8))
        sns.boxplot(x='Pangolin_GEN_Top', y='Length_GisGen', hue='Pangolin_GEN_Top', data=df,
                    palette=lineage_colors, dodge=False)
        plt.title('Sequence Length Variation for Top Pangolin Lineages')
        plt.xlabel('Pangolin Lineage')
        plt.ylabel('Sequence Length')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_folder, 'sequence_length_variation_boxplot.png'))

        # Histogram for Sequence Length Variation
        plt.figure(figsize=(14, 8))
        sns.histplot(data=df, x='Length_GisGen', hue='Pangolin_GEN_Top', element="step", bins=30, kde=True,
                     palette=lineage_colors)
        plt.title('Sequence Length Variation for Top Pangolin Lineages (and Others)')
        plt.xlabel('Sequence Length')
        plt.ylabel('Count')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_folder, 'sequence_length_variation_histogram.png'))

        # Statistical analysis (Optional, if you still want to save this as well)
        length_stats = df.groupby('Pangolin_GEN_Top')['Length_GisGen'].describe()
        length_stats.to_csv(os.path.join(self.output_folder, 'sequence_length_variation_stats.csv'))

    def submission_date_lag_analysis(self, df):
        # Convert dates to datetime objects with inferred format
        df['Collection_Date_GisGen'] = pd.to_datetime(df['Collection_Date_GisGen'], errors='coerce', format='%Y-%m-%d')
        df['Submission_data_GEN'] = pd.to_datetime(df['Submission_data_GEN'], errors='coerce')

        # Calculate lag
        df['Submission_Lag'] = (df['Submission_data_GEN'] - df['Collection_Date_GisGen']).dt.days

        # Filter out rows where datetime conversion failed (resulting in NaT)
        df_filtered = df.dropna(subset=['Collection_Date_GisGen', 'Submission_data_GEN'])

        # Lag statistics
        lag_stats = df_filtered['Submission_Lag'].describe()
        lag_stats.to_csv(os.path.join(self.output_folder, 'submission_date_lag_stats.csv'))

        # Plotting
        plt.figure(figsize=(12, 6))
        sns.histplot(df_filtered, x='Submission_Lag', bins=30, kde=True)
        plt.title('Submission Date Lag Analysis')
        plt.savefig(os.path.join(self.output_folder, 'submission_date_lag.png'))

    def amino_acid_substitutions_analysis(self, df):
        # Parse and count substitutions
        substitutions = Counter(", ".join(df['AA Substitutions_GIS'].dropna()).split(", "))
        common_substitutions = pd.Series(substitutions).sort_values(ascending=False)

        # Saving the output
        common_substitutions.to_csv(os.path.join(self.output_folder, 'amino_acid_substitutions.csv'))

    def generate_md5_parallel(self, fasta_file):
        gisaid_md5 = self.generate_md5(fasta_file)
        return os.path.basename(fasta_file).replace('.fasta', ''), gisaid_md5

    def generate_md5(self, filepath):
        """
        Generate MD5 checksum for a given file.

        Args:
            filepath (str): Path to the file for which to generate the MD5 checksum.

        Returns:
            str: The generated MD5 checksum.
        """
        try:
            hash_md5 = hashlib.md5()
            with open(filepath, "rb") as f:
                for chunk in iter(lambda: f.read(8192), b""):
                    hash_md5.update(chunk)
            return hash_md5.hexdigest()
        except FileNotFoundError:
            self.custom_print(f"File not found: {filepath}")
            return None

    def reformat_fasta(self, filepath, original_gisaid_fasta_file):
        """
        Reformat FASTA files to a specific raw format.

        Args:
            filepath (str): Path to the FASTA file to be reformatted.
            original_gisaid_fasta_file (str): Name of the original GISAID FASTA file to skip reformatting.
        """
        # Skip processing the original GISAID multi-sample FASTA file
        if os.path.basename(filepath) == original_gisaid_fasta_file:
            self.custom_print(f"Skipping reformatting of the original GISAID file: {original_gisaid_fasta_file}")
            return

        try:
            # Open the FASTA file in read mode to read the first line
            with open(filepath, 'r') as f:
                first_line = f.readline().strip()

            # Check if the first line is already in the required format; if so, log it and return
            if first_line == ">sample":
                self.custom_print(f"File {os.path.basename(filepath)} is already reformatted. Skipping.")
                return

            # Open the FASTA file in read mode to read all lines
            with open(filepath, 'r') as f:
                lines = f.readlines()

            # Open the FASTA file in write mode to overwrite it
            with open(filepath, 'w') as f:
                # Loop through each line to write it back
                for line in lines:
                    # If the line starts with '>', write '>sample' instead
                    if line.startswith(">"):
                        f.write(">sample\n")
                    else:
                        # Otherwise, write the line as is
                        f.write(line.strip())
        except Exception as e:
            self.custom_print(f"An error occurred while reformatting {filepath}: {e}")

    def fuzzy_match(self, x, choice):
        try:
            match, score, _ = process.extractOne(x, [choice], scorer=fuzz.token_set_ratio)
            # self.custom_print(f"Processing record: {x}, Best match: {match}, Score: {score}")
            if score >= 90:  # Threshold
                return match
            else:
                return None
        except FileNotFoundError as e:
            self.custom_print(f"File not found error: {e}")
        except Exception as e:
            self.custom_print(f"An error occurred: {e}")
            return None
        pass


class ProgressWindow(QWidget):
    def __init__(self):
        super().__init__()
        # Declare instance attributes here
        self.startTime = None
        self.timer = QTimer(self)
        self.layout = QVBoxLayout(self)
        self.elapsedTimeLabel = None  # Declare now, initialize in initUI
        self.main_progress_bar = None  # Declare now, initialize in initUI
        self.text_box = None  # Declare now, initialize in initUI
        # Connect the timer's timeout signal to updateElapsedTime
        self.timer.timeout.connect(self.updateElapsedTime)
        self.initUI()

    def initUI(self):
        self.setWindowTitle("Progress")

        # Elapsed time label
        self.elapsedTimeLabel = QLabel("00:00", self)
        self.layout.addWidget(self.elapsedTimeLabel)

        self.main_progress_bar = QProgressBar(self)
        self.layout.addWidget(self.main_progress_bar)

        self.text_box = QTextEdit(self)
        self.text_box.setReadOnly(True)
        self.layout.addWidget(self.text_box)

        self.setGeometry(500, 500, 600, 600)

    def startElapsedTimeTimer(self):
        self.startTime = QTime.currentTime()  # Record the start time
        self.timer.start(1000)  # Start the timer to update every second

    def updateElapsedTime(self):
        if self.startTime:
            # Calculate the elapsed time
            now = QTime.currentTime()
            elapsed_time = self.startTime.secsTo(now)

            # Convert seconds to minutes and seconds
            minutes, seconds = divmod(elapsed_time, 60)
            # Format the string to display
            time_string = "{:02d}:{:02d}".format(minutes, seconds)

            self.elapsedTimeLabel.setText(time_string)

    def update_progress(self, value):
        self.main_progress_bar.setValue(value)

    def update_message(self, message):
        # Check if the message starts with any of the specified keywords
        if message.startswith("Step"):
            # Apply HTML formatting for bold and larger font, starting with a line break
            formatted_message = f"<br><div style='font-weight:bold; font-size:14pt; margin-left: 0px;'>{message}</div>"
        elif message.startswith(("->", "Skipping", "Mismatch")):
            # For lines starting with "->", add indentation using a div for proper margin application
            formatted_message = f"<br><div style='font-size:12pt; margin-left: 35px;'>{message}</div>"
        elif message.startswith(("MD5", "Downstream")):
            formatted_message = f"<br><div style='margin-left: 0px; font-weight:bold; color: green; font-size:18pt; margin-left: 50px;'>{message}</div><br>"
        else:
            # Apply a smaller font size increase for other messages, starting with a line break
            formatted_message = f"<br><div style='font-size:12pt; margin-left: 0px;'>{message}</div>"

        # Append the formatted message to the text box
        self.text_box.append(formatted_message)


class MainApplication(QWidget):
    output_folder = ""

    def __init__(self):
        super().__init__()
        self.worker_thread = WorkerThread()
        self.main_layout = QVBoxLayout()
        self.perform_downstream_analysis = False
        self.worker_thread.resume_signal.connect(self.worker_thread.perform_heavy_tasks)
        self.gisaid_entry = QLineEdit()
        self.genbank_entry = QLineEdit()  # Assuming similar pattern for other QLineEdit attributes
        self.gisaid_fasta_folder_entry = QLineEdit()
        self.output_folder_entry = QLineEdit()
        self.filter_omicron_checkbox = QCheckBox()
        self.downstream_analysis_checkbox = QCheckBox()
        self.progress_window = None
        self.initUI()

    def initUI(self):
        self.setWindowTitle("GISAID <--> GenBanks Sample Matching")
        self.setLayout(self.main_layout)
        self.addComponents()
        self.worker_thread.pause_execution.connect(self.wait_for_user_input)

    def addComponents(self):
        self.main_layout.addLayout(self.createSubtitleLayout())
        self.main_layout.addLayout(self.createGISAIDLayout())
        self.main_layout.addLayout(self.createGenBankLayout())
        self.main_layout.addLayout(self.createOutputFolderLayout())
        self.main_layout.addLayout(self.createGISAIDFASTAFolderLayout())
        self.main_layout.addLayout(self.createFilterOmicronLayout())
        self.main_layout.addLayout(self.createDownstreamAnalysisCheckbox())
        self.main_layout.addLayout(self.createButtonsLayout())

    def createSubtitleLayout(self):
        # Create a vertical layout
        vertical_layout = QVBoxLayout()

        # Create a label for the logo and set the pixmap
        logo_label = QLabel(self)
        # Get the directory in which the script is located
        script_dir = os.path.dirname(os.path.abspath(__file__))
        image_path = os.path.join(script_dir, "images", "Final_gisgen_logo_new.jpg")
        logo_pixmap = QPixmap(image_path)  # Adjust path if necessary
        if logo_pixmap.isNull():  # Check if the QPixmap failed to load
            print("The logo image was not found. Displaying placeholder text instead.")
            logo_label.setText("Logo Not Found")  # Placeholder text
        else:
            logo_label.setPixmap(logo_pixmap.scaled(400, 400, Qt.KeepAspectRatio))

        vertical_layout.addWidget(logo_label, alignment=Qt.AlignCenter)

        # Create and configure the Info Help button
        info_button = QPushButton("Info Help \n(Step by Step Guide)\n")
        font = QFont()
        font.setPointSize(11)
        info_button.setFont(font)
        info_button.clicked.connect(self.show_instructions)
        vertical_layout.addWidget(info_button)  # Add the button below the logo

        return vertical_layout

    def createGISAIDLayout(self):
        gisaid_layout = QHBoxLayout()
        gisaid_label_text = """<b><span style="font-size:12pt;">GISAID Meta Table:</span></b><br>
            [ .tar.xz sample metatable from <br><a href="https://gisaid.org/">https://gisaid.org/</a> ]"""
        gisaid_label = QLabel(gisaid_label_text)

        gisaid_button = QPushButton("Browse")
        gisaid_button.clicked.connect(lambda: self.select_file(self.gisaid_entry))
        gisaid_layout.addWidget(gisaid_label)
        gisaid_layout.addWidget(self.gisaid_entry)
        gisaid_layout.addWidget(gisaid_button)
        return gisaid_layout

    def createGenBankLayout(self):
        genbank_layout = QHBoxLayout()
        genbank_label = QLabel(
            """<br><b><span style="font-size:12pt;">GenBank Meta Table:</span></b><br> [ CSV sample metatable from <br> <a href="https://www.ncbi.nlm.nih.gov/labs/virus">https://www.ncbi.nlm.nih.gov/labs/virus</a> ]""")
        genbank_button = QPushButton("Browse")
        genbank_button.clicked.connect(lambda: self.select_file(self.genbank_entry))
        genbank_layout.addWidget(genbank_label)
        genbank_layout.addWidget(self.genbank_entry)
        genbank_layout.addWidget(genbank_button)
        return genbank_layout

    def createOutputFolderLayout(self):
        output_folder_layout = QHBoxLayout()
        output_folder_label = QLabel(
            """<br><b><span style="font-size:12pt;">Output Folder:</span></b><br> [ Select a folder where results will be saved ]""")
        output_folder_button = QPushButton("Browse")
        output_folder_button.clicked.connect(lambda: self.select_folder(self.output_folder_entry))
        output_folder_layout.addWidget(output_folder_label)
        output_folder_layout.addWidget(self.output_folder_entry)
        output_folder_layout.addWidget(output_folder_button)
        return output_folder_layout

    def createGISAIDFASTAFolderLayout(self):
        gisaid_fasta_folder_layout = QHBoxLayout()
        gisaid_fasta_folder_label = QLabel(
            """<br><b><span style="font-size:12pt;">GISAID FASTA Folder:</span></b><br> [ Select a folder for the GISAID FASTA.<br> <span style="color:red;">This folder needs to be empty.</span> ]<br>""")
        self.gisaid_fasta_folder_entry = QLineEdit()
        gisaid_fasta_folder_button = QPushButton("Browse")
        gisaid_fasta_folder_button.clicked.connect(lambda: self.select_folder(self.gisaid_fasta_folder_entry))
        gisaid_fasta_folder_layout.addWidget(gisaid_fasta_folder_label)
        gisaid_fasta_folder_layout.addWidget(self.gisaid_fasta_folder_entry)
        gisaid_fasta_folder_layout.addWidget(gisaid_fasta_folder_button)
        return gisaid_fasta_folder_layout

    def createFilterOmicronLayout(self):
        filter_omicron_layout = QHBoxLayout()
        filter_omicron_layout.addStretch(1)
        self.filter_omicron_checkbox = QCheckBox("[ Filter by Omicron ]")  # Define the checkbox as an attribute
        font = QFont()
        font.setPointSize(11)  # Adjust the size as needed
        self.filter_omicron_checkbox.setFont(font)

        filter_omicron_layout.addWidget(self.filter_omicron_checkbox)
        filter_omicron_layout.addStretch(1)
        return filter_omicron_layout

    def createDownstreamAnalysisCheckbox(self):
        downstream_analysis_layout = QHBoxLayout()
        downstream_analysis_layout.addStretch(1)

        self.downstream_analysis_checkbox = QCheckBox("[ Enable Downstream Analysis ]")
        font = QFont()
        font.setPointSize(11)  # Adjust the size as needed
        self.downstream_analysis_checkbox.setFont(font)

        downstream_analysis_layout.addWidget(self.downstream_analysis_checkbox)
        analysis_info_button = QPushButton()
        analysis_info_button.setIcon(self.style().standardIcon(getattr(QStyle, "SP_MessageBoxInformation")))
        analysis_info_button.clicked.connect(self.show_analysis_description)  # Connect to the new method

        downstream_analysis_layout.addWidget(analysis_info_button)
        downstream_analysis_layout.addStretch(1)
        return downstream_analysis_layout

    def createButtonsLayout(self):
        buttons_layout = QHBoxLayout()
        run_button = QPushButton("Run")
        font = QFont()
        font.setPointSize(14)  # Adjust the size as needed
        font.setBold(True)  # Make the text bold
        run_button.setFont(font)

        run_button.clicked.connect(self.run_script)
        buttons_layout.addWidget(run_button)
        return buttons_layout

    def run_script(self):
        # Warning message about computational power and GUI suspension
        warning_msg = QMessageBox()
        warning_msg.setIcon(QMessageBox.Icon.Warning)
        warning_msg.setWindowTitle("Warning")
        warning_msg.setText(
            "<p><span style='font-weight:bold; font-size:13pt;'>WARNING</span></p>"
            "<p><span font-size:13pt;'>This program needs substantial computational power and may temporarily suspend the GUI to allocate more threads for processing, depending on your hardware's capabilities.</span></p>"
            "<p><span style='font-weight:bold; font-size:13pt;'>**** If the GUI seems unresponsive, please avoid closing the application. ****</span></p>")
        warning_msg.setStandardButtons(
            QMessageBox.StandardButton.Ok | QMessageBox.StandardButton.Cancel)  # Uncomment for PyQt6

        response = warning_msg.exec()

        if response == QMessageBox.StandardButton.Ok:
            if not self.worker_thread.isRunning():
                # Set up the variables for the worker thread
                self.worker_thread.gisaid_data_file = self.gisaid_entry.text()
                self.worker_thread.genbank_data_file = self.genbank_entry.text()
                self.worker_thread.filter_by_omicron = self.filter_omicron_checkbox.isChecked()
                self.worker_thread.output_folder = self.output_folder_entry.text()
                self.worker_thread.gisaid_fasta_folder_path = self.gisaid_fasta_folder_entry.text()
                self.worker_thread.perform_downstream_analysis = self.downstream_analysis_checkbox.isChecked()
                MainApplication.output_folder = self.output_folder_entry.text()
                self.worker_thread.output_folder = MainApplication.output_folder

                # Start the worker thread
                self.worker_thread.start()

                # Initialize and display the progress window
                self.progress_window = ProgressWindow()  # Instantiate ProgressWindow when needed
                self.progress_window.show()
                self.worker_thread.progress_update.connect(self.progress_window.update_progress)
                self.worker_thread.update_message.connect(self.progress_window.update_message)
                self.progress_window.startElapsedTimeTimer()

        elif response == QMessageBox.StandardButton.Cancel:
            # If the user cancels, do not start the process
            pass

    def resume_script(self):
        if self.downstream_analysis_checkbox.isChecked():
            self.worker_thread.perform_downstream_analysis = True
        self.worker_thread.resume_work()

    def select_file(self, line_edit):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File")
        if file_path:  # This checks if a file was actually selected
            line_edit.setText(file_path)

    def select_folder(self, line_edit):
        folder_path = QFileDialog.getExistingDirectory(self, "Select Folder")
        if folder_path:  # This checks if a folder was actually selected
            line_edit.setText(folder_path)

    def wait_for_user_input(self):
        # Create a modal dialog box to wait for user input
        gisaid_fasta_path = self.worker_thread.gisaid_fasta_folder_path
        file_path = self.worker_thread.sample_ids_file_path
        msgBox = QMessageBox()
        msgBox.setWindowTitle("Pause Execution")
        msgBox.setText(f"\n         **********   STOP   **********\n"
                       f"\nThe Code has PAUSED. Please complete the required steps and then press 'Continue'."
                       f"\nPlease Follow the steps below:\n"
                       f"\nStep 1:"
                       f"\n ----> Open a web browser and Go to https://www.epicov.org/ "
                       f"\n     ----> login using your GISAID credentials\n"
                       f"\nStep 2:"
                       f"\n ----> In https://www.epicov.org/ select [EpiCoV] and under EpiCoV, select [Search]"
                       f"\n     ----> At the bottom right hand side between EPI_SET and Analysis box, click [Select]"
                       f"\n         ----> Click [Choose File] and select the following text file:"
                       f"\n             ---->  {file_path}\n"
                       f"\nStep 3:"
                       f"\n ----> Download multi-fasta file for the queried sample set\n"
                       f"\nStep 4:"
                       f"\n ----> Place the single multi-fasta file in the following folder:  "
                       f"\n     ----> {gisaid_fasta_path}"
                       f"\n         ----> Make sure the single multi-fasta file is the only file in the folder\n"
                       f"\n Now you can Click 'OK' to proceed.\n")

        msgBox.setStandardButtons(QMessageBox.StandardButton.Ok)
        if msgBox.exec() == QMessageBox.StandardButton.Ok:
            # Logic after the user presses 'Continue' on the dialog box
            # This will resume the worker thread
            self.worker_thread.resume_work()

    def show_analysis_description(self):
        description = """
        Downstream Analysis Description:
        --------------------------------
        This analysis includes:
        - Temporal Analysis of Variants: 
            Tracks the evolution of virus variants over time.

        - Geographical Distribution: 
            Studies the spread of variants across different regions.

        - Sequence Length Variation: 
            Analyzes the variation in sequence lengths among different samples.

        - Amino Acid Substitutions: 
            Identifies common substitutions in the virus sequences.

        - Submission Date Lag: 
            Examines the time lag between sample collection and submission dates.

        Enable this option to include these analyses in the final report.
        """
        QMessageBox.information(self, "Downstream Analysis Description", description)

    def show_instructions(self):
        instructions = """
        1. Provide GISAID Meta Table:  
            Click the "Browse" button next to "GISAID Meta Table" and 
            select the .tar.xz file you downloaded from https://gisaid.org/.

        2. Provide GenBank Meta Table:  
            Click the "Browse" button next to "GenBank Meta Table" and 
            select the CSV file you downloaded from 
            https://www.ncbi.nlm.nih.gov/labs/virus.

        3. Select Output Folder:  
            Select a folder where the matching table and 
            other result output files will be saved.

        4. GISAID FASTA Folder:  
            Choose an empty folder for the GISAID multi-sample FASTA file. 
            Ensure it's the only file in the folder when prompted to add it.

        5. Filter Options:  
            Optionally, you can select the checkboxes to filter by Omicron variant.

        5. Enable Downstream Analysis:  
            You can select the checkboxes to perform several downstream analysis: 
            [Temporal Analysis of Variants] [Geographical Distribution] 
            [Sequence Length Variation] [Amino Acid Substitutions] 
            [Submission Date Lag]

        6. Run:  
            Click the "Run" button to execute the script.

        7. Monitor Progress:  
            A second window with a progress bar and task description will appear. 
            This program demands substantial computational power and may briefly 
            suspend the GUI to dedicate more threads to processing, depending on 
            your hardware. 
            
        ** If the GUI becomes unresponsive, please do not close the application **

        """
        QMessageBox.information(self, "Instructions", instructions)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainApplication()
    window.show()

    with cProfile.Profile() as profiler:
        app.exec()

    profile_output_path = os.path.join(MainApplication.output_folder, "profile.prof")
    profiler.dump_stats(profile_output_path)
