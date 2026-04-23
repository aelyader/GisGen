import gc
import hashlib
import os
import tarfile
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO

import matplotlib.pyplot as plt
import pandas as pd
import requests
import seaborn as sns
from Bio import SeqIO
from PySide6.QtCore import QThread, Signal
from rapidfuzz import fuzz
from requests.adapters import HTTPAdapter
from urllib3 import Retry

try:
    from memory_profiler import profile
except ImportError:
    def profile(func):
        return func



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
        # Initialize GISAID as an empty DataFrame
        GISAID = pd.DataFrame()
        # Define the data types for each column
        dtype_dict = {
            'Virus name': 'object',  # Keeping as object as it's likely a unique identifier or string
            'Last vaccinated': 'category',
            'Passage details/history': 'category',
            'Type': 'category',
            'Accession ID': 'object',  # Keeping as object (likely a unique string)
            'Collection date': 'object',  # Consider parsing as datetime if needed
            'Location': 'category',
            'Additional location information': 'category',
            'Sequence length': 'int32',
            'Host': 'category',
            'Patient age': 'category',  # If numeric, consider 'int32' or 'float32'
            'Gender': 'category',
            'Clade': 'category',
            'Pango lineage': 'category',
            'Pango version': 'category',
            'Variant': 'category',
            'AA Substitutions': 'object',  # Keeping as object if it's a complex string
            'Submission date': 'category',  # Consider parsing as datetime if needed
            'Is reference?': 'float32',  # Consider 'boolean' if it's just True/False
            'Is complete?': 'category',
            'Is high coverage?': 'category',
            'Is low coverage?': 'category',
            'N-Content': 'float32',
            'GC-Content': 'float32'
        }

        # Check if the file is a .tar.xz file
        if self.gisaid_data_file.endswith('.tar.xz'):
            self.custom_print("The .tar.xz file exceeds 700MB, indicating that the TSV table size surpasses 12GB. "
                              "(The extraction and TSV upload process will require approximately 3 to 4 minutes.)")
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
                    GISAID = pd.read_csv(extracted_file_path, sep='\t', dtype=dtype_dict, low_memory=False)
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
        # Step 6: Starting fuzzy matching
        self.custom_print("Step 6: Starting fuzzy matching")
        fuzzy_matches = [self.fuzzy_match(gis_iso, gen_iso) for gis_iso, gen_iso in zip(gisaid_isolates, genbank_isolates)]

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
        def check_network_connection(url):
            try:
                response = requests.head(url, timeout=5)
                return response.ok
            except requests.RequestException as e:
                return False

        def download_data_with_retry(url):
            session = requests.Session()
            retries = Retry(total=5, backoff_factor=1, status_forcelist=[502, 503, 504])
            session.mount('http://', HTTPAdapter(max_retries=retries))
            session.mount('https://', HTTPAdapter(max_retries=retries))

            try:
                self.custom_print("Attempting to download GenBank FASTA MD5 table...")
                response = session.get(url)
                response.raise_for_status()  # This will raise an HTTPError if the response was an error
                return response.content
            except requests.exceptions.HTTPError as errh:
                self.custom_print("Http Error:", errh)
            except requests.exceptions.ConnectionError as errc:
                self.custom_print("Error Connecting:", errc)
            except requests.exceptions.Timeout as errt:
                self.custom_print("Timeout Error:", errt)
            except requests.exceptions.RequestException as err:
                self.custom_print("Oops: Something Else", err)

        self.custom_print("Step 11: Download GenBank FASTA MD5 table")
        genbank_md5_table_url = "https://zenodo.org/record/10578584/files/updated_genbank_sequences.csv"
        if check_network_connection(genbank_md5_table_url):
            genbank_md5_table_content = download_data_with_retry(genbank_md5_table_url)
            # Assuming genbank_md5_table_content is a CSV content, you would then load it into a DataFrame as before
            genbank_md5_table = pd.read_csv(StringIO(genbank_md5_table_content.decode('utf-8')))
        else:
            self.custom_print("Network connection to GenBank MD5 table URL could not be established.")
        gc.collect()
        self.progress_update.emit(75)
        # Compare MD5 checksums
        self.custom_print("Step 12: Compare MD5 checksums between GISAID & GenBank")
        genbank_md5_lookup = dict(zip(genbank_md5_table['Accession'], genbank_md5_table['GenBank_MD5']))
        md5_rows = []

        for _, row in self.GIS_GEN_final.iterrows():
            gisaid_id = row['AccessionID_GIS']
            genbank_id = row['AccessionID_GEN']
            gisaid_md5 = md5_results.get(gisaid_id)
            genbank_md5 = genbank_md5_lookup.get(genbank_id)
            is_match = 'Yes' if gisaid_md5 == genbank_md5 else 'No'

            md5_rows.append({
                'GISAID_ID': gisaid_id,
                'GenBank_ID': genbank_id,
                'GISAID_MD5': gisaid_md5,
                'GenBank_MD5': genbank_md5,
                'Match': is_match,
            })

        md5_comparison_results = pd.DataFrame(md5_rows)

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
        variant_counts = df.groupby(['Collection_Date_GisGen', 'Variant_GIS'], observed=True).size().unstack(fill_value=0)
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
        variant_distribution = df.groupby(['Country_GIS', 'Variant_GIS'], observed=True).size().unstack(fill_value=0)
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
            if pd.isna(x) or pd.isna(choice):
                return None

            left = str(x).strip()
            right = str(choice).strip()

            if left == right:
                return right

            score = fuzz.token_set_ratio(left, right)
            return right if score >= 90 else None
        except Exception as e:
            self.custom_print(f"An error occurred: {e}")
            return None

