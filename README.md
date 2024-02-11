# GisGen V2: User Guide and README

Welcome to GisGen V2, the enhanced version of the GisGen tool designed for the harmonization of virus sequence metadata between GenBank and GISAID databases. This README provides comprehensive instructions on setting up your environment, running the **`GisGen_Run_withTimer.py`** script, and utilizing GisGen for your research needs.

## **Setting Up Your Environment**

### **Prerequisites**

- Python 3.11 or newer
- At least 32GB of RAM
- CPU equivalent or better than Intel Core i7-11700K for optimal performance

### **Using the YAML File for Environment Setup**

```

pip install biopython pandas matplotlib seaborn pyside6 rapidfuzz

```

1. **Download the YAML File**: Ensure you have the **`environment.yml`** file that was shared. This file contains all the necessary package information and Python version required for GisGen.
2. **Install Anaconda or Miniconda**: If you haven't already, download and install Anaconda or Miniconda. Both are excellent choices for managing environments and packages; Anaconda includes a large number of pre-installed packages, while Miniconda is a minimal installer for those who prefer to install only what they need.
3. **Create and Activate the Conda Environment**: Navigate to the directory containing the **`environment.yml`** file in your terminal or command prompt. Create the Conda environment by executing:
    
    ```bash
    
    conda env create -f environment.yml
    ```
    
    This command creates a new environment named as specified in the **`environment.yml`** file, installing Python and all required packages in their specified versions.
    
4. **Activate the Environment**: Once the environment is successfully created, you can activate it using:
    
    ```
    
    conda activate gisgen_env
    ```
    
    Replace **`gisgen_env`** with the actual name of your environment if it's different, as defined in the **`environment.yml`** file.
    

## **Running `GisGen_Run_withTimer.py`**

1. **Clone the Repository**: Clone the GisGen GitHub repository to your local machine:
    
    ```bash
    
    git clone https://github.com/ZooPhy/GisGen.git
    ```
    
    Navigate to the cloned directory:
    
    ```bash
    
    cd GisGen
    ```
    
2. **Execute the Script**: Run the **`GisGen_Run_withTimer.py`** script:
    
    ```
    
    python GisGen_Run_withTimer.py
    ```
    
    This script includes a timer for performance benchmarking and will provide real-time updates on the process.
    

## **Using GisGen**

### **Step 1: Data Preparation**

- Ensure you have the metadata tables from GenBank and GISAID. For GISAID data, follow the instructions provided in our conversation to download the multi-sample FASTA file.

### **Step 2: Configuration**

- Upon launching GisGen, you'll be prompted to configure settings such as the path to your data files and desired fuzzy matching parameters.

### **Step 3: Running the Tool**

- With your environment set and data prepared, initiate the tool's workflow. Expect the GUI to appear unresponsive at times, especially during heavy computational tasks like MD5 checksum generation. This is normal—GisGen is still running.

### **Step 4: Monitoring Progress**

- The tool processes up to 1,000 samples efficiently within an estimated 12 to 14 minutes, excluding the time needed to download the FASTA file from GISAID. Monitor the terminal or command prompt for real-time updates.

### **Step 5: Reviewing Outputs**

- Upon completion, GisGen generates detailed reports, including **`GIS_GEN_final.csv`** and **`md5_results.csv`**, highlighting matches, discrepancies, and MD5 checksum validations.

# GisGen Output files

The GisGen pipeline, designed for the harmonization and analysis of virus sequence metadata between GenBank and GISAID databases, generates several key files as part of its output. Each file serves a specific purpose, providing valuable insights and data for epidemiological research. Below is a detailed report of the files produced by the GisGen pipeline and the contents of each.

### **1. GIS_GEN_final.csv**

- **Description**: This is the primary output file of the GisGen pipeline. It contains the harmonized metadata for the analyzed samples, integrating data from both GenBank and GISAID sources.
- **Contents**: The file includes columns for sample identifiers, collection dates, host information, geographic location, sequence length, and any other metadata fields processed by GisGen. It also includes the results of fuzzy matching, showing how data from different sources correspond and differ.

### **2. md5_results.csv**

- **Description**: This file provides the results of the MD5 checksum validation process, which is essential for verifying the integrity and authenticity of the sequence data.
- **Contents**: It lists the MD5 checksums for each sequence file from both GISAID and GenBank, alongside a comparison result indicating whether the checksums match. This ensures the sequences are identical and have not been altered or corrupted.

### **3. Mismatch_Report.txt (or .csv)**

- **Description**: Generated during the reconciliation process, this report details the discrepancies found between the metadata entries in GenBank and GISAID.
- **Contents**: For each mismatch identified, the report lists the sample identifier, the fields that do not match, and the differing values from each database. This is crucial for identifying potential errors or inconsistencies in the metadata.

### **4. Temporal_Analysis_Variants.csv**

- **Optional Output**: If the user opts for downstream analysis, this file contains the results of the temporal analysis of virus variants.
- **Contents**: The file includes time-stamped data showing the distribution and frequency of different virus variants over time, providing insights into their evolution and spread.

### **5. Geographical_Distribution_Variants.csv**

- **Optional Output**: Another downstream analysis output, focusing on the geographic distribution of the virus variants.
- **Contents**: It maps the occurrence of virus variants across different locations, offering a spatial perspective on the data. This includes country and, if available, more specific location data such as state or city.

### **6. Sequence_Length_Variation.csv**

- **Optional Output**: This file results from analyzing the variation in sequence length among the samples.
- **Contents**: It provides statistics on the sequence lengths, such as average, median, and range, broken down by virus variant or other relevant classifications.

### **7. Amino_Acid_Substitutions.csv**

- **Optional Output**: If conducted, this analysis identifies common amino acid substitutions across the sequences.
- **Contents**: The file lists the substitutions observed, their locations in the genome, and their frequencies. This can be critical for understanding mutations and their potential impacts.

### **8. Submission_Date_Lag.csv**

- **Optional Output**: This file contains an analysis of the lag between sample collection dates and submission dates to the databases.
- **Contents**: It includes each sample's collection and submission dates, the calculated lag time, and statistical summaries of the lag across the dataset.

### **Additional Notes**

- Ensure the GISAID FASTA folder is empty before adding the newly downloaded multi-sample FASTA file to prevent any processing errors.
- Patience is required due to the tool's computational demands. The apparent unresponsiveness during heavy tasks is typical and does not indicate a crash.

For any issues or further assistance, please contact the development team or refer to the extensive support resources available on the GitHub page. This guide is designed to facilitate a smooth and effective user experience with GisGen, ensuring accurate and efficient genomic data analysis for epidemiological research.
