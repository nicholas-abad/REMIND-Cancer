# REMIND-Cancer: Folder Structure & Results JSON Setup

This script automates the setup of the **REMIND-Cancer** project’s folder structure and the creation of a JSON file to track results. It takes a **configuration file** as input and performs the following:

1. **Creates a structured folder system** for patient data.
2. **Copies WGS (Whole Genome Sequencing) VCF files** into the structured folder.
3. **Generates a results JSON file** to track processed data.

---

## **Optional: Run Validation Scripts Before This**

Before running the main pipeline, ensure everything is set up correctly by validating the required files. Navigate to the `utils` subfolder and run the following scripts:

### **Validate Configuration File**

```bash
python src/utils/validate_configuration_file.py -c path/to/configuration_file.json
```

### **Validate Metadata File**

```bash
python src/utils/validate_metadata_file.py -m path/to/metadata.csv
```

### **Validate RNA-Seq Dataframe**

```bash
python src/utils/validate_rnaseq_dataframe.py -r path/to/rnaseq.csv
```

These scripts will check for missing files, incorrect formats, and inconsistencies to prevent errors in the pipeline.

---

## **Usage**

Run the script in order to set up the pipeline structure using:

```bash
python src/pipeline_setup/create_initial_structure.py -c path/to/configuration_file.json
```

---

## **What This Script Does**

### **1. Load Configuration File**

The script begins by reading a JSON configuration file, which contains:

- The dataset name.
- The path to the **metadata CSV file**.
- The output location for patient folders.
- The output path for the results JSON file.

### **2. Create Folder Structure for Patients**

This step reads the **metadata CSV file** and generates a directory structure:

```
output_path_to_patient_folders/
 ├── patient_1_tumor/
 │   ├── sample1_original.vcf
 ├── patient_2_tumor/
 │   ├── sample2_original.vcf
```

- Each **patient’s folder** is named as `pid_tumor_origin`.
- Each **VCF file** is copied and renamed with `_original.vcf`.
- If a **VCF file is missing**, a warning is logged.

### **3. Create a Results Tracking JSON**

A JSON file is generated to track the presence of sequencing and expression data.

#### **Example Output (`results.json`)**

```json
{
    "results": {
        "original": {
            "primary_tumor_wgs": ["path/to/patient_1/sample1_original.vcf"],
            "primary_tumor_wgs_and_rnaseq": [],
            "metastasic_tumor_wgs": [],
            "metastasic_tumor_wgs_and_rnaseq": []
        }
    }
}
```

The script:

- Extracts `pid` from each patient folder.
- Checks `tumor_origin` and whether RNA-seq data exists.
- Classifies files accordingly.

## **Detailed Breakdown of Functions**

### **`create_folder_structure(metadata_path, output_folder)`**

- Reads the metadata CSV.
- Iterates through each **patient ID (`pid`)**.
- Creates a **patient subfolder** if it doesn’t exist.
- Copies the **VCF file** to the correct location.
- Logs **missing files**.

### **`create_results_json(metadata_path, patient_folders_path, results_json_path)`**

- Reads metadata to extract tumor and expression data.
- Scans patient folders for VCF files.
- Categorizes each file into one of four groups:
  - `primary_tumor_wgs`
  - `primary_tumor_wgs_and_rnaseq`
  - `metastasic_tumor_wgs`
  - `metastasic_tumor_wgs_and_rnaseq`
- Saves results as a JSON file.

### **`main(config_path)`**

- Parses the **configuration file**.
- Calls `create_folder_structure()`.
- Calls `create_results_json()`.

---

## **Expected Console Output**

When you run the script, you’ll see output similar to:

```bash
#### (1/2) Creating Folder Structure ####
Setting up folder structure: 100%|██████████████| 50/50 [00:03<00:00, 15.5it/s]
✅ Folder structure created at: /path/to/patient_folders
❗ Number of VCF files not copied: 3

#### (2/2) Creating Results JSON File ####
Generating JSON file: 100%|██████████████| 50/50 [00:02<00:00, 20.2it/s]
✅ Results JSON saved at: /path/to/results.json
```

---

## **Prerequisites**

- **Python 3.x** installed.
- Install dependencies using:
  ```bash
  pip install -r requirements.txt
  ```
- Ensure all file paths in `configuration_file.json` are correctly set.

---

## **Troubleshooting**

| Issue                                          | Solution                                                                               |
| ---------------------------------------------- | -------------------------------------------------------------------------------------- |
| `FileNotFoundError: Metadata file not found` | Check that the metadata CSV path is correct in the config file.                        |
| `VCF file does not exist`                    | Ensure the `path_to_wgs_file` column in the metadata file contains valid file paths. |
| `Results JSON is empty`                      | Check that the metadata contains valid tumor and RNA-seq data references.              |

This script automates the setup process, ensuring that data is structured properly before pipeline execution. 🚀
