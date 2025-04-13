# Validation Scripts for REMIND-Cancer

This subfolder contains three validation scripts designed to ensure the integrity and correctness of configuration, metadata, and RNA-seq dataframe files used in the **REMIND-Cancer** project.

---

## **1. `validate_configuration_file.py`**

This script validates paths in the configuration file to ensure they exist and are correctly formatted.

### **Usage:**

```bash
python src/utils/validate_configuration_file.py -c path/to/configuration_file.json
```

### **Checks Performed:**

- Ensures the configuration file exists.
- Verifies that paths (except `path_to_results`) are valid.
- Checks if `path_to_fimo_executable` is an executable file.
- Reports missing or incorrectly formatted paths.

---

## **2. `validate_metadata_file.py`**

This script validates a metadata CSV file by checking required columns and ensuring referenced files exist.

### **Usage:**

```bash
python src/utils/validate_metadata_file.py -m path/to/metadata.csv
```

### **Checks Performed:**

- Ensures the metadata file exists.
- Confirms required columns exist:`pid`, `tumor_origin`, `path_to_wgs_file`, `cohort`, `ge_data_available`.
- Verifies whether all `path_to_wgs_file` entries exist.
- Checks the validity of corresponding VCF files (expects `.vcf` filenames based on `.wgs` entries).
- Prints the number of valid WGS and VCF files.

---

## **3. `validate_rnaseq_dataframe.py`**

This script validates an RNA-seq dataframe to ensure correct formatting and column structure.

### **Usage:**

```bash
python src/utils/validate_rnaseq_dataframe.py -r path/to/rnaseq.csv
```

### **Checks Performed:**

- Ensures the RNA-seq file exists.
- Confirms that the file is comma-separated.
- Verifies the presence of the `pid` column.
- Counts and prints the number of **PIDs** and **genes**.

---

## **General Notes**

- Ensure you have **Python 3.x** installed.
- Install dependencies using:
  ```bash
  pip install -r requirements.txt
  ```
- Modify file paths accordingly before running scripts.

These validation scripts help maintain data integrity and prevent issues before running further analyses.
