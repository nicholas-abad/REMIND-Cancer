import os
import json
import argparse

def get_project_root():
    """
    Finds the project root by navigating up from the script location.

    Returns
    -------
    str
        Absolute path to the project root directory.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(script_dir, "../../../"))

def resolve_path(relative_path, is_executable=False):
    """
    Converts a relative path to an absolute path based on the project root, unless it's an executable.

    Parameters
    ----------
    relative_path : str
        The relative file path to be resolved.
    is_executable : bool, optional
        Whether the path refers to an executable file (default is False).

    Returns
    -------
    str
        Absolute path of the given relative path.
    """
    if is_executable:
        return os.path.expanduser(relative_path)  # Use home directory expansion for executables
    project_root = get_project_root()
    return os.path.abspath(os.path.join(project_root, relative_path))

def check_paths_in_config(config, parent_key=""):
    """
    Recursively checks file paths in the configuration dictionary.

    Parameters
    ----------
    config : dict
        The configuration dictionary containing file paths.
    parent_key : str, optional
        The parent key used for tracking nested dictionary paths (default is an empty string).

    Returns
    -------
    tuple
        A tuple containing two lists:
        - missing_files : list of str
            List of messages indicating missing file paths.
        - existing_files : list of str
            List of messages indicating existing file paths.
    """
    missing_files = []
    existing_files = []
    
    for key, value in config.items():
        if isinstance(value, dict):
            sub_missing, sub_existing = check_paths_in_config(value, key)
            missing_files.extend(sub_missing)
            existing_files.extend(sub_existing)
        elif isinstance(value, str) and key != "path_to_results":  # Skip path_to_results
            if "/" in value or "\\" in value:  # Likely a file path
                is_executable = key == "path_to_fimo_executable"
                abs_path = resolve_path(value, is_executable)
                
                if not os.path.exists(abs_path):
                    missing_files.append(f"Missing: {abs_path} (Key: {parent_key}.{key})")
                elif is_executable and not os.access(abs_path, os.X_OK):
                    missing_files.append(f"Not an executable: {abs_path} (Key: {parent_key}.{key})")
                else:
                    existing_files.append(f"Exists: {abs_path} (Key: {parent_key}.{key})")

    return missing_files, existing_files

def validate_config(config_path):
    """
    Loads the configuration file and checks for missing and existing paths.

    Parameters
    ----------
    config_path : str
        Path to the configuration file.
    """
    project_root = get_project_root()
    abs_config_path = os.path.abspath(config_path)

    # Ensure configuration file exists
    if not os.path.exists(abs_config_path):
        print(f"❌ Configuration file not found: {abs_config_path}")
        return

    with open(abs_config_path, "r") as file:
        config = json.load(file)

    # Check all paths except "path_to_results"
    missing_files, existing_files = check_paths_in_config(config)

    # Print results
    if existing_files:
        print("\n✅ Existing files:")
        for msg in existing_files:
            print(msg)

    if missing_files:
        print("\n🚨 Missing or incorrect files detected:")
        for msg in missing_files:
            print(msg)
    else:
        print("\n✅ All paths are valid!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validate file paths in a configuration JSON file.")
    parser.add_argument("-c", "--config", required=True, type=str, help="Path to the configuration file.")
    
    args = parser.parse_args()
    
    validate_config(args.config) 