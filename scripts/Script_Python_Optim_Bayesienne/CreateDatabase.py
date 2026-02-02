#!/usr/bin/env python3
#  Copyright (C) 2024
#  Laboratoire d'Informatique Fondamentale et Appliquée de Tours, Tours, France
#
#  DIGEP, Politecnico di Torino, Corso Duca degli Abruzzi 24, Torino, Italy
#  This file is part of bilevel-scheduling.
#
#  bilevel-scheduling is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published
#  by the Free Software Foundation, either version 3 of the License,
#  or (at your option) any later version.
#
#  bilevel-scheduling is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty
#  of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with bilevel-scheduling. If not, see <https://www.gnu.org/licenses/>.

#!/usr/bin/env python3
import os
import shutil
import re
import pandas as pd

def concatenate_csv_files(folder_path):
    """
    In the same folder, there are CSV files of the form {integer}M_results{method}.csv.
    This function concatenates all CSV files that share the same method into a single CSV file.
    If only one method group exists, the output file is named 'resData.csv'; otherwise, each group is saved as 'resData_{method}.csv'.
    It returns a dictionary mapping each method to a tuple (concatenated DataFrame, output CSV file path).
    """
    pattern = re.compile(r'^(\d+)M_results(.+)\.csv$')
    files_by_method = {}
    for file in os.listdir(folder_path):
        if file.endswith(".csv"):
            match = pattern.match(file)
            if match:
                method = match.group(2)
                files_by_method.setdefault(method, []).append(os.path.join(folder_path, file))

    concatenated_dataframes = {}
    for method, file_list in files_by_method.items():
        file_list.sort()  # sort for consistency
        df_list = []
        for file_path in file_list:
            try:
                df = pd.read_csv(file_path,sep="\t")
                df_list.append(df)
            except Exception as e:
                print(f"Error reading {file_path}: {e}")
        if df_list:
            concatenated_df = pd.concat(df_list, ignore_index=True)
            # Determine the output filename: if only one method group exists overall, use "resData.csv"
            if len(files_by_method) == 1:
                output_filename = "resData.csv"
            else:
                output_filename = f"resData_{method}.csv"
            output_path = os.path.join(folder_path, output_filename)
            concatenated_df.to_csv(output_path, index=False)
            concatenated_dataframes[method] = (concatenated_df, output_path)
    return concatenated_dataframes

def get_instance_name(file_path):
    """
    Opens the instance file and reads the first line.
    The first line is expected to be of the form "name:instanceName".
    Returns the instance name (the part after the colon) or None if not found.
    """
    with open(file_path, 'r') as f:
        first_line = f.readline().strip()
        if first_line.startswith("name:"):
            return first_line.split("name:")[1].strip()
    return None

def get_objective_and_time_from_df(df, instance_name):
    """
    Searches in the DataFrame for the row where the "InstanceName" column matches the given instance_name.
    Returns a tuple (objective, time) based on the columns "Objective" and "Time".
    Returns (None, None) if no matching row is found.
    """
    matching_rows = df[df["InstanceName"] == instance_name]
    if not matching_rows.empty:
        row = matching_rows.iloc[0]
        return row.get("Objective"), row.get("Time")
    return None, None

def main(genValidation = True):
    # Get the current working directory and define the base path for "../../instances/"
    # The os.path.abspath call converts the relative path to an absolute one.
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../instances"))
    print("Base directory:", base_dir)

    # Define the paths for the train and result databases within the base directory.
    train_dir = os.path.join(os.path.dirname(__file__), "TrainInstances")
    res_train_dir = os.path.join(os.path.dirname(__file__), "TrainResults")
    valid_dir = os.path.join(os.path.dirname(__file__), "ValidInstances")
    res_valid_dir = os.path.join(os.path.dirname(__file__), "ValidResults")

    # Create "trainDatabase" and "resDatabase" folders if they do not exist.
    os.makedirs(train_dir, exist_ok=True)
    os.makedirs(res_train_dir, exist_ok=True)
    os.makedirs(valid_dir, exist_ok=True)
    os.makedirs(res_valid_dir, exist_ok=True)

    # Loop over all subdirectories in base_dir that start with 'N' and are followed by digits.
    for folder in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder)
        # Check if the folder name starts with "N" and the rest is numeric (e.g., "N100", "N200", etc.)
        if os.path.isdir(folder_path) and folder.startswith("N") and folder[1:].isdigit():
            # Extract the numeric value after 'N' for use in the new file name.
            n_value = folder[1:]

            # In the same folder, there are CSV files of the form {integer}M_results{method}.csv.
            # Create a single CSV file by concatenating all CSV files with the same method.
            concatenated_csv_data = concatenate_csv_files(folder_path)
            # For further processing, if multiple methods exist, choose the first one (adjust as needed).
            if concatenated_csv_data:
                first_method = list(concatenated_csv_data.keys())[0]
                df, res_csv_path = concatenated_csv_data[first_method]
            else:
                df = None

            # Construct the path to the "instances" folder inside the current N* directory.
            instances_folder = os.path.join(folder_path, "instances")

            # Check if the "instances" folder exists before processing.
            if os.path.isdir(instances_folder):
                # Initialize an instance counter for the current folder.
                instance_train_number = 1
                instance_valid_number = 1

                # Loop through all files in the "instances" folder.
                for file_name in sorted(os.listdir(instances_folder)):

                    file_path = os.path.join(instances_folder, file_name)
                    # Process only files (skip directories if any)
                    if os.path.isfile(file_path):

                        # Retrieve the instance name by opening the file and reading the first line.
                        # The line is expected to be of the form "name:instanceName".
                        instance_name = get_instance_name(file_path)
                        indexInstance = int(instance_name.split('_')[0].split("instance")[1])
                        # Construct the destination file name with the format "donnees_{n_value}_{instance_number}.dat"
                        # if genValidation is true, then save one instance on train and other on validation
                        if (genValidation):
                            if (indexInstance % 2 == 0):
                                dest_file_name = f"donnees_{n_value}_{instance_valid_number}.dat"
                                dest_file_path = os.path.join(valid_dir, dest_file_name)
                            else:
                                dest_file_name = f"donnees_{n_value}_{instance_train_number}.dat"
                                dest_file_path = os.path.join(train_dir, dest_file_name)
                        else:
                            dest_file_name = f"donnees_{n_value}_{instance_train_number}.dat"
                            dest_file_path = os.path.join(train_dir, dest_file_name)

                        # Copy the instance file to the trainDatabase folder with the new name.
                        shutil.copy2(file_path, dest_file_path)

                        if instance_name and df is not None:
                            # Search in the concatenated DataFrame for the row with the matching instance name in the "InstanceName" column.
                            objective, exec_time = get_objective_and_time_from_df(df, instance_name)

                            if objective is not None and exec_time is not None:
                                # Create a file in the res_train_dir folder with the name "donnees_{n_value}_{instance_number}.dat.seq"
                                # The file contains the text "{objective} {time}".
                                # if genValidation is true, then save one instance on train and other on validation
                                if (genValidation):
                                    if (indexInstance % 2 == 0):
                                        seq_file_name = f"donnees_{n_value}_{instance_valid_number}.dat.seq"
                                        seq_file_path = os.path.join(res_valid_dir, seq_file_name)
                                    else:
                                        seq_file_name = f"donnees_{n_value}_{instance_train_number}.dat.seq"
                                        seq_file_path = os.path.join(res_train_dir, seq_file_name)
                                else:
                                    seq_file_name = f"donnees_{n_value}_{instance_train_number}.dat.seq"
                                    seq_file_path = os.path.join(res_train_dir, seq_file_name)

                                with open(seq_file_path, 'w') as seq_file:
                                    seq_file.write(f"{objective} {exec_time}")
                            else:
                                print(f"Instance {instance_name} not found in the DataFrame.")
                        else:
                            print(f"Could not retrieve instance name from {file_path} or concatenated DataFrame not available.")

                        # Increment the instance counter.
                        if (genValidation):
                            if (indexInstance % 2 == 0):
                                instance_valid_number += 1
                            else:
                                instance_train_number += 1
                        else:
                            instance_train_number += 1

                print(f"{instance_valid_number + instance_train_number - 2} files was created for instance with {n_value} jobs")
                if (genValidation):
                    print(f"{(instance_valid_number + instance_train_number - 1)//2} files was created train with {n_value} jobs")
                    print(f"{(instance_valid_number + instance_train_number - 1)//2} files was created validation with {n_value} jobs")
            else:
                print(f"Warning: The folder {instances_folder} does not exist.")

if __name__ == "__main__":
    main(True)


