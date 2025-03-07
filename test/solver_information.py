
import pandas as pd
import matplotlib.pyplot as plt
import os

def read_csv_file(file_path):
    return pd.read_csv(file_path)

def main():

    parent_folder = "/mnt/c/Users/tryfonas/Data/"
    folders = os.listdir(parent_folder)
    folders = [os.path.join(parent_folder, folder) for folder in folders if folder.startswith('hill_vortex')]
    print(folders)

    file_name = 'solve_information.csv'
    all_dfs = {}
    
    # Read all dataframes
    for folder in folders:
        print(folder)
        file_path = os.path.join(folder, file_name)
        df = pd.read_csv(file_path)
        all_dfs[folder] = df
        print(f'Read {file_path} with shape {df.shape}')
    # Get metrics from first dataframe (assuming all have same columns)
    metrics = list(all_dfs[folders[0]].columns[1:])  # All columns except 'Iteration'
    print(metrics)

    
    # Create one plot per metric, including all datasets
    for metric in metrics:
        plt.figure(figsize=(12, 8))
        
        for folder, df in all_dfs.items():
            # Create a more readable label by removing 'results_' prefix
            label = folder.replace('results_', '')
            plt.plot(df['Iteration'], df[metric], marker='o', label=label)
        
        plt.title(f'{metric} over iterations')
        plt.xlabel('Iteration')
        plt.ylabel(metric)
        plt.grid(True)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()  # Adjust layout to prevent legend cutoff
        plt.show()

if __name__ == "__main__":
    main()
