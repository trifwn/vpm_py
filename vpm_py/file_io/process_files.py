import multiprocessing as mp
from .particle_file import process_particle_file_dat, process_particle_file_h5 
from .mesh_file import process_pm_file_h5, process_pm_file_dat
import os
import re 

def process_pm_file(filename: str, folder: str | None = None):
    if filename.endswith('.h5'):
        return process_pm_file_h5(filename, folder)
    else:
        return process_pm_file_dat(filename, folder)

def process_particle_file(filename: str, folder: str | None = None):
    if filename.endswith('.h5'):
        return process_particle_file_h5(filename, folder)
    else:
        return process_particle_file_dat(filename, folder)

def process_multiple_files(files, folder, process_func):
    with mp.Pool() as pool:
        results = pool.starmap(process_func, [(f, folder) for f in files])
    return results

def process_multiple_pm_files(files, folder):
    return process_multiple_files(files, folder, process_pm_file)

def process_multiple_particle_files(files, folder):
    return process_multiple_files(files, folder, process_particle_file)

def get_latest_particle_file(
    folder : str,
    pattern: str ='particles.*',
    extension: str ='.h5',
):
    files = os.listdir(folder)
    particle_filename_pattern = pattern + extension 
    particle_filename_pattern = particle_filename_pattern.replace('.', '\.')
    particle_regex = re.compile(pattern=particle_filename_pattern)
    files = sorted([f for f in files if particle_regex.search(f)])
    if len(files) == 0:
        if not folder.endswith('results'):
            return get_latest_particle_file(
                os.path.join(folder, "results/"), 
                pattern, 
                extension
            )

        return None
    latest_file = files[-1]
    print(f"Latest file: {latest_file}")
    return *process_particle_file(filename= latest_file, folder= folder), len(files)
