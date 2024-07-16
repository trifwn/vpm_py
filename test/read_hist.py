import numpy as np
from scipy.io import FortranFile, FortranEOFError
import matplotlib.pyplot as plt

def read_hist_file(file_pattern, num_files):
    data = []
    for i in range(1, num_files+1):
        filename = f'{i:02d}{file_pattern}'
        read_meta = True
        file_data = []
        with FortranFile(filename, 'r') as f:
            print(f'Reading file {filename}... {i}')
            # Read the initial metadata if this is the first file
            while True:
                if read_meta:
                    XMIN_pm, YMIN_pm, ZMIN_pm = f.read_reals(np.float64)
                    DXpm, DYpm, DZpm = f.read_reals(np.float64)                
                    NXs_bl ,NXf_bl , j ,k = f.read_ints(np.int32)
                    print(f'Initial metadata: XMIN_pm={XMIN_pm}, YMIN_pm={YMIN_pm}, ZMIN_pm={ZMIN_pm}')
                    print(f'DXpm={DXpm}, DYpm={DYpm}, DZpm={DZpm}')
                    print(f'NXs_bl={NXs_bl}, NXf_bl={NXf_bl}, j={j}, k={k}')
                    read_meta = False
                    metadata = (
                        XMIN_pm, YMIN_pm, ZMIN_pm,
                        DXpm, DYpm, DZpm,
                        NXs_bl, NXf_bl, j, k
                    )
                    file_data.append(metadata)
                
                # Read NTIME_pm and velocity data for the current file
                try:
                    record = f.read_record(
                      '<i4',  # NTIME_pm
                      np.dtype(('<f8', (NXf_bl - NXs_bl + 1))), # velvrx_pm
                      np.dtype(('<f8', (NXf_bl - NXs_bl + 1))), # velvry_pm
                      np.dtype(('<f8', (NXf_bl - NXs_bl + 1)))  # velvrz_pm
                    )
                    NTIME_pm = record[0]
                    velvrx_pm = record[1]
                    velvry_pm = record[2]
                    velvrz_pm = record[3]
           
                    record = (NTIME_pm, velvrx_pm, velvry_pm, velvrz_pm)
                    file_data.append(record)
                    # print(
                    #     f'Read NTIME_pm={NTIME_pm} and velocity (velx {velvrx_pm.shape}, vely: {velvry_pm.shape} and velz: {velvrz_pm.shape}) data for file {filename}.'
                    # )
                # Check EOF
                except FortranEOFError:
                    print(f'EOF reached for file {filename}.')
                    print()
                    break

            # print(f'Read file {filename}: NTIME_pm={NTIME_pm}')
        data.append(file_data)
    
    return data

class DataSlice():
    def __init__(self, NTIME,j,k, velvrx, velvry, velvrz):
        self.NTIME = NTIME
        self.j = j
        self.k = k
        self.velvrx = velvrx
        self.velvry = velvry
        self.velvrz = velvrz
    
    def __str__(self):
        print(f'NTIME={self.NTIME}, j={self.j}, k={self.k}\n')

    def __repr__(self):
        return f'NTIME={self.NTIME}, j={self.j}, k={self.k}\n'


class allData():
    def __init__(self):
        self.slices = []
    
    def add_slice(self, slice):
        self.slices.append(slice)

    def add_slices_from_Parser(self, parser):
        for slice in parser.slices:
            self.slices.append(slice)

    def get_data_by_jk(self, j, k):
        # Get all slices where j and k are equal to the given values
        # And order them by NTIME in ascending order
        slices = []
        for slice in self.slices:
            if slice.j == j and slice.k == k:
                slices.append(slice)
        
        return np.array(sorted(slices, key=lambda x: x.NTIME))
    
    def get_data_by_NTIME(self, NTIME):
        # Get all slices where NTIME is equal to the given value
        # And order them by j and k
        slices = []
        for slice in self.slices:
            if slice.NTIME == NTIME:
                slices.append(slice)
        
        return np.array(sorted(slices, key=lambda x: (x.j, x.k)))

class DataParser():
    jks: list[tuple[int, int]] = []
    def __init__(self, metadata, velocity_records):
        self.metadata = metadata
        self.velocity_data = velocity_records

        # From the metadata, get the initial values
        self.XMIN_pm = metadata[0]
        self.YMIN_pm = metadata[1]
        self.ZMIN_pm = metadata[2]
        self.DXpm = metadata[3]
        self.DYpm = metadata[4]
        self.DZpm = metadata[5]
        self.NXs_bl = metadata[6]
        self.NXf_bl = metadata[7]
        self.j = metadata[8]
        self.k = metadata[9]

        DataParser.jks.append((self.j, self.k))
        self.slices = []
        for record in velocity_records:
            NTIME_pm = record[0]
            velvrx_pm = record[1]
            velvry_pm = record[2]
            velvrz_pm = record[3]
            self.slices.append(
                DataSlice(NTIME_pm,self.j,self.k, velvrx_pm, velvry_pm, velvrz_pm)
            )

# Example usage:
if __name__ == "__main__":
    file_pattern = 'hist.bin'
    num_files = 9
    data = read_hist_file(file_pattern, num_files)

    print(len(data))
    # From each data et the initial metadata and the velocity data

    Data = allData()
    for i, file_data in enumerate(data):
        print(f'File {i+1}:')
        metadata = file_data[0]
        velocity_data = file_data[1:]
        data = DataParser(metadata, velocity_data)
        Data.add_slices_from_Parser(data)

    arr = Data.get_data_by_jk(4,4)
    print(len(arr))
    print(arr.shape)
    fig, ax = plt.subplots(3,1)
    for slice in arr:
        ax[0].plot(slice.velvrx)
        ax[1].plot(slice.velvry)
        ax[2].plot(slice.velvrz)
    fig.show()
    plt.show(block= True)

    arr = Data.get_data_by_NTIME(60)
    print(len(arr))
    print(arr.shape)
    fig, ax = plt.subplots(3,1)
    for slice in arr:
        ax[0].plot(slice.velvrx, label = f'NTIME={slice.NTIME}, j = {slice.j}, k = {slice.k}')
        ax[1].plot(slice.velvry, label = f'NTIME={slice.NTIME}, j = {slice.j}, k = {slice.k}')
        ax[2].plot(slice.velvrz, label = f'NTIME={slice.NTIME}, j = {slice.j}, k = {slice.k}')

    fig.show()
    plt.show(block= True)
    print('Successfully read all files.')