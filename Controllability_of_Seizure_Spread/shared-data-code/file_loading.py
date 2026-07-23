


def get_data(zip_file_directory,txt_file):
    from zipfile import ZipFile
    import numpy as np
    with ZipFile(zip_file_directory) as myzip:
        with myzip.open(txt_file) as myfile:
            f=np.loadtxt(myfile)
    return f

