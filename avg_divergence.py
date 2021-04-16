import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.ticker as ticker


def load_2_d_array(directory, col_number_y):
    px_axis = np.loadtxt(directory, skiprows=3, usecols=(0,))

    counts = np.loadtxt(directory, skiprows=3, usecols=(col_number_y,))
    return np.transpose(np.stack((px_axis, counts), axis=0))


def mean_std(array):
    mean = np.mean(array, axis=1)
    std = np.std(array, axis=1)
    return mean, std


def label_from_header(file_name):
    return file_name[0:4]


def plot_results(array, label):
    plt.scatter(array[::, 0], array[::, 1])
    plt.xlabel('N')
    plt.ylabel('sigma 1/e divergence [mrad]')
    #plt.legend()


def get_file_list(path_picture):
    tif_files = []
    counter = 0
    for file in os.listdir(path_picture):
        print(file)
        try:
            if file.endswith(".txt"):
                tif_files.append(str(file))
                counter = counter + 1
            else:
                print("only other files found")
        except Exception as e:
            raise e
    return tif_files


directory = "z1600/"
files_1600 = get_file_list(directory)


def batch_files(files, path):
    array_all = np.zeros([15, len(files)])
    count = 0
    for x in files:
        directory = str(path) + str(x)
        print(directory)
        array = load_2_d_array(directory, 1)
        array_all[::, count] = array[::, 1]
        label = str(x[:-4])
        #plot_results(array, label)
        count = count + 1

    return array_all, array[::,0]


print(directory)

fig, ax = plt.subplots()




array_all, array_x = batch_files(files_1600, directory)
mean_array, std_array = mean_std(array_all)
ax.axes.errorbar(array_x,  mean_array, yerr=std_array, xerr=None,fmt='o', label = '1600 um')


directory = "z2000/"
files_2000 = get_file_list(directory)
array_all, array_x = batch_files(files_2000, directory)
mean_array, std_array = mean_std(array_all)
ax.axes.errorbar(array_x,  mean_array, yerr=std_array, xerr=None,fmt='s', label = '2000 um')

directory = "z2800/"
files_2800 = get_file_list(directory)
array_all, array_x = batch_files(files_2800, directory)
mean_array, std_array = mean_std(array_all)
ax.axes.errorbar(array_x,  mean_array, yerr=std_array, xerr=None,fmt='v', label = '2800 um')

directory = "z3200/"
files_3200 = get_file_list(directory)
array_all, array_x = batch_files(files_3200, directory)
mean_array, std_array = mean_std(array_all)
ax.axes.errorbar(array_x,  mean_array, yerr=std_array, xerr=None, fmt='.', label = '3200 um')

directory = "z3600/"
files_3600 = get_file_list(directory)
array_all, array_x = batch_files(files_3600, directory)
mean_array, std_array = mean_std(array_all)
ax.axes.errorbar(array_x,  mean_array, yerr=std_array, xerr=None, fmt='x',label = '3600 um')


def diffraction_limit():
    N_list = np.arange(1, 30, 1)
    N_diffraction_limit = np.zeros([29, 1])
    for x in range(15, 30 - 1):
        # halfangle
        N_list[x] = 1 + x
        N_diffraction_limit[x] = (60. *1E3/ 1500.) / (1 + x)
    return N_list, N_diffraction_limit


#plot_diffraction_limit()
array_x, diffraction_limit = diffraction_limit()
ax.axes.scatter(array_x,  diffraction_limit, label = 'ThetaL/N', alpha = 0.5)

#plt.legend()
ax.axes.set_xlabel('harmonic number /N')
ax.axes.set_ylabel('w0(z,N) divergence /mrad')
ax.axes.set_yscale('log')
#plt.axes.axes.set_yscale(1E-3)
ax.axes.set_xlim(16.5,27)
ax.axes.set_ylim(1,6)
minors = np.linspace(0, 6, 11)[1:-1]
ax.yaxis.set_minor_locator(ticker.FixedLocator(minors))
ax.axes.legend()

plt.savefig("div_over_N_best_shots" + ".png", bbox_inches="tight", dpi=1000)
plt.show()
