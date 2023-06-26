import tkinter
import os
import re
from tkinter.filedialog import askopenfilenames

tkinter.Tk().withdraw()


def find_time(input_files):
    # Use this to extract time from the .csv file generated when exporting spectra directly in the Bionavis SPR

    for file in input_files:
        head, tail = os.path.splitext(file)

        with open(file, 'r') as f:
            string = f.readline()

        string = string.replace(' L1 670nm;', '')
        string = string.replace(' L2 980nm;', '')
        string = string.replace(' L3 670nm;', '')
        string = string.replace(' L4 785nm;', '')

        #  Saving as new .csv file with time vector
        with open((head + '_time' + tail), 'w') as w:
            w.write(string)

    return 'Done!'


if __name__ == '__main__':
    print(find_time(askopenfilenames()))
