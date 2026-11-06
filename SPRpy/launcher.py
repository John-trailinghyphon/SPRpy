# Command line entry points of the installed SPRpy package (see [project.scripts] in pyproject.toml):
#   sprpy           starts SPRpy
#   sprpy-convert   converts .spr2 files to .csv (files can be given as arguments, otherwise a dialog opens)
#   sprpy-xcal      generates the X-calibration file of an instrument
#   sprpy-folder    opens the SPRpy folder (config.toml, X-calibration files, sessions) in the file explorer
#
# The SPRpy scripts read config.toml and other files from their own folder, so they are run from there.

import os
import runpy
import subprocess
import sys

PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))


def _run_script(script_name, arguments=()):
    os.chdir(PACKAGE_DIR)
    if PACKAGE_DIR not in sys.path:
        sys.path.insert(0, PACKAGE_DIR)
    sys.argv = [script_name, *arguments]
    runpy.run_path(os.path.join(PACKAGE_DIR, script_name), run_name='__main__')


def main():
    _run_script('SPRpy.py')


def convert():
    # Paths given on the command line are relative to where the command was run
    _run_script('SPRpy_spr2_to_csv.py', [os.path.abspath(path) for path in sys.argv[1:]])


def xcal():
    _run_script('SPRpy_X_cal.py')


def open_folder():
    print('SPRpy folder: ' + PACKAGE_DIR)
    if sys.platform.startswith('win'):
        os.startfile(PACKAGE_DIR)
    elif sys.platform == 'darwin':
        subprocess.run(['open', PACKAGE_DIR])
    else:
        subprocess.run(['xdg-open', PACKAGE_DIR])


if __name__ == '__main__':
    main()
