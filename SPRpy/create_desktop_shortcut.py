import platform
# Check if OS is Windows
if platform.system() != 'Windows':
    raise OSError('This script is only for Windows OS')

import sysconfig
from win32com.client import Dispatch
import os
from __about__ import version


def create_shortcut():

    # folder name
    folder_name = 'SPRpy' + ' (v' + version + ')'

    # Scripts directory (location of python packages)
    lib_dir = sysconfig.get_path('purelib')

    # Target of shortcut
    target = os.path.join(lib_dir, folder_name)

    # Name of link file
    link_name = folder_name + '.lnk'

    # Path to location of link file
    shell = Dispatch('WScript.Shell')
    shell_folder = shell.SpecialFolders('Desktop')
    shortcut = shell.CreateShortCut(shell_folder + '\\' + link_name)
    shortcut.TargetPath = target
    shortcut.WorkingDirectory = lib_dir
    shortcut.save()
