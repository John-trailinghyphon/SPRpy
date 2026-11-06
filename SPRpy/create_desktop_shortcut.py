import platform
# Check if OS is Windows
if platform.system() != 'Windows':
    raise OSError('This script is only for Windows OS')

import os
import shutil
import sys
from win32com.client import Dispatch


def create_shortcut():
    """
    Creates two desktop shortcuts: "SPRpy" starts SPRpy, and "SPRpy folder" opens the SPRpy folder (config.toml,
    X-calibration files, conversion scripts and sessions). Works for installations with uv (uv tool install SPRpy) and pip.
    """

    package_dir = os.path.dirname(os.path.abspath(__file__))
    shell = Dispatch('WScript.Shell')
    desktop = shell.SpecialFolders('Desktop')

    # Shortcut that starts SPRpy: the "sprpy" command installed next to this Python environment (or on the PATH)
    scripts_dir = os.path.dirname(sys.executable)
    candidates = [os.path.join(scripts_dir, 'sprpy.exe'), os.path.join(scripts_dir, 'Scripts', 'sprpy.exe'), shutil.which('sprpy')]
    sprpy_command = next((c for c in candidates if c and os.path.exists(c)), None)
    if sprpy_command:
        shortcut = shell.CreateShortCut(os.path.join(desktop, 'SPRpy.lnk'))
        shortcut.TargetPath = sprpy_command
        shortcut.WorkingDirectory = package_dir
        shortcut.IconLocation = os.path.join(package_dir, 'assets', 'icon.ico')
        shortcut.Description = 'Start SPRpy'
        shortcut.save()
        print('Created desktop shortcut "SPRpy" (starts SPRpy)')
    else:
        print('The "sprpy" command was not found, so only the folder shortcut was created.')

    # Shortcut to the SPRpy folder
    shortcut = shell.CreateShortCut(os.path.join(desktop, 'SPRpy folder.lnk'))
    shortcut.TargetPath = package_dir
    shortcut.WorkingDirectory = package_dir
    shortcut.save()
    print('Created desktop shortcut "SPRpy folder" (' + package_dir + ')')


if __name__ == '__main__':
    create_shortcut()
