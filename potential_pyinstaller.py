# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 09:09:45 2023

@author: lawashburn
"""

# Python To Exe
# pip install pyinstaller
import os
import sys
# convert to Console based Exe
os.system('pyinstaller --onefile --console script.py')
# convert to GUI based Exe
os.system('pyinstaller --onefile script.py')
# convert to Console based Exe with icon
os.system('pyinstaller --onefile --console --icon=icon.ico script.py')
# convert to GUI based Exe with icon
os.system('pyinstaller --onefile --icon=icon.ico script.py')
# convert to Console based Exe with icon and add a file
os.system('pyinstaller --onefile --console --icon=icon.ico --add-data="file.txt;." script.py')
# convert to GUI based Exe with icon and add a folder
os.system('pyinstaller --onefile --icon=icon.ico --add-data="folder;." script.py')