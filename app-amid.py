#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 13 13:47:18 2021

@author: mmemc
"""

from amid import AMID
import streamlit as st
import tkinter as tk
from tkinter import filedialog
import bokeh as bk
import numpy as np
from pathlib import Path


st.write("""
         # Welcome to the AMID!
         """)

fdata = st.file_uploader("Load your UHPC file here!")

if fdata is not None:
    i = 0
    for line in fdata:
        if i == 0:
            st.write("{}".format(line.decode()))
        i = i +1

# Set up tkinter
root = tk.Tk()
root.withdraw()

# Make folder picker dialog appear on top of other windows
root.wm_attributes('-topmost', 1)

# Folder picker button
st.title('Set output folder.')

st.write('Where would you like graphs and datafiles saved?')
clicked = st.button('Select DST')
if clicked:
    dstdir = st.text_input('Selected DST folder:', filedialog.askdirectory(master=root))
    
    