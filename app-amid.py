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


@st.cache(persist=True)
def read_data(bytesIO, dstdir, cell_label):
    return AMID(dstdir, "", "", cell_label, bytesIO=bytesIO)

@st.cache(persist=True, suppress_st_warning=True)
def output_dir(clicked):
    if clicked is True:
        dstdir = st.text_input('Selected DST folder:', filedialog.askdirectory(master=root))
        return dstdir
    else:
        return "./"

st.write("""
         # Welcome to the AMID!
         """)

bytesIO = st.file_uploader("Load your UHPC file here!")

# Set up tkinter
root = tk.Tk()
root.withdraw()

# Make folder picker dialog appear on top of other windows
root.wm_attributes('-topmost', 1)

# Folder picker button
st.title('Set output folder.')

st.write('Where would you like graphs and datafiles saved?')
clicked = st.button('Select DST')
dstdir = output_dir(clicked)
#dstdir = './'
#if clicked:
#    dstdir = st.text_input('Selected DST folder:', filedialog.askdirectory(master=root))
    
cell_label = st.text_input('Cell Label')

plot_protocol = st.sidebar.button("Plot Protocol")
if plot_protocol:
    save_protocol = st.sidebar.button("Save Protocol")

if bytesIO is not None:
    amid_data = read_data(bytesIO, dstdir, cell_label)
    #amid_data = AMID(dstdir, "", "", cell_label, bytesIO=fdata)
    if plot_protocol:
        pro_fig = amid_data.plot_protocol(save=False, return_fig=True)
        st.pyplot(fig=pro_fig)
    
    