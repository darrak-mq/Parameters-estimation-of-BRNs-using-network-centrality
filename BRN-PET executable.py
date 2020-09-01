# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 11:27:08 2019

@author: dmoin.msbi17rcms
"""

from tkinter import *
import tkinter as tk
import networkx as nx
import math
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import itertools
from tqdm import tqdm
from pyrthomas.network_analyser import NetworkAnalyser
import csv
import pandas as pd
import seaborn as sns

root = tk.Tk()
root.title("BRN-PET")
#using geometry attribute to change the root windows size
root.geometry("1000x680") #size of the app to be widthx500
root.resizable(0, 0) #Doesn't allow resizing in the x or y direction
root.configure(background='maroon')

G=nx.DiGraph()
pNodes = []
comb_dic = {}
dic_lst = []

############################## go1 button ###########################
def nCr(n,r):
    f = math.factorial
    return f(n) / (f(r) * f(n-r))
   
def go1_activity():
    fileName = nameVar.get()
    file = open(fileName, 'r')
    returnSepLines = file.read().split('\n') #separating lines on the basis of return key
    List_tabSepLines = []
    for num in range(len(returnSepLines)): #this for loop is adding edges between nodes
        tabSepLines = returnSepLines[num].split('\t')
        G.add_edge(tabSepLines[0],tabSepLines[2], weight=int(tabSepLines[1]))
        List_tabSepLines.append(tabSepLines)   
    f = Figure(figsize=(7,3), dpi=130)
    a = f.add_subplot(111)
    
    pos=nx.spring_layout(G)
    nx.draw(G, pos, ax = a, with_labels = True)
    
    labels = nx.get_edge_attributes(G,'weight')
    nx.draw_networkx_edge_labels(G,pos,edge_labels=labels, ax=a)
    
    f11 = tk.Frame(master = root, bg="maroon",width=1000,height=390)
    f11.grid(row=10,column=0)
    f11.grid_propagate(0)
    f11.update()
    canvas = FigureCanvasTkAgg(f, f11)
    canvas.draw()
    canvas.get_tk_widget().place(x=40, y=1000)
    canvas._tkcanvas.place(x=45, y=0)
    
    f12 = tk.Frame(master = root, bg="maroon",width=1000,height=40)
    f12.grid(row=11,column=0)
    f12.grid_propagate(0)
    f12.update()
    toolbar = NavigationToolbar2TkAgg(canvas, f12)
    toolbar.place(x=713, y=0)
    toolbar.update()

    pNodes = list(G.nodes()).copy()
    source_list = [row[0] for row in List_tabSepLines].copy()
    intr_list = [row[1] for row in List_tabSepLines].copy()
    sink_list = [row[2] for row in List_tabSepLines].copy()
    nval = 1
    val = 1
    comb_dic = {}
    for node_index, node_val in enumerate(pNodes):
        list_of_combs = []
        nval = nval * val
        comb = 0
        max_interaction_val = 0
        min_interaction_val = 0
        num_of_source_nodes = 0
        node_source_list = []
        range_list = []

        for sink_index, sink_val in enumerate(sink_list):
            if node_val == sink_val:
                num_of_source_nodes = num_of_source_nodes + 1
                for source_index, source_val in enumerate(source_list):
                    if source_index == sink_index:
                        node_source_list.append(source_val)
                for source_index, source_val in enumerate(source_list):   
                    if node_val == source_val:
                        for intr_index, intr_val in enumerate(intr_list):
                            if source_index == intr_index: #syncing index numbers of both sink and interaction lists
                                intr_val = int(intr_val)
                                if intr_val > max_interaction_val: #calculating maximum interaction value of node 
                                    max_interaction_val = intr_val
                                if intr_val <= min_interaction_val: #calculating minimum interaction value of node 
                                    min_interaction_val = intr_val
        list_of_combs = []
        if num_of_source_nodes == 0:
            list_of_combs = [[]]
        if num_of_source_nodes != 0:
            for r in range(num_of_source_nodes + 1):
                com = nCr(num_of_source_nodes, r)
                comb = comb + com
            print ("printing list of combinations")
            
            for i in range(len(node_source_list)+1):
                r_comb_list = [list(cnod) for cnod in itertools.combinations(node_source_list, i)]
                list_of_combs.extend(r_comb_list)
        comb_dic[node_val] = list_of_combs

        print ("printing list of combinations")

        for rnum in range((max_interaction_val - min_interaction_val) + 1):
            if min_interaction_val <= max_interaction_val:
                range_list.append(min_interaction_val)
                min_interaction_val = min_interaction_val + 1
    
        temp = len(range_list) 
        val = temp**comb # all possible combinations values for one node respectively             
    nval = nval * val
    file.close()

    state_graphs = NetworkAnalyser.get_possible_state_graphs(G)
    for graph in tqdm(state_graphs, total=int(nval)):
        bw_centrality = nx.betweenness_centrality(graph, normalized=False) #calculating betweenness centrality of every node
        dic_lst.append(bw_centrality)
  
    stateEntry = tk.Entry(f13, textvariable = stateVar)
    stateEntry.config(state = tk.NORMAL)
    stateEntry.place(x=115, y=6)

################################## Central button #################################

def maxCentrality():
    root1 = tk.Tk()
    root1.title("Parameters Heatmap")
    root1.geometry("1000x700") #size of the app to be widthx500
    root1.resizable(0, 0) #Doesn't allow resizing in the x or y direction
    root1.configure(background='white')
#### inputting the desired nodes of state graph for heatmap generation##########
    nodes = stateVar.get()
    bet_cen_list = [] #creating empty list to store betweenness centrality values in all state graphs of the entered node 
    for i in dic_lst: #i is the dictionary
        bet_cen_list.append(i[nodes]) #storing betweenness cen values of entered node w.r.t each state graph 
    sorted_betCen_list = sorted(([v, ind] for ind, v in enumerate(bet_cen_list)), reverse=True) #sort values(descending) meanwhile keep their original index num
    del bet_cen_list[:] 
    del sorted_betCen_list[10:]
    csvData1 = [['Betweenness Centrality', 'State Graphs']] #creating a csv file containing only the headers
    with open ('CenData.csv', 'w', newline='') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(csvData1)    
    csvFile.close()
    
    with open ('CenData.csv', 'a', newline='') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(sorted_betCen_list)
    csvFile.close()
    
    sg_list = []
    with open('CenData.csv', 'r') as csvFile:
        reader = csv.reader(csvFile)
        for row in reader:
            sg_list.append(row[1])
    
    csvFile.close()
    del sg_list[0:1]
    params = NetworkAnalyser.get_possible_parameters(G)
    
    csvData2 = [['State Graphs', 'Resources', 'Parameters']] #creating a csv file containing only the headers
    with open ('heatmapData.csv', 'w', newline='') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(csvData2)    
    csvFile.close()
    hdata = []
    fileName = nameVar.get()
    file = open(fileName, 'r')
    returnSepLines = file.read().split('\n') #separating lines on the basis of return key
    List_tabSepLines = []
    for num in range(len(returnSepLines)): #this for loop is adding edges between nodes
        tabSepLines = returnSepLines[num].split('\t')
        List_tabSepLines.append(tabSepLines)   
    
    pNodes = list(G.nodes()).copy()
    source_list = [row[0] for row in List_tabSepLines].copy()
    intr_list = [row[1] for row in List_tabSepLines].copy()
    sink_list = [row[2] for row in List_tabSepLines].copy()
    nval = 1
    val = 1
    comb_dic = {}
    for node_index, node_val in enumerate(pNodes):
        list_of_combs = []
        nval = nval * val
        comb = 0
        max_interaction_val = 0
        min_interaction_val = 0
        num_of_source_nodes = 0
        node_source_list = []
        for sink_index, sink_val in enumerate(sink_list):
            if node_val == sink_val:
                num_of_source_nodes = num_of_source_nodes + 1
                for source_index, source_val in enumerate(source_list):
                    if source_index == sink_index:
                        node_source_list.append(source_val)
                for source_index, source_val in enumerate(source_list):   
                    if node_val == source_val:
                        for intr_index, intr_val in enumerate(intr_list):
                            if source_index == intr_index: #syncing index numbers of both sink and interaction lists
                                intr_val = int(intr_val)
                                if intr_val > max_interaction_val: #calculating maximum interaction value of node 
                                    max_interaction_val = intr_val
                                if intr_val <= min_interaction_val: #calculating minimum interaction value of node 
                                    min_interaction_val = intr_val
        if num_of_source_nodes != 0:
            for r in range(num_of_source_nodes + 1):
                com = nCr(num_of_source_nodes, r)
                comb = comb + com
            list_of_combs = []
            for i in range(len(node_source_list)+1):
                r_comb_list = [list(cnod) for cnod in itertools.combinations(node_source_list, i)]
                list_of_combs.extend(r_comb_list)
            comb_dic[node_val] = list_of_combs

    for sg_ind, sg_val in enumerate(sg_list):
        j = -1
        pval_list = []
        sg_val = int(sg_val)
        par_dic = params[sg_val]
        for item in par_dic:
            for tup in par_dic[item]:
                pval_list.append(tup[1])
        for nod_ind, nod_val in enumerate(pNodes):
            for dit in comb_dic:
                if nod_val == dit:
                    clist = comb_dic[dit]
                    for cind, cval in enumerate(clist):
                        j = j + 1
                        res_string = 'k' + str(nod_val) + str(cval)
                        hdata.append([sg_val, res_string, pval_list[j]])
    with open ('heatmapData.csv', 'a', newline='') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(hdata)
    csvFile.close()

    def create_plot():
        helix = pd.read_csv('C:/Users/darrak/.spyder/heatmapData.csv')
        f, ax = plt.subplots(figsize=(11, 9))
        pivot_table = helix.pivot('Resources', 'State Graphs', 'Parameters')
        plt.xlabel('State Graphs', size = 15)
        plt.ylabel('Resources', size = 15)
        plt.title('Parameters Profiling', size = 15)
        sns.heatmap(pivot_table, annot=True, linewidths=.003, square = True, cmap = 'BuPu');
        return f

    fig = create_plot()
    nodes =str(nodes)
    sgStateString = "Heatmap for " + nodes + " (Highly Central)"
    f15 = tk.Frame(master = root1, bg="white", width=1000, height=40)
    f15.grid(row=0,column=0)
    f15.grid_propagate(0)
    f15.update()
    sgStateLabel = tk.Label(f15, text = sgStateString, font=("Times New Roman", 20, "bold"), bg = "white", fg = "maroon")
    sgStateLabel.place(x=365, y=5) 
    
    f16 = tk.Frame(master = root1, bg="white", width=1000, height=620)
    f16.grid(row=1,column=0)
    f16.grid_propagate(0)
    f16.update()
    canvas = FigureCanvasTkAgg(fig, f16)
    canvas.draw()
    canvas.get_tk_widget().place(x=40, y=1000)
    canvas._tkcanvas.place(x=100, y=0)
    
    f17 = tk.Frame(master = root1, bg="white", width=1000, height=40)
    f17.grid(row=2,column=0)
    f17.grid_propagate(0)
    f17.update()
    toolbar = NavigationToolbar2TkAgg(canvas, f17)
    toolbar.place(x=0, y=0)
    toolbar.update()

    root1.mainloop()

################## DEADLOCK BUTTON #############################################
def minCentrality():
    root2 = tk.Tk()
    root2.title("Parameters Heatmap")
    root2.geometry("1000x700") #size of the app to be widthx500
    root2.resizable(0, 0) #Doesn't allow resizing in the x or y direction
    root2.configure(background='white')
#### inputting the desired nodes of state graph for heatmap generation##########
    nodes = stateVar.get()
    bet_cen_list = [] #creating empty list to store betweenness centrality values in all state graphs of the entered node
    for i in dic_lst: #i is the dictionary
        bet_cen_list.append(i[nodes]) #storing betweenness cen values of entered node w.r.t each state graph 
    sorted_betCen_list = sorted(([v, ind] for ind, v in enumerate(bet_cen_list))) #sort values(descending) meanwhile keep their original index num
    del bet_cen_list[:] 
    del sorted_betCen_list[10:]
    csvData1 = [['Betweenness Centrality', 'State Graphs']] #creating a csv file containing only the headers
    with open ('CenData.csv', 'w', newline='') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(csvData1)    
    csvFile.close()
    
    with open ('CenData.csv', 'a', newline='') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(sorted_betCen_list)
    csvFile.close()
    
    sg_list = []
    with open('CenData.csv', 'r') as csvFile:
        reader = csv.reader(csvFile)
        for row in reader:
            sg_list.append(row[1])
    
    csvFile.close()
    del sg_list[0:1]
    params = NetworkAnalyser.get_possible_parameters(G)
    
    csvData2 = [['State Graphs', 'Resources', 'Parameters']] #creating a csv file containing only the headers
    with open ('heatmapData.csv', 'w', newline='') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(csvData2)    
    csvFile.close()
    hdata = []

    fileName = nameVar.get()
    file = open(fileName, 'r')
    returnSepLines = file.read().split('\n') #separating lines on the basis of return key
    List_tabSepLines = []
    for num in range(len(returnSepLines)): #this for loop is adding edges between nodes
        tabSepLines = returnSepLines[num].split('\t')
        List_tabSepLines.append(tabSepLines)   
    
    pNodes = list(G.nodes()).copy()
    source_list = [row[0] for row in List_tabSepLines].copy()
    intr_list = [row[1] for row in List_tabSepLines].copy()
    sink_list = [row[2] for row in List_tabSepLines].copy()
    nval = 1
    val = 1
    comb_dic = {}
    for node_index, node_val in enumerate(pNodes):
        list_of_combs = []
        nval = nval * val
        comb = 0
        max_interaction_val = 0
        min_interaction_val = 0
        num_of_source_nodes = 0
        node_source_list = []
        for sink_index, sink_val in enumerate(sink_list):
            if node_val == sink_val:
                num_of_source_nodes = num_of_source_nodes + 1
                for source_index, source_val in enumerate(source_list):
                    if source_index == sink_index:
                        node_source_list.append(source_val)
                for source_index, source_val in enumerate(source_list):   
                    if node_val == source_val:
                        for intr_index, intr_val in enumerate(intr_list):
                            if source_index == intr_index: #syncing index numbers of both sink and interaction lists
                                intr_val = int(intr_val)
                                if intr_val > max_interaction_val: #calculating maximum interaction value of node 
                                    max_interaction_val = intr_val
                                if intr_val <= min_interaction_val: #calculating minimum interaction value of node 
                                    min_interaction_val = intr_val
        if num_of_source_nodes != 0:
            for r in range(num_of_source_nodes + 1):
                com = nCr(num_of_source_nodes, r)
                comb = comb + com
            list_of_combs = []
            for i in range(len(node_source_list)+1):
                r_comb_list = [list(cnod) for cnod in itertools.combinations(node_source_list, i)]
                list_of_combs.extend(r_comb_list)
            comb_dic[node_val] = list_of_combs

    for sg_ind, sg_val in enumerate(sg_list):
        j = -1
        pval_list = []
        sg_val = int(sg_val)
        par_dic = params[sg_val]
        for item in par_dic:
            for tup in par_dic[item]:
                pval_list.append(tup[1])
        for nod_ind, nod_val in enumerate(pNodes):
            for dit in comb_dic:
                if nod_val == dit:
                    clist = comb_dic[dit]
                    for cind, cval in enumerate(clist):
                        j = j + 1
                        res_string = 'k' + str(nod_val) + str(cval)
                        hdata.append([sg_val, res_string, pval_list[j]])
    with open ('heatmapData.csv', 'a', newline='') as csvFile:
        writer = csv.writer(csvFile)
        writer.writerows(hdata)
    csvFile.close()

    def create_plot():
        helix = pd.read_csv('C:/Users/darrak/.spyder/heatmapData.csv')
        f, ax = plt.subplots(figsize=(11, 9))
        pivot_table = helix.pivot('Resources', 'State Graphs', 'Parameters')
        plt.xlabel('State Graphs', size = 15)
        plt.ylabel('Resources', size = 15)
        plt.title('Parameters Profiling', size = 15)
        sns.heatmap(pivot_table, annot=True, linewidths=.003, square = True, cmap = 'BuPu');
        return f
    fig = create_plot()
    nodes =str(nodes)
    sgStateString = "Heatmap for " + nodes + " (Least Central)"
    f15 = tk.Frame(master = root2, bg="white", width=1000, height=40)
    f15.grid(row=0,column=0)
    f15.grid_propagate(0)
    f15.update()
    sgStateLabel = tk.Label(f15, text = sgStateString, font=("Times New Roman", 20, "bold"), bg = "white", fg = "maroon")
    sgStateLabel.place(x=365, y=5) 
    
    f16 = tk.Frame(master = root2, bg="white", width=1000, height=620)
    f16.grid(row=1,column=0)
    f16.grid_propagate(0)
    f16.update()
    canvas = FigureCanvasTkAgg(fig, f16)
    canvas.draw()
    canvas.get_tk_widget().place(x=40, y=1000)
    canvas._tkcanvas.place(x=100, y=0)
    
    f17 = tk.Frame(master = root2, bg="white", width=1000, height=40)
    f17.grid(row=2,column=0)
    f17.grid_propagate(0)
    f17.update()
    toolbar = NavigationToolbar2TkAgg(canvas, f17)
    toolbar.place(x=0, y=0)
    toolbar.update()
    
    root2.mainloop()
    
f1 = tk.Frame(master = root, bg="maroon",width=1000,height=40)
f1.grid(row=0,column=0)
f1.grid_propagate(0)
f1.update()
Label1 = tk.Label(f1, text="Parameters Estimation", font=("Times New Roman", 20, "bold"), bg = "maroon", fg = "white")
Label1.place(x=365, y=5)

f2 = tk.Frame(master = root, bg="maroon",width=1000,height=25)
f2.grid(row=1,column=0)
f2.grid_propagate(0)
f2.update()
Label2 = tk.Label(f2,text="for",font=("Times New Roman", 14, "bold"), bg="maroon", fg = "white")
Label2.place(x=485, y=0)

f3 = tk.Frame(master = root, bg="maroon",width=1000,height=40)
f3.grid(row=2,column=0)
f3.grid_propagate(0)
f3.update()
Label3 = tk.Label(f3, text="Biological Regulatory Networks", font=("Times New Roman", 20, "bold"), bg = "maroon", fg = "white")
Label3.place(x=312, y=3)

f4 = tk.Frame(master = root, bg="maroon",width=1000,height=15)
f4.grid(row=3,column=0)
f4.grid_propagate(0)
f4.update()

f5 = tk.Frame(master = root, bg="white",width=1000,height=2)
f5.grid(row=4,column=0)
f5.grid_propagate(0)
f5.update()

f6 = tk.Frame(master = root, bg="maroon",width=1000,height=4)
f6.grid(row=5,column=0)
f6.grid_propagate(0)
f6.update()

f7 = tk.Frame(master = root, bg="white",width=1000,height=2)
f7.grid(row=6,column=0)
f7.grid_propagate(0)
f7.update()

f8 = tk.Frame(master = root, bg="maroon",width=1000,height=20)
f8.grid(row=7,column=0)
f8.grid_propagate(0)
f8.update()

nameVar = tk.StringVar()
f9 = tk.Frame(master = root, bg="maroon",width=1000,height=40)
f9.grid(row=8,column=0)
f9.grid_propagate(0)
f9.update()
nameLabel = tk.Label(f9, text = "File Name",font=("Times New Roman", 12, "bold"), bg = "maroon", fg = "white")
nameLabel.place(x=0, y=3)
nameEntry = tk.Entry(f9, textvariable = nameVar)
nameEntry.place(x=80, y=6)
nameButton = tk.Button(f9, text = "GO", command = go1_activity)
nameButton.place(x=210, y=3)
sampleLabel = tk.Label(f9, text = "sample: filename.sif",font=("Times New Roman", 10, "bold"), bg = "maroon", fg = "white")
sampleLabel.place(x=240, y=3)

f10 = tk.Frame(master = root, bg="maroon",width=1000,height=10)
f10.grid(row=9,column=0)
f10.grid_propagate(0)
f10.update()

f11 = tk.Frame(master = root, bg="maroon",width=1000,height=390)
f11.grid(row=10,column=0)
f11.grid_propagate(0)
f11.update()
    
f12 = tk.Frame(master = root, bg="maroon",width=1000,height=40)
f12.grid(row=11,column=0)
f12.grid_propagate(0)
f12.update()

stateVar = tk.StringVar()
f13 = tk.Frame(master = root, bg="maroon",width=1000,height=32)
f13.grid(row=12,column=0)
f13.grid_propagate(0)
f13.update()
stateLabel = tk.Label(f13, text = "State of interest",font=("Times New Roman", 12, "bold"), bg = "maroon", fg = "white")
stateLabel.place(x=0, y=3)
stateEntry = tk.Entry(f13, textvariable = stateVar)
stateEntry.config(state = tk.DISABLED)
stateEntry.place(x=115, y=6)
stateButton = tk.Button(f13, text = "Highly Central", font=("Times New Roman", 10, "bold"), command = maxCentrality)
stateButton.place(x=245, y=3)
deadlockButton = tk.Button(f13, text = "Least Central", font=("Times New Roman", 10, "bold"), command = minCentrality)
deadlockButton.place(x=343, y=3)
sampleLabel2 = tk.Label(f13, text = "sample: 0,1,1 (as per your number of entities in BRN)",font=("Times New Roman", 10, "bold"), bg = "maroon", fg = "white")
sampleLabel2.place(x=432, y=4)

root.mainloop()
