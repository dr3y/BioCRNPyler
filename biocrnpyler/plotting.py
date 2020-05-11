# pathutil.py - path utilities
# ASS, 21 Apr 2020
#
# This file contains some utility functions for plotting CRNs
#
# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import random
import networkx as nx
import statistics
from bokeh.models import (BoxSelectTool, Circle,Square, EdgesAndLinkedNodes, HoverTool,
                          MultiLine, NodesAndLinkedEdges, Plot, Range1d, TapTool,PanTool,WheelZoomTool)
from bokeh.palettes import Spectral4
from bokeh.plotting import from_networkx
from fa2 import ForceAtlas2
import numpy as np
from matplotlib import cm
import dnaplotlib as dpl
import matplotlib.pyplot as plt
from biocrnpyler.dna_construct import DNA_construct

def updateLimits(limits,xvalues):
    for value in xvalues:
        if(value < limits[0]):
            limits[0] = value
        if(value > limits[1]):
            limits[1] = value
    return limits

def makeArrows2(graph_renderer,graph,positions,headsize=3,headangle=np.pi/6):
    """this function draws an arrow shape at the end of graph lines"""
    xs,ys = [],[]
    xbounds = [0,0]
    ybounds = [0,0]
    for edge in graph.edges:
        #iterate through all the edges
        from_node = edge[0]
        to_node = edge[1]
        from_x = positions[edge[0]][0]
        from_y = positions[edge[0]][1]
        to_x = positions[edge[1]][0]
        to_y = positions[edge[1]][1]
        updateLimits(xbounds,[from_x,to_x])
        updateLimits(ybounds,[from_y,to_y])
        #above, we get all the variables
        #below, we are assuming the "to" position is the middle of the
        #coordinate space
        ydif = from_y-to_y
        xdif = from_x-to_x
        #next we calculate the angle from the destination node 
        #to the source node
        angl = np.arctan2(ydif,xdif)
        #the arrow consists of three added points, one on either side 
        #of the line and one in the middle
        p1x = to_x+headsize*np.cos(angl+headangle) #left side of the arrow
        p1y = to_y+headsize*np.sin(angl+headangle) #left side of the arrow
        p2x = to_x+headsize*np.cos(angl-headangle) #right side of the arrow
        p2y = to_y+headsize*np.sin(angl-headangle) #right side of the arrow
        p3x = to_x+headsize*.7*np.cos(angl) #middle of the arrow
        p3y = to_y+headsize*.7*np.sin(angl) #middle of the arrow
        xs.append([from_x,p3x,p1x,to_x,p2x,p3x]) #'xs' is a list of lists which represent each line from node to node
        ys.append([from_y,p3y,p1y,to_y,p2y,p3y]) #'ys' is the same thing except the y positions
    graph_renderer.edge_renderer.data_source.data['xs'] = xs #this part replaces the lines with the ones made by this function
    graph_renderer.edge_renderer.data_source.data['ys'] = ys
    return xbounds,ybounds

def graphPlot(DG,DGspecies,DGreactions,plot,layout="force"):
    """given a directed graph, plot it!"""
    if(layout=="force"):
        #below are parameters for the force directed graph visualization
        forceatlas2 = ForceAtlas2(
                            # Behavior alternatives
                            outboundAttractionDistribution=True,  # Dissuade hubs
                            linLogMode=False,  # NOT IMPLEMENTED
                            adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
                            edgeWeightInfluence=1.0,

                            # Performance
                            jitterTolerance=1.0,  # Tolerance
                            barnesHutOptimize=True,
                            barnesHutTheta=1.2,
                            multiThreaded=False,  # NOT IMPLEMENTED

                            # Tuning
                            scalingRatio=2.0,
                            strongGravityMode=False,
                            gravity=1.0,

                            # Log
                            verbose=False)

        positions = forceatlas2.forceatlas2_networkx_layout(DG, pos=None, iterations=2000) 
    elif(layout == "circle"):
        #here we would arrange the nodes in a circle
        pass
    reaction_renderer = from_networkx(DGreactions, positions, center=(0, 0))
    species_renderer = from_networkx(DGspecies, positions, center=(0, 0))
    edges_renderer = from_networkx(DG, positions, center=(0, 0))

    #edges
    edges_renderer.node_renderer.glyph = Circle(size=12,line_alpha=0,fill_alpha=0, fill_color="color")
    edges_renderer.edge_renderer.glyph = MultiLine( line_alpha=0.2, line_width=4)
    edges_renderer.edge_renderer.selection_glyph = MultiLine(line_color=Spectral4[2], line_width=5)
    edges_renderer.edge_renderer.hover_glyph = MultiLine(line_color=Spectral4[1], line_width=5)
    xbounds,ybounds = makeArrows2(edges_renderer,DG,positions,headsize=5) #make the arrows!
    
    #we want to find the middle of the graph and plot a square that is 1:1 aspect ratio
    
    #find the midpoint of the graph
    xmid = statistics.mean(xbounds)
    ymid = statistics.mean(ybounds)
    #now, subtract the middle from the edges
    xmiddlized = [a-xmid for a in xbounds]
    ymiddlized = [a-ymid for a in ybounds]
    #now, find the biggest dimension
    absdim = max([abs(a) for a in xmiddlized+ymiddlized])
    xlim = [xmid-absdim*1.05,xmid + absdim*1.05]
    ylim = [ymid-absdim*1.05,ymid + absdim*1.05]
    #now set it on the plot!
    plot.x_range = Range1d(xlim[0],xlim[1])
    plot.y_range = Range1d(ylim[0],ylim[1])

    #reactions
    reaction_renderer.node_renderer.glyph = Square(size=8, fill_color=Spectral4[0])
    reaction_renderer.node_renderer.selection_glyph = Square(size=8, fill_color=Spectral4[2])
    reaction_renderer.node_renderer.hover_glyph = Square(size=8, fill_color=Spectral4[1])

    #nodes
    species_renderer.node_renderer.glyph = Circle(size=12, fill_color="color")
    species_renderer.node_renderer.selection_glyph = Circle(size=15, fill_color=Spectral4[2])
    species_renderer.node_renderer.hover_glyph = Circle(size=15, fill_color=Spectral4[1])
    
    #this part adds the interactive elements that make it so that the lines are highlighted 
    #when you mouse over and click
    edge_hover_tool = HoverTool(tooltips= None,renderers=[edges_renderer])
    node_hover_tool = HoverTool(tooltips=[("name", "@species"), ("type", "@type")],\
                                        renderers=[reaction_renderer,species_renderer],attachment="right")
    plot.add_tools(edge_hover_tool,node_hover_tool, TapTool(), BoxSelectTool(),PanTool(),WheelZoomTool())

    edges_renderer.selection_policy = NodesAndLinkedEdges()
    edges_renderer.inspection_policy = EdgesAndLinkedNodes()

    plot.renderers.append(edges_renderer)
    plot.renderers.append(reaction_renderer)
    plot.renderers.append(species_renderer)

def generate_networkx_graph(CRN,useweights=False):
    """generates a networkx DiGraph object that represents the CRN."""
    CRNgraph = nx.DiGraph()
    allnodenum = 1 #every node has an index
    nodedict = {} #this is so that we can write out the reactions in
                #the reaction "species" field
                #it has {species:index}
    rxnlist = [] #list of numbers corresponding to only reaction nodes
    speclist = CRN.species
    nodedict["nothing"]=0
    CRNgraph.add_node(0)
    CRNgraph.nodes[0]["type"]="nothing"
    CRNgraph.nodes[0]["species"]="nothing"
    CRNgraph.nodes[0]["color"]="purple"
    for specie in CRN.species:
        mycol = "teal"
        if(specie.material_type=="complex"):
            mycol = "cyan"
        elif(specie.material_type=="protein"):
            mycol = "green"
        elif(specie.material_type=="dna"):
            mycol = "grey"
        elif(specie.material_type=="rna"):
            mycol = "orange"
        elif(specie.material_type=="phosphate"):
            mycol = "yellow"
        elif(specie.material_type=="ligand"):
            mycol = "pink"
        nodedict[specie]=allnodenum
        CRNgraph.add_node(allnodenum)
        CRNgraph.nodes[allnodenum]["type"]=str(specie.material_type)
        CRNgraph.nodes[allnodenum]["species"]=str(specie)
        CRNgraph.nodes[allnodenum]["color"]=mycol
        allnodenum +=1
    for rxn in CRN.reactions:
        CRNgraph.add_node(allnodenum)
        CRNgraph.nodes[allnodenum]["type"]=rxn.propensity_type
        mycol = "blue"
        #CRNgraph.nodes[allnodenum]
        kval = rxn.k
        if(not useweights):
            kval = 1
        krev_val = rxn.k_r
        if((krev_val > 0) and (not useweights)):
            krev_val = 1
        for reactant in rxn.inputs:
            CRNgraph.add_edge(nodedict[reactant],allnodenum,weight=kval)
            if(krev_val>0):
                CRNgraph.add_edge(allnodenum,nodedict[reactant],weight=krev_val)
        for product in rxn.outputs:
            CRNgraph.add_edge(allnodenum,nodedict[product],weight=kval)
            if(krev_val>0):
                CRNgraph.add_edge(nodedict[product],allnodenum,weight=krev_val)
        if(len(rxn.outputs)==0):
            CRNgraph.add_edge(allnodenum,0,weight=kval)
            if(krev_val>0):
                CRNgraph.add_edge(0,allnodenum,weight=krev_val)
        elif(len(rxn.inputs)==0):
            CRNgraph.add_edge(0,allnodenum,weight=kval)
            if(krev_val>0):
                CRNgraph.add_edge(allnodenum,0,weight=krev_val)
        CRNgraph.nodes[allnodenum]["color"]=mycol
        CRNgraph.nodes[allnodenum]["species"]=str(rxn)
        rxnlist += [allnodenum]
        allnodenum +=1
    CRNspeciesonly = CRNgraph.copy()
    CRNspeciesonly.remove_nodes_from(rxnlist)
    CRNreactionsonly = CRNgraph.copy()
    CRNreactionsonly.remove_nodes_from(range(rxnlist[0]))
    return CRNgraph,CRNspeciesonly,CRNreactionsonly

def make_dpl_from_construct(construct,showlabels=[]):
    outdesign = []
    cmap = cm.Set1(range(len(construct.parts_list)))
    pind = 0
    for part in construct.parts_list:
        showlabel = False
        if(part.part_type in showlabels):
            showlabel = True
        outdesign+=make_dpl_from_part(part,color=cmap[pind][:-1],color2 = random.choice(cmap)[:-1],showlabel=showlabel)
        pind+=1
    return outdesign
def make_dpl_from_part(part,direction=None,color=(1,4,2),color2=(3,2,4),showlabel=False):
    if(direction==None and part.direction != None):
        direction = part.direction=="forward"
    elif(direction==None):
        direction = True
    if(not type(part.color)==type(None)):
        color = part.color
    if(not type(part.color2)==type(None)):
        color2 = part.color2
    dnaplotlib_dict = {\
        "promoter":"Promoter",\
        "rbs":"RBS",\
        "CDS":"CDS",\
        "terminator":"Terminator",\
        "attP":"RecombinaseSite",\
        "attB":"RecombinaseSite",\
        "attL":"RecombinaseSite2",\
        "attR":"RecombinaseSite2"}
    dpl_type = dnaplotlib_dict[part.part_type]
    outdesign = [{'type':dpl_type,"name":part.name,"fwd":direction,'opts':{'color':color,'color2':color2}}]
    if(part.regulator!= None):
        outdesign += [{"type":"Operator","name":part.regulator,"fwd":direction,'opts':{'color':color,'color2':color2}}]
    if(showlabel):
        outdesign[0]["opts"].update({'label':str(part),'label_size':13,'label_y_offset':-8,})
    if(not direction):
        outdesign = outdesign[::-1]
    return outdesign

def plotConstruct(DNA_construct_obj,dna_renderer=dpl.DNARenderer(scale = 5,linewidth=3),\
                                    rna_renderer=dpl.DNARenderer(scale = 5,linewidth=3,linecolor=(1,0,0)),\
                                    plot_rnas=False,debug=False):
    """helper function for making dnaplotlib plots of a DNA_construct object. Plots the
    DNAs and the RNAs that come from that DNA, using DNA_construct.explore_txtl"""
    design = make_dpl_from_construct(DNA_construct_obj,showlabels=["attB","attP","attL","attR"])
    circular=DNA_construct_obj.circular
    plotDesign(design,circular=circular)
    if(plot_rnas):
        rnas,proteins = DNA_construct_obj.explore_txtl()
        if(debug):
            print("rna:")
            print(rnas)
            print("protein")
            print(proteins)
        for promoter in rnas:
            rnadesign = make_dpl_from_construct(DNA_construct(rnas[promoter]))
            rnacolor = rna_renderer.linecolor
            for part in rnadesign:
                if("edgecolor" not in part['opts']):
                    part['opts'].update({'edgecolor':rnacolor})
            plotDesign(rnadesign,renderer=rna_renderer)

def plotDesign(design,renderer = dpl.DNARenderer(scale = 5,linewidth=3),part_renderers=None,\
                circular=False):
    """helper function for doing dnaplotlib plots. You need to set the size and min max of the
    plot, and that's what this function does"""
    if(part_renderers==None):
        part_renderers = renderer.SBOL_part_renderers()
    fig = plt.figure(figsize=(len(design)*.75,1.1))
    ax = fig.add_axes([0,0,1,1])
    start,end = renderer.renderDNA(ax,design,part_renderers,circular=circular)
    ax.axis('off')
    addedsize=1
    ax.set_xlim([start-addedsize,end+addedsize])
    ax.set_ylim([-15,15])
    plt.show()