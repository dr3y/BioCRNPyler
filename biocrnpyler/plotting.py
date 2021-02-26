# plotting.py - plotting utilities
# ASS, 21 Apr 2020
#
# This file contains some utility functions for plotting CRNs
#
# Copyright (c) 2018, Build-A-Cell. All rights reserved.
# See LICENSE file in the project root directory for details.

import math
import random
import statistics
from warnings import warn

from .components_basic import Protein, DNA, RNA
from .dna_part_cds import CDS
from .dna_part_misc import IntegraseSite
from .dna_part_promoter import Promoter
from .dna_part_rbs import RBS
from .dna_part_terminator import Terminator
from .propensities import MassAction
from .species import ComplexSpecies, Species
from .polymer import OrderedPolymer
from .dna_construct import Construct
import io
import base64

import copy
HAVE_MATPLOTLIB = False
try:
    import matplotlib.pyplot as plt
    from matplotlib import cm

    HAVE_MATPLOTLIB = True
except ModuleNotFoundError:
    pass
PLOT_DNA = False
try:
    import dnaplotlib as dpl
    PLOT_DNA = True
except ModuleNotFoundError:
    pass

if(PLOT_DNA and not HAVE_MATPLOTLIB):
    PLOT_DNA = False

PLOT_NETWORK = False
try:
    import networkx as nx
    from bokeh.models import (BoxSelectTool, Circle, EdgesAndLinkedNodes,
                              HoverTool, MultiLine, NodesAndLinkedEdges,
                              PanTool, Plot, Range1d, Square, TapTool,
                              WheelZoomTool)
    from bokeh.plotting import from_networkx
    from bokeh.palettes import Spectral4
    from fa2 import ForceAtlas2
    PLOT_NETWORK = True
except ModuleNotFoundError:
    pass


def updateLimits(limits, xvalues):
    for value in xvalues:
        if(value < limits[0]):
            limits[0] = value
        if(value > limits[1]):
            limits[1] = value
    return limits


def makeArrows2(graph_renderer, graph, positions, headsize=3, headangle=math.pi/6):
    """this function draws an arrow shape at the end of graph lines"""
    xs, ys = [], []
    xbounds = [0, 0]
    ybounds = [0, 0]
    for edge in graph.edges:
        # iterate through all the edges
        from_node = edge[0]
        to_node = edge[1]
        from_x = positions[from_node][0]
        from_y = positions[from_node][1]
        to_x = positions[to_node][0]
        to_y = positions[to_node][1]
        updateLimits(xbounds, [from_x, to_x])
        updateLimits(ybounds, [from_y, to_y])
        # above, we get all the variables
        # below, we are assuming the "to" position is the middle of the
        # coordinate space
        ydif = from_y-to_y
        xdif = from_x-to_x
        # next we calculate the angle from the destination node
        # to the source node
        angl = math.atan2(ydif, xdif)
        # the arrow consists of three added points, one on either side
        # of the line and one in the middle
        p1x = to_x+headsize*math.cos(angl+headangle)  # left side of the arrow
        p1y = to_y+headsize*math.sin(angl+headangle)  # left side of the arrow
        p2x = to_x+headsize*math.cos(angl-headangle)  # right side of the arrow
        p2y = to_y+headsize*math.sin(angl-headangle)  # right side of the arrow
        p3x = to_x+headsize*.7*math.cos(angl)  # middle of the arrow
        p3y = to_y+headsize*.7*math.sin(angl)  # middle of the arrow
        # 'xs' is a list of lists which represent each line from node to node
        xs.append([from_x, p3x, p1x, to_x, p2x, p3x])
        # 'ys' is the same thing except the y positions
        ys.append([from_y, p3y, p1y, to_y, p2y, p3y])
    # this part replaces the lines with the ones made by this function
    graph_renderer.edge_renderer.data_source.data['xs'] = xs
    graph_renderer.edge_renderer.data_source.data['ys'] = ys
    return xbounds, ybounds


def graphPlot(DG,DGspecies,DGreactions,plot,layout="force",positions=None,\
                        posscale = 1.0,layoutfunc=None,iterations=2000,rseed=30,show_species_images=False):
    """given a directed graph, plot it!
    Inputs:
    DG: a directed graph of type DiGraph
    DGspecies: a directed graph which only contains the species nodes
    DGreactions: a directed graph which only contains the reaction nodes
    plot: a bokeh plot object
    layout: graph layout function. 
                'force' uses fa2 to push nodes apart
                'circle' plots the nodes and reactions in two overlapping circles, with the reactions on the inside of the circle
                'custom' allows user input "layoutfunc". Internally, layoutfunc is passed the three inputs (DG, DGspecies, DGreactions)
                                                        and should output a position dictionary with node {<node number>:(x,y)}

    positions: a dictionary of node names and x,y positions. this gets passed into the layout function
    posscale: multiply the scaling of the plot. This only affects the arrows because the arrows are a hack :("""
    random.seed(rseed)
    if(not PLOT_NETWORK):
        warn("network plotting disabled because some libraries are not found")
        return
    if(layout == "force"):
        # below are parameters for the force directed graph visualization
        forceatlas2 = ForceAtlas2(
            # Behavior alternatives
            outboundAttractionDistribution=True,  # Dissuade hubs
            linLogMode=False,  # NOT IMPLEMENTED
            # Prevent overlap (NOT IMPLEMENTED)
            adjustSizes=False,
            edgeWeightInfluence=1.0,

            # Performance
            jitterTolerance=1.0,  # Tolerance
            barnesHutOptimize=True,
            barnesHutTheta=1.2,
            multiThreaded=False,  # NOT IMPLEMENTED

            # Tuning
            scalingRatio=2.4*posscale,
            strongGravityMode=False,
            gravity=1.0,

            # Log
            verbose=False)

        positions = forceatlas2.forceatlas2_networkx_layout(
            DG, pos=positions, iterations=iterations)
    elif(layout == "circle"):
        positions = nx.circular_layout(DGspecies, scale=50*posscale)
        positions.update(nx.circular_layout(DGreactions, scale=35*posscale))
    elif(layout == "custom"):
        positions = layoutfunc(DG, DGspecies, DGreactions)
    reaction_renderer = from_networkx(DGreactions, positions, center=(0, 0))
    species_renderer = from_networkx(DGspecies, positions, center=(0, 0))
    edges_renderer = from_networkx(DG, positions, center=(0, 0))

    # edges
    edges_renderer.node_renderer.glyph = Circle(
        size=12, line_alpha=0, fill_alpha=0, fill_color="color")
    edges_renderer.edge_renderer.glyph = MultiLine(
        line_alpha=0.2, line_width=4, line_join="round")
    edges_renderer.edge_renderer.selection_glyph = MultiLine(
        line_color=Spectral4[2], line_width=5, line_join="round")
    edges_renderer.edge_renderer.hover_glyph = MultiLine(
        line_color=Spectral4[1], line_width=5, line_join="round")
    xbounds, ybounds = makeArrows2(
        edges_renderer, DG, positions, headsize=5)  # make the arrows!

    # we want to find the middle of the graph and plot a square that is 1:1 aspect ratio

    # find the midpoint of the graph
    xmid = statistics.mean(xbounds)
    ymid = statistics.mean(ybounds)
    # now, subtract the middle from the edges
    xmiddlized = [a-xmid for a in xbounds]
    ymiddlized = [a-ymid for a in ybounds]
    # now, find the biggest dimension
    absdim = max([abs(a) for a in xmiddlized+ymiddlized])
    xlim = [xmid-absdim*1.05, xmid + absdim*1.05]
    ylim = [ymid-absdim*1.05, ymid + absdim*1.05]
    # now set it on the plot!
    plot.x_range = Range1d(xlim[0], xlim[1])
    plot.y_range = Range1d(ylim[0], ylim[1])

    # reactions
    reaction_renderer.node_renderer.glyph = Square(
        size=8, fill_color="color")
    reaction_renderer.node_renderer.selection_glyph = Square(
        size=8, fill_color=Spectral4[2])
    reaction_renderer.node_renderer.hover_glyph = Square(
        size=8, fill_color=Spectral4[1])

    # nodes
    species_renderer.node_renderer.glyph = Circle(size=12, fill_color="color")
    species_renderer.node_renderer.selection_glyph = Circle(size=15, fill_color=Spectral4[2])
    species_renderer.node_renderer.hover_glyph = Circle(size=15, fill_color=Spectral4[1])
    
    #this part adds the interactive elements that make it so that the lines are highlighted 
    #when you mouse over and click
    edge_hover_tool = HoverTool(tooltips= None,renderers=[edges_renderer])
    if( not show_species_images):
        species_hover_tool = HoverTool(tooltips=[("name", "@species"), ("type", "@type")],\
                                        renderers=[species_renderer],attachment="right")
    else:
        species_hover_tool = HoverTool(tooltips='<div><div> <img src="data:image/png;base64,@image" style="float: left; margin: 0px 0px 0px 0px;"></img></div></div>',\
                                        renderers=[species_renderer],attachment="right")
    rxn_hover_tool = HoverTool(tooltips=[("reaction", "@species"), ("type", "@type"),("k_f","@k"),("k_r","@k_r")],\
                                        renderers=[reaction_renderer],attachment="right")
    
    plot.add_tools(edge_hover_tool,species_hover_tool,rxn_hover_tool, TapTool(), BoxSelectTool(),PanTool(),WheelZoomTool())

    edges_renderer.selection_policy = NodesAndLinkedEdges()
    edges_renderer.inspection_policy = EdgesAndLinkedNodes()

    plot.renderers.append(edges_renderer)
    plot.renderers.append(reaction_renderer)
    plot.renderers.append(species_renderer)

def generate_networkx_graph(CRN,useweights=False,use_pretty_print=False,pp_show_material=True,
                                                    pp_show_rates=True,pp_show_attributes=True, 
                                                    pp_show_compartments=True,
                                                colordict=None,imagedict = None):
    """generates a networkx DiGraph object that represents the CRN.
    input:
    ==========================
    CRN: a CRN from mixture.get_model() for example
    useweights: this will attempt to represent the reaction rates by the length of edges.
                short edges are fast rates. It doesn't look very good usually
    use_pretty_print: this uses the "pretty print" function to represent reactions and nodes a bit cleaner
    the next three parameters are pretty_print parameters
        pp_show_material: default false because this is listed in "type"
        pp_show_rates: default true because this is useful information
        pp_show_attributes
    colordict: a dictionary containing which node types are what color based upon the following keywords:
        Keywords are chosen to match species.material_type
            {"complex": "cyan",
            "protein": "green",
            "dna": "grey",
            "rna": "orange",
            "ligand": "pink",
            "phosphate": "yellow",
            "nothing":"purple"}

    When using a custom colordict, the following attributes will be checked to find colors with the first keys taking precedence:
        repr(species): "color"
        species.name: "color"
        (species.material_type, tuple(species.attributes)): "color"
        species.material_type: "color"
        tuple(species.attributes): "color"

    imagedict is a dictionary which contains species and their corresponding image representations.
    This is the output generated by CRNPlotter.renderMixture()

    output:
    ==================
    CRNgraph: the DiGraph object containing all nodes and edges
    CRNspeciesonly: a DiGraph object with only species
    CRNreactionsonly: a DiGraph object with only reactions
    """
    if(not PLOT_NETWORK):
        warn("network plotting disabled because some libraries are not found")
        return None, None, None
    if not colordict:
        colordict = {"complex": "cyan", "protein": "green",
                     "dna": "grey", "rna": "orange",
                     "ligand": "pink", "phosphate": "yellow", "nothing": "purple"}
    CRNgraph = nx.DiGraph()
    allnodenum = 1  # every node has an index
    # this starts at 1 because "nothing" is node 0
    nodedict = {}  # this is so that we can write out the reactions in
    # the reaction "species" field
    # it has {species:index}
    rxnlist = []  # list of numbers corresponding to only reaction nodes
    # sometimes reactions have no products. I want this to be represented in the graph with this
    # "nothing" node. However, usually we are making degradation reactions which yield the
    # degradation enzyme, so then it doesn't go to nothing. This means actually this node
    # isn't use for anything. But i think it's good to have just in case.
    default_species_color = "grey"
    default_reaction_color = "cornflowerblue"
    nodedict["nothing"] = 0
    CRNgraph.add_node(0)
    CRNgraph.nodes[0]["type"] = "nothing"
    CRNgraph.nodes[0]["species"] = "nothing"
    if("nothing" in colordict):
        CRNgraph.nodes[0]["color"] = colordict["nothing"]
    for species in CRN.species:
        # add all species first

        if repr(species) in colordict:
            species_color = colordict[repr(species)]
        elif species.name in colordict:
            species_color = colordict[species.name]
        elif (species.material_type, tuple(species.attributes)) in colordict:
            species_color = colordict[(species.material_type,
                               tuple(species.attributes))]
        elif(species.material_type in colordict):
            species_color = colordict[species.material_type]
        elif tuple(species.attributes) in colordict:
            species_color = colordict[tuple(species.attributes)]
        else:
            species_color = default_species_color

        nodedict[species] = allnodenum
        CRNgraph.add_node(allnodenum)
        CRNgraph.nodes[allnodenum]["type"]=str(species.material_type)
        if((imagedict is not None) and (species in imagedict)):
            CRNgraph.nodes[allnodenum]["image"]= imagedict[species].decode()
        if(not use_pretty_print):
            CRNgraph.nodes[allnodenum]["species"] = str(species)
        else:
            spectxt = species.pretty_print(
                show_material=pp_show_material, show_compartment=pp_show_compartments)
            CRNgraph.nodes[allnodenum]["species"] = spectxt
        CRNgraph.nodes[allnodenum]["color"] = species_color
        allnodenum += 1
    # reactions follow, allnodenum is not reset between these two loops
    for rxn in CRN.reactions:
        CRNgraph.add_node(allnodenum)
        CRNgraph.nodes[allnodenum]["type"] = str(rxn.propensity_type)
        if isinstance(rxn.propensity_type, MassAction):
            CRNgraph.nodes[allnodenum]["k"] = str(
                rxn.propensity_type.k_forward)
            CRNgraph.nodes[allnodenum]["k_r"] = str(
                rxn.propensity_type.k_reverse)
        else:
            CRNgraph.nodes[allnodenum]["k"] = str(rxn.propensity_type.k)
            CRNgraph.nodes[allnodenum]["k_r"] = ''

        reaction_color = None
        if repr(rxn) in colordict:
            reaction_color = colordict[repr(rxn)]
        for k,p in rxn.propensity_type.propensity_dict["parameters"].items():
            if(hasattr(p,"search_key")):
                mech_str = repr(p.search_key.mechanism).strip('\'\"')
                partid_str = repr(p.search_key.part_id).strip('\'\"')
                name_str = repr(p.search_key.name).strip('\'\"')
                new_color = reaction_color
                if(name_str in colordict):
                    new_color = colordict[name_str] #name of the mechanism that made the reaction
                if(partid_str in colordict):
                    new_color = colordict[partid_str] #partid used to make the reaction
                if(mech_str in colordict):
                    new_color = colordict[mech_str] #the type of mechanism that made the reaction
                if(reaction_color is not None and reaction_color != new_color):
                        #if you change the color of a reaction, the final color will be unexpected,
                        #so there is a warning that that happened
                        warn(f"reaction color was {reaction_color} but now you want it to be {colordict[name_str]}")
                reaction_color = new_color


        # CRNgraph.nodes[allnodenum]
        if isinstance(rxn.propensity_type, MassAction):
            kval = rxn.propensity_type.k_forward
            CRNgraph.nodes[allnodenum]["k"] = str(kval)
        else:
            kval = rxn.k
            CRNgraph.nodes[allnodenum]["k"] = str(rxn.propensity_type.k)

        if(not useweights):
            kval = 1
        if isinstance(rxn.propensity_type, MassAction):
            krev_val = rxn.propensity_type.k_reverse
        else:
            krev_val = None
        if((krev_val is not None) and (not useweights)):
            krev_val = 1
        for reactant in rxn.inputs:
            
            CRNgraph.add_edge(nodedict[reactant.species],allnodenum,weight=kval)
            if(krev_val is not None):
                # if the k is 0 then the node does not exist, right?
                CRNgraph.add_edge(
                    allnodenum, nodedict[reactant.species], weight=krev_val)
        for product in rxn.outputs:
            #TODO species cannot find another species in the nodedict????
            CRNgraph.add_edge(allnodenum,nodedict[product.species],weight=kval)
            if(krev_val is not None):
                CRNgraph.add_edge(
                    nodedict[product.species], allnodenum, weight=krev_val)
        if(len(rxn.outputs) == 0):
            # this adds an edge to the "nothing" node we made in the beginning
            CRNgraph.add_edge(allnodenum, 0, weight=kval)
            if(krev_val is not None):
                CRNgraph.add_edge(0, allnodenum, weight=krev_val)
        elif(len(rxn.inputs) == 0):
            # this adds an edge from the "nothing" node we made in the beginning
            CRNgraph.add_edge(0, allnodenum, weight=kval)
            if(krev_val is not None):
                CRNgraph.add_edge(allnodenum, 0, weight=krev_val)
        if(reaction_color is None):
            reaction_color = default_reaction_color
        CRNgraph.nodes[allnodenum]["color"] = reaction_color
        if(not use_pretty_print):
            CRNgraph.nodes[allnodenum]["species"] = str(rxn)
        else:
            rxntxt = rxn.pretty_print(
                show_material=pp_show_material, show_rates=pp_show_rates, show_attributes=pp_show_attributes)
            # this will show up as "reaction" in the tool tip
            CRNgraph.nodes[allnodenum]["species"] = rxntxt
        # the name of the reaction is the string representation
        rxnlist += [allnodenum]
        allnodenum += 1
    CRNspeciesonly = CRNgraph.copy()
    CRNspeciesonly.remove_nodes_from(rxnlist)
    CRNreactionsonly = CRNgraph.copy()
    CRNreactionsonly.remove_nodes_from(range(rxnlist[0]))
    return CRNgraph, CRNspeciesonly, CRNreactionsonly


class CRNPlotter:
    class MultiPart:
        def __init__(self,name,parts_list,bound = None):
            """multiple simple parts which are treated as one"""
            self.name = name
            self.parts_list = parts_list
            self.bound = None
        def get_directed(self,direction,bound=None,non_binder=None):

            if(non_binder is None):
                non_binder = ["Promoter"]
            new_multipart = copy.deepcopy(self)
            bound_for_distribution = None
            if(bound is not None):
                bound_for_distribution = copy.copy(bound)
            elif(self.bound is not None):
                bound_for_distribution = copy.copy(self.bound)
            if(bound_for_distribution is not None):
                recursion = 10
                while(len(bound_for_distribution)>0 and recursion > 0):
                    for part in new_multipart.parts_list:
                        if(part.dpl_type not in non_binder):
                            if(part.bound is None):
                                part.bound = [bound_for_distribution.pop(0)]
                            else:
                                part.bound += [bound_for_distribution.pop(0)]
                            if(len(bound_for_distribution) == 0):
                                break
                    recursion -= 1
                if(recursion == 0):
                    raise ValueError(f"reached maximum recursion when trying to populate multipart {self}")
                
            new_multipart.parts_list = [a.get_directed(direction) for a in new_multipart.parts_list]
            if(direction=="reverse"):
                new_multipart.parts_list = new_multipart.parts_list[::-1]
                return new_multipart
            else:
                return new_multipart
        def get_dpl(self):
            outlist = []
            for part in self.parts_list:
                outlist+= part.get_dpl()
            return outlist
        def __repr__(self):
            return "MultiPart("+",".join([str(a) for a in self.parts_list])+")"
          
    class SimpleConstruct:
        def __init__(self,name,parts_list,circular=False,material_type="dna",label_size=13,added_opts = None):
            self.name = name
            self.parts_list = parts_list
            self.circular = circular
            self.material_type = material_type
            self.label_size = label_size
            if(added_opts is None):
                self.added_opts = {}
        def get_dpl(self):
            outlist = []
            for part in self.parts_list:
                part_dpl = part.get_dpl()
                for subpart_dpl in part_dpl:
                    subpart_dpl["opts"].update(self.added_opts)
                outlist+= part_dpl
            return outlist
        def get_dpl_binders(self):
            my_dpl_output = self.get_dpl()
            out_regs = []
            for design in my_dpl_output:
                if('make_binders' in design):
                    for binder in design['make_binders']:
                        linecolor = 'blue'
                        if(binder.material_type=='dna'):
                            linecolor = 'black'
                        elif(binder.material_type=='rna'):
                            linecolor = 'red'
                        if(hasattr(binder,"color")):
                            bindcolor = binder.color
                        else:
                            bindcolor = (0.5,0.5,0.5)
                        out_reg = {'type':'Binding', 'from_part':design, 'to_part':design,
                                            'opts':{'label':binder.name,'label_size':self.label_size*.7,\
                                                    'color':linecolor, 'label_x_offset':-1,'y_offset':10,\
                                                    'face_color':bindcolor}
                                  }
                        out_reg['opts'].update(binder.added_opts)
                        out_regs+=[out_reg]
            return out_regs,my_dpl_output
        def renderDNA(self,  dna_renderer,ax=None,plot_backbone=True):
            
            part_renderers = dna_renderer.SBOL_part_renderers()
            reg_renderers = dna_renderer.std_reg_renderers()
            if(ax is None):
                figsize = (1,1)
                #fig,ax = plt.subplots(constrained_layout=True)
                fig = plt.figure(figsize=figsize)
                ax = fig.add_axes([0,0,1,1])
                plt.tight_layout(pad=0.0001)
            
            my_regs,my_designs = self.get_dpl_binders()
            for part in my_designs:
                part['opts'].update({'edgecolor':dna_renderer.linecolor})
            start,end = dna_renderer.renderDNA(ax,my_designs,part_renderers,circular=self.circular,\
                                    regs=my_regs, reg_renderers=reg_renderers,plot_backbone=plot_backbone)
            
            fig = ax.get_figure()
            ax.axis('off')
            ''
            ylimits = [None,None]
            xlimits = [None,None]
            relevant_stuff = ax.patches+ax.texts
            for patch in relevant_stuff:
                bbox = patch.get_window_extent(renderer=fig.canvas.get_renderer())
                if(ylimits == [None,None]):
                    ylimits = [bbox.ymin,bbox.ymax]
                if(bbox.ymax > ylimits[1]):
                    ylimits[1] = bbox.ymax
                if(bbox.ymin < ylimits[0]):
                    ylimits[0] = bbox.ymin
                if(xlimits == [None,None]):
                    xlimits = [bbox.xmin,bbox.xmax]
                if(bbox.xmax > xlimits[1]):
                    xlimits[1] = bbox.xmax
                if(bbox.xmin < xlimits[0]):
                    xlimits[0] = bbox.xmin
            xlimits[0],ylimits[0] = fig.transFigure.inverted().transform((xlimits[0], ylimits[0]))

            xlimits[1],ylimits[1] = fig.transFigure.inverted().transform((xlimits[1], ylimits[1]))
            ax.relim()
            yheight = ylimits[1]-ylimits[0]
            xheight = xlimits[1]-xlimits[0]
            fig.set_size_inches(xheight/24,yheight/24)
            ax.set_aspect('equal')
            ax.autoscale_view()
            #ax.set_xlim((xlimits[0],xlimits[1]))
            #ax.set_ylim((ylimits[0],ylimits[1]))
            return ax
    
    class SimplePart:
        def __init__(self,name,dpl_type,direction='forward',bound=None,color=None,\
                        color2=None,show_label=True,label_size = 13,label_y_offset = -8,added_opts=None,material_type=None):
            """a simple 'part' for the sole purpose of rendering using DNAplotlib"""
            self.name = name
            self.color = color
            self.color2 = color2 
            self.dpl_type = dpl_type #this is a string which dnaplotlib knows about
            self.direction=direction #"forward" or "reverse"
            self.bound = bound #this should be a list of SimpleParts
            if(added_opts is None):
                self.added_opts = {}
            else:
                self.added_opts = added_opts #dictionary of keywords for dnaplotlib
            self.show_label =show_label #if the label should be added to the 'opts' upon output
            self.label_size = label_size #font size of the label
            self.label_y_offset = label_y_offset
            self.material_type = material_type
        def get_directed(self,direction,bound=None):
            copied_part = copy.copy(self)
            copied_part.direction = direction
            if(bound is not None):
                copied_part.bound=bound
            elif(self.bound is not None):
                copied_part.bound = self.bound
            return copied_part
        def get_dpl(self,bound=None):
            direction = True
            if(self.direction == "reverse"):
                direction = False
            dpl_out = {'type':self.dpl_type, 'name':self.name, 'fwd':direction}
            opts = {'color':self.color,'color2':self.color2}
            if(self.show_label):
                opts['label']=self.name
                opts['label_size']=self.label_size
                opts['label_y_offset']=self.label_y_offset
            opts.update(self.added_opts)
            dpl_out['opts']=opts
            if(bound is not None):
                dpl_out['make_binders']=bound
            elif(self.bound is not None):
                dpl_out['make_binders']=self.bound

            return [dpl_out]
        def __repr__(self):
            return "SimplePart("+str(self.name)+"-"+str(self.direction)[0]+")"

    def __init__(self,dna_renderer=None,rna_renderer=None,cmap = "Set3"):
        if(dna_renderer is None):
            self.dna_renderer=dpl.DNARenderer(scale = 5,linewidth=3)
        else:
            self.dna_renderer = dna_renderer
        if(rna_renderer is None):
            self.rna_renderer=dpl.DNARenderer(scale = 5,linewidth=3,linecolor=(1,0,0))
        else:
            self.rna_renderer = rna_renderer
        self.cmap = plt.get_cmap(cmap).colors
        self.color_counter = 0
        self.clear_dicts()
    def renderMixture(self,mixture,crn = None,rna_renderer=None,dna_renderer=None,store=True):
        """creates dnaplotlib images for all relevant species in a mixture"""
        if(crn is None):
            mycrn = mixture.compile_crn()
        else:
            mycrn = crn
        if(rna_renderer is None and self.rna_renderer is not None):
            rna_renderer = self.rna_renderer
        else:
            raise ValueError("rna_renderer cannot be None")
        if(dna_renderer is None and self.dna_renderer is not None):
            dna_renderer = self.dna_renderer
        else:
            raise ValueError("dna_renderer cannot be None")
        self.clear_dicts()
        for component in mixture.component_enumeration():
            if(isinstance(component,Construct)):
                a = self.make_dpls_from_construct(component)
        plt.ioff()
        for species in mycrn.species:
            a = self.make_dpl_from_species(species)
            if(a.material_type is not None):
                plot_bb = True
                if(not isinstance(a,self.SimpleConstruct)):
                    plot_bb = False
                    newcon = self.SimpleConstruct(a.name,[a],material_type=a.material_type)
                else:
                    newcon = a
                
                if(newcon.material_type=='rna'):
                    ax = newcon.renderDNA(rna_renderer,plot_backbone=plot_bb)
                else:
                    ax = newcon.renderDNA(dna_renderer,plot_backbone=plot_bb)
                if(store):
                    imagestream = io.BytesIO()
                    fig = ax.get_figure()
                    fig.savefig(imagestream,bbox_inches='tight')
                    png_str = base64.b64encode(imagestream.getvalue())
                    self.species_image_dict[species]= png_str
        return self.species_image_dict
    def renderConstruct(self,construct_obj,rna_renderer=None,dna_renderer=None,showlabels=True, render_rna=True):
        if(rna_renderer is None and self.rna_renderer is not None):
            rna_renderer = self.rna_renderer
        else:
            raise ValueError("rna_renderer cannot be None")
        if(dna_renderer is None and self.dna_renderer is not None):
            dna_renderer = self.dna_renderer
        else:
            raise ValueError("dna_renderer cannot be None")
        a = self.make_dpls_from_construct(construct_obj)
        outaxs = [a.renderDNA(dna_renderer)]
        if(render_rna):
            components = construct_obj.enumerate_constructs()
            for component in components:
                if(isinstance(component,RNA)):
                    a = self.make_dpls_from_construct(component)
                    outaxs += [a.renderDNA(rna_renderer)]
        return outaxs
    def clear_dicts(self):
        self.part_dpl_dict = {}
        self.construct_dpl_dict = {}
        self.species_dpl_dict = {}
        self.species_image_dict ={}
    def get_color(self):
        if(self.cmap is not None):
            out_color = self.cmap[self.color_counter]
            self.color_counter+=1
            if(self.color_counter>=len(self.cmap)):
                self.color_counter = 0
            return out_color
        else:
            raise ValueError("No colormap set")
    def make_dpls_from_construct(self,construct,save_in_dict=True):
        if(construct in self.construct_dpl_dict):
            return self.construct_dpl_dict[construct]
        else:
            new_parts_list = []
            for part in construct:
                new_parts_list += [self.make_dpl_from_part(part,save_in_dict=save_in_dict)]
            
            if(isinstance(construct,DNA)):
                mat_type = "dna"
            elif(isinstance(construct,RNA)):
                mat_type = "rna"
            simple_construct = self.SimpleConstruct(name = construct.name,\
                                                    parts_list=new_parts_list,\
                                                    circular=construct.circular,\
                                                    material_type = mat_type)
            if(save_in_dict):
                self.construct_dpl_dict[construct]=simple_construct
            return simple_construct
    
    def make_dpl_from_species(self,species):
        if(species in self.species_dpl_dict):
            return self.species_dpl_dict[species]
        else:
            if(isinstance(species,OrderedPolymer)):
                #we are dealing with a polymer
                polylist = [] #accumulating the list of SimpleParts
                for monomer in species:
                    #going through the monomers
                    removed_monomer = monomer.get_removed()
                    if(removed_monomer in self.species_dpl_dict):
                        #if we already know about this monomer, just use that
                        polylist += [self.species_dpl_dict[removed_monomer].get_directed(monomer.direction)]
                    elif(isinstance(monomer,ComplexSpecies)):
                        #if the monomer is a complex that means we have to make a bound simplepart
                        binders = []
                        base_simplepart = None #we must figure out who the base is. This is how we do that
                        for specie in monomer.get_species(recursive=True):
                            if(isinstance(specie,ComplexSpecies)):
                                continue
                            if(specie.material_type=='part'):
                                #this material type is only made by dna constructs
                                base_simplepart = copy.copy(self.make_dpl_from_species(specie)) #copy it because now we make it bound to stuff
                                #ideally you already ran make_dpl_from_construct and so this will be populated
                            else:
                                binders += [self.make_dpl_from_species(specie)] #other stuff gets added as a binder
                        base_simplepart.bound = binders
                        self.species_dpl_dict[removed_monomer]=base_simplepart
                        polylist += [base_simplepart.get_directed(monomer.direction)]
                    elif(isinstance(monomer,Species)):
                        #in this case we couldnt find the monomer in our dictionary, but it is
                        #at the base level
                        base_simplepart = self.SimplePart(monomer.name,'UserDefined',color=self.get_color()) #new simplepart
                        self.species_dpl_dict[removed_monomer]=base_simplepart #adding it to the dictionary for the future
                        polylist += [base_simplepart.get_directed(monomer.direction)] #add it to the polymer, with a direction
                out_dpl = self.SimpleConstruct(name = species.name,parts_list = polylist,circular=species.circular,material_type = species.material_type)
            elif(isinstance(species,ComplexSpecies)):
                #if we have a complex but it is not a polymer, we just do the "binding" determination of the polymer part
                base_simplepart = None
                binders = []
                for specie in species.species:
                    #i don't really care who is the base of this thing
                    if(base_simplepart is None):
                        base_simplepart = self.make_dpl_from_species(specie)
                    else:
                        binders += [self.make_dpl_from_species(specie)]
                
                out_dpl = base_simplepart.get_directed(None,binders)
            elif(isinstance(species,Species)):
                #this happens if we have another type of species
                dpl_type = 'Origin' #default is 'origin', which is just a circle
                if(species.material_type=='dna'):
                    dpl_type = 'UserDefined'
                elif(species.material_type=='rna'):
                    dpl_type = 'NCRNA'
                out_dpl = self.SimplePart(name = species.name,dpl_type = dpl_type,color=self.get_color(),material_type=species.material_type)
            
            self.species_dpl_dict[species] = out_dpl
            return out_dpl
    def make_dpl_from_part(self,part,set_color=None,save_in_dict=True,set_color2=None):
        removed_part = part.get_removed()
        if(removed_part in self.part_dpl_dict):
            
            return self.part_dpl_dict[removed_part].get_directed(part.direction)
        else:
            dpl_type= "UserDefined"
            needs_color2 = False
            regs = None
            if(isinstance(part,Promoter)):
                dpl_type = "Promoter"
                if(hasattr(part,"regulators")):
                    regs = part.regulators
            elif(isinstance(part,RBS)):
                dpl_type = "RBS"
            elif(isinstance(part,CDS)):
                dpl_type = "CDS"
            elif(isinstance(part,Protein)):
                dpl_type = "CDS"
            elif(isinstance(part,Terminator)):
                dpl_type = "Terminator"
            elif(isinstance(part,IntegraseSite)):
                if(part.site_type in ["attP","attB"]):
                    dpl_type = "RecombinaseSite"
                elif(part.site_type in ["attL","attR"]):
                    dpl_type = "RecombinaseSite2"
                    needs_color2 = True
                for key_part in self.part_dpl_dict:
                    stored_part = self.part_dpl_dict[key_part]
                    if(isinstance(key_part,IntegraseSite) and \
                        key_part.integrase==part.integrase and \
                            key_part.dinucleotide==part.dinucleotide):
                        types_list = [key_part.site_type,part.site_type]
                        if(types_list in [["attP","attL"],["attL","attR"],["attR","attL"], ["attB","attR"]]):
                            if(set_color2 is None):
                                set_color2 = stored_part.color
                        if(types_list in [["attB","attL"],["attP","attR"],["attL","attB"],\
                                                                                ["attR","attP"]]):
                            if(set_color is None):
                                set_color = stored_part.color
                        elif(types_list in [["attL","attR"],["attR","attL"],\
                                                                ["attL","attP"],["attR","attB"]]):
                            if(set_color is None):
                                set_color = stored_part.color2
                    if(set_color is not None and set_color2 is not None):
                        #we found all the needed colors so give up
                        break
            if(set_color is None):
                color = self.get_color()
            else:
                color = set_color
            color2 = None
            if(set_color2 is not None):
                color2 = set_color2
            elif(set_color2 is None and needs_color2):
                color2 = self.get_color()
            
            
            outpart = self.SimplePart(name=part.name,\
                                      dpl_type=dpl_type,
                                      color=color,
                                      color2=color2)
            if(regs is not None):
                regparts = []
                for reg in regs:
                    regpart = self.SimplePart(name=reg,dpl_type="Operator",color=color,show_label=False)
                    regparts += [regpart]
                retpart = self.MultiPart(name=part.name,parts_list =[outpart]+regparts)
            else:
                retpart = outpart
            if(save_in_dict):
                self.part_dpl_dict[removed_part] = retpart
                self.species_dpl_dict[part.dna_species] = retpart
            return retpart.get_directed(part.direction)
