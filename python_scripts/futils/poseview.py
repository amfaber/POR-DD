import nglview as nv
import ipywidgets
from ipywidgets import widgets
import plotly.graph_objects as go
from functools import partial
from plotly.subplots import make_subplots

def create_widget(biopython_protein, ligs, cnndata = None, title = None):
    wid = nv.NGLWidget()
    wid.layout.width = "750px"
    wid.layout.height = "500px"
    wid.add_structure(nv.BiopythonStructure(biopython_protein))

    ligcomponent = None
    def update_lig(change):
        nonlocal ligcomponent
        idx = change.new
        if ligcomponent is not None:
            wid.remove_component(ligcomponent)
        ori = wid._camera_orientation
        ligcomponent = wid.add_component(ligs[idx], no_zoom = True)
        wid._set_camera_orientation(ori)
    
    def increment_widget(but, wid, amount = 1):
        wid.value += amount
    
    def update_fig(change, fig):
        idx = change.new
        fig.data[0].y = [cnndata.loc[idx, "CNNscore"]]
        fig.data[1].y = [cnndata.loc[idx, "CNNaffinity"]]
    
    i = widgets.IntText(value = 1, layout = widgets.Layout(width = "100px"))

    if cnndata is not None:
        trace = go.Bar(width = 0.4)
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        fig = go.FigureWidget(data = fig, layout_margin = dict(r = 20, l = 20, b = 20, t = 20), 
                            layout_width = 300, layout_height = 500)
        
        trace.x = ["CNNscore"]
        trace.marker.color = "#1f77b4"
        fig.add_trace(trace)
        trace.x = ["CNNaffinity"]
        trace.marker.color = "#ff7f0e"
        fig.add_trace(trace, secondary_y=True)
        fig.update_yaxes(range = [0, 1], secondary_y=False)
        fig.update_yaxes(range = [2, 7], secondary_y=True)
        fig.update_layout(showlegend = False)
        i.observe(partial(update_fig, fig = fig), names = "value")
    i.observe(update_lig, names = "value")
    i.value = 0
    nextbutton = widgets.Button(description = "Next")
    nextbutton.on_click(partial(increment_widget, wid = i))
    prevbutton = widgets.Button(description = "Previous")
    prevbutton.on_click(partial(increment_widget, wid = i, amount = -1))

    controls = ipywidgets.HBox([i, prevbutton, nextbutton])
    molview = ipywidgets.VBox([controls, wid])
    if cnndata is not None:
        box = ipywidgets.HBox([molview, fig])
    else:
        box = molview
    if title is not None:
        title = widgets.HTML(f'<p style="font-size:2.5em">{title}</p>')
        titlebox = widgets.VBox([box, title])
        return titlebox
    return box