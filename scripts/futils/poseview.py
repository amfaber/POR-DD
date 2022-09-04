import nglview as nv
import ipywidgets
from ipywidgets import widgets
import plotly.graph_objects as go
from functools import partial
from plotly.subplots import make_subplots
from matplotlib.colors import TABLEAU_COLORS

class Poseviewer:
    def __init__(self, biopython_protein, df, lig_cols, datacols = None, title = None, colors = None):
        self.widget = self.create_widget(biopython_protein, df, lig_cols, datacols = datacols, title = title, colors = colors)

    def _repr_pretty_(self):
        return self.widget

    def create_widget(self, biopython_protein, df, lig_cols, datacols = None, title = None, colors = None):
        self.wid = nv.NGLWidget()
        self.wid.layout.width = "750px"
        self.wid.layout.height = "500px"
        self.wid.add_structure(nv.BiopythonStructure(biopython_protein))

        if colors is None:
            self.colors = list(TABLEAU_COLORS.values())
        else:
            self.colors = colors

        self.ligcomponents = None
        def update_lig(change):
            idx = change.new
            if self.ligcomponents is not None:
                for ligcomponent in self.ligcomponents:
                    self.wid.remove_component(ligcomponent)
            ori = self.wid._camera_orientation
            self.ligcomponents = []
            for i, lig_col in enumerate(lig_cols):
                self.ligcomponents.append(self.wid.add_component(df.loc[idx, lig_col]))
                self.ligcomponents[i].clear()
                self.ligcomponents[i].add_ball_and_stick(colorValue = self.colors[i%len(self.colors)])
            self.wid._set_camera_orientation(ori)
        
        def increment_widget(but, wid, amount = 1):
            wid.value += amount
        
        def update_fig(change, fig):
            idx = change.new
            fig.data[0].y = [df.loc[idx, "CNNscore"]]
            fig.data[1].y = [df.loc[idx, "CNNaffinity"]]
        
        i = widgets.IntText(value = 1, layout = widgets.Layout(width = "100px"))

        if datacols is not None:
            trace = go.Bar(width = 0.4)
            self.fig = make_subplots(specs=[[{"secondary_y": True}]])
            self.fig = go.FigureWidget(data = self.fig, layout_margin = dict(r = 20, l = 20, b = 20, t = 20), 
                                layout_width = 300, layout_height = 500)
            
            trace.x = ["CNNscore"]
            trace.marker.color = "#1f77b4"
            self.fig.add_trace(trace)
            trace.x = ["CNNaffinity"]
            trace.marker.color = "#ff7f0e"
            self.fig.add_trace(trace, secondary_y=True)
            self.fig.update_yaxes(range = [0, 1], secondary_y=False)
            self.fig.update_yaxes(range = [2, 8], secondary_y=True)
            self.fig.update_layout(showlegend = False)
            i.observe(partial(update_fig, fig = self.fig), names = "value")
        i.observe(update_lig, names = "value")
        i.value = 0
        nextbutton = widgets.Button(description = "Next")
        nextbutton.on_click(partial(increment_widget, wid = i))
        prevbutton = widgets.Button(description = "Previous")
        prevbutton.on_click(partial(increment_widget, wid = i, amount = -1))

        controls = ipywidgets.HBox([i, prevbutton, nextbutton])
        molview = ipywidgets.VBox([controls, self.wid])
        if datacols is not None:
            box = ipywidgets.HBox([molview, self.fig])
        else:
            box = molview
        if title is not None:
            title = widgets.HTML(f'<p style="font-size:2.5em">{title}</p>')
            titlebox = widgets.VBox([box, title])
            return titlebox
        return box