import urllib.request as urlreq
import dash
from dash.dependencies import Input, Output
import dash_bio as dashbio
from dash import html, dcc
from dash_bio.utils import xyz_reader

app = dash.Dash(__name__)


data = urlreq.urlopen(
    'https://git.io/speck_methane.xyz'
).read().decode('utf-8')

data = xyz_reader.read_xyz(datapath_or_datastring=data, is_datafile=False)

app.layout = html.Div([
    dcc.Dropdown(
        id='default-speck-preset-views',
        options=[
            {'label': 'Default', 'value': 'default'},
            {'label': 'Ball and stick', 'value': 'stickball'}
        ],
        value='default'
    ),
    dashbio.Speck(
        id='default-speck',
        data=data
    ),
])

@app.callback(
    Output('default-speck', 'presetView'),
    Input('default-speck-preset-views', 'value')
)
def update_preset_view(preset_name):
    return preset_name


if __name__ == '__main__':
    app.run_server(debug=True)
