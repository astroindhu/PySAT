from bokeh.plotting import figure, show, output_file
from bokeh.palettes import d3
from bokeh.layouts import column
from bokeh.models import Legend


def bokehplot(xy,keylist, ptitle, xlabel, ylabel,lstyle ='solid',lwidth=5):
    '''customised bokeh plot supporting psl files'''
    fig = figure(title= ptitle, x_axis_label=xlabel, y_axis_label=ylabel, plot_width=1000, plot_height=700,toolbar_location='above')
    fig.title.text_font = "times"
    fig.title.text_font_style = "bold"
    fig.title.text_font_size = '24pt'
    fig.axis.axis_label_text_font_size = '24pt'
    fig.axis.major_label_text_font_size = '24pt'
    fig.axis.major_label_text_font_style = 'bold'
    fig.axis.axis_label_text_font_style = 'bold'
    fig.xgrid.grid_line_color = None
    fig.ygrid.grid_line_color = None
    fig.outline_line_width = 1
    fig.outline_line_color = "black"

    legend_it = []
    clibrary = d3["Category10"][10]+['#000000','#440154', '#30678D','#30678D','#DD4968']
    clist = clibrary[:len(keylist)]
    cl = 0
    for l in sorted(keylist): 
        plot = fig.line(xy[l][:,0], xy[l][:,1], color=clist[cl],line_width=lwidth, line_dash=lstyle, muted_color=clist[cl],muted_alpha=0.3)
        cl+=1
        legend_it.append((l,[plot]))
    legend = Legend(items=legend_it, location=(5, -60))
    legend.click_policy="mute"
    legend.label_text_font_size = '20pt'
    fig.add_layout(legend, 'right')
    # fig.legend.background_fill_alpha = 0
    return fig