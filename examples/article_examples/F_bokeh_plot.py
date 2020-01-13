from bokeh.embed import file_html
from bokeh.resources import CDN
from A_linear_plot import CustomTranslator

translator = CustomTranslator()
graphic_record = translator.translate_record("plasmid.gb")
bokeh_plot = graphic_record.plot_with_bokeh(figure_width=10, figure_height=2)
html = file_html(bokeh_plot, CDN, "my plot")
with open("F_bokeh_plot.html", "w") as f:
    f.write(html)
