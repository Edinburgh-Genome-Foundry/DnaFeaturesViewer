from A_linear_plot import CustomTranslator

translator = CustomTranslator()
graphic_record = translator.translate_record("./plasmid.gb")
graphic_record.plot_on_multiple_pages(
    "multiline_plot.pdf",
    nucl_per_line=80,
    lines_per_page=10,
    plot_sequence=True
)
