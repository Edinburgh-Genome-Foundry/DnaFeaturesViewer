from A_linear_plot import CustomTranslator

translator = CustomTranslator()
graphic_record = translator.translate_record("plasmid.gb")
cropped_record = graphic_record.crop((6320, 6350))
cropped_record.ticks_resolution = 5  # One tick every 5 nucleotides
ax, _ = cropped_record.plot(figure_width=6)
cropped_record.plot_sequence(ax, guides_intensity=0.2)
cropped_record.plot_translation(
    ax=ax, location=(6335, 6350, -1), fontdict={"weight": "bold"}
)
ax.figure.savefig("B_detail_plot.svg", bbox_inches="tight")
