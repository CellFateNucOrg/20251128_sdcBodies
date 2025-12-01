import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.image as mpimg

# === CONFIG ===
table_file = "/Volumes/external.data/MeisterLab/jsemple/20251118_sdcBodies/ppts/rnaqc.xlsx"   # Replace with your table file (CSV/Excel)
plots_dir = "/Volumes/external.data/MeisterLab/jsemple/20251118_sdcBodies/NovogeneQC_sdcBody_coh1/qc_reports_novogene/SEU20251021117_X208SC25105402-Z01_NVDE2025100219-CH-BernUni-Semple-42-Celegans-mRNASeq-6Gb-WOBI-__20251031102809543/src/Caliper"        # Directory containing .png files
output_pdf = "RNAqcPlots.pdf"  # Final multipage PDF

# === LOAD TABLE ===
# If Excel: use pd.read_excel("table.xlsx")
if table_file.endswith('.csv'):
    df = pd.read_csv(table_file)
elif table_file.endswith('.tsv'):
    df = pd.read_csv(table_file, sep='\t')
elif table_file.endswith('.xlsx') or table_file.endswith('.xls'):
    df = pd.read_excel(table_file)
else:
    raise ValueError("Unsupported table file format. Use .csv or .xlsx/.xls")

filename_col="Nucleic Acid ID"  # Column with IDs matching filenames
filetype=".PNG"               # File extension of plots
title_col="Sample Name"        # Column with titles for each plot
num_rows = 4
num_cols = 2

# === CREATE PDF ===
with PdfPages(output_pdf) as pdf:
    # Process in chunks of 8 (4x2 grid)
    for i in range(0, len(df), num_rows * num_cols):
        chunk = df.iloc[i:i+(num_rows * num_cols)]
        
        fig, axes = plt.subplots(num_rows, num_cols, figsize=(11, 14))  # A4-ish size
        axes = axes.flatten()
        
        for ax, (_, row) in zip(axes, chunk.iterrows()):
            nuc_id = str(row[filename_col])
            title = str(row[title_col])
            
            # Find matching file in plots_dir
            matched_files = [f for f in os.listdir(plots_dir) if nuc_id in f and f.endswith(filetype)]
            
            if matched_files:
                img_path = os.path.join(plots_dir, matched_files[0])
                img = mpimg.imread(img_path)
                ax.imshow(img)
                ax.set_title(title, fontsize=10)
                ax.axis("off")
            else:
                ax.text(0.5, 0.5, f"No file for {nuc_id}", ha="center", va="center")
                ax.set_title(title, fontsize=10)
                ax.axis("off")
        
        # Hide unused subplots if chunk < 8
        for j in range(len(chunk), num_rows * num_cols):
            axes[j].axis("off")
        
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

print(f"✅ PDF created: {output_pdf}")
