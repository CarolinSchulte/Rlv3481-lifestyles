
import cobra
import riptide
import csv


model = cobra.io.read_sbml_model('rhizo_riptide.xml')

transcriptome_rhizo = riptide.read_transcription_file('Data/rhizosphere_WT_rpkm.txt',norm=False)

with open('Data/rhizosphereGenes.csv', newline='',encoding='utf-8-sig') as f:
    reader = csv.reader(f)
    data = list(reader)

[rhizo_genes] = data


riptide_maxfit_WT_full = riptide.maxfit_contextualize(model,transcriptome=transcriptome_rhizo,\
                                                      tasks=rhizo_genes,frac_min=0.5,frac_max=0.95)

riptide.save_output(riptide_maxfit_WT_full,path='Results/rhizosphere')
